MODULE TASK_B
        IMPLICIT NONE

        INTEGER, PRIVATE :: i, j
        REAL(KIND=8), PRIVATE :: pi = 4d0 * DATAN(1d0)
        COMPLEX(KIND=8), PRIVATE :: im = (0d0, 1d0)
         
        CONTAINS
                FUNCTION FFT(x, V, N) 
                        ! Computes the discrete Fourier transform of the function V
                        ! through the FFTW routine (http://www.fftw.org)

                        IMPLICIT NONE

                        INCLUDE 'fftw3.f'

                        ! INPUT
                        INTEGER, INTENT(IN) :: N ! dimension of V
                        REAL(KIND=8), DIMENSION(N), INTENT(IN) :: V, x
                      
                        ! OUTPUT
                        COMPLEX(KIND=8), DIMENSION(N/2 + 1) :: FFT ! Fourier transfrom of V

                        ! ROUTINE PARAMETERS
                        INTEGER(KIND=8) :: plan
                        
                        REAL(KIND=8) :: k(N/2 +1), L
                        REAL(KIND=8), DIMENSION(N) :: IFT   ! Inverse Fourier Transform

                        CALL DFFTW_PLAN_DFT_R2C_1D(plan, N, V, FFT, FFTW_ESTIMATE)
                        CALL DFFTW_EXECUTE_DFT_R2C(plan, V, FFT)
                        CALL DFFTW_DESTROY_PLAN(plan)

                        FFT = FFT / DBLE(N)

                        CALL DFFTW_PLAN_DFT_C2R_1D(plan, N, FFT, IFT, FFTW_PRESERVE_INPUT, FFTW_ESTIMATE)
                        CALL DFFTW_EXECUTE_DFT_C2R(plan, FFT, IFT)
                        CALL DFFTW_DESTROY_PLAN(plan)

                        L = x(N) - x(1)
                        k(:) = (/ ((2d0 * i) * pi / L, i=0,N/2) /)

                        OPEN(UNIT=10, FILE='fourier.txt')
                          WRITE(10, '(6F15.7)') (x(i), V(i), k(i), FFT(i), IFT(i), i=1,N/2+1)
                          WRITE(10, '(6F15.7)') (x(i), V(i), k(N/2+1), FFT(N/2+1), IFT(i), i = N/2+2, N)
                        CLOSE(10)

                END FUNCTION FFT

                FUNCTION HAMILTONIAN(V, N, k) RESULT(H)
                        ! Returns the Hamiltonian (the upper triangle) obtained 
                        ! after expanding the wave function on the plane 
                        ! waves basis set.

                        IMPLICIT NONE

                        ! INPUT
                        INTEGER, INTENT(IN) :: N ! order of H
                        REAL(KIND=8), DIMENSION(N), INTENT(IN) :: k   ! wave vectors
                        COMPLEX(KIND=8), DIMENSION(:), INTENT(IN) :: V ! Fourier transform of the potential

                        ! OUTPUT
                        COMPLEX(KIND=8), DIMENSION(N, N) :: H 

                        H(:,:) = 0d0
                        DO i = 1,N
                          H(i,i) = V(1) + 0.5 * k(i)**2
                          H(i,i+1:) = V(2:N+1-i)
                        END DO

                END FUNCTION HAMILTONIAN
          
                SUBROUTINE SOLVE_EIGH(x, N, V, eigval_only, eigh_range, IL, IU, VL, VU, store)
                        ! Computes the eigenvalues and the eigenvectors of an
                        ! upper triangle hermitian matrix H by calling the LAPACK 
                        ! routine ZHEEVX (http://www.netlib.org/lapack)

                        IMPLICIT NONE

                        INCLUDE 'fftw3.f'

                        ! INPUT
                        INTEGER, INTENT(IN) :: N   
                        REAL(KIND=8), INTENT(IN) :: x(N)

                        REAL(KIND=8), INTENT(IN) :: V(N)   ! potential V(x)
                        
                        CHARACTER(LEN=1), INTENT(IN), OPTIONAL :: eigh_range  ! 'A' to compute all eigenvalues, 
                                                                              ! 'I' to compute from the IL-th to the
                                                                              ! IU-th eigenvalue and 'V' to compute all
                                                                              ! the eigenvalues in the half-open 
                                                                              ! interval (VL, VU]. Default is 'A'.

                        INTEGER, INTENT(IN), OPTIONAL :: IL, IU
                        REAL(KIND=8), INTENT(IN), OPTIONAL :: VL, VU

                        LOGICAL, INTENT(IN), OPTIONAL :: eigval_only, store

                        ! ROUTINE 
                        REAL(KIND=8) :: t_start, t_end

                        INTEGER :: d                              ! number of plane waves considered
                        COMPLEX(KIND=8) :: V_fourier(N/2 + 1)     ! Fourier transform of potential V
                        REAL(KIND=8) :: L                         ! interval length
                        REAL(KIND=8), ALLOCATABLE :: k(:)         ! wave vectors

                        COMPLEX(KIND=8), ALLOCATABLE :: H(:,:)

                        CHARACTER(LEN=1) :: JOBZ = 'V', RANGE = 'A', UPLO = 'U'   ! DSTEVR routine parameters
                        INTEGER :: INFO, M, LWORK                                 ! 
                        INTEGER, ALLOCATABLE :: IWORK(:), IFAIL(:)                !
                        REAL(KIND=8) :: ABSTOL, DLAMCH                            !
                        REAL(KIND=8), ALLOCATABLE :: W(:), RWORK(:)               !
                        COMPLEX(KIND=8), ALLOCATABLE :: Z(:,:), WORK(:)           !

                        CHARACTER(LEN=15) :: fmt
      
                        OPEN(UNIT=10, FILE='Input/data_sheet.dat', STATUS='Old', ACTION='Read')
                          DO i = 1, 15
                            READ(10, *)
                          ENDDO

                          READ(10, *) d
                        CLOSE(10)
                        
                        V_Fourier = FFT(x, V, N)   ! Fourier transform of V(x)
                        
                        ALLOCATE(k(d))
                        L = DABS(x(N) - x(1))
                        k(:) = (/ ((2d0 * i - d * 1d0 - 1d0) * pi / L, i = 1,d) /) ! wave vectors

                        ALLOCATE(H(d,d))
                        H = Hamiltonian(V_fourier, d, k)   ! Hamiltonian

                        ! EIGENVALUES CALCULUS

                        IF (PRESENT(eigval_only)) THEN
                                IF (eigval_only .eqv. .TRUE.) JOBZ = 'N'
                        ENDIF
                        IF (PRESENT(eigh_range)) RANGE = eigh_range

                        LWORK = 2*d
                        ALLOCATE(W(d), Z(d,d))
                        ALLOCATE(IWORK(5*d), IFAIL(d), RWORK(7*d), WORK(LWORK))

                        ABSTOL = 2*DLAMCH('S')
                        
                        CALL CPU_TIME(t_start)
                            CALL ZHEEVX(JOBZ, RANGE, UPLO, d, H, d, VL, VU, IL, IU, ABSTOL, M, W, &
                                    Z, d, WORK, LWORK, RWORK, IWORK, IFAIL, INFO)
                        CALL CPU_TIME(t_end)

                        IF (INFO==0) THEN
                                PRINT '(A, I0, A, F15.7, A)', 'Found ', M, ' eigenvalues in', t_end - t_start, ' seconds'

                                IF (PRESENT(store)) THEN
                                        IF (store .eqv. .TRUE.) THEN
                                                PRINT *, 'Writing eigenvalues to "eigenvalues_k.txt"...'

                                                OPEN(UNIT=10, FILE='eigenvalues_k.txt', ACTION='Write')
                                                  WRITE(10, '(I5, F15.7)') (i, W(i), i=1,M)
                                                CLOSE(10)
                                        
                                                PRINT *, ''
                                        
                                                IF (JOBZ == 'V') THEN
                                                        PRINT *, 'Writing eigenvectors to "eigenvectors_k.txt"...'
                                                
                                                        WRITE(fmt, '("(",I5,"F15.7)")') 2*M+1
                                                        OPEN(UNIT=10, FILE='eigenvectors_k.txt', ACTION='Write')
                                                          WRITE(10, fmt) (k(i), Z(i,:M), i=1,d)
                                                        CLOSE(10)

                                                        PRINT *, ''
                                                ENDIF
                                        ELSE
                                                PRINT '(I5, F15.7)', (i, W(i), i=1, MIN(10,M))
                                                PRINT *, ''
                                        ENDIF

                                ELSE
                                        PRINT '(I5, F15.7)', (i, W(i), i=1, MIN(10,M))
                                        PRINT *, ''
                                ENDIF


                        ELSE
                                PRINT *, 'Something went wrong durig the computation of the eigenvalues'
                        ENDIF
                END SUBROUTINE SOLVE_EIGH
                        
END MODULE


