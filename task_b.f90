MODULE TASK_B
        USE read_var
        CONTAINS
                FUNCTION FFT(x, V, N) 
                        ! Computes the discrete Fourier transform of the function V
                        ! through the FFTW routine (http://www.fftw.org)

                        IMPLICIT NONE

                        INCLUDE 'fftw3.f'

                        ! INPUT
                        REAL(KIND=8), INTENT(IN) :: x(:), V(:)
                        INTEGER, INTENT(IN) :: N
                      
                        ! OUTPUT
                        COMPLEX(KIND=8):: FFT(N/2+1) ! Fourier transfrom of V
                        COMPLEX(KIND=8) :: FFT2(N/2+1)

                        ! ROUTINE PARAMETERS
                        INTEGER(KIND=8) :: plan
                        REAL(KIND=8) :: k(N/2+1), IFT(N)
                        REAL(KIND=8) :: L
                
                        CALL DFFTW_PLAN_DFT_R2C_1D(plan, N, V, FFT, FFTW_ESTIMATE)
                        CALL DFFTW_EXECUTE_DFT_R2C(plan, V, FFT)
                        CALL DFFTW_DESTROY_PLAN(plan)

                        FFT = FFT / DBLE(N)
                        FFT2 = FFT

                        CALL DFFTW_PLAN_DFT_C2R_1D(plan, N, FFT2, IFT, FFTW_ESTIMATE)
                        CALL DFFTW_EXECUTE_DFT_C2R(plan, FFT2, IFT)
                        CALL DFFTW_DESTROY_PLAN(plan)


                        L = DABS(x(N) - x(1))
                        k(:) = (/ ((2d0 * i) * pi / L, i=0,N/2) /)

                        OPEN(UNIT=10, FILE='Output/fourier.txt')
                          WRITE(10, '(6F15.7)') (x(i), V(i), k(i), FFT(i), IFT(i), i=1,N/2+1)
                          WRITE(10, '(6F15.7)') (x(i), V(i), k(N/2+1), FFT(N/2+1), IFT(i), i = N/2+2, N)
                        CLOSE(10)

                END FUNCTION FFT

                FUNCTION HAMILTONIAN(V, k, d) RESULT(H)
                        ! Returns the Hamiltonian (the upper triangle) obtained 
                        ! after expanding the wave function on the plane 
                        ! waves basis set.

                        IMPLICIT NONE

                        ! INPUT
                        REAL(KIND=8), INTENT(IN) :: k(:)   ! wave vectors
                        COMPLEX(KIND=8), INTENT(IN) :: V(:) ! Fourier transform of the potential
                        INTEGER, INTENT(IN) :: d

                        ! OUTPUT
                        COMPLEX(KIND=8) ::  H(d,d)
                        
                        H(:,:) = 0.d0
                        DO i = 1,d      
                          H(i,i) = V(1) + 0.5 * k(i)**2
                          H(i,i+1:) = V(2:d+1-i)
                        END DO

                END FUNCTION HAMILTONIAN
          
                SUBROUTINE SOLVE_EIGH_K(x, V, N, d, err_eva, index)
                        ! Computes the eigenvalues and the eigenvectors of an
                        ! upper triangle hermitian matrix H by calling the LAPACK 
                        ! routine ZHEEVX (http://www.netlib.org/lapack)

                        IMPLICIT NONE

                        INCLUDE 'fftw3.f'

                        ! INPUT
                        REAL(KIND=8), INTENT(IN) :: x(:)
                        REAL(KIND=8), INTENT(IN) :: V(:)   ! potential V(x)
                        INTEGER, INTENT(IN), OPTIONAL::  index
                        INTEGER, INTENT(IN) :: d, N 

                        ! ROUTINE 
                        REAL(KIND=8) :: t_start, t_end
                        COMPLEX(KIND=8) :: V_fourier(N/2+1)     ! Fourier transform of potential V  
                        REAL(KIND=8) :: L                         ! interval length
                        REAL(KIND=8) :: k(d)                      ! wave vectors (d)
                        COMPLEX(KIND=8):: H(d,d)                       !(d,d)
                
                        CHARACTER(LEN=1) :: UPLO = 'U'                    ! DSTEVR routine parameters
                        INTEGER :: INFO, M                     ! 
                        INTEGER :: IWORK(5*d), IFAIL(d)        !   (5*d), (d)
                        REAL(KIND=8) :: ABSTOL, DLAMCH 
                        REAL(KIND=8) :: NORM       !
                        REAL(KIND=8) :: W(d), RWORK(7*d)       !    (d), (7*d)
                        COMPLEX(KIND=8):: Z(d,d), WORK(2*d)   !    (d,d), (2*d)
                        CHARACTER(LEN=15) :: fmt

                        !OUTPUT
                        REAL(KIND=8), INTENT(OUT), OPTIONAL :: err_eva

                        V_Fourier = FFT(x, V, N)   ! Fourier transform of V(x)

                        L = DABS(x(N) - x(1))      

                        ABSTOL = 2*DLAMCH('S')

                        k(:) = (/ ((2d0 * i - d * 1d0 - 1d0) * pi / L, i = 1,d) /)
                        H = Hamiltonian(V_fourier, k, d)   ! Hamiltonian: complex, sparse, hermitian

                        IF (PRESENT(index)) THEN
                                CALL ZHEEVX(JOBZ, RANGE, UPLO, d, H, d, VL, VU, index, index, ABSTOL, M, W, &
                                Z, d, WORK, 2*d, RWORK, IWORK, IFAIL, INFO)
                                err_eva = W(1)    ! During the error analysis only the eigenvalue of intrest is returned
                                RETURN
                        ELSE 
                                CALL CPU_TIME(t_start)
                                CALL ZHEEVX(JOBZ, RANGE, UPLO, d, H, d, VL, VU, IL, IU, ABSTOL, M, W, &
                                Z, d, WORK, 2*d, RWORK, IWORK, IFAIL, INFO)
                                CALL CPU_TIME(t_end)
                        END IF
                        
                        IF (INFO==0) THEN
                                PRINT '(A, I0, A, F15.7, A)', 'Found ', M, ' eigenvalues in', t_end - t_start, ' seconds'
                                        IF (store .eqv. .TRUE.) THEN
                                                PRINT *, 'Writing eigenvalues to "Output/eigenvalues_k.txt"...'

                                                OPEN(UNIT=10, FILE='Output/eigenvalues_k.txt', ACTION='Write')
                                                  WRITE(10, '(I5, F15.7)') (i, W(i), i=1,M)
                                                CLOSE(10)
                                        
                                                PRINT *, ''
                                        
                                                IF (JOBZ == 'V') THEN
                                                        PRINT *, 'Writing eigenvectors to "eigenvectors_k.txt"...'
                                                
                                                        WRITE(fmt, '("(",I5,"F15.7)")') 2*M+1
                                                        OPEN(UNIT=10, FILE='Output/eigenvectors_k.txt', ACTION='Write')
                                                          WRITE(10, fmt) (k(i), Z(i,:M), i=1,d)
                                                        CLOSE(10)

                                                        
                                                        DO j = 1,3
                                                                NORM = 0.d0  
                                                                DO i = 1,d
                                                                NORM = NORM + ABS(Z(i,j)*Z(i,j))
                                                                END DO 
                                                                WRITE(*,'(A15,I2,A40,F20.10)') "Squared sum of the",j,"-th &
                                                                                & eigefunction Fourier's coefficient:", NORM
                                                                                        
                                                        END DO 

                                                        PRINT *, ''
                                                ENDIF
                                        ELSE
                                                PRINT '(I5, F15.7)', (i, W(i), i=1, MIN(10,M))
                                                PRINT *, ''
                                        ENDIF
                        ELSE
                                PRINT *, 'Something went wrong during the computation of the eigenvalues'
                        ENDIF
                END SUBROUTINE SOLVE_EIGH_K
                        
END MODULE

