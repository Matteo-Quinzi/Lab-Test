MODULE task_c 
    USE potential
    USE task_a
    IMPLICIT NONE 

    INTEGER, PRIVATE :: i,j 

    CONTAINS 
            !
            !
            !
            FUNCTION QUADRATURE_1D(N, dx, f)
                ! Compute the integral of the function f, sampled 
                ! on intervals of width dx, using the trapezoidal formula.

                ! Input 
                INTEGER, INTENT(IN) :: N 
                REAL(KIND=8), INTENT(IN) :: dx 
                REAL(KIND=8), DIMENSION(N), INTENT(IN) :: f 

                ! Output
                REAL(KIND=8) :: QUADRATURE_1D

                QUADRATURE_1D = 0.5 * (f(1) + f(N)) + SUM(f(2:N-1))
                QUADRATURE_1D = dx * QUADRATURE_1D
            END FUNCTION
            !
            !
            !
            FUNCTION HAMILTONIAN_HO(M, N, dx, eigenvalues, eigenstates, V_eff)
                ! Defines the Hamiltonian in the space spanned by the first 
                ! M eigenstates of the harmonic oscillator.
                ! Has the Hamiltonian will be real and symmetric, only the upper
                ! triangle is computed.

                ! Input 
                INTEGER, INTENT(IN) :: N, M 
                REAL(KIND=8), INTENT(IN) :: dx 
                REAL(KIND=8), INTENT(IN) :: eigenvalues(M), eigenstates(N,M), V_eff(N)
                
                !Output
                REAL(KIND=8), DIMENSION(M,M) :: HAMILTONIAN_HO

                !Internal variables
                REAL(KIND=8), DIMENSION(N) :: temp_func ! temporary function useful during integration
                INTEGER :: i, j, k 

                ! Initializing H with elements on the diagonal 

                DO i = 1, M
                    DO j = i, M
                        ! The matrix element H(i,j) is defined by the relation
                        ! < i | V_eff | j >
                        temp_func = (/ ( eigenstates(k,i) * V_eff(k) * eigenstates(k,j) , k = 1,N ) /)
                        HAMILTONIAN_HO(i,j) = QUADRATURE_1D(N, dx, temp_func)
                        IF (j==i) THEN
                            HAMILTONIAN_HO(i,j) = HAMILTONIAN_HO(i,j) + eigenvalues(i)
                        END IF
                    END DO 
                END DO 

            END FUNCTION HAMILTONIAN_HO
            !
            !
            !
            SUBROUTINE SOLVE_HAMILTONIAN_REAL(N, H, eigval_only, eigh_range, IL, IU, VL, VU, store)
                ! Solve the iegenvalues problem for the real symmetric matrix H.
                ! Lapack DSYEVX routine is used (http://www.netlib.org/lapack)

                ! Input 
                INTEGER, INTENT(INOUT):: N ! Order of matrix H and number of wanted eigenvalues
                REAL(KIND=8), INTENT(INOUT) :: H(N,N)
                LOGICAL, INTENT(IN), OPTIONAL :: eigval_only, store
                CHARACTER(1), INTENT(IN), OPTIONAL:: eigh_range
                INTEGER, INTENT(IN), OPTIONAL :: IL, IU        ! defines which eigenvalues to compute
                REAL(KIND=8), INTENT(IN), OPTIONAL :: VL, VU   ! defines the range in which looking for eigenvalues 

                ! time variables
                REAL(KIND=8) :: t_start, t_end

                !INTERNAL DSYEVZ VARIABLES
                CHARACTER(1) :: JOBZ = 'V', UPLO = 'U', RANGE = 'A'       !
                INTEGER :: INFO, M, LWORK, IL_in=0, IU_in=0               ! 
                INTEGER, ALLOCATABLE :: IWORK(:), IFAIL(:)                !
                REAL(KIND=8) :: ABSTOL, DLAMCH, VL_in=0., VU_in=0.        !
                REAL(KIND=8), ALLOCATABLE :: W(:)                         !
                REAL(KIND=8), ALLOCATABLE :: Z(:,:), WORK(:)              !

                ! Reading optional input parameters
                IF (PRESENT(eigval_only)) THEN
                    IF (eigval_only .eqv. .TRUE.) JOBZ = 'N'
                ENDIF
                IF (PRESENT(eigh_range)) RANGE = eigh_range
                IF (PRESENT(IL)) IL_in = IL
                IF (PRESENT(IU)) IU_in = IU
                IF (PRESENT(VL)) VL_in = VL 
                IF (PRESENT(VU)) VU_in = VU

                ! Allocating memory
                LWORK = 8*N 
                ALLOCATE(W(N), Z(N,N))
                ALLOCATE(IWORK(5*N), IFAIL(N), WORK(LWORK))

                
                !Defining optimal precision
                ABSTOL = 2*DLAMCH('S')

                
                !Calling DSYEVX routine
                CALL CPU_TIME(t_start)
                            CALL DSYEVX(JOBZ, RANGE, UPLO, N, H, N, VL_in, VU_in, IL_in, IU_in, ABSTOL, M, W, &
                                    Z, N+1, WORK, LWORK, IWORK, IFAIL, INFO)
                CALL CPU_TIME(t_end)
                
                IF (INFO==0) THEN
                    PRINT '(A, I0, A, F15.7, A)', 'Found ', M, ' eigenvalues in', t_end - t_start, ' seconds'
                ELSE
                    PRINT *, 'Something went wrong durig the computation of the eigenvalues'
                ENDIF
                WRITE(*,*)
                !
                IF (PRESENT(store)) THEN
                    IF (store .eqv. .TRUE.) THEN
                            PRINT *, 'Writing eigenvalues to "Output/eigenvalues_ho.txt"...'

                            OPEN(UNIT=10, FILE='Output/eigenvalues_ho.txt', ACTION='Write')
                              WRITE(10, '(I5, F15.7)') (i, W(i), i=1,M)
                            CLOSE(10)
                    
                            PRINT *, ''
                    
                            IF (JOBZ == 'V') THEN
                                    PRINT *, 'Writing eigenvectors to "Output/eigenvectors_ho.txt"...'
                            
                                    OPEN(UNIT=10, FILE='Output/eigenvectors_ho.txt', ACTION='Write')
                                      DO i=1,M
                                          WRITE(10, *) Z(i,:)
                                      END DO 
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
                STOP 'Bye'
            END SUBROUTINE SOLVE_HAMILTONIAN_REAL
            !
            !
            !
            SUBROUTINE SOLVE_EIGH_HO(x, N, V, eigval_only, eigh_range, IL, IU, VL, VU, store)
                ! Solve eigenvalues problem in the space
                ! spanned by the harmonic oscillator eigenstates

                !Input 
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

                !Internal Variables
                INTEGER :: i, M = 10, j
                INTEGER:: io_unit = 12
                REAL(KIND=8) :: V_ho(N), V_eff(N), dx, temp
                REAL(KIND=8), ALLOCATABLE :: eva_ho(:), eve_ho(:,:), H_ho(:,:) 

                IF (PRESENT(IU))  M = IU 
                ALLOCATE(eva_ho(M), eve_ho(N,M), H_ho(M,M))
                

                ! Approximating Morse potential with harmonic potential
                V_ho =  SHIFTED_HARMONIC(x, N)
                OPEN(unit = io_unit, file='Output/pot_ho.txt')
                    WRITE(io_unit, fmt='(2F15.7)') ( x(i), V_ho(i), i=1,N)
                CLOSE(io_unit)

                ! Solving eigenvalues problem for the harmonic potential
                CALL SOLVE_EIGH_X(x, N, V_ho, eigval_only = .FALSE., store = .TRUE., &
                eigh_range = 'I', IL = 1, IU = IU)

                ! Reading eigenvalues and eigenvectors 
                OPEN(unit = io_unit, file='Output/eigenvalues.txt', STATUS='old')
                   DO i=1, M 
                    READ(io_unit,*) j, eva_ho(i)
                   END DO
                CLOSE(io_unit)

                OPEN(unit = io_unit, file='Output/eigenvectors.txt', STATUS='old')
                    DO i = 1,N
                        READ(io_unit, *) temp, eve_ho(i,:) 
                    END DO
                CLOSE(io_unit)

                
                ! Defining effective potential
                V_eff = (/ ( V(i) - V_ho(i) , i = 1,N ) /)

                !Building Hamiltonian in HO space 
                dx = x(2) - x(1)
                H_ho = HAMILTONIAN_HO(M, N, dx, eva_ho, eve_ho, V_eff)


                !Solving eigenvalues problem in HO space 
                ! Computing all eigenvalues
                CALL SOLVE_HAMILTONIAN_REAL(M, H_ho, eigval_only=.FALSE.,eigh_range='A', store=store)

            END SUBROUTINE SOLVE_EIGH_HO 
            !
            !
            !

END MODULE task_c