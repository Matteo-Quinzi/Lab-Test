MODULE TASK_A
        USE read_var
        CONTAINS 
                SUBROUTINE HAMILTONIAN(x, V, dg, e)
                        ! Builds a real tridiagonal matrix starting from the discrete
                        ! potential V(x). Here dg is the diagonal array while 
                        ! e is the sub-diagonal array. 

                        IMPLICIT NONE
                        
                        ! INPUT
                        REAL(KIND=8), INTENT(IN) :: x(N), V(N)  

                        ! OUTPUT
                        REAL(KIND=8) :: dg(N), e(N-1)   ! diagonal and sub-diagonal array
                        
                        ! ROUTINE
                        REAL(KIND=8) :: h   ! sampling width
                        
                        h = DABS(x(2) - x(1))

                        dg(:) = 1d0/h**2 + V(:)
                        e(:) = -1d0/(2d0 * h**2)
                END SUBROUTINE HAMILTONIAN

                SUBROUTINE SOLVE_EIGH_X(x, V)
                        ! Solves the eigenvalue problem in the reciprocal space using the LAPACK
                        ! routine DSTEVR (http://www.netlib.org/lapack/).

                        IMPLICIT NONE

                        ! INPUT
                        REAL(KIND=8), DIMENSION(N) :: x, V

                        ! ROUTINE
                        REAL(KIND=8) :: t_start, t_end

                        REAL(KIND=8) :: dg(N), e(N-1)   ! diagonal and sub-diagonal array

                        INTEGER :: M, ISUPPZ(2*N), IWORK(10*N), INFO   !
                        REAL(KIND=8) :: ABSTOL, DLAMCH, WORK(20*N)     !
                        REAL(KIND=8) :: W(N), Z(N,N)                   ! Eigenvalues and eigenvectors array

                        CHARACTER(LEN=15) :: fmt
                        REAL(KIND=8) :: dx
                        
                        CALL HAMILTONIAN(x, V, dg, e)   ! builds the tridiagonal Hamiltonian 
                        
                        ! EIGENVALUES CALCULUS
                        ABSTOL = 2d0 * DLAMCH('S')
                        
                        CALL CPU_TIME(t_start)
                            CALL DSTEVR(JOBZ, RANGE, N, dg, e, VL, VU, IL, IU, ABSTOL, M, &
                                    W, Z, N, ISUPPZ, WORK, 20*N, IWORK, 10*N, INFO)
                        CALL CPU_TIME(t_end)

                        dx = x(2) - x(1)
                        Z(:,:) = Z(:,:) * dx**(-0.5)
                        
                        IF (INFO==0) THEN
                                PRINT '(A, I0, A, F15.7, A)', 'Found ', M, ' eigenvalues in', t_end - t_start, ' seconds'

                                        IF (store .eqv. .TRUE.) THEN
                                                PRINT *, 'Writing eigenvalues to "eigenvalues.txt"...'

                                                OPEN(UNIT=10, FILE='Output/eigenvalues.txt', ACTION='Write')
                                                  WRITE(10, '(I5, F15.7)') (i, W(i), i=1,M)
                                                CLOSE(10)
                                        
                                                PRINT *, ''
                                        
                                                IF (JOBZ == 'V') THEN
                                                        PRINT *, 'Writing eigenvectors to "eigenvectors.txt"...'
                                                
                                                        WRITE(fmt, '("(",I5,"F15.7)")') M+1
                                                        OPEN(UNIT=10, FILE='Output/eigenvectors.txt', ACTION='Write')
                                                          WRITE(10, fmt) (x(i), Z(i,:M), i=1,N)
                                                        CLOSE(10)

                                                        PRINT *, ''
                                                ENDIF
                                        ELSE
                                                PRINT '(I5, F15.7)', (i, W(i), i=1, MIN(10,M))
                                                PRINT *, ''
                                        ENDIF

                        ELSE
                                PRINT *, 'Something went wrong durig the computation of the eigenvalues'
                        ENDIF
                END SUBROUTINE SOLVE_EIGH_X
                !
                !
                !
        END MODULE 
