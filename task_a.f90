MODULE TASK_A
        USE read_var
        CONTAINS 
                SUBROUTINE HAMILTONIAN(x, V, dg, e)
                        ! Builds a real tridiagonal matrix starting from the discrete
                        ! potential V(x). Here dg is the diagonal array while 
                        ! e is the sub-diagonal array. 

                        IMPLICIT NONE
                        
                        ! INPUT
                        REAL(KIND=8), INTENT(IN) :: x(:), V(:)  

                        ! OUTPUT
                        REAL(KIND=8) ::  dg(:), e(:)   ! diagonal and sub-diagonal array
                        
                        ! ROUTINE
                        REAL(KIND=8) :: h   ! sampling width
                        
                        h = DABS(x(2) - x(1))

                        dg(:) = 1d0/h**2 + V(:)
                        e(:) = -1d0/(2d0 * h**2)
                END SUBROUTINE HAMILTONIAN

                SUBROUTINE SOLVE_EIGH_X(x, V, N, err_eva, index)
                        ! Solves the eigenvalue problem in the direct space using the LAPACK
                        ! routine DSTEVR (http://www.netlib.org/lapack/).

                        IMPLICIT NONE

                        ! INPUT
                        REAL(KIND=8), INTENT(IN):: x(:), V(:)
                        INTEGER, INTENT(IN) :: N 
                        INTEGER, INTENT(IN), OPTIONAL :: index

                        ! ROUTINE
                        REAL(KIND=8) :: t_start, t_end
                        REAL(KIND=8) :: dg(N), e(N-1)   ! diagonal and sub-diagonal array
                        INTEGER :: M, INFO   
                        INTEGER :: ISUPPZ(2*N), IWORK(10*N)  ! ISUPPZ(2*N), IWORK(10*N),WORK(20*N) 
                        REAL(KIND=8) :: ABSTOL, DLAMCH, WORK(20*N), W(N)       
                        REAL(KIND=8), ALLOCATABLE :: Z(:,:)  ! Eigenvalues and eigenvectors arrays
                        CHARACTER(LEN=15) :: fmt
                        REAL(KIND=8) :: NORM

                        ! OUTPUT
                        REAL(KIND=8), INTENT(OUT), OPTIONAL :: err_eva
  
                        CALL HAMILTONIAN(x, V, dg, e)   ! builds the tridiagonal Hamiltonian 
                        
                        ! EIGENVALUES CALCULUS
                        ABSTOL = 2d0 * DLAMCH('S')

                        IF (PRESENT(index)) THEN  
                                ALLOCATE(Z(N,1))                        
                                CALL DSTEVR(JOBZ, RANGE, N, dg, e, VL, VU, index, index, ABSTOL, M, &
                                            W, Z, N, ISUPPZ, WORK, 20*N, IWORK, 10*N, INFO)
      
                                err_eva = W(1)
                                RETURN 
                        ELSE 
                                ALLOCATE(Z(N,IU-IL+1))   
                                CALL CPU_TIME(t_start)
                                CALL DSTEVR(JOBZ, RANGE, N, dg, e, VL, VU, IL, IU, ABSTOL, M, &
                                W, Z, N, ISUPPZ, WORK, 20*N, IWORK, 10*N, INFO)
                                CALL CPU_TIME(t_end)
                        END IF
                        
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
                                                          WRITE(10, fmt) (x(i), Z(i,1:M), i=1,N)
                                                        CLOSE(10)

                                                        DO j = 1,3
                                                                NORM = 0.d0  
                                                                DO i = 1,N
                                                                NORM = NORM + ABS(Z(i,j)*Z(i,j))
                                                                END DO 
                                                                WRITE(*,'(I2,A24,F20.10)') j,"-th &
                                                                                & eigenfunction norm:", NORM              
                                                        END DO 

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
