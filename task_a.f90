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
                        REAL(KIND=8), ALLOCATABLE :: dg(:), e(:)   ! diagonal and sub-diagonal array
                        
                        ! ROUTINE
                        REAL(KIND=8) :: h   ! sampling width
                        INTEGER :: Npoints
                        !
                        Npoints = size(x)        
                        !
                        h = DABS(x(2) - x(1))

                        ALLOCATE(dg(Npoints),e(Npoints-1))

                        dg(:) = 1d0/h**2 + V(:)
                        e(:) = -1d0/(2d0 * h**2)
                END SUBROUTINE HAMILTONIAN

                SUBROUTINE SOLVE_EIGH_X(x, V, err_eva, index)
                        ! Solves the eigenvalue problem in the direct space using the LAPACK
                        ! routine DSTEVR (http://www.netlib.org/lapack/).

                        IMPLICIT NONE

                        ! INPUT
                        REAL(KIND=8), INTENT(IN):: x(:), V(:)
                        INTEGER, INTENT(IN), OPTIONAL:: index

                        ! ROUTINE
                        REAL(KIND=8) :: t_start, t_end

                        REAL(KIND=8), ALLOCATABLE :: dg(:), e(:)   ! diagonal and sub-diagonal array

                        INTEGER :: M, INFO   
                        INTEGER, ALLOCATABLE :: ISUPPZ(:), IWORK(:)  ! ISUPPZ(2*N), IWORK(10*N),WORK(20*N) 
                        REAL(KIND=8) :: ABSTOL, DLAMCH   
                        real(kind=8), ALLOCATABLE :: WORK(:)  ! WORK(20*N)        

                        REAL(KIND=8), ALLOCATABLE :: W(:), Z(:,:)  ! Eigenvalues and eigenvectors arrays

                        CHARACTER(LEN=15) :: fmt
                        REAL(KIND=8) :: dx
                        INTEGER :: Npoints, ierr

                        ! OUTPUT
                        REAL(KIND=8), INTENT(OUT), OPTIONAL :: err_eva

                        !
                        Npoints = size(x)
                        !
                        
                        CALL HAMILTONIAN(x, V, dg, e)   ! builds the tridiagonal Hamiltonian 
                        
                        ! EIGENVALUES CALCULUS
                        ABSTOL = 2d0 * DLAMCH('S')

                        IF (PRESENT(index)) THEN  
                                IL = index
                                IU = index 
                        END IF 

                        ALLOCATE(ISUPPZ(2*Npoints),IWORK(10*Npoints),WORK(20*Npoints),stat=ierr)
                        if (ierr /= 0) print*, 'Memory error 2!'
                        ALLOCATE(W(Npoints),Z(Npoints,IU-IL+1),stat=ierr)
                        if (ierr /= 0) print*, 'Memory error 3!'
                        
                        CALL CPU_TIME(t_start)
                            CALL DSTEVR(JOBZ, RANGE, Npoints, dg, e, VL, VU, IL, IU, ABSTOL, M, &
                                    W, Z, Npoints, ISUPPZ, WORK, 20*Npoints, IWORK, 10*Npoints, INFO)
                        CALL CPU_TIME(t_end)

                        IF (PRESENT(err_eva))THEN
                                err_eva = W(1)
                                RETURN 
                        END IF 

                        DEALLOCATE(ISUPPZ,IWORK,WORK,dg,e)

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
                                                          WRITE(10, fmt) (x(i), Z(i,1:M), i=1,Npoints)
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
