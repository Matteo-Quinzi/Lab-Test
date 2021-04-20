MODULE error 
    !
    ! performs convergence tests to evaluate computational errors
    !
    USE task_a
    USE potential
    

        IMPLICIT NONE 
    !
        CONTAINS


        SUBROUTINE L_CONV(toll, alpha, x0, percent, index, Nconv, L)
            IMPLICIT NONE 

            !! INPUT
            REAL(KIND=8), INTENT(IN) :: x0, alpha, toll
            INTEGER, INTENT(IN) :: percent, index
           
            !! OUTPUT 
            INTEGER, INTENT(OUT) :: Nconv
            REAL(KIND=8), INTENT(OUT) :: L
            !
            L = x0 - (1/alpha)*DLOG(1.d0-sqrt(1-percent/100.d0))
            !
           ! L = 20.d0
            CALL POINTS_CONVERGENCE(toll, 0.d0, L, index, Nconv)
            RETURN 
        END SUBROUTINE L_CONV

            

        SUBROUTINE POINTS_CONVERGENCE(toll, a, b, index,Nend)
            IMPLICIT NONE 

            !! INPUT
            REAL(KIND=8), INTENT(IN) :: toll
            INTEGER, INTENT(IN) :: index
            REAL(KIND=8), INTENT(IN) :: a,b

            !! ROUTINE
            REAL(KIND=8) :: eva_1, eva_2
            INTEGER :: Npoints  
            REAL(KIND=8) :: abs_diff
            INTEGER :: maxit, count = 1

            !! OUTPUT
            INTEGER, INTENT(OUT) :: Nend

            maxit = 1000    !safe exit condition: maximum number of iterations 

            Npoints = 1000  ! minimum number of points to be considered
            !
            eva_2 = SINGLE_EVA(Npoints, a, b, index)
            !        
            abs_diff = toll
            DO WHILE (toll<=abs_diff)
                !
                count = count + 1
                IF (count>=maxit) THEN
                    PRINT*, "Maximum number of iterations reached! Try to lift up the tolerance..."
                    RETURN 
                END IF 
                !
                eva_1 = eva_2
                !
                Npoints = 1000*count
                eva_2 = SINGLE_EVA(Npoints, a, b, index)
                !print*, eva_2, Npoints, count
        
                abs_diff = abs(eva_1-eva_2)

            END DO

            Nend = Npoints
            
            RETURN 
        END SUBROUTINE POINTS_CONVERGENCE



        FUNCTION SINGLE_EVA (Npoints, a, b, index)
            IMPLICIT NONE 

            !! INPUT 
            INTEGER :: Npoints
            REAL(KIND=8), INTENT(IN) :: a , b
            INTEGER :: index , ierr

            !! ROUTINE 
            INTEGER :: i 
            REAL(KIND=8), ALLOCATABLE :: x(:), V(:)

            !! OUTPUT 
            REAL(KIND=8) :: SINGLE_EVA

            ALLOCATE(x(Npoints), V(Npoints),stat=ierr)
            if (ierr /= 0) print*, 'Memory error 1!'
            !
            x(:) = (/ (a + (b-a) * (i*1d0) / (Npoints-1), i = 0, Npoints-1) /)
            V = MORSE(x)
            !
            CALL SOLVE_EIGH_X(x, V, SINGLE_EVA,index)            
            !
            DEALLOCATE(x,V)
            
            RETURN 
        END FUNCTION SINGLE_EVA


END MODULE


    


