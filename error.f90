MODULE error 
    !
    ! performs convergence tests to evaluate computational errors
    !
    USE task_a
    USE task_b
    USE potential
    USE read_var

        IMPLICIT NONE 

        INTEGER, PARAMETER, PRIVATE:: x_step = 1000, k_step = 1

        CONTAINS

        SUBROUTINE ERROR_ANALYSIS(toll, space, percent, max_index)
            IMPLICIT NONE
            
            !! INPUT
            REAL(KIND=8), INTENT(IN) :: toll
            CHARACTER(len=1), INTENT(IN) :: space 
            INTEGER, INTENT(IN) :: percent, max_index

            !!ROUTINE
            INTEGER :: ios, i
            INTEGER :: N, N_default
            REAL(KIND=8) :: L
            INTEGER, ALLOCATABLE :: temp(:,:)

            PRINT*, "DATA ARE PRINTED FOR DEBUGGING..."

            IF (space== "X") THEN 
                PRINT*, "Computing the points of convergence as a function of eigenvalue index in real space."
                PRINT*, "L is so that V(L)=-",percent,"%De"
                PRINT*, "Writing data in ""eva_variation_Ltxt""."

                OPEN(UNIT=15,file="Output/eva_variation_L.txt")  ! evaluation with "analitical" L 
                DO i = 1, max_index
                    CALL WIDTH_CONVERGENCE(toll,percent,i,N,L,"P","X")
                    IF (N == -1) THEN
                        PRINT*, "Maximum number of iterations reached! Try to increse the tolerance, De or decrease the maximum&
                                 &eigenvalue index "
                        EXIT 
                    END IF 
                    CALL WIDTH_CONVERGENCE(toll,percent,i,N_default,L,"D","X")
                    IF (N == -1) THEN
                        PRINT*, "Maximum number of iterations reached! Try to lift up the tolerance..."
                        EXIT 
                    END IF 
                    WRITE(15,'(3I7)') i, N, N_default
                    WRITE(*,'(3I7)') i, N, N_default
                END DO 
                CLOSE(15)
                PRINT*, "computing the points to converge as a function of L in real space."
        
                OPEN(unit=15,FILE="Output/L_variation_perc.txt")
                    DO i = 1,20
                        CALL WIDTH_CONVERGENCE(toll,i,max_index, N,L,"P","X")
                        IF (N == -1) THEN
                            PRINT*, "Maximum number of iterations reached! Try to lift up the tolerance..."
                            EXIT 
                        END IF 
                        WRITE(15,*) i, L, N    
                        print*, i, L, N
                    END DO
                CLOSE(15)
                RETURN
              
            ELSEIF (space == "K") THEN 
                    PRINT*, "Computing the points of convergence as a function of eigenvalue index in reciprocal space."
                    PRINT*, "L is so that V(L)=-",percent,"%De"

                    OPEN(UNIT=15,file="Output/eva_variation_L_K.txt") !evaluation with "analitical" L
                    DO i = 1, max_index
                        CALL WIDTH_CONVERGENCE(toll,percent,i,N,L,"P","K")
                        IF (N==-1) THEN 
                            PRINT*, "Convergence test has been stopped."
                            PRINT*, "Eigenvalue of index",i,"cannot converge!" 
                            PRINT*, "Consider reducing the tolerance or increasing De."
                            EXIT 
                        END IF 
                        CALL WIDTH_CONVERGENCE(toll,percent,i,N_default,L,"D","K")
                        WRITE(15,'(3I7)') i, N, N_default
                        IF (N==-1) THEN
                            PRINT*, "Convergence test has been stopped."
                            PRINT*, "Eigenvalue of index",i,"can not converge!" 
                            PRINT*, "Consider reducing the tolerance or increasing De."
                            EXIT 
                        END IF 
                        print*, i, N, N_default
                    END DO 
                    CLOSE(15)

                    PRINT*, "computing the points to converge as a function of L in reciprocal space."
                    OPEN(unit=15,FILE="Output/L_variation_perc_K.txt")

                    DO i = 1,20
                        CALL WIDTH_CONVERGENCE(toll,i,max_index, N,L,"P","K")
                        WRITE(15,*) i, L, N
                        IF (N==-1) THEN 
                            PRINT*, "Convergence test has been stopped."
                            PRINT*, "L is so that V(L)=-",percent,"%De can not converge!" 
                            PRINT*, "Consider reducing the tolerance, the eigenvalue index or increasing De."
                            EXIT 
                        END IF 
                        print*, i, L, N
                    END DO 
                    CLOSE(15)
                    RETURN

            ELSEIF(space == "H") THEN 
                ALLOCATE(temp(10,20))
                PRINT*, "Computing the points to converge as a function of N and L in real space."
                PRINT*, "Please, wait. it may takes a while..."

                DO i = 1, 10
                    DO j = 1, 20 
                        CALL WIDTH_CONVERGENCE(toll,j,i,N,L,"P","X")
                        temp(i,j) = N
                    END DO 
                    print*, i 
                END DO

                OPEN(UNIT=15,file="Output/heat_map.txt", iostat = ios)
                PRINT*, "Writing the converging points in ""heat_map.txt""..."         
                DO i=1,10
                    WRITE(15,*) temp(i,:)
                END DO
                CLOSE(15)
                RETURN

            ELSE
                PRINT*, "Error raised in subroutine real_variation, in module ""error""."
                PRINT*, """mansion"" received an illegal input value!"
                PRINT*, "mansion can be either ""L"", ""E"" or ""H""."
                RETURN
            END IF  
        END SUBROUTINE ERROR_ANALYSIS
        

        SUBROUTINE WIDTH_CONVERGENCE(toll, percent, index, N, L, mansion, space)
            IMPLICIT NONE 

            !! INPUT
            REAL(KIND=8), INTENT(IN) :: toll
            INTEGER, INTENT(IN) :: percent, index
            CHARACTER(len=1), INTENT(IN) :: mansion, space ! mansion can be either P or D(percent or default)
           
            !! OUTPUT 
            INTEGER, INTENT(OUT) :: N
            REAL(KIND=8), INTENT(OUT) :: L

            IF (mansion == "P") THEN 
            L = x0 - (1/alpha)*DLOG(1.d0-sqrt(1-percent/100.d0))
            ELSEIF (mansion == "D") THEN 
            L = b 
            END IF 

            IF (space == "X") THEN 
                CALL N_d_CONV(toll, L, index, N, "X")
            ELSE IF (space =="K") THEN 
                CALL N_d_CONV(toll, L, index, N, "K")
            ELSE 
                PRINT*, "Error raised in L_conv"
                PRINT*, """space"" received an illegal input value!"
                PRINT*, "mansion can be either X or K"
                RETURN
            END IF

            RETURN 
        END SUBROUTINE WIDTH_CONVERGENCE


        SUBROUTINE N_d_CONV(toll, L, index, N, mansion)
            IMPLICIT NONE 

            !! INPUT
            REAL(KIND=8), INTENT(IN) :: toll
            INTEGER, INTENT(IN) :: index
            REAL(KIND=8), INTENT(IN) :: L
            CHARACTER (LEN=1) :: mansion

            !! ROUTINE
            REAL(KIND=8) :: eva_1, eva_2
            REAL(KIND=8) :: abs_diff
            INTEGER :: maxit, counter,  N_points

            !! OUTPUT
            INTEGER, INTENT(OUT) :: N

            maxit = 200

            IF (MANSION == "X") THEN              
                N_points = x_step
                
                IF (N_points<=INDEX) THEN 
                    PRINT*, "error raised in N_d_conv, task b. x_step was to small."
                    PRINT*, "Provide a number of x_step greater than the eigenvalue index:", index
                    RETURN 
                END IF 
                eva_2 = SINGLE_EVA(N_points, L, index,"X")                
                abs_diff = toll
                counter = 1
                DO WHILE (toll<=abs_diff)
                    counter = counter + 1

                    IF (counter>=maxit) THEN
                        N = -1 
                        RETURN 
                    END IF 

                    N_points = x_step*counter
                    eva_1 = eva_2
                    eva_2 = SINGLE_EVA(N_points, L, index, "X")
                    abs_diff = abs(eva_1-eva_2)
                END DO
                N = N_points
                RETURN

            ELSEIF (MANSION == "K") THEN 
                N_points = 1000
                d = 50

                eva_2 = SINGLE_EVA(N_points, L, index, "K", d)
                abs_diff = toll
                counter = 1
                DO WHILE (toll<=abs_diff)
                    counter = counter + 1
                    
                    IF (counter>=maxit) THEN
                        PRINT*, "Maximum number of iterations reached! Try to lift up the tolerance..."
                        RETURN 
                    END IF 
                    
                    eva_1 = eva_2
                    d = k_step*counter

                    IF (d>=N_points/2) THEN 
                        N = -1 
                        RETURN 
                    END IF 
                    
                    IF (d<=index) THEN 
                        CYCLE 
                    END IF  
        
                    eva_2 = SINGLE_EVA(N_points, L, index, "K",d)
                    abs_diff = abs(eva_1-eva_2)                    
                END DO
                N = d
            END IF 
            RETURN 
        END SUBROUTINE N_d_CONV


        FUNCTION SINGLE_EVA (N, L, index, space, d) RESULT(w_1)
            IMPLICIT NONE 

            !! INPUT 
            INTEGER :: N, index
            REAL(KIND=8) :: L
            INTEGER, OPTIONAL :: d
            CHARACTER(LEN=1) :: space

            !! ROUTINE 
            INTEGER :: i 
            REAL(KIND=8), ALLOCATABLE :: x(:), V(:)
            
            !! OUTPUT 
            REAL(KIND=8) :: w_1
            
            ALLOCATE(x(N), V(N))
            
            x(:) = (/ (L* (i*1d0) / (N-1), i = 0, N-1) /)
            V = MORSE(x,N)
           
            IF (space == "X") THEN
                CALL SOLVE_EIGH_X(x, V,N,w_1,index)       
            ELSEIF (space == "K") THEN 
                CALL SOLVE_EIGH_K(x, V, N, d, w_1,index)
            END IF 

            DEALLOCATE(x,V)
            RETURN 
        END FUNCTION SINGLE_EVA

END MODULE


    


