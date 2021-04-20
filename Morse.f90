PROGRAM MAIN

        !! MODULE

        USE read_var
        USE potential   ! defines the potential function
        
        USE task_a     ! the computation is performed in the coordinates space
        USE task_b    ! the computation is performed in the reciprocal space
        USE task_c
        USE error

        
        !! VARIABLE DECLARATIONS       

        IMPLICIT NONE 
 
        REAL(KIND=8), ALLOCATABLE :: x(:)   ! Points of the interval
        REAL(KIND=8), ALLOCATABLE :: V(:)   ! Potential V(x)
        
        !error analysis parameters
        INTEGER :: Nend, index
        REAL(KIND=8) :: toll, L_end

        !! CODE
        ! DATA READING
        CALL READ_DATA()
        WRITE(*,*) 'You asked for task -->  ', task
        
        ALLOCATE(x(N), V(N))
        x(:) = (/ (a + (b-a) * (i*1d0) / (N-1), i = 0, N-1) /)   ! discrete interval [a, b]
       
        IF (V_str == 'M') V = MORSE(x)              ! Morse potential V(x)
        IF (V_str == 'H') V = HARMONIC(x)           ! harmonic potential V(x)
        IF (V_str == 'X') V = SHIFTED_HARMONIC(x)   ! 
        
        OPEN(UNIT=10, FILE='Output/pot.txt')
          WRITE(10, '(F15.7,F20.10)') (x(i), V(i), i=1,N)
        CLOSE(10)

        IF (task == 'a') CALL SOLVE_EIGH_X(x, V)                                  
        IF (task == 'b') CALL SOLVE_EIGH_K(x, V)
        !IF (task == 'c') CALL SOLVE_EIGH_HO(x, V)

        !!!!!


        
        toll = 1.d-7
        index = 20
        CALL L_CONV(toll, alpha, x0, 5, index, Nend, L_end) !toll, alpha, x0, percent, index, Nconv, L
        print*, Nend, L_end


                                                
END PROGRAM
