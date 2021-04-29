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
        REAL(KIND=8), ALLOCATABLE  :: V(:)   ! Potential V(x)
        
        !error analysis parameters
        REAL(KIND=8) :: toll
  
        !! CODE
        ! DATA READING
        CALL READ_DATA()

        WRITE(*,*) 'You asked for task -->  ', task
        
        ALLOCATE(x(N), V(N))
        x(:) = (/ (a + (b-a) * (i*1d0) / (N-1), i = 0, N-1) /)   ! discrete interval [a, b]
       
        IF (V_str == 'M') V = MORSE(x,N)              ! Morse potential V(x)
        IF (V_str == 'H') V = HARMONIC(x,N)           ! harmonic potential V(x)
        IF (V_str == 'X') V = SHIFTED_HARMONIC(x,N)   ! 
        
        OPEN(UNIT=10, FILE='Output/pot.txt')
          WRITE(10, '(F15.7,F20.10)') (x(i), V(i), i=1,N)
        CLOSE(10)

        IF (task == 'a') CALL SOLVE_EIGH_X(x, V, N)                                  
        IF (task == 'b') CALL SOLVE_EIGH_K(x, V, N, d)
        IF (task == 'c') CALL SOLVE_EIGH_HO(x, V)

        
        toll = 1.d-7
        CALL ERROR_ANALYSIS(toll,"X",10,15)    

                                                
END PROGRAM
