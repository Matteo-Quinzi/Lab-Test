PROGRAM MAIN

        !! MODULE

        USE potential   ! defines the potential function
        
        USE task_a     ! the computation is performed in the coordinates space
        USE task_b    ! the computation is performed in the reciprocal space
        USE task_c
        
        !! VARIABLE DECLARATIONS       

        IMPLICIT NONE 
 
        INTEGER :: N, i
        CHARACTER(1) :: task
        REAL(KIND=8) :: a, b   ! Interval [a,b] made up of N points

        REAL(KIND=8), ALLOCATABLE :: x(:)   ! Points of the interval
        REAL(KIND=8), ALLOCATABLE :: V(:)   ! Potential V(x)

        !! CODE

        ! Reading Task 
        OPEN(UNIT = 10, FILE='Input/data_sheet.dat', STATUS = 'old', ACTION = 'read')
            DO i = 1,18
              READ(10,*)
            END DO

            READ(10,*) task
            WRITE(*,*) 'You asked for task -->  ', task
        CLOSE(10)

        ! DATA READING
        OPEN(UNIT = 10, FILE = 'Input/data_sheet.dat', STATUS = 'old', ACTION = 'read')
          READ(10, *)
          READ(10, *)      
          READ(10, *) a, b, N
            
        CLOSE(10)

        ALLOCATE(x(N), V(N))
        x(:) = (/ (a + (b-a) * (i*1d0) / (N-1), i = 0, N-1) /)   ! discrete interval [a, b]
       
        !V = HARMONIC(x, N)   ! harmonic potential V(x)
        V = MORSE(x, N)     ! Morse potential V(x)

        OPEN(UNIT=10, FILE='Output/pot.txt')
          WRITE(10, '(2F15.7)') (x(i), V(i), i=1,N)
        CLOSE(10)

        IF (task == 'a') THEN
            CALL SOLVE_EIGH_X(x, N, V, eigval_only = .FALSE., store = .TRUE., &   ! store is an optional logical variable: .TRUE. to 
                      eigh_range='I', IL = 1, IU = 10)                            ! write the eigenvalues and the eigenvectors to a 
                                                                                  ! text file, .FALSE. to print out the first eigenvalues.
        ELSE IF (task == 'b') THEN                                                ! store is an optional logical variable: .TRUE. to
            CALL SOLVE_EIGH_K(x, N, V, eigval_only = .FALSE., store = .TRUE., &   ! Default is .FALSE.. 
                      eigh_range='I', IL = 1, IU = 10)                            ! eigh_range (optional) controls the range of the
                                                                                  !'I' (from the IL-th to the IU-th) or 'V' (all 
        ELSE IF (task == 'c') THEN                                                ! eigenvalues in the half-open interval (VL,VU]).
            CALL SOLVE_EIGH_HO(x, N, V, eigval_only = .FALSE., store = .TRUE., &  ! Default is 'A'.   
            eigh_range='I', IL = 1, IU = 20)
        ELSE 
          STOP 'Invalid task !!!'
        END IF
                                                                                                                                                  
END PROGRAM
