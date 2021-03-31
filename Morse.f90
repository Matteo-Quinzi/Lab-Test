PROGRAM MAIN

        !! MODULE

        USE potential   ! defines the potential function
        
        !USE task_a     ! the computation is performed in the coordinates space
        USE task_b    ! the computation is performed in the reciprocal space
        
        !! VARIABLE DECLARATIONS       

        IMPLICIT NONE 
 
        INTEGER :: N   
        REAL(KIND=8) :: a, b   ! Interval [a,b] made up of N points

        REAL(KIND=8), ALLOCATABLE :: x(:)   ! Points of the interval
        REAL(KIND=8), ALLOCATABLE :: V(:)   ! Potential V(x)

        !! CODE

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

        OPEN(UNIT=10, FILE='pot.txt')
          WRITE(10, '(2F15.7)') (x(i), V(i), i=1,N)
        CLOSE(10)

        CALL SOLVE_EIGH(x, N, V, eigval_only = .FALSE., store = .TRUE., &   ! store is an optional logical variable: .TRUE. to 
                eigh_range='I', IL = 1, IU = 10)                            ! write the eigenvalues and the eigenvectors to a 
                                                                            ! text file, .FALSE. to print out the first eigenvalues.
                                                                            ! Default is .FALSE..
                                                                            !
                                                                            ! eigh_range (optional) controls the range of the 
                                                                            ! eigenvalues to compute. It's 'A' (all eigenvalues), 
                                                                            !'I' (from the IL-th to the IU-th) or 'V' (all 
                                                                            ! eigenvalues in the half-open interval (VL,VU]).
                                                                            ! Default is 'A'. 
END PROGRAM
