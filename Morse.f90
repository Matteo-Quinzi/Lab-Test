PROGRAM MAIN

! "Morse" evaluates selected eigenvalues and, optionally, eigenfunctions 
!  of the a particle in presence of the Morse potential. 
!  Computation can be performed either in coordinate or momentum representation.
!  Convergences tests can be executed to find the input parameters that miaximize stability.


        USE read_var    ! read from "data_sheet.dat" and allocate the input paramaters 
        USE potential   ! defines the potential function
        USE task_a     ! perform the computation in the coordinates space
        USE task_b    ! perform the computation in the reciprocal space
        USE error

        IMPLICIT NONE 

        REAL(KIND=8), ALLOCATABLE :: x(:)   ! Interval [a,b] equapartition
        REAL(KIND=8), ALLOCATABLE  :: V(:)   ! Potential V(x)

        REAL(KIND=8) :: toll ! Tolerance of the convergence tests
        INTEGER :: maxInd, maxPotPerc  !Error analysis parameters


        ! DATA READING
        CALL READ_DATA()
        WRITE(*,*) 'You asked for task -->  ', task
        
        ALLOCATE(x(N), V(N))
        x(:) = (/ (a + (b-a) * (i*1d0) / (N-1), i = 0, N-1) /)   ! discrete interval [a, b]
        IF (V_str == 'M') V = MORSE(x,N)              ! Morse potential V(x)
        IF (V_str == 'H') V = HARMONIC(x,N)           ! harmonic potential V(x)
        IF (V_str == 'X') V = SHIFTED_HARMONIC(x,N)   ! harmonic potential shifted to overlap with Morse potential
        
        OPEN(UNIT=10, FILE='Output/pot.txt')
          WRITE(10, '(F15.7,F20.10)') (x(i), V(i), i=1,N)
        CLOSE(10)

        IF (task == 'a') CALL SOLVE_EIGH_X(x, V, N)   ! computation is performed in real space   
        IF (task == 'b') CALL SOLVE_EIGH_K(x, V, N, d)  ! computation is perfomed in reciprocal space (wave plane basis)
        

        ! ERROR ANALYSIS PARAMETERS

        toll = 1.d-7  ! Convergence tolerance: lower limit of absolute distance of two eigenvalues
        maxInd = 10     ! Maximum index have two different meaning:
                        !     - convergence VS eigenvalue index: highest index for which the test is executed.
                        !     - convergence VS length: index of actual eigenvalue used for the computation.
        maxPotPerc = 10! Maximum potential percentage: smallest value of L_min used for the computation. 
                        ! L_min satifies the relation: V(L_min) = -De*(maxPotPerc/100)

        CALL ERROR_ANALYSIS(toll,"K",maxPotPerc,maxInd)    
                                               
END PROGRAM
