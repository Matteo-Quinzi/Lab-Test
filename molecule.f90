PROGRAM molecule 
    !
    USE eigen
    ! 
    IMPLICIT NONE 
    !
    REAL(double) :: L, time, alpha, x_0, x_min, D_morse
    INTEGER :: N, M, i, j 
    CHARACTER(1) :: jobz
    CHARACTER(50) :: fmt
    REAL(double), ALLOCATABLE :: eigenvals(:), eigenvecs(:,:), x(:), d(:), e(:), potential(:)
    INTEGER :: io_unit = 12
    !
    x_0 = 2.0d0
    x_min = 0.d0
    D_morse = 50.0
    !
    OPEN(io_unit, file='input.txt', STATUS='old')
       WRITE(*,*) "Reading from 'input.txt'. "
       ! Read the interval length L, the number of wanted eigenvalues M, the number of points 
       ! on the grid N, and the morse potential parameter alpha
       READ(io_unit,*) L, M, N, jobz, alpha         
    CLOSE(io_unit)
    !
    ALLOCATE(x(N), d(N), e(N-1), eigenvals(M), eigenvecs(N, M), potential(N))        
    !
    CALL linspace(x_min, L, N, x)  ! define points on x axis
    !
    CALL Morse(x, x_0, D_morse, alpha, potential) ! define potential on x space
    ! 
    CALL Hamiltonian(x, d, e, potential) ! define the hamiltonian
    !
    CALL solve_eigens(M, d, e, jobz, eigenvals, eigenvecs, time)  ! solve eigenproblem for the hamiltonian
    !
    CALL rescale_eigens(eigenvals, eigenvecs)  ! rescale the results
    !
    WRITE(*,*) 'Calculation in x space completed !'
    !
    OPEN(io_unit, file='eigenvals.txt')
        WRITE(*,*) " Saving eigenvalues in 'eigenvals.txt'. "
        DO i = 1, M 
            WRITE(io_unit,'(I0, F18.10)') i, eigenvals(i)
        END DO
    CLOSE(io_unit)
    !
    IF ( jobz == 'V' ) THEN
        OPEN(io_unit, file='waves.txt')
        WRITE(*,*) " Saving eigenvectors in 'waves.txt'. "
            DO i = 1, N
                WRITE(io_unit, *) x(i), potential(i), eigenvecs(i, :)
            END DO 
        CLOSE(io_unit)
    END IF
    !
    !TASK B
    !
    
END PROGRAM molecule