program tester 
    use eigen
    use mirror 
    !
    implicit none 
    !
    integer, parameter :: io_unit = 12
    real(kind=8) :: pi = acos(-1.d0)
    integer :: i, j, M = 10
    ! Real space variables 
    integer :: N_x
    real(kind=8) :: L, x_0, alpha, D_e
    real(kind=8), allocatable :: x(:), dx, d(:), e(:), V_re_x(:), eigenvals(:), eigenvecs(:,:), absolute_squares(:,:)
    complex(kind=8), allocatable :: V_x(:), at_V_x(:), fd2(:)
    complex(kind=8), allocatable :: eve_test(:)
    !
    ! Reciprocal space variables 
    integer :: N_k 
    real(kind=8) :: dk
    real(kind=8), allocatable :: k(:), eigenvals_k(:)
    complex(kind=8), allocatable :: V_K(:), H_k(:,:)
    !
    !
    ! Acquiring input data
    open(unit = io_unit, file='input.dat', status='old')
        read( io_unit , *) 
        read( io_unit , fmt=*) L, x_0, D_e, alpha, N_x, N_k 
    close(io_unit) 
    M = N_k
    !
    !
    ! Real space eigenenergies and eigenfunctions computation 
    ! 
    ! Allocating memory ... 
    allocate(x(N_x), d(N_x), e(N_x - 1), V_x(N_x), V_re_x(N_x))
    !
    !Defining x space potential 
    dx = L / dble(N_x)
    call linspace(0.d0, L, N_x, x)
    !do i = 1,N_x
    !    if ( x(i) <= 10.d0 ) then
    !        V_re_x(i) = 0.d0
    !    else if ( x(i) <= 11.d0) then
    !        V_re_x(i) = + 5.d0 
    !    else 
    !        V_re_x(i) = 0.d0
    !    end if
    !end do
    call Morse(x, x_0, D_e, alpha, V_re_x)
    !
    ! Defining x space Hamiltonian and solving eigenvalues problem
    ! (Only the first N_k eigenvalues are computed)
    allocate(eigenvals(M), eigenvecs(N_x, M))
    call Hamiltonian_x(x, d, e, V_re_x)
    call solve_eigens(M, d, e, 'v', eigenvals, eigenvecs)
    eigenvecs = eigenvecs / sqrt(dx)
    !
    ! Saving eigenvalues and eigenvectors 
    write(*,*) ' Computation in x space completed '
    write(*,*) ' Saving eigenvalues in eva_x.dat '
    open( io_unit, file='eva_x.dat')
        do i = 1, M
            write(io_unit, *) i, eigenvals(i)  
        end do
    close(io_unit)
    !
    write(*,*) ' Saving eigenvectors in eve_x.dat '
    open( io_unit, file='eve_x.dat')
    do i = 1, N_x 
        write(io_unit, *) i, x(i), V_re_x(i),  eigenvecs(i,:)
    end do
    close( io_unit )
    !
    !
    write(*,*) ' Starting computation in k space'
    !
    ! Defining complex V_x potential 
    do i = 1, N_x
        V_x(i) = complex( V_re_x(i) , 0.d0 ) 
        !if ( x(i) <= (x_0 - log(2.d0)/alpha) ) then 
        !    V_x(i) = 0.d0 
        !else 
        !    V_x(i) = D_e * ( ( 1 - exp(- alpha*(x(i)-x_0) )) ** 2.d0 - 1.d0 )
        !end if 
        !V_x(i) = 0.5 * x(i) ** 2.d0
        !
    end do
    !
    ! Fourier Transforming V_x potential
    ! Anti transform check is added 
    allocate(at_V_x(N_x), fd2(N_x))
    call do_fft(V_x, fd2) 
    call do_fft(fd2, at_V_x, 0)
    !
    !
    write(*,*) ' Saving Fourier transform results in four.dat'
    open(io_unit, file='four.dat')
        do i = 1, N_x
            write(12, fmt='(8e15.7)') x(i), V_x(i), 2*pi/L *i , L*fd2(i), at_V_x(i)  
        end do
    close(io_unit)
    !
    !
    ! Defining k wavevectors (only positive are recquired)
    allocate(k(N_k))
    dk = 2*pi/L
    do i = 1,N_k
        k(i) = dk * (i - 1 - N_k/2)
    end do
    !
    ! Extracting Transformed potential
    allocate(V_k(N_k)) 
    V_k(1 : N_k) = fd2( 1 : N_k)
    !
    ! Defining Hamiltonian in k space 
    allocate(H_k(N_k, N_k), eigenvals_k(N_k))
    call hamiltonian_k(k, V_k, H_k)
    !call print_matrix(H_k, N_k)
    call solve_hamiltonian_k(N_k, H_k, eigenvals_k)
    !
    ! rebuilding eigenstates in x space for comparison
    allocate(eve_test(N_x), absolute_squares(N_x,2))
    call rebuild_wave(x, k, H_k(:,5), eve_test)
    call absolute_square(eve_test, absolute_squares(:,1))
    open(io_unit, file='eve_test.dat')
        do i = 1,N_x
            write(io_unit, fmt='(2e15.7)') absolute_squares(i,1)/L, eigenvecs(i,5)**2.d0
        end do
    close(io_unit)
    !
    ! 
    write(*,*) ' Computation in k space completed '
    !
    !
    write(*,*) ' Saving eigenvalues in eva_k.dat'
    open(io_unit, file='eva_k.dat')
        do i=1, N_k
            write(io_unit, *) i, eigenvals_k(i)
        end do
    close(io_unit)
    !
    !
    write(*,*) ' Saving eigenvectors in eve_k.dat'
    open(io_unit, file='eve_k.dat')
        do i=1,N_k
            write(io_unit,*) i, H_k(i,:)
        end do
    close(io_unit)
    !
    !
    write(*,*) ' Execution terminated'
    write(*,*) ' See you soon ;) '
end program tester 