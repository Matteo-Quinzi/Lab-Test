program ho_basis
    use eigen  
    use mirror
    !
    implicit none
    !
    integer, parameter :: N = 5000, M = 25, io_unit=12
    integer :: i
    real(kind=8), parameter :: L = 15.d0, D_e = 50.d0, alpha = 0.5d0, x_0 = 2.5d0
    real(kind=8) :: x(N), V_morse(N), V_harmonic(N), d_morse(N), e_morse(N-1), d_harmonic(N), e_harmonic(N-1) 
    real(kind=8) :: eigenvals_m(M), eigenvecs_m(N,M), eigenvals_h(M), eigenvecs_h(N,M)
    real(kind=8) :: dx, integral, test_fun(N), V_eff(N), H_ho(M,M)
    real(kind=8) :: eigenvals_ho(M), summ
    dx = L / (N+1) 
    !
    ! defining discretization in x space and potentials
    call linspace(0.d0, L, N, x)
    call Morse(x, x_0, D_e, alpha, V_morse)
    call m_to_h(x, x_0, D_e, alpha, V_harmonic)
    !
    ! defining Hamiltonians for Morse and harmonic potentials 
    call Hamiltonian_x(x, d_morse, e_morse, V_morse)
    call Hamiltonian_x(x, d_harmonic, e_harmonic, V_harmonic)
    !
    ! solve eigenvalues problem for both potentials 
    call solve_eigens(M, d_morse, e_morse, 'v', eigenvals_m, eigenvecs_m)
    call solve_eigens(M, d_harmonic, e_harmonic, 'v', eigenvals_h, eigenvecs_h)
    !
    ! rescaling eigenvecs 
    eigenvecs_m = 1/sqrt(dx) * eigenvecs_m
    eigenvecs_h = 1/sqrt(dx) * eigenvecs_h
    !
    ! quadrature test 
    test_fun = (/ (eigenvecs_m(i,1) * eigenvecs_m(i,1), i=1,N) /) 
    call quadrature_1d(dx, test_fun, integral)
    write(*,*) 'Integral is ', integral
    !
    !
    !
    write(*,*) 'Saving potentials in potentials.dat '
    open(io_unit, file='potentials.dat')
        do i = 1,N
            write(io_unit, fmt=*) x(i), V_morse(i), V_harmonic(i)  
        end do 
    close(io_unit)
    !
    !
    write(*,*) 'Saving eigenvalues in eva_mh.dat '
    open(io_unit, file='eva_mh.dat')
        do i = 1, M
            write(io_unit, fmt=*) eigenvals_m(i), eigenvals_h(i)
        end do
    close(io_unit)
    !
    !
    write(*,*) 'Saving Morse eigenvectors in eve_morse.dat'
    open(io_unit, file='eve_morse.dat')
        do i = 1, N
            write(io_unit, fmt=*) eigenvecs_m(i,:)
        end do
    close(io_unit)
    !
    !
    write(*,*) 'Saving Harmonic eigenvectors in eve_harmonic.dat'
    open(io_unit, file='eve_harmonic.dat')
        do i = 1, N
            write(io_unit, fmt=*) eigenvecs_h(i,:)
        end do
    close(io_unit)
    !
    ! the first M eigenvecs will form the basis in ho_space
    !
    ! Defining V_eff and Hamiltonian in ho space
    V_eff = (/ ( V_morse(i) - V_harmonic(i) , i =1,N )  /)
    call hamiltonian_ho_space(dx, eigenvals_h, eigenvecs_h, V_eff, H_ho)
    !
    ! Solving Hamiltonian in ho space 
    call solve_hamiltonian_real(H_ho, eigenvals_ho)
    !
    ! Saving results ...
    !
    !
    write(*,*) 'Saving eigenvalues and eigenvectors of ho space in eva_eve_ho.dat'
    open(io_unit, file='eva_eve_ho.dat')
        do i = 1,M
            write(io_unit, fmt=*) eigenvals_ho(i), H_ho(i,:)
        end do
    close(io_unit)
    !
    !
    summ = 0
    do i = 1, M 
        summ = summ + H_ho(i,1)*H_ho(i,1)
    end do 
    print*, summ

end program ho_basis