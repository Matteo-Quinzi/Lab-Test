module mirror 
    !
    use, intrinsic :: iso_c_binding
    !
    implicit none
    !
    include 'fftw3.f03'
    !
    real(kind=8), private :: pi = acos(-1.d0)
    integer(kind=8), private :: i, j, k
    !
    contains
    !
    !
    !
    subroutine do_fft(fd1, fd2, forward)
        ! Compute the fast fourier transform of a 1D function.
        ! This subroutine works with the FFTW library.
        !
        ! Arguments :
        ! ------------
        ! fd1(double complex, dimension N ) : input function 
        ! forward (optional) : species (1) forward or (0) backward transformation
        ! 
        ! Outputs :
        ! ------------
        ! fd2 (double complex, dimension N ) : (anti)transformed function
        !
        complex(kind=8), dimension(:) :: fd1
        integer, optional, intent(in) :: forward   
        integer :: flow, N 
        complex(kind=8), dimension(:) :: fd2
        type(C_PTR) :: plan 
        !
        N = size(fd1)
        flow = 1
        if (present(forward)) flow = forward 
        if (flow == 1) then
            plan = fftw_plan_dft_1d(N, fd1, fd2, FFTW_forward, FFTW_estimate)
            call fftw_execute_dft(plan, fd1, fd2)
            call fftw_destroy_plan(plan)
            fd2 = fd2 / dble(N)            ! output is correctly renormalized
        else if (flow == 0) then 
            plan = fftw_plan_dft_1d(N, fd1, fd2, FFTW_backward, FFTW_estimate)
            call fftw_execute_dft(plan, fd1, fd2)
            call fftw_destroy_plan(plan)
            ! output here does not need renormalization
        else 
            stop 'Value Error : Illegal value of optional parameter forward'
        end if
    end subroutine do_fft
    !
    !
    !
    subroutine gaussian(x, std, mu, y)
        ! Compute a gaussian at points given by x with
        ! mean value mu and standard deviation std
        real(kind=8), dimension(:), intent(in) :: x 
        real(kind=8), intent(in) :: std, mu 
        complex(kind=8), dimension(:), intent(out) :: y
        integer :: N, i  
        N = size(x)
        do i = 1, N
            y(i) = 1 / ( std * sqrt(2.d0 * pi) ) * exp( - (x(i) - mu)**2.d0 / (2.d0* std**2.d0) ) 
        end do
    end subroutine gaussian
    !
    !
    !
    subroutine hamiltonian_k(k, V_k, H_k)
        ! Define the Hamiltonian in k space
        !
        ! Arguments:
        !------------
        ! k(double, dimension N) : array with wave momenta
        ! V_k(complex, dimension N) : array with potential in k_space
        !
        ! Outputs:
        ! H_k(complex, dimension N x N) : hermitian matrix describing the 
        !                                 k-soace hamiltonian (only the upper triangle
        !                                 is built up)
        real(kind=8), dimension(:), intent(in) :: k
        complex(kind=8), dimension(:), intent(in) :: V_k 
        complex(kind=8), dimension(:,:), intent(out) :: H_k 
        integer :: i,j, N 
        N = size(V_k) 
        do i = 1, N 
            do j = i, N 
                H_k(i,j) = V_k(j-i+1)
            end do
        end do
        do i = 1,N
            H_K(i,i) =  0.5 * k(i) ** 2.d0 + H_k(i,i) 
        end do
    end subroutine hamiltonian_k
    !
    !
    !
    subroutine solve_hamiltonian_k(N, H_k, eva_k)
        ! solve hamiltonian in k space using routine zheevr from 
        ! LAPACK library 
        ! 
        ! Arguments:
        ! -----------
        ! H_k (complex, N x N) : Hamiltonian in k_space
        ! N (integer) : matrix H order
        !
        ! Outputs:
        ! eva_k (real, M) : array with eigenvalues 
        integer, intent(in) :: N
        complex(kind=8), intent(inout) :: H_k(:,:) ! on exit the lower triangle is destroyed
        real(kind=8), intent(out) :: eva_k(:)
        !
        ! zheevr in parameters 
        character(1) :: jobz = 'v', uplo = 'u'
        integer :: lda , lwork, lrwork
        ! zheevr out parameters
        integer :: info
        real(kind=8), allocatable:: w(:), rwork(:)
        complex(kind=8), allocatable:: work(:)

        lda = N
        lwork = 2*N-1
        lrwork = 3*N-2
        !
        allocate( w(N), rwork(2*lrwork), work(2*lwork) )
        !
        !
        call zheev(jobz, uplo, N, H_k, lda, eva_k, work, lwork,&
        rwork, info)
        !
        !
        ! diagnostic printout 
        if (info == 0) then 
            !write(*,*) 'Subroutine zheevr : successful computation !'
        else if (info < 0 ) then 
            write(*,*) 'Subroutine zheevr : element i = ', -info, ' had an illegal value |:('
        else 
            write(*,*) 'Subroutine zheevr : algorithm failed to converge :( '
        end if 
    end subroutine solve_hamiltonian_k
    !
    !
    !
    subroutine m_to_h(x, x_0, D_e, alpha, h_potential)
        ! Given a Morse potential defines the harmonic potential 
        ! obtained by Taylor expanding near the energy minimum
        real(kind=8), dimension(:), intent(in) :: x
        real(kind=8), intent(in) :: x_0, D_e, alpha ! Morse potential parameters
        real(kind=8), dimension(:), intent(out) :: h_potential
        real(kind=8) :: mw2 ! harmonic constant 
        integer :: N 
        N = size(h_potential)
        mw2 = D_e * alpha ** 2.d0 
        h_potential = (/ (mw2 * (x(i) - x_0)**2.d0 - D_e , i=1,N )/)
    end subroutine m_to_h
    !
    !
    !
    subroutine quadrature_1d(dx, f, integral)
        ! Perform quadrature in 1d using trapezoidal formula
        real(kind=8), dimension(:), intent(in) :: f(:)
        real(kind=8), intent(in) :: dx
        real(kind=8), intent(out) :: integral 
        integer :: N 
        N = size(f)
        integral = 0.5 * (f(1) + f(N)) 
        do i = 2,N-1
            integral = integral + f(i)
        end do
        integral = dx*integral
    end subroutine
    !
    !
    !
    subroutine hamiltonian_ho_space(dx, eigenvalues, eigenstates, V_eff, H)
        ! Defines the hamiltonian on a subset of the eigenstates of the
        ! harmonic oscillator
        !integer :: M ! number of base states
        real(kind=8), intent(in) :: dx 
        real(kind=8), dimension(:,:), intent(in) :: eigenstates
        real(kind=8), dimension(:), intent(in) :: V_eff, eigenvalues
        real(kind=8), dimension(:,:), intent(out) :: H 
        integer :: M, N, i, j, k 
        real(kind=8), dimension(:), allocatable :: temp_intg ! temporary array for integration 
        M = size(eigenstates(1,:))
        N = size(eigenstates(:,1))
        allocate(temp_intg(N))
        do i = 1,M
            do j=i,M
                ! define integral as <i| V_eff |j>
                temp_intg = (/ ( eigenstates(k,i) * V_eff(k) * eigenstates(k,j) , k=1,N)   /)
                call quadrature_1d(dx, temp_intg, H(i,j))
            end do 
        end do
        !
        ! adding values on the diagonal
        do i = 1,M
            H(i,i) = H(i,i) + eigenvalues(i)
        end do
    end subroutine hamiltonian_ho_space
    !
    !
    !
    subroutine solve_hamiltonian_real(H, eigenvalues)
        ! Solve eigenvalues problem for a real hamiltonian
        real(kind=8), dimension(:,:), intent(inout) :: H 
        real(kind=8), dimension(:), intent(out) :: eigenvalues 
        integer :: N 
        !
        ! internal dsyev variables 
        character(1) :: jobz='v', uplo='u'
        integer :: lda 
        real(kind=8), dimension(:), allocatable :: work 
        integer :: lwork, info 
        N = size(H(:,1))
        lda = max(1,N) + 1 
        lwork = max(1, 3*N - 1) + 1 
        allocate(work(lwork))
        ! 
        ! calling dsyev 
        call dsyev(jobz, uplo, N, H, N, eigenvalues, work, lwork, info)
        !
        !
        ! diagnostic printout 
        if (info == 0) then 
            write(*,*) 'Subroutine dsyev : successful computation !'
        else if (info < 0 ) then 
            write(*,*) 'Subroutine dsyev : element i = ', -info, ' had an illegal value |:('
        else 
            write(*,*) 'Subroutine dsyev : algorithm failed to converge :( '
        end if 
    end subroutine solve_hamiltonian_real
    !
    !
    !
    subroutine rebuild_wave(x, k, A_k, f)
        ! The aim of this routine is to rebuild the eigenstates in
        ! x_space to compare them with the eigenstates directly built in
        ! that space
        !
        ! Arguments:
        !------------
        ! x(double) (dimension N) : array of coordinates
        ! k(double) (dimension N_k) : array of wavevectors, determining which waves to use
        ! A_k(complex) (dimension N_k) : array of the waves coefficients
        !
        ! Output:
        ! --------
        ! f(complex) (dimension N) : array with the rebuilt function in x space
        real(kind=8), dimension(:), intent(in) :: x, k
        complex(kind=8), dimension(:), intent(in) :: A_k 
        complex(kind=8), dimension(:), intent(out) :: f 
        complex(kind=8), dimension(:,:), allocatable :: T ! auxiliary matrix 
        integer :: N, N_k, i, j
        N = size(x)
        N_k = size(k)
        allocate( T(N_k, N) )
        ! Building up a matrix which has a plane wave on each row 
        ! weighted by the coefficient A_k
        do i = 1, N_k
            T(i,:) = (/ ( A_k(i) * zexp(-complex(0.d0,k(i)*x(j))) , j=1,N ) /)
        end do
        !
        ! Summing up all the waves pointwise
        f = (/ (sum(T(:,i)) , i=1,N) /)
    end subroutine
    !
    !
    !
    subroutine absolute_square(V, A)
        ! Given an array with complex entries, the 
        ! routine builds up an array with the absolute square 
        ! of each entry
        complex(kind=8), dimension(:), intent(in) :: V
        real(kind=8), dimension(:), intent(out) :: A 
        integer :: N 
        N = size(V)
        A = (/ (real(V(i))**2.d0 + aimag(V(i))**2.d0 , i=1,N)   /)
    end subroutine absolute_square
    !
    !
    !
    subroutine print_matrix(M, N)
        !diagnostic subroutine
        !print NxN matrix in a nice format
        complex(kind=8), dimension(:,:), intent(in) :: M 
        integer, intent(in) :: N 
        integer :: i,j
        do i = 1,N
                write(*,*) M(i,:) 
        end do
    end subroutine
end module mirror