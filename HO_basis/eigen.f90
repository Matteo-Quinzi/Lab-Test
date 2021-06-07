MODULE eigen 
    !
    ! Contains subroutines to produce the 1-D quantum harmonic
    ! oscillator's Hamiltonian in coordinates representation
    ! and solve the eigenvalues problem.
    !
    IMPLICIT NONE
    INTEGER, PRIVATE :: i,j             ! loop variables
    REAL(KIND=8), PRIVATE :: h = 1.d0    ! distance between two successive 
                                         ! coordinates values
    SAVE
    !
    INTEGER, PARAMETER :: single = SELECTED_REAL_KIND(p=6, r=37)
    INTEGER, PARAMETER :: double = SELECTED_REAL_KIND(p=13, r=200)
    !
    !
    !
    !
    CONTAINS 
    !
    !
    !
    SUBROUTINE linspace(start, end, N, x, close_in)
        !
        ! Build up an array of N equispaced values in
        ! the interval [start, end]
        !
        ! Arguments:
        ! start (real)(double precision) : first value in the array
        ! end   (real)(double precision) : last value in the array
        ! N     (integer) : number of points in the interval
        ! close_in (logical)(optional) : if .TRUE. extremes are included  
        ! 
        ! Outputs :
        ! x (real)(double)(dimension(N)) : array with equispaced values
        !                                  in the chosen interval 
        !
        REAL(double), INTENT(IN) :: start, end
        INTEGER, INTENT(IN) :: N 
        LOGICAL, OPTIONAL, INTENT(IN) :: close_in
        REAL(double), DIMENSION(:), INTENT(OUT) :: x
        LOGICAL :: close = .FALSE.
        !
        IF ( SIZE(x) /= N ) THEN 
          STOP "Error raised in subroutine linspace: &
          & x's dimension differs from number of points"
        END IF
        IF (PRESENT(close_in)) close = close_in
        IF (close .EQV. .FALSE.) THEN
          h = (end - start)/(N + 1)
          DO i = 1, N
            x(i) = start + i*h
          END DO
        ELSE 
          h = (end - start)/(N - 1)
          DO i = 0, N-1
            x(i+1) = start + i*h 
          END DO
        END IF
    END SUBROUTINE linspace 
    !
    !
    !
    SUBROUTINE Morse(x , x_0, D_morse, alpha, potential)
        !
        ! Compute Morse potential on selected x points
        !
        ! Arguments :
        ! x (double)(dimension N) : array of coordinates points
        ! x_0 (double) : coordinate of minimum
        ! D (double) : potential well depth
        ! alpha (double) : Morse parameter 
        !
        ! Outputs :
        ! potential(double)(dimension N): array with potential values
        ! 
        REAL(double), DIMENSION(:), INTENT(IN) :: x
        REAL(double), INTENT(IN) :: x_0, alpha, D_morse
        REAL(double), DIMENSION(:), INTENT(OUT) :: potential
        INTEGER :: s 
        s = SIZE(x)
        DO i=1,s
            potential(i) = D_morse * ( ( 1 - EXP( - alpha * ( x(i) - x_0 ) ) )**2.d0 - 1.d0 )
        END DO 
    END SUBROUTINE    
    !
    !
    !
    SUBROUTINE Hamiltonian_x(x, d, e, potential)
        !
        ! Build up the tridiagonal Hamiltonian of
        ! the 1-D quantum harmonic oscillator in 
        ! coordinates representation
        !
        ! Arguments :
        ! x (double)(dimension N) : array with points in 
        !                           the chosen interval
        ! potential(double)(dimension N) : array with potential 
        !                                  computed on x values
        !
        ! Outputs :
        ! d(double)(dimension N) : array with values on the 
        !                           Hamiltonian diagonal
        ! e(double)(dimension N-1) : array with values out of
        !                            the Hamiltonian diagonal
        !
        REAL(double), INTENT(IN) :: x(:), potential(:)
        REAL(double), INTENT(OUT) :: d(:), e(:)
        INTEGER :: N 
        !
        IF (SIZE(d) /= SIZE(x)) THEN 
            STOP "Error raised in subroutine hamiltonian: &
            &d's dimension must be equal to x's"
        ELSE IF (SIZE(e) /= SIZE(d)-1) THEN
            STOP "Error raised in subroutine hamiltonian: &
            &e's dimension must be 1 less then d's"
        END IF
        N = SIZE(d)
        DO i = 1, N
            d(i) = 1.d0 / h**2.d0 + potential(i) 
        END DO
        DO i = 1, N-1
            e(i) = - 1/(2.d0 * h**2.d0) 
        END DO
    END SUBROUTINE Hamiltonian_x
    !
    !
    !
    subroutine print_matrix(M, N)
        !diagnostic subroutine
        !print NxN matrix in a nice format
        complex(kind=8), dimension(:,:), intent(in) :: M 
        integer, intent(in) :: N 
        integer :: i
        do i = 1,N
                write(*,*) M(i,:) 
        end do
    end subroutine
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
    SUBROUTINE solve_eigens(M, d, e, jobz, eigenvals, eigenvecs, time)
        !
        ! Solve the eigenvalues problem of a tridiagonal
        ! matrix using lapack library's dstevr() subroutine.
        ! Returns only the first M eigenvalues and, optionally,
        ! the relative eigenvectors.
        !
        ! Arguments :
        ! M (integer) : number of eigenvalues
        ! D (real)(double)(dimension(N)) : array with the diagonal values
        ! E (real)(double)(dimension(N-1)) : array with out of diagonal elements
        ! jobz (character(1)) : if jobz='V' computes also eigenvectors, else 
        !                       if jobz='N' computes only eigenvalues
        !
        ! Outputs :
        ! eigenvals (real)(double)(dimension(M)) : contains the first M
        !                                          computed eigenvalues
        ! eigenvecs (real)(double)(dimension(N,M)) : if jobz='V' contains
        !                         the first M eigenvectors on each column 
        ! time (real)(double)(optional) : time passed in dstevr subroutine
        !
        ! For a detailed description of the dstevr() subroutine parameters
        ! look at the official manual page on netlib.org
        ! http://www.netlib.org/lapack/explore-html/dc/dd2/group__double_o_t_h_e_reigen_ga323734560b8bd052fbc474dc2f0b5605.html
        !
        INTEGER, INTENT(IN) :: M 
        REAL(double), INTENT(IN) :: d(:), e(:)
        CHARACTER(1), INTENT(IN) :: jobz 
        REAL(double), INTENT(OUT) :: eigenvals(:), eigenvecs(:,:)
        REAL(double), OPTIONAL, INTENT(OUT) :: time
        CHARACTER(1) :: range
        INTEGER :: N, il, iu, m_out, info, liwork, lwork
        INTEGER, DIMENSION(:), ALLOCATABLE :: iwork, isuppz 
        REAL(double) :: abstol, vl, vu, tic, toc, DLAMCH
        REAL(double), DIMENSION(:), ALLOCATABLE :: w, work
        !
        IF (SIZE(eigenvals) /= M) THEN
            STOP "Error raised in subroutine solve_eigens: &
            &eigenvals's size must be equal to M"
        ELSE IF (SIZE(eigenvecs,2) /= M) THEN
            STOP "Error raised in subroutine solve_eigens: &
            &eigenvecs must have M columns"
        ELSE IF (SIZE(eigenvecs,1) /= SIZE(d)) THEN
            STOP "Error raised in subroutine solve_eigens: &
            &eigenvecs must have as rows as d"
        END IF 
        !
        N = SIZE(d)
        range = 'I'             ! compute eigenvals in an indexed interval
        il = 1                  ! first index
        iu = M                  ! last index
        abstol = 2*DLAMCH('S')  ! absolute tolerance
        liwork = 10*N + 1 
        lwork = 20*N + 1
        !
        ALLOCATE(w(N), work(lwork), iwork(liwork), isuppz(lwork))
        !
        CALL CPU_TIME(tic)
        CALL dstevr(jobz, range, N, d, e, vl, vu, il, iu, abstol, m_out, &
          w, eigenvecs, N, isuppz, work, lwork, iwork, liwork, info)
        CALL CPU_TIME(toc)
        IF (PRESENT(time)) time = toc - tic
        !
        IF (info > 0) THEN
            STOP "Error raised in dstevr subroutine: &
            &internal error"
        ELSE IF (info < 0) THEN
            WRITE(*,'(A,I0)') 'i = ', -info
            STOP "Error raised in dstevr subroutine: &
            &there was an illegal value for index i"
        END IF
        !
        DO i = 1, M
            eigenvals(i) = w(i)
        END DO
    END SUBROUTINE solve_eigens 
    !
    !
    !
    SUBROUTINE rescale_eigens(eigenvals, eigenvecs)
        !
        ! Rescale eigenvalues and, optionally, eigenvectors
        ! with appropriate factors
        !
        REAL(double), INTENT(INOUT) :: eigenvals(:)
        REAL(double), OPTIONAL, INTENT(INOUT) :: eigenvecs(:,:)
        INTEGER :: N, M 
        M = SIZE(eigenvals)
        DO i = 1, M
            eigenvals(i) = 0.5d0 * eigenvals(i) 
        END DO
        IF (PRESENT(eigenvecs)) THEN 
            N = SIZE(eigenvecs, 1)
            DO i = 1, N 
                DO j = 1, M
                    eigenvecs(i, j) = 1/SQRT(h) * eigenvecs(i, j)
                END DO
            END DO
        END IF
    END SUBROUTINE rescale_eigens
    !
    !
END MODULE 
