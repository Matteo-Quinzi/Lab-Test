PROGRAM MAIN
        USE READ_VAR
        USE TASK
        IMPLICIT NONE

        REAL(KIND=8), ALLOCATABLE :: x(:), V(:), k(:)
        COMPLEX(KIND=8), ALLOCATABLE :: F_X(:), F_K(:), F_evol(:,:)
 
        CALL READ_DATA()
        ALLOCATE(x(N), V(N), k(N))
        ALLOCATE(F_X(N), F_K(N), F_evol(N,save_wave + 1))
        
        x(:) = (/ (i * L / (N - 1), i = 0, N-1) /)
        !k(:) = (/ (2.D0 * pi * i / L, i = 0, N-1) /)
        k(:) = (/ ((2 * i - N) * pi / L, i = 1, N) /)

        V = POT(x)
        !V(:) = 0.D0
        !V(6000:) = 1D10

        F_X = PHI(x)
        F_evol(:,1) = F_X(:)

        save_timestep = M / save_wave 
        DO i = 1, M
            F_X(:) = F_X(:) * EXP(-im * V(:) * dt/2.D0)
            
            F_K = FFT(F_X, 'F') 
            F_K(:) = F_K(:) * EXP(-im * k(:)**2 * dt/2.D0) / N

            F_X = FFT(F_K, 'B')
            F_X(:) = F_X(:) * EXP(-im * V(:) * dt/2.D0)

            IF (MOD(i, save_timestep) == 0) F_evol(:, i / save_timestep + 1) = F_X(:)
        END DO
        
        OPEN(UNIT=10, FILE='moving.txt')
        WRITE(10, '(203F15.7)') (x(i), F_evol(i,:), i=1,N)
        CLOSE(10)

        OPEN(UNIT=10, FILE='pot.txt')
        WRITE(10, '(2F20.7)') (x(i), V(i), i=1,N)
        CLOSE(10)

END PROGRAM MAIN


