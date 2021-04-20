PROGRAM MAIN
        USE READ_VAR
        USE TASK
        IMPLICIT NONE

        REAL(KIND=8), ALLOCATABLE :: x(:), V(:), k(:)
        COMPLEX(KIND=8), ALLOCATABLE :: F_X(:), F_K(:)
        REAL(KIND=8), ALLOCATABLE :: F_evol(:,:)
 
        CALL READ_DATA()
        ALLOCATE(x(N), V(N), k(N))
        ALLOCATE(F_X(N), F_K(N), F_evol(N,save_wave + 1))
        
        x(:) = (/ (i * L / (N - 1), i = 0, N-1) /)
        
        k( : N/2) = (/ (2.D0 * pi * i / L, i = 0, N/2-1) /)
        k(N/2 + 1 : ) = (/ (-2.D0 * pi / L * (N/2.D0 - i), i=0, N/2-1) /) 

        V = POT(x)
        
        F_X = PHI(x)
        F_evol(:,1) = REAL(F_X(:))**2 + AIMAG(F_X(:))**2

        save_timestep = M / save_wave 
        DO i = 1, M
            F_X(:) = F_X(:) * EXP(-im * V(:) * dt/2.D0)
            
            F_K = FFT(F_X, 'F') 
            F_K(:) = F_K(:) * EXP(-im * k(:)**2 * dt/2.D0) / N

            F_X = FFT(F_K, 'B')
            F_X(:) = F_X(:) * EXP(-im * V(:) * dt/2.D0)

            IF (MOD(i, save_timestep) == 0) F_evol(:, i / save_timestep + 1) = REAL(F_X(:))**2 + AIMAG(F_X(:))**2
        END DO
        
        OPEN(UNIT=10, FILE='Output/moving.txt')
        WRITE(10, '(202F15.7)') (x(i), F_evol(i,:), i=1,N)
        CLOSE(10)

        OPEN(UNIT=10, FILE='Output/pot.txt')
        WRITE(10, '(2F20.7)') (x(i), V(i), i=1,N)
        CLOSE(10)

END PROGRAM MAIN

