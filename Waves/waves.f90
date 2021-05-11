PROGRAM MAIN
        USE READ_VAR
        USE TASK
        USE OMP_LIB
        IMPLICIT NONE

        INTEGER :: l2, l1, cc=0, cd=0, i_ent, i_exi, j
        REAL(KIND=8) :: peak, x_peak, ent, exi, t_ent, t_exi, time_hartman, T_B(15,4)

        REAL(KIND=8), ALLOCATABLE :: x(:), V(:), k(:), F_EVOL(:,:), F_SQMOD(:)
        COMPLEX(KIND=8), ALLOCATABLE :: V_OP(:), K_OP(:), F_X(:), F_K(:)

        CALL READ_DATA()
        ALLOCATE(x(N), V(N), k(N))
        ALLOCATE(V_OP(N), K_OP(N))
        ALLOCATE(F_X(N), F_K(N), F_SQMOD(N), F_EVOL(N, save_wave+1))
        
        x(:) = (/ (i * L / (N - 1), i = 0, N-1) /)
        
        k( : N/2) = (/ (2.D0 * pi * i / L, i = 0, N/2-1) /)
        k(N/2 + 1 : ) = (/ (-2.D0 * pi / L * (N/2.D0 - i), i=0, N/2-1) /) 

        DO j = 1,15
        V = POT(x)

        V_OP(:) = EXP(-im * V(:) * dt / 2.D0)
        K_OP(:) = EXP(-im * k(:)**2 * dt / 2.D0)

        !$OMP PARALLEL PRIVATE(F_X, F_K) SHARED(F_EVOL)

        F_X = PHI(x)

        F_EVOL(:,1) = REAL(F_X(:))**2 + AIMAG(F_X(:))**2

        PRINT *, 'INITIAL NORM:', SUM(F_EVOL(:,1)) * (x(2)-x(1))

        save_timestep = M / save_wave

        l1 = INT(N * l_1 / (2 * L))
        l2 = INT(N * l_2 / (2 * L))

        ent = 0.
        exi = 0.

        cc = 0
        cd = 0
        !$OMP DO
        DO i = 1, M
            F_X(:) = F_X(:) * EXP(-im * V(:) * dt/2.D0)
            
            F_K = FFT(F_X, 'F') 
            F_K(:) = F_K(:) * EXP(-im * k(:)**2 * dt/2.D0) / N

            F_X = FFT(F_K, 'B')
            F_X(:) = F_X(:) * EXP(-im * V(:) * dt/2.D0)

            F_SQMOD(:) = REAL(F_X(:))**2 + AIMAG(F_X(:))**2

            IF (F_SQMOD(N/2-l1+1) > ent) THEN 
                    ent = F_SQMOD(N/2-l1+1)
                    t_ent = i*dt
                    i_ent = i
            ENDIF

            IF (F_SQMOD(N/2+l1+1) > exi) THEN 
                    exi = F_SQMOD(N/2+l1+1)
                    t_exi = i*dt
                    i_exi = i
            ENDIF

            IF (MOD(i, save_timestep) == 0) F_EVOL(:, i / save_timestep + 1) = F_SQMOD(:)

            IF (cc == 0) THEN
                    IF (F_SQMOD(N/2-l1+1) > 0.0001) THEN 
                            time_hartman = i*dt
                            cc = 1
                            PRINT *, '1',i
                    ENDIF
            ENDIF
                 
            IF (cd == 0) THEN
                    IF (F_SQMOD(N/2+l1-1) > 0.0001) THEN
                            time_hartman = i*dt - time_hartman
                            cd = 1
                            PRINT *, '2', i
                    ENDIF
            ENDIF
 
        END DO
        !$OMP END DO
        !$OMP END PARALLEL
        
        l_1 = l_1 + 1.
        x0 = x0 - 0.5

        T_B(j,1) = l_1
        T_B(j,2) = time_hartman
        T_B(j,3) = t_ent
        T_B(j,4) = t_exi

        END DO

        OPEN(UNIT=10, FILE='hartman.txt')
                WRITE(10, '(4F15.7)') (T_B(i,:), i=1,15)
        CLOSE(10)

        PRINT *, 'FINAL NORM:  ', SUM(F_EVOL(:,save_wave +1)) * (x(2)-x(1))
        PRINT *, ''
        PRINT *, 'TIME DELAY:', t_exi - t_ent
        PRINT *, 'ENTERED AT', i_ent, ', EXITED AT', i_exi

        CALL WRITE_INTO_FILE()
        
        PRINT '(A, F7.5, A)', 'R  = ', SUM(F_EVOL(           : N/2 - l1   , save_wave+1)) * (x(2)-x(1)), ' %'
        PRINT '(A, F7.5, A)', 'A1 = ', SUM(F_EVOL(N/2 - l1   : N/2 + l1   , save_wave+1)) * (x(2)-x(1)), ' %'
        PRINT '(A, F7.5, A)', 'A12= ', SUM(F_EVOL(N/2 + l1   : 2*N/3 - l2 , save_wave+1)) * (x(2)-x(1)), ' %'
        PRINT '(A, F7.5, A)', 'A2 = ', SUM(F_EVOL(2*N/3 - l2 : 2*N/3 + l2 , save_wave+1)) * (x(2)-x(1)), ' %'
        PRINT '(A, F7.5, A)', 'T  = ', SUM(F_EVOL(2*N/3 + l2 :            , save_wave+1)) * (x(2)-x(1)), ' %'

        peak = 0.D0
        x_peak = 0.D0

        DO i = 1, N
          IF (i >= N/2 + l2) THEN
                  IF (F_EVOL(i, save_wave+1) > peak) THEN
                          peak = F_EVOL(i, save_wave+1)
                          x_peak = x(i)
                  END IF
         END IF
        END DO

        PRINT *, 'Found peak at', x_peak

        PRINT *, 'Hartman =', time_hartman


          
        CONTAINS 
                SUBROUTINE WRITE_INTO_FILE()
                        CHARACTER(20) :: FMT

                        WRITE(FMT, '("(", I0, "F20.10)")') save_wave + 1
                        OPEN(UNIT = 10, FILE = 'Output/moving.txt')
                          WRITE(10, FMT) (F_evol(i,:), i=1,N)
                        CLOSE(10)

                        OPEN(UNIT = 10, FILE = 'Output/final_wave.txt')
                          WRITE(10, '(2F20.10)') (F_X(i), i = 1, N)
                        CLOSE(10)
                END SUBROUTINE WRITE_INTO_FILE

END PROGRAM MAIN

