PROGRAM MAIN
        USE READ_VAR
        USE TASK
     
        IMPLICIT NONE

        INTEGER :: l2, l1
        REAL(KIND=8), ALLOCATABLE :: x(:), V(:), k(:), PHI_EVOL(:,:), PHI_SQMOD(:)
        COMPLEX(KIND=8), ALLOCATABLE :: V_OP(:), K_OP(:), PHI_X(:), PHI_K(:)
        
        CALL TITLE()

        CALL READ_DATA()
        ALLOCATE(x(N), V(N), k(N))
        ALLOCATE(V_OP(N), K_OP(N))
        ALLOCATE(PHI_X(N), PHI_K(N), PHI_SQMOD(N), PHI_EVOL(N, save_wave+1))
        
        x(:) = (/ (i * L / (N - 1), i = 0, N-1) /)
        
        k( : N/2) = (/ (2.D0 * pi * i / L, i = 0, N/2-1) /)
        k(N/2 + 1 : ) = (/ (-2.D0 * pi / L * (N/2.D0 - i), i=0, N/2-1) /) 
        K_OP(:) = EXP(-im * k(:)**2 * dt / 2.D0)
       
        V = POT()
        V_OP(:) = EXP(-im * V(:) * dt / 2.D0)
        l1 = INT(N * l_1 / (2 * L))
        l2 = INT(N * l_2 / (2 * L)) 
 
        PHI_X = PHI(x)
        PHI_EVOL(:,1) = REAL(PHI_X(:))**2 + AIMAG(PHI_X(:))**2

        save_timestep = M / save_wave

        PRINT *, 'Computing wavepacket time evolution...'
        
        DO i = 1, M
            PHI_X(:) = PHI_X(:) * V_OP(:)
                        
            PHI_K = FFT(PHI_X, 'F') 
            PHI_K(:) = PHI_K(:) * K_OP(:) / N

            PHI_X = FFT(PHI_K, 'B')
            PHI_X(:) = PHI_X(:) * V_OP(:)

            PHI_SQMOD(:) = REAL(PHI_X(:))**2 + AIMAG(PHI_X(:))**2

            IF (MOD(i, save_timestep) == 0) PHI_EVOL(:, i / save_timestep + 1) = PHI_SQMOD(:) 
        END DO
        
        PRINT *, 'computation done.'
        PRINT *, ''

        CALL WRITE_INTO_FILE()

        PRINT *, 'The program has successfully completed the computation.'
        PRINT *, 'Good bye.'
        PRINT *, ''
        PRINT *, ''
        
        CONTAINS 
                SUBROUTINE TITLE()
                        PRINT *, ''
                        PRINT *, 'NUMERICAL SOLUTION OF TIME DEPENDENT SCHROEDINGER EQUATION'
                        PRINT *, '  QUANTUM WAVEPACKET TIME EVOLUTION'
                        PRINT *, ''
                        PRINT *, 'Code written by: Cuoghi A., Quinzi M., Rizzi G.'
                        PRINT *, ''
                        PRINT *, 'Dipartimento di Scienze Fisiche, Informatiche e Matematiche,'
                        PRINT *, 'Università di Modena e Reggio Emilia'
                        PRINT *, ''
                        PRINT *, 'May 2021'
                        PRINT *, ''
                        PRINT *, ''
                END SUBROUTINE TITLE
                
                SUBROUTINE WRITE_INTO_FILE()
                        CHARACTER(50) :: FMT
                        CHARACTER(14) :: f_pot  = 'Output/pot.txt'
                        CHARACTER(17) :: f_evol = 'Output/moving.txt'
                        CHARACTER(21) :: f_last = 'Output/final_wave.txt'
                        CHARACTER(16) :: f_coef = 'Output/coeff.txt'

                        PRINT *, 'Writing potential energy in "', f_pot, '"...'
                        OPEN(UNIT=10, FILE=f_pot)
                          WRITE(10, '(2F20.10)') (x(i), V(i), i = 1, N)
                        CLOSE(10)
                        PRINT *, ''

                        PRINT *, 'Saving wavepacket evolution in "', f_evol, '"...'
                        WRITE(FMT, '("(", I0, "F20.10)")') save_wave + 1
                        OPEN(UNIT=10, FILE=f_evol)
                          WRITE(10, FMT) (PHI_EVOL(i,:), i=1,N)
                        CLOSE(10)
                        PRINT *, ''

                        PRINT *, 'Saving last computed wavefunction in "', f_last, '"...'
                        OPEN(UNIT=10, FILE=f_last)
                          WRITE(10, '(2F20.10)') (PHI_X(i), i = 1, N)
                        CLOSE(10)
                        PRINT *, ''

                        PRINT *, 'Saving evolution parameters and computed coefficients in "', f_coef, '"...'
                        OPEN(UNIT=10, FILE=f_coef, POSITION='APPEND')
                         WRITE(10, '("L=", F20.10)') L
                         WRITE(10, '("E1=", F20.10, ", l1=", F20.10, ", E2=", F20.10, ", l2=", F20.10)') E_1, l_1, E_2, l_2
                         WRITE(10, '("s=", F20.10, ", x0=", F20.10, ", k0=", F20.10)') sigma, x0, k0
                         WRITE(10, '("M=", I10, ", dt=", F20.10)') M, dt
                         WRITE(10, *) 
                         WRITE(10, '(A, F15.10, A)') 'R  = ', SUM(PHI_EVOL(: N/2-l1, save_wave+1)) * (x(2)-x(1)) * 100.D0, ' %'
                         WRITE(10, '(A, F15.10, A)') 'A1 = ', SUM(PHI_EVOL(N/2-l1 : N/2+l1, save_wave+1)) * (x(2)-x(1)) * 100.D0, ' %'
                         WRITE(10, '(A, F15.10, A)') 'A12= ', SUM(PHI_EVOL(N/2+l1 : 2*N/3-l2, save_wave+1)) * (x(2)-x(1)) * 100.D0, ' %'
                         WRITE(10, '(A, F15.10, A)') 'A2 = ', SUM(PHI_EVOL(2*N/3-l2 : 2*N/3+l2, save_wave+1)) * (x(2)-x(1)) * 100.D0, ' %'
                         WRITE(10, '(A, F15.10, A)') 'T  = ', SUM(PHI_EVOL(2*N/3+l2 : , save_wave+1)) * (x(2)-x(1)) * 100.D0, ' %'
                         WRITE(10, *)
                         WRITE(10, '(A, F12.10, A)') 'Initial norm = ', SUM(PHI_EVOL(:, 1)) * (x(2) - x(1))
                         WRITE(10, '(A, F12.10, A)') 'Final norm = ', SUM(PHI_EVOL(:, save_wave+1)) * (x(2) - x(1))
                         WRITE(10, *)
                         WRITE(10, '(A)') '####################################'
                         WRITE(10, *)
                       CLOSE(10)
                       PRINT *, ''
                END SUBROUTINE WRITE_INTO_FILE
END PROGRAM MAIN

