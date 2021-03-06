MODULE READ_VAR
        ! Declares the initial variables.
        ! Some of them are read from an external file.
        IMPLICIT NONE

        ! Global values accessible from the other programs:
        ! constant values and routine parameters
        REAL(KIND=8) :: pi = 4.D0 * DATAN(1.D0)  ! pi
        COMPLEX(KIND=8) :: im = (0.D0, 1.D0)  ! imaginary unit
        INTEGER :: i, save_timestep, IOS
        
        INTEGER :: l1, l2   ! barriers semilengths converted in array points unit

        ! values to be read from an external file
        INTEGER :: N, M, save_wave
        REAL(KIND=8) :: L, E_1, l_1, E_2, l_2
        LOGICAL :: load_f_inp
        CHARACTER(25) :: f_inp
        REAL(KIND=8) :: sigma, x0, k0, dt
        
        CONTAINS
                SUBROUTINE READ_DATA()
                        ! Read values from an external file.
                        IMPLICIT NONE

                        CHARACTER(20) :: f_data = 'Input/data_sheet.dat'  ! file name

                        PRINT *, 'Reading values from "', f_data, '"...'  
                        OPEN(UNIT=10, FILE=f_data, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
                        IF (IOS==0) THEN
                                READ(10, *) 
                                READ(10, '(3X, I10)') N  ! number of points of the interval
                                READ(10, '(3X, F20.10)') L  ! length of the interval
                                READ(10, *) 
                                READ(10, *) 
                                READ(10, '(5X, F20.10)') E_1  ! first barrier height
                                READ(10, '(5X, F20.10)') l_1  ! first barrier length
                                IF (l_1 > L) THEN
                                        PRINT *, 'Invalid value: l_1 is too large. Program stopping...'
                                        STOP
                                ELSEIF (l_1 < L/N) THEN
                                        PRINT *, 'WARNING: l_1 is less than the spatial points density L/N.'
                                        PRINT *, 'To turn off the barrier, please set E_1 = 0.'
                                ENDIF
                                READ(10, *) 
                                READ(10, '(5X, F20.10)') E_2  ! second barrier height
                                READ(10, '(5X, F20.10)') l_2  ! second barrier length
                                IF (l_2 > 2*L/3) THEN
                                        PRINT *, 'Invalid value: l_2 is too large. Program stopping...'
                                        STOP
                                ELSEIF (l_2 < L/N) THEN
                                        PRINT *, 'WARNING: l_2 is less than the spatial points density L/N.'
                                        PRINT *, 'To turn off the barrier, please set E_2 = 0.'
                                ENDIF
                                READ(10, *) 
                                READ(10, *)
                                READ(10, *)
                                READ(10, '(12X, L1)') load_f_inp  ! boolean for wavepacket shape
                                READ(10, *)
                                READ(10, *)
                                READ(10, '(7X, A)') f_inp  ! wavpeacket external file
                                READ(10, *)
                                READ(10, *)
                                READ(10, '(7X, F20.10)') sigma  ! gaussian wavepacket width
                                READ(10, '(4X, F20.10)') x0  ! gaussian wavepacket position
                                IF (x0 < 0 .OR. x0 > L) THEN
                                        PRINT *, 'Invalid value: x_0 has to be given a value between 0 and L. Program stopping...'
                                        STOP
                                ENDIF
                                READ(10, '(4X, F20.10)') k0  ! gaussian wavepacket velocity
                                READ(10, *)
                                READ(10, *)
                                READ(10, '(3X, I10)') M  ! number of timesteps
                                READ(10, '(4X, F20.10)') dt  ! timestep width
                                READ(10, '(11X, I5)') save_wave  ! number of wavepacket to save
                                IF (save_wave > M) THEN
                                        PRINT *, 'Invalid value: save_wave is too large, it has to be < M. Program stopping...'
                                        STOP
                                ENDIF

                                l1 = INT(N * l_1 / (2 * L))   ! semilength of the first barrier in array points unit
                                l2 = INT(N * l_2 / (2 * L))

                        ELSE
                                PRINT *, 'Error in opening "', f_data, '". Program stopping...'
                                STOP
                        ENDIF
                        CLOSE(10)
                        PRINT *, ''
                END SUBROUTINE READ_DATA
END MODULE READ_VAR

