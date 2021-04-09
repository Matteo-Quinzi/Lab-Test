MODULE READ_VAR
        IMPLICIT NONE

        REAL(KIND=8) :: pi = 4.D0 * DATAN(1.D0)
        COMPLEX(KIND=8) :: im = (0.D0, 1.D0)

        INTEGER :: i, save_timestep

        INTEGER :: N, M, save_wave
        REAL(KIND=8) :: L, E_1, l_1, E_2, l_2
        REAL(KIND=8) :: sigma, x0, k0
        REAL(KIND=8) :: dt
        
        CONTAINS
                SUBROUTINE READ_DATA()
                        IMPLICIT NONE

                        OPEN(UNIT=10, FILE='../Input/data_sheet.dat')
                          READ(10, *)
                          READ(10, '(3X, I10)') N
                          READ(10, '(3X, F15.7)') L
                          READ(10, *) 
                          READ(10, *) 
                          READ(10, '(5X, F15.7)') E_1
                          READ(10, '(5X, F15.7)') l_1
                          READ(10, *) 
                          READ(10, '(5X, F15.7)') E_2
                          READ(10, '(5X, F15.7)') l_2
                          READ(10, *) 
                          READ(10, *)
                          READ(10, '(7X, F15.7)') sigma
                          READ(10, '(4X, F15.7)') x0
                          READ(10, '(4X, F15.7)') k0
                          READ(10, *)
                          READ(10, *)
                          READ(10, '(3X, I10)') M
                          READ(10, '(4X, F15.10)') dt
                          READ(10, '(11X, I5)') save_wave
                        CLOSE(10)

                END SUBROUTINE READ_DATA
END MODULE READ_VAR


        
