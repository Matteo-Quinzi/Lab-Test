MODULE READ_VAR
    IMPLICIT NONE

    REAL(KIND=8) :: pi = 4.D0 * DATAN(1.D0)
    COMPLEX(KIND=8) :: im = (0.D0, 1.D0)

   
    INTEGER :: N, d
    REAL(KIND=8) :: a, b
    REAL(KIND=8) :: k_el
    REAL(KIND=8) :: D_e, alpha, x0
    CHARACTER(1) :: task, V_str, JOBZ, RANGE
    INTEGER :: IL, IU
    REAL(KIND=8) :: VL, VU
    LOGICAL :: store
    INTEGER :: i, j

    
    
    CONTAINS
            SUBROUTINE READ_DATA()
                    IMPLICIT NONE

                    OPEN(UNIT=10, FILE='Input/data_sheet.dat')
                      READ(10, *)
                      READ(10, '(4X, F15.7)') a
                      READ(10, '(4X, F15.7)') b
                      READ(10, '(4X, I10)') N
                      READ(10, *)
                      READ(10, *) 
                      READ(10, '(6X, A1)') V_str
                      READ(10, *)
                      READ(10, *) 
                      READ(10, '(7X, F15.7)') k_el
                      READ(10, *) 
                      READ(10, *) 
                      READ(10, '(6X, F15.7)') D_e
                      READ(10, '(8X, F15.7)') alpha
                      READ(10, '(5X, F15.7)') x0
                      READ(10, *) 
                      READ(10, *)
                      READ(10, '(4X, I10)') d
                      READ(10, *) 
                      READ(10, *)
                      READ(10, '(7X, A1)') task
                      READ(10, *) 
                      READ(10, '(7X, A1)') JOBZ
                      READ(10, '(8X, A1)') RANGE
                      READ(10, *) 
                      READ(10, '(5X, I10)') IL
                      READ(10, '(5X, I10)') IU
                      READ(10, *)
                      READ(10, '(5X, F15.7)') VL
                      READ(10, '(5X, F15.7)') VU
                      READ(10, *)
                      READ(10, '(8X, L1)') store
                    CLOSE(10)

            END SUBROUTINE READ_DATA
END MODULE READ_VAR
