MODULE READ_VAR
    IMPLICIT NONE

    REAL(KIND=8) :: pi = 4.D0 * DATAN(1.D0)
    COMPLEX(KIND=8) :: im = (0.D0, 1.D0)

   
    INTEGER :: N, d
    REAL(KIND=8) :: a, b, temp 
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
                      IF (a>b) THEN 
                        temp = a 
                        a = b 
                        b = temp
                      ELSE IF (a == b) THEN 
                        PRINT*, "Invalid input: a has to be different from b. Program's stopping."
                        STOP
                      END IF 
                      READ(10, '(4X, I10)') N
                      IF (N<=0) THEN 
                        PRINT*, "Invalid input: N has to be a positive integer number. Program's stopping."
                        STOP 
                      END IF 
                      READ(10, *)
                      READ(10, *) 
                      READ(10, '(6X, A1)') V_str
                      IF (V_str.ne."X".and.V_str.ne."M".and.V_str.ne."H") THEN 
                        PRINT*, "Invalid input: V_str has to be X, M or H. Program's stopping."
                        STOP
                      END IF 
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
                      IF (d<=0) THEN 
                        PRINT*, "Invalid input: d has to be a positive integer number. Program's stopping."
                        STOP 
                      END IF
                      READ(10, *) 
                      READ(10, *)
                      READ(10, '(7X, A1)') task
                      IF (task.ne."a".and.task.ne."b") THEN 
                        PRINT*, "Invalid input: task has to be either a or b. Program's stopping."
                        STOP
                      END IF 
                      READ(10, *) 
                      READ(10, '(7X, A1)') JOBZ
                      IF (JOBZ.ne."V".and.JOBZ.ne."M") THEN 
                        PRINT*, "Invalid input: JOBZ has to be either V or N. Program's stopping."
                        STOP
                      END IF 
                      READ(10, '(8X, A1)') RANGE
                      IF (RANGE.ne."V".and.RANGE.ne."I".and.RANGE.ne."A") THEN 
                        PRINT*, "Invalid input: JOBZ has to be A, V or I. Program's stopping."
                        STOP
                      END IF
                      READ(10, *) 
                      READ(10, '(5X, I10)') IL
                      READ(10, '(5X, I10)') IU
                      IF (RANGE.eq."I") THEN 
                        IF (task.eq."a") THEN 
                          IF (IU-IL>N-1) THEN 
                            PRINT*, "Invalid input: IU-IL has to be smaller than N-1. Program's stopping."
                          STOP 
                          END IF 
                        ELSE IF (task.eq."b") THEN 
                          IF (IU-IL>d-1) THEN 
                            PRINT*, "Invalid input: IU-IL has to be smaller than d-1. Program's stopping."
                          STOP 
                          END IF
                        END IF 
                      END IF 
                      READ(10, *)
                      READ(10, '(5X, F15.7)') VL
                      READ(10, '(5X, F15.7)') VU
                      READ(10, *)
                      READ(10, '(8X, L1)') store
                    CLOSE(10)

            END SUBROUTINE READ_DATA
END MODULE READ_VAR
