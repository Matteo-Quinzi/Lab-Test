MODULE POTENTIAL
        IMPLICIT NONE
        INTEGER :: i

        CONTAINS
                FUNCTION HARMONIC(x, N) RESULT(V)
                        ! Returns the harmonic potential V(x)

                        IMPLICIT NONE

                        ! INPUT
                        INTEGER, INTENT(IN) :: N
                        REAL(KIND=8), INTENT(IN) :: x(N)   
                        
                        ! OUTPUT
                        REAL(KIND=8) :: V(N)

                        ! ROUTINE
                        REAL(KIND=8) :: k   ! spring constant

                        OPEN(UNIT=10, FILE='Input/data_sheet.dat', STATUS='old', ACTION='read')
                          DO i = 1,7
                            READ(10, *)
                          END DO

                          READ(10, *) k
                        CLOSE(10)

                        V(:) = 0.5d0 * k * x(:)**2
                END FUNCTION HARMONIC

                FUNCTION MORSE(x, N) RESULT(V)
                        ! Returns the Morse potential V(x)

                        IMPLICIT NONE

                        ! INPUT 
                        INTEGER, INTENT(IN) :: N
                        REAL(KIND=8), INTENT(IN) :: x(N)

                        ! OUTPUT
                        REAL(KIND=8) :: V(N)

                        ! ROUTINE
                        REAL(KIND=8) :: D_e, a, x0   ! D_e is the well depth
                                                     ! a controls the width of the well
                                                     ! x0 is the equilibrium bond distance

                        OPEN(UNIT=10, FILE='Input/data_sheet.dat', STATUS='old', ACTION='read')
                          DO i = 1,11
                            READ(10, *)
                          END DO

                          READ(10, *) D_e, a, x0
                        CLOSE(10)

                        V(:) = D_e * ((1d0 - DEXP(-a * (x(:) - x0)))**2 - 1d0)
                END FUNCTION MORSE
END MODULE
