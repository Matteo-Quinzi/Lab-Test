MODULE POTENTIAL
        IMPLICIT NONE
        INTEGER, PRIVATE:: i
                                                       ! D_e is the well depth
                                                       ! a controls the width of the well
        REAL(KIND=8), SAVE:: D_e, alpha, x_0                                                 ! x0 is the equilibrium bond distance

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

                        OPEN(UNIT=10, FILE='Input/data_sheet.dat', STATUS='old', ACTION='read')
                          DO i = 1,11
                            READ(10, *)
                          END DO

                          READ(10, *) D_e, alpha, x_0
                        CLOSE(10)
                        WRITE(*,*) D_e, alpha
                        V(:) = D_e * ((1d0 - DEXP(-alpha * (x(:) - x_0)))**2 - 1d0)
                END FUNCTION MORSE
                !
                !
                !
                FUNCTION SHIFTED_HARMONIC(x, N) RESULT(V)
                  ! Defines an harmonic potential which approximate the 
                  ! Morse potential around the minimum 

                  ! Input 
                  INTEGER, INTENT(IN) :: N
                  REAL(KIND=8), INTENT(IN) :: x(N)

                  ! Output 
                  REAL(KIND=8) :: V(N)
                  WRITE(*,*) D_e, alpha
                  V(:) = (/ ( D_e * alpha**2.d0 * (x(i) - x_0)**2.d0 - D_e, i = 1,N )   /)
                END FUNCTION SHIFTED_HARMONIC
END MODULE
