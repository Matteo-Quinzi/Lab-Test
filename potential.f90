MODULE POTENTIAL
        USE read_var

        CONTAINS
                FUNCTION HARMONIC(x) RESULT(V)
                        ! Returns the harmonic potential V(x)
                        IMPLICIT NONE

                        ! INPUT
                        REAL(KIND=8), INTENT(IN) :: x(:)   
                        
                        ! OUTPUT
                        REAL(KIND=8), ALLOCATABLE :: V(:)

                        INTEGER :: Npoints
                        !
                        Npoints = size(x)
                        ALLOCATE(V(Npoints))
                        V(:) = 0.5d0 * k_el * x(:)**2
                        !
                END FUNCTION HARMONIC

                FUNCTION MORSE(x) RESULT(V)
                        ! Returns the Morse potential V(x)
                        IMPLICIT NONE

                        ! INPUT 
                        REAL(KIND=8), INTENT(IN) :: x(:)

                        ! OUTPUT
                        REAL(KIND=8), ALLOCATABLE :: V(:)
                        !        
                        INTEGER :: Npoints 
                        !

                        Npoints = size(x)
                        ALLOCATE(V(Npoints))
                        V(:) = D_e * ((1d0 - DEXP(-alpha * (x(:) - x0)))**2 - 1d0)
                        !
                END FUNCTION MORSE
                !
                !
                !
                FUNCTION SHIFTED_HARMONIC(x) RESULT(V)
                  ! Defines an harmonic potential which approximate the 
                  ! Morse potential around the minimum 

                  ! Input 
                  REAL(KIND=8), INTENT(IN) :: x(:)

                  ! Output 
                  REAL(KIND=8), ALLOCATABLE :: V(:)
                  INTEGER :: Npoints 
                  !
                  Npoints = size(x)
                  ALLOCATE(V(Npoints))
                  V(:) = D_e * alpha**2.d0 * (x(:) - x0)**2.d0 - D_e
                  !
                END FUNCTION SHIFTED_HARMONIC
END MODULE
