MODULE TASK
        USE READ_VAR 
        CONTAINS
                FUNCTION POT()
                        IMPLICIT NONE
                        
                        ! OUTPUT 
                        REAL(KIND=8) :: POT(N)

                        ! ROUTINE
                        INTEGER(KIND=8) :: l1, l2
                        
                        l1 = INT(N * l_1 / (2 * L))
                        l2 = INT(N * l_2 / (2 * L))

                        POT(:) = 0.D0
                        POT(  N/2 - l1 :   N/2 + l1) = E_1
                        POT(2*N/3 - l2 : 2*N/3 + l2) = E_2
                END FUNCTION POT

                FUNCTION PHI(x) 
                        IMPLICIT NONE

                        ! INPUT
                        REAL(KIND=8), INTENT(IN) :: x(N)
                        
                        ! OUTPUT
                        COMPLEX(KIND=8) :: PHI(N)

                        IF (load_f_inp .eqv. .TRUE.) THEN
                                OPEN(UNIT=10, FILE=f_inp, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
                                IF (IOS == 0) THEN
                                        READ(10, '(2F20.10)') (PHI(i), i = 1,N)

                                ELSE
                                        PRINT *, 'Error in opening "', f_inp, '". Program stopping...'
                                        STOP
                                ENDIF

                        ELSE
                                PHI(:) = (2.D0 * pi * sigma**2)**(-0.25) * EXP(im * k0 * x(:)) * EXP(-(x(:) - x0)**2 / (4.D0 * sigma**2))
                        ENDIF

                END FUNCTION PHI

                FUNCTION FFT(PHI, DIR)
                        IMPLICIT NONE
                        INCLUDE 'fftw3.f'

                        ! INPUT
                        COMPLEX(KIND=8), INTENT(IN) :: PHI(N)
                        CHARACTER(LEN=1), INTENT(IN) :: DIR

                        ! OUTPUT
                        COMPLEX(KIND=8) :: FFT(N)

                        ! ROUTINE
                        INTEGER(KIND=8) :: plan

                        IF (DIR == 'F') THEN   ! forward 
                                CALL DFFTW_PLAN_DFT_1D(plan, N, PHI, FFT, FFTW_FORWARD, FFTW_ESTIMATE)

                        ELSEIF (DIR == 'B') THEN   ! backward
                                CALL DFFTW_PLAN_DFT_1D(plan, N, PHI, FFT, FFTW_BACKWARD, FFTW_ESTIMATE)

                        ENDIF

                        CALL DFFTW_EXECUTE_DFT(plan, PHI, FFT)
                        CALL DFFTW_DESTROY_PLAN(plan)

                END FUNCTION FFT
END MODULE TASK
