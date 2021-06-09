MODULE TASK
        ! Defines the routines called by the main program to perform the computation.
        USE READ_VAR   ! Gains access to the variables declared in READ_VAR module.
        CONTAINS
                FUNCTION POT()
                        ! Builds the potential energy array.
                        IMPLICIT NONE
                        
                        ! OUTPUT
                        REAL(KIND=8) :: POT(N)   ! potential energy array

                        ! The potential energy is everywhere zero except for two regions where it has a constant values of E_1 and E_2.
                        POT(:) = 0.D0   
                        POT(  N/2 - l1 :   N/2 + l1) = E_1   ! first barrier
                        POT(2*N/3 - l2 : 2*N/3 + l2) = E_2   ! second barrier
                END FUNCTION POT

                FUNCTION PHI(x) 
                        ! Returns the initial wavepacket array, computed over the points of the spatial interval array x. 
                        IMPLICIT NONE

                        ! INPUT
                        REAL(KIND=8), INTENT(IN) :: x(N)   ! spatial interval array 
                        
                        ! OUTPUT
                        COMPLEX(KIND=8) :: PHI(N)   ! wavepacket array

                        ! If the boolean parameter load_f_inp (read from the external file) is .TRUE. the wavepacket is read from 
                        ! the external file f_inp. If .FALSE., it is initialized as a Gaussian wavepacket.
                        IF (load_f_inp .eqv. .TRUE.) THEN
                                OPEN(UNIT=10, FILE=f_inp, STATUS='OLD', ACTION='READ', IOSTAT=IOS)
                                IF (IOS == 0) THEN
                                        READ(10, '(2F20.10)') (PHI(i), i = 1,N)
                                ELSE
                                        PRINT *, 'Error in opening "', f_inp, '". Program stopping...'
                                        STOP
                                ENDIF

                        ELSE  
                                PHI(:) = (2.D0 * pi * sigma**2)**(-0.25) * EXP(im * k0 * x(:)) * EXP(-(x(:) - x0)**2 / (4.D0 * sigma**2))   ! Gaussian wavepacket
                        ENDIF
                END FUNCTION PHI

                FUNCTION FFT(PHI, DIR)
                        ! Computes the discrete Fourier transform of an input array PHI. The computation is performed with the help
                        ! of the external subroutine library FFTW (http://www.fftw.org).
                        IMPLICIT NONE
                        INCLUDE 'fftw3.f'  ! FFTW interface

                        ! INPUT
                        COMPLEX(KIND=8), INTENT(IN) :: PHI(N)  ! input array to transform
                        CHARACTER(LEN=1), INTENT(IN) :: DIR   ! transform direction

                        ! OUTPUT
                        COMPLEX(KIND=8) :: FFT(N)   ! tranformed array

                        ! ROUTINE
                        INTEGER(KIND=8) :: plan   ! integer storing the information about the transform

                        IF (DIR == 'F') THEN  
                                CALL DFFTW_PLAN_DFT_1D(plan, N, PHI, FFT, FFTW_FORWARD, FFTW_ESTIMATE)   ! forward transform

                        ELSEIF (DIR == 'B') THEN   
                                CALL DFFTW_PLAN_DFT_1D(plan, N, PHI, FFT, FFTW_BACKWARD, FFTW_ESTIMATE)   ! backward transform

                        ENDIF

                        CALL DFFTW_EXECUTE_DFT(plan, PHI, FFT)   ! subroutine call to execute the calculus
                        CALL DFFTW_DESTROY_PLAN(plan)
                END FUNCTION FFT
END MODULE TASK
