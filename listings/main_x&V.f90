ALLOCATE(x(N), V(N))
x(:) = (/ (a + (b-a) * (i*1d0) / (N-1), i = 0, N-1) /)   ! discrete interval [a, b]
IF (V_str == 'M') V = MORSE(x,N)              ! Morse potential V(x)
IF (V_str == 'H') V = HARMONIC(x,N)           ! harmonic potential V(x)
IF (V_str == 'X') V = SHIFTED_HARMONIC(x,N)   ! harmonic potential shifted to overlap with Morse potential