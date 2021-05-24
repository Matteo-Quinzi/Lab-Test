REAL(KIND=8), ALLOCATABLE :: x(:)   ! Interval [a,b] equapartition
REAL(KIND=8), ALLOCATABLE  :: V(:)   ! Potential V(x)
REAL(KIND=8) :: toll ! Tolerance of the convergence tests
INTEGER :: maxInd, maxPotPerc  !Error analysis parameters