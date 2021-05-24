PROGRAM main 

USE read_var    ! read from "data_sheet.dat" and allocate the input paramaters 
USE potential   ! defines the potential function
USE task_a      ! perform the computation in the coordinates space
USE task_b      ! perform the computation in the reciprocal space
USE error       ! conduct a series of convergence tests to identify an efficient set of input parameters.
