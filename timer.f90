FUNCTION timer()

  USE mpi

  IMPLICIT none

  REAL(KIND=8) :: timer

  timer=mpi_WTIME()

END FUNCTION

