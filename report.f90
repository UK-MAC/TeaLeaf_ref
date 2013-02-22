MODULE report_module

  USE data_module
  USE clover_module
 
CONTAINS

SUBROUTINE report_error(location, error)

  IMPLICIT NONE

  CHARACTER(LEN=*)  :: location, error

  WRITE(*    ,*)
  WRITE(*    ,*)  'Error from ',location,':'
  WRITE(*    ,*)  error
  WRITE(g_out,*)
  WRITE(g_out,*)  'Error from ',location,':'
  WRITE(g_out,*)  error
  WRITE(0    ,*)
  WRITE(0    ,*)  'Error from ',location,':'
  WRITE(0    ,*)  error
  WRITE(*    ,*)
  WRITE(g_out,*)
  WRITE(0    ,*)
  WRITE(*    ,*) 'CLOVER is terminating.'
  WRITE(*,    *)
  WRITE(g_out,*) 'CLOVER is terminating.'
  WRITE(g_out,*)
  WRITE(0    ,*) 'CLOVER is terminating.'
  WRITE(0    ,*)

  CALL clover_abort

END SUBROUTINE report_error

END MODULE report_module
