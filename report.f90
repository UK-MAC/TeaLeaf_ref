!Crown Copyright 2014 AWE.
!
! This file is part of TeaLeaf.
!
! TeaLeaf is free software: you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the 
! Free Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! TeaLeaf is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details.
!
! You should have received a copy of the GNU General Public License along with 
! TeaLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Controls error reporting 
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Outputs error messages and aborts the calculation.

MODULE report_module

  USE data_module
  USE tea_module
 
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
  WRITE(*    ,*) 'TEA is terminating.'
  WRITE(*,    *)
  WRITE(g_out,*) 'TEA is terminating.'
  WRITE(g_out,*)
  WRITE(0    ,*) 'TEA is terminating.'
  WRITE(0    ,*)

  CALL tea_abort

END SUBROUTINE report_error

END MODULE report_module
