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

!>  @brief Driver for the timestep kernels
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Invokes the user specified timestep kernel.

MODULE calc_dt_module

CONTAINS

SUBROUTINE calc_dt(chunk,local_dt)

  USE tea_module

  IMPLICIT NONE

  INTEGER          :: chunk
  REAL(KIND=8)     :: local_dt

  local_dt = dtinit

END SUBROUTINE calc_dt

END MODULE calc_dt_module
