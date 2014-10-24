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

!>  @brief Fortran timestep kernel
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Returns the fixed timestep

MODULE calc_dt_kernel_module

CONTAINS

SUBROUTINE calc_dt_kernel(dt)

  IMPLICIT NONE

  REAL(KIND=8) :: dt

END SUBROUTINE calc_dt_kernel

END MODULE calc_dt_kernel_module

