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

!>  @brief set filed driver
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Invokes the user specified set field kernel.

MODULE set_field_module

CONTAINS

SUBROUTINE set_field()

  USE tea_module
  USE set_field_kernel_module

  IMPLICIT NONE

  INTEGER :: t

  REAL(KIND=8) :: kernel_time,timer

  IF(profiler_on) kernel_time=timer()
  DO t=1,tiles_per_task
    IF(use_fortran_kernels)THEN
      CALL set_field_kernel(chunk%tiles(t)%field%x_min,     &
                            chunk%tiles(t)%field%x_max,     &
                            chunk%tiles(t)%field%y_min,     &
                            chunk%tiles(t)%field%y_max,     &
                            halo_exchange_depth,     &
                            chunk%tiles(t)%field%energy0,   &
                            chunk%tiles(t)%field%energy1)
    ENDIF
  ENDDO
  IF(profiler_on) profiler%set_field=profiler%set_field+(timer()-kernel_time)

END SUBROUTINE set_field

END MODULE set_field_module
