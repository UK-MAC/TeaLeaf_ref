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

!>  @brief Calculate the minimum timestep for all mesh chunks.
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Invokes the kernels needed to calculate the timestep.

MODULE timestep_module

CONTAINS

SUBROUTINE timestep()

  USE tea_module
  USE report_module
  USE update_halo_module
  USE calc_dt_module
  USE definitions_module

  IMPLICIT NONE

  INTEGER :: t

  REAL(KIND=8)    :: dtlp

  REAL(KIND=8)    :: kernel_time,timer


  IF(profiler_on) kernel_time=timer()

  DO t=1,tiles_per_task
    CALL calc_dt(dtlp)

    IF(dtlp.LE.dt) THEN
      dt=dtlp
    ENDIF
  END DO

  CALL tea_min(dt)

  IF(profiler_on) profiler%timestep=profiler%timestep+(timer()-kernel_time)

  IF (parallel%boss) THEN
      WRITE(g_out,"(' Step ', i7,' time ', f11.7,' timestep  ',1pe9.2,i8)") step,time,dt
      WRITE(0,    "(' Step ', i7,' time ', f11.7,' timestep  ',1pe9.2,i8)") step,time,dt
  ENDIF

END SUBROUTINE timestep

END MODULE timestep_module
