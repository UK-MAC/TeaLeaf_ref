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

!>  @brief Driver for the halo updates
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Invokes the kernels for the internal and external halo cells for
!>  the fields specified.

MODULE update_halo_module

  USE tea_module
  USE report_module
  USE update_halo_kernel_module
  USE update_internal_halo_kernel_module
  USE caliscope_module
 
CONTAINS

SUBROUTINE update_halo(fields,depth)

  IMPLICIT NONE

  INTEGER :: fields(NUM_FIELDS),depth
  REAL(KIND=8) :: timer,halo_time

  TYPE(SCOPE_TYPE):: caliprof

  CALL caliprof%create("update_halo")

  IF (profiler_on) halo_time=timer()
  CALL tea_exchange( fields, depth)
  IF (profiler_on) profiler%halo_exchange = profiler%halo_exchange + (timer() - halo_time)

  CALL update_boundary( fields, depth)

  CALL update_tile_boundary( fields, depth)

END SUBROUTINE update_halo

SUBROUTINE update_boundary( fields,depth)

  IMPLICIT NONE

  INTEGER :: t,fields(NUM_FIELDS),depth
  REAL(KIND=8) :: timer,halo_time

  TYPE(SCOPE_TYPE):: caliprof

  CALL caliprof%create("update_boundary")

  IF (profiler_on) halo_time=timer()

  IF (reflective_boundary .EQV. .TRUE. .AND. ANY(chunk%chunk_neighbours .EQ. EXTERNAL_FACE)) THEN
    IF (use_fortran_kernels)THEN
!$OMP PARALLEL
!$OMP DO
      DO t=1,tiles_per_task
        CALL update_halo_kernel(chunk%tiles(t)%field%x_min,          &
                                chunk%tiles(t)%field%x_max,          &
                                chunk%tiles(t)%field%y_min,          &
                                chunk%tiles(t)%field%y_max,          &
                                chunk%halo_exchange_depth,           &
                                chunk%chunk_neighbours,     &
                                chunk%tiles(t)%tile_neighbours,     &
                                chunk%tiles(t)%field%density,        &
                                chunk%tiles(t)%field%energy0,        &
                                chunk%tiles(t)%field%energy1,        &
                                chunk%tiles(t)%field%u,              &
                                chunk%tiles(t)%field%vector_p,       &
                                chunk%tiles(t)%field%vector_sd,      &
                                chunk%tiles(t)%field%vector_rtemp,      &
                                chunk%tiles(t)%field%vector_z,      &
                                chunk%tiles(t)%field%vector_kx,     &
                                chunk%tiles(t)%field%vector_ky,     &
                                chunk%tiles(t)%field%vector_di,     &
                                fields,                         &
                                depth                           )
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    ENDIF
  ENDIF

  IF (profiler_on) profiler%halo_update = profiler%halo_update + (timer() - halo_time)

END SUBROUTINE update_boundary

SUBROUTINE update_tile_boundary( fields, depth)

  IMPLICIT NONE

  INTEGER :: t,fields(NUM_FIELDS),depth, right_idx, up_idx
  REAL(KIND=8) :: timer,halo_time

  TYPE(SCOPE_TYPE):: caliprof

  CALL caliprof%create("update_tile_boundary")

  IF (profiler_on) halo_time=timer()

  IF (tiles_per_task .GT. 1) THEN
    IF (use_fortran_kernels)THEN
!$OMP PARALLEL PRIVATE(right_idx, up_idx)
!$OMP DO
      DO t=1,tiles_per_task
        right_idx = chunk%tiles(t)%tile_neighbours(CHUNK_RIGHT)

        IF (right_idx .NE. EXTERNAL_FACE) THEN
          CALL update_internal_halo_left_right_kernel(                &
                                  chunk%tiles(t)%field%x_min,          &
                                  chunk%tiles(t)%field%x_max,          &
                                  chunk%tiles(t)%field%y_min,          &
                                  chunk%tiles(t)%field%y_max,          &
                                  chunk%tiles(t)%field%density,        &
                                  chunk%tiles(t)%field%energy0,        &
                                  chunk%tiles(t)%field%energy1,        &
                                  chunk%tiles(t)%field%u,              &
                                  chunk%tiles(t)%field%vector_p,       &
                                  chunk%tiles(t)%field%vector_sd,      &
                                  chunk%tiles(t)%field%vector_rtemp,       &
                                  chunk%tiles(t)%field%vector_z,       &
                                  chunk%tiles(t)%field%vector_kx,      &
                                  chunk%tiles(t)%field%vector_ky,      &
                                  chunk%tiles(t)%field%vector_di,      &
                                  chunk%tiles(right_idx)%field%x_min,          &
                                  chunk%tiles(right_idx)%field%x_max,          &
                                  chunk%tiles(right_idx)%field%y_min,          &
                                  chunk%tiles(right_idx)%field%y_max,          &
                                  chunk%tiles(right_idx)%field%density,        &
                                  chunk%tiles(right_idx)%field%energy0,        &
                                  chunk%tiles(right_idx)%field%energy1,        &
                                  chunk%tiles(right_idx)%field%u,              &
                                  chunk%tiles(right_idx)%field%vector_p,       &
                                  chunk%tiles(right_idx)%field%vector_sd,      &
                                  chunk%tiles(right_idx)%field%vector_rtemp,       &
                                  chunk%tiles(right_idx)%field%vector_z,       &
                                  chunk%tiles(right_idx)%field%vector_kx,      &
                                  chunk%tiles(right_idx)%field%vector_ky,      &
                                  chunk%tiles(right_idx)%field%vector_di,      &
                                  chunk%halo_exchange_depth,                   &
                                  fields,                         &
                                  depth                           )
        ENDIF
      ENDDO
!$OMP END DO NOWAIT

!$  IF (depth .GT. 1) THEN
!$OMP BARRIER
!$  ENDIF

!$OMP DO
      DO t=1,tiles_per_task
        up_idx = chunk%tiles(t)%tile_neighbours(CHUNK_TOP)

        IF (up_idx .NE. EXTERNAL_FACE) THEN
          CALL update_internal_halo_bottom_top_kernel(                &
                                  chunk%tiles(t)%field%x_min,          &
                                  chunk%tiles(t)%field%x_max,          &
                                  chunk%tiles(t)%field%y_min,          &
                                  chunk%tiles(t)%field%y_max,          &
                                  chunk%tiles(t)%field%density,        &
                                  chunk%tiles(t)%field%energy0,        &
                                  chunk%tiles(t)%field%energy1,        &
                                  chunk%tiles(t)%field%u,              &
                                  chunk%tiles(t)%field%vector_p,       &
                                  chunk%tiles(t)%field%vector_sd,      &
                                  chunk%tiles(t)%field%vector_rtemp,      &
                                  chunk%tiles(t)%field%vector_z,       &
                                  chunk%tiles(t)%field%vector_kx,      &
                                  chunk%tiles(t)%field%vector_ky,      &
                                  chunk%tiles(t)%field%vector_di,      &
                                  chunk%tiles(up_idx)%field%x_min,          &
                                  chunk%tiles(up_idx)%field%x_max,          &
                                  chunk%tiles(up_idx)%field%y_min,          &
                                  chunk%tiles(up_idx)%field%y_max,          &
                                  chunk%tiles(up_idx)%field%density,        &
                                  chunk%tiles(up_idx)%field%energy0,        &
                                  chunk%tiles(up_idx)%field%energy1,        &
                                  chunk%tiles(up_idx)%field%u,              &
                                  chunk%tiles(up_idx)%field%vector_p,       &
                                  chunk%tiles(up_idx)%field%vector_sd,      &
                                  chunk%tiles(up_idx)%field%vector_rtemp,      &
                                  chunk%tiles(up_idx)%field%vector_z,       &
                                  chunk%tiles(up_idx)%field%vector_kx,      &
                                  chunk%tiles(up_idx)%field%vector_ky,      &
                                  chunk%tiles(up_idx)%field%vector_di,      &
                                  chunk%halo_exchange_depth,                &
                                  fields,                         &
                                  depth                           )
        ENDIF
      ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
    ENDIF
  ENDIF

  IF (profiler_on) profiler%internal_halo_update = profiler%internal_halo_update + (timer() - halo_time)

END SUBROUTINE update_tile_boundary

END MODULE update_halo_module

