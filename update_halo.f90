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
  USE update_halo_kernel_module

CONTAINS

SUBROUTINE update_halo(fields,depth)

  IMPLICIT NONE

  INTEGER :: t,fields(NUM_FIELDS),depth
  REAL(KIND=8) :: timer,halo_time

  IF (profiler_on) halo_time=timer()
  CALL tea_exchange(fields,depth)
  IF (profiler_on) profiler%halo_exchange = profiler%halo_exchange + (timer() - halo_time)

  CALL update_boundary(fields, depth)

END SUBROUTINE update_halo

SUBROUTINE update_boundary(fields,depth)

  IMPLICIT NONE

  INTEGER :: t,fields(NUM_FIELDS),depth, right_idx, up_idx
  REAL(KIND=8) :: timer,halo_time

  IF (profiler_on) halo_time=timer()
  IF (reflective_boundary .EQV. .TRUE.) THEN
    IF(use_fortran_kernels)THEN
      DO t=1,tiles_per_task
        CALL update_halo_kernel(chunk%tiles(t)%field%x_min,          &
                                chunk%tiles(t)%field%x_max,          &
                                chunk%tiles(t)%field%y_min,          &
                                chunk%tiles(t)%field%y_max,          &
                                halo_exchange_depth,          &
                                chunk%tiles(t)%tile_neighbours,     &
                                chunk%tiles(t)%field%density,        &
                                chunk%tiles(t)%field%energy0,        &
                                chunk%tiles(t)%field%energy1,        &
                                chunk%tiles(t)%field%u,              &
                                chunk%tiles(t)%field%vector_p,       &
                                chunk%tiles(t)%field%vector_sd,      &
                                fields,                         &
                                depth                           )
      ENDDO
    ENDIF
  ENDIF
  IF (profiler_on) profiler%halo_update = profiler%halo_update + (timer() - halo_time)

  ! TODO profiling info
  IF(use_fortran_kernels)THEN
    DO t=1,tiles_per_task
      IF (chunk%tiles(t)%tile_neighbours(CHUNK_RIGHT) .NE. EXTERNAL_FACE) THEN
        right_idx = chunk%tiles(t)%tile_neighbours(CHUNK_RIGHT)

        IF (chunk%tiles(t)%field%y_cells .NE. chunk%tiles(right_idx)%field%y_cells) THEN
          CALL report_error("update_halo", "Tried to exchange between tow tiles which ahd different sizes")
        ENDIF

        CALL update_internal_halo_left_right_kernel(                &
                                chunk%tiles(t)%field%x_min,          &
                                chunk%tiles(t)%field%x_max,          &
                                chunk%tiles(t)%field%y_min,          &
                                chunk%tiles(t)%field%y_max,          &
                                chunk%tiles(t)%tile_neighbours,     &
                                chunk%tiles(t)%field%density,        &
                                chunk%tiles(t)%field%energy0,        &
                                chunk%tiles(t)%field%energy1,        &
                                chunk%tiles(t)%field%u,              &
                                chunk%tiles(t)%field%vector_p,       &
                                chunk%tiles(t)%field%vector_sd,      &
                                chunk%tiles(right_idx)%field%x_min,          &
                                chunk%tiles(right_idx)%field%x_max,          &
                                chunk%tiles(right_idx)%field%y_min,          &
                                chunk%tiles(right_idx)%field%y_max,          &
                                chunk%tiles(right_idx)%tile_neighbours,     &
                                chunk%tiles(right_idx)%field%density,        &
                                chunk%tiles(right_idx)%field%energy0,        &
                                chunk%tiles(right_idx)%field%energy1,        &
                                chunk%tiles(right_idx)%field%u,              &
                                chunk%tiles(right_idx)%field%vector_p,       &
                                chunk%tiles(right_idx)%field%vector_sd,      &
                                halo_exchange_depth,          &
                                fields,                         &
                                depth                           )
      ENDIF
    ENDDO
  ENDIF

END SUBROUTINE update_boundary

END MODULE update_halo_module
