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

!>  @brief Mesh chunk generation driver
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Invoked the users specified chunk generator.

SUBROUTINE generate_chunk()

  USE tea_module
  USE generate_chunk_kernel_module

  IMPLICIT NONE

  INTEGER         :: t

  INTEGER         :: state
  REAL(KIND=8), DIMENSION(number_of_states) :: state_density,state_energy
  REAL(KIND=8), DIMENSION(number_of_states) :: state_xmin,state_xmax,state_ymin,state_ymax,state_radius
  INTEGER,      DIMENSION(number_of_states) :: state_geometry

  DO state=1,number_of_states
   state_density(state)=states(state)%density
   state_energy(state)=states(state)%energy
   state_xmin(state)=states(state)%xmin
   state_xmax(state)=states(state)%xmax
   state_ymin(state)=states(state)%ymin
   state_ymax(state)=states(state)%ymax
   state_radius(state)=states(state)%radius
   state_geometry(state)=states(state)%geometry
  ENDDO

  IF(use_fortran_kernels) THEN
!$OMP PARALLEL NUM_THREADS(tiles_per_task)
!$OMP DO
    DO t=1,tiles_per_task
      CALL generate_chunk_kernel(chunk%tiles(t)%field%x_min,             &
                                 chunk%tiles(t)%field%x_max,             &
                                 chunk%tiles(t)%field%y_min,             &
                                 chunk%tiles(t)%field%y_max,halo_exchange_depth,             &
                                 chunk%tiles(t)%field%vertexx,           &
                                 chunk%tiles(t)%field%vertexy,           &
                                 chunk%tiles(t)%field%cellx,             &
                                 chunk%tiles(t)%field%celly,             &
                                 chunk%tiles(t)%field%density,           &
                                 chunk%tiles(t)%field%energy0,           &
                                 chunk%tiles(t)%field%u,                 &
                                 number_of_states,                      &
                                 state_density,                         &
                                 state_energy,                          &
                                 state_xmin,                            &
                                 state_xmax,                            &
                                 state_ymin,                            &
                                 state_ymax,                            &
                                 state_radius,                          &
                                 state_geometry,                        &
                                 g_rect,                                &
                                 g_circ,                                &
                                 g_point)
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE generate_chunk
