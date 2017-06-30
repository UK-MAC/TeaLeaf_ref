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

!>  @brief Driver for chunk initialisation.
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Invokes the user specified chunk initialisation kernel.

SUBROUTINE initialise_chunk()

  USE definitions_module
  USE tea_module
  USE initialise_chunk_kernel_module

  IMPLICIT NONE

  INTEGER      :: t

  REAL(KIND=8) :: xmin,ymin,dx,dy

  dx=(grid%xmax - grid%xmin)/REAL(grid%x_cells)
  dy=(grid%ymax - grid%ymin)/REAL(grid%y_cells)

  IF(use_fortran_kernels) THEN
!$OMP PARALLEL PRIVATE(xmin, ymin)
!$OMP DO
    DO t=1,tiles_per_task
      xmin=grid%xmin + dx*REAL(chunk%tiles(t)%left-1)

      ymin=grid%ymin + dy*REAL(chunk%tiles(t)%bottom-1)

      CALL initialise_chunk_kernel(chunk%tiles(t)%field%x_min,    &
                                   chunk%tiles(t)%field%x_max,    &
                                   chunk%tiles(t)%field%y_min,    &
                                   chunk%tiles(t)%field%y_max,    &
                                   xmin,ymin,dx,dy,              &
                                   chunk%tiles(t)%field%vertexx,  &
                                   chunk%tiles(t)%field%vertexdx, &
                                   chunk%tiles(t)%field%vertexy,  &
                                   chunk%tiles(t)%field%vertexdy, &
                                   chunk%tiles(t)%field%cellx,    &
                                   chunk%tiles(t)%field%celldx,   &
                                   chunk%tiles(t)%field%celly,    &
                                   chunk%tiles(t)%field%celldy,   &
                                   chunk%tiles(t)%field%volume,   &
                                   chunk%tiles(t)%field%xarea,    &
                                   chunk%tiles(t)%field%yarea     )
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL
  ENDIF

END SUBROUTINE initialise_chunk
