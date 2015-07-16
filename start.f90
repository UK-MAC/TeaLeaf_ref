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

!>  @brief Main set up routine
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Invokes the mesh decomposer and sets up chunk connectivity. It then
!>  allocates the communication buffers and call the chunk initialisation and
!>  generation routines and primes the halo cells and writes an initial field summary.

SUBROUTINE start

  USE tea_module
  USE parse_module
  USE update_halo_module
  USE set_field_module

  IMPLICIT NONE

  INTEGER :: t

  INTEGER :: fields(NUM_FIELDS)

  LOGICAL :: profiler_original

  ! Do no profile the start up costs otherwise the total times will not add up
  ! at the end
  profiler_original=profiler_on
  profiler_on=.FALSE.

  IF (parallel%boss)THEN
    WRITE(g_out,*) 'Setting up initial geometry'
    WRITE(g_out,*)
  ENDIF

  time  = 0.0
  step  = 0
  dt    = dtinit

  CALL tea_barrier()

  CALL tea_decompose(grid%x_cells, grid%y_cells)

  ALLOCATE(chunk%tiles(tiles_per_task))

  chunk%x_cells = chunk%right -chunk%left  +1
  chunk%y_cells = chunk%top   -chunk%bottom+1

  chunk%chunk_x_min = 1
  chunk%chunk_y_min = 1
  chunk%chunk_x_max = chunk%x_cells
  chunk%chunk_y_max = chunk%y_cells

  CALL tea_decompose_tiles(chunk%x_cells, chunk%y_cells)

  DO t=1,tiles_per_task
    chunk%tiles(t)%x_cells = chunk%tiles(t)%right -chunk%tiles(t)%left  +1
    chunk%tiles(t)%y_cells = chunk%tiles(t)%top   -chunk%tiles(t)%bottom+1

    chunk%tiles(t)%field%x_min = 1
    chunk%tiles(t)%field%y_min = 1
    chunk%tiles(t)%field%x_max = chunk%tiles(t)%x_cells
    chunk%tiles(t)%field%y_max = chunk%tiles(t)%y_cells
  ENDDO

  IF (parallel%boss)THEN
    WRITE(g_out,*)"Tile size ",chunk%tiles(1)%x_cells," by ",chunk%tiles(1)%y_cells," cells"
  ENDIF

  CALL build_field()

  CALL tea_allocate_buffers()

  CALL initialise_chunk()

  IF (parallel%boss)THEN
    WRITE(g_out,*) 'Generating chunk'
  ENDIF

  CALL generate_chunk()

  ! Prime all halo data for the first step
  fields=0
  fields(FIELD_DENSITY)=1
  fields(FIELD_ENERGY0)=1
  fields(FIELD_ENERGY1)=1

  CALL update_halo(fields,halo_exchange_depth)

  IF (parallel%boss)THEN
    WRITE(g_out,*)
    WRITE(g_out,*) 'Problem initialised and generated'
  ENDIF

  ! copy time level 0 to time level 1 before the first print
  CALL set_field()

  CALL field_summary()

  IF (visit_frequency.NE.0) CALL visit()

  CALL tea_barrier()

  profiler_on=profiler_original

END SUBROUTINE start

