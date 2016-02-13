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

  INTEGER :: level,t

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

  DO level=1,2
  !write(6,*) "level=",level," tiles_per_task=",tiles_per_task
  CALL tea_decompose(level, grid(level)%x_cells, grid(level)%y_cells)

  ALLOCATE(chunk(level)%tiles(tiles_per_task))

  chunk(level)%x_cells = chunk(level)%right -chunk(level)%left  +1
  chunk(level)%y_cells = chunk(level)%top   -chunk(level)%bottom+1

  chunk(level)%chunk_x_min = 1
  chunk(level)%chunk_y_min = 1
  chunk(level)%chunk_x_max = chunk(level)%x_cells
  chunk(level)%chunk_y_max = chunk(level)%y_cells

  !write(6,*) "calling tea_decompose_tiles"
  CALL tea_decompose_tiles(level, chunk(level)%x_cells, chunk(level)%y_cells)
  !write(6,*) "after   tea_decompose_tiles"

  DO t=1,tiles_per_task
    chunk(level)%tiles(t)%x_cells = chunk(level)%tiles(t)%right -chunk(level)%tiles(t)%left  +1
    chunk(level)%tiles(t)%y_cells = chunk(level)%tiles(t)%top   -chunk(level)%tiles(t)%bottom+1

    chunk(level)%tiles(t)%field%x_min = 1
    chunk(level)%tiles(t)%field%y_min = 1
    chunk(level)%tiles(t)%field%x_max = chunk(level)%tiles(t)%x_cells
    chunk(level)%tiles(t)%field%y_max = chunk(level)%tiles(t)%y_cells
  ENDDO

  IF (parallel%boss)THEN
    WRITE(g_out,*)"Tile size ",chunk(level)%tiles(1)%x_cells," by ",chunk(level)%tiles(1)%y_cells," cells"
  ENDIF

  IF (parallel%boss .AND. level < 2)THEN
    WRITE(g_out,*)"Sub-tile size ranges from ",floor  (chunk(level)%tiles(tiles_per_task)%x_cells/ &
                                               real(chunk(level)%sub_tile_dims(1)))," by ", &
                                               floor  (chunk(level)%tiles(tiles_per_task)%y_cells/ &
                                               real(chunk(level)%sub_tile_dims(2)))," cells", &
                                        " to ",ceiling(chunk(level)%tiles(1)             %x_cells/ &
                                               real(chunk(level)%sub_tile_dims(1)))," by ", &
                                               ceiling(chunk(level)%tiles(1)             %y_cells/ &
                                               real(chunk(level)%sub_tile_dims(2)))," cells"
  ENDIF

  CALL build_field(level)

  CALL tea_allocate_buffers(level)

  CALL initialise_chunk(level)

  IF (parallel%boss)THEN
    WRITE(g_out,*) 'Generating chunk on level ',level
  ENDIF

  IF (level < 2) THEN
    grid(level+1)%x_cells=mpi_dims(1)*chunk(level)%tile_dims(1)*chunk(level)%sub_tile_dims(1)
    grid(level+1)%y_cells=mpi_dims(2)*chunk(level)%tile_dims(2)*chunk(level)%sub_tile_dims(2)
  ENDIF
  ENDDO

  level=1
  CALL generate_chunk(level)

  ! Prime all halo data for the first step
  fields=0
  fields(FIELD_DENSITY)=1
  fields(FIELD_ENERGY0)=1
  fields(FIELD_ENERGY1)=1

  CALL update_halo(level,fields,chunk(level)%halo_exchange_depth)

  IF (parallel%boss)THEN
    WRITE(g_out,*)
    WRITE(g_out,*) 'Problem initialised and generated'
  ENDIF

  ! copy time level 0 to time level 1 before the first print
  CALL set_field(level)

  CALL field_summary()

  IF (visit_frequency.NE.0) CALL visit()

  CALL tea_barrier()

  profiler_on=profiler_original

END SUBROUTINE start

