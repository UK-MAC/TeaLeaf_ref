
MODULE pack_module

  USE definitions_module
  USE pack_kernel_module
  USE report_module

CONTAINS

SUBROUTINE tea_pack_buffers(fields, depth, face, mpi_buffer, offsets)

  IMPLICIT NONE

  INTEGER      :: fields(:),depth
  INTEGER      :: offsets(:)
  REAL(KIND=8) :: mpi_buffer(:)
  INTEGER       :: face

  CALL call_packing_functions(fields, depth, face, .TRUE., mpi_buffer, offsets)

END SUBROUTINE

SUBROUTINE tea_unpack_buffers(fields, depth, face, mpi_buffer, offsets)

  IMPLICIT NONE

  INTEGER      :: fields(:),depth
  INTEGER      :: offsets(:)
  REAL(KIND=8) :: mpi_buffer(:)
  INTEGER       :: face

  CALL call_packing_functions(fields, depth, face, .FALSE., mpi_buffer, offsets)

END SUBROUTINE

SUBROUTINE call_packing_functions(fields, depth, face, packing, mpi_buffer, offsets)

  IMPLICIT NONE

  INTEGER      :: fields(:),depth
  INTEGER      :: offsets(:)
  REAL(KIND=8) :: mpi_buffer(:)
  INTEGER      :: face,t,tile_offset
  LOGICAL      :: packing

  DO t=1,tiles_per_task
    SELECT CASE (face)
    CASE (CHUNK_LEFT, CHUNK_RIGHT)
      tile_offset = (chunk%tiles(t)%bottom - chunk%bottom)*depth
    CASE (CHUNK_BOTTOM, CHUNK_TOP)
      tile_offset = (chunk%tiles(t)%left - chunk%left)*depth
    CASE DEFAULT
      CALL report_error("pack.f90","Invalid face pased to buffer packing")
    END SELECT

    IF (chunk%tiles(t)%tile_neighbours(face) .NE. EXTERNAL_FACE) THEN
      CYCLE
    ENDIF

    CALL pack_all(chunk%tiles(t)%field%x_min,                    &
                  chunk%tiles(t)%field%x_max,                    &
                  chunk%tiles(t)%field%y_min,                    &
                  chunk%tiles(t)%field%y_max,                    &
                  halo_exchange_depth,                    &
                  chunk%tiles(t)%tile_neighbours,     &
                  chunk%tiles(t)%field%density,        &
                  chunk%tiles(t)%field%energy0,        &
                  chunk%tiles(t)%field%energy1,        &
                  chunk%tiles(t)%field%u,              &
                  chunk%tiles(t)%field%vector_p,       &
                  chunk%tiles(t)%field%vector_sd,      &
                  chunk%tiles(t)%field%vector_r,      &
                  chunk%tiles(t)%field%vector_rtemp,      &                  
                  fields, &
                  depth, &
                  face, &
                  packing, &
                  mpi_buffer,                &
                  offsets, &
                  tile_offset)
  ENDDO

END SUBROUTINE

END MODULE

