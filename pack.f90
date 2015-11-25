
MODULE pack_module

  USE definitions_module
  USE pack_kernel_module
  USE report_module

CONTAINS

SUBROUTINE tea_pack_buffers(level, fields, depth, face, mpi_buffer, offsets)

  IMPLICIT NONE

  INTEGER      :: level,fields(:),depth
  INTEGER      :: offsets(:)
  REAL(KIND=8) :: mpi_buffer(:)
  INTEGER       :: face

  CALL call_packing_functions(level, fields, depth, face, .TRUE., mpi_buffer, offsets)

END SUBROUTINE

SUBROUTINE tea_unpack_buffers(level, fields, depth, face, mpi_buffer, offsets)

  IMPLICIT NONE

  INTEGER      :: level,fields(:),depth
  INTEGER      :: offsets(:)
  REAL(KIND=8) :: mpi_buffer(:)
  INTEGER       :: face

  CALL call_packing_functions(level, fields, depth, face, .FALSE., mpi_buffer, offsets)

END SUBROUTINE

SUBROUTINE call_packing_functions(level, fields, depth, face, packing, mpi_buffer, offsets)

  IMPLICIT NONE

  INTEGER      :: level,fields(:),depth
  INTEGER      :: offsets(:)
  REAL(KIND=8) :: mpi_buffer(:)
  INTEGER      :: face,t,tile_offset
  LOGICAL      :: packing

!$OMP PARALLEL PRIVATE(tile_offset)
!$OMP DO
  DO t=1,tiles_per_task
    SELECT CASE (face)
    CASE (CHUNK_LEFT, CHUNK_RIGHT)
      tile_offset = (chunk(level)%tiles(t)%bottom - chunk(level)%bottom)*depth
    CASE (CHUNK_BOTTOM, CHUNK_TOP)
      tile_offset = (chunk(level)%tiles(t)%left - chunk(level)%left)*depth
      IF (tile_offset .NE. 0) THEN
        tile_offset = tile_offset + depth*depth
      ENDIF
    CASE DEFAULT
      CALL report_error("pack.f90","Invalid face pased to buffer packing")
    END SELECT

    IF (chunk(level)%tiles(t)%tile_neighbours(face) .NE. EXTERNAL_FACE) THEN
      CYCLE
    ENDIF

    CALL pack_all(chunk(level)%tiles(t)%field%x_min,                    &
                  chunk(level)%tiles(t)%field%x_max,                    &
                  chunk(level)%tiles(t)%field%y_min,                    &
                  chunk(level)%tiles(t)%field%y_max,                    &
                  halo_exchange_depth,                    &
                  chunk(level)%tiles(t)%tile_neighbours,     &
                  chunk(level)%tiles(t)%field%density,        &
                  chunk(level)%tiles(t)%field%energy0,        &
                  chunk(level)%tiles(t)%field%energy1,        &
                  chunk(level)%tiles(t)%field%u,              &
                  chunk(level)%tiles(t)%field%vector_p,       &
                  chunk(level)%tiles(t)%field%vector_sd,      &
                  chunk(level)%tiles(t)%field%vector_r,      &
                  chunk(level)%tiles(t)%field%vector_z,      &
                  chunk(level)%tiles(t)%field%vector_kx,     &
                  chunk(level)%tiles(t)%field%vector_ky,     &
                  chunk(level)%tiles(t)%field%vector_di,     &
                  fields, &
                  depth, &
                  face, &
                  packing, &
                  mpi_buffer,                &
                  offsets, &
                  tile_offset)
  ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE

END MODULE

