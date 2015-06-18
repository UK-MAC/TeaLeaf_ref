
MODULE pack_module

  USE definitions_module
  USE pack_kernel_module

CONTAINS

FUNCTION yincs(field_type) RESULT(y_inc)
  integer :: field_type, y_inc

  IF (field_type.EQ.CELL_DATA) THEN
    y_inc=0
  ELSEIF (field_type.EQ.VERTEX_DATA) THEN
    y_inc=1
  ELSEIF (field_type.EQ.X_FACE_DATA) THEN
    y_inc=0
  ELSEIF (field_type.EQ.Y_FACE_DATA) THEN
    y_inc=1
  ENDIF
END FUNCTION

FUNCTION xincs(field_type) RESULT(x_inc)
  integer :: field_type, x_inc

  IF (field_type.EQ.CELL_DATA) THEN
    x_inc=0
  ELSEIF (field_type.EQ.VERTEX_DATA) THEN
    x_inc=1
  ELSEIF (field_type.EQ.X_FACE_DATA) THEN
    x_inc=1
  ELSEIF (field_type.EQ.Y_FACE_DATA) THEN
    x_inc=0
  ENDIF
END FUNCTION

SUBROUTINE tea_pack_buffers(fields, depth, face, mpi_buffer, offsets)
  IMPLICIT NONE

  INTEGER      :: fields(:),depth
  INTEGER      :: offsets(:)
  REAL(KIND=8) :: mpi_buffer(:)
  INTEGER       :: face
  LOGICAL       :: packing

  CALL call_packing_functions(fields, depth, face, .TRUE., mpi_buffer, offsets)

END SUBROUTINE

SUBROUTINE tea_unpack_buffers(fields, depth, face, mpi_buffer, offsets)
  IMPLICIT NONE

  INTEGER      :: fields(:),depth
  INTEGER      :: offsets(:)
  REAL(KIND=8) :: mpi_buffer(:)
  INTEGER       :: face
  LOGICAL       :: packing

  CALL call_packing_functions(fields, depth, face, .FALSE., mpi_buffer, offsets)

END SUBROUTINE

SUBROUTINE call_packing_functions(fields, depth, face, packing, mpi_buffer, offsets)

  IMPLICIT NONE

  INTERFACE
    SUBROUTINE pack_or_unpack(x_min,x_max,y_min,y_max,halo_exchange_depth,    &
                              field, mpi_buffer,          &
                              depth, x_inc, y_inc,        &
                              buffer_offset)

      IMPLICIT NONE

      INTEGER      :: depth,x_min,x_max,y_min,y_max,buffer_offset, x_inc, y_inc,halo_exchange_depth
      REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)
      REAL(KIND=8) :: mpi_buffer(:)
    END SUBROUTINE
  END INTERFACE

  INTEGER      :: fields(:),depth
  INTEGER      :: offsets(:)
  REAL(KIND=8) :: mpi_buffer(:)
  INTEGER      :: face,t
  LOGICAL      :: packing

  PROCEDURE(pack_or_unpack), POINTER :: pack_func => NULL()

  IF (packing .EQV. .TRUE.) THEN
    SELECT CASE (face)
    CASE (CHUNK_LEFT)
      pack_func => tea_pack_message_left
    CASE (CHUNK_RIGHT)
      pack_func => tea_pack_message_right
    CASE (CHUNK_BOTTOM)
      pack_func => tea_pack_message_bottom
    CASE (CHUNK_TOP)
      pack_func => tea_pack_message_top
    CASE DEFAULT
      !call report_error("pack.f90","Invalid face pased to buffer packing")
    END SELECT
  ELSE
    SELECT CASE (face)
    CASE (CHUNK_LEFT)
      pack_func => tea_unpack_message_left
    CASE (CHUNK_RIGHT)
      pack_func => tea_unpack_message_right
    CASE (CHUNK_BOTTOM)
      pack_func => tea_unpack_message_bottom
    CASE (CHUNK_TOP)
      pack_func => tea_unpack_message_top
    CASE DEFAULT
      !call report_error("pack.f90","Invalid face pased to buffer packing")
    END SELECT
  ENDIF

  ! FIXME make these pack correctly

!$OMP PARALLEL
!$OMP DO
  DO t=1,tiles_per_task
    IF (fields(FIELD_DENSITY).EQ.1) THEN
        CALL pack_func(chunk%tiles(t)%field%x_min,                    &
                       chunk%tiles(t)%field%x_max,                    &
                       chunk%tiles(t)%field%y_min,                    &
                       chunk%tiles(t)%field%y_max,                    &
                       halo_exchange_depth,                    &
                       chunk%tiles(t)%field%density,                 &
                       mpi_buffer,                &
                       depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                       offsets(FIELD_DENSITY))
    ENDIF
    IF (fields(FIELD_ENERGY0).EQ.1) THEN
        CALL pack_func(chunk%tiles(t)%field%x_min,                    &
                       chunk%tiles(t)%field%x_max,                    &
                       chunk%tiles(t)%field%y_min,                    &
                       chunk%tiles(t)%field%y_max,                    &
                       halo_exchange_depth,                    &
                       chunk%tiles(t)%field%energy0,                  &
                       mpi_buffer,                &
                       depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                       offsets(FIELD_ENERGY0))
    ENDIF
    IF (fields(FIELD_ENERGY1).EQ.1) THEN
        CALL pack_func(chunk%tiles(t)%field%x_min,                    &
                       chunk%tiles(t)%field%x_max,                    &
                       chunk%tiles(t)%field%y_min,                    &
                       chunk%tiles(t)%field%y_max,                    &
                       halo_exchange_depth,                    &
                       chunk%tiles(t)%field%energy1,                  &
                       mpi_buffer,                &
                       depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                       offsets(FIELD_ENERGY1))
    ENDIF
    IF (fields(FIELD_P).EQ.1) THEN
        CALL pack_func(chunk%tiles(t)%field%x_min,                    &
                       chunk%tiles(t)%field%x_max,                    &
                       chunk%tiles(t)%field%y_min,                    &
                       chunk%tiles(t)%field%y_max,                    &
                       halo_exchange_depth,                    &
                       chunk%tiles(t)%field%vector_p,                  &
                       mpi_buffer,                &
                       depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                       offsets(FIELD_P))
    ENDIF
    IF (fields(FIELD_U).EQ.1) THEN
        CALL pack_func(chunk%tiles(t)%field%x_min,                    &
                       chunk%tiles(t)%field%x_max,                    &
                       chunk%tiles(t)%field%y_min,                    &
                       chunk%tiles(t)%field%y_max,                    &
                       halo_exchange_depth,                    &
                       chunk%tiles(t)%field%u,                  &
                       mpi_buffer,                &
                       depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                       offsets(FIELD_U))
    ENDIF
    IF (fields(FIELD_SD).EQ.1) THEN
        CALL pack_func(chunk%tiles(t)%field%x_min,                    &
                       chunk%tiles(t)%field%x_max,                    &
                       chunk%tiles(t)%field%y_min,                    &
                       chunk%tiles(t)%field%y_max,                    &
                       halo_exchange_depth,                    &
                       chunk%tiles(t)%field%vector_sd,                  &
                       mpi_buffer,                &
                       depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                       offsets(FIELD_SD))
    ENDIF
    IF (fields(FIELD_R).EQ.1) THEN
        CALL pack_func(chunk%tiles(t)%field%x_min,                    &
                       chunk%tiles(t)%field%x_max,                    &
                       chunk%tiles(t)%field%y_min,                    &
                       chunk%tiles(t)%field%y_max,                    &
                       halo_exchange_depth,                    &
                       chunk%tiles(t)%field%vector_r,                  &
                       mpi_buffer,                &
                       depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                       offsets(FIELD_R))
    ENDIF
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

END MODULE

