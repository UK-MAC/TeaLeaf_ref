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

!>  @brief Fortran mpi buffer packing kernel
!>  @author Wayne Gaudin
!>  @details Packs/unpacks mpi send and receive buffers

MODULE pack_kernel_module

  ! These need to be kept consistent with the data module to avoid use statement
  INTEGER,private,PARAMETER :: CHUNK_LEFT   =1    &
                              ,CHUNK_RIGHT  =2    &
                              ,CHUNK_BOTTOM =3    &
                              ,CHUNK_TOP    =4    &
                              ,EXTERNAL_FACE=-1

  INTEGER,private,PARAMETER   :: FIELD_DENSITY    = 1         &
                                ,FIELD_ENERGY0    = 2         &
                                ,FIELD_ENERGY1    = 3         &
                                ,FIELD_U          = 4         &
                                ,FIELD_P          = 5         &
                                ,FIELD_SD         = 6         &
                                ,FIELD_R          = 7         &
                                ,FIELD_Z          = 8         &
                                ,FIELD_KX         = 9         &
                                ,FIELD_KY         = 10        &
                                ,FIELD_DI         = 11        &
                                ,NUM_FIELDS       = 11

   INTEGER,         PARAMETER :: CELL_DATA     = 1,        &
                                 VERTEX_DATA   = 2,        &
                                 X_FACE_DATA   = 3,        &
                                 y_FACE_DATA   = 4

CONTAINS

FUNCTION yincs(field_type) RESULT(y_inc)

  integer :: field_type, y_inc

  y_inc = 0

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

  x_inc = 0

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

SUBROUTINE pack_all(x_min, x_max, y_min, y_max, halo_exchange_depth, &
    tile_neighbours, &
    density,                                                    &
    energy0,                                                    &
    energy1,                                                    &
    u,                                                          &
    p,                                                          &
    sd,                                                         &
    r,                                                          &
    z,                                                          &
    kx,                                                         &
    ky,                                                         &
    di,                                                         &
    fields, depth, face, packing, mpi_buffer, offsets, tile_offset)

  IMPLICIT NONE

  INTERFACE
    SUBROUTINE pack_or_unpack(x_min,x_max,y_min,y_max,halo_exchange_depth,    &
                              field, mpi_buffer,          &
                              depth, x_inc, y_inc,        &
                              buffer_offset, edge_minus, edge_plus)

      IMPLICIT NONE

      INTEGER      :: depth,x_min,x_max,y_min,y_max,buffer_offset, x_inc, y_inc,halo_exchange_depth, edge_minus, edge_plus
      REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)
      REAL(KIND=8) :: mpi_buffer(:)
    END SUBROUTINE
  END INTERFACE

  INTEGER      :: fields(:)
  INTEGER      :: offsets(:)
  REAL(KIND=8) :: mpi_buffer(:)
  INTEGER      :: face,tile_offset
  LOGICAL      :: packing
  INTEGER      :: depth,x_min,x_max,y_min,y_max, halo_exchange_depth, edge_minus, edge_plus
  INTEGER, DIMENSION(4) :: tile_neighbours

  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                        :: density,energy0,energy1, u, sd, p, r, z, kx, ky, di

  PROCEDURE(pack_or_unpack), POINTER :: pack_func => NULL()

  SELECT CASE (face)
  CASE (CHUNK_LEFT, CHUNK_RIGHT)
    IF (tile_neighbours(CHUNK_BOTTOM) .EQ. EXTERNAL_FACE) THEN
      edge_minus = depth
    ELSE
      edge_minus = 0
    ENDIF

    IF (tile_neighbours(CHUNK_TOP) .EQ. EXTERNAL_FACE) THEN
      edge_plus = depth
    ELSE
      edge_plus = 0
    ENDIF
  CASE (CHUNK_BOTTOM, CHUNK_TOP)
    IF (tile_neighbours(CHUNK_LEFT) .EQ. EXTERNAL_FACE) THEN
      edge_minus = depth
    ELSE
      edge_minus = 0
    ENDIF

    IF (tile_neighbours(CHUNK_RIGHT) .EQ. EXTERNAL_FACE) THEN
      edge_plus = depth
    ELSE
      edge_plus = 0
    ENDIF
  END SELECT

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
    END SELECT
  ENDIF

!$OMP PARALLEL
  IF (fields(FIELD_DENSITY).EQ.1) THEN
      CALL pack_func(x_min,                    &
                     x_max,                    &
                     y_min,                    &
                     y_max,                    &
                     halo_exchange_depth,                    &
                     density,                 &
                     mpi_buffer,                &
                     depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                     tile_offset + offsets(FIELD_DENSITY),   &
                     edge_minus, edge_plus)
  ENDIF
  IF (fields(FIELD_ENERGY0).EQ.1) THEN
      CALL pack_func(x_min,                    &
                     x_max,                    &
                     y_min,                    &
                     y_max,                    &
                     halo_exchange_depth,                    &
                     energy0,                  &
                     mpi_buffer,                &
                     depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                     tile_offset + offsets(FIELD_ENERGY0),   &
                     edge_minus, edge_plus)
  ENDIF
  IF (fields(FIELD_ENERGY1).EQ.1) THEN
      CALL pack_func(x_min,                    &
                     x_max,                    &
                     y_min,                    &
                     y_max,                    &
                     halo_exchange_depth,                    &
                     energy1,                  &
                     mpi_buffer,                &
                     depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                     tile_offset + offsets(FIELD_ENERGY1),   &
                     edge_minus, edge_plus)
  ENDIF
  IF (fields(FIELD_U).EQ.1) THEN
      CALL pack_func(x_min,                    &
                     x_max,                    &
                     y_min,                    &
                     y_max,                    &
                     halo_exchange_depth,                    &
                     u,                  &
                     mpi_buffer,                &
                     depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                     tile_offset + offsets(FIELD_U),   &
                     edge_minus, edge_plus)
  ENDIF
  IF (fields(FIELD_P).EQ.1) THEN
      CALL pack_func(x_min,                    &
                     x_max,                    &
                     y_min,                    &
                     y_max,                    &
                     halo_exchange_depth,                    &
                     p,                  &
                     mpi_buffer,                &
                     depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                     tile_offset + offsets(FIELD_P),   &
                     edge_minus, edge_plus)
  ENDIF
  IF (fields(FIELD_SD).EQ.1) THEN
      CALL pack_func(x_min,                    &
                     x_max,                    &
                     y_min,                    &
                     y_max,                    &
                     halo_exchange_depth,                    &
                     sd,                  &
                     mpi_buffer,                &
                     depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                     tile_offset + offsets(FIELD_SD),   &
                     edge_minus, edge_plus)
  ENDIF
  IF (fields(FIELD_R).EQ.1) THEN
      CALL pack_func(x_min,                    &
                     x_max,                    &
                     y_min,                    &
                     y_max,                    &
                     halo_exchange_depth,                    &
                     r,                  &
                     mpi_buffer,                &
                     depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                     tile_offset + offsets(FIELD_R),   &
                     edge_minus, edge_plus)
  ENDIF
  IF (fields(FIELD_z).EQ.1) THEN
      CALL pack_func(x_min,                    &
                     x_max,                    &
                     y_min,                    &
                     y_max,                    &
                     halo_exchange_depth,                    &
                     z,                  &
                     mpi_buffer,                &
                     depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                     tile_offset + offsets(FIELD_z),   &
                     edge_minus, edge_plus)
  ENDIF
  IF (fields(FIELD_kx).EQ.1) THEN
      CALL pack_func(x_min,                    &
                     x_max,                    &
                     y_min,                    &
                     y_max,                    &
                     halo_exchange_depth,                    &
                     kx,                  &
                     mpi_buffer,                &
                     depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                     tile_offset + offsets(FIELD_kx),   &
                     edge_minus, edge_plus)
                     !depth, xincs(X_FACE_DATA), yincs(X_FACE_DATA),   &
  ENDIF
  IF (fields(FIELD_ky).EQ.1) THEN
      CALL pack_func(x_min,                    &
                     x_max,                    &
                     y_min,                    &
                     y_max,                    &
                     halo_exchange_depth,                    &
                     ky,                  &
                     mpi_buffer,                &
                     depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                     tile_offset + offsets(FIELD_ky),   &
                     edge_minus, edge_plus)
                     !depth, xincs(Y_FACE_DATA), yincs(Y_FACE_DATA),   &
  ENDIF
  IF (fields(FIELD_di).EQ.1) THEN
      CALL pack_func(x_min,                    &
                     x_max,                    &
                     y_min,                    &
                     y_max,                    &
                     halo_exchange_depth,                    &
                     di,                  &
                     mpi_buffer,                &
                     depth, xincs(CELL_DATA), yincs(CELL_DATA),   &
                     tile_offset + offsets(FIELD_di),   &
                     edge_minus, edge_plus)
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE tea_pack_message_left(x_min,x_max,y_min,y_max,halo_exchange_depth,field,                &
                                 left_snd_buffer,                              &
                                 depth,x_inc, y_inc,                             &
                                 buffer_offset, edge_minus, edge_plus)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,halo_exchange_depth
  INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset, edge_minus, edge_plus

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)
  REAL(KIND=8) :: left_snd_buffer(:)

  ! Pack

!$OMP DO
  DO k=y_min-edge_minus,y_max+y_inc+edge_plus
    DO j=1,depth
      index=buffer_offset + j + (k + depth - 1)*depth
      left_snd_buffer(index)=field(x_min+x_inc-1+j,k)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE tea_pack_message_left

SUBROUTINE tea_unpack_message_left(x_min,x_max,y_min,y_max,halo_exchange_depth,field,                &
                                   left_rcv_buffer,                              &
                                   depth,x_inc, y_inc,                             &
                                   buffer_offset, edge_minus, edge_plus)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,halo_exchange_depth
  INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset, edge_minus, edge_plus

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)
  REAL(KIND=8) :: left_rcv_buffer(:)

  ! Unpack

!$OMP DO
  DO k=y_min-edge_minus,y_max+y_inc+edge_plus
    DO j=1,depth
      index= buffer_offset + j+(k+depth-1)*depth
      field(x_min-j,k)=left_rcv_buffer(index)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE tea_unpack_message_left

SUBROUTINE tea_pack_message_right(x_min,x_max,y_min,y_max,halo_exchange_depth,field,                &
                                  right_snd_buffer,                             &
                                   depth,x_inc, y_inc,                             &
                                  buffer_offset, edge_minus, edge_plus)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,halo_exchange_depth
  INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset, edge_minus, edge_plus

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)
  REAL(KIND=8) :: right_snd_buffer(:)

  ! Pack

!$OMP DO
  DO k=y_min-edge_minus,y_max+y_inc+edge_plus
    DO j=1,depth
      index= buffer_offset + j+(k+depth-1)*depth
      right_snd_buffer(index)=field(x_max+1-j,k)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE tea_pack_message_right

SUBROUTINE tea_unpack_message_right(x_min,x_max,y_min,y_max,halo_exchange_depth,field,                &
                                    right_rcv_buffer,                             &
                                   depth,x_inc, y_inc,                             &
                                    buffer_offset, edge_minus, edge_plus)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,halo_exchange_depth
  INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset, edge_minus, edge_plus

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)
  REAL(KIND=8) :: right_rcv_buffer(:)

  ! Unpack

!$OMP DO
  DO k=y_min-edge_minus,y_max+y_inc+edge_plus
    DO j=1,depth
      index= buffer_offset + j+(k+depth-1)*depth
      field(x_max+x_inc+j,k)=right_rcv_buffer(index)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE tea_unpack_message_right

SUBROUTINE tea_pack_message_top(x_min,x_max,y_min,y_max,halo_exchange_depth,field,                &
                                top_snd_buffer,                               &
                                   depth,x_inc, y_inc,                             &
                                buffer_offset, edge_minus, edge_plus)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,halo_exchange_depth
  INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset, edge_minus, edge_plus

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)
  REAL(KIND=8) :: top_snd_buffer(:)

  ! Pack

!$OMP DO
  DO k=1,depth
    DO j=x_min-edge_minus,x_max+x_inc+edge_plus
      index= buffer_offset + j + edge_minus + (k - 1)*(x_max + x_inc + (edge_plus+edge_minus))
      top_snd_buffer(index)=field(j,y_max+1-k)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE tea_pack_message_top

SUBROUTINE tea_unpack_message_top(x_min,x_max,y_min,y_max,halo_exchange_depth,field,                &
                                  top_rcv_buffer,                               &
                                   depth,x_inc, y_inc,                             &
                                  buffer_offset, edge_minus, edge_plus)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,halo_exchange_depth
  INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset, edge_minus, edge_plus

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)
  REAL(KIND=8) :: top_rcv_buffer(:)

  ! Unpack

!$OMP DO
  DO k=1,depth
    DO j=x_min-edge_minus,x_max+x_inc+edge_plus
      index= buffer_offset + j + edge_minus + (k - 1)*(x_max + x_inc + (edge_plus+edge_minus))
      field(j,y_max+y_inc+k)=top_rcv_buffer(index)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE tea_unpack_message_top

SUBROUTINE tea_pack_message_bottom(x_min,x_max,y_min,y_max,halo_exchange_depth,field,                &
                                   bottom_snd_buffer,                            &
                                   depth,x_inc, y_inc,                             &
                                   buffer_offset, edge_minus, edge_plus)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,halo_exchange_depth
  INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset, edge_minus, edge_plus

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)
  REAL(KIND=8) :: bottom_snd_buffer(:)

  ! Pack

!$OMP DO
  DO k=1,depth
    DO j=x_min-edge_minus,x_max+x_inc+edge_plus
      index= buffer_offset + j + edge_minus + (k - 1)*(x_max + x_inc + (edge_plus+edge_minus))
      bottom_snd_buffer(index)=field(j,y_min+y_inc-1+k)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE tea_pack_message_bottom

SUBROUTINE tea_unpack_message_bottom(x_min,x_max,y_min,y_max,halo_exchange_depth,field,                &
                                     bottom_rcv_buffer,                            &
                                   depth,x_inc, y_inc,                             &
                                     buffer_offset, edge_minus, edge_plus)

  IMPLICIT NONE

  INTEGER      :: depth,x_min,x_max,y_min,y_max,halo_exchange_depth
  INTEGER      :: j,k,x_inc,y_inc,index,buffer_offset, edge_minus, edge_plus

  REAL(KIND=8) :: field(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)
  REAL(KIND=8) :: bottom_rcv_buffer(:)

  ! Unpack

!$OMP DO
  DO k=1,depth
    DO j=x_min-edge_minus,x_max+x_inc+edge_plus
      index= buffer_offset + j + edge_minus + (k - 1)*(x_max + x_inc + (edge_plus+edge_minus))
      field(j,y_min-k)=bottom_rcv_buffer(index)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT

END SUBROUTINE tea_unpack_message_bottom

END MODULE pack_kernel_module
