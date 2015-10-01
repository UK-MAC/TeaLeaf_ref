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

!>  @brief Communication Utilities
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Contains all utilities required to run TeaLeaf in a distributed
!>  environment, including initialisation, mesh decompostion, reductions and
!>  halo exchange using explicit buffers.
!>
!>  Note the halo exchange is currently coded as simply as possible and no
!>  optimisations have been implemented, such as post receives before sends or packing
!>  buffers with multiple data fields. This is intentional so the effect of these
!>  optimisations can be measured on large systems, as and when they are added.
!>
!>  Even without these modifications TeaLeaf weak scales well on moderately sized
!>  systems of the order of 10K cores.

MODULE tea_module

  USE definitions_module
  USE pack_module
  USE global_mpi_module
  USE report_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_init_comms

  IMPLICIT NONE

  INTEGER :: err,rank,size
  INTEGER, dimension(2)  :: periodic
  ! not periodic
  DATA periodic/0, 0/

  mpi_dims = 0

  rank=0
  size=1

  CALL MPI_INIT(err)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,err)

  ! Create comm and get coords
  CALL MPI_DIMS_CREATE(size, 2, mpi_dims, err)
  CALL MPI_CART_CREATE(MPI_COMM_WORLD, 2, mpi_dims, periodic, 1, mpi_cart_comm, err)

  CALL MPI_COMM_RANK(mpi_cart_comm,rank,err)
  CALL MPI_COMM_SIZE(mpi_cart_comm,size,err)
  CALL MPI_CART_COORDS(mpi_cart_comm, rank, 2, mpi_coords, err)

  IF (rank.EQ.0) THEN
    parallel%boss=.TRUE.
  ENDIF

  parallel%task = rank
  parallel%boss_task=0
  parallel%max_task=size

END SUBROUTINE tea_init_comms

SUBROUTINE tea_finalize

  INTEGER :: err

  CLOSE(g_out)
  CALL FLUSH(0)
  CALL FLUSH(6)
  CALL FLUSH(g_out)
  CALL MPI_FINALIZE(err)

END SUBROUTINE tea_finalize

SUBROUTINE tea_decompose(x_cells,y_cells)

  ! This decomposes the mesh into a number of chunks.

  IMPLICIT NONE

  INTEGER :: x_cells,y_cells
  INTEGER :: delta_x,delta_y

  INTEGER  :: chunk_x,chunk_y,mod_x,mod_y

  INTEGER  :: err

  ! Get destinations/sources
  CALL mpi_cart_shift(mpi_cart_comm, 1, 1,      &
    chunk%chunk_neighbours(CHUNK_BOTTOM),   &
    chunk%chunk_neighbours(CHUNK_TOP),      &
    err)
  CALL mpi_cart_shift(mpi_cart_comm, 0, 1,      &
    chunk%chunk_neighbours(CHUNK_LEFT),     &
    chunk%chunk_neighbours(CHUNK_RIGHT),    &
    err)

  WHERE (chunk%chunk_neighbours .EQ. MPI_PROC_NULL)
    chunk%chunk_neighbours = EXTERNAL_FACE
  END WHERE

  chunk_x = mpi_dims(1)
  chunk_y = mpi_dims(2)

  delta_x=x_cells/chunk_x
  delta_y=y_cells/chunk_y
  mod_x=MOD(x_cells,chunk_x)
  mod_y=MOD(y_cells,chunk_y)

  chunk%left = mpi_coords(1)*delta_x + 1
  if (mpi_coords(1) .le. mod_x) then
    chunk%left = chunk%left + mpi_coords(1)
  else
    chunk%left = chunk%left + mod_x
  endif
  chunk%right = chunk%left+delta_x - 1
  if (mpi_coords(1) .lt. mod_x) then
    chunk%right = chunk%right + 1
  endif

  chunk%bottom = mpi_coords(2)*delta_y + 1
  if (mpi_coords(2) .le. mod_y) then
    chunk%bottom = chunk%bottom + mpi_coords(2)
  else
    chunk%bottom = chunk%bottom + mod_y
  endif
  chunk%top = chunk%bottom+delta_y - 1
  if (mpi_coords(2) .lt. mod_y) then
    chunk%top = chunk%top + 1
  endif

  IF (parallel%boss)THEN
    WRITE(g_out,*)
    WRITE(g_out,*)"Mesh ratio of ",REAL(x_cells)/REAL(y_cells)
    WRITE(g_out,*)"Decomposing the mesh into ",chunk_x," by ",chunk_y," chunks"
    WRITE(g_out,*)
  ENDIF

END SUBROUTINE tea_decompose

SUBROUTINE tea_decompose_tiles(x_cells, y_cells)

  IMPLICIT NONE

  INTEGER :: x_cells, y_cells

  INTEGER :: delta_x,delta_y
  INTEGER  :: tiles_x,tiles_y,mod_x,mod_y

  ! create a fake communicator to easily decompose the tiles as before
  INTEGER :: err, j, k, t
  INTEGER, DIMENSION(2)     :: tile_dims, tile_coords

  tile_dims = 0

  ! TODO input parameter to say how to split it - rows, columns, etc
  !tile_dims(1) = 1

  ! get good split for tiles
  CALL MPI_DIMS_CREATE(tiles_per_task, 2, tile_dims, err)

  tiles_x = tile_dims(2)
  tiles_y = tile_dims(1)

  delta_x=x_cells/tiles_x
  delta_y=y_cells/tiles_y
  mod_x=MOD(x_cells,tiles_x)
  mod_y=MOD(y_cells,tiles_y)

  DO k=0,tile_dims(2)-1
    DO j=0,tile_dims(1)-1
      t = k*tile_dims(1) + j + 1

      tile_coords(1) = j
      tile_coords(2) = k

      chunk%tiles(t)%tile_coords(1) = j+1
      chunk%tiles(t)%tile_coords(2) = k+1

      chunk%tiles(t)%def_tile_coords(1) = mpi_coords(1)*tile_dims(1) + (j + 1)
      chunk%tiles(t)%def_tile_coords(2) = mpi_coords(2)*tile_dims(2) + (k + 1)

      chunk%tiles(t)%left = chunk%left + tile_coords(1)*delta_x
      if (tile_coords(1) .le. mod_x) then
        chunk%tiles(t)%left = chunk%tiles(t)%left + tile_coords(1)
      else
        chunk%tiles(t)%left = chunk%tiles(t)%left + mod_x
      endif

      chunk%tiles(t)%right = chunk%tiles(t)%left+delta_x - 1
      if (tile_coords(1) .lt. mod_x) then
        chunk%tiles(t)%right = chunk%tiles(t)%right + 1
      endif

      chunk%tiles(t)%bottom = chunk%bottom + tile_coords(2)*delta_y
      if (tile_coords(2) .le. mod_y) then
        chunk%tiles(t)%bottom = chunk%tiles(t)%bottom + tile_coords(2)
      else
        chunk%tiles(t)%bottom = chunk%tiles(t)%bottom + mod_y
      endif

      chunk%tiles(t)%top = chunk%tiles(t)%bottom+delta_y - 1
      if (tile_coords(2) .lt. mod_y) then
        chunk%tiles(t)%top = chunk%tiles(t)%top + 1
      endif

      chunk%tiles(t)%tile_neighbours = EXTERNAL_FACE

      IF (j .GT. 0) THEN
        chunk%tiles(t)%tile_neighbours(CHUNK_LEFT) = (k+0)*tile_dims(1) + (j-1) + 1
      ENDIF

      IF (j .LT. tile_dims(1)-1) THEN
        chunk%tiles(t)%tile_neighbours(CHUNK_RIGHT) = (k+0)*tile_dims(1) + (j+1) + 1
      ENDIF

      IF (k .GT. 0) THEN
        chunk%tiles(t)%tile_neighbours(CHUNK_BOTTOM) = (k-1)*tile_dims(1) + (j+0) + 1
      ENDIF

      IF (k .LT. tile_dims(2)-1) THEN
        chunk%tiles(t)%tile_neighbours(CHUNK_TOP) = (k+1)*tile_dims(1) + (j+0) + 1
      ENDIF

    ENDDO
  ENDDO

  chunk%tile_dims = tile_dims

  IF (parallel%boss)THEN
    WRITE(g_out,*)
    WRITE(g_out,*)"Decomposing each chunk into ",tiles_x," by ",tiles_y," tiles"
    WRITE(g_out,*)
  ENDIF

END SUBROUTINE tea_decompose_tiles

SUBROUTINE tea_allocate_buffers()

  IMPLICIT NONE

  INTEGER      :: bt_size, lr_size
  INTEGER,PARAMETER :: num_buffered=NUM_FIELDS

  INTEGER           :: allocate_extra_size

  allocate_extra_size = max(2, halo_exchange_depth)

  lr_size = num_buffered*(chunk%y_cells + 2*allocate_extra_size)*halo_exchange_depth
  bt_size = num_buffered*(chunk%x_cells + 2*allocate_extra_size)*halo_exchange_depth

  ! Unallocated buffers for external boundaries caused issues on some systems so they are now
  !  all allocated
  !IF (chunk%chunk_neighbours(CHUNK_LEFT).NE.EXTERNAL_FACE) THEN
    ALLOCATE(chunk%left_snd_buffer(lr_size))
    ALLOCATE(chunk%left_rcv_buffer(lr_size))
  !ENDIF
  !IF (chunk%chunk_neighbours(CHUNK_RIGHT).NE.EXTERNAL_FACE) THEN
    ALLOCATE(chunk%right_snd_buffer(lr_size))
    ALLOCATE(chunk%right_rcv_buffer(lr_size))
  !ENDIF
  !IF (chunk%chunk_neighbours(CHUNK_BOTTOM).NE.EXTERNAL_FACE) THEN
    ALLOCATE(chunk%bottom_snd_buffer(bt_size))
    ALLOCATE(chunk%bottom_rcv_buffer(bt_size))
  !ENDIF
  !IF (chunk%chunk_neighbours(CHUNK_TOP).NE.EXTERNAL_FACE) THEN
    ALLOCATE(chunk%top_snd_buffer(bt_size))
    ALLOCATE(chunk%top_rcv_buffer(bt_size))
  !ENDIF

  chunk%left_snd_buffer = 0
  chunk%right_snd_buffer = 0
  chunk%bottom_snd_buffer = 0
  chunk%top_snd_buffer = 0

  chunk%left_rcv_buffer = 0
  chunk%right_rcv_buffer = 0
  chunk%bottom_rcv_buffer = 0
  chunk%top_rcv_buffer = 0

END SUBROUTINE tea_allocate_buffers

SUBROUTINE tea_exchange(fields,depth)

  IMPLICIT NONE

    INTEGER      :: fields(NUM_FIELDS),depth, err
    INTEGER      :: left_right_offset(NUM_FIELDS),bottom_top_offset(NUM_FIELDS)
    INTEGER      :: end_pack_index_left_right, end_pack_index_bottom_top,field
    INTEGER      :: message_count_lr, message_count_ud
    INTEGER      :: exchange_size_lr, exchange_size_ud
    INTEGER, DIMENSION(4)                 :: request_lr, request_ud
    INTEGER, DIMENSION(MPI_STATUS_SIZE,4) :: status_lr, status_ud
    LOGICAL :: test_complete

    IF (ALL(chunk%chunk_neighbours .eq. EXTERNAL_FACE)) return

    exchange_size_lr = depth*(chunk%y_cells+2*depth)
    exchange_size_ud = depth*(chunk%x_cells+2*depth)

    request_lr = 0
    message_count_lr = 0
    request_ud = 0
    message_count_ud = 0

    end_pack_index_left_right=0
    end_pack_index_bottom_top=0
    left_right_offset = 0
    bottom_top_offset = 0

    DO field=1,NUM_FIELDS
      IF (fields(field).EQ.1) THEN
        left_right_offset(field)=end_pack_index_left_right
        bottom_top_offset(field)=end_pack_index_bottom_top
        end_pack_index_left_right=end_pack_index_left_right + exchange_size_lr
        end_pack_index_bottom_top=end_pack_index_bottom_top + exchange_size_ud
      ENDIF
    ENDDO

    IF (chunk%chunk_neighbours(CHUNK_LEFT).NE.EXTERNAL_FACE) THEN
      ! do left exchanges
      CALL tea_pack_buffers(fields, depth, CHUNK_LEFT, &
        chunk%left_snd_buffer, left_right_offset)

      !send and recv messagse to the left
      CALL tea_send_recv_message_left(chunk%left_snd_buffer,                      &
                                         chunk%left_rcv_buffer,                      &
                                         end_pack_index_left_right,                    &
                                         1, 2,                                               &
                                         request_lr(message_count_lr+1), request_lr(message_count_lr+2))
      message_count_lr = message_count_lr + 2
    ENDIF

    IF (chunk%chunk_neighbours(CHUNK_RIGHT).NE.EXTERNAL_FACE) THEN
      ! do right exchanges
      CALL tea_pack_buffers(fields, depth, CHUNK_RIGHT, &
        chunk%right_snd_buffer, left_right_offset)

      !send message to the right
      CALL tea_send_recv_message_right(chunk%right_snd_buffer,                     &
                                          chunk%right_rcv_buffer,                     &
                                          end_pack_index_left_right,                    &
                                          2, 1,                                               &
                                          request_lr(message_count_lr+1), request_lr(message_count_lr+2))
      message_count_lr = message_count_lr + 2
    ENDIF

    IF (depth .EQ. 1) THEN
      test_complete = .false.
      ! don't have to transfer now
      CALL MPI_TESTALL(message_count_lr, request_lr, test_complete, status_lr, err)
    ELSE
      test_complete = .true.
      !make a call to wait / sync
      CALL MPI_WAITALL(message_count_lr,request_lr,status_lr,err)
    ENDIF

    IF (test_complete .EQV. .TRUE.) THEN
      !unpack in left direction
      IF (chunk%chunk_neighbours(CHUNK_LEFT).NE.EXTERNAL_FACE) THEN
        CALL tea_unpack_buffers(fields, depth, CHUNK_LEFT, &
          chunk%left_rcv_buffer, left_right_offset)
      ENDIF

      !unpack in right direction
      IF (chunk%chunk_neighbours(CHUNK_RIGHT).NE.EXTERNAL_FACE) THEN
        CALL tea_unpack_buffers(fields, depth, CHUNK_RIGHT, &
          chunk%right_rcv_buffer, left_right_offset)
      ENDIF
    ENDIF

    IF (chunk%chunk_neighbours(CHUNK_BOTTOM).NE.EXTERNAL_FACE) THEN
      ! do bottom exchanges
      CALL tea_pack_buffers(fields, depth, CHUNK_BOTTOM, &
        chunk%bottom_snd_buffer, bottom_top_offset)

      !send message downwards
      CALL tea_send_recv_message_bottom(chunk%bottom_snd_buffer,                     &
                                           chunk%bottom_rcv_buffer,                     &
                                           end_pack_index_bottom_top,                     &
                                           3, 4,                                                &
                                           request_ud(message_count_ud+1), request_ud(message_count_ud+2))
      message_count_ud = message_count_ud + 2
    ENDIF

    IF (chunk%chunk_neighbours(CHUNK_TOP).NE.EXTERNAL_FACE) THEN
      ! do top exchanges
      CALL tea_pack_buffers(fields, depth, CHUNK_TOP, &
        chunk%top_snd_buffer, bottom_top_offset)

      !send message upwards
      CALL tea_send_recv_message_top(chunk%top_snd_buffer,                           &
                                        chunk%top_rcv_buffer,                           &
                                        end_pack_index_bottom_top,                        &
                                        4, 3,                                                   &
                                        request_ud(message_count_ud+1), request_ud(message_count_ud+2))
      message_count_ud = message_count_ud + 2
    ENDIF

    IF (test_complete .EQV. .FALSE.) THEN
      !make a call to wait / sync
      CALL MPI_WAITALL(message_count_lr,request_lr,status_lr,err)

      !unpack in left direction
      IF (chunk%chunk_neighbours(CHUNK_LEFT).NE.EXTERNAL_FACE) THEN
        CALL tea_unpack_buffers(fields, depth, CHUNK_LEFT, &
          chunk%left_rcv_buffer, left_right_offset)
      ENDIF

      !unpack in right direction
      IF (chunk%chunk_neighbours(CHUNK_RIGHT).NE.EXTERNAL_FACE) THEN
        CALL tea_unpack_buffers(fields, depth, CHUNK_RIGHT, &
          chunk%right_rcv_buffer, left_right_offset)
      ENDIF
    ENDIF

    !need to make a call to wait / sync
    CALL MPI_WAITALL(message_count_ud,request_ud,status_ud,err)

    !unpack in top direction
    IF (chunk%chunk_neighbours(CHUNK_TOP).NE.EXTERNAL_FACE) THEN
      CALL tea_unpack_buffers(fields, depth, CHUNK_TOP, &
        chunk%top_rcv_buffer, bottom_top_offset)
    ENDIF

    !unpack in bottom direction
    IF (chunk%chunk_neighbours(CHUNK_BOTTOM).NE.EXTERNAL_FACE) THEN
      CALL tea_unpack_buffers(fields, depth, CHUNK_BOTTOM, &
        chunk%bottom_rcv_buffer, bottom_top_offset)
    ENDIF

END SUBROUTINE tea_exchange

SUBROUTINE tea_send_recv_message_left(left_snd_buffer, left_rcv_buffer,      &
                                         total_size,                     &
                                         tag_send, tag_recv,                    &
                                         req_send, req_recv)

  REAL(KIND=8)    :: left_snd_buffer(:), left_rcv_buffer(:)
  INTEGER         :: left_task
  INTEGER         :: total_size, tag_send, tag_recv, err
  INTEGER         :: req_send, req_recv

  left_task =chunk%chunk_neighbours(CHUNK_LEFT)

  CALL MPI_ISEND(left_snd_buffer,total_size,MPI_DOUBLE_PRECISION,left_task,tag_send &
                ,mpi_cart_comm,req_send,err)

  CALL MPI_IRECV(left_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,left_task,tag_recv &
                ,mpi_cart_comm,req_recv,err)

END SUBROUTINE tea_send_recv_message_left

SUBROUTINE tea_send_recv_message_right(right_snd_buffer, right_rcv_buffer,   &
                                          total_size,                    &
                                          tag_send, tag_recv,                   &
                                          req_send, req_recv)

  IMPLICIT NONE

  REAL(KIND=8) :: right_snd_buffer(:), right_rcv_buffer(:)
  INTEGER      :: right_task
  INTEGER      :: total_size, tag_send, tag_recv, err
  INTEGER      :: req_send, req_recv

  right_task=chunk%chunk_neighbours(CHUNK_RIGHT)

  CALL MPI_ISEND(right_snd_buffer,total_size,MPI_DOUBLE_PRECISION,right_task,tag_send, &
                 mpi_cart_comm,req_send,err)

  CALL MPI_IRECV(right_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,right_task,tag_recv, &
                 mpi_cart_comm,req_recv,err)

END SUBROUTINE tea_send_recv_message_right

SUBROUTINE tea_send_recv_message_top(top_snd_buffer, top_rcv_buffer,     &
                                        total_size,                  &
                                        tag_send, tag_recv,                 &
                                        req_send, req_recv)

    IMPLICIT NONE

    REAL(KIND=8) :: top_snd_buffer(:), top_rcv_buffer(:)
    INTEGER      :: top_task
    INTEGER      :: total_size, tag_send, tag_recv, err
    INTEGER      :: req_send, req_recv

    top_task=chunk%chunk_neighbours(CHUNK_TOP)

    CALL MPI_ISEND(top_snd_buffer,total_size,MPI_DOUBLE_PRECISION,top_task,tag_send, &
                   mpi_cart_comm,req_send,err)

    CALL MPI_IRECV(top_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,top_task,tag_recv, &
                   mpi_cart_comm,req_recv,err)

END SUBROUTINE tea_send_recv_message_top

SUBROUTINE tea_send_recv_message_bottom(bottom_snd_buffer, bottom_rcv_buffer,        &
                                           total_size,                           &
                                           tag_send, tag_recv,                          &
                                           req_send, req_recv)

  IMPLICIT NONE

  REAL(KIND=8) :: bottom_snd_buffer(:), bottom_rcv_buffer(:)
  INTEGER      :: bottom_task
  INTEGER      :: total_size, tag_send, tag_recv, err
  INTEGER      :: req_send, req_recv

  bottom_task=chunk%chunk_neighbours(CHUNK_BOTTOM)

  CALL MPI_ISEND(bottom_snd_buffer,total_size,MPI_DOUBLE_PRECISION,bottom_task,tag_send &
                ,mpi_cart_comm,req_send,err)

  CALL MPI_IRECV(bottom_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,bottom_task,tag_recv &
                ,mpi_cart_comm,req_recv,err)

END SUBROUTINE tea_send_recv_message_bottom

END MODULE tea_module

