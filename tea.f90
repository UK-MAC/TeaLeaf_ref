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

SUBROUTINE tea_decompose(level,x_cells,y_cells)

  ! This decomposes the mesh into a number of chunks.

  IMPLICIT NONE

  INTEGER :: level,x_cells,y_cells
  INTEGER :: delta_x,delta_y

  INTEGER  :: chunk_x,chunk_y,mod_x,mod_y

  INTEGER  :: err

  ! Get destinations/sources
  CALL mpi_cart_shift(mpi_cart_comm, 1, 1,      &
    chunk(level)%chunk_neighbours(CHUNK_BOTTOM),   &
    chunk(level)%chunk_neighbours(CHUNK_TOP),      &
    err)
  CALL mpi_cart_shift(mpi_cart_comm, 0, 1,      &
    chunk(level)%chunk_neighbours(CHUNK_LEFT),     &
    chunk(level)%chunk_neighbours(CHUNK_RIGHT),    &
    err)

  WHERE (chunk(level)%chunk_neighbours .EQ. MPI_PROC_NULL)
    chunk(level)%chunk_neighbours = EXTERNAL_FACE
  END WHERE

  chunk_x = mpi_dims(1)
  chunk_y = mpi_dims(2)

  delta_x=x_cells/chunk_x
  delta_y=y_cells/chunk_y
  mod_x=MOD(x_cells,chunk_x)
  mod_y=MOD(y_cells,chunk_y)

  chunk(level)%left = mpi_coords(1)*delta_x + 1
  if (mpi_coords(1) .le. mod_x) then
    chunk(level)%left = chunk(level)%left + mpi_coords(1)
  else
    chunk(level)%left = chunk(level)%left + mod_x
  endif
  chunk(level)%right = chunk(level)%left+delta_x - 1
  if (mpi_coords(1) .lt. mod_x) then
    chunk(level)%right = chunk(level)%right + 1
  endif

  chunk(level)%bottom = mpi_coords(2)*delta_y + 1
  if (mpi_coords(2) .le. mod_y) then
    chunk(level)%bottom = chunk(level)%bottom + mpi_coords(2)
  else
    chunk(level)%bottom = chunk(level)%bottom + mod_y
  endif
  chunk(level)%top = chunk(level)%bottom+delta_y - 1
  if (mpi_coords(2) .lt. mod_y) then
    chunk(level)%top = chunk(level)%top + 1
  endif

  IF (parallel%boss)THEN
    WRITE(g_out,*)
    WRITE(g_out,*)"Mesh ratio of ",REAL(x_cells)/REAL(y_cells)
    WRITE(g_out,*)"Decomposing the mesh into ",chunk_x," by ",chunk_y," chunks"
    WRITE(g_out,*)
  ENDIF

END SUBROUTINE tea_decompose

SUBROUTINE tea_decompose_tiles(level,x_cells, y_cells)

  IMPLICIT NONE

  INTEGER :: level,x_cells,y_cells,xs_cells,ys_cells

  INTEGER :: delta_x,delta_y
  INTEGER  :: tiles_x,tiles_y,mod_x,mod_y

  INTEGER :: err, j, k, t
  INTEGER,PARAMETER :: sub_tile_nx=2, sub_tile_ny=2

  INTEGER :: best_fit_i,i
  REAL(KIND=8) :: best_fit_v,fit_v

  chunk(level)%tile_dims = 0

  ! TODO input parameter to say how to split it - rows, columns, etc
  !chunk(level)%tile_dims(1) = 1

  ! get good split for tiles
  !write(6,*) "MPI_DIMS_CREATE:tile input :",level, tiles_per_task, 2
  !CALL MPI_DIMS_CREATE(tiles_per_task, 2, chunk(level)%tile_dims, err)
  !write(6,*) "MPI_DIMS_CREATE:tile output:",level, chunk(level)%tile_dims, err
  best_fit_v=0.0_8
  best_fit_i=0
  !write(6,*) x_cells,y_cells,tiles_per_task
  IF (MOD(x_cells*y_cells,tiles_per_task) /= 0) THEN
    DO i=1,tiles_per_task-1
      IF (MOD(x_cells*y_cells,tiles_per_task+i) == 0) THEN
        tiles_per_task=tiles_per_task+i
        EXIT
      ENDIF
      IF (MOD(x_cells*y_cells,tiles_per_task-i) == 0) THEN
        tiles_per_task=tiles_per_task-i
        EXIT
      ENDIF
    ENDDO
    WRITE(6,*) "modified task_per_tile to match domains size:",tiles_per_task
  ENDIF
  DO i=1,tiles_per_task
    IF (mod(tiles_per_task,i) /= 0) CYCLE
    j=tiles_per_task/i
    IF (mod(x_cells,i) /= 0) CYCLE
    IF (mod(y_cells,j) /= 0) CYCLE
    fit_v=real(min(x_cells/i,y_cells/j),8)/real(max(x_cells/i,y_cells/j),8)
    IF (fit_v > best_fit_v) THEN
      best_fit_v=fit_v
      best_fit_i=i
    ENDIF
  ENDDO
  IF (best_fit_i == 0) STOP
  chunk(level)%tile_dims(1)=best_fit_i
  chunk(level)%tile_dims(2)=tiles_per_task/chunk(level)%tile_dims(1)
!  WRITE(6,*) "tiles_per_task optimisation    :",tiles_per_task,x_cells,y_cells, &
!    " best fit:",best_fit_i,tiles_per_task/chunk(level)%tile_dims(1),fit_v

  ! get good split for sub-tiles
  IF (level < 2) THEN
!    write(6,*) "MPI_DIMS_CREATE:sub_tile input :",level, sub_tiles_per_tile, 2
!    CALL MPI_DIMS_CREATE(sub_tiles_per_tile, 2, chunk(level)%sub_tile_dims, err)
!    write(6,*) "MPI_DIMS_CREATE:sub_tile output:",level, chunk(level)%sub_tile_dims, err
  best_fit_v=0.0_8
  best_fit_i=0
  xs_cells=x_cells/chunk(level)%tile_dims(1); ys_cells=y_cells/chunk(level)%tile_dims(2)
  IF (MOD(xs_cells*ys_cells,sub_tiles_per_tile) /= 0) THEN
    DO i=1,sub_tiles_per_tile-1
      IF (MOD(xs_cells*ys_cells,sub_tiles_per_tile+i) == 0) THEN
        sub_tiles_per_tile=sub_tiles_per_tile+i
        EXIT
      ENDIF
      IF (MOD(xs_cells*ys_cells,sub_tiles_per_tile-i) == 0) THEN
        sub_tiles_per_tile=sub_tiles_per_tile-i
        EXIT
      ENDIF
    ENDDO
    WRITE(6,*) "modified sub_tiles_per_tile to match domains size:",sub_tiles_per_tile
  ENDIF
  DO i=1,sub_tiles_per_tile
    IF (mod(sub_tiles_per_tile,i) /= 0) CYCLE
    j=sub_tiles_per_tile/i
    IF (mod(x_cells/chunk(level)%tile_dims(1),i) /= 0) CYCLE
    IF (mod(y_cells/chunk(level)%tile_dims(2),j) /= 0) CYCLE
    fit_v=real(min(x_cells/chunk(level)%tile_dims(1)/i, &
                   y_cells/chunk(level)%tile_dims(2)/j),8)/ &
          real(max(x_cells/chunk(level)%tile_dims(1)/i, &
                   y_cells/chunk(level)%tile_dims(2)/j),8)
    IF (fit_v > best_fit_v) THEN
      best_fit_v=fit_v
      best_fit_i=i
    ENDIF
  ENDDO
  IF (best_fit_i == 0) STOP
  chunk(level)%sub_tile_dims(1)=best_fit_i
  chunk(level)%sub_tile_dims(2)=sub_tiles_per_tile/chunk(level)%sub_tile_dims(1)
!  WRITE(6,*) "sub_tiles_per_tile optimisation:",sub_tiles_per_tile,xs_cells,ys_cells, &
!    " best fit:",best_fit_i,sub_tiles_per_tile/chunk(level)%sub_tile_dims(1),fit_v
  ENDIF

  tiles_x = chunk(level)%tile_dims(1)
  tiles_y = chunk(level)%tile_dims(2)

  delta_x=x_cells/tiles_x
  delta_y=y_cells/tiles_y
  mod_x=MOD(x_cells,tiles_x)
  mod_y=MOD(y_cells,tiles_y)

  DO j=0,chunk(level)%tile_dims(1)-1
    DO k=0,chunk(level)%tile_dims(2)-1
      t = j*chunk(level)%tile_dims(2) + k + 1

      ! start off with 0-indexed for figuring out where in the grid it is
      chunk(level)%tiles(t)%tile_coords(1) = j
      chunk(level)%tiles(t)%tile_coords(2) = k

      chunk(level)%tiles(t)%left = chunk(level)%left + chunk(level)%tiles(t)%tile_coords(1)*delta_x
      if (chunk(level)%tiles(t)%tile_coords(1) .le. mod_x) then
        chunk(level)%tiles(t)%left = chunk(level)%tiles(t)%left + chunk(level)%tiles(t)%tile_coords(1)
      else
        chunk(level)%tiles(t)%left = chunk(level)%tiles(t)%left + mod_x
      endif

      chunk(level)%tiles(t)%right = chunk(level)%tiles(t)%left+delta_x - 1
      if (chunk(level)%tiles(t)%tile_coords(1) .lt. mod_x) then
        chunk(level)%tiles(t)%right = chunk(level)%tiles(t)%right + 1
      endif

      chunk(level)%tiles(t)%bottom = chunk(level)%bottom + chunk(level)%tiles(t)%tile_coords(2)*delta_y
      if (chunk(level)%tiles(t)%tile_coords(2) .le. mod_y) then
        chunk(level)%tiles(t)%bottom = chunk(level)%tiles(t)%bottom + chunk(level)%tiles(t)%tile_coords(2)
      else
        chunk(level)%tiles(t)%bottom = chunk(level)%tiles(t)%bottom + mod_y
      endif

      chunk(level)%tiles(t)%top = chunk(level)%tiles(t)%bottom+delta_y - 1
      if (chunk(level)%tiles(t)%tile_coords(2) .lt. mod_y) then
        chunk(level)%tiles(t)%top = chunk(level)%tiles(t)%top + 1
      endif

      ! add one to make it into 1 indexed
      chunk(level)%tiles(t)%tile_coords = chunk(level)%tiles(t)%tile_coords + 1

      ! absolute position of tile compared to all other tiles in grid
      chunk(level)%tiles(t)%def_tile_coords(1) = (mpi_coords(1)*chunk(level)%tile_dims(1) + j)*chunk(level)%sub_tile_dims(1) + 1
      chunk(level)%tiles(t)%def_tile_coords(2) = (mpi_coords(2)*chunk(level)%tile_dims(2) + k)*chunk(level)%sub_tile_dims(2) + 1

      chunk(level)%tiles(t)%def_tile_idx = (chunk(level)%tiles(t)%def_tile_coords(2) - 1)*mpi_dims(1)*chunk(level)%tile_dims(1) + &
                                            chunk(level)%tiles(t)%def_tile_coords(1)

      chunk(level)%tiles(t)%tile_neighbours = EXTERNAL_FACE

      IF (j .GT. 0) THEN
        chunk(level)%tiles(t)%tile_neighbours(CHUNK_LEFT) = (j-1)*chunk(level)%tile_dims(2) + (k+0) + 1
      ENDIF

      IF (j .LT. chunk(level)%tile_dims(1)-1) THEN
        chunk(level)%tiles(t)%tile_neighbours(CHUNK_RIGHT) = (j+1)*chunk(level)%tile_dims(2) + (k+0) + 1
      ENDIF

      IF (k .GT. 0) THEN
        chunk(level)%tiles(t)%tile_neighbours(CHUNK_BOTTOM) = (j+0)*chunk(level)%tile_dims(2) + (k-1) + 1
      ENDIF

      IF (k .LT. chunk(level)%tile_dims(2)-1) THEN
        chunk(level)%tiles(t)%tile_neighbours(CHUNK_TOP) = (j+0)*chunk(level)%tile_dims(2) + (k+1) + 1
      ENDIF
    ENDDO
  ENDDO

  IF (parallel%boss)THEN
    WRITE(g_out,*)
    WRITE(g_out,*)"Decomposing each chunk into ",tiles_x," by ",tiles_y," tiles"
    WRITE(g_out,*)
  ENDIF

END SUBROUTINE tea_decompose_tiles

SUBROUTINE tea_allocate_buffers(level)

  IMPLICIT NONE

  INTEGER           :: level

  INTEGER           :: bt_size, lr_size
  INTEGER,PARAMETER :: num_buffered=NUM_FIELDS

  INTEGER           :: allocate_extra_size

  allocate_extra_size = max(2, chunk(level)%halo_exchange_depth)

  lr_size = num_buffered*(chunk(level)%y_cells + 2*allocate_extra_size)*chunk(level)%halo_exchange_depth
  bt_size = num_buffered*(chunk(level)%x_cells + 2*allocate_extra_size)*chunk(level)%halo_exchange_depth
  !write(6,*) lr_size,bt_size

  ! Unallocated buffers for external boundaries caused issues on some systems so they are now
  !  all allocated
  !IF (chunk(level)%chunk_neighbours(CHUNK_LEFT).NE.EXTERNAL_FACE) THEN
    ALLOCATE(chunk(level)%left_snd_buffer(lr_size))
    ALLOCATE(chunk(level)%left_rcv_buffer(lr_size))
  !ENDIF
  !IF (chunk(level)%chunk_neighbours(CHUNK_RIGHT).NE.EXTERNAL_FACE) THEN
    ALLOCATE(chunk(level)%right_snd_buffer(lr_size))
    ALLOCATE(chunk(level)%right_rcv_buffer(lr_size))
  !ENDIF
  !IF (chunk(level)%chunk_neighbours(CHUNK_BOTTOM).NE.EXTERNAL_FACE) THEN
    ALLOCATE(chunk(level)%bottom_snd_buffer(bt_size))
    ALLOCATE(chunk(level)%bottom_rcv_buffer(bt_size))
  !ENDIF
  !IF (chunk(level)%chunk_neighbours(CHUNK_TOP).NE.EXTERNAL_FACE) THEN
    ALLOCATE(chunk(level)%top_snd_buffer(bt_size))
    ALLOCATE(chunk(level)%top_rcv_buffer(bt_size))
  !ENDIF

  chunk(level)%left_snd_buffer = 0
  chunk(level)%right_snd_buffer = 0
  chunk(level)%bottom_snd_buffer = 0
  chunk(level)%top_snd_buffer = 0

  chunk(level)%left_rcv_buffer = 0
  chunk(level)%right_rcv_buffer = 0
  chunk(level)%bottom_rcv_buffer = 0
  chunk(level)%top_rcv_buffer = 0

END SUBROUTINE tea_allocate_buffers

SUBROUTINE tea_exchange(level,fields,depth)

  IMPLICIT NONE

    INTEGER         :: level
    INTEGER      :: fields(NUM_FIELDS),depth, err
    INTEGER      :: left_right_offset(NUM_FIELDS),bottom_top_offset(NUM_FIELDS)
    INTEGER      :: end_pack_index_left_right, end_pack_index_bottom_top,field
    INTEGER      :: message_count_lr, message_count_ud
    INTEGER      :: exchange_size_lr, exchange_size_ud
    INTEGER, DIMENSION(4)                 :: request_lr, request_ud
    INTEGER, DIMENSION(MPI_STATUS_SIZE,4) :: status_lr, status_ud
    LOGICAL :: test_complete

    IF (ALL(chunk(level)%chunk_neighbours .eq. EXTERNAL_FACE)) return

    exchange_size_lr = depth*(chunk(level)%y_cells + 2*depth)
    exchange_size_ud = depth*(chunk(level)%x_cells + 2*depth)

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

    IF (chunk(level)%chunk_neighbours(CHUNK_LEFT).NE.EXTERNAL_FACE) THEN
      ! do left exchanges
      CALL tea_pack_buffers(level, fields, depth, CHUNK_LEFT, &
        chunk(level)%left_snd_buffer, left_right_offset)

      !send and recv messagse to the left
      CALL tea_send_recv_message_left(level, chunk(level)%left_snd_buffer,                      &
                                         chunk(level)%left_rcv_buffer,                      &
                                         end_pack_index_left_right,                    &
                                         1, 2,                                               &
                                         request_lr(message_count_lr+1), request_lr(message_count_lr+2))
      message_count_lr = message_count_lr + 2
    ENDIF

    IF (chunk(level)%chunk_neighbours(CHUNK_RIGHT).NE.EXTERNAL_FACE) THEN
      ! do right exchanges
      CALL tea_pack_buffers(level, fields, depth, CHUNK_RIGHT, &
        chunk(level)%right_snd_buffer, left_right_offset)

      !send message to the right
      CALL tea_send_recv_message_right(level, chunk(level)%right_snd_buffer,                     &
                                          chunk(level)%right_rcv_buffer,                     &
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
      IF (chunk(level)%chunk_neighbours(CHUNK_LEFT).NE.EXTERNAL_FACE) THEN
        CALL tea_unpack_buffers(level, fields, depth, CHUNK_LEFT, &
          chunk(level)%left_rcv_buffer, left_right_offset)
      ENDIF

      !unpack in right direction
      IF (chunk(level)%chunk_neighbours(CHUNK_RIGHT).NE.EXTERNAL_FACE) THEN
        CALL tea_unpack_buffers(level, fields, depth, CHUNK_RIGHT, &
          chunk(level)%right_rcv_buffer, left_right_offset)
      ENDIF
    ENDIF

    IF (chunk(level)%chunk_neighbours(CHUNK_BOTTOM).NE.EXTERNAL_FACE) THEN
      ! do bottom exchanges
      CALL tea_pack_buffers(level, fields, depth, CHUNK_BOTTOM, &
        chunk(level)%bottom_snd_buffer, bottom_top_offset)

      !send message downwards
      CALL tea_send_recv_message_bottom(level, chunk(level)%bottom_snd_buffer,                     &
                                           chunk(level)%bottom_rcv_buffer,                     &
                                           end_pack_index_bottom_top,                     &
                                           3, 4,                                                &
                                           request_ud(message_count_ud+1), request_ud(message_count_ud+2))
      message_count_ud = message_count_ud + 2
    ENDIF

    IF (chunk(level)%chunk_neighbours(CHUNK_TOP).NE.EXTERNAL_FACE) THEN
      ! do top exchanges
      CALL tea_pack_buffers(level, fields, depth, CHUNK_TOP, &
        chunk(level)%top_snd_buffer, bottom_top_offset)

      !send message upwards
      CALL tea_send_recv_message_top(level, chunk(level)%top_snd_buffer,                           &
                                        chunk(level)%top_rcv_buffer,                           &
                                        end_pack_index_bottom_top,                        &
                                        4, 3,                                                   &
                                        request_ud(message_count_ud+1), request_ud(message_count_ud+2))
      message_count_ud = message_count_ud + 2
    ENDIF

    IF (test_complete .EQV. .FALSE.) THEN
      !make a call to wait / sync
      CALL MPI_WAITALL(message_count_lr,request_lr,status_lr,err)

      !unpack in left direction
      IF (chunk(level)%chunk_neighbours(CHUNK_LEFT).NE.EXTERNAL_FACE) THEN
        CALL tea_unpack_buffers(level, fields, depth, CHUNK_LEFT, &
          chunk(level)%left_rcv_buffer, left_right_offset)
      ENDIF

      !unpack in right direction
      IF (chunk(level)%chunk_neighbours(CHUNK_RIGHT).NE.EXTERNAL_FACE) THEN
        CALL tea_unpack_buffers(level, fields, depth, CHUNK_RIGHT, &
          chunk(level)%right_rcv_buffer, left_right_offset)
      ENDIF
    ENDIF

    !need to make a call to wait / sync
    CALL MPI_WAITALL(message_count_ud,request_ud,status_ud,err)

    !unpack in top direction
    IF (chunk(level)%chunk_neighbours(CHUNK_TOP).NE.EXTERNAL_FACE) THEN
      CALL tea_unpack_buffers(level, fields, depth, CHUNK_TOP, &
        chunk(level)%top_rcv_buffer, bottom_top_offset)
    ENDIF

    !unpack in bottom direction
    IF (chunk(level)%chunk_neighbours(CHUNK_BOTTOM).NE.EXTERNAL_FACE) THEN
      CALL tea_unpack_buffers(level, fields, depth, CHUNK_BOTTOM, &
        chunk(level)%bottom_rcv_buffer, bottom_top_offset)
    ENDIF

END SUBROUTINE tea_exchange

SUBROUTINE tea_send_recv_message_left(level,left_snd_buffer, left_rcv_buffer,      &
                                         total_size,                     &
                                         tag_send, tag_recv,                    &
                                         req_send, req_recv)

  INTEGER         :: level
  REAL(KIND=8)    :: left_snd_buffer(:), left_rcv_buffer(:)
  INTEGER         :: left_task
  INTEGER         :: total_size, tag_send, tag_recv, err
  INTEGER         :: req_send, req_recv

  left_task =chunk(level)%chunk_neighbours(CHUNK_LEFT)

  CALL MPI_ISEND(left_snd_buffer,total_size,MPI_DOUBLE_PRECISION,left_task,tag_send &
                ,mpi_cart_comm,req_send,err)

  CALL MPI_IRECV(left_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,left_task,tag_recv &
                ,mpi_cart_comm,req_recv,err)

END SUBROUTINE tea_send_recv_message_left

SUBROUTINE tea_send_recv_message_right(level, right_snd_buffer, right_rcv_buffer,   &
                                          total_size,                    &
                                          tag_send, tag_recv,                   &
                                          req_send, req_recv)

  IMPLICIT NONE

  INTEGER         :: level
  REAL(KIND=8) :: right_snd_buffer(:), right_rcv_buffer(:)
  INTEGER      :: right_task
  INTEGER      :: total_size, tag_send, tag_recv, err
  INTEGER      :: req_send, req_recv

  right_task=chunk(level)%chunk_neighbours(CHUNK_RIGHT)

  CALL MPI_ISEND(right_snd_buffer,total_size,MPI_DOUBLE_PRECISION,right_task,tag_send, &
                 mpi_cart_comm,req_send,err)

  CALL MPI_IRECV(right_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,right_task,tag_recv, &
                 mpi_cart_comm,req_recv,err)

END SUBROUTINE tea_send_recv_message_right

SUBROUTINE tea_send_recv_message_top(level, top_snd_buffer, top_rcv_buffer,     &
                                        total_size,                  &
                                        tag_send, tag_recv,                 &
                                        req_send, req_recv)

    IMPLICIT NONE

    INTEGER         :: level
    REAL(KIND=8) :: top_snd_buffer(:), top_rcv_buffer(:)
    INTEGER      :: top_task
    INTEGER      :: total_size, tag_send, tag_recv, err
    INTEGER      :: req_send, req_recv

    top_task=chunk(level)%chunk_neighbours(CHUNK_TOP)

    CALL MPI_ISEND(top_snd_buffer,total_size,MPI_DOUBLE_PRECISION,top_task,tag_send, &
                   mpi_cart_comm,req_send,err)

    CALL MPI_IRECV(top_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,top_task,tag_recv, &
                   mpi_cart_comm,req_recv,err)

END SUBROUTINE tea_send_recv_message_top

SUBROUTINE tea_send_recv_message_bottom(level, bottom_snd_buffer, bottom_rcv_buffer,        &
                                           total_size,                           &
                                           tag_send, tag_recv,                          &
                                           req_send, req_recv)

  IMPLICIT NONE

  INTEGER         :: level
  REAL(KIND=8) :: bottom_snd_buffer(:), bottom_rcv_buffer(:)
  INTEGER      :: bottom_task
  INTEGER      :: total_size, tag_send, tag_recv, err
  INTEGER      :: req_send, req_recv

  bottom_task=chunk(level)%chunk_neighbours(CHUNK_BOTTOM)

  CALL MPI_ISEND(bottom_snd_buffer,total_size,MPI_DOUBLE_PRECISION,bottom_task,tag_send &
                ,mpi_cart_comm,req_send,err)

  CALL MPI_IRECV(bottom_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,bottom_task,tag_recv &
                ,mpi_cart_comm,req_recv,err)

END SUBROUTINE tea_send_recv_message_bottom

END MODULE tea_module

