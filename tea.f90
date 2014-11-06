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

  USE data_module
  USE definitions_module
  USE MPI

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_barrier

  INTEGER :: err

  CALL MPI_BARRIER(MPI_COMM_WORLD,err)

END SUBROUTINE tea_barrier

SUBROUTINE tea_abort

  INTEGER :: ierr,err

  CALL MPI_ABORT(MPI_COMM_WORLD,ierr,err)

END SUBROUTINE tea_abort

SUBROUTINE tea_finalize

  INTEGER :: err

  CLOSE(g_out)
  CALL FLUSH(0)
  CALL FLUSH(6)
  CALL FLUSH(g_out)
  CALL MPI_FINALIZE(err)

END SUBROUTINE tea_finalize

SUBROUTINE tea_init_comms

  IMPLICIT NONE

  INTEGER :: err,rank,size

  rank=0
  size=1

  CALL MPI_INIT(err)

  CALL MPI_COMM_RANK(MPI_COMM_WORLD,rank,err)
  CALL MPI_COMM_SIZE(MPI_COMM_WORLD,size,err)

  parallel%parallel=.TRUE.
  parallel%task=rank

  IF(rank.EQ.0) THEN
    parallel%boss=.TRUE.
  ENDIF

  parallel%boss_task=0
  parallel%max_task=size

END SUBROUTINE tea_init_comms

SUBROUTINE tea_get_num_chunks(count)

  IMPLICIT NONE

  INTEGER :: count

! Should be changed so there can be more than one chunk per mpi task

  count=parallel%max_task

END SUBROUTINE tea_get_num_chunks

SUBROUTINE tea_decompose(x_cells,y_cells,left,right,bottom,top)

  ! This decomposes the mesh into a number of chunks.
  ! The number of chunks may be a multiple of the number of mpi tasks
  ! Doesn't always return the best split if there are few factors
  ! All factors need to be stored and the best picked. But its ok for now

  IMPLICIT NONE

  INTEGER :: x_cells,y_cells,left(:),right(:),top(:),bottom(:)
  INTEGER :: c,delta_x,delta_y

  REAL(KIND=8) :: mesh_ratio,factor_x,factor_y
  INTEGER  :: chunk_x,chunk_y,mod_x,mod_y,split_found

  INTEGER  :: cx,cy,chunk,add_x,add_y,add_x_prev,add_y_prev

  ! 2D Decomposition of the mesh

  mesh_ratio=real(x_cells)/real(y_cells)

  chunk_x=number_of_chunks
  chunk_y=1

  split_found=0 ! Used to detect 1D decomposition
  DO c=1,number_of_chunks
    IF (MOD(number_of_chunks,c).EQ.0) THEN
      factor_x=number_of_chunks/real(c)
      factor_y=c
      !Compare the factor ratio with the mesh ratio
      IF(factor_x/factor_y.LE.mesh_ratio) THEN
        chunk_y=c
        chunk_x=number_of_chunks/c
        split_found=1
        EXIT
      ENDIF
    ENDIF
  ENDDO

  IF(split_found.EQ.0.OR.chunk_y.EQ.number_of_chunks) THEN ! Prime number or 1D decomp detected
    IF(mesh_ratio.GE.1.0) THEN
      chunk_x=number_of_chunks
      chunk_y=1
    ELSE
      chunk_x=1
      chunk_y=number_of_chunks
    ENDIF
  ENDIF

  delta_x=x_cells/chunk_x
  delta_y=y_cells/chunk_y
  mod_x=MOD(x_cells,chunk_x)
  mod_y=MOD(y_cells,chunk_y)

  ! Set up chunk mesh ranges and chunk connectivity

    add_x_prev=0
    add_y_prev=0
    chunk=1
    DO cy=1,chunk_y
        DO cx=1,chunk_x
            add_x=0
            add_y=0
            IF(cx.LE.mod_x)add_x=1
            IF(cy.LE.mod_y)add_y=1

            IF (chunk .EQ. parallel%task+1) THEN
                left(1)   = (cx-1)*delta_x+1+add_x_prev
                right(1)  = left(1)+delta_x-1+add_x
                bottom(1) = (cy-1)*delta_y+1+add_y_prev
                top(1)    = bottom(1)+delta_y-1+add_y

                chunks(1)%chunk_neighbours(chunk_left)=chunk_x*(cy-1)+cx-1
                chunks(1)%chunk_neighbours(chunk_right)=chunk_x*(cy-1)+cx+1
                chunks(1)%chunk_neighbours(chunk_bottom)=chunk_x*(cy-2)+cx
                chunks(1)%chunk_neighbours(chunk_top)=chunk_x*(cy)+cx

                IF(cx.EQ.1)       chunks(1)%chunk_neighbours(chunk_left)=external_face
                IF(cx.EQ.chunk_x) chunks(1)%chunk_neighbours(chunk_right)=external_face
                IF(cy.EQ.1)       chunks(1)%chunk_neighbours(chunk_bottom)=external_face
                IF(cy.EQ.chunk_y) chunks(1)%chunk_neighbours(chunk_top)=external_face
            ENDIF

            IF(cx.LE.mod_x)add_x_prev=add_x_prev+1
            chunk=chunk+1
        ENDDO
        add_x_prev=0
        IF(cy.LE.mod_y)add_y_prev=add_y_prev+1
    ENDDO

  IF(parallel%boss)THEN
    WRITE(g_out,*)
    WRITE(g_out,*)"Mesh ratio of ",mesh_ratio
    WRITE(g_out,*)"Decomposing the mesh into ",chunk_x," by ",chunk_y," chunks"
    WRITE(g_out,*)
  ENDIF

END SUBROUTINE tea_decompose

SUBROUTINE tea_allocate_buffers(chunk)

  IMPLICIT NONE

  INTEGER      :: chunk

  lr_pack_buffer_size = (chunks(chunk)%field%y_max+5)
  bt_pack_buffer_size = (chunks(chunk)%field%x_max+5)

  ! Unallocated buffers for external boundaries caused issues on some systems so they are now
  !  all allocated
  IF(parallel%task.EQ.chunks(chunk)%task)THEN
    !IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%left_snd_buffer(6*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%left_rcv_buffer(6*(chunks(chunk)%field%y_max+5)))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%right_snd_buffer(6*(chunks(chunk)%field%y_max+5)))
      ALLOCATE(chunks(chunk)%right_rcv_buffer(6*(chunks(chunk)%field%y_max+5)))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%bottom_snd_buffer(6*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%bottom_rcv_buffer(6*(chunks(chunk)%field%x_max+5)))
    !ENDIF
    !IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
      ALLOCATE(chunks(chunk)%top_snd_buffer(6*(chunks(chunk)%field%x_max+5)))
      ALLOCATE(chunks(chunk)%top_rcv_buffer(6*(chunks(chunk)%field%x_max+5)))
    !ENDIF
  ENDIF

END SUBROUTINE tea_allocate_buffers

SUBROUTINE tea_exchange(fields,depth)

  IMPLICIT NONE

    INTEGER      :: fields(NUM_FIELDS),depth, chunk
    INTEGER      :: left_right_offset(NUM_FIELDS),bottom_top_offset(NUM_FIELDS)
    INTEGER      :: request(4)
    INTEGER      :: message_count,err
    INTEGER      :: status(MPI_STATUS_SIZE,4)
    INTEGER      :: end_pack_index_left_right, end_pack_index_bottom_top,field

    ! Assuming 1 patch per task, this will be changed

    request=0
    message_count=0

    chunk = 1

    end_pack_index_left_right=0
    end_pack_index_bottom_top=0
    left_right_offset = 0
    bottom_top_offset = 0
    DO field=1,NUM_FIELDS
      IF(fields(field).EQ.1) THEN
        left_right_offset(field)=end_pack_index_left_right
        bottom_top_offset(field)=end_pack_index_bottom_top
        end_pack_index_left_right=end_pack_index_left_right + depth*lr_pack_buffer_size
        end_pack_index_bottom_top=end_pack_index_bottom_top + depth*bt_pack_buffer_size
      ENDIF
    ENDDO

    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
      ! do left exchanges
      CALL tea_pack_left(chunk, fields, depth, left_right_offset)

      !send and recv messagse to the left
      CALL tea_send_recv_message_left(chunks(chunk)%left_snd_buffer,                      &
                                         chunks(chunk)%left_rcv_buffer,                      &
                                         chunk,end_pack_index_left_right,                    &
                                         1, 2,                                               &
                                         request(message_count+1), request(message_count+2))
      message_count = message_count + 2
    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
      ! do right exchanges
      CALL tea_pack_right(chunk, fields, depth, left_right_offset)

      !send message to the right
      CALL tea_send_recv_message_right(chunks(chunk)%right_snd_buffer,                     &
                                          chunks(chunk)%right_rcv_buffer,                     &
                                          chunk,end_pack_index_left_right,                    &
                                          2, 1,                                               &
                                          request(message_count+1), request(message_count+2))
      message_count = message_count + 2
    ENDIF

    !make a call to wait / sync
    CALL MPI_WAITALL(message_count,request,status,err)

    !unpack in left direction
    IF(chunks(chunk)%chunk_neighbours(chunk_left).NE.external_face) THEN
        CALL tea_unpack_left(fields, chunk, depth,                      &
                                chunks(chunk)%left_rcv_buffer,                      &
                                left_right_offset)
    ENDIF


    !unpack in right direction
    IF(chunks(chunk)%chunk_neighbours(chunk_right).NE.external_face) THEN
        CALL tea_unpack_right(fields, chunk, depth,                     &
                                 chunks(chunk)%right_rcv_buffer,                      &
                                 left_right_offset)
    ENDIF

    message_count = 0
    request = 0

    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
      ! do bottom exchanges
      CALL tea_pack_bottom(chunk, fields, depth, bottom_top_offset)

      !send message downwards
      CALL tea_send_recv_message_bottom(chunks(chunk)%bottom_snd_buffer,                     &
                                           chunks(chunk)%bottom_rcv_buffer,                     &
                                           chunk,end_pack_index_bottom_top,                     &
                                           3, 4,                                                &
                                           request(message_count+1), request(message_count+2))
      message_count = message_count + 2
    ENDIF

    IF(chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face) THEN
      ! do top exchanges
      CALL tea_pack_top(chunk, fields, depth, bottom_top_offset)

      !send message upwards
      CALL tea_send_recv_message_top(chunks(chunk)%top_snd_buffer,                           &
                                        chunks(chunk)%top_rcv_buffer,                           &
                                        chunk,end_pack_index_bottom_top,                        &
                                        4, 3,                                                   &
                                        request(message_count+1), request(message_count+2))
      message_count = message_count + 2
    ENDIF

    !need to make a call to wait / sync
    CALL MPI_WAITALL(message_count,request,status,err)

    !unpack in top direction
    IF( chunks(chunk)%chunk_neighbours(chunk_top).NE.external_face ) THEN
        CALL tea_unpack_top(fields, chunk, depth,                       &
                               chunks(chunk)%top_rcv_buffer,                      &
                               bottom_top_offset)
    ENDIF

    !unpack in bottom direction
    IF(chunks(chunk)%chunk_neighbours(chunk_bottom).NE.external_face) THEN
        CALL tea_unpack_bottom(fields, chunk, depth,                   &
                                 chunks(chunk)%bottom_rcv_buffer,                      &
                                 bottom_top_offset)
    ENDIF

END SUBROUTINE tea_exchange

SUBROUTINE tea_pack_left(chunk, fields, depth, left_right_offset)

  USE pack_kernel_module

  IMPLICIT NONE

  INTEGER      :: fields(:),depth, chunk
  INTEGER      :: left_right_offset(:)

  IF(fields(FIELD_DENSITY).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%density,                 &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, CELL_DATA,                             &
                                    left_right_offset(FIELD_DENSITY))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%energy0,                  &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, CELL_DATA,                             &
                                    left_right_offset(FIELD_ENERGY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%energy1,                  &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, CELL_DATA,                             &
                                    left_right_offset(FIELD_ENERGY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_P).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%vector_p,                  &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, CELL_DATA,                             &
                                    left_right_offset(FIELD_P))
    ENDIF
  ENDIF
  IF(fields(FIELD_U).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%u,                  &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, CELL_DATA,                             &
                                    left_right_offset(FIELD_U))
    ENDIF
  ENDIF
  IF(fields(FIELD_SD).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_left(chunks(chunk)%field%x_min,                    &
                                    chunks(chunk)%field%x_max,                    &
                                    chunks(chunk)%field%y_min,                    &
                                    chunks(chunk)%field%y_max,                    &
                                    chunks(chunk)%field%vector_sd,                  &
                                    chunks(chunk)%left_snd_buffer,                &
                                    CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                    depth, CELL_DATA,                             &
                                    left_right_offset(FIELD_SD))
    ENDIF
  ENDIF

END SUBROUTINE tea_pack_left

SUBROUTINE tea_send_recv_message_left(left_snd_buffer, left_rcv_buffer,      &
                                         chunk, total_size,                     &
                                         tag_send, tag_recv,                    &
                                         req_send, req_recv)

  REAL(KIND=8)    :: left_snd_buffer(:), left_rcv_buffer(:)
  INTEGER         :: left_task
  INTEGER         :: chunk
  INTEGER         :: total_size, tag_send, tag_recv, err
  INTEGER         :: req_send, req_recv

  left_task =chunks(chunk)%chunk_neighbours(chunk_left) - 1

  CALL MPI_ISEND(left_snd_buffer,total_size,MPI_DOUBLE_PRECISION,left_task,tag_send &
                ,MPI_COMM_WORLD,req_send,err)

  CALL MPI_IRECV(left_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,left_task,tag_recv &
                ,MPI_COMM_WORLD,req_recv,err)

END SUBROUTINE tea_send_recv_message_left

SUBROUTINE tea_unpack_left(fields, chunk, depth,                         &
                              left_rcv_buffer,                              &
                              left_right_offset)

  USE pack_kernel_module

  IMPLICIT NONE

  INTEGER         :: fields(:), chunk, depth
  INTEGER         :: left_right_offset(:)
  REAL(KIND=8)    :: left_rcv_buffer(:)


  IF(fields(FIELD_DENSITY).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%density,                 &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_DENSITY))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%energy0,                  &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_ENERGY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%energy1,                  &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_ENERGY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_P).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%vector_p,                  &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_P))
    ENDIF
  ENDIF
  IF(fields(FIELD_U).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%u,                  &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_U))
    ENDIF
  ENDIF
  IF(fields(FIELD_SD).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_left(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%vector_sd,                  &
                                      chunks(chunk)%left_rcv_buffer,                &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      left_right_offset(FIELD_SD))
    ENDIF
  ENDIF

END SUBROUTINE tea_unpack_left

SUBROUTINE tea_pack_right(chunk, fields, depth, left_right_offset)

  USE pack_kernel_module

  IMPLICIT NONE

  INTEGER        :: chunk, fields(:), depth, tot_packr, left_right_offset(:)

  IF(fields(FIELD_DENSITY).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%density,                 &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     left_right_offset(FIELD_DENSITY))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%energy0,                  &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     left_right_offset(FIELD_ENERGY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%energy1,                  &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     left_right_offset(FIELD_ENERGY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_P).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%vector_p,                  &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     left_right_offset(FIELD_P))
    ENDIF
  ENDIF
  IF(fields(FIELD_U).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%u,                  &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     left_right_offset(FIELD_U))
    ENDIF
  ENDIF
  IF(fields(FIELD_SD).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_right(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%vector_sd,                  &
                                     chunks(chunk)%right_snd_buffer,               &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     left_right_offset(FIELD_SD))
    ENDIF
  ENDIF

END SUBROUTINE tea_pack_right

SUBROUTINE tea_send_recv_message_right(right_snd_buffer, right_rcv_buffer,   &
                                          chunk, total_size,                    &
                                          tag_send, tag_recv,                   &
                                          req_send, req_recv)

  IMPLICIT NONE

  REAL(KIND=8) :: right_snd_buffer(:), right_rcv_buffer(:)
  INTEGER      :: right_task
  INTEGER      :: chunk
  INTEGER      :: total_size, tag_send, tag_recv, err
  INTEGER      :: req_send, req_recv

  right_task=chunks(chunk)%chunk_neighbours(chunk_right) - 1

  CALL MPI_ISEND(right_snd_buffer,total_size,MPI_DOUBLE_PRECISION,right_task,tag_send, &
                 MPI_COMM_WORLD,req_send,err)

  CALL MPI_IRECV(right_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,right_task,tag_recv, &
                 MPI_COMM_WORLD,req_recv,err)

END SUBROUTINE tea_send_recv_message_right

SUBROUTINE tea_unpack_right(fields, chunk, depth,                          &
                               right_rcv_buffer,                              &
                               left_right_offset)

  USE pack_kernel_module

  IMPLICIT NONE

  INTEGER         :: fields(:), chunk, total_in_right_buff, depth, left_right_offset(:)
  REAL(KIND=8)    :: right_rcv_buffer(:)

  IF(fields(FIELD_DENSITY).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%density,                 &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_DENSITY))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%energy0,                  &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_ENERGY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%energy1,                  &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_ENERGY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_P).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%vector_p,                  &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_P))
    ENDIF
  ENDIF
  IF(fields(FIELD_U).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%u,                  &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_U))
    ENDIF
  ENDIF
  IF(fields(FIELD_SD).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_right(chunks(chunk)%field%x_min,                    &
                                       chunks(chunk)%field%x_max,                    &
                                       chunks(chunk)%field%y_min,                    &
                                       chunks(chunk)%field%y_max,                    &
                                       chunks(chunk)%field%vector_sd,                  &
                                       chunks(chunk)%right_rcv_buffer,               &
                                       CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                       depth, CELL_DATA,                             &
                                       left_right_offset(FIELD_SD))
    ENDIF
  ENDIF

END SUBROUTINE tea_unpack_right

SUBROUTINE tea_pack_top(chunk, fields, depth, bottom_top_offset)

  USE pack_kernel_module

  IMPLICIT NONE

  INTEGER        :: chunk, fields(:), depth, bottom_top_offset(:)

  IF(fields(FIELD_DENSITY).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%density,                 &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, CELL_DATA,                             &
                                   bottom_top_offset(FIELD_DENSITY))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%energy0,                  &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, CELL_DATA,                             &
                                   bottom_top_offset(FIELD_ENERGY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%energy1,                  &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, CELL_DATA,                             &
                                   bottom_top_offset(FIELD_ENERGY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_P).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%vector_p,                  &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, CELL_DATA,                             &
                                   bottom_top_offset(FIELD_P))
    ENDIF
  ENDIF
  IF(fields(FIELD_U).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%u,                  &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, CELL_DATA,                             &
                                   bottom_top_offset(FIELD_U))
    ENDIF
  ENDIF
  IF(fields(FIELD_SD).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_top(chunks(chunk)%field%x_min,                    &
                                   chunks(chunk)%field%x_max,                    &
                                   chunks(chunk)%field%y_min,                    &
                                   chunks(chunk)%field%y_max,                    &
                                   chunks(chunk)%field%vector_sd,                  &
                                   chunks(chunk)%top_snd_buffer,                 &
                                   CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                   depth, CELL_DATA,                             &
                                   bottom_top_offset(FIELD_SD))
    ENDIF
  ENDIF

END SUBROUTINE tea_pack_top

SUBROUTINE tea_send_recv_message_top(top_snd_buffer, top_rcv_buffer,     &
                                        chunk, total_size,                  &
                                        tag_send, tag_recv,                 &
                                        req_send, req_recv)

    IMPLICIT NONE

    REAL(KIND=8) :: top_snd_buffer(:), top_rcv_buffer(:)
    INTEGER      :: top_task
    INTEGER      :: chunk
    INTEGER      :: total_size, tag_send, tag_recv, err
    INTEGER      :: req_send, req_recv

    top_task=chunks(chunk)%chunk_neighbours(chunk_top) - 1

    CALL MPI_ISEND(top_snd_buffer,total_size,MPI_DOUBLE_PRECISION,top_task,tag_send, &
                   MPI_COMM_WORLD,req_send,err)

    CALL MPI_IRECV(top_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,top_task,tag_recv, &
                   MPI_COMM_WORLD,req_recv,err)

END SUBROUTINE tea_send_recv_message_top

SUBROUTINE tea_unpack_top(fields, chunk, depth,                        &
                             top_rcv_buffer,                              &
                             bottom_top_offset)

  USE pack_kernel_module

  IMPLICIT NONE

  INTEGER         :: fields(:), chunk, total_in_top_buff, depth, bottom_top_offset(:)
  REAL(KIND=8)    :: top_rcv_buffer(:)


  IF(fields(FIELD_DENSITY).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%density,                 &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_DENSITY))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%energy0,                  &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_ENERGY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%energy1,                  &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_ENERGY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_P).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%vector_p,                  &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_P))
    ENDIF
  ENDIF
  IF(fields(FIELD_U).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%u,                  &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_U))
    ENDIF
  ENDIF
  IF(fields(FIELD_SD).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_top(chunks(chunk)%field%x_min,                    &
                                     chunks(chunk)%field%x_max,                    &
                                     chunks(chunk)%field%y_min,                    &
                                     chunks(chunk)%field%y_max,                    &
                                     chunks(chunk)%field%vector_sd,                  &
                                     chunks(chunk)%top_rcv_buffer,                 &
                                     CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                     depth, CELL_DATA,                             &
                                     bottom_top_offset(FIELD_SD))
    ENDIF
  ENDIF

END SUBROUTINE tea_unpack_top

SUBROUTINE tea_pack_bottom(chunk, fields, depth, bottom_top_offset)

  USE pack_kernel_module

  IMPLICIT NONE

  INTEGER        :: chunk, fields(:), depth, tot_packb, bottom_top_offset(:)

  IF(fields(FIELD_DENSITY).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%density,                 &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      bottom_top_offset(FIELD_DENSITY))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%energy0,                  &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      bottom_top_offset(FIELD_ENERGY0))
    ENDIF
  ENDIF
  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%energy1,                  &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      bottom_top_offset(FIELD_ENERGY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_P).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%vector_p,                  &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      bottom_top_offset(FIELD_P))
    ENDIF
  ENDIF
  IF(fields(FIELD_U).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%u,                  &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      bottom_top_offset(FIELD_U))
    ENDIF
  ENDIF
  IF(fields(FIELD_SD).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_pack_message_bottom(chunks(chunk)%field%x_min,                    &
                                      chunks(chunk)%field%x_max,                    &
                                      chunks(chunk)%field%y_min,                    &
                                      chunks(chunk)%field%y_max,                    &
                                      chunks(chunk)%field%vector_sd,                  &
                                      chunks(chunk)%bottom_snd_buffer,              &
                                      CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                      depth, CELL_DATA,                             &
                                      bottom_top_offset(FIELD_SD))
    ENDIF
  ENDIF

END SUBROUTINE tea_pack_bottom

SUBROUTINE tea_send_recv_message_bottom(bottom_snd_buffer, bottom_rcv_buffer,        &
                                           chunk, total_size,                           &
                                           tag_send, tag_recv,                          &
                                           req_send, req_recv)

  IMPLICIT NONE

  REAL(KIND=8) :: bottom_snd_buffer(:), bottom_rcv_buffer(:)
  INTEGER      :: bottom_task
  INTEGER      :: chunk
  INTEGER      :: total_size, tag_send, tag_recv, err
  INTEGER      :: req_send, req_recv

  bottom_task=chunks(chunk)%chunk_neighbours(chunk_bottom) - 1

  CALL MPI_ISEND(bottom_snd_buffer,total_size,MPI_DOUBLE_PRECISION,bottom_task,tag_send &
                ,MPI_COMM_WORLD,req_send,err)

  CALL MPI_IRECV(bottom_rcv_buffer,total_size,MPI_DOUBLE_PRECISION,bottom_task,tag_recv &
                ,MPI_COMM_WORLD,req_recv,err)

END SUBROUTINE tea_send_recv_message_bottom

SUBROUTINE tea_unpack_bottom(fields, chunk, depth,                        &
                             bottom_rcv_buffer,                              &
                             bottom_top_offset)

  USE pack_kernel_module

  IMPLICIT NONE

  INTEGER         :: fields(:), chunk, depth, bottom_top_offset(:)
  REAL(KIND=8)    :: bottom_rcv_buffer(:)

  IF(fields(FIELD_DENSITY).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%density,                 &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_DENSITY))
    ENDIF
  ENDIF

  IF(fields(FIELD_ENERGY0).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%energy0,                  &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_ENERGY0))
    ENDIF
  ENDIF

  IF(fields(FIELD_ENERGY1).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%energy1,                  &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_ENERGY1))
    ENDIF
  ENDIF
  IF(fields(FIELD_P).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%vector_p,                  &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_P))
    ENDIF
  ENDIF
  IF(fields(FIELD_U).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%u,                  &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_U))
    ENDIF
  ENDIF
  IF(fields(FIELD_SD).EQ.1) THEN
    IF(use_fortran_kernels) THEN
      CALL tea_unpack_message_bottom(chunks(chunk)%field%x_min,                    &
                                        chunks(chunk)%field%x_max,                    &
                                        chunks(chunk)%field%y_min,                    &
                                        chunks(chunk)%field%y_max,                    &
                                        chunks(chunk)%field%vector_sd,                  &
                                        chunks(chunk)%bottom_rcv_buffer,              &
                                        CELL_DATA,VERTEX_DATA,X_FACE_DATA,Y_FACE_DATA,&
                                        depth, CELL_DATA,                             &
                                        bottom_top_offset(FIELD_SD))
    ENDIF
  ENDIF

END SUBROUTINE tea_unpack_bottom

SUBROUTINE tea_sum(value)

  ! Only sums to the master

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: total

  INTEGER :: err

  total=value

  CALL MPI_REDUCE(value,total,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,err)

  value=total

END SUBROUTINE tea_sum

SUBROUTINE tea_allsum(value)

  ! Global reduction for CG solver

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: total

  INTEGER :: err

  total=value

  CALL MPI_ALLREDUCE(value,total,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,err)

  value=total

END SUBROUTINE tea_allsum

SUBROUTINE tea_min(value)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: minimum

  INTEGER :: err

  minimum=value

  CALL MPI_ALLREDUCE(value,minimum,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,err)

  value=minimum

END SUBROUTINE tea_min

SUBROUTINE tea_max(value)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: maximum

  INTEGER :: err

  maximum=value

  CALL MPI_ALLREDUCE(value,maximum,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,err)

  value=maximum

END SUBROUTINE tea_max

SUBROUTINE tea_allgather(value,values)

  IMPLICIT NONE

  REAL(KIND=8) :: value

  REAL(KIND=8) :: values(parallel%max_task)

  INTEGER :: err

  values(1)=value ! Just to ensure it will work in serial

  CALL MPI_ALLGATHER(value,1,MPI_DOUBLE_PRECISION,values,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,err)

END SUBROUTINE tea_allgather

SUBROUTINE tea_check_error(error)

  IMPLICIT NONE

  INTEGER :: error

  INTEGER :: maximum

  INTEGER :: err

  maximum=error

  CALL MPI_ALLREDUCE(error,maximum,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,err)

  error=maximum

END SUBROUTINE tea_check_error


END MODULE tea_module
