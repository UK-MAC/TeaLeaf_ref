SUBROUTINE start

  USE clover_module
  USE parse_module
  USE update_halo_module
  USE ideal_gas_module

  IMPLICIT NONE

  INTEGER :: c,j,k

  INTEGER :: x_cells,y_cells
  INTEGER, ALLOCATABLE :: right(:),left(:),top(:),bottom(:)

  INTEGER :: fields(NUM_FIELDS)

  IF(parallel%boss)THEN
     WRITE(g_out,*) 'Setting up initial geometry'
     WRITE(g_out,*)
  ENDIF

  time  = 0.0
  step  = 0
  dtold = dtinit
  dt    = dtinit

  CALL clover_barrier

  CALL clover_get_num_chunks(number_of_chunks)

  ALLOCATE(chunks(1:number_of_chunks))
  ALLOCATE(left(1:number_of_chunks))
  ALLOCATE(right(1:number_of_chunks))
  ALLOCATE(bottom(1:number_of_chunks))
  ALLOCATE(top(1:number_of_chunks))

  CALL clover_decompose(grid%x_cells,grid%y_cells,left,right,bottom,top)

  DO c=1,number_of_chunks
      
    ! Needs changing so there can be more than 1 chunk per task
    chunks(c)%task = c-1

    x_cells = right(c) -left(c)  +1
    y_cells = top(c)   -bottom(c)+1
      
    IF(chunks(c)%task.EQ.parallel%task)THEN
      CALL build_field(c,x_cells,y_cells)
    ENDIF
    chunks(c)%field%left    = left(c)
    chunks(c)%field%bottom  = bottom(c)
    chunks(c)%field%right   = right(c)
    chunks(c)%field%top     = top(c)
    chunks(c)%field%left_boundary   = 1
    chunks(c)%field%bottom_boundary = 1
    chunks(c)%field%right_boundary  = grid%x_cells
    chunks(c)%field%top_boundary    = grid%y_cells
    chunks(c)%field%x_min = 1
    chunks(c)%field%y_min = 1
    chunks(c)%field%x_max = right(c)-left(c)+1
    chunks(c)%field%y_max = top(c)-bottom(c)+1

  ENDDO

  DEALLOCATE(left,right,bottom,top)

  CALL clover_barrier

  DO c=1,number_of_chunks
    IF(chunks(c)%task.EQ.parallel%task)THEN
      CALL clover_allocate_buffers(c)
    ENDIF
  ENDDO

  DO c=1,number_of_chunks
    IF(chunks(c)%task.EQ.parallel%task)THEN
      CALL initialise_chunk(c)
    ENDIF
  ENDDO

  IF(parallel%boss)THEN
     WRITE(g_out,*) 'Generating chunks'
  ENDIF

  DO c=1,number_of_chunks
    IF(chunks(c)%task.EQ.parallel%task)THEN
      CALL generate_chunk(c)
    ENDIF
  ENDDO

  advect_x=.TRUE.

  CALL clover_barrier

  DO c = 1, number_of_chunks
    CALL ideal_gas(c,.FALSE.)
  END DO

  ! Prime all halo data for the first step
  fields=0
  fields(FIELD_DENSITY0)=1
  fields(FIELD_ENERGY0)=1
  fields(FIELD_PRESSURE)=1
  fields(FIELD_VISCOSITY)=1
  fields(FIELD_DENSITY1)=1
  fields(FIELD_ENERGY1)=1
  fields(FIELD_XVEL0)=1
  fields(FIELD_YVEL0)=1
  fields(FIELD_XVEL1)=1
  fields(FIELD_YVEL1)=1

  CALL update_halo(fields,2)

  IF(parallel%boss)THEN
     WRITE(g_out,*)
     WRITE(g_out,*) 'Problem initialised and generated'
  ENDIF

  CALL field_summary()

  IF(visit_frequency.NE.0) CALL visit()

  CALL clover_barrier

END SUBROUTINE start
