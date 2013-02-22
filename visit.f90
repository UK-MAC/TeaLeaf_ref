SUBROUTINE visit

  USE clover_module
  USE update_halo_module
  USE viscosity_module
  USE ideal_gas_module

  IMPLICIT NONE

  INTEGER :: i,j,k,c,l,err,get_unit,u,dummy
  INTEGER :: nxc,nyc,nxv,nyv,nblocks
  REAL(KIND=8)    :: t0,vx,vy,temp_var

  CHARACTER(len=80)           :: name
  CHARACTER(len=10)           :: chunk_name,step_name
  CHARACTER(len=90)           :: filename

  LOGICAL, SAVE :: first_call=.TRUE.

  INTEGER :: fields(NUM_FIELDS)

  name = 'clover'

  IF(first_call) THEN

    nblocks=number_of_chunks
    filename = "clover.visit"
    u=get_unit(dummy)
    OPEN(UNIT=u,FILE=filename,STATUS='UNKNOWN',IOSTAT=err)
    WRITE(u,'(a,i5)')'!NBLOCKS ',nblocks
    CLOSE(u)

    first_call=.FALSE.

  ENDIF

  DO c=1,number_of_chunks
    CALL ideal_gas(c,.FALSE.)
  ENDDO

  fields=0
  fields(FIELD_PRESSURE)=1
  fields(FIELD_XVEL0)=1
  fields(FIELD_YVEL0)=1
  CALL update_halo(fields,1)

  CALL viscosity()

  IF ( parallel%boss ) THEN

    filename = "clover.visit"
    u=get_unit(dummy)
    OPEN(UNIT=u,FILE=filename,STATUS='UNKNOWN',POSITION='APPEND',IOSTAT=err)

    DO c = 1, number_of_chunks
      WRITE(chunk_name, '(i6)') c+100000
      chunk_name(1:1) = "."
      WRITE(step_name, '(i6)') step+100000
      step_name(1:1) = "."
      filename = trim(trim(name) //trim(chunk_name)//trim(step_name))//".vtk"
      WRITE(u,'(a)')TRIM(filename)
    ENDDO
    CLOSE(u)

  ENDIF

  DO c = 1, number_of_chunks
    IF(chunks(c)%task.EQ.parallel%task) THEN
      nxc=chunks(c)%field%x_max-chunks(c)%field%x_min+1
      nyc=chunks(c)%field%y_max-chunks(c)%field%y_min+1
      nxv=nxc+1
      nyv=nyc+1
      WRITE(chunk_name, '(i6)') c+100000
      chunk_name(1:1) = "."
      WRITE(step_name, '(i6)') step+100000
      step_name(1:1) = "."
      filename = trim(trim(name) //trim(chunk_name)//trim(step_name))//".vtk"
      u=get_unit(dummy)
      OPEN(UNIT=u,FILE=filename,STATUS='UNKNOWN',IOSTAT=err)
      WRITE(u,'(a)')'# vtk DataFile Version 3.0'
      WRITE(u,'(a)')'vtk output'
      WRITE(u,'(a)')'ASCII'
      WRITE(u,'(a)')'DATASET RECTILINEAR_GRID'
      WRITE(u,'(a,2i12,a)')'DIMENSIONS',nxv,nyv,' 1'
      WRITE(u,'(a,i5,a)')'X_COORDINATES ',nxv,' double'
      DO j=chunks(c)%field%x_min,chunks(c)%field%x_max+1
        WRITE(u,'(e12.4)')chunks(c)%field%vertexx(j)
      ENDDO
      WRITE(u,'(a,i5,a)')'Y_COORDINATES ',nyv,' double'
      DO k=chunks(c)%field%y_min,chunks(c)%field%y_max+1
        WRITE(u,'(e12.4)')chunks(c)%field%vertexy(k)
      ENDDO
      WRITE(u,'(a)')'Z_COORDINATES 1 double'
      WRITE(u,'(a)')'0'
      WRITE(u,'(a,i20)')'CELL_DATA ',nxc*nyc
      WRITE(u,'(a)')'FIELD FieldData 4'
      WRITE(u,'(a,i20,a)')'density 1 ',nxc*nyc,' double'
      DO k=chunks(c)%field%y_min,chunks(c)%field%y_max
        WRITE(u,'(e12.4)')(chunks(c)%field%density0(j,k),j=chunks(c)%field%x_min,chunks(c)%field%x_max)
      ENDDO
      WRITE(u,'(a,i20,a)')'energy 1 ',nxc*nyc,' double'
      DO k=chunks(c)%field%y_min,chunks(c)%field%y_max
        WRITE(u,'(e12.4)')(chunks(c)%field%energy0(j,k),j=chunks(c)%field%x_min,chunks(c)%field%x_max)
      ENDDO
      WRITE(u,'(a,i20,a)')'pressure 1 ',nxc*nyc,' double'
      DO k=chunks(c)%field%y_min,chunks(c)%field%y_max
        WRITE(u,'(e12.4)')(chunks(c)%field%pressure(j,k),j=chunks(c)%field%x_min,chunks(c)%field%x_max)
      ENDDO
      WRITE(u,'(a,i20,a)')'viscosity 1 ',nxc*nyc,' double'
      DO k=chunks(c)%field%y_min,chunks(c)%field%y_max
        DO j=chunks(c)%field%x_min,chunks(c)%field%x_max
          temp_var=0.0
          IF(chunks(c)%field%viscosity(j,k).GT.0.00000001) temp_var=chunks(c)%field%viscosity(j,k)
          WRITE(u,'(e12.4)') temp_var
        ENDDO
      ENDDO
      WRITE(u,'(a,i20)')'POINT_DATA ',nxv*nyv
      WRITE(u,'(a)')'FIELD FieldData 2'
      WRITE(u,'(a,i20,a)')'x_vel 1 ',nxv*nyv,' double'
      DO k=chunks(c)%field%y_min,chunks(c)%field%y_max+1
        DO j=chunks(c)%field%x_min,chunks(c)%field%x_max+1
          temp_var=0.0
          IF(ABS(chunks(c)%field%xvel0(j,k)).GT.0.00000001) temp_var=chunks(c)%field%xvel0(j,k)
          WRITE(u,'(e12.4)') temp_var
        ENDDO
      ENDDO
      WRITE(u,'(a,i20,a)')'y_vel 1 ',nxv*nyv,' double'
      DO k=chunks(c)%field%y_min,chunks(c)%field%y_max+1
        DO j=chunks(c)%field%x_min,chunks(c)%field%x_max+1
          temp_var=0.0
          IF(ABS(chunks(c)%field%yvel0(j,k)).GT.0.00000001) temp_var=chunks(c)%field%yvel0(j,k)
          WRITE(u,'(e12.4)') temp_var
        ENDDO
      ENDDO
      CLOSE(u)
    ENDIF
  ENDDO

END SUBROUTINE visit
