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

!>  @brief Generates graphics output files.
!>  @author David Beckingsale, Wayne Gaudin
!>  @details The field data over all mesh chunks is written to a .vtk files and
!>  the .visit file is written that defines the time for each set of vtk files.

SUBROUTINE visit

  USE tea_module
  USE update_halo_module

  IMPLICIT NONE

  INTEGER :: j,k,c,err,get_unit,u,dummy
  INTEGER :: nxc,nyc,nxv,nyv,nblocks
  REAL(KIND=8)    :: temp_var

  CHARACTER(len=80)           :: name
  CHARACTER(len=10)           :: chunk_name,step_name
  CHARACTER(len=90)           :: filename

  LOGICAL, SAVE :: first_call=.TRUE.

  INTEGER :: fields(NUM_FIELDS)

  REAL(KIND=8) :: kernel_time,timer

  name = 'tea'

  IF(first_call) THEN

    nblocks=number_of_chunks
    filename = "tea.visit"
    u=get_unit(dummy)
    OPEN(UNIT=u,FILE=filename,STATUS='UNKNOWN',IOSTAT=err)
    WRITE(u,'(a,i5)')'!NBLOCKS ',nblocks
    CLOSE(u)

    first_call=.FALSE.

  ENDIF

  IF ( parallel%boss ) THEN

    filename = "tea.visit"
    u=get_unit(dummy)
    OPEN(UNIT=u,FILE=filename,STATUS='UNKNOWN',POSITION='APPEND',IOSTAT=err)

    DO c = 1, chunks_per_task
      WRITE(chunk_name, '(i6)') c+100000
      chunk_name(1:1) = "."
      WRITE(step_name, '(i6)') step+100000
      step_name(1:1) = "."
      filename = trim(trim(name) //trim(chunk_name)//trim(step_name))//".vtk"
      WRITE(u,'(a)')TRIM(filename)
    ENDDO
    CLOSE(u)

  ENDIF

  IF(profiler_on) kernel_time=timer()
  DO c = 1, chunks_per_task
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
      WRITE(u,'(a)')'FIELD FieldData 3'
      WRITE(u,'(a,i20,a)')'density 1 ',nxc*nyc,' double'
      DO k=chunks(c)%field%y_min,chunks(c)%field%y_max
        WRITE(u,'(e12.4)')(chunks(c)%field%density(j,k),j=chunks(c)%field%x_min,chunks(c)%field%x_max)
      ENDDO
      WRITE(u,'(a,i20,a)')'energy 1 ',nxc*nyc,' double'
      DO k=chunks(c)%field%y_min,chunks(c)%field%y_max
        WRITE(u,'(e12.4)')(chunks(c)%field%energy0(j,k),j=chunks(c)%field%x_min,chunks(c)%field%x_max)
      ENDDO
      WRITE(u,'(a,i20,a)')'temperature 1 ',nxc*nyc,' double'
      DO k=chunks(c)%field%y_min,chunks(c)%field%y_max
        WRITE(u,'(e12.4)')(chunks(c)%field%u(j,k),j=chunks(c)%field%x_min,chunks(c)%field%x_max)
      ENDDO
      CLOSE(u)
    ENDIF
  ENDDO
  IF(profiler_on) profiler%visit=profiler%visit+(timer()-kernel_time)

END SUBROUTINE visit
