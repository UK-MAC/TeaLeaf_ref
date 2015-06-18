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

!>  @brief Top level initialisation routine
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Checks for the user input and either invokes the input reader or
!>  switches to the internal test problem. It processes the input and strips
!>  comments before writing a final input file.
!>  It then calls the start routine.

SUBROUTINE initialise

  USE tea_module
  USE parse_module
  USE report_module

  IMPLICIT NONE

  INTEGER :: ios
  INTEGER :: get_unit,stat,uin,out_unit
!$ INTEGER :: OMP_GET_THREAD_NUM,OMP_GET_NUM_THREADS
  CHARACTER(LEN=g_len_max) :: ltmp

  IF(parallel%boss)THEN
    g_out=get_unit(dummy)

    OPEN(FILE='tea.out',ACTION='WRITE',UNIT=g_out,IOSTAT=ios)
    IF(ios.NE.0) CALL report_error('initialise','Error opening tea.out file.')

  ELSE
    g_out=6
  ENDIF

!$OMP PARALLEL
  IF(parallel%boss)THEN
!$  IF(OMP_GET_THREAD_NUM().EQ.0) THEN
      WRITE(g_out,*)
      WRITE(g_out,'(a15,f8.3)') 'Tea Version ',g_version
      WRITE(g_out,'(a18)') 'MPI Version'
!$    WRITE(g_out,'(a18)') 'OpenMP Version'
      WRITE(g_out,'(a14,i6)') 'Task Count ',parallel%max_task !MPI
!$    WRITE(g_out,'(a15,i5)') 'Thread Count: ',OMP_GET_NUM_THREADS()
      WRITE(g_out,*)
      WRITE(*,*)'Output file tea.out opened. All output will go there.'
!$  ENDIF
  ENDIF
!$OMP END PARALLEL

  CALL tea_barrier

  IF(parallel%boss)THEN
    WRITE(g_out,*) 'Tea will run from the following input:-'
    WRITE(g_out,*)
  ENDIF

  IF(parallel%boss)THEN
    uin=get_unit(dummy)

    OPEN(FILE='tea.in',ACTION='READ',STATUS='OLD',UNIT=uin,IOSTAT=ios)
    IF(ios.NE.0) THEN
      out_unit=get_unit(dummy)
      OPEN(FILE='tea.in',UNIT=out_unit,STATUS='REPLACE',ACTION='WRITE',IOSTAT=ios)
      WRITE(out_unit,'(A)')'*tea'
      WRITE(out_unit,'(A)')'state 1 density=100.0 energy=0.0001'
      WRITE(out_unit,'(A)')'state 2 density=0.1 energy=25.0 geometry=rectangle xmin=0.0 xmax=1.0 ymin=1.0 ymax=2.0'
      WRITE(out_unit,'(A)')'state 3 density=0.1 energy=0.1 geometry=rectangle xmin=1.0 xmax=6.0 ymin=1.0 ymax=2.0'
      WRITE(out_unit,'(A)')'state 4 density=0.1 energy=0.1 geometry=rectangle xmin=5.0 xmax=6.0 ymin=1.0 ymax=8.0'
      WRITE(out_unit,'(A)')'state 5 density=0.1 energy=0.1 geometry=rectangle xmin=5.0 xmax=10.0 ymin=7.0 ymax=8.0'
      WRITE(out_unit,'(A)')'x_cells=10'
      WRITE(out_unit,'(A)')'y_cells=10'
      WRITE(out_unit,'(A)')'xmin=0.0'
      WRITE(out_unit,'(A)')'ymin=0.0'
      WRITE(out_unit,'(A)')'xmax=10.0'
      WRITE(out_unit,'(A)')'ymax=10.0'
      WRITE(out_unit,'(A)')'initial_timestep=0.004'
      WRITE(out_unit,'(A)')'end_step=10'
      WRITE(out_unit,'(A)')'tl_max_iters=1000'
      WRITE(out_unit,'(A)')' test_problem 1'
      WRITE(out_unit,'(A)')'tl_use_jacobi'
      WRITE(out_unit,'(A)')'tl_eps=1.0e-15'
      WRITE(out_unit,'(A)')'*endtea'
      CLOSE(out_unit)
      uin=get_unit(dummy)
      OPEN(FILE='tea.in',ACTION='READ',STATUS='OLD',UNIT=uin,IOSTAT=ios)
    ENDIF
    IF(ios.NE.0) CALL report_error('initialise','Error opening tea.in')

    out_unit=get_unit(dummy)
    OPEN(FILE='tea.in.tmp',UNIT=out_unit,STATUS='REPLACE',ACTION='WRITE',IOSTAT=ios)
    IF(ios.NE.0) CALL  report_error('initialise','Error opening tea.in.tmp file')
    stat=parse_init(uin,'')
    DO
       stat=parse_getline(-1_4)
       IF(stat.NE.0)EXIT
       WRITE(out_unit,'(A)') line
    ENDDO
    CLOSE(out_unit)
  ENDIF

  CALL tea_barrier

  g_in=get_unit(dummy)
  OPEN(FILE='tea.in.tmp',ACTION='READ',STATUS='OLD',UNIT=g_in,IOSTAT=ios)

  IF(ios.NE.0) CALL report_error('initialise','Error opening tea.in.tmp file')

  CALL tea_barrier

  IF(parallel%boss)THEN
     REWIND(uin)
     DO
        READ(UNIT=uin,IOSTAT=ios,FMT='(a100)') ltmp ! Read in next line.
        IF(ios.NE.0)EXIT
        WRITE(g_out,FMT='(a100)') ltmp
     ENDDO
  ENDIF

  IF(parallel%boss)THEN
     WRITE(g_out,*)
     WRITE(g_out,*) 'Initialising and generating'
     WRITE(g_out,*)
  ENDIF

  CALL read_input()

  CALL tea_barrier

  step=0

  CALL start()

  CALL tea_barrier

  IF(parallel%boss)THEN
     WRITE(g_out,*) 'Starting the calculation'
  ENDIF

  CLOSE(g_in)

END SUBROUTINE initialise

FUNCTION get_unit(dummy)
  INTEGER :: get_unit,dummy

  INTEGER :: u
  LOGICAL :: used

  DO u=7,99
     INQUIRE(UNIT=u,OPENED=used)
     IF(.NOT.used)THEN
        EXIT
     ENDIF
  ENDDO

  get_unit=u

END FUNCTION get_unit
