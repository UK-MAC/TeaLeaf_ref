SUBROUTINE initialise

  USE clover_module
  USE parse_module
  USE report_module

  IMPLICIT NONE

  INTEGER :: ios
  INTEGER :: get_unit,stat,n,uin,out_unit
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
      WRITE(g_out,'(a15,f8.3)') 'Clover Version ',g_version
      WRITE(g_out,'(a18)') 'MPI Version'
!$    WRITE(g_out,'(a18)') 'OpenMP Version'
      WRITE(g_out,'(a14,i6)') 'Task Count ',parallel%max_task !MPI
!$    WRITE(g_out,'(a15,i5)') 'Thread Count: ',OMP_GET_NUM_THREADS()
      WRITE(g_out,*)
      WRITE(*,*)'Output file tea.out opened. All output will go there.'
!$  ENDIF
  ENDIF
!$OMP END PARALLEL

  CALL clover_barrier

  IF(parallel%boss)THEN
    WRITE(g_out,*) 'Clover will run from the following input:-'
    WRITE(g_out,*)
  ENDIF

  IF(parallel%boss)THEN
    uin=get_unit(dummy)

    OPEN(FILE='tea.in',ACTION='READ',STATUS='OLD',UNIT=uin,IOSTAT=ios)
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

  CALL clover_barrier

  g_in=get_unit(dummy)
  OPEN(FILE='tea.in.tmp',ACTION='READ',STATUS='OLD',UNIT=g_in,IOSTAT=ios)

  IF(ios.NE.0) CALL report_error('initialise','Error opening tea.in.tmp file')

  CALL clover_barrier

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

  CALL clover_barrier

  step=0

  CALL start

  CALL clover_barrier

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
