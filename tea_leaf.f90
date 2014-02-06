PROGRAM tea_leaf

  USE clover_module

  IMPLICIT NONE

!$ INTEGER :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM

  CALL clover_init_comms()

!$OMP PARALLEL
  IF(parallel%boss)THEN
!$  IF(OMP_GET_THREAD_NUM().EQ.0) THEN
      WRITE(*,*)
      WRITE(*,'(a15,f8.3)') 'Tea Version ',g_version
      WRITE(*,'(a18)') 'MPI Version'
!$    WRITE(*,'(a18)') 'OpenMP Version'
      WRITE(*,'(a14,i6)') 'Task Count ',parallel%max_task !MPI
!$    WRITE(*,'(a15,i5)') 'Thread Count: ',OMP_GET_NUM_THREADS()
      WRITE(*,*)
      WRITE(0,*)
      WRITE(0,'(a15,f8.3)') 'Tea Version ',g_version
      WRITE(0,'(a18)') 'MPI Version'
!$    WRITE(0,'(a18)') 'OpenMP Version'
      WRITE(0,'(a14,i6)') 'Task Count ',parallel%max_task !MPI
!$    WRITE(0,'(a15,i5)') 'Thread Count: ',OMP_GET_NUM_THREADS()
      WRITE(0,*)
!$  ENDIF
  ENDIF
!$OMP END PARALLEL

  CALL initialise

  CALL hydro
  
  ! Deallocate everything
  
END PROGRAM tea_leaf

