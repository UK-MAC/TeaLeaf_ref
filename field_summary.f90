SUBROUTINE field_summary()

  USE clover_module
  USE ideal_gas_module
  USE field_summary_kernel_module

  IMPLICIT NONE

  REAL(KIND=8) :: vol,mass,ie,ke,press

!$ INTEGER :: OMP_GET_THREAD_NUM

  INTEGER      :: c

  IF(parallel%boss)THEN
    WRITE(g_out,*)
    WRITE(g_out,*) 'Time ',time
    WRITE(g_out,'(a13,7a16)')'           ','Volume','Mass','Density','Pressure','Internal Energy','Kinetic Energy','Total Energy'
  ENDIF

  DO c=1,number_of_chunks
    CALL ideal_gas(c,.FALSE.)
  ENDDO

  DO c=1,number_of_chunks
    IF(chunks(c)%task.EQ.parallel%task) THEN
      CALL field_summary_kernel(chunks(c)%field%x_min,                   &
                                chunks(c)%field%x_max,                   &
                                chunks(c)%field%y_min,                   &
                                chunks(c)%field%y_max,                   &
                                chunks(c)%field%volume,                  &
                                chunks(c)%field%density0,                &
                                chunks(c)%field%energy0,                 &
                                chunks(c)%field%pressure,                &
                                chunks(c)%field%xvel0,                   &
                                chunks(c)%field%yvel0,                   &
                                vol,mass,ie,ke,press                     )
    ENDIF
  ENDDO

  ! For mpi I need a reduction here
  CALL clover_sum(vol)
  CALL clover_sum(mass)
  CALL clover_sum(press)
  CALL clover_sum(ie)
  CALL clover_sum(ke)

  IF(parallel%boss) THEN
!$  IF(OMP_GET_THREAD_NUM().EQ.0) THEN
      WRITE(g_out,'(a6,i7,7e16.4)')' step:',step,vol,mass,mass/vol,press/vol,ie,ke,ie+ke
      WRITE(g_out,*)
!$  ENDIF
   ENDIF

END SUBROUTINE field_summary
