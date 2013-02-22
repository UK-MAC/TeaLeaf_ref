MODULE timestep_module

CONTAINS

SUBROUTINE timestep()

  USE clover_module
  USE report_module
  USE update_halo_module
  USE viscosity_module
  USE calc_dt_module
  USE ideal_gas_module

  IMPLICIT NONE

  INTEGER :: c
  INTEGER :: jldt,kldt

  REAL(KIND=8)    :: dtlp
  REAL(KIND=8)    :: x_pos,y_pos,xl_pos,yl_pos

  CHARACTER(LEN=8) :: dt_control,dtl_control

  INTEGER :: small

  INTEGER :: fields(NUM_FIELDS)

!$ INTEGER :: OMP_GET_THREAD_NUM

  dt    = g_big
  small=0

  DO c = 1, number_of_chunks
    CALL ideal_gas(c,.FALSE.)
  END DO

  fields=0
  fields(FIELD_PRESSURE)=1
  fields(FIELD_ENERGY0)=1
  fields(FIELD_DENSITY0)=1
  fields(FIELD_XVEL0)=1
  fields(FIELD_YVEL0)=1
  CALL update_halo(fields,1)

  CALL viscosity()

  fields=0
  fields(FIELD_VISCOSITY)=1
  CALL update_halo(fields,1)

  DO c = 1, number_of_chunks
    CALL calc_dt(c,dtlp,dtl_control,xl_pos,yl_pos,jldt,kldt)

    IF(dtlp.LE.dt) THEN
      dt=dtlp
      dt_control=dtl_control
      x_pos=xl_pos
      y_pos=yl_pos
      jdt=jldt
      kdt=kldt
    ENDIF
  END DO

  dt = MIN(dt, (dtold * dtrise), dtmax)

  CALL clover_min(dt)

  IF(dt.LT.dtmin) small=1

  IF (parallel%boss) THEN
!$  IF(OMP_GET_THREAD_NUM().EQ.0) THEN
      WRITE(g_out,"(' Step ', i7,' time ', f11.7,' control ',a11,' timestep  ',1pe9.2,i8,',',i8,' x ',1pe9.2,' y ',1pe9.2)") &
                      step,time,dt_control,dt,jdt,kdt,x_pos,y_pos
      WRITE(0,"(' Step ', i7,' time ', f11.7,' control ',a11,' timestep  ',1pe9.2,i8,',',i8,' x ',1pe9.2,' y ',1pe9.2)") &
                      step,time,dt_control,dt,jdt,kdt,x_pos,y_pos
!$  ENDIF
  ENDIF

  IF(small.EQ.1) THEN
    CALL report_error('timestep','small timestep')
  ENDIF

  dtold = dt

END SUBROUTINE timestep

END MODULE timestep_module
