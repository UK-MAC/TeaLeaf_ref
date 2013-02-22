MODULE calc_dt_module

CONTAINS

SUBROUTINE calc_dt(chunk,local_dt,local_control,xl_pos,yl_pos,jldt,kldt)

  USE clover_module
  USE calc_dt_kernel_module

  IMPLICIT NONE

  INTEGER          :: chunk
  REAL(KIND=8)     :: local_dt
  CHARACTER(LEN=8) :: local_control
  REAL(KIND=8)     :: xl_pos,yl_pos
  INTEGER          :: jldt,kldt

  INTEGER          :: l_control
  INTEGER          :: small

  local_dt=g_big

  IF(chunks(chunk)%task.NE.parallel%task) RETURN

  small = 0

  CALL calc_dt_kernel(chunks(chunk)%field%x_min,     &
                      chunks(chunk)%field%x_max,     &
                      chunks(chunk)%field%y_min,     &
                      chunks(chunk)%field%y_max,     &
                      g_small,                       &
                      g_big,                         &
                      dtmin,                         &
                      dtc_safe,                      &
                      dtu_safe,                      &
                      dtv_safe,                      &
                      dtdiv_safe,                    &
                      chunks(chunk)%field%xarea,     &
                      chunks(chunk)%field%yarea,     &
                      chunks(chunk)%field%cellx,     &
                      chunks(chunk)%field%celly,     &
                      chunks(chunk)%field%celldx,    &
                      chunks(chunk)%field%celldy,    &
                      chunks(chunk)%field%volume,    &
                      chunks(chunk)%field%density0,  &
                      chunks(chunk)%field%energy0,   &
                      chunks(chunk)%field%pressure,  &
                      chunks(chunk)%field%viscosity, &
                      chunks(chunk)%field%soundspeed,&
                      chunks(chunk)%field%xvel0,     &
                      chunks(chunk)%field%yvel0,     &
                      chunks(chunk)%field%work_array1,&
                      local_dt,                      &
                      l_control,                     &
                      xl_pos,                        &
                      yl_pos,                        &
                      jldt,                          &
                      kldt,                          &
                      small                          )

  IF(l_control.EQ.1) local_control='sound'
  IF(l_control.EQ.2) local_control='xvel'
  IF(l_control.EQ.3) local_control='yvel'
  IF(l_control.EQ.4) local_control='div'

END SUBROUTINE calc_dt

END MODULE calc_dt_module
