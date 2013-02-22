MODULE ideal_gas_module

CONTAINS

SUBROUTINE ideal_gas(chunk,predict)

  USE clover_module
  USE ideal_gas_kernel_module

  IMPLICIT NONE

  INTEGER :: chunk

  LOGICAl :: predict

  IF(chunks(chunk)%task .EQ. parallel%task) THEN

    IF(.NOT.predict) THEN
      IF(use_fortran_kernels)THEN
        CALL ideal_gas_kernel(chunks(chunk)%field%x_min,    &
                            chunks(chunk)%field%x_max,      &
                            chunks(chunk)%field%y_min,      &
                            chunks(chunk)%field%y_max,      &
                            chunks(chunk)%field%density0,   &
                            chunks(chunk)%field%energy0,    &
                            chunks(chunk)%field%pressure,   &
                            chunks(chunk)%field%soundspeed  )
      ELSEIF(use_C_kernels)THEN
        CALL ideal_gas_kernel_c(chunks(chunk)%field%x_min,  &
                            chunks(chunk)%field%x_max,      &
                            chunks(chunk)%field%y_min,      &
                            chunks(chunk)%field%y_max,      &
                            chunks(chunk)%field%density0,   &
                            chunks(chunk)%field%energy0,    &
                            chunks(chunk)%field%pressure,   &
                            chunks(chunk)%field%soundspeed  )
      ENDIF
    ELSE
      IF(use_fortran_kernels)THEN
        CALL ideal_gas_kernel(chunks(chunk)%field%x_min,    &
                            chunks(chunk)%field%x_max,      &
                            chunks(chunk)%field%y_min,      &
                            chunks(chunk)%field%y_max,      &
                            chunks(chunk)%field%density1,   &
                            chunks(chunk)%field%energy1,    &
                            chunks(chunk)%field%pressure,   &
                            chunks(chunk)%field%soundspeed  )
      ELSEIF(use_C_kernels)THEN
        CALL ideal_gas_kernel_c(chunks(chunk)%field%x_min,  &
                            chunks(chunk)%field%x_max,      &
                            chunks(chunk)%field%y_min,      &
                            chunks(chunk)%field%y_max,      &
                            chunks(chunk)%field%density1,   &
                            chunks(chunk)%field%energy1,    &
                            chunks(chunk)%field%pressure,   &
                            chunks(chunk)%field%soundspeed  )
      ENDIF
    ENDIF

  ENDIF

END SUBROUTINE ideal_gas

END MODULE ideal_gas_module
