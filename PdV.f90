MODULE PdV_module

CONTAINS

SUBROUTINE PdV(predict)

  USE clover_module
  USE report_module
  USE PdV_kernel_module
  USE revert_module
  USE update_halo_module
  USE ideal_gas_module

  IMPLICIT NONE

  LOGICAL :: predict

  INTEGER :: prdct

  INTEGER :: c
  INTEGER :: fields(NUM_FIELDS)

  error_condition=0 ! Not used yet due to issue with OpenA reduction

  DO c=1,number_of_chunks

    IF(chunks(c)%task.EQ.parallel%task) THEN

      IF(use_fortran_kernels)THEN
        CALL PdV_kernel(predict,                  &
                      chunks(c)%field%x_min,      &
                      chunks(c)%field%x_max,      &
                      chunks(c)%field%y_min,      &
                      chunks(c)%field%y_max,      &
                      dt,                         &
                      chunks(c)%field%xarea,      &
                      chunks(c)%field%yarea,      &
                      chunks(c)%field%volume ,    &
                      chunks(c)%field%density0,   &
                      chunks(c)%field%density1,   &
                      chunks(c)%field%energy0,    &
                      chunks(c)%field%energy1,    &
                      chunks(c)%field%pressure,   &
                      chunks(c)%field%viscosity,  &
                      chunks(c)%field%xvel0,      &
                      chunks(c)%field%xvel1,      &
                      chunks(c)%field%yvel0,      &
                      chunks(c)%field%yvel1,      &
                      chunks(c)%field%work_array1 )
      ELSEIF(use_C_kernels)THEN

        IF(predict) THEN
          prdct=0
        ELSE
          prdct=1
        ENDIF

        CALL PdV_kernel_c(prdct,                  &
                      chunks(c)%field%x_min,      &
                      chunks(c)%field%x_max,      &
                      chunks(c)%field%y_min,      &
                      chunks(c)%field%y_max,      &
                      dt,                         &
                      chunks(c)%field%xarea,      &
                      chunks(c)%field%yarea,      &
                      chunks(c)%field%volume ,    &
                      chunks(c)%field%density0,   &
                      chunks(c)%field%density1,   &
                      chunks(c)%field%energy0,    &
                      chunks(c)%field%energy1,    &
                      chunks(c)%field%pressure,   &
                      chunks(c)%field%viscosity,  &
                      chunks(c)%field%xvel0,      &
                      chunks(c)%field%xvel1,      &
                      chunks(c)%field%yvel0,      &
                      chunks(c)%field%yvel1,      &
                      chunks(c)%field%work_array1 )
      ENDIF
    ENDIF

  ENDDO

  CALL clover_check_error(error_condition)

  IF(error_condition.EQ.1) THEN
    CALL report_error('PdV','error in PdV')
  ENDIF

  IF(predict)THEN
    DO c=1,number_of_chunks
      CALL ideal_gas(c,.TRUE.)
    ENDDO
    fields=0
    fields(FIELD_PRESSURE)=1
    CALL update_halo(fields,1)
  ENDIF

  IF ( predict ) THEN
    CALL revert()
  ENDIF

END SUBROUTINE PdV

END MODULE PdV_module
