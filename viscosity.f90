MODULE viscosity_module

CONTAINS

SUBROUTINE viscosity()

  USE clover_module
  USE viscosity_kernel_module
  
  IMPLICIT NONE

  INTEGER :: c

  DO c=1,number_of_chunks

    IF(chunks(c)%task.EQ.parallel%task) THEN

      IF(use_fortran_kernels)THEN
        CALL viscosity_kernel(chunks(c)%field%x_min,                   &
                            chunks(c)%field%x_max,                     &
                            chunks(c)%field%y_min,                     &
                            chunks(c)%field%y_max,                     &
                            chunks(c)%field%celldx,                    &
                            chunks(c)%field%celldy,                    &
                            chunks(c)%field%density0,                  &
                            chunks(c)%field%pressure,                  &
                            chunks(c)%field%viscosity,                 &
                            chunks(c)%field%xvel0,                     &
                            chunks(c)%field%yvel0                      )
      ELSEIF(use_C_kernels)THEN
        CALL viscosity_kernel_c(chunks(c)%field%x_min,                 &
                            chunks(c)%field%x_max,                     &
                            chunks(c)%field%y_min,                     &
                            chunks(c)%field%y_max,                     &
                            chunks(c)%field%celldx,                    &
                            chunks(c)%field%celldy,                    &
                            chunks(c)%field%density0,                  &
                            chunks(c)%field%pressure,                  &
                            chunks(c)%field%viscosity,                 &
                            chunks(c)%field%xvel0,                     &
                            chunks(c)%field%yvel0                      )
      ENDIF

    ENDIF

  ENDDO

END SUBROUTINE viscosity

END MODULE viscosity_module
