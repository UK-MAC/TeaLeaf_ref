MODULE accelerate_module

CONTAINS

SUBROUTINE accelerate()

  USE clover_module
  USE accelerate_kernel_module

  IMPLICIT NONE

  INTEGER :: c

  DO c=1,number_of_chunks

    IF(chunks(c)%task.EQ.parallel%task) THEN

      IF(use_fortran_kernels) THEN
        CALL accelerate_kernel(chunks(c)%field%x_min,                &
                             chunks(c)%field%x_max,                  &
                             chunks(c)%field%y_min,                  &
                             chunks(c)%field%y_max,                  &
                             dt,                                     &
                             chunks(c)%field%xarea,                  &
                             chunks(c)%field%yarea,                  &
                             chunks(c)%field%volume,                 &
                             chunks(c)%field%density0,               &
                             chunks(c)%field%pressure,               &
                             chunks(c)%field%viscosity,              &
                             chunks(c)%field%xvel0,                  &
                             chunks(c)%field%yvel0,                  &
                             chunks(c)%field%xvel1,                  &
                             chunks(c)%field%yvel1,                  &
                             chunks(c)%field%work_array1             )
      ELSEIF(use_C_kernels)THEN
        CALL accelerate_kernel_c(chunks(c)%field%x_min,              &
                             chunks(c)%field%x_max,                  &
                             chunks(c)%field%y_min,                  &
                             chunks(c)%field%y_max,                  &
                             dt,                                     &
                             chunks(c)%field%xarea,                  &
                             chunks(c)%field%yarea,                  &
                             chunks(c)%field%volume,                 &
                             chunks(c)%field%density0,               &
                             chunks(c)%field%pressure,               &
                             chunks(c)%field%viscosity,              &
                             chunks(c)%field%xvel0,                  &
                             chunks(c)%field%yvel0,                  &
                             chunks(c)%field%xvel1,                  &
                             chunks(c)%field%yvel1,                  &
                             chunks(c)%field%work_array1             )
      ENDIF

    ENDIF

  ENDDO

END SUBROUTINE accelerate

END MODULE accelerate_module
