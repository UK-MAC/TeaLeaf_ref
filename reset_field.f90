MODULE reset_field_module

CONTAINS

SUBROUTINE reset_field()

  USE clover_module
  USE reset_field_kernel_module

  IMPLICIT NONE

  INTEGER :: c

  DO c=1,number_of_chunks

    IF(chunks(c)%task.EQ.parallel%task) THEN

      IF(use_fortran_kernels)THEN
        CALL reset_field_kernel(chunks(c)%field%x_min,   &
                              chunks(c)%field%x_max,     &
                              chunks(c)%field%y_min,     &
                              chunks(c)%field%y_max,     &
                              chunks(c)%field%density0,  &
                              chunks(c)%field%density1,  &
                              chunks(c)%field%energy0,   &
                              chunks(c)%field%energy1,   &
                              chunks(c)%field%xvel0,     &
                              chunks(c)%field%xvel1,     &
                              chunks(c)%field%yvel0,     &
                              chunks(c)%field%yvel1      )
      ELSEIF(use_C_kernels)THEN
        CALL reset_field_kernel_c(chunks(c)%field%x_min, &
                              chunks(c)%field%x_max,     &
                              chunks(c)%field%y_min,     &
                              chunks(c)%field%y_max,     &
                              chunks(c)%field%density0,  &
                              chunks(c)%field%density1,  &
                              chunks(c)%field%energy0,   &
                              chunks(c)%field%energy1,   &
                              chunks(c)%field%xvel0,     &
                              chunks(c)%field%xvel1,     &
                              chunks(c)%field%yvel0,     &
                              chunks(c)%field%yvel1      )
      ENDIF

    ENDIF

  ENDDO

END SUBROUTINE reset_field

END MODULE reset_field_module
