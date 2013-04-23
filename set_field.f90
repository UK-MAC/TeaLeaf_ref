MODULE set_field_module

CONTAINS

SUBROUTINE set_field()

  USE clover_module
  USE set_field_kernel_module

  IMPLICIT NONE

  INTEGER :: c

  DO c=1,number_of_chunks

    IF(chunks(c)%task.EQ.parallel%task) THEN

        CALL set_field_kernel(chunks(c)%field%x_min,   &
                              chunks(c)%field%x_max,     &
                              chunks(c)%field%y_min,     &
                              chunks(c)%field%y_max,     &
                              chunks(c)%field%density0,  &
                              chunks(c)%field%density1,  &
                              chunks(c)%field%energy0,   &
                              chunks(c)%field%energy1)

    ENDIF

  ENDDO

END SUBROUTINE set_field

END MODULE set_field_module
