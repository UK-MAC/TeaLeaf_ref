MODULE revert_module

CONTAINS

SUBROUTINE revert()

  USE clover_module
  USE revert_kernel_module

  IMPLICIT NONE

  INTEGER :: c

  DO c=1,number_of_chunks

    IF(chunks(c)%task.EQ.parallel%task) THEN

      IF(use_fortran_kernels)THEN
        CALL revert_kernel(chunks(c)%field%x_min,   &
                         chunks(c)%field%x_max,     &
                         chunks(c)%field%y_min,     &
                         chunks(c)%field%y_max,     &
                         chunks(c)%field%density0,  &
                         chunks(c)%field%density1,  &
                         chunks(c)%field%energy0,   &
                         chunks(c)%field%energy1    )
      ELSEIF(use_C_kernels)THEN
        CALL revert_kernel_c(chunks(c)%field%x_min, &
                         chunks(c)%field%x_max,     &
                         chunks(c)%field%y_min,     &
                         chunks(c)%field%y_max,     &
                         chunks(c)%field%density0,  &
                         chunks(c)%field%density1,  &
                         chunks(c)%field%energy0,   &
                         chunks(c)%field%energy1    )
      ENDIF

    ENDIF

  ENDDO

END SUBROUTINE revert

END MODULE revert_module
