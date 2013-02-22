MODULE flux_calc_module

CONTAINS

SUBROUTINE flux_calc()

  USE clover_module
  USE flux_calc_kernel_module

  IMPLICIT NONE

  INTEGER :: c

  DO c=1,number_of_chunks

    IF(chunks(c)%task.EQ.parallel%task) THEN

      IF(use_fortran_kernels)THEN
        CALL flux_calc_kernel(chunks(c)%field%x_min,         &
                            chunks(c)%field%x_max,           &
                            chunks(c)%field%y_min,           &
                            chunks(c)%field%y_max,           &
                            dt,                              &
                            chunks(c)%field%xarea,           &
                            chunks(c)%field%yarea,           &
                            chunks(c)%field%xvel0,           &
                            chunks(c)%field%yvel0,           &
                            chunks(c)%field%xvel1,           &
                            chunks(c)%field%yvel1,           &
                            chunks(c)%field%vol_flux_x,      &
                            chunks(c)%field%vol_flux_y       )
      ELSEIF(use_C_kernels)THEN
        CALL flux_calc_kernel_c(chunks(c)%field%x_min,       &
                            chunks(c)%field%x_max,           &
                            chunks(c)%field%y_min,           &
                            chunks(c)%field%y_max,           &
                            dt,                              &
                            chunks(c)%field%xarea,           &
                            chunks(c)%field%yarea,           &
                            chunks(c)%field%xvel0,           &
                            chunks(c)%field%yvel0,           &
                            chunks(c)%field%xvel1,           &
                            chunks(c)%field%yvel1,           &
                            chunks(c)%field%vol_flux_x,      &
                            chunks(c)%field%vol_flux_y       )
      ENDIF

    ENDIF

  ENDDO

END SUBROUTINE flux_calc

END MODULE flux_calc_module
