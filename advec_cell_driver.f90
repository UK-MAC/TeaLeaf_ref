MODULE  advec_cell_driver_module

CONTAINS

SUBROUTINE advec_cell_driver(chunk,sweep_number,dir)

  USE clover_module
  USE advec_cell_kernel_module

  IMPLICIT NONE

  INTEGER :: chunk,sweep_number,dir,vector

  IF(chunks(chunk)%task.EQ.parallel%task) THEN

    IF(use_fortran_kernels)THEN
      CALL advec_cell_kernel(chunks(chunk)%field%x_min,               &
                           chunks(chunk)%field%x_max,                 &
                           chunks(chunk)%field%y_min,                 &
                           chunks(chunk)%field%y_max,                 &
                           dir,                                       &
                           sweep_number,                              &
                           use_vector_loops,                          &
                           chunks(chunk)%field%vertexdx,              &
                           chunks(chunk)%field%vertexdy,              &
                           chunks(chunk)%field%volume,                &
                           chunks(chunk)%field%density1,              &
                           chunks(chunk)%field%energy1,               &
                           chunks(chunk)%field%mass_flux_x,           &
                           chunks(chunk)%field%vol_flux_x,            &
                           chunks(chunk)%field%mass_flux_y,           &
                           chunks(chunk)%field%vol_flux_y,            &
                           chunks(chunk)%field%work_array1,           &
                           chunks(chunk)%field%work_array2,           &
                           chunks(chunk)%field%work_array3,           &
                           chunks(chunk)%field%work_array4,           &
                           chunks(chunk)%field%work_array5,           &
                           chunks(chunk)%field%work_array6,           &
                           chunks(chunk)%field%work_array7            )
    ELSEIF(use_C_kernels)THEN
      IF(use_vector_loops) THEN
        vector=1
      ELSE
        vector=0
      ENDIF
      CALL advec_cell_kernel_c(chunks(chunk)%field%x_min,             &
                           chunks(chunk)%field%x_max,                 &
                           chunks(chunk)%field%y_min,                 &
                           chunks(chunk)%field%y_max,                 &
                           dir,                                       &
                           sweep_number,                              &
                           vector,                                    &
                           chunks(chunk)%field%vertexdx,              &
                           chunks(chunk)%field%vertexdy,              &
                           chunks(chunk)%field%volume,                &
                           chunks(chunk)%field%density1,              &
                           chunks(chunk)%field%energy1,               &
                           chunks(chunk)%field%mass_flux_x,           &
                           chunks(chunk)%field%vol_flux_x,            &
                           chunks(chunk)%field%mass_flux_y,           &
                           chunks(chunk)%field%vol_flux_y,            &
                           chunks(chunk)%field%work_array1,           &
                           chunks(chunk)%field%work_array2,           &
                           chunks(chunk)%field%work_array3,           &
                           chunks(chunk)%field%work_array4,           &
                           chunks(chunk)%field%work_array5,           &
                           chunks(chunk)%field%work_array6,           &
                           chunks(chunk)%field%work_array7            )
    ENDIF

  ENDIF

END SUBROUTINE advec_cell_driver

END MODULE  advec_cell_driver_module

