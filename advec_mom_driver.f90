MODULE advec_mom_driver_module

CONTAINS

SUBROUTINE advec_mom_driver(chunk,which_vel,direction,sweep_number)

  USE clover_module
  USE advec_mom_kernel_mod

  IMPLICIT NONE

  INTEGER :: chunk,which_vel,direction,sweep_number,vector

  IF(chunks(chunk)%task.EQ.parallel%task) THEN

    IF(use_fortran_kernels)THEN
      CALL advec_mom_kernel(chunks(chunk)%field%x_min,            &
                          chunks(chunk)%field%x_max,              &
                          chunks(chunk)%field%y_min,              &
                          chunks(chunk)%field%y_max,              &
                          chunks(chunk)%field%xvel1,              &
                          chunks(chunk)%field%yvel1,              &
                          chunks(chunk)%field%mass_flux_x,        &
                          chunks(chunk)%field%vol_flux_x,         &
                          chunks(chunk)%field%mass_flux_y,        &
                          chunks(chunk)%field%vol_flux_y,         &
                          chunks(chunk)%field%volume,             &
                          chunks(chunk)%field%density1,           &
                          chunks(chunk)%field%work_array1,        &
                          chunks(chunk)%field%work_array2,        &
                          chunks(chunk)%field%work_array3,        &
                          chunks(chunk)%field%work_array4,        &
                          chunks(chunk)%field%work_array5,        &
                          chunks(chunk)%field%work_array6,        &
                          chunks(chunk)%field%work_array7,        &
                          chunks(chunk)%field%celldx,             &
                          chunks(chunk)%field%celldy,             &
                          which_vel,                              &
                          sweep_number,                           &
                          direction,                              &
                          use_vector_loops                        )
    ELSEIF(use_C_kernels)THEN
      IF(use_vector_loops) THEN
        vector=1
      ELSE
        vector=0
      ENDIF
      CALL advec_mom_kernel_c(chunks(chunk)%field%x_min,          &
                          chunks(chunk)%field%x_max,              &
                          chunks(chunk)%field%y_min,              &
                          chunks(chunk)%field%y_max,              &
                          chunks(chunk)%field%xvel1,              &
                          chunks(chunk)%field%yvel1,              &
                          chunks(chunk)%field%mass_flux_x,        &
                          chunks(chunk)%field%vol_flux_x,         &
                          chunks(chunk)%field%mass_flux_y,        &
                          chunks(chunk)%field%vol_flux_y,         &
                          chunks(chunk)%field%volume,             &
                          chunks(chunk)%field%density1,           &
                          chunks(chunk)%field%work_array1,        &
                          chunks(chunk)%field%work_array2,        &
                          chunks(chunk)%field%work_array3,        &
                          chunks(chunk)%field%work_array4,        &
                          chunks(chunk)%field%work_array5,        &
                          chunks(chunk)%field%work_array6,        &
                          chunks(chunk)%field%work_array7,        &
                          chunks(chunk)%field%celldx,             &
                          chunks(chunk)%field%celldy,             &
                          which_vel,                              &
                          sweep_number,                           &
                          direction,                              &
                          vector                                  )
    ENDIF

  ENDIF

END SUBROUTINE advec_mom_driver

END MODULE advec_mom_driver_module
