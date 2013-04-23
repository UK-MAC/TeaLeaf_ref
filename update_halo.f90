MODULE update_halo_module

CONTAINS

SUBROUTINE update_halo(fields,depth)

  USE clover_module
  USE update_halo_kernel_module

  IMPLICIT NONE

  INTEGER :: c,fields(NUM_FIELDS),depth

  CALL clover_exchange(fields,depth)

  DO c=1,number_of_chunks

    IF(chunks(c)%task.EQ.parallel%task) THEN

      CALL update_halo_kernel(chunks(c)%field%x_min,          &
                              chunks(c)%field%x_max,          &
                              chunks(c)%field%y_min,          &
                              chunks(c)%field%y_max,          &
                              chunks(c)%field%left,           &
                              chunks(c)%field%bottom,         &
                              chunks(c)%field%right,          &
                              chunks(c)%field%top,            &
                              chunks(c)%field%left_boundary,  &
                              chunks(c)%field%bottom_boundary,&
                              chunks(c)%field%right_boundary, &
                              chunks(c)%field%top_boundary,   &
                              chunks(c)%chunk_neighbours,     &
                              chunks(c)%field%density0,       &
                              chunks(c)%field%energy0,        &
                              chunks(c)%field%pressure,       &
                              chunks(c)%field%viscosity,      &
                              chunks(c)%field%soundspeed,     &
                              chunks(c)%field%density1,       &
                              chunks(c)%field%energy1,        &
                              chunks(c)%field%xvel0,          &
                              chunks(c)%field%yvel0,          &
                              chunks(c)%field%xvel1,          &
                              chunks(c)%field%yvel1,          &
                              chunks(c)%field%vol_flux_x,     &
                              chunks(c)%field%vol_flux_y,     &
                              chunks(c)%field%mass_flux_x,    &
                              chunks(c)%field%mass_flux_y,    &
                              chunks(c)%field%u,              &
                              fields,                         &
                              depth                           )

    ENDIF

  ENDDO

END SUBROUTINE update_halo

END MODULE update_halo_module
