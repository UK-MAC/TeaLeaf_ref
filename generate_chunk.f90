SUBROUTINE generate_chunk(chunk)

  USE clover_module
  USE generate_chunk_kernel_module

  IMPLICIT NONE

  INTEGER         :: chunk

  INTEGER         :: state
  REAL(KIND=8), DIMENSION(number_of_states) :: state_density,state_energy,state_xvel,state_yvel
  REAL(KIND=8), DIMENSION(number_of_states) :: state_xmin,state_xmax,state_ymin,state_ymax,state_radius
  INTEGER,      DIMENSION(number_of_states) :: state_geometry

  DO state=1,number_of_states 
   state_density(state)=states(state)%density
   state_energy(state)=states(state)%energy
   state_xvel(state)=states(state)%xvel
   state_yvel(state)=states(state)%yvel
   state_xmin(state)=states(state)%xmin
   state_xmax(state)=states(state)%xmax
   state_ymin(state)=states(state)%ymin
   state_ymax(state)=states(state)%ymax
   state_radius(state)=states(state)%radius
   state_geometry(state)=states(state)%geometry
  ENDDO

  CALL generate_chunk_kernel(chunks(chunk)%field%x_min,             &
                             chunks(chunk)%field%x_max,             &
                             chunks(chunk)%field%y_min,             &
                             chunks(chunk)%field%y_max,             &
                             chunks(chunk)%field%vertexx,           &
                             chunks(chunk)%field%vertexy,           &
                             chunks(chunk)%field%cellx,             &
                             chunks(chunk)%field%celly,             &
                             chunks(chunk)%field%density0,          &
                             chunks(chunk)%field%energy0,           &
                             chunks(chunk)%field%xvel0,             &
                             chunks(chunk)%field%yvel0,             &
                             number_of_states,                      &
                             state_density,                         &
                             state_energy,                          &
                             state_xvel,                            &
                             state_yvel,                            &
                             state_xmin,                            &
                             state_xmax,                            &
                             state_ymin,                            &
                             state_ymax,                            &
                             state_radius,                          &
                             state_geometry,                        &
                             g_rect,                                &
                             g_circ                                 )

END SUBROUTINE generate_chunk
