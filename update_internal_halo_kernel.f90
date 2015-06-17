
MODULE update_internal_halo_kernel_module

  ! These need to be kept consistent with the data module to avoid use statement
  INTEGER,private,PARAMETER :: CHUNK_LEFT   =1    &
                            ,CHUNK_RIGHT  =2    &
                            ,CHUNK_BOTTOM =3    &
                            ,CHUNK_TOP    =4    &
                            ,EXTERNAL_FACE=-1

  INTEGER,private,PARAMETER :: FIELD_DENSITY    = 1         &
                            ,FIELD_ENERGY0    = 2         &
                            ,FIELD_ENERGY1    = 3         &
                            ,FIELD_U          = 4         &
                            ,FIELD_P          = 5         &
                            ,FIELD_SD         = 6         &
                            ,FIELD_R          = 7         &
                            ,NUM_FIELDS       = 7

CONTAINS

  SUBROUTINE update_internal_halo_left_right_kernel(                                &
                        x_min,x_max,y_min,y_max,                                    &
                        density,                                                    &
                        energy0,                                                    &
                        energy1,                                                    &
                        u,                                                          &
                        p,                                                          &
                        sd,                                                         &
                        x_min_right,x_max_right,y_min_right,y_max_right,            &
                        density_right,                                              &
                        energy0_right,                                              &
                        energy1_right,                                              &
                        u_right,                                                    &
                        p_right,                                                    &
                        sd_right,                                                   &
                        halo_exchange_depth,                                        &
                        fields,                                                     &
                        depth                                                       )
  IMPLICIT NONE

  INTEGER :: halo_exchange_depth
  INTEGER :: fields(NUM_FIELDS),depth

  INTEGER :: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: density,energy0,energy1, u, sd, p

  INTEGER :: x_min_right,x_max_right,y_min_right,y_max_right
  REAL(KIND=8), DIMENSION(x_min_right-halo_exchange_depth:x_max_right+halo_exchange_depth,y_min_right-halo_exchange_depth:y_max_right+halo_exchange_depth) :: density_right,energy0_right,energy1_right, u_right, sd_right, p_right

  END SUBROUTINE

END MODULE update_internal_halo_kernel_module

