SUBROUTINE build_field(chunk,x_cells,y_cells)

   USE clover_module

   IMPLICIT NONE

   INTEGER :: chunk,x_cells,y_cells

   chunks(chunk)%field%x_min=1
   chunks(chunk)%field%y_min=1

   chunks(chunk)%field%x_max=x_cells
   chunks(chunk)%field%y_max=y_cells

   ALLOCATE(chunks(chunk)%field%density0  (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                   chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%density1  (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                   chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%energy0   (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                   chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%energy1   (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                   chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%pressure  (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                   chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%viscosity (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                   chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%soundspeed(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                   chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))

   ALLOCATE(chunks(chunk)%field%xvel0(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                      chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%xvel1(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                      chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%yvel0(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                      chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%yvel1(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                      chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))

   ALLOCATE(chunks(chunk)%field%vol_flux_x (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%mass_flux_x(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%vol_flux_y (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%mass_flux_y(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))

   ALLOCATE(chunks(chunk)%field%u         (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                   chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))

   ALLOCATE(chunks(chunk)%field%work_array1(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%work_array2(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%work_array3(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%work_array4(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%work_array5(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%work_array6(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%work_array7(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                            chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))

   ALLOCATE(chunks(chunk)%field%cellx   (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2))
   ALLOCATE(chunks(chunk)%field%celly   (chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%vertexx (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3))
   ALLOCATE(chunks(chunk)%field%vertexy (chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%celldx  (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2))
   ALLOCATE(chunks(chunk)%field%celldy  (chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%vertexdx(chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3))
   ALLOCATE(chunks(chunk)%field%vertexdy(chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))
   ALLOCATE(chunks(chunk)%field%volume  (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                                         chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%xarea   (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+3, &
                                         chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+2))
   ALLOCATE(chunks(chunk)%field%yarea   (chunks(chunk)%field%x_min-2:chunks(chunk)%field%x_max+2, &
                                         chunks(chunk)%field%y_min-2:chunks(chunk)%field%y_max+3))

END SUBROUTINE build_field
