SUBROUTINE initialise_chunk(chunk)

  USE clover_module
  USE initialise_chunk_kernel_module

  IMPLICIT NONE

  INTEGER :: chunk

  REAL(KIND=8) :: xmin,ymin,dx,dy

  dx=(grid%xmax-grid%xmin)/float(grid%x_cells)
  dy=(grid%ymax-grid%ymin)/float(grid%y_cells)

  xmin=grid%xmin+dx*float(chunks(chunk)%field%left-1)

  ymin=grid%ymin+dy*float(chunks(chunk)%field%bottom-1)

  CALL initialise_chunk_kernel(chunks(chunk)%field%x_min,    &
                               chunks(chunk)%field%x_max,    &
                               chunks(chunk)%field%y_min,    &
                               chunks(chunk)%field%y_max,    &
                               xmin,ymin,dx,dy,              &
                               chunks(chunk)%field%vertexx,  &
                               chunks(chunk)%field%vertexdx, &
                               chunks(chunk)%field%vertexy,  &
                               chunks(chunk)%field%vertexdy, &
                               chunks(chunk)%field%cellx,    &
                               chunks(chunk)%field%celldx,   &
                               chunks(chunk)%field%celly,    &
                               chunks(chunk)%field%celldy,   &
                               chunks(chunk)%field%volume,   &
                               chunks(chunk)%field%xarea,    &
                               chunks(chunk)%field%yarea     )


END SUBROUTINE initialise_chunk
