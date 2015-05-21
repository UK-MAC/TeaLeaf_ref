!Crown Copyright 2014 AWE.
!
! This file is part of TeaLeaf.
!
! TeaLeaf is free software: you can redistribute it and/or modify it under 
! the terms of the GNU General Public License as published by the 
! Free Software Foundation, either version 3 of the License, or (at your option) 
! any later version.
!
! TeaLeaf is distributed in the hope that it will be useful, but 
! WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
! FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
! details.
!
! You should have received a copy of the GNU General Public License along with 
! TeaLeaf. If not, see http://www.gnu.org/licenses/.

!>  @brief Driver for the halo updates
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Invokes the kernels for the internal and external halo cells for
!>  the fields specified.

MODULE update_halo_module

CONTAINS

SUBROUTINE update_halo(fields,depth)

  USE tea_module
  USE update_halo_kernel_module

  IMPLICIT NONE

  INTEGER :: c,fields(NUM_FIELDS),depth
  REAL(KIND=8) :: timer,halo_time

  IF (profiler_on) halo_time=timer()
  CALL tea_exchange(fields,depth)
  IF (profiler_on) profiler%halo_exchange = profiler%halo_exchange + (timer() - halo_time)

  IF (profiler_on) halo_time=timer()
  DO c=1,chunks_per_task

  IF(chunks(c)%task.EQ.parallel%task) THEN
    IF(use_fortran_kernels)THEN
      CALL update_halo_kernel(chunks(c)%field%x_min,          &
                              chunks(c)%field%x_max,          &
                              chunks(c)%field%y_min,          &
                              chunks(c)%field%y_max, halo_exchange_depth,          &
                              chunks(c)%chunk_neighbours,     &
                              chunks(c)%field%density,        &
                              chunks(c)%field%energy0,        &
                              chunks(c)%field%energy1,        &
                              chunks(c)%field%u,              &
                              chunks(c)%field%vector_p,       &
                              chunks(c)%field%vector_sd,      &
                              fields,                         &
                              reflective_boundary,            &
                              depth                           )
    ENDIF
  ENDIF

  ENDDO
  IF (profiler_on) profiler%halo_update = profiler%halo_update + (timer() - halo_time)

END SUBROUTINE update_halo

END MODULE update_halo_module
