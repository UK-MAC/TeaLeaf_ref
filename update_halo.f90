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

  CALL tea_exchange(fields,depth)

  DO c=1,chunks_per_task

    IF(chunks(c)%task.EQ.parallel%task) THEN
      IF(use_fortran_kernels)THEN
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
                                chunks(c)%field%density,        &
                                chunks(c)%field%energy0,        &
                                chunks(c)%field%energy1,        &
                                chunks(c)%field%u,              &
                                chunks(c)%field%vector_p,       &
                                chunks(c)%field%vector_sd,      &
                                fields,                         &
                                depth                           )
      ENDIF
    ENDIF

  ENDDO

END SUBROUTINE update_halo

END MODULE update_halo_module
