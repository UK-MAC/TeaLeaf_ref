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

!>  @brief Fortran set field kernel.
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Copies all of the final start of step filed data to the end of
!>  step data.

MODULE set_field_kernel_module

CONTAINS

SUBROUTINE set_field_kernel(x_min,x_max,y_min,y_max,    &
                            energy0,                    &
                            energy1)                    

  IMPLICIT NONE

  INTEGER :: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy1

  INTEGER :: j,k

!$OMP PARALLEL

!$OMP DO
  DO k=y_min,y_max
     DO j=x_min,x_max
        energy1(j,k)=energy0(j,k)
     ENDDO
  ENDDO
!$OMP END DO

!$OMP END PARALLEL

END SUBROUTINE set_field_kernel

END MODULE set_field_kernel_module
