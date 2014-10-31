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

!>  @brief Fortran field summary kernel
!>  @author David Beckingsale, Wayne Gaudin
!>  @details The total mass, internal energy, temperature are calculated

MODULE field_summary_kernel_module

CONTAINS

SUBROUTINE field_summary_kernel(x_min,x_max,y_min,y_max, &
                                volume,                  &
                                density,                 &
                                energy0,                 &
                                u,                       &
                                vol,mass,ie,temp         )

  IMPLICIT NONE

  INTEGER      :: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: volume
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density,energy0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u
  REAL(KIND=8) :: vol,mass,ie,temp

  INTEGER      :: j,k
  REAL(KIND=8) :: cell_vol,cell_mass

  vol=0.0
  mass=0.0
  ie=0.0
  temp=0.0

!$OMP PARALLEL
!$OMP DO PRIVATE(cell_vol,cell_mass) REDUCTION(+ : vol,mass,ie,temp)
  DO k=y_min,y_max
    DO j=x_min,x_max
      cell_vol=volume(j,k)
      cell_mass=cell_vol*density(j,k)
      vol=vol+cell_vol
      mass=mass+cell_mass
      ie=ie+cell_mass*energy0(j,k)
      temp=temp+cell_mass*u(j,k)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE field_summary_kernel

END MODULE field_summary_kernel_module
