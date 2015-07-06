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

!>  @brief  Allocates the data for each mesh chunk
!>  @author David Beckingsale, Wayne Gaudin
!>  @details The data fields for the mesh chunk are allocated based on the mesh
!>  size.

SUBROUTINE build_field()

  USE tea_module

  IMPLICIT NONE

  INTEGER :: j,k
  INTEGER :: t

  DO t=1,tiles_per_task
    chunk%tiles(t)%field%x_min=1
    chunk%tiles(t)%field%y_min=1

    chunk%tiles(t)%field%x_max=chunk%tiles(t)%x_cells
    chunk%tiles(t)%field%y_max=chunk%tiles(t)%y_cells

    ALLOCATE(chunk%tiles(t)%field%density  (chunk%tiles(t)%field%x_min-halo_exchange_depth:chunk%tiles(t)%field%x_max+halo_exchange_depth, &
         chunk%tiles(t)%field%y_min-halo_exchange_depth:chunk%tiles(t)%field%y_max+halo_exchange_depth))
    ALLOCATE(chunk%tiles(t)%field%energy0  (chunk%tiles(t)%field%x_min-halo_exchange_depth:chunk%tiles(t)%field%x_max+halo_exchange_depth, &
         chunk%tiles(t)%field%y_min-halo_exchange_depth:chunk%tiles(t)%field%y_max+halo_exchange_depth))
    ALLOCATE(chunk%tiles(t)%field%energy1  (chunk%tiles(t)%field%x_min-halo_exchange_depth:chunk%tiles(t)%field%x_max+halo_exchange_depth, &
         chunk%tiles(t)%field%y_min-halo_exchange_depth:chunk%tiles(t)%field%y_max+halo_exchange_depth))
    ALLOCATE(chunk%tiles(t)%field%u        (chunk%tiles(t)%field%x_min-halo_exchange_depth:chunk%tiles(t)%field%x_max+halo_exchange_depth, &
         chunk%tiles(t)%field%y_min-halo_exchange_depth:chunk%tiles(t)%field%y_max+halo_exchange_depth))
    ALLOCATE(chunk%tiles(t)%field%u0       (chunk%tiles(t)%field%x_min-halo_exchange_depth:chunk%tiles(t)%field%x_max+halo_exchange_depth, &
         chunk%tiles(t)%field%y_min-halo_exchange_depth:chunk%tiles(t)%field%y_max+halo_exchange_depth))
    ALLOCATE(chunk%tiles(t)%field%vector_r (chunk%tiles(t)%field%x_min-halo_exchange_depth:chunk%tiles(t)%field%x_max+halo_exchange_depth, &
         chunk%tiles(t)%field%y_min-halo_exchange_depth:chunk%tiles(t)%field%y_max+halo_exchange_depth))
    ALLOCATE(chunk%tiles(t)%field%vector_Mi(chunk%tiles(t)%field%x_min-halo_exchange_depth:chunk%tiles(t)%field%x_max+halo_exchange_depth, &
         chunk%tiles(t)%field%y_min-halo_exchange_depth:chunk%tiles(t)%field%y_max+halo_exchange_depth))
    ALLOCATE(chunk%tiles(t)%field%vector_w (chunk%tiles(t)%field%x_min-halo_exchange_depth:chunk%tiles(t)%field%x_max+halo_exchange_depth, &
         chunk%tiles(t)%field%y_min-halo_exchange_depth:chunk%tiles(t)%field%y_max+halo_exchange_depth))
    ALLOCATE(chunk%tiles(t)%field%vector_z (chunk%tiles(t)%field%x_min-halo_exchange_depth:chunk%tiles(t)%field%x_max+halo_exchange_depth, &
         chunk%tiles(t)%field%y_min-halo_exchange_depth:chunk%tiles(t)%field%y_max+halo_exchange_depth))
    ALLOCATE(chunk%tiles(t)%field%vector_Kx(chunk%tiles(t)%field%x_min-halo_exchange_depth:chunk%tiles(t)%field%x_max+halo_exchange_depth, &
         chunk%tiles(t)%field%y_min-halo_exchange_depth:chunk%tiles(t)%field%y_max+halo_exchange_depth))
    ALLOCATE(chunk%tiles(t)%field%vector_Ky(chunk%tiles(t)%field%x_min-halo_exchange_depth:chunk%tiles(t)%field%x_max+halo_exchange_depth, &
         chunk%tiles(t)%field%y_min-halo_exchange_depth:chunk%tiles(t)%field%y_max+halo_exchange_depth))
    ALLOCATE(chunk%tiles(t)%field%vector_p (chunk%tiles(t)%field%x_min-halo_exchange_depth:chunk%tiles(t)%field%x_max+halo_exchange_depth, &
         chunk%tiles(t)%field%y_min-halo_exchange_depth:chunk%tiles(t)%field%y_max+halo_exchange_depth))
    ALLOCATE(chunk%tiles(t)%field%vector_sd(chunk%tiles(t)%field%x_min-halo_exchange_depth:chunk%tiles(t)%field%x_max+halo_exchange_depth, &
         chunk%tiles(t)%field%y_min-halo_exchange_depth:chunk%tiles(t)%field%y_max+halo_exchange_depth))

    ALLOCATE(chunk%tiles(t)%field%tri_cp(chunk%tiles(t)%field%x_min:chunk%tiles(t)%field%x_max, &
         chunk%tiles(t)%field%y_min:chunk%tiles(t)%field%y_max))
    ALLOCATE(chunk%tiles(t)%field%tri_bfp(chunk%tiles(t)%field%x_min:chunk%tiles(t)%field%x_max, &
         chunk%tiles(t)%field%y_min:chunk%tiles(t)%field%y_max))

    ALLOCATE(chunk%tiles(t)%field%cellx   (chunk%tiles(t)%field%x_min-2:chunk%tiles(t)%field%x_max+2))
    ALLOCATE(chunk%tiles(t)%field%celly   (chunk%tiles(t)%field%y_min-2:chunk%tiles(t)%field%y_max+2))
    ALLOCATE(chunk%tiles(t)%field%vertexx (chunk%tiles(t)%field%x_min-2:chunk%tiles(t)%field%x_max+3))
    ALLOCATE(chunk%tiles(t)%field%vertexy (chunk%tiles(t)%field%y_min-2:chunk%tiles(t)%field%y_max+3))
    ALLOCATE(chunk%tiles(t)%field%celldx  (chunk%tiles(t)%field%x_min-2:chunk%tiles(t)%field%x_max+2))
    ALLOCATE(chunk%tiles(t)%field%celldy  (chunk%tiles(t)%field%y_min-2:chunk%tiles(t)%field%y_max+2))
    ALLOCATE(chunk%tiles(t)%field%vertexdx(chunk%tiles(t)%field%x_min-2:chunk%tiles(t)%field%x_max+3))
    ALLOCATE(chunk%tiles(t)%field%vertexdy(chunk%tiles(t)%field%y_min-2:chunk%tiles(t)%field%y_max+3))
    ALLOCATE(chunk%tiles(t)%field%volume  (chunk%tiles(t)%field%x_min-2:chunk%tiles(t)%field%x_max+2, &
         chunk%tiles(t)%field%y_min-2:chunk%tiles(t)%field%y_max+2))
    ALLOCATE(chunk%tiles(t)%field%xarea   (chunk%tiles(t)%field%x_min-2:chunk%tiles(t)%field%x_max+3, &
         chunk%tiles(t)%field%y_min-2:chunk%tiles(t)%field%y_max+2))
    ALLOCATE(chunk%tiles(t)%field%yarea   (chunk%tiles(t)%field%x_min-2:chunk%tiles(t)%field%x_max+2, &
         chunk%tiles(t)%field%y_min-2:chunk%tiles(t)%field%y_max+3))
  ENDDO

  ! Zeroing isn't strictly neccessary but it ensures physical pages
  ! are allocated. This prevents first touch overheads in the main code
  ! cycle which can skew timings in the first step
  ! Explicit loop limits ensures correct NUMA access, which array syntax does
  ! not
  DO t=1,tiles_per_task
!$OMP PARALLEL
!$OMP DO
  DO k=chunk%tiles(t)%field%y_min-halo_exchange_depth,chunk%tiles(t)%field%y_max+halo_exchange_depth
    DO j=chunk%tiles(t)%field%x_min-halo_exchange_depth,chunk%tiles(t)%field%x_max+halo_exchange_depth
      chunk%tiles(t)%field%density(j,k)=0.0
      chunk%tiles(t)%field%energy0(j,k)=0.0
      chunk%tiles(t)%field%energy1(j,k)=0.0
      chunk%tiles(t)%field%u(j,k)=0.0
      chunk%tiles(t)%field%u0(j,k)=0.0

      chunk%tiles(t)%field%vector_r(j,k)=0.0
      chunk%tiles(t)%field%vector_Mi(j,k)=0.0
      chunk%tiles(t)%field%vector_w(j,k)=0.0
      chunk%tiles(t)%field%vector_z(j,k)=0.0
      chunk%tiles(t)%field%vector_Kx(j,k)=0.0
      chunk%tiles(t)%field%vector_Ky(j,k)=0.0
      chunk%tiles(t)%field%vector_p(j,k)=0.0
      chunk%tiles(t)%field%vector_sd(j,k)=0.0
    ENDDO
  ENDDO
!$OMP ENDDO
!$OMP DO
  DO k=chunk%tiles(t)%field%y_min,chunk%tiles(t)%field%y_max
    DO j=chunk%tiles(t)%field%x_min,chunk%tiles(t)%field%x_max
      chunk%tiles(t)%field%tri_cp(j,k)=0.0
      chunk%tiles(t)%field%tri_bfp(j,k)=0.0
    ENDDO
  ENDDO
!$OMP ENDDO
!$OMP DO
  DO k=chunk%tiles(t)%field%y_min-2,chunk%tiles(t)%field%y_max+2
    DO j=chunk%tiles(t)%field%x_min-2,chunk%tiles(t)%field%x_max+2
      chunk%tiles(t)%field%volume(j,k)=0.0
    ENDDO
  ENDDO
!$OMP ENDDO
!$OMP DO
  DO k=chunk%tiles(t)%field%y_min-2,chunk%tiles(t)%field%y_max+2
      DO j=chunk%tiles(t)%field%x_min-2,chunk%tiles(t)%field%x_max+3
          chunk%tiles(t)%field%xarea(j,k)=0.0
      ENDDO
  ENDDO
!$OMP END DO
!$OMP DO
  DO k=chunk%tiles(t)%field%y_min-2,chunk%tiles(t)%field%y_max+3
      DO j=chunk%tiles(t)%field%x_min-2,chunk%tiles(t)%field%x_max+2
          chunk%tiles(t)%field%yarea(j,k)=0.0
      ENDDO
  ENDDO
!$OMP END DO
!$OMP DO
  DO j=chunk%tiles(t)%field%x_min-2,chunk%tiles(t)%field%x_max+2
      chunk%tiles(t)%field%cellx(j)=0.0
      chunk%tiles(t)%field%celldx(j)=0.0
  ENDDO
!$OMP END DO
!$OMP DO
  DO k=chunk%tiles(t)%field%y_min-2,chunk%tiles(t)%field%y_max+2
      chunk%tiles(t)%field%celly(k)=0.0
      chunk%tiles(t)%field%celldy(k)=0.0
  ENDDO
!$OMP END DO
!$OMP DO
  DO j=chunk%tiles(t)%field%x_min-2,chunk%tiles(t)%field%x_max+3
      chunk%tiles(t)%field%vertexx(j)=0.0
      chunk%tiles(t)%field%vertexdx(j)=0.0
  ENDDO
!$OMP END DO
!$OMP DO
  DO k=chunk%tiles(t)%field%y_min-2,chunk%tiles(t)%field%y_max+3
      chunk%tiles(t)%field%vertexy(k)=0.0
      chunk%tiles(t)%field%vertexdy(k)=0.0
  ENDDO
!$OMP END DO
!$OMP END PARALLEL
  ENDDO

END SUBROUTINE build_field
