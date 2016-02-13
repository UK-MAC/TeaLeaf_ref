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

SUBROUTINE build_field(level)

  USE tea_module

  IMPLICIT NONE

  INTEGER :: level

  INTEGER :: j,k
  INTEGER :: t

  !write(6,*) "build_field:",level,tiles_per_task,size(chunk(level)%tiles)

!$OMP PARALLEL
!$OMP DO
  DO t=1,tiles_per_task
    chunk(level)%tiles(t)%field%x_min=1
    chunk(level)%tiles(t)%field%y_min=1

    chunk(level)%tiles(t)%field%x_max=chunk(level)%tiles(t)%x_cells
    chunk(level)%tiles(t)%field%y_max=chunk(level)%tiles(t)%y_cells

    ! TODO only allocate extra halo on external tiles

    !write(6,*) "build_field:",level,t,chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth, &
    !     chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
    !     chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth, &
    !     chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth
    ALLOCATE(chunk(level)%tiles(t)%field%density &
        (chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth))
    ALLOCATE(chunk(level)%tiles(t)%field%energy0 &
        (chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth))
    ALLOCATE(chunk(level)%tiles(t)%field%energy1 &
        (chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth))
    ALLOCATE(chunk(level)%tiles(t)%field%u       &
        (chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth))
    ALLOCATE(chunk(level)%tiles(t)%field%u0      &
        (chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth))
    ALLOCATE(chunk(level)%tiles(t)%field%vector_r &
        (chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth))
    ALLOCATE(chunk(level)%tiles(t)%field%vector_r_m1 &
        (chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth))
    ALLOCATE(chunk(level)%tiles(t)%field%vector_Mi &
        (chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth))
    ALLOCATE(chunk(level)%tiles(t)%field%vector_w &
        (chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth))
    ALLOCATE(chunk(level)%tiles(t)%field%vector_z &
        (chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth))
    ALLOCATE(chunk(level)%tiles(t)%field%vector_y &
        (chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth))
    ALLOCATE(chunk(level)%tiles(t)%field%vector_r1 &
        (chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth))
    ALLOCATE(chunk(level)%tiles(t)%field%vector_Di &
        (chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth))
    ALLOCATE(chunk(level)%tiles(t)%field%vector_Kx &
        (chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth))
    ALLOCATE(chunk(level)%tiles(t)%field%vector_Ky &
        (chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth))
    ALLOCATE(chunk(level)%tiles(t)%field%vector_p &
        (chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth))
    ALLOCATE(chunk(level)%tiles(t)%field%vector_sd &
        (chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth))
    ALLOCATE(chunk(level)%tiles(t)%field%row_sums &
        (chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth: &
         chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth))

    ALLOCATE(chunk(level)%tiles(t)%field%tri_cp &
        (chunk(level)%tiles(t)%field%x_min:chunk(level)%tiles(t)%field%x_max, &
         chunk(level)%tiles(t)%field%y_min:chunk(level)%tiles(t)%field%y_max))
    ALLOCATE(chunk(level)%tiles(t)%field%tri_bfp &
        (chunk(level)%tiles(t)%field%x_min:chunk(level)%tiles(t)%field%x_max, &
         chunk(level)%tiles(t)%field%y_min:chunk(level)%tiles(t)%field%y_max))

    ALLOCATE(chunk(level)%tiles(t)%field%cellx   (chunk(level)%tiles(t)%field%x_min-2:chunk(level)%tiles(t)%field%x_max+2))
    ALLOCATE(chunk(level)%tiles(t)%field%celly   (chunk(level)%tiles(t)%field%y_min-2:chunk(level)%tiles(t)%field%y_max+2))
    ALLOCATE(chunk(level)%tiles(t)%field%vertexx (chunk(level)%tiles(t)%field%x_min-2:chunk(level)%tiles(t)%field%x_max+3))
    ALLOCATE(chunk(level)%tiles(t)%field%vertexy (chunk(level)%tiles(t)%field%y_min-2:chunk(level)%tiles(t)%field%y_max+3))
    ALLOCATE(chunk(level)%tiles(t)%field%celldx  (chunk(level)%tiles(t)%field%x_min-2:chunk(level)%tiles(t)%field%x_max+2))
    ALLOCATE(chunk(level)%tiles(t)%field%celldy  (chunk(level)%tiles(t)%field%y_min-2:chunk(level)%tiles(t)%field%y_max+2))
    ALLOCATE(chunk(level)%tiles(t)%field%vertexdx(chunk(level)%tiles(t)%field%x_min-2:chunk(level)%tiles(t)%field%x_max+3))
    ALLOCATE(chunk(level)%tiles(t)%field%vertexdy(chunk(level)%tiles(t)%field%y_min-2:chunk(level)%tiles(t)%field%y_max+3))
    ALLOCATE(chunk(level)%tiles(t)%field%volume  (chunk(level)%tiles(t)%field%x_min-2:chunk(level)%tiles(t)%field%x_max+2, &
         chunk(level)%tiles(t)%field%y_min-2:chunk(level)%tiles(t)%field%y_max+2))
    ALLOCATE(chunk(level)%tiles(t)%field%xarea   (chunk(level)%tiles(t)%field%x_min-2:chunk(level)%tiles(t)%field%x_max+3, &
         chunk(level)%tiles(t)%field%y_min-2:chunk(level)%tiles(t)%field%y_max+2))
    ALLOCATE(chunk(level)%tiles(t)%field%yarea   (chunk(level)%tiles(t)%field%x_min-2:chunk(level)%tiles(t)%field%x_max+2, &
         chunk(level)%tiles(t)%field%y_min-2:chunk(level)%tiles(t)%field%y_max+3))

  ! TODO nested parallelism would require nested zero allocation...?

  ! Zeroing isn't strictly neccessary but it ensures physical pages
  ! are allocated. This prevents first touch overheads in the main code
  ! cycle which can skew timings in the first step
  ! Explicit loop limits ensures correct NUMA access, which array syntax does
  ! not
!$OMP PARALLEL
!$OMP DO
  DO k=chunk(level)%tiles(t)%field%y_min-chunk(level)%halo_exchange_depth, &
       chunk(level)%tiles(t)%field%y_max+chunk(level)%halo_exchange_depth
    DO j=chunk(level)%tiles(t)%field%x_min-chunk(level)%halo_exchange_depth, &
         chunk(level)%tiles(t)%field%x_max+chunk(level)%halo_exchange_depth
      chunk(level)%tiles(t)%field%density(j,k)=0.0
      chunk(level)%tiles(t)%field%energy0(j,k)=0.0
      chunk(level)%tiles(t)%field%energy1(j,k)=0.0
      chunk(level)%tiles(t)%field%u(j,k)=0.0
      chunk(level)%tiles(t)%field%u0(j,k)=0.0

      chunk(level)%tiles(t)%field%vector_r(j,k)=0.0
      chunk(level)%tiles(t)%field%vector_r_m1(j,k)=0.0
      chunk(level)%tiles(t)%field%vector_Mi(j,k)=0.0
      chunk(level)%tiles(t)%field%vector_w(j,k)=0.0
      chunk(level)%tiles(t)%field%vector_z(j,k)=0.0
      chunk(level)%tiles(t)%field%vector_y(j,k)=0.0
      chunk(level)%tiles(t)%field%vector_r1(j,k)=0.0
      chunk(level)%tiles(t)%field%vector_Di(j,k)=0.0
      chunk(level)%tiles(t)%field%vector_Kx(j,k)=0.0
      chunk(level)%tiles(t)%field%vector_Ky(j,k)=0.0
      chunk(level)%tiles(t)%field%vector_p(j,k)=0.0
      chunk(level)%tiles(t)%field%vector_sd(j,k)=0.0
      chunk(level)%tiles(t)%field%row_sums(j,k)=0.0
    ENDDO
  ENDDO
!$OMP ENDDO
!$OMP DO
  DO k=chunk(level)%tiles(t)%field%y_min,chunk(level)%tiles(t)%field%y_max
    DO j=chunk(level)%tiles(t)%field%x_min,chunk(level)%tiles(t)%field%x_max
      chunk(level)%tiles(t)%field%tri_cp(j,k)=0.0
      chunk(level)%tiles(t)%field%tri_bfp(j,k)=0.0
    ENDDO
  ENDDO
!$OMP ENDDO
!$OMP DO
  DO k=chunk(level)%tiles(t)%field%y_min-2,chunk(level)%tiles(t)%field%y_max+2
    DO j=chunk(level)%tiles(t)%field%x_min-2,chunk(level)%tiles(t)%field%x_max+2
      chunk(level)%tiles(t)%field%volume(j,k)=0.0
    ENDDO
  ENDDO
!$OMP ENDDO
!$OMP DO
  DO k=chunk(level)%tiles(t)%field%y_min-2,chunk(level)%tiles(t)%field%y_max+2
      DO j=chunk(level)%tiles(t)%field%x_min-2,chunk(level)%tiles(t)%field%x_max+3
          chunk(level)%tiles(t)%field%xarea(j,k)=0.0
      ENDDO
  ENDDO
!$OMP END DO
!$OMP DO
  DO k=chunk(level)%tiles(t)%field%y_min-2,chunk(level)%tiles(t)%field%y_max+3
      DO j=chunk(level)%tiles(t)%field%x_min-2,chunk(level)%tiles(t)%field%x_max+2
          chunk(level)%tiles(t)%field%yarea(j,k)=0.0
      ENDDO
  ENDDO
!$OMP END DO
!$OMP DO
  DO j=chunk(level)%tiles(t)%field%x_min-2,chunk(level)%tiles(t)%field%x_max+2
      chunk(level)%tiles(t)%field%cellx(j)=0.0
      chunk(level)%tiles(t)%field%celldx(j)=0.0
  ENDDO
!$OMP END DO
!$OMP DO
  DO k=chunk(level)%tiles(t)%field%y_min-2,chunk(level)%tiles(t)%field%y_max+2
      chunk(level)%tiles(t)%field%celly(k)=0.0
      chunk(level)%tiles(t)%field%celldy(k)=0.0
  ENDDO
!$OMP END DO
!$OMP DO
  DO j=chunk(level)%tiles(t)%field%x_min-2,chunk(level)%tiles(t)%field%x_max+3
      chunk(level)%tiles(t)%field%vertexx(j)=0.0
      chunk(level)%tiles(t)%field%vertexdx(j)=0.0
  ENDDO
!$OMP END DO
!$OMP DO
  DO k=chunk(level)%tiles(t)%field%y_min-2,chunk(level)%tiles(t)%field%y_max+3
      chunk(level)%tiles(t)%field%vertexy(k)=0.0
      chunk(level)%tiles(t)%field%vertexdy(k)=0.0
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE build_field

