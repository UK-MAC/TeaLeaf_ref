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

!>  @brief Fortran heat conduction kernel
!>  @author Michael Boulton, Wayne Gaudin, Douglas Shanks
!>  @details Implicitly calculates the change in temperature using the PPCG method

MODULE tea_leaf_ppcg_kernel_module

USE definitions_module, only: tl_ppcg_active ! added for ppcg init
USE tea_leaf_common_kernel_module

IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_kernel_ppcg_init_sd(x_min,             &
                                        x_max,             &
                                        y_min,             &
                                        y_max,             &
                                        halo_exchange_depth,             &
                                        r,                 &
                                        kx,                 &
                                        ky,                 &
                                        Di,                 &
                                        sd,                &
                                        z,                &
                                        utemp,                &
                                        rtemp,                &
                                        cp,                &
                                        bfp,                &
                                        Mi,                &
                                        rx, ry,             &
                                        theta,             &
                                        preconditioner_type)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: r, sd, kx, ky , z, Mi, Di, utemp, rtemp
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  REAL(KIND=8) :: theta, theta_r, rx, ry

  INTEGER :: j,k

  theta_r = 1.0_8/theta

!$OMP PARALLEL
  IF (preconditioner_type .NE. TL_PREC_NONE) THEN
!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        sd(j, k) =  z(j, k)*theta_r
        rtemp(j, k) =  r(j, k)
        utemp (j, k) = sd(j, k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ELSE
!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        sd(j, k) =  r(j, k)*theta_r
        rtemp(j, k) =  r(j, k)
        utemp (j, k) = sd(j, k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_ppcg_init_sd

SUBROUTINE tea_leaf_ppcg_init_kernel(x_min,  &
                           x_max,                  &
                           y_min,                  &
                           y_max,                  &
                           halo_exchange_depth,                  &
                           p,                      &
                           r,                      &
                           Mi,                     &
                           z,                      &
                           Kx,                     &
                           Ky,                     &
                           Di,                     &
                           cp,                     &
                           bfp,                    &
                           rx,                     &
                           ry,                     &
                           rro,                    &
                           preconditioner_type,    &
                           ppcg_inner_iters,       &
                           ch_alphas,              &
                           ch_betas,               &
                           theta,                  &
                           solve_time,             &
                           step)

  IMPLICIT NONE

  INTEGER :: preconditioner_type,step
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth,ppcg_inner_iters
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: r, Kx, Ky, Di, z, Mi, p
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  REAL(KIND=8) :: theta,solve_time
  REAL(KIND=8), DIMENSION(:) :: ch_alphas,ch_betas

  INTEGER(KIND=4) :: j,k

  REAL(kind=8) :: rro
  REAL(KIND=8) :: rx, ry

  ! We divide the algorithm up into steps
  ! step = 1, is a CG step
  ! step = 2 or 3, is a partial step of the PPCG algorithm
  
!$OMP PARALLEL 
  IF (step == 1 .or. step ==3) rro = 0.0_8

  IF (step == 1) then
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      p(j, k) = 0.0_8
      z(j, k) = 0.0_8
    ENDDO
  ENDDO
!$OMP END DO
  ELSEIF (step == 3) then

!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      p(j, k) = 0.0_8
    ENDDO
  ENDDO
!$OMP END DO
  ENDIF

  IF (preconditioner_type .NE. TL_PREC_NONE .or. (tl_ppcg_active .and. step == 3)) THEN
  
    ! We don't apply the preconditioner on the final application of the polynomial acceleration   
    IF (step == 1 .or. step == 2) then
      IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
        CALL tea_block_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,&
                             r, z, cp, bfp, Kx, Ky, Di, rx, ry)
      ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
        CALL tea_diag_solve(x_min, x_max, y_min, y_max, halo_exchange_depth, 0,&
                             r, z, Mi)
      ENDIF
   ENDIF

    IF (step == 1 .or. step == 3) then
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            p(j, k) = z(j, k)
        ENDDO
    ENDDO
!$OMP END DO NOWAIT
    ENDIF
  ELSE
    IF (step == 1) then
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            p(j, k) = r(j, k)
        ENDDO
    ENDDO
!$OMP END DO NOWAIT
    ENDIF
  ENDIF
  IF (step == 1 .or. step ==3) then
!$OMP DO REDUCTION(+:rro)
  DO k=y_min,y_max
    DO j=x_min,x_max
      rro = rro + r(j, k)*p(j, k)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_ppcg_init_kernel


SUBROUTINE tea_leaf_kernel_ppcg_inner(x_min,             &
                                      x_max,             &
                                      y_min,             &
                                      y_max,             &
                                      halo_exchange_depth,             &
                                      x_min_bound,      &
                                      x_max_bound,      &
                                      y_min_bound,      &
                                      y_max_bound,      &
                                      alpha,             &
                                      beta,              &
                                      rx, ry,            &
                                      inner_step,       &
                                      u,                 &
                                      r,                 &
                                      Kx,                &
                                      Ky,                &
                                      Di,                &
                                      sd,                &
                                      z,                &
                                      utemp,                &
                                      rtemp,                &
                                      cp,                &
                                      bfp,                &
                                      Mi,                &
                                      preconditioner_type)
  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: u, r, Kx, Ky, sd , z, Mi, Di, utemp, rtemp
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  INTEGER(KIND=4) :: j,k
  REAL(KIND=8), DIMENSION(:) :: alpha, beta
  REAL(KIND=8) :: smvp, rx, ry

  INTEGER(KIND=4) :: x_min_bound, x_max_bound, y_min_bound, y_max_bound, inner_step

!$OMP PARALLEL PRIVATE(smvp)
      IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
!$OMP DO
        DO k=y_min_bound,y_max_bound
          DO j=x_min_bound,x_max_bound
            smvp = Di(j,k)*sd(j, k)                                 &
              - ry*(Ky(j, k+1)*sd(j, k+1) + Ky(j, k)*sd(j, k-1))  &
              - rx*(Kx(j+1, k)*sd(j+1, k) + Kx(j, k)*sd(j-1, k))
!don't change r or u
            rtemp(j, k) = rtemp(j, k) - smvp
          ENDDO
        ENDDO
!$OMP END DO
!$OMP DO
        DO k=y_min_bound,y_max_bound
          DO j=x_min_bound,x_max_bound
            sd(j, k) = alpha(inner_step)*sd(j, k) + beta(inner_step)*Mi(j, k)*rtemp(j, k)
            utemp(j, k) = utemp(j, k) + sd(j, k)
          ENDDO
        ENDDO
!$OMP END DO
      ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
!$OMP DO
        DO k=y_min_bound,y_max_bound
          DO j=x_min_bound,x_max_bound
            smvp = Di(j,k)*sd(j, k)                                 &
              - ry*(Ky(j, k+1)*sd(j, k+1) + Ky(j, k)*sd(j, k-1))  &
              - rx*(Kx(j+1, k)*sd(j+1, k) + Kx(j, k)*sd(j-1, k))
!don't change r or u
            rtemp(j, k) = rtemp(j, k) - smvp
          ENDDO
        ENDDO
!$OMP END DO
        CALL tea_block_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                               rtemp, z, cp, bfp, Kx, Ky, Di, rx, ry)
  
!$OMP DO
        DO k=y_min_bound,y_max_bound
          DO j=x_min_bound,x_max_bound
            sd(j, k) = alpha(inner_step)*sd(j, k) + beta(inner_step)*z(j, k)
            utemp(j, k) = utemp(j, k) + sd(j, k)
          ENDDO
        ENDDO
!$OMP END DO
      ELSE
!$OMP DO
        DO k=y_min_bound,y_max_bound
          DO j=x_min_bound,x_max_bound
            smvp = Di(j,k)*sd(j, k)                                 &
              - ry*(Ky(j, k+1)*sd(j, k+1) + Ky(j, k)*sd(j, k-1))  &
              - rx*(Kx(j+1, k)*sd(j+1, k) + Kx(j, k)*sd(j-1, k))
!don't change r or u
            rtemp(j, k) = rtemp(j, k) - smvp
          ENDDO
        ENDDO
!$OMP END DO
!$OMP DO
        DO k=y_min_bound,y_max_bound
          DO j=x_min_bound,x_max_bound
              sd(j, k) = alpha(inner_step)*sd(j, k) + beta(inner_step)*rtemp(j, k)
              utemp(j, k) = utemp(j, k) + sd(j, k)
          ENDDO
        ENDDO
!$OMP END DO
      ENDIF
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_ppcg_inner

SUBROUTINE tea_leaf_kernel_ppcg_inner_norxy(x_min,             &
                                      x_max,             &
                                      y_min,             &
                                      y_max,             &
                                      halo_exchange_depth,             &
                                      x_min_bound,      &
                                      x_max_bound,      &
                                      y_min_bound,      &
                                      y_max_bound,      &
                                      alpha,             &
                                      beta,              &
                                      inner_step,       &
                                      u,                 &
                                      r,                 &
                                      Kx,                &
                                      Ky,                &
                                      Di,                &
                                      sd,                &
                                      z,                &
                                      utemp,                &
                                      rtemp,                &
                                      cp,                &
                                      bfp,                &
                                      Mi,                &
                                      preconditioner_type)
  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: u, r, Kx, Ky, sd , z, Mi, Di, utemp, rtemp
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  INTEGER(KIND=4) :: j,k
  REAL(KIND=8), DIMENSION(:) :: alpha, beta
  REAL(KIND=8) :: smvp, rx, ry

  INTEGER(KIND=4) :: x_min_bound, x_max_bound, y_min_bound, y_max_bound, inner_step

!$OMP PARALLEL PRIVATE(smvp)
      IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
!$OMP DO
        DO k=y_min_bound,y_max_bound
          DO j=x_min_bound,x_max_bound
            smvp = Di(j,k)*sd(j, k)                                 &
              - (Ky(j, k+1)*sd(j, k+1) + Ky(j, k)*sd(j, k-1))  &
              - (Kx(j+1, k)*sd(j+1, k) + Kx(j, k)*sd(j-1, k))

            rtemp(j, k) = rtemp(j, k) - smvp
          ENDDO
        ENDDO
!$OMP END DO
!$OMP DO
        DO k=y_min_bound,y_max_bound
          DO j=x_min_bound,x_max_bound
            sd(j, k) = alpha(inner_step)*sd(j, k) + beta(inner_step)*Mi(j, k)*rtemp(j, k)
            utemp(j, k) = utemp(j, k) + sd(j, k)
          ENDDO
        ENDDO
!$OMP END DO NOWAIT
      ELSE
!$OMP DO
    DO k=y_min_bound,y_max_bound
        DO j=x_min_bound,x_max_bound
            smvp = Di(j,k)*sd(j, k)                                 &
                - (Ky(j, k+1)*sd(j, k+1) + Ky(j, k)*sd(j, k-1))  &
                - (Kx(j+1, k)*sd(j+1, k) + Kx(j, k)*sd(j-1, k))

            rtemp(j, k) = rtemp(j, k) - smvp

        ENDDO
    ENDDO
!$OMP END DO

      IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
        CALL tea_block_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                               rtemp, z, cp, bfp, Kx, Ky, Di, rx, ry)
!$OMP DO
        DO k=y_min_bound,y_max_bound
            DO j=x_min_bound,x_max_bound
                sd(j, k) = alpha(inner_step)*sd(j, k) + beta(inner_step)*z(j, k)
                utemp(j, k) = utemp(j, k) + sd(j, k)
            ENDDO
        ENDDO
!$OMP END DO NOWAIT
      ELSE
!$OMP DO
        DO k=y_min_bound,y_max_bound
          DO j=x_min_bound,x_max_bound
            sd(j, k) = alpha(inner_step)*sd(j, k) + beta(inner_step)*rtemp(j, k)
            utemp(j, k) = utemp(j, k) + sd(j, k)
          ENDDO
        ENDDO
!$OMP END DO NOWAIT
      ENDIF
    ENDIF
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_ppcg_inner_norxy

SUBROUTINE tea_leaf_ppcg_calc_zrnorm_kernel(x_min, &
                          x_max,             &
                          y_min,             &
                          y_max,             &
                          halo_exchange_depth,             &
                          z, r,               &
                          preconditioner_type,    &
                          norm)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: r, z
  REAL(KIND=8) :: norm
  integer :: j, k

  norm = 0.0_8

!$OMP PARALLEL
  IF (preconditioner_type .NE. TL_PREC_NONE .or. tl_ppcg_active) THEN
!$OMP DO REDUCTION(+:norm)
    DO k=y_min,y_max
        DO j=x_min,x_max
            norm = norm + z(j, k)*r(j, k)
        ENDDO
    ENDDO
!$OMP END DO
  ELSE
!$OMP DO REDUCTION(+:norm)
    DO k=y_min,y_max
        DO j=x_min,x_max
            norm = norm + r(j, k)*r(j, k)
        ENDDO
    ENDDO
!$OMP END DO
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_ppcg_calc_zrnorm_kernel

SUBROUTINE tea_leaf_ppcg_update_z_kernel(x_min, &
                          x_max,             &
                          y_min,             &
                          y_max,             &
                          halo_exchange_depth,             &
                          z, utemp)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: z, utemp
  integer :: j, k

!$OMP PARALLEL
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            z(j, k) = utemp(j, k)
        ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_ppcg_update_z_kernel

SUBROUTINE tea_leaf_ppcg_pupdate_kernel(x_min, &
                          x_max,             &
                          y_min,             &
                          y_max,             &
                          halo_exchange_depth,             &
                          z, p)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: z, p
  integer :: j, k

!$OMP PARALLEL
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            p(j, k) = z(j, k)
        ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_ppcg_pupdate_kernel

! New routine to store the previous residual

SUBROUTINE tea_leaf_ppcg_store_r_kernel(x_min, &
                                        x_max, &
                                        y_min, &
                                        y_max, &
                                        halo_exchange_depth, &
                                        r, r_store )

  IMPLICIT NONE
  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: x_min, x_max, y_min, y_max, halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: r, r_store

!$OMP PARALLEL
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      r_store(j, k) = r(j, k)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_ppcg_store_r_kernel

SUBROUTINE tea_leaf_ppcg_calc_rrn_kernel(x_min, &
                                         x_max, &
                                         y_min, &
                                         y_max, &
                                         halo_exchange_depth, &
                                         r, &
                                         r_store, &
                                         z, &
                                         rrn )

  IMPLICIT NONE
  INTEGER(KIND=4) :: j,k
  INTEGER(KIND=4) :: x_min, x_max, y_min, y_max, halo_exchange_depth
  REAL(KIND=8) :: rrn
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                          y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: r, r_store, z

  rrn = 0.0_8

!$OMP PARALLEL
!$OMP DO REDUCTION (+:rrn)
  DO k=y_min,y_max
    DO j=x_min,x_max
      rrn = rrn + (r(j, k) - r_store(j, k))*z(j, k)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_ppcg_calc_rrn_kernel

END MODULE

