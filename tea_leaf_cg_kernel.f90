!cROWn Copyright 2014 AWE.
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
!>  @author Michael Boulton, Wayne Gaudin
!>  @details Implicitly calculates the change in temperature using CG method

MODULE tea_leaf_cg_kernel_module

  USE definitions_module, only: tl_ppcg_active, tl_ppcg_inner_steps
  USE tea_leaf_common_kernel_module
  USE tea_leaf_ppcg_module, only: tea_leaf_run_ppcg_inner_steps

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_cg_init_kernel(x_min,  &
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
                           level,                  &
                           ppcg_inner_iters,       &
                           ch_alphas,              &
                           ch_betas,               &
                           theta,                  &
                           solve_time,             &
                           step)

  IMPLICIT NONE

  INTEGER :: preconditioner_type,level,step
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth,ppcg_inner_iters
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: r, Kx, Ky, Di, z, Mi, p
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  REAL(KIND=8) :: theta,solve_time
  REAL(KIND=8), DIMENSION(:) :: ch_alphas,ch_betas

  INTEGER(KIND=4) :: j,k

  REAL(kind=8) :: rro
  REAL(KIND=8) :: rx, ry

! step 1 is a CG step, whereas steps 2 and 3 are partial steps from the PPCG algorithm to allow a middle step for the PP application

  !print*,step," r2=",sum(r**2)
  if (step == 1 .or. step ==3) rro = 0.0_8

!$OMP PARALLEL REDUCTION(+:rro)
  if (step == 1) then
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      p(j, k) = 0.0_8
      z(j, k) = 0.0_8
    ENDDO
  ENDDO
!$OMP END DO
  elseif (step == 3) then
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      p(j, k) = 0.0_8
    ENDDO
  ENDDO
!$OMP END DO
  endif

  !IF (preconditioner_type .NE. TL_PREC_NONE .or. tl_ppcg_active) THEN
  IF (preconditioner_type .NE. TL_PREC_NONE .or. (tl_ppcg_active .and. step == 3)) THEN

    if (step == 1 .or. step == 2) then
    IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
      CALL tea_block_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                             r, z, cp, bfp, Kx, Ky, Di, rx, ry)
    ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
      CALL tea_diag_solve(x_min, x_max, y_min, y_max, halo_exchange_depth, 0,           &
                             r, z, Mi)
    ENDIF
    !print*,step," z2=",sum(z**2)
    endif

    !IF (tl_ppcg_active) THEN
    !  CALL tea_leaf_run_ppcg_inner_steps(level, ch_alphas, ch_betas, theta, &
    !      tl_ppcg_inner_steps, solve_time)
    !  ppcg_inner_iters = ppcg_inner_iters + tl_ppcg_inner_steps
    !ENDIF

    if (step == 1 .or. step == 3) then
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            p(j, k) = z(j, k)
        ENDDO
    ENDDO
    !print*,step," z2=",sum(z**2)," r.z=",sum(r*z)
!$OMP END DO NOWAIT
    endif
  ELSE
    if (step == 1) then
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            p(j, k) = r(j, k)
        ENDDO
    ENDDO
!$OMP END DO NOWAIT
    endif
  ENDIF
  if (step == 1 .or. step ==3) then
!$OMP DO
  DO k=y_min,y_max
    DO j=x_min,x_max
      rro = rro + r(j, k)*p(j, k)
    ENDDO
  ENDDO
!$OMP END DO NOWAIT
  endif
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_cg_init_kernel

SUBROUTINE tea_leaf_cg_calc_w_kernel(x_min,             &
                                                   x_max,             &
                                                   y_min,             &
                                                   y_max,             &
                                                   halo_exchange_depth,             &
                                                   p,                 &
                                                   w,                 &
                                                   Kx,                &
                                                   Ky,                &
                                                   Di,                &
                                                   rx,                &
                                                   ry,                &
                                                   pw                 )

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: w, Kx, Ky, p, Di

    REAL(KIND=8) ::  rx, ry

    INTEGER(KIND=4) :: j,k
    REAL(kind=8) :: pw

  pw = 0.0_8

!$OMP PARALLEL REDUCTION(+:pw)
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            w(j, k) = Di(j,k)*p(j, k)                             &
                - ry*(Ky(j, k+1)*p(j, k+1) + Ky(j, k)*p(j, k-1))  &
                - rx*(Kx(j+1, k)*p(j+1, k) + Kx(j, k)*p(j-1, k))
        ENDDO
        DO j=x_min,x_max
            pw = pw + w(j, k)*p(j, k)
        ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_cg_calc_w_kernel

SUBROUTINE tea_leaf_cg_calc_w_kernel_norxy(x_min,             &
                                                   x_max,             &
                                                   y_min,             &
                                                   y_max,             &
                                                   halo_exchange_depth,             &
                                                   p,                 &
                                                   w,                 &
                                                   Kx,                &
                                                   Ky,                &
                                                   Di,                &
                                                   pw                 )

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: w, Kx, Ky, p, Di

    INTEGER(KIND=4) :: j,k
    REAL(kind=8) :: pw

  pw = 0.0_8

!$OMP PARALLEL REDUCTION(+:pw)
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            w(j, k) = Di(j,k)*p(j, k)                             &
                - (Ky(j, k+1)*p(j, k+1) + Ky(j, k)*p(j, k-1))  &
                - (Kx(j+1, k)*p(j+1, k) + Kx(j, k)*p(j-1, k))
        ENDDO
        DO j=x_min,x_max
            pw = pw + w(j, k)*p(j, k)
        ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_cg_calc_w_kernel_norxy

SUBROUTINE tea_leaf_cg_calc_ur_kernel(x_min,             &
                                                    x_max,             &
                                                    y_min,             &
                                                    y_max,             &
                                                    halo_exchange_depth,             &
                                                    u,                 &
                                                    p,                 &
                                                    r,                 &
                                                    Mi,                &
                                                    w,                 &
                                                    z,                 &
                                                    cp,                     &
                                                    bfp,                     &
                                                    Kx, &
                                                    Ky, &
                                                    Di, &
                                                    rx, &
                                                    ry, &
                                                    alpha,             &
                                                    rrn,               &
                                                    preconditioner_type)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: u, r, Mi, w, z, Kx, Ky, p, Di
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  REAL(KIND=8) :: rx, ry

  INTEGER(KIND=4) :: j,k
  REAL(kind=8) :: alpha, rrn

  rrn = 0.0_8

!$OMP PARALLEL REDUCTION(+:rrn)
  IF (preconditioner_type .NE. TL_PREC_NONE) THEN

    IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN

!$OMP DO
      DO k=y_min,y_max
        DO j=x_min,x_max
          u(j, k) =  u(j, k) + alpha   *p(j, k)
        ENDDO
      ENDDO
!$OMP END DO NOWAIT
!$OMP DO
      DO k=y_min,y_max
        DO j=x_min,x_max
          r(j, k) =  r(j, k) - alpha   *w(j, k)
          z(j, k) =            Mi(j, k)*r(j, k)
          rrn     = rrn      +  r(j, k)*z(j, k)
        ENDDO
      ENDDO
!$OMP END DO NOWAIT

    ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN

!$OMP DO
      DO k=y_min,y_max
        DO j=x_min,x_max
          u(j, k) = u(j, k) + alpha*p(j, k)
        ENDDO
      ENDDO
!$OMP END DO
!$OMP DO
      DO k=y_min,y_max
        DO j=x_min,x_max
          r(j, k) = r(j, k) - alpha*w(j, k)
        ENDDO
      ENDDO
!$OMP END DO

      CALL tea_block_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                             r, z, cp, bfp, Kx, Ky, Di, rx, ry)

!$OMP DO
      DO k=y_min,y_max
        DO j=x_min,x_max
          rrn = rrn + r(j, k)*z(j, k)
        ENDDO
      ENDDO
!$OMP END DO NOWAIT

    ENDIF
  ELSE
!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        u(j, k) = u(j, k) + alpha*p(j, k)
      ENDDO
    ENDDO
!$OMP END DO
!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        r(j, k) = r(j, k) - alpha*w(j, k)
        rrn     = rrn + r(j, k)*r(j, k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_cg_calc_ur_kernel

SUBROUTINE tea_leaf_cg_calc_p_kernel(x_min,             &
                                                   x_max,             &
                                                   y_min,             &
                                                   y_max,             &
                                                   halo_exchange_depth,             &
                                                   p,                 &
                                                   r,                 &
                                                   z,                 &
                                                   beta,              &
                                                   preconditioner_type)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: z, r, p

  INTEGER(KIND=4) :: j,k
  REAL(kind=8) :: beta

!$OMP PARALLEL
  IF (preconditioner_type .NE. TL_PREC_NONE .or. tl_ppcg_active) THEN
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            p(j, k) = z(j, k) + beta*p(j, k)
        ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ELSE
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            p(j, k) = r(j, k) + beta*p(j, k)
        ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_cg_calc_p_kernel

END MODULE

