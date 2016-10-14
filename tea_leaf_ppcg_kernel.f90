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

USE tea_leaf_common_kernel_module
USE definitions_module, only: tl_ppcg_inner_steps, tl_ppcg_active ! added for ppcg init

IMPLICIT NONE

CONTAINS


SUBROUTINE tea_leaf_kernel_ppcg_init_sd(x_min,             &
                                        x_max,             &
                                        y_min,             &
                                        y_max,             &
                                        halo_exchange_depth,             &
                                        r,                 &
                                        rtemp,		&                                        
                                        kx,                 &
                                        ky,                 &
                                        sd,                &
                                        z,                &
                                        utemp,		&                                        
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
                          :: r, sd, kx, ky , z, Mi, rtemp, utemp
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  REAL(KIND=8) :: theta, theta_r, rx, ry

  INTEGER :: j,k

  theta_r = 1.0_8/theta

!$OMP PARALLEL

  IF (preconditioner_type .NE. TL_PREC_NONE) THEN

!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        sd(j, k) = z(j, k)*theta_r
        rtemp(j,k) = r(j,k)
        utemp(j,k) = sd(j,k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ELSE
!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        sd(j, k) = r(j, k)*theta_r
        rtemp(j,k) = r(j,k)
        utemp(j,k) = sd(j,k) 
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
                           cp,                     &
                           bfp,                    &
                           rx,                     &
                           ry,                     &
                           rro,                    &
                           preconditioner_type,	&
                           ppcg_inner_iters,	&
                           ch_alphas,	&
                           ch_betas,	&
                           theta,	&
                           solve_time,step)

  IMPLICIT NONE

  INTEGER :: preconditioner_type,step
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth, ppcg_inner_iters
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: r, Kx, Ky, z, Mi, p
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  REAL(KIND=8), DIMENSION(:) :: ch_alphas, ch_betas
  INTEGER(KIND=4) :: j,k

  REAL(KIND=8) ::  rro,rx, ry, theta, solve_time

  ! We divide the algorithm up into steps
  ! step = 1, is a CG step
  ! step = 2 or 3, is a partial step of the PPCG algorithm
  
!$OMP PARALLEL  
  IF (step == 1 .OR. step == 3) rro = 0.0_8

  IF (step == 1) THEN

!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
	p(j, k) = 0.0_8
	z(j, k) = 0.0_8
      ENDDO
    ENDDO
!$OMP END DO
  
  ELSEIF (step == 3) THEN
   
!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
	p(j, k) = 0.0_8
      ENDDO
    ENDDO
!$OMP END DO
  
  ENDIF

  IF (preconditioner_type .NE. TL_PREC_NONE .OR. (tl_ppcg_active .AND. step == 3)) THEN
  
    ! We don't apply the preconditioner on the final application of the polynomial acceleration   
    IF (step == 1 .OR. step == 2) THEN
      IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
        CALL tea_block_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                             r, z, cp, bfp, Kx, Ky, rx, ry)
      ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
        CALL tea_diag_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                             r, z, Mi, Kx, Ky, rx, ry)
      ENDIF        
    ENDIF
    
    IF ( step == 1 .OR. step ==3 ) THEN  
!$OMP DO
      DO k=y_min,y_max
          DO j=x_min,x_max
              p(j, k) = z(j, k)
          ENDDO
      ENDDO
!$OMP END DO NOWAIT
    ENDIF    
  ELSE 
    IF (step == 1) THEN  
!$OMP DO
      DO k=y_min,y_max
          DO j=x_min,x_max
              p(j, k) = r(j, k)
          ENDDO
      ENDDO
!$OMP END DO NOWAIT
    ENDIF
  ENDIF  
  IF (step == 1 .OR. step == 3) THEN
!$OMP DO REDUCTION(+:rro)
    DO k=y_min,y_max
      DO j=x_min,x_max
          rro = rro + r(j, k)*p(j, k);
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
                                      rtemp,          &
                                      Kx,                &
                                      Ky,                &
                                      sd,                &
                                      z,                &
                                      utemp,          &
                                      cp,                &
                                      bfp,                &
                                      Mi,                &
                                      preconditioner_type)
  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: u, r, Kx, Ky, sd , z, Mi, rtemp, utemp
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  INTEGER(KIND=4) :: j,k
  REAL(KIND=8), DIMENSION(:) :: alpha, beta
  REAL(KIND=8) :: smvp, rx, ry

  INTEGER(KIND=4) :: x_min_bound, x_max_bound, y_min_bound, y_max_bound, inner_step

!$OMP PARALLEL PRIVATE(smvp)
!$OMP DO
      DO k=y_min_bound,y_max_bound
          DO j=x_min_bound,x_max_bound        
                  smvp = (1.0_8                                           &
                  + ry*(Ky(j, k+1) + Ky(j, k))                        &
                  + rx*(Kx(j+1, k) + Kx(j, k)))*sd(j, k)              &
                  - ry*(Ky(j, k+1)*sd(j, k+1) + Ky(j, k)*sd(j, k-1))  &
                  - rx*(Kx(j+1, k)*sd(j+1, k) + Kx(j, k)*sd(j-1, k))

              rtemp(j, k) = rtemp(j, k) - smvp
          ENDDO
      ENDDO
!$OMP END DO

    IF (preconditioner_type .NE. TL_PREC_NONE) THEN
 
      IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
        CALL tea_block_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                               rtemp, z, cp, bfp, Kx, Ky, rx, ry)
      ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
        CALL tea_diag_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                               rtemp, z, Mi, Kx, Ky, rx, ry)
      ENDIF
  
!$OMP DO
      DO k=y_min_bound,y_max_bound
          DO j=x_min_bound,x_max_bound
              sd(j, k) = alpha(inner_step)*sd(j, k) + beta(inner_step)*z(j, k)
              utemp(j,k) = utemp(j,k) + sd(j,k)             
          ENDDO
      ENDDO
!$OMP END DO
    ELSE
!$OMP DO
      DO k=y_min_bound,y_max_bound
          DO j=x_min_bound,x_max_bound        
              sd(j, k) = alpha(inner_step)*sd(j, k) + beta(inner_step)*rtemp(j, k)
              utemp(j,k) = utemp(j,k) + sd(j,k)  
          ENDDO
      ENDDO
!$OMP END DO
    ENDIF
!$OMP END PARALLEL

END SUBROUTINE

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
  IF (preconditioner_type .NE. TL_PREC_NONE .OR. tl_ppcg_active) THEN
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

end SUBROUTINE tea_leaf_ppcg_calc_zrnorm_kernel

! New routine to update z

SUBROUTINE tea_leaf_ppcg_update_z_kernel(x_min, &
                          x_max,             &
                          y_min,             &
                          y_max,             &
                          halo_exchange_depth, &
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
            z(j,k) = utemp(j, k)
        ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_ppcg_update_z_kernel

! New routine to store the previous residual

SUBROUTINE tea_leaf_ppcg_store_r_kernel(x_min, &
                          x_max,             &
                          y_min,             &
                          y_max,             &
                          halo_exchange_depth,             &
                          r, rstore)

  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: r, rstore
  integer :: j, k

!$OMP PARALLEL
!$OMP DO 
    DO k=y_min,y_max
        DO j=x_min,x_max     
            rstore(j,k) = r(j, k)   
        ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

end SUBROUTINE tea_leaf_ppcg_store_r_kernel

! This does FCG(1) residual compute to minimise rounding error in ppcg

SUBROUTINE tea_leaf_ppcg_calc_rrn_kernel(x_min,              &
                                         x_max,              &
                                         y_min,              &
                                         y_max,              &
                                         halo_exchange_depth,&
                                         r,                  &
                                         rstore,             &
                                         z,                  &
                                         rrn)

  IMPLICIT NONE
  
  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: r, rstore, z
  REAL(KIND=8) :: rrn
  integer :: j, k

  rrn = 0.0_8

!$OMP PARALLEL
!$OMP DO REDUCTION(+:rrn)
    DO k=y_min,y_max
        DO j=x_min,x_max
            rrn = rrn + (r(j, k)- rstore(j,k))*z(j, k)  
        ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_ppcg_calc_rrn_kernel

END MODULE


