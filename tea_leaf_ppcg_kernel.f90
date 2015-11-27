MODULE tea_leaf_ppcg_kernel_module

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
                          :: r, sd, kx, ky , z, Mi, Di
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  REAL(KIND=8) :: theta, theta_r, rx, ry

  INTEGER :: j,k

  theta_r = 1.0_8/theta

!$OMP PARALLEL

  IF (preconditioner_type .NE. TL_PREC_NONE) THEN

    !IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
    !  CALL tea_block_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
    !                         r, z, cp, bfp, Kx, Ky, Di, rx, ry)
    !ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
    !  CALL tea_diag_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
    !                         r, z, Mi)
    !ENDIF

!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        sd(j, k) = z(j, k)*theta_r
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ELSE
!$OMP DO
    DO k=y_min,y_max
      DO j=x_min,x_max
        sd(j, k) = r(j, k)*theta_r
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
  ENDIF
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_ppcg_init_sd

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
                                      cp,                &
                                      bfp,                &
                                      Mi,                &
                                      preconditioner_type)
  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: u, r, Kx, Ky, sd , z, Mi, Di
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  INTEGER(KIND=4) :: j,k,depth
  REAL(KIND=8), DIMENSION(:) :: alpha, beta
  REAL(KIND=8) :: smvp, rx, ry

  INTEGER(KIND=4) :: x_min_bound, x_max_bound, y_min_bound, y_max_bound, inner_step

    !write(6,*) "sd:"; write(6,"(9f10.6)") sd
    !if (any(kx /= kx)) then; write(6,*) "Kx:"; write(6,"(9f10.6)") rx*Kx; stop; endif
    !if (any(ky /= ky)) then; write(6,*) "Ky:"; write(6,"(9f10.6)") ry*Ky; stop; endif
    !if (any(di /= di)) then; write(6,*) "Di:"; write(6,"(9f10.6)") Di   ; stop; endif
    !if (any(sd /= sd)) then; write(6,*) "sd:"; write(6,"(9f10.6)") sd   ; stop; endif

!$OMP PARALLEL PRIVATE(smvp)
!$OMP DO
    DO k=y_min_bound,y_max_bound
        DO j=x_min_bound,x_max_bound
            smvp = Di(j,k)*sd(j, k)                                 &
                - ry*(Ky(j, k+1)*sd(j, k+1) + Ky(j, k)*sd(j, k-1))  &
                - rx*(Kx(j+1, k)*sd(j+1, k) + Kx(j, k)*sd(j-1, k))

            r(j, k) = r(j, k) - smvp
            u(j, k) = u(j, k) + sd(j, k)
            !if (j == 1 .and. k == 1) then
            !  write(6,"('1,1:t:',5es12.5)") Di(j,k)*sd(j, k),           &
            !    - ry*Ky(j, k+1)*sd(j, k+1), -ry*Ky(j, k)*sd(j, k-1),  &
            !    - rx*Kx(j+1, k)*sd(j+1, k), -rx*Kx(j, k)*sd(j-1, k)
            !  write(6,"('1,1:c:',5es12.5)") Di(j,k),           &
            !    - ry*Ky(j, k+1), -ry*Ky(j, k),  &
            !    - rx*Kx(j+1, k), -rx*Kx(j, k)
            !  write(6,"('1,1:v:',5es12.5)") sd(j, k),           &
            !    sd(j, k+1), sd(j, k-1),  &
            !    sd(j+1, k), sd(j-1, k)
            !endif
        ENDDO
    ENDDO
!$OMP END DO

    depth=max(x_min-x_min_bound,y_min-y_min_bound,x_max_bound-x_max,y_max_bound-y_max)
    IF (preconditioner_type .NE. TL_PREC_NONE) THEN
      IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
        CALL tea_block_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                               r, z, cp, bfp, Kx, Ky, Di, rx, ry)
      ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
        CALL tea_diag_solve(x_min, x_max, y_min, y_max, halo_exchange_depth, depth,       &
                               r, z, Mi)
      ENDIF
  
!$OMP DO
      DO k=y_min_bound,y_max_bound
          DO j=x_min_bound,x_max_bound
              sd(j, k) = alpha(inner_step)*sd(j, k) + beta(inner_step)*z(j, k)
          ENDDO
      ENDDO
!$OMP END DO
    ELSE
!$OMP DO
      DO k=y_min_bound,y_max_bound
          DO j=x_min_bound,x_max_bound
              sd(j, k) = alpha(inner_step)*sd(j, k) + beta(inner_step)*r(j, k)
          ENDDO
      ENDDO
!$OMP END DO
    ENDIF
!$OMP END PARALLEL

    !write(6,*) "sd:"; write(6,"(9f10.6)") sd
    !if (any(kx /= kx)) then; write(6,*) "Kx:"; write(6,"(9f10.6)") rx*Kx; stop; endif
    !if (any(ky /= ky)) then; write(6,*) "Ky:"; write(6,"(9f10.6)") ry*Ky; stop; endif
    !if (any(di /= di)) then; write(6,*) "Di:"; write(6,"(9f10.6)") Di   ; stop; endif
    !if (any(sd /= sd)) then; write(6,*) "sd:"; write(6,"(9f10.6)") sd   ; stop; endif

    !write(6,*) inner_step," rnorm =",sqrt(sum(r (x_min:x_max,y_min:y_max)**2)),sqrt(sum(r (x_min+1:x_max-1,y_min+1:y_max-1)**2))
    !write(6,*) inner_step," unorm =",sqrt(sum(u (x_min:x_max,y_min:y_max)**2)),sqrt(sum(u (x_min+1:x_max-1,y_min+1:y_max-1)**2))
    !write(6,*) inner_step," sdnorm=",sqrt(sum(sd(x_min:x_max,y_min:y_max)**2)),sqrt(sum(sd(x_min+1:x_max-1,y_min+1:y_max-1)**2))

END SUBROUTINE

SUBROUTINE tea_leaf_kernel_ppcg_inner_nouup(x_min,             &
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
                                      cp,                &
                                      bfp,                &
                                      Mi,                &
                                      preconditioner_type)
  IMPLICIT NONE

  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                          :: u, r, Kx, Ky, sd , z, Mi, Di
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  INTEGER(KIND=4) :: j,k,depth
  REAL(KIND=8), DIMENSION(:) :: alpha, beta
  REAL(KIND=8) :: smvp, rx, ry

  INTEGER(KIND=4) :: x_min_bound, x_max_bound, y_min_bound, y_max_bound, inner_step

    !write(6,*) "sd:"; write(6,"(9f10.6)") sd
    !if (any(kx /= kx)) then; write(6,*) "Kx:"; write(6,"(9f10.6)") rx*Kx; stop; endif
    !if (any(ky /= ky)) then; write(6,*) "Ky:"; write(6,"(9f10.6)") ry*Ky; stop; endif
    !if (any(di /= di)) then; write(6,*) "Di:"; write(6,"(9f10.6)") Di   ; stop; endif
    !if (any(sd /= sd)) then; write(6,*) "sd:"; write(6,"(9f10.6)") sd   ; stop; endif

!$OMP PARALLEL PRIVATE(smvp)
!$OMP DO
    DO k=y_min_bound,y_max_bound
        DO j=x_min_bound,x_max_bound
            smvp = Di(j,k)*sd(j, k)                                 &
                - ry*(Ky(j, k+1)*sd(j, k+1) + Ky(j, k)*sd(j, k-1))  &
                - rx*(Kx(j+1, k)*sd(j+1, k) + Kx(j, k)*sd(j-1, k))

            r(j, k) = r(j, k) - smvp
            !u(j, k) = u(j, k) + sd(j, k)
            !if (j == 1 .and. k == 1) then
            !  write(6,"('1,1:t:',5es12.5)") Di(j,k)*sd(j, k),           &
            !    - ry*Ky(j, k+1)*sd(j, k+1), -ry*Ky(j, k)*sd(j, k-1),  &
            !    - rx*Kx(j+1, k)*sd(j+1, k), -rx*Kx(j, k)*sd(j-1, k)
            !  write(6,"('1,1:c:',5es12.5)") Di(j,k),           &
            !    - ry*Ky(j, k+1), -ry*Ky(j, k),  &
            !    - rx*Kx(j+1, k), -rx*Kx(j, k)
            !  write(6,"('1,1:v:',5es12.5)") sd(j, k),           &
            !    sd(j, k+1), sd(j, k-1),  &
            !    sd(j+1, k), sd(j-1, k)
            !endif
        ENDDO
    ENDDO
!$OMP END DO

    depth=max(x_min-x_min_bound,y_min-y_min_bound,x_max_bound-x_max,y_max_bound-y_max)
    IF (preconditioner_type .NE. TL_PREC_NONE) THEN
      IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
        CALL tea_block_solve(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                               r, z, cp, bfp, Kx, Ky, Di, rx, ry)
      ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
        CALL tea_diag_solve(x_min, x_max, y_min, y_max, halo_exchange_depth, depth,       &
                               r, z, Mi)
      ENDIF
  
!$OMP DO
      DO k=y_min_bound,y_max_bound
          DO j=x_min_bound,x_max_bound
              sd(j, k) = alpha(inner_step)*sd(j, k) + beta(inner_step)*z(j, k)
          ENDDO
      ENDDO
!$OMP END DO
    ELSE
!$OMP DO
      DO k=y_min_bound,y_max_bound
          DO j=x_min_bound,x_max_bound
              sd(j, k) = alpha(inner_step)*sd(j, k) + beta(inner_step)*r(j, k)
          ENDDO
      ENDDO
!$OMP END DO
    ENDIF
!$OMP END PARALLEL

    !write(6,*) "sd:"; write(6,"(9f10.6)") sd
    !if (any(kx /= kx)) then; write(6,*) "Kx:"; write(6,"(9f10.6)") rx*Kx; stop; endif
    !if (any(ky /= ky)) then; write(6,*) "Ky:"; write(6,"(9f10.6)") ry*Ky; stop; endif
    !if (any(di /= di)) then; write(6,*) "Di:"; write(6,"(9f10.6)") Di   ; stop; endif
    !if (any(sd /= sd)) then; write(6,*) "sd:"; write(6,"(9f10.6)") sd   ; stop; endif

    !write(6,*) inner_step," rnorm =",sqrt(sum(r (x_min:x_max,y_min:y_max)**2)),sqrt(sum(r (x_min+1:x_max-1,y_min+1:y_max-1)**2))
    !write(6,*) inner_step," unorm =",sqrt(sum(u (x_min:x_max,y_min:y_max)**2)),sqrt(sum(u (x_min+1:x_max-1,y_min+1:y_max-1)**2))
    !write(6,*) inner_step," sdnorm=",sqrt(sum(sd(x_min:x_max,y_min:y_max)**2)),sqrt(sum(sd(x_min+1:x_max-1,y_min+1:y_max-1)**2))

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
  IF (preconditioner_type .NE. TL_PREC_NONE) THEN
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

END MODULE

