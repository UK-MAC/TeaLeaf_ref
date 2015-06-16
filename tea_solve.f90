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

!>  @brief Driver for the heat conduction kernel
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Invokes the user specified kernel for the heat conduction

MODULE tea_leaf_module

  USE report_module
  USE data_module
  USE tea_leaf_kernel_common_module
  USE tea_leaf_kernel_jacobi_module
  USE tea_leaf_kernel_cg_module
  USE tea_leaf_kernel_ppcg_module
  USE tea_leaf_kernel_cheby_module
  USE update_halo_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf()

  IMPLICIT NONE

!$ INTEGER :: OMP_GET_THREAD_NUM
  INTEGER :: n, t
  REAL(KIND=8) :: ry,rx,old_error,error,exact_error,initial_residual

  INTEGER :: fields(NUM_FIELDS)

  REAL(KIND=8) :: timer,halo_time,solve_time,init_time,reset_time,dot_product_time

  ! For CG solver
  REAL(KIND=8) :: rro, pw, rrn, alpha, beta

  ! For chebyshev solver
  REAL(KIND=8), DIMENSION(max_iters) :: cg_alphas, cg_betas
  REAL(KIND=8), DIMENSION(max_iters) :: ch_alphas, ch_betas
  REAL(KIND=8),SAVE :: eigmin, eigmax, theta, cn
  INTEGER :: est_itc, cheby_calc_steps, max_cheby_iters, info, ppcg_inner_iters
  LOGICAL :: ch_switch_check
  LOGICAL, SAVE :: first=.TRUE.

  INTEGER :: cg_calc_steps

  REAL(KIND=8) :: cg_time, ch_time, total_solve_time, ch_per_it, cg_per_it, iteration_time

  cg_time = 0.0_8
  ch_time = 0.0_8
  cg_calc_steps = 0
  ppcg_inner_iters = 0
  ch_switch_check = .false.

  total_solve_time = 0.0_8
  init_time = 0.0_8
  halo_time = 0.0_8
  solve_time = 0.0_8
  dot_product_time = 0.0_8

  IF (coefficient .NE. RECIP_CONDUCTIVITY .AND. coefficient .NE. conductivity) THEN
    CALL report_error('tea_leaf', 'unknown coefficient option')
  ENDIF

  cheby_calc_steps = 0
  cg_calc_steps = 0

  initial_residual = 0.0_8

  total_solve_time = timer()

  ! INIT
  IF (profiler_on) init_time=timer()

  fields=0
  fields(FIELD_ENERGY1) = 1
  fields(FIELD_DENSITY) = 1

  IF (profiler_on) halo_time=timer()
  CALL update_halo(fields,halo_exchange_depth)
  IF (profiler_on) init_time = init_time + (timer()-halo_time)

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      rx = dt/(chunk%tiles(t)%field%celldx(chunk%tiles(t)%field%x_min)**2)
      ry = dt/(chunk%tiles(t)%field%celldy(chunk%tiles(t)%field%y_min)**2)

      CALL tea_leaf_kernel_init_common(chunk%tiles(t)%field%x_min, &
          chunk%tiles(t)%field%x_max,                                  &
          chunk%tiles(t)%field%y_min,                                  &
          chunk%tiles(t)%field%y_max,                                  &
          halo_exchange_depth,                                  &
          chunk%chunk_neighbours,                             &
          reflective_boundary,                                    &
          chunk%tiles(t)%field%density,                                &
          chunk%tiles(t)%field%energy1,                                &
          chunk%tiles(t)%field%u,                                      &
          chunk%tiles(t)%field%u0,                                      &
          chunk%tiles(t)%field%vector_r,                               &
          chunk%tiles(t)%field%vector_w,                               &
          chunk%tiles(t)%field%vector_Kx,                              &
          chunk%tiles(t)%field%vector_Ky,                              &
          chunk%tiles(t)%field%tri_cp,   &
          chunk%tiles(t)%field%tri_bfp,    &
          chunk%tiles(t)%field%vector_Mi,    &
          rx, ry, tl_preconditioner_type, coefficient)
    ENDDO
  ENDIF

  fields=0
  fields(FIELD_U) = 1

  IF (profiler_on) halo_time=timer()
  CALL update_halo(fields,1)
  IF (profiler_on) init_time = init_time + (timer()-halo_time)

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      CALL tea_leaf_calc_residual(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                        &
          chunk%tiles(t)%field%y_min,                        &
          chunk%tiles(t)%field%y_max,                        &
          halo_exchange_depth,                        &
          chunk%tiles(t)%field%u,                            &
          chunk%tiles(t)%field%u0,                           &
          chunk%tiles(t)%field%vector_r,                     &
          chunk%tiles(t)%field%vector_Kx,                    &
          chunk%tiles(t)%field%vector_Ky,                    &
          rx, ry)
      CALL tea_leaf_calc_2norm_kernel(chunk%tiles(t)%field%x_min,        &
          chunk%tiles(t)%field%x_max,                                    &
          chunk%tiles(t)%field%y_min,                                    &
          chunk%tiles(t)%field%y_max,                                    &
          halo_exchange_depth,                                    &
          chunk%tiles(t)%field%vector_r,                                 &
          initial_residual)
    ENDDO
  ENDIF

  IF (profiler_on) dot_product_time=timer()
  CALL tea_allsum(initial_residual)
  IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
  IF (profiler_on) init_time = init_time + (timer()-dot_product_time)

  old_error = initial_residual

  initial_residual=SQRT(initial_residual)

  IF (parallel%boss.AND.verbose_on) THEN
!$  IF (OMP_GET_THREAD_NUM().EQ.0) THEN
      WRITE(g_out,*)"Initial residual ",initial_residual
!$  ENDIF
  ENDIF

  IF (tl_use_cg .OR. tl_use_chebyshev .OR. tl_use_ppcg) THEN
    ! All 3 of these solvers use the CG kernels
    IF (use_fortran_kernels) THEN
      DO t=1,tiles_per_task
        CALL tea_leaf_kernel_init_cg_fortran(chunk%tiles(t)%field%x_min, &
            chunk%tiles(t)%field%x_max,                                  &
            chunk%tiles(t)%field%y_min,                                  &
            chunk%tiles(t)%field%y_max,                                  &
            halo_exchange_depth,                                  &
            chunk%tiles(t)%field%density,                                &
            chunk%tiles(t)%field%energy1,                                &
            chunk%tiles(t)%field%u,                                      &
            chunk%tiles(t)%field%vector_p,                               &
            chunk%tiles(t)%field%vector_r,                               &
            chunk%tiles(t)%field%vector_Mi,                              &
            chunk%tiles(t)%field%vector_w,                               &
            chunk%tiles(t)%field%vector_z,                               &
            chunk%tiles(t)%field%vector_Kx,                              &
            chunk%tiles(t)%field%vector_Ky,                              &
            chunk%tiles(t)%field%tri_cp,   &
            chunk%tiles(t)%field%tri_bfp,    &
            rx, ry, rro, tl_preconditioner_type)
      ENDDO
    ENDIF

    ! and globally sum rro
    IF (profiler_on) dot_product_time=timer()
    CALL tea_allsum(rro)
    IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
    IF (profiler_on) init_time = init_time + (timer()-dot_product_time)

    ! need to update p when using CG due to matrix/vector multiplication
    fields=0
    fields(FIELD_U) = 1
    fields(FIELD_P) = 1

    IF (profiler_on) halo_time=timer()
    CALL update_halo(fields,1)
    IF (profiler_on) init_time=init_time+(timer()-halo_time)

    fields=0
    fields(FIELD_P) = 1
  ELSEIF (tl_use_jacobi) THEN
    fields=0
    fields(FIELD_U) = 1
  ENDIF

  IF (profiler_on) profiler%tea_init = profiler%tea_init + (timer() - init_time)

  IF (profiler_on) solve_time = timer()

  DO n=1,max_iters

    iteration_time = timer()

    IF (ch_switch_check .EQV. .FALSE.) THEN
      IF ((cheby_calc_steps .GT. 0)) THEN
        ! already started or already have good guesses for eigenvalues
        ch_switch_check = .TRUE.
      ELSE IF ((first .EQV. .FALSE.) .AND. tl_use_ppcg .AND. n .GT. 1) THEN
        ! If using PPCG, it can start almost immediately
        ch_switch_check = .TRUE.
      ELSE IF ((ABS(old_error) .LE. tl_ch_cg_epslim) .AND. (n .GE. tl_ch_cg_presteps)) THEN
        ! Error is less than set limit, and enough steps have passed to get a good eigenvalue guess
        ch_switch_check = .TRUE.
      ELSE
        ! keep doing CG (or jacobi)
        ch_switch_check = .FALSE.
      ENDIF
    ENDIF

    IF ((tl_use_chebyshev .OR. tl_use_ppcg) .AND. ch_switch_check) THEN
      ! on the first chebyshev steps, find the eigenvalues, coefficients,
      ! and expected number of iterations
      IF (cheby_calc_steps .EQ. 0) THEN
        ! maximum number of iterations in chebyshev solver
        max_cheby_iters = max_iters - n + 2

        IF (first) THEN
          ! calculate eigenvalues
          CALL tea_calc_eigenvalues(cg_alphas, cg_betas, eigmin, eigmax, &
              max_iters, n-1, info)
          first=.FALSE.
          IF (info .NE. 0) CALL report_error('tea_leaf', 'Error in calculating eigenvalues')
          eigmin = eigmin * 0.95
          eigmax = eigmax * 1.05
        ENDIF

        IF (tl_use_chebyshev) THEN
          ! calculate chebyshev coefficients
          CALL tea_calc_ch_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
              theta, max_cheby_iters)

          ! don't need to update p any more
          fields = 0
          fields(FIELD_U) = 1
        ELSE IF (tl_use_ppcg) THEN
          ! To avoid some irritating bounds checking in ppcg inner iterations
          max_cheby_iters = max_cheby_iters + halo_exchange_depth

          ! currently also calculate chebyshev coefficients
          CALL tea_calc_ls_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
              theta, tl_ppcg_inner_steps)
        ENDIF

        cn = eigmax/eigmin

        IF (parallel%boss) THEN
!$        IF (OMP_GET_THREAD_NUM().EQ.0) THEN
100 FORMAT("Eigen min",e14.6," Eigen max",e14.6," Condition number",f14.6," Error",e14.6)
            WRITE(g_out,'(a,i3,a,e15.7)') "Switching after ",n," CG its, error ",rro
            WRITE(g_out, 100) eigmin,eigmax,cn,old_error
            WRITE(0,'(a,i3,a,e15.7)') "Switching after ",n," CG its, error ",rro
            WRITE(0, 100) eigmin,eigmax,cn,old_error
!$        ENDIF
        ENDIF
      ENDIF

      IF (tl_use_chebyshev) THEN
        IF (cheby_calc_steps .EQ. 0) THEN
          CALL tea_leaf_cheby_first_step(ch_alphas, ch_betas, fields, &
              old_error, rx, ry, theta, cn, max_cheby_iters, est_itc, solve_time)

          cheby_calc_steps = 1
        ELSE
          IF (use_fortran_kernels) THEN
            DO t=1,tiles_per_task
              CALL tea_leaf_kernel_cheby_iterate(chunk%tiles(t)%field%x_min,&
                          chunk%tiles(t)%field%x_max,                       &
                          chunk%tiles(t)%field%y_min,                       &
                          chunk%tiles(t)%field%y_max,                       &
                          halo_exchange_depth,                       &
                          chunk%tiles(t)%field%u,                           &
                          chunk%tiles(t)%field%u0,                          &
                          chunk%tiles(t)%field%vector_p,                    &
                          chunk%tiles(t)%field%vector_r,                    &
                          chunk%tiles(t)%field%vector_Mi,                   &
                          chunk%tiles(t)%field%vector_w,                    &
                          chunk%tiles(t)%field%vector_z,                    &
                          chunk%tiles(t)%field%vector_Kx,                   &
                          chunk%tiles(t)%field%vector_Ky,                   &
                          chunk%tiles(t)%field%tri_cp,   &
                          chunk%tiles(t)%field%tri_bfp,    &
                          ch_alphas, ch_betas, max_cheby_iters,        &
                          rx, ry, cheby_calc_steps, tl_preconditioner_type)
            ENDDO
          ENDIF

          ! after estimated number of iterations has passed, calc resid.
          ! Leaving 10 iterations between each global reduction won't affect
          ! total time spent much if at all (number of steps spent in
          ! chebyshev is typically O(300+)) but will greatly reduce global
          ! synchronisations needed
          IF ((n .GE. est_itc) .AND. (MOD(n, 10) .eq. 0)) THEN
            IF (use_fortran_kernels) THEN
              DO t=1,tiles_per_task
                CALL tea_leaf_calc_2norm_kernel(chunk%tiles(t)%field%x_min,&
                      chunk%tiles(t)%field%x_max,                          &
                      chunk%tiles(t)%field%y_min,                          &
                      chunk%tiles(t)%field%y_max,                          &
                      halo_exchange_depth,                          &
                      chunk%tiles(t)%field%vector_r,                       &
                      error                                           )
              ENDDO
            ENDIF

            IF (profiler_on) dot_product_time=timer()
            CALL tea_allsum(error)
            IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
            IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)
          ENDIF
        ENDIF
      ELSE IF (tl_use_ppcg) THEN
        IF (use_fortran_kernels) THEN
          DO t=1,tiles_per_task
            CALL tea_leaf_kernel_solve_cg_fortran_calc_w(chunk%tiles(t)%field%x_min,&
                chunk%tiles(t)%field%x_max,                                         &
                chunk%tiles(t)%field%y_min,                                         &
                chunk%tiles(t)%field%y_max,                                         &
                halo_exchange_depth,                                         &
                chunk%tiles(t)%field%vector_p,                                      &
                chunk%tiles(t)%field%vector_w,                                      &
                chunk%tiles(t)%field%vector_Kx,                                     &
                chunk%tiles(t)%field%vector_Ky,                                     &
                rx, ry, pw)
          ENDDO
        ENDIF

        IF (profiler_on) dot_product_time=timer()
        CALL tea_allsum(pw)
        IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
        IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

        alpha = rro/pw

        IF (use_fortran_kernels) THEN
          DO t=1,tiles_per_task
            CALL tea_leaf_kernel_solve_cg_fortran_calc_ur(chunk%tiles(t)%field%x_min,&
                chunk%tiles(t)%field%x_max,                                          &
                chunk%tiles(t)%field%y_min,                                          &
                chunk%tiles(t)%field%y_max,                                          &
                halo_exchange_depth,                                          &
                chunk%tiles(t)%field%u,                                              &
                chunk%tiles(t)%field%vector_p,                                       &
                chunk%tiles(t)%field%vector_r,                                       &
                chunk%tiles(t)%field%vector_Mi,                                      &
                chunk%tiles(t)%field%vector_w,                                       &
                chunk%tiles(t)%field%vector_z,                                       &
                chunk%tiles(t)%field%tri_cp,   &
                chunk%tiles(t)%field%tri_bfp,    &
                chunk%tiles(t)%field%vector_Kx,                              &
                chunk%tiles(t)%field%vector_Ky,                              &
                rx, ry, &
                alpha, rrn, tl_preconditioner_type)
          ENDDO
        ENDIF

        ! not using rrn, so don't do a tea_allsum

        CALL tea_leaf_run_ppcg_inner_steps(ch_alphas, ch_betas, theta, &
            rx, ry, tl_ppcg_inner_steps, solve_time)
        ppcg_inner_iters = ppcg_inner_iters + tl_ppcg_inner_steps

        IF (use_fortran_kernels) THEN
          DO t=1,tiles_per_task
            CALL tea_leaf_ppcg_calc_zrnorm_kernel(chunk%tiles(t)%field%x_min, &
                  chunk%tiles(t)%field%x_max,                           &
                  chunk%tiles(t)%field%y_min,                           &
                  chunk%tiles(t)%field%y_max,                           &
                  halo_exchange_depth,                           &
                  chunk%tiles(t)%field%vector_z,                        &
                  chunk%tiles(t)%field%vector_r,                        &
                  tl_preconditioner_type, rrn)
          ENDDO
        ENDIF

        IF (profiler_on) dot_product_time=timer()
        CALL tea_allsum(rrn)
        IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
        IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

        beta = rrn/rro

        IF (use_fortran_kernels) THEN
          DO t=1,tiles_per_task
            CALL tea_leaf_kernel_solve_cg_fortran_calc_p(chunk%tiles(t)%field%x_min,&
                chunk%tiles(t)%field%x_max,                                         &
                chunk%tiles(t)%field%y_min,                                         &
                chunk%tiles(t)%field%y_max,                                         &
                halo_exchange_depth,                                         &
                chunk%tiles(t)%field%vector_p,                                      &
                chunk%tiles(t)%field%vector_r,                                      &
                chunk%tiles(t)%field%vector_z,                                      &
                beta, tl_preconditioner_type)
          ENDDO
        ENDIF

        error = rrn
        rro = rrn
      ENDIF

      cheby_calc_steps = cheby_calc_steps + 1
    ELSEIF (tl_use_cg .OR. tl_use_chebyshev .OR. tl_use_ppcg) THEN
      fields(FIELD_P) = 1
      cg_calc_steps = cg_calc_steps + 1

      pw = 0.0_08

      IF (use_fortran_kernels) THEN
        DO t=1,tiles_per_task
          CALL tea_leaf_kernel_solve_cg_fortran_calc_w(chunk%tiles(t)%field%x_min,&
              chunk%tiles(t)%field%x_max,                                         &
              chunk%tiles(t)%field%y_min,                                         &
              chunk%tiles(t)%field%y_max,                                         &
              halo_exchange_depth,                                         &
              chunk%tiles(t)%field%vector_p,                                      &
              chunk%tiles(t)%field%vector_w,                                      &
              chunk%tiles(t)%field%vector_Kx,                                     &
              chunk%tiles(t)%field%vector_Ky,                                     &
              rx, ry, pw)
        ENDDO
      ENDIF

      IF (profiler_on) dot_product_time=timer()
      CALL tea_allsum(pw)
      IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
      IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

      alpha = rro/pw
      cg_alphas(n) = alpha

      rrn = 0.0_8

      IF (use_fortran_kernels) THEN
        DO t=1,tiles_per_task
          CALL tea_leaf_kernel_solve_cg_fortran_calc_ur(chunk%tiles(t)%field%x_min,&
              chunk%tiles(t)%field%x_max,                                          &
              chunk%tiles(t)%field%y_min,                                          &
              chunk%tiles(t)%field%y_max,                                          &
              halo_exchange_depth,                                          &
              chunk%tiles(t)%field%u,                                              &
              chunk%tiles(t)%field%vector_p,                                       &
              chunk%tiles(t)%field%vector_r,                                       &
              chunk%tiles(t)%field%vector_Mi,                                      &
              chunk%tiles(t)%field%vector_w,                                       &
              chunk%tiles(t)%field%vector_z,                                       &
              chunk%tiles(t)%field%tri_cp,   &
              chunk%tiles(t)%field%tri_bfp,    &
              chunk%tiles(t)%field%vector_Kx,                              &
              chunk%tiles(t)%field%vector_Ky,                              &
              rx, ry, &
              alpha, rrn, tl_preconditioner_type)
        ENDDO
      ENDIF

      IF (profiler_on) dot_product_time=timer()
      CALL tea_allsum(rrn)
      IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
      IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

      beta = rrn/rro
      cg_betas(n) = beta

      IF (use_fortran_kernels) THEN
        DO t=1,tiles_per_task
          CALL tea_leaf_kernel_solve_cg_fortran_calc_p(chunk%tiles(t)%field%x_min,&
              chunk%tiles(t)%field%x_max,                                         &
              chunk%tiles(t)%field%y_min,                                         &
              chunk%tiles(t)%field%y_max,                                         &
              halo_exchange_depth,                                         &
              chunk%tiles(t)%field%vector_p,                                      &
              chunk%tiles(t)%field%vector_r,                                      &
              chunk%tiles(t)%field%vector_z,                                      &
              beta, tl_preconditioner_type)
        ENDDO
      ENDIF

      error = rrn
      rro = rrn
    ELSEIF (tl_use_jacobi) THEN
      IF (use_fortran_kernels) THEN
        DO t=1,tiles_per_task
          CALL tea_leaf_kernel_jacobi_solve(chunk%tiles(t)%field%x_min,&
              chunk%tiles(t)%field%x_max,                       &
              chunk%tiles(t)%field%y_min,                       &
              chunk%tiles(t)%field%y_max,                       &
              halo_exchange_depth,                       &
              rx,                                          &
              ry,                                          &
              chunk%tiles(t)%field%vector_Kx,                   &
              chunk%tiles(t)%field%vector_Ky,                   &
              error,                                       &
              chunk%tiles(t)%field%u0,                          &
              chunk%tiles(t)%field%u,                           &
              chunk%tiles(t)%field%vector_r)
        ENDDO
      ENDIF

      IF (profiler_on) dot_product_time=timer()
      CALL tea_allsum(error)
      IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
      IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)
    ENDIF

    ! updates u and possibly p
    IF (profiler_on) halo_time = timer()
    CALL update_halo(fields,1)
    IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

    IF (profiler_on) THEN
      IF (tl_use_chebyshev .AND. ch_switch_check) THEN
        ch_time=ch_time+(timer()-iteration_time)
      ELSE
        cg_time=cg_time+(timer()-iteration_time)
      ENDIF
    ENDIF

    error=SQRT(error)

    IF (parallel%boss.AND.verbose_on) THEN
!$    IF (OMP_GET_THREAD_NUM().EQ.0) THEN
        WRITE(g_out,*)"Residual ",error
!$    ENDIF
    ENDIF

    IF (ABS(error) .LT. eps*initial_residual) EXIT

    old_error = error

  ENDDO

  IF (tl_check_result) THEN
    fields = 0
    fields(FIELD_U) = 1

    IF (profiler_on) halo_time = timer()
    CALL update_halo(fields,1)
    IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

    IF (use_fortran_kernels) THEN
      DO t=1,tiles_per_task
        CALL tea_leaf_calc_residual(chunk%tiles(t)%field%x_min,&
            chunk%tiles(t)%field%x_max,                        &
            chunk%tiles(t)%field%y_min,                        &
            chunk%tiles(t)%field%y_max,                        &
            halo_exchange_depth,                        &
            chunk%tiles(t)%field%u,                            &
            chunk%tiles(t)%field%u0,                           &
            chunk%tiles(t)%field%vector_r,                     &
            chunk%tiles(t)%field%vector_Kx,                    &
            chunk%tiles(t)%field%vector_Ky,                    &
            rx, ry)
        CALL tea_leaf_calc_2norm_kernel(chunk%tiles(t)%field%x_min,        &
            chunk%tiles(t)%field%x_max,                                    &
            chunk%tiles(t)%field%y_min,                                    &
            chunk%tiles(t)%field%y_max,                                    &
            halo_exchange_depth,                                    &
            chunk%tiles(t)%field%vector_r,                                 &
            exact_error)
      ENDDO
    ENDIF

    IF (profiler_on) dot_product_time=timer()
    CALL tea_allsum(exact_error)
    IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
    IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

    exact_error = SQRT(exact_error)
  ENDIF

  IF (profiler_on) profiler%tea_solve = profiler%tea_solve + (timer() - solve_time)

  IF (parallel%boss) THEN
!$  IF (OMP_GET_THREAD_NUM().EQ.0) THEN

102 FORMAT('Conduction error ',e14.7)
      WRITE(g_out,102) error/initial_residual
      WRITE(0,102) error/initial_residual

      IF (tl_check_result) THEN
101 FORMAT('EXACT error calculated as', e14.7)
        WRITE(0, 101) exact_error/initial_residual
        WRITE(g_out, 101) exact_error/initial_residual
      ENDIF

      WRITE(g_out,"('Iteration count ',i8)") n-1
      WRITE(0,"('Iteration count ', i8)") n-1
      IF (tl_use_ppcg) THEN
103 FORMAT('PPCG Iteration count', i8, ' (Total ',i8,')')
        WRITE(g_out,103) ppcg_inner_iters, ppcg_inner_iters + n-1
        WRITE(0,103) ppcg_inner_iters, ppcg_inner_iters + n-1
      ENDIF
!$  ENDIF
  ENDIF

  ! RESET
  IF (profiler_on) reset_time=timer()

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      CALL tea_leaf_kernel_finalise(chunk%tiles(t)%field%x_min, &
          chunk%tiles(t)%field%x_max,                           &
          chunk%tiles(t)%field%y_min,                           &
          chunk%tiles(t)%field%y_max,                           &
          halo_exchange_depth,                           &
          chunk%tiles(t)%field%energy1,                         &
          chunk%tiles(t)%field%density,                         &
          chunk%tiles(t)%field%u)
    ENDDO
  ENDIF

  fields=0
  fields(FIELD_ENERGY1) = 1

  IF (profiler_on) halo_time=timer()
  CALL update_halo(fields,1)
  IF (profiler_on) reset_time = reset_time + (timer()-halo_time)

  IF (profiler_on) profiler%tea_reset = profiler%tea_reset + (timer() - reset_time)

  IF (profiler_on .AND. parallel%boss) THEN
    total_solve_time = (timer() - total_solve_time)
    WRITE(0, "(a16,f16.10,a7,i7,a16,f16.10)") "Solve Time",total_solve_time,"Its",n,"Time Per It",total_solve_time/n
    WRITE(g_out, "(a16,f16.10,a7,i7,a16,f16.10)") "Solve Time",total_solve_time,"Its",n,"Time Per It",total_solve_time/n
  ENDIF

  IF (profiler_on .AND. tl_use_chebyshev) THEN
    CALL tea_sum(ch_time)
    CALL tea_sum(cg_time)
    IF (parallel%boss) THEN
      cg_per_it = MERGE((cg_time/cg_calc_steps)/parallel%max_task, 0.0_8, cg_calc_steps .GT. 0)
      ch_per_it = MERGE((ch_time/cheby_calc_steps)/parallel%max_task, 0.0_8, cheby_calc_steps .GT. 0)

      WRITE(0, "(a3, a16, a7, a16, a7)") "", "Time", "Its", "Per it", "Ratio"
      WRITE(0, "(a3, f16.10, i7, f16.10, f7.2)") &
          "CG", cg_time + 0.0_8, cg_calc_steps, cg_per_it, 1.0_8
      WRITE(0, "(a3, f16.10, i7, f16.10, f7.2)") "CH", ch_time + 0.0_8, cheby_calc_steps, &
          ch_per_it, MERGE(ch_per_it/cg_per_it, 0.0_8, cheby_calc_steps .GT. 0)
      WRITE(0, "('Chebyshev actually took ', i6, ' (' i6, ' off guess)')") &
          cheby_calc_steps, cheby_calc_steps-est_itc

      WRITE(g_out, "(a3, a16, a7, a16, a7)") "", "Time", "Its", "Per it", "Ratio"
      WRITE(g_out, "(a3, f16.10, i7, f16.10, f7.2)") &
          "CG", cg_time + 0.0_8, cg_calc_steps, cg_per_it, 1.0_8
      WRITE(g_out, "(a3, f16.10, i7, f16.10, f7.2)") "CH", ch_time + 0.0_8, cheby_calc_steps, &
          ch_per_it, MERGE(ch_per_it/cg_per_it, 0.0_8, cheby_calc_steps .GT. 0)
      WRITE(g_out, "('Chebyshev actually took ', i6, ' (' i6, ' off guess)')") &
          cheby_calc_steps, cheby_calc_steps-est_itc
    ENDIF
  ENDIF

END SUBROUTINE tea_leaf

SUBROUTINE tea_leaf_run_ppcg_inner_steps(ch_alphas, ch_betas, theta, &
    rx, ry, tl_ppcg_inner_steps, solve_time)

  INTEGER :: fields(NUM_FIELDS)
  INTEGER :: t, tl_ppcg_inner_steps, ppcg_cur_step
  REAL(KIND=8) :: rx, ry, theta
  REAL(KIND=8) :: halo_time, timer, solve_time
  REAL(KIND=8), DIMENSION(max_iters) :: ch_alphas, ch_betas

  INTEGER(KIND=4) :: x_min_bound, x_max_bound, y_min_bound, y_max_bound, inner_step, bounds_extra

  fields = 0
  fields(FIELD_U) = 1

  IF (profiler_on) halo_time=timer()
  CALL update_halo(fields,1)
  IF (profiler_on) solve_time = solve_time + (timer() - halo_time)

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      CALL tea_leaf_kernel_ppcg_init_sd(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                              &
          chunk%tiles(t)%field%y_min,                              &
          chunk%tiles(t)%field%y_max,                              &
          halo_exchange_depth,                              &
          chunk%tiles(t)%field%vector_r,                           &
          chunk%tiles(t)%field%vector_Kx,                        &
          chunk%tiles(t)%field%vector_Ky,                        &
          chunk%tiles(t)%field%vector_sd,                          &
          chunk%tiles(t)%field%vector_z,                          &
          chunk%tiles(t)%field%tri_cp,                          &
          chunk%tiles(t)%field%tri_bfp,                          &
          chunk%tiles(t)%field%vector_Mi,                          &
          rx, ry,                          &
          theta, tl_preconditioner_type)
    ENDDO
  ENDIF

  ! inner steps
  DO ppcg_cur_step=1,tl_ppcg_inner_steps,halo_exchange_depth

    fields = 0
    fields(FIELD_SD) = 1
    fields(FIELD_R) = 1

    IF (profiler_on) halo_time = timer()
    CALL update_halo(fields,halo_exchange_depth)
    IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

    inner_step = ppcg_cur_step

    DO bounds_extra = halo_exchange_depth-1, 0, -1
      IF (chunk%chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
        x_min_bound = chunk%tiles(t)%field%x_min
      ELSE
        x_min_bound = chunk%tiles(t)%field%x_min - bounds_extra
      ENDIF

      IF (chunk%chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
        x_max_bound = chunk%tiles(t)%field%x_max
      ELSE
        x_max_bound = chunk%tiles(t)%field%x_max + bounds_extra
      ENDIF

      IF (chunk%chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
        y_min_bound = chunk%tiles(t)%field%y_min
      ELSE
        y_min_bound = chunk%tiles(t)%field%y_min - bounds_extra
      ENDIF

      IF (chunk%chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
        y_max_bound = chunk%tiles(t)%field%y_max
      ELSE
        y_max_bound = chunk%tiles(t)%field%y_max + bounds_extra
      ENDIF

      IF (use_fortran_kernels) THEN
        DO t=1,tiles_per_task
          CALL tea_leaf_kernel_ppcg_inner(chunk%tiles(t)%field%x_min,&
              chunk%tiles(t)%field%x_max,                            &
              chunk%tiles(t)%field%y_min,                            &
              chunk%tiles(t)%field%y_max,                            &
              halo_exchange_depth,                            &
              x_min_bound,                                    &
              x_max_bound,                                    &
              y_min_bound,                                    &
              y_max_bound,                                    &
              ch_alphas, ch_betas,                              &
              rx, ry,                                           &
              inner_step,                                     &
              chunk%tiles(t)%field%u,                                &
              chunk%tiles(t)%field%vector_r,                         &
              chunk%tiles(t)%field%vector_Kx,                        &
              chunk%tiles(t)%field%vector_Ky,                        &
              chunk%tiles(t)%field%vector_sd,                        &
              chunk%tiles(t)%field%vector_z,                          &
              chunk%tiles(t)%field%tri_cp,                          &
              chunk%tiles(t)%field%tri_bfp,                          &
              chunk%tiles(t)%field%vector_Mi,                          &
              tl_preconditioner_type)
        ENDDO
      ENDIF

      fields = 0
      fields(FIELD_SD) = 1

      IF (profiler_on) halo_time = timer()
      CALL update_boundary(fields, bounds_extra)
      IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

      inner_step = inner_step + 1
      IF (inner_step .gt. tl_ppcg_inner_steps) EXIT
    ENDDO
  ENDDO

  fields = 0
  fields(FIELD_P) = 1

END SUBROUTINE tea_leaf_run_ppcg_inner_steps

SUBROUTINE tea_leaf_cheby_first_step(ch_alphas, ch_betas, fields, &
    error, rx, ry, theta, cn, max_cheby_iters, est_itc, solve_time)

  IMPLICIT NONE

  integer :: t, est_itc, max_cheby_iters
  integer, dimension(:) :: fields
  REAL(KIND=8) :: it_alpha, cn, gamm, bb, error, rx, ry, theta
  REAL(KIND=8), DIMENSION(:) :: ch_alphas, ch_betas
  REAL(KIND=8) :: halo_time, timer, dot_product_time, solve_time

  ! calculate 2 norm of u0
  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      CALL tea_leaf_calc_2norm_kernel(chunk%tiles(t)%field%x_min,&
            chunk%tiles(t)%field%x_max,                          &
            chunk%tiles(t)%field%y_min,                          &
            chunk%tiles(t)%field%y_max,                          &
            halo_exchange_depth,                          &
            chunk%tiles(t)%field%u0,                             &
            bb)
    ENDDO
  ENDIF

  IF (profiler_on) dot_product_time=timer()
  CALL tea_allsum(bb)
  IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
  IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

  ! initialise 'p' array
  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      CALL tea_leaf_kernel_cheby_init(chunk%tiles(t)%field%x_min,&
            chunk%tiles(t)%field%x_max,                          &
            chunk%tiles(t)%field%y_min,                          &
            chunk%tiles(t)%field%y_max,                          &
            halo_exchange_depth,                          &
            chunk%tiles(t)%field%u,                              &
            chunk%tiles(t)%field%u0,                             &
            chunk%tiles(t)%field%vector_p,                       &
            chunk%tiles(t)%field%vector_r,                       &
            chunk%tiles(t)%field%vector_Mi,                      &
            chunk%tiles(t)%field%vector_w,                       &
            chunk%tiles(t)%field%vector_z,                       &
            chunk%tiles(t)%field%vector_Kx,                      &
            chunk%tiles(t)%field%vector_Ky,                      &
            chunk%tiles(t)%field%tri_cp,   &
            chunk%tiles(t)%field%tri_bfp,    &
            ch_alphas, ch_betas, max_cheby_iters,           &
            rx, ry, theta, error, tl_preconditioner_type)
    ENDDO
  ENDIF

  IF (profiler_on) halo_time = timer()
  CALL update_halo(fields,1)
  IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      CALL tea_leaf_kernel_cheby_iterate(chunk%tiles(t)%field%x_min,&
          chunk%tiles(t)%field%x_max,                               &
          chunk%tiles(t)%field%y_min,                               &
          chunk%tiles(t)%field%y_max,                               &
          halo_exchange_depth,                               &
          chunk%tiles(t)%field%u,                                   &
          chunk%tiles(t)%field%u0,                                  &
          chunk%tiles(t)%field%vector_p,                            &
          chunk%tiles(t)%field%vector_r,                            &
          chunk%tiles(t)%field%vector_Mi,                           &
          chunk%tiles(t)%field%vector_w,                            &
          chunk%tiles(t)%field%vector_z,                            &
          chunk%tiles(t)%field%vector_Kx,                           &
          chunk%tiles(t)%field%vector_Ky,                           &
          chunk%tiles(t)%field%tri_cp,   &
          chunk%tiles(t)%field%tri_bfp,    &
          ch_alphas, ch_betas, max_cheby_iters,                &
          rx, ry, 1, tl_preconditioner_type)
    ENDDO
  ENDIF

  IF (use_fortran_kernels) THEN
    DO t=1,tiles_per_task
      CALL tea_leaf_calc_2norm_kernel(chunk%tiles(t)%field%x_min,&
            chunk%tiles(t)%field%x_max,                          &
            chunk%tiles(t)%field%y_min,                          &
            chunk%tiles(t)%field%y_max,                          &
            halo_exchange_depth,                          &
            chunk%tiles(t)%field%vector_r,                       &
            error)
    ENDDO
  ENDIF

  IF (profiler_on) dot_product_time=timer()
  CALL tea_allsum(error)
  IF (profiler_on) profiler%dot_product= profiler%dot_product+ (timer() - dot_product_time)
  IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

  it_alpha = eps*bb/(4.0_8*error)
  gamm = (SQRT(cn) - 1.0_8)/(SQRT(cn) + 1.0_8)
  est_itc = NINT(LOG(it_alpha)/(2.0_8*LOG(gamm)))

  IF (parallel%boss) THEN
      WRITE(g_out,'(a11)')"est itc"
      WRITE(g_out,'(11i11)')est_itc
      WRITE(0,'(a11)')"est itc"
      WRITE(0,'(11i11)')est_itc
  ENDIF

END SUBROUTINE tea_leaf_cheby_first_step

END MODULE tea_leaf_module
