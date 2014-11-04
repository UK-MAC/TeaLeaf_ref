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
  USE tea_leaf_kernel_module
  USE tea_leaf_kernel_cg_module
  USE tea_leaf_kernel_ppcg_module
  USE tea_leaf_kernel_cheby_module
  USE update_halo_module

  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf()

  IMPLICIT NONE

!$ INTEGER :: OMP_GET_THREAD_NUM
  INTEGER :: c, n
  REAL(KIND=8) :: ry,rx, error, exact_error

  INTEGER :: fields(NUM_FIELDS)

  REAL(KIND=8) :: timer,halo_time, solve_time, init_time

  ! For CG solver
  REAL(KIND=8) :: rro, pw, rrn, alpha, beta

  ! For chebyshev solver
  REAL(KIND=8), DIMENSION(max_iters) :: cg_alphas, cg_betas
  REAL(KIND=8), DIMENSION(max_iters) :: ch_alphas, ch_betas
  REAL(KIND=8) :: eigmin, eigmax, theta, cn
  INTEGER :: est_itc, cheby_calc_steps, max_cheby_iters, info, switch_step
  LOGICAL :: ch_switch_check

  INTEGER :: cg_calc_steps

  REAL(KIND=8) :: cg_time, ch_time, total_solve_time, ch_per_it, cg_per_it

  cg_time = 0.0_8
  ch_time = 0.0_8
  cg_calc_steps = 0

  total_solve_time = 0.0_8
  init_time = 0.0_8
  halo_time = 0.0_8
  solve_time = 0.0_8

  IF(coefficient .NE. RECIP_CONDUCTIVITY .AND. coefficient .NE. conductivity) THEN
    CALL report_error('tea_leaf', 'unknown coefficient option')
  ENDIF

  error = 1e10
  cheby_calc_steps = 0

  DO c=1,chunks_per_task

    IF(chunks(c)%task.EQ.parallel%task) THEN

      ! INIT
      IF (profiler_on) init_time=timer()

      fields=0
      fields(FIELD_ENERGY1) = 1
      fields(FIELD_DENSITY) = 1
      CALL update_halo(fields,2)

      IF (use_fortran_kernels .OR. use_c_kernels) THEN
        rx = dt/(chunks(c)%field%celldx(chunks(c)%field%x_min)**2)
        ry = dt/(chunks(c)%field%celldy(chunks(c)%field%y_min)**2)
      ENDIF

      IF(tl_use_cg .OR. tl_use_chebyshev .OR. tl_use_ppcg) THEN
        ! All 3 of these solvers use the CG kernels
        IF(use_fortran_kernels) THEN
          CALL tea_leaf_kernel_init_cg_fortran(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                       &
              chunks(c)%field%y_min,                       &
              chunks(c)%field%y_max,                       &
              chunks(c)%field%density,                     &
              chunks(c)%field%energy1,                     &
              chunks(c)%field%u,                           &
              chunks(c)%field%vector_p,                 &
              chunks(c)%field%vector_r,                 &
              chunks(c)%field%vector_Mi,                 &
              chunks(c)%field%vector_w,                 &
              chunks(c)%field%vector_z,                 &
              chunks(c)%field%vector_Kx,                 &
              chunks(c)%field%vector_Ky,                 &
              rx, ry, rro, coefficient)
        ELSEIF(use_C_kernels) THEN
          CALL tea_leaf_kernel_init_cg_c(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                       &
              chunks(c)%field%y_min,                       &
              chunks(c)%field%y_max,                       &
              chunks(c)%field%density,                     &
              chunks(c)%field%energy1,                     &
              chunks(c)%field%u,                           &
              chunks(c)%field%vector_p,                 &
              chunks(c)%field%vector_r,                 &
              chunks(c)%field%vector_Mi,                 &
              chunks(c)%field%vector_w,                 &
              chunks(c)%field%vector_z,                 &
              chunks(c)%field%vector_Kx,                 &
              chunks(c)%field%vector_Ky,                 &
              rx, ry, rro, coefficient)
        ENDIF

        ! need to update p when using CG due to matrix/vector multiplication
        fields=0
        fields(FIELD_U) = 1
        fields(FIELD_P) = 1
        CALL update_halo(fields,1)

        ! and globally sum rro
        CALL tea_allsum(rro)
      ELSEIF(tl_use_jacobi) THEN
        IF (use_fortran_kernels) THEN
          CALL tea_leaf_kernel_init(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                       &
              chunks(c)%field%y_min,                       &
              chunks(c)%field%y_max,                       &
              chunks(c)%field%celldx,                      &
              chunks(c)%field%celldy,                      &
              chunks(c)%field%volume,                      &
              chunks(c)%field%density,                     &
              chunks(c)%field%energy1,                     &
              chunks(c)%field%u0,                          & !u0
              chunks(c)%field%u,                           & !u1
              chunks(c)%field%vector_r,                 & !un
              chunks(c)%field%vector_w,                 & !Kx temp
              chunks(c)%field%vector_z,                 & !Ky temp
              chunks(c)%field%vector_Kx,                 & !Kx
              chunks(c)%field%vector_Ky,                 & !Ky
              coefficient)
        ELSEIF(use_C_kernels) THEN
          CALL tea_leaf_kernel_init_c(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                         &
              chunks(c)%field%y_min,                         &
              chunks(c)%field%y_max,                         &
              chunks(c)%field%celldx,                        &
              chunks(c)%field%celldy,                        &
              chunks(c)%field%volume,                        &
              chunks(c)%field%density,                       &
              chunks(c)%field%energy1,                       &
              chunks(c)%field%u0,                            & !u0
              chunks(c)%field%u,                             & !u1
              chunks(c)%field%vector_r,                   & !un
              chunks(c)%field%vector_w,                   & !Kx temp
              chunks(c)%field%vector_z,                   & !Ky temp
              chunks(c)%field%vector_Kx,                   & !Kx
              chunks(c)%field%vector_Ky,                   & !Ky
              coefficient)
        ENDIF

        fields=0
        fields(FIELD_U) = 1
      ENDIF

      IF(use_fortran_kernels) THEN
        CALL tea_leaf_kernel_cheby_copy_u(chunks(c)%field%x_min,&
          chunks(c)%field%x_max,                       &
          chunks(c)%field%y_min,                       &
          chunks(c)%field%y_max,                       &
          chunks(c)%field%u0,                &
          chunks(c)%field%u)
      ENDIF

      IF (profiler_on) profiler%tea_init = profiler%tea_init + (timer() - init_time)
      IF (profiler_on) solve_time = timer()

      DO n=1,max_iters

        IF (tl_ch_cg_errswitch) THEN
            ! either the error has got below tolerance, or it's already going
            ch_switch_check = (cheby_calc_steps .GT. 0) .OR. (error .LE. tl_ch_cg_epslim)
        ELSE
            ! enough steps have passed and error < 1, otherwise it's nowhere near converging on eigenvalues
            ch_switch_check = (n .GE. tl_ch_cg_presteps) .AND. (error .le. 1.0_8)
        ENDIF

        IF ((tl_use_chebyshev .OR. tl_use_ppcg) .AND. ch_switch_check) THEN
          ! on the first chebyshev steps, find the eigenvalues, coefficients,
          ! and expected number of iterations
          IF (cheby_calc_steps .EQ. 0) THEN
            ! maximum number of iterations in chebyshev solver
            max_cheby_iters = max_iters - n + 2
            rro = error

            ! calculate eigenvalues
            CALL tea_calc_eigenvalues(cg_alphas, cg_betas, eigmin, eigmax, &
                max_iters, n-1, info)

            IF (info .NE. 0) CALL report_error('tea_leaf', 'Error in calculating eigenvalues')

            eigmin = eigmin*0.95
            eigmax = eigmax*1.05

            IF (tl_use_chebyshev) THEN
              ! calculate chebyshev coefficients
              CALL tea_calc_ch_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
                  theta, max_cheby_iters)
            ELSE IF (tl_use_ppcg) THEN
              ! currently also calculate chebyshev coefficients
              ! TODO least squares
              CALL tea_calc_ch_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
                  theta, tl_ppcg_inner_steps)
            ENDIF

            cn = eigmax/eigmin

            IF (parallel%boss) THEN
              WRITE(g_out,'(a,i3,a,e15.7)') "Switching after ",n," steps, error ",rro
              WRITE(g_out,'(4a11)')"eigmin", "eigmax", "cn", "error"
              WRITE(g_out,'(2f11.5,2e11.4)')eigmin, eigmax, cn, error
              WRITE(0,'(a,i3,a,e15.7)') "Switching after ",n," steps, error ",rro
              WRITE(0,'(4a11)')"eigmin", "eigmax", "cn", "error"
              WRITE(0,'(2f11.5,2e11.4)')eigmin, eigmax, cn, error
            ENDIF
          ENDIF

          IF (tl_use_chebyshev) THEN
              ! don't need to update p any more
              fields(FIELD_P) = 0

              IF (cheby_calc_steps .EQ. 0) THEN
                CALL tea_leaf_cheby_first_step(c, ch_alphas, ch_betas, fields, &
                    error, rx, ry, theta, cn, max_cheby_iters, est_itc)

                cheby_calc_steps = 2

                switch_step = n
              ELSE
                  IF(use_fortran_kernels) THEN
                      CALL tea_leaf_kernel_cheby_iterate(chunks(c)%field%x_min,&
                          chunks(c)%field%x_max,                       &
                          chunks(c)%field%y_min,                       &
                          chunks(c)%field%y_max,                       &
                          chunks(c)%field%u,                           &
                          chunks(c)%field%u0,                          &
                          chunks(c)%field%vector_p,                 &
                          chunks(c)%field%vector_r,                 &
                          chunks(c)%field%vector_Mi,                 &
                          chunks(c)%field%vector_w,                 &
                          chunks(c)%field%vector_z,                 &
                          chunks(c)%field%vector_Kx,                 &
                          chunks(c)%field%vector_Ky,                 &
                          ch_alphas, ch_betas, max_cheby_iters,        &
                          rx, ry, cheby_calc_steps)
                  ENDIF

                  ! after estimated number of iterations has passed, calc resid.
                  ! Leaving 10 iterations between each global reduction won't affect
                  ! total time spent much if at all (number of steps spent in
                  ! chebyshev is typically O(300+)) but will greatly reduce global
                  ! synchronisations needed
                  IF ((n .GE. est_itc) .AND. (mod(n, 10) .eq. 0)) THEN
                    IF(use_fortran_kernels) THEN
                      CALL tea_leaf_calc_2norm_kernel(chunks(c)%field%x_min,        &
                            chunks(c)%field%x_max,                       &
                            chunks(c)%field%y_min,                       &
                            chunks(c)%field%y_max,                       &
                            chunks(c)%field%vector_r,                 &
                            error)
                    ENDIF

                    CALL tea_allsum(error)
                  ENDIF
              ENDIF
          ELSE IF (tl_use_ppcg) THEN
            IF (cheby_calc_steps .EQ. 0) THEN
              cheby_calc_steps = 1

              IF(use_fortran_kernels) THEN
                CALL tea_leaf_calc_residual(chunks(c)%field%x_min,&
                    chunks(c)%field%x_max,                       &
                    chunks(c)%field%y_min,                       &
                    chunks(c)%field%y_max,                       &
                    chunks(c)%field%u,                           &
                    chunks(c)%field%u0,                 &
                    chunks(c)%field%vector_r,                 &
                    chunks(c)%field%vector_Kx,                 &
                    chunks(c)%field%vector_Ky,                 &
                    rx, ry)
              ENDIF

              IF (profiler_on) halo_time = timer()
              ! update p
              CALL update_halo(fields,1)
              IF (profiler_on) profiler%halo_exchange = profiler%halo_exchange + (timer() - halo_time)
              IF (profiler_on) solve_time = solve_time - halo_time

              CALL tea_allsum(rro)
            ENDIF

            IF(use_fortran_kernels) THEN
              CALL tea_leaf_kernel_solve_cg_fortran_calc_w(chunks(c)%field%x_min,&
                  chunks(c)%field%x_max,                       &
                  chunks(c)%field%y_min,                       &
                  chunks(c)%field%y_max,                       &
                  chunks(c)%field%vector_p,                 &
                  chunks(c)%field%vector_w,                 &
                  chunks(c)%field%vector_Kx,                 &
                  chunks(c)%field%vector_Ky,                 &
                  rx, ry, pw)
            ENDIF

            CALL tea_allsum(pw)
            alpha = rro/pw

            IF(use_fortran_kernels) THEN
              CALL tea_leaf_kernel_solve_cg_fortran_calc_ur(chunks(c)%field%x_min,&
                  chunks(c)%field%x_max,                       &
                  chunks(c)%field%y_min,                       &
                  chunks(c)%field%y_max,                       &
                  chunks(c)%field%u,                           &
                  chunks(c)%field%vector_p,                 &
                  chunks(c)%field%vector_r,                 &
                  chunks(c)%field%vector_Mi,                 &
                  chunks(c)%field%vector_w,                 &
                  chunks(c)%field%vector_z,                 &
                  alpha, rrn)
            ENDIF

            ! not using rrn, so don't do a tea_allsum

            CALL tea_leaf_run_ppcg_inner_steps(ch_alphas, ch_betas, theta, &
                rx, ry, tl_ppcg_inner_steps, c)

            IF(use_fortran_kernels) THEN
              CALL tea_leaf_calc_2norm_kernel(chunks(c)%field%x_min,        &
                    chunks(c)%field%x_max,                       &
                    chunks(c)%field%y_min,                       &
                    chunks(c)%field%y_max,                       &
                    chunks(c)%field%vector_r,                 &
                    rrn)
            ENDIF

            CALL tea_allsum(rrn)

            beta = rrn/rro

            IF(use_fortran_kernels) THEN
              CALL tea_leaf_kernel_solve_cg_fortran_calc_p(chunks(c)%field%x_min,&
                  chunks(c)%field%x_max,                       &
                  chunks(c)%field%y_min,                       &
                  chunks(c)%field%y_max,                       &
                  chunks(c)%field%vector_p,                 &
                  chunks(c)%field%vector_r,                 &
                  chunks(c)%field%vector_z,                 &
                  beta)
            ENDIF

            error = rrn
            rro = rrn
          ENDIF

          cheby_calc_steps = cheby_calc_steps + 1
        ELSEIF(tl_use_cg .OR. tl_use_chebyshev .OR. tl_use_ppcg) THEN
          fields(FIELD_P) = 1
          cg_calc_steps = cg_calc_steps + 1

          IF(use_fortran_kernels) THEN
            CALL tea_leaf_kernel_solve_cg_fortran_calc_w(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                                         &
                chunks(c)%field%y_min,                                         &
                chunks(c)%field%y_max,                                         &
                chunks(c)%field%vector_p,                                   &
                chunks(c)%field%vector_w,                                   &
                chunks(c)%field%vector_Kx,                                   &
                chunks(c)%field%vector_Ky,                                   &
                rx, ry, pw)
          ELSEIF(use_c_kernels) THEN
            CALL tea_leaf_kernel_solve_cg_c_calc_w(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                                   &
                chunks(c)%field%y_min,                                   &
                chunks(c)%field%y_max,                                   &
                chunks(c)%field%vector_p,                             &
                chunks(c)%field%vector_w,                             &
                chunks(c)%field%vector_Kx,                             &
                chunks(c)%field%vector_Ky,                             &
                rx, ry, pw)
          ENDIF

          CALL tea_allsum(pw)
          alpha = rro/pw
          IF(tl_use_chebyshev .OR. tl_use_ppcg) cg_alphas(n) = alpha

          IF(use_fortran_kernels) THEN
            CALL tea_leaf_kernel_solve_cg_fortran_calc_ur(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                       &
                chunks(c)%field%y_min,                       &
                chunks(c)%field%y_max,                       &
                chunks(c)%field%u,                           &
                chunks(c)%field%vector_p,                 &
                chunks(c)%field%vector_r,                 &
                chunks(c)%field%vector_Mi,                 &
                chunks(c)%field%vector_w,                 &
                chunks(c)%field%vector_z,                 &
                alpha, rrn)
          ELSEIF(use_c_kernels) THEN
            CALL tea_leaf_kernel_solve_cg_c_calc_ur(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                                    &
                chunks(c)%field%y_min,                                    &
                chunks(c)%field%y_max,                                    &
                chunks(c)%field%u,                                        &
                chunks(c)%field%vector_p,                              &
                chunks(c)%field%vector_r,                              &
                chunks(c)%field%vector_Mi,                              &
                chunks(c)%field%vector_w,                              &
                chunks(c)%field%vector_z,                              &
                alpha, rrn)
          ENDIF

          CALL tea_allsum(rrn)
          beta = rrn/rro
          IF(tl_use_chebyshev .OR. tl_use_ppcg) cg_betas(n) = beta

          IF(use_fortran_kernels) THEN
            CALL tea_leaf_kernel_solve_cg_fortran_calc_p(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                                         &
                chunks(c)%field%y_min,                                         &
                chunks(c)%field%y_max,                                         &
                chunks(c)%field%vector_p,                                   &
                chunks(c)%field%vector_r,                                   &
                chunks(c)%field%vector_z,                                   &
                beta)
          ELSEIF(use_c_kernels) THEN
            CALL tea_leaf_kernel_solve_cg_c_calc_p(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                                   &
                chunks(c)%field%y_min,                                   &
                chunks(c)%field%y_max,                                   &
                chunks(c)%field%vector_p,                             &
                chunks(c)%field%vector_r,                             &
                chunks(c)%field%vector_z,                             &
                beta)
          ENDIF

          error = rrn
          rro = rrn
        ELSEIF(tl_use_jacobi) THEN
          IF(use_fortran_kernels) THEN
            CALL tea_leaf_kernel_solve(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                       &
                chunks(c)%field%y_min,                       &
                chunks(c)%field%y_max,                       &
                rx,                                          &
                ry,                                          &
                chunks(c)%field%vector_Kx,                 &
                chunks(c)%field%vector_Ky,                 &
                error,                                       &
                chunks(c)%field%u0,                          &
                chunks(c)%field%u,                           &
                chunks(c)%field%vector_r)
          ELSEIF(use_C_kernels) THEN
            CALL tea_leaf_kernel_solve_c(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                         &
                chunks(c)%field%y_min,                         &
                chunks(c)%field%y_max,                         &
                rx,                                            &
                ry,                                            &
                chunks(c)%field%vector_Kx,                   &
                chunks(c)%field%vector_Ky,                   &
                error,                                         &
                chunks(c)%field%u0,                            &
                chunks(c)%field%u,                             &
                chunks(c)%field%vector_r)
          ENDIF

          CALL tea_allsum(error)
        ENDIF

        ! updates u and possibly p
        IF (profiler_on) halo_time = timer()
        CALL update_halo(fields,1)
        IF (profiler_on) profiler%halo_exchange = profiler%halo_exchange + (timer() - halo_time)
        IF (profiler_on) solve_time = solve_time - halo_time

        IF (profiler_on) THEN
          IF (tl_use_chebyshev .AND. ch_switch_check) THEN
              ch_time=ch_time+(timer()-profiler%tea_solve)
          ELSE
              cg_time=cg_time+(timer()-profiler%tea_solve)
          ENDIF
          total_solve_time = total_solve_time + (timer() - profiler%tea_solve)
        ENDIF

        IF (abs(error) .LT. eps) EXIT

      ENDDO

      IF (profiler_on) profiler%tea_solve = profiler%tea_solve + (timer() - solve_time)

      IF (tl_check_result) THEN
        IF(use_fortran_kernels) THEN
          CALL tea_leaf_calc_residual(chunks(c)%field%x_min,&
              chunks(c)%field%x_max,                       &
              chunks(c)%field%y_min,                       &
              chunks(c)%field%y_max,                       &
              chunks(c)%field%u,                           &
              chunks(c)%field%u0,                 &
              chunks(c)%field%vector_r,                 &
              chunks(c)%field%vector_Kx,                 &
              chunks(c)%field%vector_Ky,                 &
              rx, ry)
          CALL tea_leaf_calc_2norm_kernel(chunks(c)%field%x_min,        &
              chunks(c)%field%x_max,                       &
              chunks(c)%field%y_min,                       &
              chunks(c)%field%y_max,                       &
              chunks(c)%field%vector_r,                 &
              exact_error)
        ENDIF

        CALL tea_allsum(exact_error)
      ENDIF

      IF (parallel%boss) THEN
!$      IF(OMP_GET_THREAD_NUM().EQ.0) THEN
          WRITE(g_out,"('Conduction error ',e14.7)") error
          WRITE(0,"('Conduction error ',e14.7)") error

          IF (tl_check_result) THEN
            WRITE(0,"('EXACT error calculated as', e14.7)") exact_error
            WRITE(g_out,"('EXACT error calculated as', e14.7)") exact_error
          ENDIF

          WRITE(g_out,"('Iteration count ',i8)") n-1
          WRITE(0,"('Iteration count ', i8)") n-1
!$      ENDIF
      ENDIF

      ! RESET
      IF(use_fortran_kernels) THEN
          CALL tea_leaf_kernel_finalise(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                           &
              chunks(c)%field%y_min,                           &
              chunks(c)%field%y_max,                           &
              chunks(c)%field%energy1,                         &
              chunks(c)%field%density,                         &
              chunks(c)%field%u)
      ELSEIF(use_C_kernels) THEN
          CALL tea_leaf_kernel_finalise_c(chunks(c)%field%x_min,&
              chunks(c)%field%x_max,                            &
              chunks(c)%field%y_min,                            &
              chunks(c)%field%y_max,                            &
              chunks(c)%field%energy1,                          &
              chunks(c)%field%density,                          &
              chunks(c)%field%u)
      ENDIF

      fields=0
      fields(FIELD_ENERGY1) = 1
      CALL update_halo(fields,1)

    ENDIF

  ENDDO

  IF (profiler_on .AND. parallel%boss) THEN
    WRITE(0, "(a16, a7, a16)") "Time", "Steps", "Per it"
    WRITE(0, "(f16.10, i7, f16.10, f7.2)") total_solve_time, n, total_solve_time/n
    WRITE(g_out, "(a16, a7, a16)") "Time", "Steps", "Per it"
    WRITE(g_out, "(f16.10, i7, f16.10, f7.2)") total_solve_time, n, total_solve_time/n
  ENDIF

  IF (profiler_on .AND. tl_use_chebyshev) THEN
    CALL tea_sum(ch_time)
    CALL tea_sum(cg_time)
    IF (parallel%boss) THEN
      cg_per_it = MERGE((cg_time/cg_calc_steps)/parallel%max_task, 0.0_8, cg_calc_steps .GT. 0)
      ch_per_it = MERGE((ch_time/cheby_calc_steps)/parallel%max_task, 0.0_8, cheby_calc_steps .GT. 0)

      WRITE(0, "(a3, a16, a7, a16, a7)") "", "Time", "Steps", "Per it", "Ratio"
      WRITE(0, "(a3, f16.10, i7, f16.10, f7.2)") &
          "CG", cg_time + 0.0_8, cg_calc_steps, cg_per_it, 1.0_8
      WRITE(0, "(a3, f16.10, i7, f16.10, f7.2)") "CH", ch_time + 0.0_8, cheby_calc_steps, &
          ch_per_it, MERGE(ch_per_it/cg_per_it, 0.0_8, cheby_calc_steps .GT. 0)
      WRITE(0, "('Chebyshev actually took ', i6, ' (' i6, ' off guess)')") &
          cheby_calc_steps, cheby_calc_steps-est_itc

      WRITE(g_out, "(a3, a16, a7, a16, a7)") "", "Time", "Steps", "Per it", "Ratio"
      WRITE(g_out, "(a3, f16.10, i7, f16.10, f7.2)") &
          "CG", cg_time + 0.0_8, cg_calc_steps, cg_per_it, 1.0_8
      WRITE(g_out, "(a3, f16.10, i7, f16.10, f7.2)") "CH", ch_time + 0.0_8, cheby_calc_steps, &
          ch_per_it, MERGE(ch_per_it/cg_per_it, 0.0_8, cheby_calc_steps .GT. 0)
      WRITE(g_out, "('Chebyshev actually took ', i6, ' (' i6, ' off guess)')") &
          cheby_calc_steps, cheby_calc_steps-est_itc
    ENDIF
  ENDIF

END SUBROUTINE tea_leaf

SUBROUTINE tea_leaF_run_ppcg_inner_steps(ch_alphas, ch_betas, theta, &
    rx, ry, tl_ppcg_inner_steps, c)

  INTEGER :: fields(NUM_FIELDS)
  INTEGER :: c, tl_ppcg_inner_steps, ppcg_cur_step
  REAL(KIND=8) :: rx, ry, theta
  REAL(KIND=8), DIMENSION(max_iters) :: ch_alphas, ch_betas

  IF(use_fortran_kernels) THEN
    CALL tea_leaf_kernel_ppcg_init_sd(chunks(c)%field%x_min,&
        chunks(c)%field%x_max,                       &
        chunks(c)%field%y_min,                       &
        chunks(c)%field%y_max,                       &
        chunks(c)%field%vector_r,                 &
        chunks(c)%field%vector_sd,                 &
        theta)
  ENDIF

  fields = 0
  fields(FIELD_SD) = 1

  ! inner steps
  DO ppcg_cur_step=1,tl_ppcg_inner_steps
    CALL update_halo(fields,1)

    IF(use_fortran_kernels) THEN
      CALL tea_leaf_kernel_ppcg_inner(chunks(c)%field%x_min,&
          chunks(c)%field%x_max,                       &
          chunks(c)%field%y_min,                       &
          chunks(c)%field%y_max,                       &
          ppcg_cur_step, &
          ch_alphas, ch_betas, &
          rx, ry, &
          chunks(c)%field%u,                           &
          chunks(c)%field%vector_r,                 &
          chunks(c)%field%vector_Kx,                 &
          chunks(c)%field%vector_Ky,                 &
          chunks(c)%field%vector_sd)
    ENDIF
  ENDDO

  fields = 0
  fields(FIELD_P) = 1
END SUBROUTINE

SUBROUTINE tea_leaf_cheby_first_step(c, ch_alphas, ch_betas, fields, &
    error, rx, ry, theta, cn, max_cheby_iters, est_itc)

  IMPLICIT NONE

  integer :: c, est_itc, max_cheby_iters
  integer, dimension(:) :: fields
  REAL(KIND=8) :: it_alpha, cn, gamm, bb, error, rx, ry, theta
  REAL(KIND=8), DIMENSION(:) :: ch_alphas, ch_betas

  ! calculate 2 norm of u0
  IF(use_fortran_kernels) THEN
    CALL tea_leaf_calc_2norm_kernel(chunks(c)%field%x_min,        &
          chunks(c)%field%x_max,                       &
          chunks(c)%field%y_min,                       &
          chunks(c)%field%y_max,                       &
          chunks(c)%field%u0,                 &
          bb)
  ENDIF

  CALL tea_allsum(bb)

  ! initialise 'p' array
  IF(use_fortran_kernels) THEN
    CALL tea_leaf_kernel_cheby_init(chunks(c)%field%x_min,&
          chunks(c)%field%x_max,                       &
          chunks(c)%field%y_min,                       &
          chunks(c)%field%y_max,                       &
          chunks(c)%field%u,                           &
          chunks(c)%field%u0,                 &
          chunks(c)%field%vector_p,                 &
          chunks(c)%field%vector_r,                 &
          chunks(c)%field%vector_Mi,                 &
          chunks(c)%field%vector_w,                 &
          chunks(c)%field%vector_z,                 &
          chunks(c)%field%vector_Kx,                 &
          chunks(c)%field%vector_Ky,                 &
          ch_alphas, ch_betas, max_cheby_iters, &
          rx, ry, theta, error)
  ENDIF

  CALL update_halo(fields,1)

  IF(use_fortran_kernels) THEN
      CALL tea_leaf_kernel_cheby_iterate(chunks(c)%field%x_min,&
          chunks(c)%field%x_max,                       &
          chunks(c)%field%y_min,                       &
          chunks(c)%field%y_max,                       &
          chunks(c)%field%u,                           &
          chunks(c)%field%u0,                          &
          chunks(c)%field%vector_p,                 &
          chunks(c)%field%vector_r,                 &
          chunks(c)%field%vector_Mi,                 &
          chunks(c)%field%vector_w,                 &
          chunks(c)%field%vector_z,                 &
          chunks(c)%field%vector_Kx,                 &
          chunks(c)%field%vector_Ky,                 &
          ch_alphas, ch_betas, max_cheby_iters,        &
          rx, ry, 1)
  ENDIF

  IF(use_fortran_kernels) THEN
    CALL tea_leaf_calc_2norm_kernel(chunks(c)%field%x_min,        &
          chunks(c)%field%x_max,                       &
          chunks(c)%field%y_min,                       &
          chunks(c)%field%y_max,                       &
          chunks(c)%field%vector_r,                 &
          error)
  ENDIF

  CALL tea_allsum(error)

  it_alpha = EPSILON(1.0_8)*bb/(4.0_8*error)
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
