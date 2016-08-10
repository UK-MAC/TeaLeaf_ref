
MODULE tea_leaf_cheby_module

  USE tea_leaf_cheby_kernel_module
  USE definitions_module
  USE tea_leaf_common_module
  USE update_halo_module
  
  IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_cheby_init(theta)

  IMPLICIT NONE

  INTEGER :: t
  REAL(KIND=8) :: theta

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
            chunk%tiles(t)%field%rx,    &
            chunk%tiles(t)%field%ry,    &
            theta, tl_preconditioner_type)
    ENDDO
  ENDIF

END SUBROUTINE tea_leaf_cheby_init

SUBROUTINE tea_leaf_cheby_iterate(ch_alphas, ch_betas, max_cheby_iters, cheby_calc_steps)

  IMPLICIT NONE

  INTEGER :: t, cheby_calc_steps, max_cheby_iters
  REAL(KIND=8), DIMENSION(:) :: ch_alphas, ch_betas

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
                  chunk%tiles(t)%field%rx,  &
                  chunk%tiles(t)%field%ry,  &
                  cheby_calc_steps, tl_preconditioner_type)
    ENDDO
  ENDIF

END SUBROUTINE tea_leaf_cheby_iterate

SUBROUTINE tqli(d,e,n, info)
    ! http://physics.sharif.edu/~jafari/fortran-codes/lanczos/tqli.f90
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(n) :: d,e
    INTEGER :: i,iter,l,m,n,info
    REAL(KIND=8) :: b,c,dd,f,g,p,r,s
    e(:)=eoshift(e(:),1)
    info = 0
    DO l=1,n
      iter=0
      iterate: DO
        DO m=l,n-1
          dd=ABS(d(m))+ABS(d(m+1))
          IF (ABS(e(m))+dd == dd) EXIT
        ENDDO
        IF (m == l) EXIT iterate
        IF (iter == 30) THEN
          info=1
          RETURN
        ENDIF
        iter=iter+1
        g=(d(l+1)-d(l))/(2.0_8*e(l))
        r=SQRT(g**2.0_8+1.0_8**2.0_8)
        g=d(m)-d(l)+e(l)/(g+SIGN(r,g))
        s=1.0_8
        c=1.0_8
        p=0.0_8
        DO i=m-1,l,-1
          f=s*e(i)
          b=c*e(i)
          r=SQRT(f**2.0_8+g**2.0_8)
          e(i+1)=r
          IF (r == 0.0_8) THEN
            d(i+1)=d(i+1)-p
            e(m)=0.0_8
            CYCLE iterate
          ENDIF
          s=f/r
          c=g/r
          g=d(i+1)-p
          r=(d(i)-g)*s+2.0_8*c*b
          p=s*r
          d(i+1)=g+p
          g=c*r-b
        ENDDO
        d(l)=d(l)-p
        e(l)=g
        e(m)=0.0_8
      END DO iterate
    END DO
END SUBROUTINE tqli

SUBROUTINE tea_calc_eigenvalues(cg_alphas, cg_betas, eigmin, eigmax, &
                                max_iters, tl_ch_cg_presteps, info)

  INTEGER :: tl_ch_cg_presteps, max_iters
  REAL(KIND=8), DIMENSION(max_iters) :: cg_alphas, cg_betas
  REAL(KIND=8), DIMENSION(tl_ch_cg_presteps) :: diag, offdiag
  ! z not used for this
  REAL(KIND=8) :: eigmin, eigmax, tmp
  INTEGER :: n, info
  LOGICAL :: swapped

  diag = 0
  offdiag = 0

  DO n=1,tl_ch_cg_presteps
    diag(n) = 1.0_8/cg_alphas(n)
    IF (n .GT. 1) diag(n) = diag(n) + cg_betas(n-1)/cg_alphas(n-1)
    IF (n .LT. tl_ch_cg_presteps) offdiag(n+1) = SQRT(cg_betas(n))/cg_alphas(n)
  ENDDO

  CALL tqli(diag, offdiag, tl_ch_cg_presteps, info)

  ! could just call this instead
  !offdiag(:)=eoshift(offdiag(:),1)
  !CALL dsterf(tl_ch_cg_presteps, diag, offdiag, info)

  IF (info .NE. 0) RETURN

  ! bubble sort eigenvalues
  DO
    DO n=1,tl_ch_cg_presteps-1
      IF (diag(n) .GE. diag(n+1)) THEN
        tmp = diag(n)
        diag(n) = diag(n+1)
        diag(n+1) = tmp
        swapped = .TRUE.
      ENDIF
    ENDDO
    IF (.NOT. swapped) EXIT
    swapped = .FALSE.
  ENDDO

  eigmin = diag(1)
  eigmax = diag(tl_ch_cg_presteps)

  IF (eigmin .LT. 0.0_8 .OR. eigmax .LT. 0.0_8) info = 1

END SUBROUTINE tea_calc_eigenvalues

SUBROUTINE tea_calc_ch_coefs(ch_alphas, ch_betas, eigmin, eigmax, &
    theta, max_cheby_iters)

  INTEGER :: n, max_cheby_iters
  REAL(KIND=8), DIMENSION(max_cheby_iters) :: ch_alphas, ch_betas
  REAL(KIND=8) :: eigmin, eigmax

  REAL(KIND=8) :: theta, delta, sigma, rho_old, rho_new, cur_alpha, cur_beta

  theta = (eigmax + eigmin)/2
  delta = (eigmax - eigmin)/2
  sigma = theta/delta

  rho_old = 1.0_8/sigma

  DO n=1,max_cheby_iters
    rho_new = 1.0_8/(2.0_8*sigma - rho_old)

    cur_alpha = rho_new*rho_old
    cur_beta = 2.0_8*rho_new/delta

    ch_alphas(n) = cur_alpha
    ch_betas(n) = cur_beta

    rho_old = rho_new
  ENDDO

END SUBROUTINE tea_calc_ch_coefs

! Moved from tea_solve.f90 
SUBROUTINE tea_leaf_cheby_first_step(ch_alphas, ch_betas, fields, &
    error, theta, cn, max_cheby_iters, est_itc, solve_time)

  IMPLICIT NONE

  INTEGER :: est_itc, max_cheby_iters
  INTEGER, DIMENSION(:) :: fields
  REAL(KIND=8) :: it_alpha, cn, gamm, bb, error, theta
  REAL(KIND=8), DIMENSION(:) :: ch_alphas, ch_betas
  REAL(KIND=8) :: halo_time, timer, dot_product_time, solve_time

  ! calculate 2 norm of u0
  CALL tea_leaf_calc_2norm(0, bb)

  IF (profiler_on) dot_product_time=timer()
  CALL tea_allsum(bb)
  IF (profiler_on) solve_time = solve_time + (timer()-dot_product_time)

  ! initialise 'p' array
  CALL tea_leaf_cheby_init(theta)

  IF (profiler_on) halo_time = timer()
  CALL update_halo(fields,1)
  IF (profiler_on) solve_time = solve_time + (timer()-halo_time)

  CALL tea_leaf_cheby_iterate(ch_alphas, ch_betas, max_cheby_iters, 1)

  CALL tea_leaf_calc_2norm(1, error)

  IF (profiler_on) dot_product_time=timer()
  CALL tea_allsum(error)
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

END MODULE tea_leaf_cheby_module

