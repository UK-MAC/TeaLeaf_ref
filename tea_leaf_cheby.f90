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
!>  @details Implicitly calculates the change in temperature using accelerated Chebyshev method

MODULE tea_leaf_kernel_cheby_module

IMPLICIT NONE

CONTAINS

SUBROUTINE tea_leaf_calc_2norm_kernel(x_min, &
                          x_max,             &
                          y_min,             &
                          y_max,             &
                          arr,               &
                          norm)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: arr
  REAL(KIND=8) :: norm
  integer :: j, k

  norm = 0.0_8

!$OMP PARALLEL
!$OMP DO REDUCTION(+:norm)
    DO k=y_min,y_max
        DO j=x_min,x_max
            norm = norm + arr(j, k)*arr(j, k)
        ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

end SUBROUTINE tea_leaf_calc_2norm_kernel

SUBROUTINE tea_leaf_kernel_cheby_init(x_min,             &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           u,                &
                           u0,                &
                           p,                &
                           r,            &
                           Mi,            &
                           w,     &
                           z,            &
                           Kx,                &
                           Ky,  &
                           ch_alphas, &
                           ch_betas, &
                           max_cheby_iters, &
                           rx, &
                           ry, &
                           theta, &
                           error)
  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: w
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: p, r, Mi, z
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: Kx, Ky

  INTEGER :: j,k, max_cheby_iters
  REAL(KIND=8) ::  rx, ry, error, theta
  REAL(KIND=8), DIMENSION(max_cheby_iters) :: ch_alphas, ch_betas

!$OMP PARALLEL
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            w(j, k) = (1.0_8                                      &
                + ry*(Ky(j, k+1) + Ky(j, k))                      &
                + rx*(Kx(j+1, k) + Kx(j, k)))*u(j, k)             &
                - ry*(Ky(j, k+1)*u(j, k+1) + Ky(j, k)*u(j, k-1))  &
                - rx*(Kx(j+1, k)*u(j+1, k) + Kx(j, k)*u(j-1, k))
            r(j, k) = u0(j, k) - w(j, k)
            !z(j, k) = Mi(j, k)*r(j, k)
          p(j, k) = (Mi(j, k)*r(j, k))/theta
        ENDDO
    ENDDO
!$OMP END DO
!$OMP DO
  DO k=y_min,y_max
      DO j=x_min,x_max
          u(j, k) = u(j, k) + p(j, k)
      ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE

SUBROUTINE tea_leaf_kernel_cheby_iterate(x_min,             &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           u,                &
                           u0,                &
                           p,                &
                           r,            &
                           Mi,            &
                           w,     &
                           z,            &
                           Kx,                &
                           Ky,  &
                           ch_alphas, &
                           ch_betas, &
                           max_cheby_iters, &
                           rx, &
                           ry, &
                           step)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: w
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: p, r, Mi, z
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: Kx, Ky

  INTEGER :: j,k

    REAL(KIND=8) ::  rx, ry

    INTEGER :: step, max_cheby_iters
    REAL(KIND=8), DIMENSION(max_cheby_iters) :: ch_alphas, ch_betas

!$OMP PARALLEL
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            w(j, k) = (1.0_8                                      &
                + ry*(Ky(j, k+1) + Ky(j, k))                      &
                + rx*(Kx(j+1, k) + Kx(j, k)))*u(j, k)             &
                - ry*(Ky(j, k+1)*u(j, k+1) + Ky(j, k)*u(j, k-1))  &
                - rx*(Kx(j+1, k)*u(j+1, k) + Kx(j, k)*u(j-1, k))
            r(j, k) = u0(j, k) - w(j, k)
            !z(j, k) = Mi(j, k)*r(j, k)
            p(j, k) = ch_alphas(step)*p(j, k) + ch_betas(step)*Mi(j, k)*r(j, k)
        ENDDO
    ENDDO
!$OMP END DO
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            u(j, k) = u(j, k) + p(j, k)
        ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_cheby_iterate

SUBROUTINE tea_leaf_kernel_cheby_copy_u(x_min,             &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           u0, u)
  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u, u0
  INTEGER(KIND=4) :: j,k

!$OMP PARALLEL
!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            u0(j, k) = u(j, k)
        ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

end SUBROUTINE

SUBROUTINE tea_leaf_kernel_cheby_reset_Mi(x_min,             &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           p,           & ! 1
                           r,           & ! 2
                           Mi,          & ! 3
                           z,           & ! 5
                           rro)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: p
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: r
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: Mi
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: z

  INTEGER(KIND=4) :: j,k

  REAL(kind=8) :: rro

!$OMP PARALLEL
!$OMP DO REDUCTION(+:rro)
    DO k=y_min,y_max
        DO j=x_min,x_max
            Mi(j, k) = 1.0_8

            z(j, k) = r(j, k)
            p(j, k) = r(j, k)

            rro = rro + r(j, k)*r(j, k);
        ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

end SUBROUTINE

SUBROUTINE tqli(d,e,n, info)
    ! http://physics.sharif.edu/~jafari/fortran-codes/lanczos/tqli.f90
    IMPLICIT NONE
    REAL(KIND=8), DIMENSION(n) :: d,e
    INTEGER :: i,iter,l,m,n,info
    REAL(KIND=8) :: b,c,dd,f,g,p,r,s
    e(:)=eoshift(e(:),1)
    info = 0
    do l=1,n
        iter=0
        iterate: do
            do m=l,n-1
                dd=abs(d(m))+abs(d(m+1))
                if (abs(e(m))+dd == dd) exit
            end do
            if (m == l) exit iterate
            if (iter == 30) then
                info=1
                return
            endif
            iter=iter+1
            g=(d(l+1)-d(l))/(2.0_8*e(l))
            r=hypot(g,1.0_8)
            g=d(m)-d(l)+e(l)/(g+sign(r,g))
            s=1.0_8
            c=1.0_8
            p=0.0_8
            do i=m-1,l,-1
                f=s*e(i)
                b=c*e(i)
                r=hypot(f,g)
                e(i+1)=r
                if (r == 0.0_8) then
                    d(i+1)=d(i+1)-p
                    e(m)=0.0_8
                    cycle iterate
                end if
                s=f/r
                c=g/r
                g=d(i+1)-p
                r=(d(i)-g)*s+2.0_8*c*b
                p=s*r
                d(i+1)=g+p
                g=c*r-b
            end do
            d(l)=d(l)-p
            e(l)=g
            e(m)=0.0_8
        end do iterate
    end do
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

  do n=1,tl_ch_cg_presteps
    diag(n) = 1.0_8/cg_alphas(n)
    if (n .gt. 1) diag(n) = diag(n) + cg_betas(n-1)/cg_alphas(n-1)
    if (n .lt. tl_ch_cg_presteps) offdiag(n+1) = sqrt(cg_betas(n))/cg_alphas(n)
  enddo

  CALL tqli(diag, offdiag, tl_ch_cg_presteps, info)

  ! could just call this instead
  !offdiag(:)=eoshift(offdiag(:),1)
  !CALL dsterf(tl_ch_cg_presteps, diag, offdiag, info)

  if (info .ne. 0) return

  ! bubble sort eigenvalues
  do
    do n=1,tl_ch_cg_presteps-1
      if (diag(n) .ge. diag(n+1)) then
        tmp = diag(n)
        diag(n) = diag(n+1)
        diag(n+1) = tmp
        swapped = .true.
      endif
    enddo
    if (.not. swapped) exit
    swapped = .false.
  enddo

  eigmin = diag(1)
  eigmax = diag(tl_ch_cg_presteps)

  if (eigmin .lt. 0.0_8 .or. eigmax .lt. 0.0_8) info = 1

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

  do n=1,max_cheby_iters
    rho_new = 1.0_8/(2.0_8*sigma - rho_old)
    cur_alpha = rho_new*rho_old
    cur_beta = 2.0_8*rho_new/delta

    ch_alphas(n) = cur_alpha
    ch_betas(n) = cur_beta

    rho_old = rho_new
  enddo

END SUBROUTINE tea_calc_ch_coefs

end module

