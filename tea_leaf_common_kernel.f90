MODULE tea_leaf_common_kernel_module

  IMPLICIT NONE

  ! These need to be kept consistent with the data module to avoid use statement
  INTEGER,private,PARAMETER :: CHUNK_LEFT   =1    &
                            ,CHUNK_RIGHT  =2    &
                            ,CHUNK_BOTTOM =3    &
                            ,CHUNK_TOP    =4    &
                            ,EXTERNAL_FACE=-1

   INTEGER,PARAMETER        ::   TL_PREC_NONE       = 1 &
                                ,TL_PREC_JAC_DIAG   = 2 &
                                ,TL_PREC_JAC_BLOCK  = 3

   INTEGER,PRIVATE         ::    CONDUCTIVITY        = 1 &
                                ,RECIP_CONDUCTIVITY  = 2

  integer, private, parameter:: jac_block_size = 4

  INTEGER(KIND=4), parameter :: block_size=1
  INTEGER(KIND=4), parameter :: kstep = block_size*jac_block_size

CONTAINS

SUBROUTINE tea_leaf_common_init_kernel(x_min,  &
                           x_max,                  &
                           y_min,                  &
                           y_max,                  &
                           halo_exchange_depth,                  &
                           chunk_neighbours,       &
                           zero_boundary,       &
                           reflective_boundary,    &
                           density,                &
                           energy,                 &
                           u,                      &
                           u0,                     &
                           r,                      &
                           w,                      &
                           Kx,                     &
                           Ky,                     &
                           cp,                     &
                           bfp,                    &
                           Mi,                     &
                           rx,                     &
                           ry,                     &
                           preconditioner_type,      &
                           coef)

  IMPLICIT NONE

  LOGICAL :: reflective_boundary
  INTEGER :: preconditioner_type
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  INTEGER, DIMENSION(4) :: chunk_neighbours, zero_boundary
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: density, energy, u, r, w, Kx, Ky, Mi, u0
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp

  INTEGER(KIND=4) :: coef
  INTEGER(KIND=4) :: j,k

  REAL(KIND=8) ::  rx, ry

!$OMP PARALLEL NUM_THREADS(INNER_NUM_THREADS)
!$OMP DO
  DO k=y_min, y_max
    DO j=x_min, x_max
      u(j,k) = energy(j,k)*density(j,k)
      u0(j,k) = energy(j,k)*density(j,k)
    ENDDO
  ENDDO
!$OMP END DO

  IF (coef .EQ. RECIP_CONDUCTIVITY) THEN
!$OMP DO
    ! use w as temp val
    DO k=y_min-halo_exchange_depth,y_max+halo_exchange_depth
      DO j=x_min-halo_exchange_depth,x_max+halo_exchange_depth
         w(j  ,k  )=1.0_8/density(j  ,k  )
      ENDDO
    ENDDO
!$OMP END DO
  ELSE IF (coef .EQ. CONDUCTIVITY) THEN
!$OMP DO
    DO k=y_min-halo_exchange_depth,y_max+halo_exchange_depth
      DO j=x_min-halo_exchange_depth,x_max+halo_exchange_depth
         w(j  ,k  )=density(j  ,k  )
      ENDDO
    ENDDO
!$OMP END DO
  ENDIF

!$OMP DO
   DO k=y_min-halo_exchange_depth + 1,y_max+halo_exchange_depth
     DO j=x_min-halo_exchange_depth + 1,x_max+halo_exchange_depth
          Kx(j,k)=(w(j-1,k  ) + w(j,k))/(2.0_8*w(j-1,k  )*w(j,k))
          Ky(j,k)=(w(j  ,k-1) + w(j,k))/(2.0_8*w(j  ,k-1)*w(j,k))
     ENDDO
   ENDDO
!$OMP END DO

! Whether to apply reflective boundary conditions to all external faces
  IF (reflective_boundary .EQV. .FALSE.) THEN
    IF (chunk_neighbours(CHUNK_LEFT).EQ.EXTERNAL_FACE .AND. zero_boundary(CHUNK_LEFT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-halo_exchange_depth,y_max+halo_exchange_depth
        DO j=x_min-halo_exchange_depth,x_min
          Kx(j,k)=0.0_8
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF (chunk_neighbours(CHUNK_RIGHT).EQ.EXTERNAL_FACE .AND. zero_boundary(CHUNK_RIGHT).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-halo_exchange_depth,y_max+halo_exchange_depth
        DO j=x_max + 1,x_max+halo_exchange_depth
          Kx(j,k)=0.0_8
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF (chunk_neighbours(CHUNK_BOTTOM).EQ.EXTERNAL_FACE .AND. zero_boundary(CHUNK_BOTTOM).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_min-halo_exchange_depth,y_min
        DO j=x_min-halo_exchange_depth,x_max+halo_exchange_depth
          Ky(j,k)=0.0_8
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
    IF (chunk_neighbours(CHUNK_TOP).EQ.EXTERNAL_FACE .AND. zero_boundary(CHUNK_TOP).EQ.EXTERNAL_FACE) THEN
!$OMP DO
      DO k=y_max + 1,y_max+halo_exchange_depth
        DO j=x_min-halo_exchange_depth,x_max+halo_exchange_depth
          Ky(j,k)=0.0_8
        ENDDO
      ENDDO
!$OMP END DO
    ENDIF
  ENDIF

  IF (preconditioner_type .EQ. TL_PREC_JAC_BLOCK) THEN
    CALL tea_block_init(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                           cp, bfp, Kx, Ky, rx, ry)
  ELSE IF (preconditioner_type .EQ. TL_PREC_JAC_DIAG) THEN
    CALL tea_diag_init(x_min, x_max, y_min, y_max, halo_exchange_depth,             &
                           Mi, Kx, Ky, rx, ry)
  ENDIF

!$OMP DO
    DO k=y_min,y_max
        DO j=x_min,x_max
            w(j, k) = (1.0_8                                      &
                + ry*(Ky(j, k+1) + Ky(j, k))                      &
                + rx*(Kx(j+1, k) + Kx(j, k)))*u(j, k)             &
                - ry*(Ky(j, k+1)*u(j, k+1) + Ky(j, k)*u(j, k-1))  &
                - rx*(Kx(j+1, k)*u(j+1, k) + Kx(j, k)*u(j-1, k))

            r(j, k) = u(j, k) - w(j, k)
            !r(j, k) = u(j, k)! This is required to make a zero initial guess to match petsc errant behaviour
                              ! Only works one timestep is run
        ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_common_init_kernel

SUBROUTINE tea_leaf_kernel_finalise(x_min,    &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           halo_exchange_depth,             &
                           energy,            &
                           density,           &
                           u)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: u, energy, density

  INTEGER(KIND=4) :: j,k

!$OMP PARALLEL NUM_THREADS(INNER_NUM_THREADS)
!$OMP DO
  DO k=y_min, y_max
    DO j=x_min, x_max
      energy(j,k) = u(j,k) / density(j,k)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_finalise

SUBROUTINE tea_leaf_calc_residual_kernel(x_min,       &
                                  x_max,       &
                                  y_min,       &
                                  y_max,       &
                                  halo_exchange_depth,       &
                                  u ,          &
                                  u0,          &
                                  r,           &
                                  Kx,          &
                                  Ky,          &
                                  rx, ry       )

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: Kx, u, r, Ky, u0

  REAL(KIND=8) :: smvp, rx, ry

  INTEGER(KIND=4) :: j,k

!$OMP PARALLEL NUM_THREADS(INNER_NUM_THREADS) PRIVATE(smvp)
!$OMP DO
    DO k=y_min, y_max
      DO j=x_min, x_max
        smvp = (1.0_8                                         &
            + ry*(Ky(j, k+1) + Ky(j, k))                      &
            + rx*(Kx(j+1, k) + Kx(j, k)))*u(j, k)             &
            - ry*(Ky(j, k+1)*u(j, k+1) + Ky(j, k)*u(j, k-1))  &
            - rx*(Kx(j+1, k)*u(j+1, k) + Kx(j, k)*u(j-1, k))
        r(j, k) = u0(j, k) - smvp
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_calc_residual_kernel

SUBROUTINE tea_leaf_calc_2norm_kernel(x_min, &
                          x_max,             &
                          y_min,             &
                          y_max,             &
                          halo_exchange_depth,             &
                          arr,               &
                          norm)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: arr
  REAL(KIND=8) :: norm
  integer :: j, k

  norm = 0.0_8

!$OMP PARALLEL NUM_THREADS(INNER_NUM_THREADS)
!$OMP DO REDUCTION(+:norm)
    DO k=y_min,y_max
        DO j=x_min,x_max
            norm = norm + arr(j, k)*arr(j, k)
        ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_calc_2norm_kernel

SUBROUTINE tea_diag_init(x_min,             &
                         x_max,             &
                         y_min,             &
                         y_max,             &
                         halo_exchange_depth,             &
                         Mi,                &
                         Kx, Ky, rx, ry)

  IMPLICIT NONE

  INTEGER(KIND=4):: j, k
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: Kx, Ky, Mi
  REAL(KIND=8) :: rx, ry

!$OMP DO
    DO k=y_min-halo_exchange_depth,y_max+halo_exchange_depth
      DO j=x_min-halo_exchange_depth,x_max+halo_exchange_depth
        Mi(j, k) = 1.0_8/(1.0_8                 &
                + ry*(Ky(j, k+1) + Ky(j, k))    &
                + rx*(Kx(j+1, k) + Kx(j, k)))
      ENDDO
    ENDDO
!$OMP END DO

END SUBROUTINE

SUBROUTINE tea_diag_solve(x_min,             &
                         x_max,             &
                         y_min,             &
                         y_max,             &
                         halo_exchange_depth,             &
                         r,                 &
                         z,                 &
                         Mi,                &
                         Kx, Ky, rx, ry)

  IMPLICIT NONE

  INTEGER(KIND=4):: j, k
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: Kx, Ky, r, z, Mi
  REAL(KIND=8) :: rx, ry

!$OMP DO
    DO k=y_min-halo_exchange_depth,y_max+halo_exchange_depth
      DO j=x_min-halo_exchange_depth,x_max+halo_exchange_depth
        z(j, k) = Mi(j, k)*r(j, k)
      ENDDO
    ENDDO
!$OMP END DO

END SUBROUTINE

#define COEF_A (-Ky(j, k)*ry)
#define COEF_B (1.0_8 + ry*(Ky(j, k+1) + Ky(j, k)) + rx*(Kx(j+1, k) + Kx(j, k)))
#define COEF_C (-Ky(j, k+1)*ry)

SUBROUTINE tea_block_init(x_min,             &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           halo_exchange_depth,             &
                           cp,                     &
                           bfp,                     &
                           Kx, Ky, rx, ry)

  IMPLICIT NONE

  INTEGER(KIND=4):: j, ko, k, bottom, top
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: Kx, Ky
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  REAL(KIND=8) :: rx, ry

!$OMP DO
    DO ko=y_min,y_max,jac_block_size

      bottom = ko
      top = MIN(ko + jac_block_size - 1, y_max)

#if defined(WITH_OMP4)
!$OMP SIMD
#endif
      DO j=x_min, x_max
        k = bottom
        cp(j,k) = COEF_C/COEF_B

        DO k=bottom+1,top
            bfp(j, k) = 1.0_8/(COEF_B - COEF_A*cp(j, k-1))
            cp(j, k) = COEF_C*bfp(j, k)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO

END SUBROUTINE

SUBROUTINE tea_block_solve(x_min,             &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           halo_exchange_depth,             &
                           r,                 &
                           z,                 &
                           cp,                     &
                           bfp,                     &
                           Kx, Ky, rx, ry)

  IMPLICIT NONE

  INTEGER(KIND=4):: j, ko, k, bottom, top, ki, upper_k, k_extra
  INTEGER(KIND=4):: x_min,x_max,y_min,y_max,halo_exchange_depth
  REAL(KIND=8), DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth) :: Kx, Ky, r, z
  REAL(KIND=8), DIMENSION(x_min:x_max,y_min:y_max) :: cp, bfp
  REAL(KIND=8) :: rx, ry
  REAL(KIND=8), dimension(0:jac_block_size-1) :: dp_l, z_l

  k_extra = y_max - MOD(y_max, kstep)

!$OMP DO
    DO ko=y_min, k_extra, kstep
      upper_k = ko+kstep - jac_block_size

      DO ki=ko,upper_k,jac_block_size
        bottom = ki
        top = ki+jac_block_size - 1

#if defined(WITH_OMP4)
!$OMP SIMD PRIVATE(dp_l, z_l)
#endif
        DO j=x_min,x_max
          k = bottom
          dp_l(k-bottom) = r(j, k)/COEF_B

          DO k=bottom+1,top
            dp_l(k-bottom) = (r(j, k) - COEF_A*dp_l(k-bottom-1))*bfp(j, k)
          ENDDO

          k=top
          z_l(k-bottom) = dp_l(k-bottom)

          DO k=top-1, bottom, -1
            z_l(k-bottom) = dp_l(k-bottom) - cp(j, k)*z_l(k-bottom+1)
          ENDDO

          DO k=bottom,top
            z(j, k) = z_l(k-bottom)
          ENDDO
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO NOWAIT

!$OMP DO
    DO ki=k_extra+1, y_max, jac_block_size
      bottom = MIN(ki, y_max)
      top = MIN(ki+jac_block_size-1, y_max)

#if defined(WITH_OMP4)
!$OMP SIMD PRIVATE(dp_l, z_l)
#endif
      DO j=x_min,x_max
        k = bottom
        dp_l(k-bottom) = r(j, k)/COEF_B

        DO k=bottom+1,top
          dp_l(k-bottom) = (r(j, k) - COEF_A*dp_l(k-bottom-1))*bfp(j, k)
        ENDDO

        k=top
        z_l(k-bottom) = dp_l(k-bottom)

        DO k=top-1, bottom, -1
          z_l(k-bottom) = dp_l(k-bottom) - cp(j, k)*z_l(k-bottom+1)
        ENDDO

        DO k=bottom,top
          z(j, k) = z_l(k-bottom)
        ENDDO
      ENDDO
    ENDDO
!$OMP END DO

END SUBROUTINE

END MODULE tea_leaf_common_kernel_module
