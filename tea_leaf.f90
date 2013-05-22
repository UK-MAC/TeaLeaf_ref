MODULE tea_leaf_kernel_module

CONTAINS

SUBROUTINE tea_leaf_kernel_init(x_min,             &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           celldx,            &
                           celldy,            &
                           volume,            &
                           density,           &
                           energy,            &
                           u0,                &
                           u1,                &
                           un,                &
                           heat_capacity,     &
                           Kx_tmp,             &
                           Ky_tmp,             &
                           Kx,                &
                           Ky)

  USE clover_module
  USE report_module

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: celldx
  REAL(KIND=8), DIMENSION(y_min-2:y_max+2) :: celldy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: volume
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: u0
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: heat_capacity
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Kx_tmp
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Ky_tmp
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Kx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Ky
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: un

  INTEGER(KIND=4), DIMENSION(2) :: error_loc

  INTEGER(KIND=4) :: j,k,n


! CALC DIFFUSION COEFFICIENT
!$OMP PARALLEL
  IF(coefficient .EQ. RECIP_CONDUCTIVITY) THEN
!$OMP DO 
    DO k=y_min-1,y_max+1
      DO j=x_min-1,x_max+1
         Kx_tmp(j  ,k  )=1.0_8/density(j  ,k  )
         Kx_tmp(j+1,k  )=1.0_8/density(j+1,k  )
         Ky_tmp(j  ,k  )=1.0_8/density(j  ,k  )
         Ky_tmp(j  ,k+1)=1.0_8/density(j  ,k+1)
      ENDDO
    ENDDO
!$OMP END DO
  ELSE IF(coefficient .EQ. CONDUCTIVITY) THEN
!$OMP DO
    DO k=y_min-1,y_max+1
      DO j=x_min-1,x_max+1
         Kx_tmp(j  ,k  )=density(j  ,k  )
         Kx_tmp(j+1,k  )=density(j+1,k  )
         Ky_tmp(j  ,k  )=density(j  ,k  )
         Ky_tmp(j  ,k+1)=density(j  ,k+1)
      ENDDO
    ENDDO
!$OMP END DO
  ELSE
    CALL report_error('tea_leaf', 'unknown coefficient option')
  ENDIF

!$OMP DO
  DO k=y_min-1,y_max+1
    DO j=x_min-1,x_max+1
         Kx(j,k)=(Kx_tmp(j  ,k  )+Kx_tmp(j+1,k  ))/(2.0*Kx_tmp(j  ,k  )*Kx_tmp(j+1,k  ))
         Ky(j,k)=(Ky_tmp(j  ,k  )+Ky_tmp(j  ,k+1))/(2.0*Ky_tmp(j  ,k  )*Ky_tmp(j  ,k+1))
    ENDDO
  ENDDO
!$OMP END DO

!$OMP DO 
  DO k=y_min-1, y_max+1
    DO j=x_min-1, x_max+1
      u0(j,k) =  energy(j,k) * density(j,k)
    ENDDO
  ENDDO
!$OMP END DO

  ! INITIAL GUESS
!$OMP DO
  DO k=y_min-1, y_max+1
    DO j=x_min-1, x_max+1
      u1(j,k) = u0(j,k)
    ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL


END SUBROUTINE tea_leaf_kernel_init

SUBROUTINE tea_leaf_kernel_solve(x_min,       &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           rx,                &
                           ry,                &
                           Kx,                &
                           Ky,                &
                           error,             &
                           u0,                &
                           u1,                &
                           un)


  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: u0, un
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u1
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Kx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Ky

  REAL(KIND=8) :: ry,rx, error

  INTEGER(KIND=4) :: j,k

  error = 0.0_8

!$OMP PARALLEL
!$OMP DO
    DO k=y_min-1, y_max+1
      DO j=x_min-1, x_max+1
        un(j,k) = u1(j,k)
      ENDDO
    ENDDO
!$OMP END DO

!$OMP DO REDUCTION(MAX:error)
    DO k=y_min, y_max
      DO j=x_min, x_max
        u1(j,k) = (u0(j,k) + Kx(j+1,k)*rx*un(j+1,k) + Kx(j,k)*rx*un(j-1,k) &
                           + Ky(j,k+1)*ry*un(j,k+1) + Ky(j,k)*ry*un(j,k-1)) &
                             /(1+2.0_8*rx+2.0_8*ry)

        error = MAX(error, ABS(u1(j,k)-u0(j,k)))
      ENDDO
    ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE tea_leaf_kernel_solve

SUBROUTINE tea_leaf_kernel_finalise(x_min,    &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           energy,            &
                           density,           &
                           u)

  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: u
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density

  INTEGER(KIND=4) :: j,k

!$OMP PARALLEL DO 
  DO k=y_min, y_max
    DO j=x_min, x_max
      energy(j,k) = u(j,k) / density(j,k)
    ENDDO
  ENDDO
!$OMP END PARALLEL DO

END SUBROUTINE tea_leaf_kernel_finalise

END MODULE tea_leaf_kernel_module

MODULE tea_leaf_module

CONTAINS

SUBROUTINE tea_leaf()
 
  USE clover_module
  USE tea_leaf_kernel_module
  USE update_halo_module

  IMPLICIT NONE

  INTEGER :: c, n, j,k
  REAL(KIND=8) :: ry,rx, error

  INTEGER :: fields(NUM_FIELDS)


  DO c=1,number_of_chunks

    IF(chunks(c)%task.EQ.parallel%task) THEN

      IF(use_fortran_kernels) THEN

            fields=0
            fields(FIELD_ENERGY1) = 1
            fields(FIELD_DENSITY1) = 1
            CALL update_halo(fields,2)

          ! INIT
          CALL tea_leaf_kernel_init(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                       &
              chunks(c)%field%y_min,                       &
              chunks(c)%field%y_max,                       &
              chunks(c)%field%celldx,                      &
              chunks(c)%field%celldy,                      &
              chunks(c)%field%volume,                      &
              chunks(c)%field%density1,                    &
              chunks(c)%field%energy1,                     &
              chunks(c)%field%work_array1,                 &
              chunks(c)%field%u,                           &
              chunks(c)%field%work_array3,                 &
              chunks(c)%field%work_array4,                 &
              chunks(c)%field%work_array5,                 &
              chunks(c)%field%work_array6,                 &
              chunks(c)%field%work_array7,                 &
              chunks(c)%field%work_array8)


          ! JACOBI

          rx = dt/(chunks(c)%field%celldx(chunks(c)%field%x_min)**2);
          ry = dt/(chunks(c)%field%celldy(chunks(c)%field%y_min)**2);

          DO n=1,max_iters

            CALL tea_leaf_kernel_solve(chunks(c)%field%x_min,&
                chunks(c)%field%x_max,                       &
                chunks(c)%field%y_min,                       &
                chunks(c)%field%y_max,                       &
                rx,                                          &
                ry,                                          &
                chunks(c)%field%work_array7,                 &
                chunks(c)%field%work_array8,                 &
                error,                                       &
                chunks(c)%field%work_array1,                 &
                chunks(c)%field%u,                           &
                chunks(c)%field%work_array3)

            ! CALL update_halo
            fields=0
            fields(FIELD_U) = 1
            CALL update_halo(fields,2)

            CALL clover_max(error)

            !PRINT *, 'ERR: ', error
            IF (error .LT. eps) EXIT
          ENDDO

          ! RESET
          CALL tea_leaf_kernel_finalise(chunks(c)%field%x_min, &
              chunks(c)%field%x_max,                           &
              chunks(c)%field%y_min,                           &
              chunks(c)%field%y_max,                           &
              chunks(c)%field%energy1,                         &
              chunks(c)%field%density1,                        &
              chunks(c)%field%u)

            fields=0
            fields(FIELD_ENERGY1) = 1
            CALL update_halo(fields,1)

      ENDIF

    ENDIF

  ENDDO

END SUBROUTINE tea_leaf

END MODULE tea_leaf_module
