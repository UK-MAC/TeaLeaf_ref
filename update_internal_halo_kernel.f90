
MODULE update_internal_halo_kernel_module

  ! These need to be kept consistent with the data module to avoid use statement
  INTEGER,private,PARAMETER :: CHUNK_LEFT   =1    &
                              ,CHUNK_RIGHT  =2    &
                              ,CHUNK_BOTTOM =3    &
                              ,CHUNK_TOP    =4    &
                              ,EXTERNAL_FACE=-1

  INTEGER,private,PARAMETER   :: FIELD_DENSITY    = 1         &
                                ,FIELD_ENERGY0    = 2         &
                                ,FIELD_ENERGY1    = 3         &
                                ,FIELD_U          = 4         &
                                ,FIELD_P          = 5         &
                                ,FIELD_SD         = 6         &
                                ,FIELD_R          = 7         &
                                ,FIELD_Z          = 8         &
                                ,FIELD_KX         = 9         &
                                ,FIELD_KY         = 10        &
                                ,FIELD_DI         = 11        &
                                ,NUM_FIELDS       = 11

CONTAINS

  SUBROUTINE update_internal_halo_left_right_kernel(                                &
                          x_min,x_max,y_min,y_max,                                    &
                          density,                                                    &
                          energy0,                                                    &
                          energy1,                                                    &
                          u,                                                          &
                          p,                                                          &
                          sd,                                                         &
                          r,                                                          &
                          z,                                                          &
                          kx,                                                         &
                          ky,                                                         &
                          di,                                                         &
                          x_min_right,x_max_right,y_min_right,y_max_right,            &
                          density_right,                                              &
                          energy0_right,                                              &
                          energy1_right,                                              &
                          u_right,                                                    &
                          p_right,                                                    &
                          sd_right,                                                   &
                          r_right,                                                    &
                          z_right,                                                    &
                          kx_right,                                                   &
                          ky_right,                                                   &
                          di_right,                                                   &
                          halo_exchange_depth,                                        &
                          fields,                                                     &
                          depth                                                       )
    IMPLICIT NONE

    INTEGER :: halo_exchange_depth
    INTEGER :: fields(NUM_FIELDS),depth

    INTEGER :: x_min,x_max,y_min,y_max
    REAL(KIND=8),DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,&
                           y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                           :: density,energy0,energy1, u, sd, p, r, z, kx, ky, di

    INTEGER :: x_min_right,x_max_right,y_min_right,y_max_right
    REAL(KIND=8),DIMENSION(x_min_right-halo_exchange_depth:x_max_right+halo_exchange_depth,y_min_right-halo_exchange_depth:&
                           y_max_right+halo_exchange_depth) :: density_right,energy0_right,energy1_right, &
                                                               u_right,sd_right,p_right,r_right,z_right, &
                                                               kx_right,ky_right,di_right

!$OMP PARALLEL
    IF (fields(FIELD_DENSITY).EQ.1) THEN
      CALL update_internal_halo_cell_left_right(x_min, x_max, y_min, y_max, density, &
        x_min_right, x_max_right, y_min_right, y_max_right, density_right, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_ENERGY0).EQ.1) THEN
      CALL update_internal_halo_cell_left_right(x_min, x_max, y_min, y_max, energy0, &
        x_min_right, x_max_right, y_min_right, y_max_right, energy0_right, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_ENERGY1).EQ.1) THEN
      CALL update_internal_halo_cell_left_right(x_min, x_max, y_min, y_max, energy1, &
        x_min_right, x_max_right, y_min_right, y_max_right, energy1_right, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_U).EQ.1) THEN
      CALL update_internal_halo_cell_left_right(x_min, x_max, y_min, y_max, u, &
        x_min_right, x_max_right, y_min_right, y_max_right, u_right, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_p).EQ.1) THEN
      CALL update_internal_halo_cell_left_right(x_min, x_max, y_min, y_max, p, &
        x_min_right, x_max_right, y_min_right, y_max_right, p_right, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_sd).EQ.1) THEN
      CALL update_internal_halo_cell_left_right(x_min, x_max, y_min, y_max, sd, &
        x_min_right, x_max_right, y_min_right, y_max_right, sd_right, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_r).EQ.1) THEN
      CALL update_internal_halo_cell_left_right(x_min, x_max, y_min, y_max, r, &
        x_min_right, x_max_right, y_min_right, y_max_right, r_right, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_z).EQ.1) THEN
      CALL update_internal_halo_cell_left_right(x_min, x_max, y_min, y_max, z, &
        x_min_right, x_max_right, y_min_right, y_max_right, z_right, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_kx).EQ.1) THEN
      CALL update_internal_halo_cell_left_right(x_min, x_max, y_min, y_max, kx, &
        x_min_right, x_max_right, y_min_right, y_max_right, kx_right, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_ky).EQ.1) THEN
      CALL update_internal_halo_cell_left_right(x_min, x_max, y_min, y_max, ky, &
        x_min_right, x_max_right, y_min_right, y_max_right, ky_right, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_di).EQ.1) THEN
      CALL update_internal_halo_cell_left_right(x_min, x_max, y_min, y_max, di, &
        x_min_right, x_max_right, y_min_right, y_max_right, di_right, &
        halo_exchange_depth, depth)
    ENDIF
!$OMP END PARALLEL

  END SUBROUTINE

  SUBROUTINE update_internal_halo_cell_left_right(x_min_left,x_max_left,y_min_left,y_max_left,  &
                          mesh_left,   &
                          x_min_right,x_max_right,y_min_right,y_max_right,            &
                          mesh_right,                           &
                          halo_exchange_depth, &
                          depth                           )
    IMPLICIT NONE

    INTEGER :: halo_exchange_depth, depth

    INTEGER :: x_min_left,x_max_left,y_min_left,y_max_left
    REAL(KIND=8), DIMENSION(x_min_left-halo_exchange_depth:x_max_left+halo_exchange_depth,y_min_left-halo_exchange_depth:&
                            y_max_left+halo_exchange_depth) :: mesh_left

    INTEGER :: x_min_right,x_max_right,y_min_right,y_max_right
    REAL(KIND=8), DIMENSION(x_min_right-halo_exchange_depth:x_max_right+halo_exchange_depth,y_min_right-halo_exchange_depth:&
    y_max_right+halo_exchange_depth) :: mesh_right

    INTEGER :: j,k

!$OMP DO
    DO k=y_min_left-depth,y_max_left+depth
      DO j=1,depth
        mesh_right(1-j,k)=mesh_left(x_max_left+1-j,k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP DO
    DO k=y_min_left-depth,y_max_left+depth
      DO j=1,depth
        mesh_left(x_max_left+j,k)=mesh_right(0+j,k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT

  END SUBROUTINE update_internal_halo_cell_left_right

  SUBROUTINE update_internal_halo_bottom_top_kernel(                                &
                          x_min,x_max,y_min,y_max,                                    &
                          density,                                                    &
                          energy0,                                                    &
                          energy1,                                                    &
                          u,                                                          &
                          p,                                                          &
                          sd,                                                         &
                          r,                                                          &
                          z,                                                          &
                          kx,                                                         &
                          ky,                                                         &
                          di,                                                         &
                          x_min_top,x_max_top,y_min_top,y_max_top,            &
                          density_top,                                              &
                          energy0_top,                                              &
                          energy1_top,                                              &
                          u_top,                                                    &
                          p_top,                                                    &
                          sd_top,                                                   &
                          r_top,                                                    &
                          z_top,                                                    &
                          kx_top,                                                   &
                          ky_top,                                                   &
                          di_top,                                                   &
                          halo_exchange_depth,                                        &
                          fields,                                                     &
                          depth                                                       )
    IMPLICIT NONE

    INTEGER :: halo_exchange_depth
    INTEGER :: fields(NUM_FIELDS),depth

    INTEGER :: x_min,x_max,y_min,y_max
    REAL(KIND=8),DIMENSION(x_min-halo_exchange_depth:x_max+halo_exchange_depth,y_min-halo_exchange_depth:y_max+halo_exchange_depth)&
                           :: density,energy0,energy1, u, sd, p, r, z, kx, ky, di

    INTEGER :: x_min_top,x_max_top,y_min_top,y_max_top
    REAL(KIND=8), DIMENSION(x_min_top-halo_exchange_depth:x_max_top+halo_exchange_depth,y_min_top-halo_exchange_depth:&
                            y_max_top+halo_exchange_depth) :: density_top,energy0_top,energy1_top,&
                                                              u_top, sd_top, p_top, r_top, z_top, kx_top, ky_top, di_top

!$OMP PARALLEL
    IF (fields(FIELD_DENSITY).EQ.1) THEN
      CALL update_internal_halo_cell_bottom_top(x_min, x_max, y_min, y_max, density, &
        x_min_top, x_max_top, y_min_top, y_max_top, density_top, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_ENERGY0).EQ.1) THEN
      CALL update_internal_halo_cell_bottom_top(x_min, x_max, y_min, y_max, energy0, &
        x_min_top, x_max_top, y_min_top, y_max_top, energy0_top, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_ENERGY1).EQ.1) THEN
      CALL update_internal_halo_cell_bottom_top(x_min, x_max, y_min, y_max, energy1, &
        x_min_top, x_max_top, y_min_top, y_max_top, energy1_top, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_U).EQ.1) THEN
      CALL update_internal_halo_cell_bottom_top(x_min, x_max, y_min, y_max, u, &
        x_min_top, x_max_top, y_min_top, y_max_top, u_top, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_p).EQ.1) THEN
      CALL update_internal_halo_cell_bottom_top(x_min, x_max, y_min, y_max, p, &
        x_min_top, x_max_top, y_min_top, y_max_top, p_top, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_sd).EQ.1) THEN
      CALL update_internal_halo_cell_bottom_top(x_min, x_max, y_min, y_max, sd, &
        x_min_top, x_max_top, y_min_top, y_max_top, sd_top, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_r).EQ.1) THEN
      CALL update_internal_halo_cell_bottom_top(x_min, x_max, y_min, y_max, r, &
        x_min_top, x_max_top, y_min_top, y_max_top, r_top, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_z).EQ.1) THEN
      CALL update_internal_halo_cell_bottom_top(x_min, x_max, y_min, y_max, z, &
        x_min_top, x_max_top, y_min_top, y_max_top, z_top, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_kx).EQ.1) THEN
      CALL update_internal_halo_cell_bottom_top(x_min, x_max, y_min, y_max, kx, &
        x_min_top, x_max_top, y_min_top, y_max_top, kx_top, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_ky).EQ.1) THEN
      CALL update_internal_halo_cell_bottom_top(x_min, x_max, y_min, y_max, ky, &
        x_min_top, x_max_top, y_min_top, y_max_top, ky_top, &
        halo_exchange_depth, depth)
    ENDIF

    IF (fields(FIELD_di).EQ.1) THEN
      CALL update_internal_halo_cell_bottom_top(x_min, x_max, y_min, y_max, di, &
        x_min_top, x_max_top, y_min_top, y_max_top, di_top, &
        halo_exchange_depth, depth)
    ENDIF
!$OMP END PARALLEL

  END SUBROUTINE update_internal_halo_bottom_top_kernel

  SUBROUTINE update_internal_halo_cell_bottom_top(x_min_bottom,x_max_bottom,y_min_bottom,y_max_bottom,  &
                          mesh_bottom,   &
                          x_min_top,x_max_top,y_min_top,y_max_top,            &
                          mesh_top,                           &
                          halo_exchange_depth, &
                          depth                           )
    IMPLICIT NONE

    INTEGER :: halo_exchange_depth, depth

    INTEGER :: x_min_bottom,x_max_bottom,y_min_bottom,y_max_bottom
    REAL(KIND=8), DIMENSION(x_min_bottom-halo_exchange_depth:x_max_bottom+halo_exchange_depth,y_min_bottom-halo_exchange_depth:&
                            y_max_bottom+halo_exchange_depth) :: mesh_bottom
    INTEGER :: x_min_top,x_max_top,y_min_top,y_max_top
    REAL(KIND=8), DIMENSION(x_min_top-halo_exchange_depth:x_max_top+halo_exchange_depth,y_min_top-halo_exchange_depth:&
                            y_max_top+halo_exchange_depth) :: mesh_top
    INTEGER :: j,k

!$OMP DO
    DO k=1,depth
      DO j=x_min_bottom-depth,x_max_bottom+depth
        mesh_top(j,1-k)=mesh_bottom(j,y_max_bottom+1-k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT
!$OMP DO
    DO k=1,depth
      DO j=x_min_bottom-depth,x_max_bottom+depth
        mesh_bottom(j,y_max_bottom+k)=mesh_top(j,0+k)
      ENDDO
    ENDDO
!$OMP END DO NOWAIT

  END SUBROUTINE update_internal_halo_cell_bottom_top

END MODULE update_internal_halo_kernel_module

