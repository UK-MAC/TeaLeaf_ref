MODULE initialise_chunk_kernel_module

CONTAINS

SUBROUTINE initialise_chunk_kernel(x_min,x_max,y_min,y_max,       &
                                   xmin,ymin,dx,dy,               &
                                   vertexx,                       &
                                   vertexdx,                      &
                                   vertexy,                       &
                                   vertexdy,                      &
                                   cellx,                         &
                                   celldx,                        &
                                   celly,                         &
                                   celldy,                        &
                                   volume,                        &
                                   xarea,                         &
                                   yarea                          )

  IMPLICIT NONE

  INTEGER      :: x_min,x_max,y_min,y_max
  REAL(KIND=8) :: xmin,ymin,dx,dy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3) :: vertexx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3) :: vertexdx
  REAL(KIND=8), DIMENSION(y_min-2:y_max+3) :: vertexy
  REAL(KIND=8), DIMENSION(y_min-2:y_max+3) :: vertexdy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: cellx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: celldx
  REAL(KIND=8), DIMENSION(y_min-2:y_max+2) :: celly
  REAL(KIND=8), DIMENSION(y_min-2:y_max+2) :: celldy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+2) :: volume
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3 ,y_min-2:y_max+2) :: xarea
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2 ,y_min-2:y_max+3) :: yarea

  INTEGER      :: j,k

!$OMP PARALLEL
!$OMP DO
  DO j=x_min-2,x_max+3
     vertexx(j)=xmin+dx*float(j-x_min)
  ENDDO
!$OMP END DO

!$OMP DO
  DO j=x_min-2,x_max+3
    vertexdx(j)=dx
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min-2,y_max+3
     vertexy(k)=ymin+dy*float(k-y_min)
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min-2,y_max+3
    vertexdy(k)=dy
  ENDDO
!$OMP END DO

!$OMP DO
  DO j=x_min-2,x_max+2
     cellx(j)=0.5*(vertexx(j)+vertexx(j+1))
  ENDDO
!$OMP END DO

!$OMP DO
  DO j=x_min-2,x_max+2
     celldx(j)=dx
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min-2,y_max+2
     celly(k)=0.5*(vertexy(k)+vertexy(k+1))
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min-2,y_max+2
     celldy(k)=dy
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min-2,y_max+2
    DO j=x_min-2,x_max+2
        volume(j,k)=dx*dy
     ENDDO
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min-2,y_max+2
    DO j=x_min-2,x_max+2
        xarea(j,k)=celldy(k)
     ENDDO
  ENDDO
!$OMP END DO

!$OMP DO
  DO k=y_min-2,y_max+2
    DO j=x_min-2,x_max+2
        yarea(j,k)=celldx(j)
     ENDDO
  ENDDO
!$OMP END DO
!$OMP END PARALLEL

END SUBROUTINE initialise_chunk_kernel

END MODULE initialise_chunk_kernel_module
