MODULE tea_leaf_kernel_module

CONTAINS

SUBROUTINE tea_leaf_kernel(x_min,             &
                           x_max,             &
                           y_min,             &
                           y_max,             &
                           dt,                &
                           celldx,            &
                           celldy,            &
                           volume,            &
                           density,           &
                           energy,            &
                           temperature,       &
                           heat_capacity,     &
                           xtemp,             &
                           ytemp,             &
                           Kx,                &
                           Ky                 )


  IMPLICIT NONE

  INTEGER(KIND=4):: x_min,x_max,y_min,y_max
  REAL(KIND=8) :: dt
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: celldx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: celldy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: volume
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: temperature
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: heat_capacity
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xtemp
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: ytemp
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Kx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Ky

  REAL(KIND=8),ALLOCATABLE    :: matrix(:,:)
  INTEGER(KIND=4),ALLOCATABLE :: columns(:,:)
  REAL(KIND=8),ALLOCATABLE    :: rhs(:)
  REAL(KIND=8),ALLOCATABLE    :: solution(:)
  REAL(KIND=8),ALLOCATABLE    :: initial_temperature(:)

  INTEGER(KIND=4):: max_task,task

  INTEGER(KIND=4):: failed

  ALLOCATE(matrix(5,(x_max-x_min+1)*(y_max-y_min+1)))
  ALLOCATE(columns(5,(x_max-x_min+1)*(y_max-y_min+1)))
  ALLOCATE(rhs((x_max-x_min+1)*(y_max-y_min+1)))
  ALLOCATE(solution((x_max-x_min+1)*(y_max-y_min+1)))
  ALLOCATE(initial_temperature((x_max-x_min+1)*(y_max-y_min+1)))

  ! Calculate heat conduction coefficients
  CALL tl_calc_K(Kx,Ky,xtemp,ytemp,density,temperature,heat_capacity,celldx,celldy,x_min,x_max,y_min,y_max)

  ! Calculate sources
  !CALL tl_calc_source()

  ! Construct the matrix
  !CALL tl_calc_matrix(max_task,matrix,columns,rhs,solution,initial_temperature, &
  !                  Kx,Ky,heat_capacity,volume,dt,             &
  !                  x_min,             &
  !                  x_max,             &
  !                  y_min,             &
  !                  y_max              )


  ! Inverts the matrix and solves
  !CALL tl_solve(max_task,task,x_min,x_max,y_min,y_max,failed,columns,matrix,rhs,temperature,solution)

  !CALL tl_clean()

  !CALL clover_check_error(failed) ! Need to be called outside
  !IF(failed.GT.0) THEN
  !  WRITE(0,*)"TeaLeaf failure"
  !ENDIF

  DEALLOCATE(matrix)
  DEALLOCATE(columns)
  DEALLOCATE(rhs)
  DEALLOCATE(solution)
  DEALLOCATE(initial_temperature)

END SUBROUTINE tea_leaf_kernel

SUBROUTINE tl_calc_K(Kx,Ky,xtemp,ytemp,density,temperature,heat_capacity,celldx,celldy,x_min,x_max,y_min,y_max)

  INTEGER(KIND=4) :: x_min,x_max,y_min,y_max
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: celldx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2) :: celldy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: volume
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: density
  REAL(KIND=8), DIMENSION(x_min-2:x_max+2,y_min-2:y_max+2) :: energy
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: temperature
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: heat_capacity
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: xtemp
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: ytemp
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Kx
  REAL(KIND=8), DIMENSION(x_min-2:x_max+3,y_min-2:y_max+3) :: Ky

  REAL(KIND=8):: xarea(x_min-2:x_max+1,y_min:y_max)
  REAL(KIND=8):: yarea(x_min-2:x_max,y_min:y_max+1)

  INTEGER(KIND=4) :: j,k

  DO k=y_min,y_max
    DO j=x_min,x_max
       heat_capacity(j,k)=0.1*density(j,k)
    ENDDO
  ENDDO

  DO k=y_min-1,y_max+1
    DO j=x_min-1,x_max+1
      xtemp(j,k) = (temperature(j,k)+temperature(j+1,k))*0.5
    ENDDO
  ENDDO

  !DO k=y_min-1,y_max+1
  !  DO j=x_min-1,x_max+1
  !    ytemp(j,k) = (temperature(j,k)+temperature(j,k+1))*0.5
  !  ENDDO
  !ENDDO

  !DO k=y_min-1,y_max+1
  !  DO j=x_min-1,x_max+1
  !     Kx(j  ,k  )=xtemp(j  ,k  )*density(j  ,k  )
  !     Kx(j+1,k  )=xtemp(j+1,k  )*density(j+1,k  )
  !     Ky(j  ,k  )=ytemp(j  ,k  )*density(j  ,k  )
  !     Ky(j  ,k+1)=ytemp(j  ,k+1)*density(j  ,k+1)
  !  ENDDO
  !ENDDO

  !DO k=y_min-1,y_max+1
  !  DO j=x_min-1,x_max+1
  !       Kx(j,k)=(Kx(j  ,k  )+Kx(j+1,k  ))/(2.0*Kx(j  ,k  )*Kx(j+1,k  ))
  !       Ky(j,k)=(Ky(j  ,k  )+Ky(j  ,k+1))/(2.0*Ky(j  ,k  )*Ky(j  ,k+1))
  !  ENDDO
  !ENDDO

  !DO k=y_min-1,y_max+1
  !  DO j=x_min-1,x_max+1
  !    dx=celldx(j)
  !    dy=celldy(k+1)

  !    IF(j.NE.jmx+1)Kx(j,k) = xarea(j+1,k)/(dx*Kx(j,k))
  !    IF(k.NE.kmx+1)Ky(j,k) = yarea(j,k+1)/(dy*Ky(j,k))
  !  ENDDO
  !ENDDO

END SUBROUTINE tl_calc_k

SUBROUTINE tl_calc_matrix(max_task,matrix,columns,rhs,solution,initial_temperature, &
                    Kx,Ky,heat_capacity,volume,dt,             &
                    x_min,             &
                    x_max,             &
                    y_min,             &
                    y_max              )

  ! Using a five point stencil, calculate the matrix coefficients

  IMPLICIT NONE

  INTEGER(4) :: max_task,x_min,x_max,y_min,y_max
  REAL(8) :: matrix(5,(x_max-x_min+1)*(y_max-y_min+1))
  INTEGER(4) :: columns(5,(x_max-x_min+1)*(y_max-y_min+1))
  REAL(8) :: rhs((x_max-x_min+1)*(y_max-y_min+1))
  REAL(8) :: solution((x_max-x_min+1)*(y_max-y_min+1))
  REAL(8) :: initial_temperature((x_max-x_min+1)*(y_max-y_min+1))
  REAL(8) :: Kx(:,:),Ky(:,:),heat_capacity(:,:),volume(:,:),dt

  REAL(8) :: temperature(x_min-3:x_max+3,y_min-3:y_max+3)

  INTEGER(4) :: j,k,n,nx,ny
  INTEGER(4) :: ncells(0:max_task),col_start,col_add,delta_col

  REAL(8) :: leftone,rightone,upone,downone,diag

  matrix=0.0
  columns=0

  n=1
  nx=x_max-x_min+1
  ny=y_max-y_min+1
  !col_start=jw+(ks*(jdome+1))
  !delta_col=jdome-(jmx-3)
  col_add=0

  DO k=y_min,y_max
    DO j=x_min,x_max

      leftone =Kx(j-1,k  )
      rightone=Kx(j  ,k  )
      downone =Ky(j  ,k-1)
      upone   =Ky(j  ,k  )

      ! Contribution to diagonal
      diag=(volume(j,k)/dt)*(heat_capacity(j,k))

      !IF(j.EQ.jmn.AND.jw.EQ.0) leftone=0.0
      !IF(k.EQ.kmn.AND.ks.EQ.0) downone=0.0
      !IF(j.EQ.jmx.AND.je.EQ.jdome) rightone=0.0
      !IF(k.EQ.kmx.AND.kn.EQ.kdomn) upone=0.0

      ! Form the matrix
      matrix(1,n) = -downone
      matrix(2,n) = -leftone
      matrix(3,n) =  diag+upone+leftone+rightone+downone
      matrix(4,n) = -rightone
      matrix(5,n) = -upone
      columns(1,n) = n-nx
      columns(2,n) = n-1
      columns(3,n) = n
      columns(4,n) = n+1
      columns(5,n) = n+nx

      columns(:,n)=columns(:,n)+col_add+col_start

      !IF(k.EQ.kmn.AND.ks.EQ.0)columns(1,n) = 0
      !IF(j.EQ.jmn.AND.jw.EQ.0)columns(2,n) = 0
      !IF(j.EQ.jmx.AND.je.EQ.jdome)columns(4,n) = 0
      !IF(k.EQ.kmx.AND.kn.EQ.kdomn)columns(5,n) = 0

      n=n+1

    ENDDO
    col_add=col_add+delta_col
  ENDDO

  n=1
  DO k=y_min,y_max
    DO j=x_min,x_max

      leftone=0.0
      rightone=0.0
      upone=0.0
      downone=0.0

      !IF(j.EQ.jmn) leftone =temperature(j-1,k  )*Kx(j-1,k  )
      !IF(j.EQ.jmx) rightone=temperature(j+1,k  )*Kx(j  ,k  )
      !IF(k.EQ.kmn) downone =temperature(j  ,k-1)*Ky(j  ,k-1)
      !IF(k.EQ.kmx) upone   =temperature(j  ,k+1)*Ky(j  ,k  )

      !IF(j.EQ.jmn.AND.jw.eq.0) leftone=0.0
      !IF(k.EQ.kmn.AND.ks.EQ.0) downone=0.0
      !IF(j.EQ.jmx.AND.je.EQ.jdome) rightone=0.0
      !IF(k.EQ.kmx.AND.kn.EQ.kdomn) upone=0.0

      diag=(volume(j,k)/dt)*heat_capacity(j,k)

      !rhs(n)=temperature(j,k)*diag+(upone+leftone+rightone+downone)+source(j,k)*volume(j,k)*density(k,k)
      initial_temperature(n)=temperature(j,k)
      n=n+1

    ENDDO
  ENDDO

  n=1
  DO k=y_min,y_max
    DO j=x_min,x_max
      solution(n)=temperature(j,k)
      n=n+1
    ENDDO
  ENDDO
 
END SUBROUTINE tl_calc_matrix

SUBROUTINE tl_solve(max_task,task,x_min,x_max,y_min,y_max,failed,columns,matrix,rhs,temperature,solution)

  IMPLICIT NONE

  INTEGER(4) :: max_task,task,x_min,x_max,y_min,y_max
  REAL(8) :: matrix(5,(x_max-x_min+1)*(y_max-y_min+1))
  INTEGER(4) :: columns(5,(x_max-x_min+1)*(y_max-y_min+1))
  REAL(8) :: rhs((x_max-x_min+1)*(y_max-y_min+1))
  REAL(8) :: solution((x_max-x_min+1)*(y_max-y_min+1))
  REAL(8) :: initial_temperature((x_max-x_min+1)*(y_max-y_min+1))
  REAL(8) :: temperature(:,:)

  INTEGER(4) :: j,k,n,procs,rank,order,local_size,delta_row,row_start,nx,ny
  INTEGER(4) :: i,l
  LOGICAL :: global,write_matrix
  INTEGER(4) :: failed

  procs=max_task+1
  rank=task
  local_size=(x_max-x_min+1)*(y_max-y_min+1)
  order= 100 !x_cells*y_cells
  global=.TRUE.

  !row_start=jw+(ks*(jdome+1))
  !delta_row=jdome-(jmx-3)
  !nx=jmx-2
  !ny=kmx-2

  !CALL tea_solve(order, local_size, columns, matrix, rhs, solution, procs, rank, row_start, delta_row, nx, ny, global,write_matrix)

  n=1
  DO k=y_min,y_max
    DO j=x_min,x_max
      IF (solution(n).LT.0.0) THEN
        failed=1
        WRITE(0,*)"Negative solution in temperature ",j,k
      ENDIF
      n=n+1
    ENDDO
  ENDDO

  n=1
  DO k=y_min,y_max
    DO j=x_min,x_max
      temperature(j,k)=solution(n)
      n=n+1
    ENDDO
  ENDDO

END SUBROUTINE tl_solve

SUBROUTINE tl_calc_source()

  IMPLICIT NONE

  INTEGER(KIND=4) :: j,k

END SUBROUTINE tl_calc_source

SUBROUTINE tl_clean()

  IMPLICIT NONE

END SUBROUTINE tl_clean

END MODULE tea_leaf_kernel_module

MODULE tea_leaf_module

CONTAINS

SUBROUTINE tea_leaf()
 
  USE clover_module
  USE tea_leaf_kernel_module

  IMPLICIT NONE

  INTEGER :: c

  DO c=1,number_of_chunks

    IF(chunks(c)%task.EQ.parallel%task) THEN

      IF(use_fortran_kernels) THEN
        CALL tea_leaf_kernel(chunks(c)%field%x_min,                  &
                             chunks(c)%field%x_max,                  &
                             chunks(c)%field%y_min,                  &
                             chunks(c)%field%y_max,                  &
                             dt,                                     &
                             chunks(c)%field%celldx,                 &
                             chunks(c)%field%celldy,                 &
                             chunks(c)%field%volume,                 &
                             chunks(c)%field%density0,               &
                             chunks(c)%field%energy0,                &
                             chunks(c)%field%work_array1,            &
                             chunks(c)%field%work_array2,            &
                             chunks(c)%field%work_array3,            &
                             chunks(c)%field%work_array4,            &
                             chunks(c)%field%work_array5,            &
                             chunks(c)%field%work_array6             )
      ENDIF

    ENDIF

  ENDDO

END SUBROUTINE tea_leaf

END MODULE tea_leaf_module


