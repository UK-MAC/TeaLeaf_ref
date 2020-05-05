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

!>  @brief TeaLeaf top level program: Invokes the main cycle
!>  @author David Beckingsale, Wayne Gaudin
!>  @details TeaLeaf in a proxy-app that solves the linear heat conduction
!>  equations using an implicit finite volume method on a Cartesian grid.
!>  The grid is staggered with internal energy, density, and temperature at cell
!>  centres.
!>
!>  It can be run in distributed mode using MPI.
!>
!>  It can use OpenMP, OpenACC on a compute device.
!>
!>  NOTE: that the proxy-app uses uniformly spaced mesh. The actual method will
!>  work on a mesh with varying spacing to keep it relevant to it's parent code.
!>  For this reason, optimisations should only be carried out on the software
!>  that do not change the underlying numerical method. For example, the
!>  volume, though constant for all cells, should remain array and not be
!>  converted to a scalar.
PROGRAM tea_leaf

  USE tea_module

  IMPLICIT NONE

!$ INTEGER :: OMP_GET_NUM_THREADS,OMP_GET_THREAD_NUM

  CALL tea_init_comms()

!$OMP PARALLEL
  IF(parallel%boss)THEN
!$  IF(OMP_GET_THREAD_NUM().EQ.0) THEN
      WRITE(*,*)
      WRITE(*,'(a15,f8.3)') 'Tea Version ',g_version
      WRITE(*,'(a18)') 'MPI Version'
!$    WRITE(*,'(a18)') 'OpenMP Version'
      WRITE(*,'(a14,i6)') 'Task Count ',parallel%max_task !MPI
!$    WRITE(*,'(a15,i5)') 'Thread Count: ',OMP_GET_NUM_THREADS()
      WRITE(*,*)
      WRITE(0,*)
      WRITE(0,'(a15,f8.3)') 'Tea Version ',g_version
      WRITE(0,'(a18)') 'MPI Version'
!$    WRITE(0,'(a18)') 'OpenMP Version'
      WRITE(0,'(a14,i6)') 'Task Count ',parallel%max_task !MPI
!$    WRITE(0,'(a15,i5)') 'Thread Count: ',OMP_GET_NUM_THREADS()
      WRITE(0,*)
!$  ENDIF
  ENDIF
!$OMP END PARALLEL

  CALL initialise

  CALL diffuse

  ! Deallocate everything

END PROGRAM tea_leaf

