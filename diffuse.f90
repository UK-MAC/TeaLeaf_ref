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

!>  @brief Controls the main diffusion cycle.
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Controls the top level cycle, invoking all the drivers and checks
!>  for outputs and completion.

SUBROUTINE diffuse

  USE tea_module
  USE timestep_module
  USE tea_leaf_module
  USE set_field_module

  IMPLICIT NONE

  INTEGER         :: loc(1)
  REAL(KIND=8)    :: timer,timerstart,wall_clock,step_clock
  
  REAL(KIND=8)    :: grind_time,cells,rstep
  REAL(KIND=8)    :: step_time,step_grind
  REAL(KIND=8)    :: first_step,second_step
  REAL(KIND=8)    :: kernel_total,totals(parallel%max_task)

  timerstart = timer()

  second_step=0.0 ! In order to prevent unused error

  ! copy time level 0 to time level 1
  CALL set_field()

  DO

    step_time = timer()

    step = step + 1

    CALL timestep()

    CALL tea_leaf()
    
    time = time + dt

    IF(summary_frequency.NE.0) THEN
      IF(MOD(step, summary_frequency).EQ.0) CALL field_summary()
    ENDIF
    IF(visit_frequency.NE.0) THEN
      IF(MOD(step, visit_frequency).EQ.0) CALL visit()
    ENDIF

    ! Sometimes there can be a significant start up cost that appears in the first step.
    ! Sometimes it is due to the number of MPI tasks, or OpenCL kernel compilation.
    ! On the short test runs, this can skew the results, so should be taken into account
    !  in recorded run times.
    IF(step.EQ.1) first_step=(timer() - step_time)
    IF(step.EQ.2) second_step=(timer() - step_time)

    IF (parallel%boss) THEN
      wall_clock=timer()-timerstart
      step_clock=timer()-step_time
      WRITE(g_out,*)"Wall clock ",wall_clock
      WRITE(0    ,*)"Wall clock ",wall_clock
      cells = grid%x_cells * grid%y_cells
      rstep = step
      grind_time   = wall_clock/(rstep * cells)
      step_grind   = step_clock/cells
      WRITE(0    ,*)"Average time per cell ",grind_time
      WRITE(g_out,*)"Average time per cell ",grind_time
      WRITE(0    ,*)"Step time per cell    ",step_grind
      WRITE(g_out,*)"Step time per cell    ",step_grind

    END IF

    IF(time+g_small.GT.end_time.OR.step.GE.end_step) THEN

      complete=.TRUE.
      CALL field_summary()
      IF(visit_frequency.NE.0) CALL visit()

      wall_clock=timer() - timerstart
      IF ( parallel%boss ) THEN
        WRITE(g_out,*)
        WRITE(g_out,*) 'Calculation complete'
        WRITE(g_out,*) 'Tea is finishing'
        WRITE(g_out,*) 'Wall clock ', wall_clock
        WRITE(g_out,*) 'First step overhead', first_step-second_step
        WRITE(    0,*) 'Wall clock ', wall_clock
        WRITE(    0,*) 'First step overhead', first_step-second_step
      ENDIF

      EXIT

    ENDIF
  ENDDO

  IF ( profiler_on ) THEN
    ! First we need to find the maximum kernel time for each task. This
    ! seems to work better than finding the maximum time for each kernel and
    ! adding it up, which always gives over 100%. I think this is because it
    ! does not take into account compute overlaps before syncronisations
    ! caused by halo exhanges.
    kernel_total=profiler%timestep+profiler%halo_exchange+profiler%summary+profiler%visit+&
        profiler%tea_init+profiler%set_field+profiler%tea_solve+profiler%tea_reset
    CALL tea_allgather(kernel_total,totals)
    ! So then what I do is use the individual kernel times for the
    ! maximum kernel time task for the profile print
    loc=MAXLOC(totals)
    kernel_total=totals(loc(1))
    CALL tea_allgather(profiler%timestep,totals)
    profiler%timestep=totals(loc(1))
    CALL tea_allgather(profiler%halo_exchange,totals)
    profiler%halo_exchange=totals(loc(1))
    CALL tea_allgather(profiler%summary,totals)
    profiler%summary=totals(loc(1))
    CALL tea_allgather(profiler%visit,totals)
    profiler%visit=totals(loc(1))
    CALL tea_allgather(profiler%tea_init,totals)
    profiler%tea_init=totals(loc(1))
    CALL tea_allgather(profiler%tea_solve,totals)
    profiler%tea_solve=totals(loc(1))
    CALL tea_allgather(profiler%tea_reset,totals)
    profiler%tea_reset=totals(loc(1))
    CALL tea_allgather(profiler%set_field,totals)
    profiler%set_field=totals(loc(1))

    IF ( parallel%boss ) THEN
      WRITE(g_out,*)
      WRITE(g_out,'(a58,2f16.4)')"Profiler Output                 Time            Percentage"
      WRITE(g_out,'(a23,2f16.4)')"Timestep              :",profiler%timestep,100.0*(profiler%timestep/wall_clock)
      WRITE(g_out,'(a23,2f16.4)')"Halo Exchange         :",profiler%halo_exchange,100.0*(profiler%halo_exchange/wall_clock)
      WRITE(g_out,'(a23,2f16.4)')"Summary               :",profiler%summary,100.0*(profiler%summary/wall_clock)
      WRITE(g_out,'(a23,2f16.4)')"Visit                 :",profiler%visit,100.0*(profiler%visit/wall_clock)
      WRITE(g_out,'(a23,2f16.4)')"Tea Init              :",profiler%tea_init,100.0*(profiler%tea_init/wall_clock)
      WRITE(g_out,'(a23,2f16.4)')"Tea Solve             :",profiler%tea_solve,100.0*(profiler%tea_solve/wall_clock)
      WRITE(g_out,'(a23,2f16.4)')"Tea Reset             :",profiler%tea_reset,100.0*(profiler%tea_reset/wall_clock)
      WRITE(g_out,'(a23,2f16.4)')"Set Field             :",profiler%set_field,100.0*(profiler%set_field/wall_clock)
      WRITE(g_out,'(a23,2f16.4)')"Total                 :",kernel_total,100.0*(kernel_total/wall_clock)
      WRITE(g_out,'(a23,2f16.4)')"The Rest              :",wall_clock-kernel_total,100.0*(wall_clock-kernel_total)/wall_clock
    ENDIF

  ENDIF

  CALL tea_finalize

END SUBROUTINE diffuse
