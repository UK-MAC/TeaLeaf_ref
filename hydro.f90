SUBROUTINE hydro

  USE clover_module
  USE timestep_module
  USE viscosity_module
  USE PdV_module
  USE accelerate_module
  USE flux_calc_module
  USE advection_module
  USE tea_leaf_module
  USE reset_field_module
  USE set_field_module

  IMPLICIT NONE

  INTEGER         :: cells
  REAL(KIND=8)    :: timer,timerstart,timerend
  
  REAL(KIND=8)    :: grind_time
  REAL(KIND=8)    :: step_time,step_grind

  timerstart = timer()

  DO

    step_time = timer()

    step = step + 1

    CALL timestep()

    IF (use_Hydro) THEN
    CALL PdV(.TRUE.)

    CALL accelerate()

    CALL PdV(.FALSE.)

    CALL flux_calc()

    CALL advection()
    ENDIF

    IF(use_TeaLeaf) THEN
        IF(.NOT. use_Hydro) THEN
        ! copy tl0 to tl1
        CALL set_field()
        ENDIF
        
        CALL tea_leaf()
    ENDIF
    
    CALL reset_field()

    advect_x = .NOT. advect_x
  
    time = time + dt

    IF(summary_frequency.NE.0) THEN
      IF(MOD(step, summary_frequency).EQ.0) CALL field_summary()
    ENDIF
    IF(visit_frequency.NE.0) THEN
      IF(MOD(step, visit_frequency).EQ.0) CALL visit()
    ENDIF

    IF(time+g_small.GT.end_time.OR.step.GE.end_step) THEN

      CALL field_summary()
      IF(visit_frequency.NE.0) CALL visit()

      IF ( parallel%boss ) THEN
        WRITE(g_out,*)
        WRITE(g_out,*) 'Calculation complete'
        WRITE(g_out,*) 'Clover is finishing'
        WRITE(g_out,*) 'Wall clock ', timer() - timerstart
        WRITE(    0,*) 'Wall clock ', timer() - timerstart
      ENDIF

      CALL clover_finalize

      EXIT

    END IF

    IF (parallel%boss) THEN
      WRITE(g_out,*)"Wall clock ",timer()-timerstart
      WRITE(0    ,*)"Wall clock ",timer()-timerstart
      cells = grid%x_cells * grid%y_cells
      grind_time   = (timer() - timerstart) / (step * cells)
      step_grind   = (timer() - step_time)/cells
      WRITE(0    ,*)"Average time per cell ",grind_time
      WRITE(g_out,*)"Average time per cell ",grind_time
      WRITE(0    ,*)"Step time per cell    ",step_grind
      WRITE(g_out,*)"Step time per cell    ",step_grind

     END IF

  END DO

END SUBROUTINE hydro
