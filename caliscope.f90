!Crown Copyright 2021 AWE.
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

!>  @brief Caliper routine
!>  @author Wei Liu
!>  @details Integrate external LLNL performance profiling library  

MODULE caliscope_module

   IMPLICIT NONE

   TYPE scope_type
       CHARACTER(LEN=:), ALLOCATABLE, PRIVATE :: name_
       CONTAINS
           PROCEDURE, PUBLIC :: create=>scope_create
           FINAL :: scope_finalise
   END TYPE scope_type

CONTAINS

    SUBROUTINE scope_create(this,name)

        USE caliper_mod, ONLY: cali_begin_region

        CLASS(scope_type), INTENT(INOUT) :: this
        CHARACTER(LEN=*), INTENT(IN) :: name

        IF(.not.ALLOCATED(this%name_)) THEN
            this%name_=name
            CALL cali_begin_region(name)
        ENDIF

    END SUBROUTINE scope_create

    SUBROUTINE scope_finalise(this)

        USE caliper_mod, ONLY: cali_end_region

        TYPE(scope_type), INTENT(INOUT) :: this

        IF(ALLOCATED(this%name_)) THEN
            CALL cali_end_region(this%name_)
            DEALLOCATE(this%name_)
        ENDIF

    END SUBROUTINE scope_finalise

END MODULE caliscope_module

