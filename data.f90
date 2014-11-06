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

!>  @brief Holds parameters definitions
!>  @author David Beckingsale, Wayne Gaudin
!>  @details Parameters used in the TeaLeaf are defined here.

MODULE data_module

   IMPLICIT NONE

   REAL(KIND=8), PARAMETER :: g_version=1.0

   INTEGER,      PARAMETER :: g_ibig=640000

   REAL(KIND=8), PARAMETER :: g_small=1.0e-16  &
                             ,g_big  =1.0e+21

   INTEGER,      PARAMETER :: g_name_len_max=255 &
                             ,g_xdir=1           &
                             ,g_ydir=2

   ! These two need to be kept consistent with update_halo
   INTEGER,      PARAMETER :: CHUNK_LEFT   =1    &
                             ,CHUNK_RIGHT  =2    &
                             ,CHUNK_BOTTOM =3    &
                             ,CHUNK_TOP    =4    &
                             ,EXTERNAL_FACE=-1

   INTEGER,         PARAMETER :: FIELD_DENSITY    = 1         &
                                ,FIELD_ENERGY0    = 2         &
                                ,FIELD_ENERGY1    = 3         &
                                ,FIELD_U          = 4         &
                                ,FIELD_P          = 5         &
                                ,FIELD_SD         = 6         &
                                ,NUM_FIELDS       = 6

   INTEGER,         PARAMETER :: CELL_DATA     = 1,        &
                                 VERTEX_DATA   = 2,        &
                                 X_FACE_DATA   = 3,        &
                                 y_FACE_DATA   = 4


   ! Time step control constants
   INTEGER,        PARAMETER ::  FIXED = 1

   INTEGER,                      PARAMETER :: g_rect=1 &
                                             ,g_circ=2 &
                                             ,g_point=3

   INTEGER         ::            g_in           & ! File for input data.
                                ,g_out

   INTEGER         ::            CONDUCTIVITY        = 1 &
                                ,RECIP_CONDUCTIVITY  = 2

   TYPE parallel_type
      LOGICAL           ::      parallel &
                               ,boss
      INTEGER         ::        max_task &
                               ,task     &
                               ,boss_task

   END TYPE parallel_type

   TYPE(parallel_type) :: parallel

   INTEGER,        PARAMETER ::g_len_max=500
   INTEGER,        PARAMETER ::chunks_per_task=1

   INTEGER                   ::lr_pack_buffer_size,bt_pack_buffer_size

END MODULE data_module
