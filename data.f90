MODULE data_module

   IMPLICIT NONE

   REAL(KIND=8), PARAMETER :: g_version=0.0

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

   INTEGER,         PARAMETER :: FIELD_DENSITY0   = 1         &
                                ,FIELD_DENSITY1   = 2         &
                                ,FIELD_ENERGY0    = 3         &
                                ,FIELD_ENERGY1    = 4         &
                                ,FIELD_PRESSURE   = 5         &
                                ,FIELD_VISCOSITY  = 6         &
                                ,FIELD_SOUNDSPEED = 7         &
                                ,FIELD_XVEL0      = 8         &
                                ,FIELD_XVEL1      = 9         &
                                ,FIELD_YVEL0      =10         &
                                ,FIELD_YVEL1      =11         &
                                ,FIELD_VOL_FLUX_X =12         &
                                ,FIELD_VOL_FLUX_Y =13         &
                                ,FIELD_MASS_FLUX_X=14         &
                                ,FIELD_MASS_FLUX_Y=15         &
                                ,FIELD_U          =16         &
                                ,NUM_FIELDS       =16

   INTEGER,         PARAMETER :: CELL_DATA     = 1,        &
                                 VERTEX_DATA   = 3,        &
                                 X_FACE_DATA   = 3,        &
                                 y_FACE_DATA   = 4


   ! Time step control constants
   INTEGER,        PARAMETER ::  SOUND = 1     &
                                ,X_VEL = 2     &
                                ,Y_VEL = 3     &
                                ,DIVERG= 4

   INTEGER,                      PARAMETER :: g_rect=1 &
                                ,g_circ=2

   INTEGER         ::            g_in           & ! File for input data.
                                ,g_out


   TYPE parallel_type
      LOGICAL           ::      parallel &
                               ,boss
      INTEGER         ::        max_task &
                               ,task     &
                               ,boss_task

   END TYPE parallel_type
   
   TYPE(parallel_type) :: parallel
     
   INTEGER,        PARAMETER ::g_len_max=500

END MODULE data_module
