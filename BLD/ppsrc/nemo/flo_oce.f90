# 1 "flo_oce.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "flo_oce.F90"
MODULE flo_oce
   !!======================================================================
   !!                     ***  MODULE flo_oce  ***
   !! lagrangian floats :   define in memory all floats parameters and variables
   !!======================================================================
   !! History :   OPA  ! 1999-10  (CLIPPER projet)
   !!   NEMO      1.0  ! 2002-11  (G. Madec, A. Bozec)  F90: Free form and module
   !!----------------------------------------------------------------------
# 71 "flo_oce.F90"
   !!----------------------------------------------------------------------
   !!   Default option :                                 NO drifting floats
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_floats = .FALSE.   !: float flag


   !!======================================================================
END MODULE flo_oce
