# 1 "floats.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "floats.F90"
MODULE floats
   !!======================================================================
   !!                       ***  MODULE  floats  ***
   !! Ocean floats : floats
   !!======================================================================
   !! History :  OPA  !          (CLIPPER)   original Code
   !!   NEMO     1.0  ! 2002-06  (A. Bozec)  F90, Free form and module
   !!----------------------------------------------------------------------
# 139 "floats.F90"
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE flo_stp( kt )          ! Empty routine
      IMPLICIT NONE
      INTEGER, INTENT( in ) :: kt
      WRITE(*,*) 'flo_stp: You should not have seen this print! error?', kt
   END SUBROUTINE flo_stp
   SUBROUTINE flo_init          ! Empty routine
      IMPLICIT NONE
   END SUBROUTINE flo_init


   !!======================================================================
 END MODULE floats
