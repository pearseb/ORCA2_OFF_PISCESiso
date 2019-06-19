# 1 "diahth.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "diahth.F90"
MODULE diahth
   !!======================================================================
   !!                       ***  MODULE  diahth  ***
   !! Ocean diagnostics: thermocline and 20 degree depth
   !!======================================================================
   !! History :  OPA  !  1994-09  (J.-P. Boulanger)  Original code
   !!                 !  1996-11  (E. Guilyardi)  OPA8 
   !!                 !  1997-08  (G. Madec)  optimization
   !!                 !  1999-07  (E. Guilyardi)  hd28 + heat content 
   !!   NEMO     1.0  !  2002-06  (G. Madec)  F90: Free form and module
   !!            3.2  !  2009-07  (S. Masson) hc300 bugfix + cleaning + add new diag
   !!----------------------------------------------------------------------
# 339 "diahth.F90"
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------
   LOGICAL , PUBLIC, PARAMETER ::   lk_diahth = .FALSE.  !: thermocline-20d depths flag
CONTAINS
   SUBROUTINE dia_hth( kt )         ! Empty routine
      IMPLICIT NONE
      INTEGER, INTENT( in ) :: kt
      WRITE(*,*) 'dia_hth: You should not have seen this print! error?', kt
   END SUBROUTINE dia_hth


   !!======================================================================
END MODULE diahth
