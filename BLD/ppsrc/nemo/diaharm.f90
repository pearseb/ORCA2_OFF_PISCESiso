# 1 "diaharm.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "diaharm.F90"
MODULE diaharm 
   !!======================================================================
   !!                       ***  MODULE  diaharm  ***
   !! Harmonic analysis of tidal constituents 
   !!======================================================================
   !! History :  3.1  !  2007  (O. Le Galloudec, J. Chanut)  Original code
   !!----------------------------------------------------------------------
# 518 "diaharm.F90"
   !!----------------------------------------------------------------------
   !!   Default case :   Empty module
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_diaharm = .FALSE.
CONTAINS
   SUBROUTINE dia_harm ( kt )     ! Empty routine
      INTEGER, INTENT( IN ) :: kt  
      WRITE(*,*) 'dia_harm: you should not have seen this print'
   END SUBROUTINE dia_harm


   !!======================================================================
END MODULE diaharm
