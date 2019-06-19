# 1 "dyncor_c1d.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dyncor_c1d.F90"
MODULE dyncor_c1d
   !!======================================================================
   !!                     ***  MODULE  dyncor_c1d  ***
   !! Ocean Dynamics :   Coriolis term in 1D configuration
   !!=====================================================================
   !! History :  2.0  !  2004-09  (C. Ethe)  Original code
   !!            3.0  !  2008-04  (G. Madec)  style only
   !!----------------------------------------------------------------------
# 103 "dyncor_c1d.F90"
   !!----------------------------------------------------------------------
   !!   Default key                                     NO 1D Configuration
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE cor_c1d              ! Empty routine
       IMPLICIT NONE
   END SUBROUTINE cor_c1d   
   SUBROUTINE dyn_cor_c1d ( kt )      ! Empty routine
      IMPLICIT NONE
      INTEGER, INTENT( in ) :: kt
      WRITE(*,*) 'dyn_cor_c1d: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_cor_c1d


   !!=====================================================================
END MODULE dyncor_c1d
