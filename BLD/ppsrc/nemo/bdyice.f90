# 1 "bdyice.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "bdyice.F90"
MODULE bdyice
   !!======================================================================
   !!                       ***  MODULE  bdyice  ***
   !! Unstructured Open Boundary Cond. :  Open boundary conditions for sea-ice (SI3)
   !!======================================================================
   !!  History :  3.3  !  2010-09 (D. Storkey)  Original code
   !!             3.4  !  2012-01 (C. Rousset)  add new sea ice model 
   !!             4.0  !  2018    (C. Rousset)  SI3 compatibility 
   !!----------------------------------------------------------------------
# 360 "bdyice.F90"
   !!---------------------------------------------------------------------------------
   !!   Default option                                                    Empty module
   !!---------------------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_ice( kt )      ! Empty routine
      IMPLICIT NONE
      INTEGER, INTENT( in ) :: kt
      WRITE(*,*) 'bdy_ice: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_ice


   !!=================================================================================
END MODULE bdyice
