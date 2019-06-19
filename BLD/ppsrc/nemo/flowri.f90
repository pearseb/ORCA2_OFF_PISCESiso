# 1 "flowri.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "flowri.F90"
MODULE flowri
   !!======================================================================
   !!                       ***  MODULE  flowri  ***
   !!
   !! Ocean floats: write floats trajectory in ascii                    ln_flo_ascii = T
   !!                                    or in netcdf ( IOM or IOSPSL ) ln_flo_ascii = F           
   !!======================================================================
   !!  History :  OPA  !  1999-09  (Y. Drillet)    : Original code
   !!              -   !  2000-06  (J.-M. Molines) : Profiling floats for CLS 
   !!   NEMO      1.0  !  2002-10  (A. Bozec)  F90 : Free form and module
   !!             3.2  !  2010-08  (slaw, cbricaud): netcdf outputs and others 
   !!----------------------------------------------------------------------
# 281 "flowri.F90"
   !!----------------------------------------------------------------------
   !!   Default option                                         Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE flo_wri                 ! Empty routine
   END SUBROUTINE flo_wri


   !!=======================================================================
END MODULE flowri
