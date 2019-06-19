# 1 "diadct.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "diadct.F90"
MODULE diadct
   !!======================================================================
   !!                       ***  MODULE  diadct  ***
   !! Ocean diagnostics: Compute the transport trough a sec.
   !!======================================================================
   !! History :  OPA  ! 02/1999 (Y Drillet)  original code
   !!                 ! 10/2001 (Y Drillet, R Bourdalle Badie)
   !!   NEMO     1.0  ! 10/2005 (M Laborie) F90
   !!            3.0  ! 04/2007 (G Garric) Ice sections
   !!             -   ! 04/2007 (C Bricaud) test on sec%nb_point, initialisation of ztransp1,ztransp2,...
   !!            3.4  ! 09/2011 (C Bricaud)
   !!----------------------------------------------------------------------
# 1242 "diadct.F90"
   !!----------------------------------------------------------------------
   !!   Default option :                                       Dummy module
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_diadct = .FALSE.    !: diamht flag
   PUBLIC 
   !! $Id: diadct.F90 10425 2018-12-19 21:54:16Z smasson $
CONTAINS

   SUBROUTINE dia_dct_init          ! Dummy routine
      IMPLICIT NONE
      WRITE(*,*) 'dia_dct_init: You should not have seen this print! error?'
   END SUBROUTINE dia_dct_init

   SUBROUTINE dia_dct( kt )         ! Dummy routine
      IMPLICIT NONE
      INTEGER, INTENT( in ) :: kt   ! ocean time-step index
      WRITE(*,*) 'dia_dct: You should not have seen this print! error?', kt
   END SUBROUTINE dia_dct


   !!======================================================================
END MODULE diadct
