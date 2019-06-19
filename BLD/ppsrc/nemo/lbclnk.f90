# 1 "lbclnk.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "lbclnk.F90"
MODULE lbclnk
   !!======================================================================
   !!                       ***  MODULE  lbclnk  ***
   !! NEMO        : lateral boundary conditions
   !!=====================================================================
   !! History :  OPA  ! 1997-06  (G. Madec)  Original code
   !!   NEMO     1.0  ! 2002-09  (G. Madec)  F90: Free form and module
   !!            3.2  ! 2009-03  (R. Benshila)  External north fold treatment  
   !!            3.5  ! 2012     (S.Mocavero, I. Epicoco)  optimization of BDY comm. via lbc_bdy_lnk and lbc_obc_lnk
   !!            3.4  ! 2012-12  (R. Bourdalle-Badie, G. Reffray)  add a C1D case  
   !!            3.6  ! 2015-06  (O. Tintó and M. Castrillo)  add lbc_lnk_multi  
   !!            4.0  ! 2017-03  (G. Madec) automatique allocation of array size (use with any 3rd dim size)
   !!             -   ! 2017-04  (G. Madec) remove duplicated routines (lbc_lnk_2d_9, lbc_lnk_2d_multiple, lbc_lnk_3d_gather)
   !!             -   ! 2017-05  (G. Madec) create generic.h90 files to generate all lbc and north fold routines
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_mpp_mpi'             MPI massively parallel processing library
   !!----------------------------------------------------------------------
   !!           define the generic interfaces of lib_mpp routines
   !!----------------------------------------------------------------------
   !!   lbc_lnk       : generic interface for mpp_lnk_3d and mpp_lnk_2d routines defined in lib_mpp
   !!   lbc_bdy_lnk   : generic interface for mpp_lnk_bdy_2d and mpp_lnk_bdy_3d routines defined in lib_mpp
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean dynamics and tracers   
   USE lib_mpp        ! distributed memory computing library
   USE lbcnfd         ! north fold

   INTERFACE lbc_lnk
      MODULE PROCEDURE   mpp_lnk_2d      , mpp_lnk_3d      , mpp_lnk_4d
   END INTERFACE
   INTERFACE lbc_lnk_ptr
      MODULE PROCEDURE   mpp_lnk_2d_ptr  , mpp_lnk_3d_ptr  , mpp_lnk_4d_ptr
   END INTERFACE
   INTERFACE lbc_lnk_multi
      MODULE PROCEDURE   lbc_lnk_2d_multi, lbc_lnk_3d_multi, lbc_lnk_4d_multi
   END INTERFACE
   !
   INTERFACE lbc_bdy_lnk
      MODULE PROCEDURE mpp_lnk_bdy_2d, mpp_lnk_bdy_3d, mpp_lnk_bdy_4d
   END INTERFACE
   !
   INTERFACE lbc_lnk_icb
      MODULE PROCEDURE mpp_lnk_2d_icb
   END INTERFACE

   PUBLIC   lbc_lnk       ! ocean/ice lateral boundary conditions
   PUBLIC   lbc_lnk_multi ! modified ocean/ice lateral boundary conditions
   PUBLIC   lbc_bdy_lnk   ! ocean lateral BDY boundary conditions
   PUBLIC   lbc_lnk_icb   ! iceberg lateral boundary conditions

   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: lbclnk.F90 10425 2018-12-19 21:54:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

# 231 "lbclnk.F90"

   !!======================================================================
   !!   identical routines in both distributed and shared memory computing
   !!======================================================================

   !!----------------------------------------------------------------------
   !!                   ***   load_ptr_(2,3,4)d   ***
   !!
   !!   * Dummy Argument :
   !!       in    ==>   ptab       ! array to be loaded (2D, 3D or 4D)
   !!                   cd_nat     ! nature of pt2d array grid-points
   !!                   psgn       ! sign used across the north fold boundary
   !!       inout <=>   ptab_ptr   ! array of 2D, 3D or 4D pointers
   !!                   cdna_ptr   ! nature of ptab array grid-points
   !!                   psgn_ptr   ! sign used across the north fold boundary
   !!                   kfld       ! number of elements that has been attributed
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!                  ***   lbc_lnk_(2,3,4)d_multi   ***
   !!                     ***   load_ptr_(2,3,4)d   ***
   !!
   !!   * Argument : dummy argument use in lbc_lnk_multi_... routines
   !!
   !!----------------------------------------------------------------------





# 1 "lbc_lnk_multi_generic.h90" 1
# 16 "lbc_lnk_multi_generic.h90"
   SUBROUTINE lbc_lnk_2d_multi( cdname                                                    &
      &                    , pt1, cdna1, psgn1, pt2, cdna2, psgn2, pt3, cdna3, psgn3   &
      &                    , pt4, cdna4, psgn4, pt5, cdna5, psgn5, pt6, cdna6, psgn6   &
      &                    , pt7, cdna7, psgn7, pt8, cdna8, psgn8, pt9, cdna9, psgn9, cd_mpp, pval)
      !!---------------------------------------------------------------------
      CHARACTER(len=*)   ,                   INTENT(in   ) ::   cdname  ! name of the calling subroutine
      REAL(wp), DIMENSION(:,:)          , TARGET, INTENT(inout) ::   pt1     ! arrays on which the lbc is applied
      REAL(wp), DIMENSION(:,:), OPTIONAL, TARGET, INTENT(inout) ::   pt2  ,  pt3  , pt4  , pt5  , pt6  , pt7  , pt8  , pt9
      CHARACTER(len=1)                     , INTENT(in   ) ::   cdna1   ! nature of pt2D. array grid-points
      CHARACTER(len=1)   , OPTIONAL        , INTENT(in   ) ::   cdna2,  cdna3, cdna4, cdna5, cdna6, cdna7, cdna8, cdna9
      REAL(wp)                             , INTENT(in   ) ::   psgn1   ! sign used across the north fold
      REAL(wp)           , OPTIONAL        , INTENT(in   ) ::   psgn2,  psgn3, psgn4, psgn5, psgn6, psgn7, psgn8, psgn9   
      CHARACTER(len=3)   , OPTIONAL        , INTENT(in   ) ::   cd_mpp  ! fill the overlap area only
      REAL(wp)           , OPTIONAL        , INTENT(in   ) ::   pval    ! background value (used at closed boundaries)
      !!
      INTEGER                         ::   kfld        ! number of elements that will be attributed
      TYPE(PTR_2D)         , DIMENSION(9) ::   ptab_ptr    ! pointer array
      CHARACTER(len=1) , DIMENSION(9) ::   cdna_ptr    ! nature of ptab_ptr grid-points
      REAL(wp)         , DIMENSION(9) ::   psgn_ptr    ! sign used across the north fold boundary
      !!---------------------------------------------------------------------
      !
      kfld = 0          ! initial array of pointer size
      !
      !                 ! Load the first array
      CALL load_ptr_2d( pt1, cdna1, psgn1, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !
      !                 ! Look if more arrays are added
      IF( PRESENT(psgn2) )   CALL load_ptr_2d( pt2, cdna2, psgn2, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn3) )   CALL load_ptr_2d( pt3, cdna3, psgn3, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn4) )   CALL load_ptr_2d( pt4, cdna4, psgn4, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn5) )   CALL load_ptr_2d( pt5, cdna5, psgn5, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn6) )   CALL load_ptr_2d( pt6, cdna6, psgn6, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn7) )   CALL load_ptr_2d( pt7, cdna7, psgn7, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn8) )   CALL load_ptr_2d( pt8, cdna8, psgn8, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn9) )   CALL load_ptr_2d( pt9, cdna9, psgn9, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !
      CALL lbc_lnk_ptr( cdname, ptab_ptr, cdna_ptr, psgn_ptr, kfld, cd_mpp, pval )
      !
   END SUBROUTINE lbc_lnk_2d_multi


   SUBROUTINE load_ptr_2d( ptab, cdna, psgn, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:)   , TARGET, INTENT(inout) ::   ptab       ! arrays on which the lbc is applied
      CHARACTER(len=1)              , INTENT(in   ) ::   cdna       ! nature of pt2d array grid-points
      REAL(wp)                      , INTENT(in   ) ::   psgn       ! sign used across the north fold boundary
      TYPE(PTR_2D)        , DIMENSION(:), INTENT(inout) ::   ptab_ptr   ! array of pointers
      CHARACTER(len=1), DIMENSION(:), INTENT(inout) ::   cdna_ptr   ! nature of pt2d_array array grid-points
      REAL(wp)        , DIMENSION(:), INTENT(inout) ::   psgn_ptr   ! sign used across the north fold boundary
      INTEGER                       , INTENT(inout) ::   kfld       ! number of elements that has been attributed
      !!---------------------------------------------------------------------
      !
      kfld                    =  kfld + 1
      ptab_ptr(kfld)%pt2d => ptab
      cdna_ptr(kfld)          =  cdna
      psgn_ptr(kfld)          =  psgn
      !
   END SUBROUTINE load_ptr_2d
# 261 "lbclnk.F90" 2









# 1 "lbc_lnk_multi_generic.h90" 1
# 16 "lbc_lnk_multi_generic.h90"
   SUBROUTINE lbc_lnk_3d_multi( cdname                                                    &
      &                    , pt1, cdna1, psgn1, pt2, cdna2, psgn2, pt3, cdna3, psgn3   &
      &                    , pt4, cdna4, psgn4, pt5, cdna5, psgn5, pt6, cdna6, psgn6   &
      &                    , pt7, cdna7, psgn7, pt8, cdna8, psgn8, pt9, cdna9, psgn9, cd_mpp, pval)
      !!---------------------------------------------------------------------
      CHARACTER(len=*)   ,                   INTENT(in   ) ::   cdname  ! name of the calling subroutine
      REAL(wp), DIMENSION(:,:,:)          , TARGET, INTENT(inout) ::   pt1     ! arrays on which the lbc is applied
      REAL(wp), DIMENSION(:,:,:), OPTIONAL, TARGET, INTENT(inout) ::   pt2  ,  pt3  , pt4  , pt5  , pt6  , pt7  , pt8  , pt9
      CHARACTER(len=1)                     , INTENT(in   ) ::   cdna1   ! nature of pt2D. array grid-points
      CHARACTER(len=1)   , OPTIONAL        , INTENT(in   ) ::   cdna2,  cdna3, cdna4, cdna5, cdna6, cdna7, cdna8, cdna9
      REAL(wp)                             , INTENT(in   ) ::   psgn1   ! sign used across the north fold
      REAL(wp)           , OPTIONAL        , INTENT(in   ) ::   psgn2,  psgn3, psgn4, psgn5, psgn6, psgn7, psgn8, psgn9   
      CHARACTER(len=3)   , OPTIONAL        , INTENT(in   ) ::   cd_mpp  ! fill the overlap area only
      REAL(wp)           , OPTIONAL        , INTENT(in   ) ::   pval    ! background value (used at closed boundaries)
      !!
      INTEGER                         ::   kfld        ! number of elements that will be attributed
      TYPE(PTR_3D)         , DIMENSION(9) ::   ptab_ptr    ! pointer array
      CHARACTER(len=1) , DIMENSION(9) ::   cdna_ptr    ! nature of ptab_ptr grid-points
      REAL(wp)         , DIMENSION(9) ::   psgn_ptr    ! sign used across the north fold boundary
      !!---------------------------------------------------------------------
      !
      kfld = 0          ! initial array of pointer size
      !
      !                 ! Load the first array
      CALL load_ptr_3d( pt1, cdna1, psgn1, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !
      !                 ! Look if more arrays are added
      IF( PRESENT(psgn2) )   CALL load_ptr_3d( pt2, cdna2, psgn2, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn3) )   CALL load_ptr_3d( pt3, cdna3, psgn3, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn4) )   CALL load_ptr_3d( pt4, cdna4, psgn4, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn5) )   CALL load_ptr_3d( pt5, cdna5, psgn5, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn6) )   CALL load_ptr_3d( pt6, cdna6, psgn6, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn7) )   CALL load_ptr_3d( pt7, cdna7, psgn7, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn8) )   CALL load_ptr_3d( pt8, cdna8, psgn8, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn9) )   CALL load_ptr_3d( pt9, cdna9, psgn9, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !
      CALL lbc_lnk_ptr( cdname, ptab_ptr, cdna_ptr, psgn_ptr, kfld, cd_mpp, pval )
      !
   END SUBROUTINE lbc_lnk_3d_multi


   SUBROUTINE load_ptr_3d( ptab, cdna, psgn, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:)   , TARGET, INTENT(inout) ::   ptab       ! arrays on which the lbc is applied
      CHARACTER(len=1)              , INTENT(in   ) ::   cdna       ! nature of pt2d array grid-points
      REAL(wp)                      , INTENT(in   ) ::   psgn       ! sign used across the north fold boundary
      TYPE(PTR_3D)        , DIMENSION(:), INTENT(inout) ::   ptab_ptr   ! array of pointers
      CHARACTER(len=1), DIMENSION(:), INTENT(inout) ::   cdna_ptr   ! nature of pt2d_array array grid-points
      REAL(wp)        , DIMENSION(:), INTENT(inout) ::   psgn_ptr   ! sign used across the north fold boundary
      INTEGER                       , INTENT(inout) ::   kfld       ! number of elements that has been attributed
      !!---------------------------------------------------------------------
      !
      kfld                    =  kfld + 1
      ptab_ptr(kfld)%pt3d => ptab
      cdna_ptr(kfld)          =  cdna
      psgn_ptr(kfld)          =  psgn
      !
   END SUBROUTINE load_ptr_3d
# 270 "lbclnk.F90" 2









# 1 "lbc_lnk_multi_generic.h90" 1
# 16 "lbc_lnk_multi_generic.h90"
   SUBROUTINE lbc_lnk_4d_multi( cdname                                                    &
      &                    , pt1, cdna1, psgn1, pt2, cdna2, psgn2, pt3, cdna3, psgn3   &
      &                    , pt4, cdna4, psgn4, pt5, cdna5, psgn5, pt6, cdna6, psgn6   &
      &                    , pt7, cdna7, psgn7, pt8, cdna8, psgn8, pt9, cdna9, psgn9, cd_mpp, pval)
      !!---------------------------------------------------------------------
      CHARACTER(len=*)   ,                   INTENT(in   ) ::   cdname  ! name of the calling subroutine
      REAL(wp), DIMENSION(:,:,:,:)          , TARGET, INTENT(inout) ::   pt1     ! arrays on which the lbc is applied
      REAL(wp), DIMENSION(:,:,:,:), OPTIONAL, TARGET, INTENT(inout) ::   pt2  ,  pt3  , pt4  , pt5  , pt6  , pt7  , pt8  , pt9
      CHARACTER(len=1)                     , INTENT(in   ) ::   cdna1   ! nature of pt2D. array grid-points
      CHARACTER(len=1)   , OPTIONAL        , INTENT(in   ) ::   cdna2,  cdna3, cdna4, cdna5, cdna6, cdna7, cdna8, cdna9
      REAL(wp)                             , INTENT(in   ) ::   psgn1   ! sign used across the north fold
      REAL(wp)           , OPTIONAL        , INTENT(in   ) ::   psgn2,  psgn3, psgn4, psgn5, psgn6, psgn7, psgn8, psgn9   
      CHARACTER(len=3)   , OPTIONAL        , INTENT(in   ) ::   cd_mpp  ! fill the overlap area only
      REAL(wp)           , OPTIONAL        , INTENT(in   ) ::   pval    ! background value (used at closed boundaries)
      !!
      INTEGER                         ::   kfld        ! number of elements that will be attributed
      TYPE(PTR_4D)         , DIMENSION(9) ::   ptab_ptr    ! pointer array
      CHARACTER(len=1) , DIMENSION(9) ::   cdna_ptr    ! nature of ptab_ptr grid-points
      REAL(wp)         , DIMENSION(9) ::   psgn_ptr    ! sign used across the north fold boundary
      !!---------------------------------------------------------------------
      !
      kfld = 0          ! initial array of pointer size
      !
      !                 ! Load the first array
      CALL load_ptr_4d( pt1, cdna1, psgn1, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !
      !                 ! Look if more arrays are added
      IF( PRESENT(psgn2) )   CALL load_ptr_4d( pt2, cdna2, psgn2, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn3) )   CALL load_ptr_4d( pt3, cdna3, psgn3, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn4) )   CALL load_ptr_4d( pt4, cdna4, psgn4, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn5) )   CALL load_ptr_4d( pt5, cdna5, psgn5, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn6) )   CALL load_ptr_4d( pt6, cdna6, psgn6, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn7) )   CALL load_ptr_4d( pt7, cdna7, psgn7, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn8) )   CALL load_ptr_4d( pt8, cdna8, psgn8, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      IF( PRESENT(psgn9) )   CALL load_ptr_4d( pt9, cdna9, psgn9, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !
      CALL lbc_lnk_ptr( cdname, ptab_ptr, cdna_ptr, psgn_ptr, kfld, cd_mpp, pval )
      !
   END SUBROUTINE lbc_lnk_4d_multi


   SUBROUTINE load_ptr_4d( ptab, cdna, psgn, ptab_ptr, cdna_ptr, psgn_ptr, kfld )
      !!---------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:,:)   , TARGET, INTENT(inout) ::   ptab       ! arrays on which the lbc is applied
      CHARACTER(len=1)              , INTENT(in   ) ::   cdna       ! nature of pt2d array grid-points
      REAL(wp)                      , INTENT(in   ) ::   psgn       ! sign used across the north fold boundary
      TYPE(PTR_4D)        , DIMENSION(:), INTENT(inout) ::   ptab_ptr   ! array of pointers
      CHARACTER(len=1), DIMENSION(:), INTENT(inout) ::   cdna_ptr   ! nature of pt2d_array array grid-points
      REAL(wp)        , DIMENSION(:), INTENT(inout) ::   psgn_ptr   ! sign used across the north fold boundary
      INTEGER                       , INTENT(inout) ::   kfld       ! number of elements that has been attributed
      !!---------------------------------------------------------------------
      !
      kfld                    =  kfld + 1
      ptab_ptr(kfld)%pt4d => ptab
      cdna_ptr(kfld)          =  cdna
      psgn_ptr(kfld)          =  psgn
      !
   END SUBROUTINE load_ptr_4d
# 279 "lbclnk.F90" 2




   !!======================================================================
END MODULE lbclnk

