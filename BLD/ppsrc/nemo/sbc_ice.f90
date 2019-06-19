# 1 "sbc_ice.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "sbc_ice.F90"
MODULE sbc_ice
   !!======================================================================
   !!                 ***  MODULE  sbc_ice  ***
   !! Surface module - SI3 & CICE: parameters & variables defined in memory
   !!======================================================================
   !! History :  3.0   !  2006-08  (G. Madec)        Surface module
   !!            3.2   !  2009-06  (S. Masson)       merge with ice_oce
   !!            3.3.1 !  2011-01  (A. R. Porter, STFC Daresbury) dynamical allocation
   !!            3.4   !  2011-11  (C. Harris)       CICE added as an option
   !!            4.0   !  2018     (many people)     SI3 compatibility
   !!----------------------------------------------------------------------
# 157 "sbc_ice.F90"
   !!----------------------------------------------------------------------
   !!   Default option                      NO SI3 or CICE sea-ice model
   !!----------------------------------------------------------------------
   USE lib_mpp          ! MPP library
   USE in_out_manager   ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_ice_alloc   !

   LOGICAL         , PUBLIC, PARAMETER ::   lk_si3     = .FALSE.  !: no SI3 ice model
   LOGICAL         , PUBLIC, PARAMETER ::   lk_cice    = .FALSE.  !: no CICE ice model
   REAL(wp)        , PUBLIC, PARAMETER ::   cldf_ice = 0.81       !: cloud fraction over sea ice, summer CLIO value   [-]
   INTEGER         , PUBLIC, PARAMETER ::   jpl = 1 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   u_ice, v_ice                        ! jpi, jpj
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tn_ice, alb_ice, qns_ice, dqns_ice  ! (jpi,jpj,jpl)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   a_i
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   emp_ice
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qsr_ice
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   h_i, h_s
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   topmelt, botmelt
   !
   !! arrays related to embedding ice in the ocean. 
   !! These arrays need to be declared even if no ice model is required. 
   !! In the no ice model or traditional levitating ice cases they contain only zeros
   !! ---------------------
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   snwice_mass        !: mass of snow and ice at current  ice time step   [Kg/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   snwice_mass_b      !: mass of snow and ice at previous ice time step   [Kg/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   snwice_fmass       !: time evolution of mass of snow+ice               [Kg/m2/s]
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_ice_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  FUNCTION sbc_ice_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: ierr(1)
      !!----------------------------------------------------------------------
      ierr(:) = 0
      ALLOCATE( snwice_mass(jpi,jpj) , snwice_mass_b(jpi,jpj), snwice_fmass(jpi,jpj) , STAT=ierr(1) )
      sbc_ice_alloc = MAXVAL( ierr )
      CALL mpp_sum ( 'sbc_ice', sbc_ice_alloc )
      IF( sbc_ice_alloc > 0 )   CALL ctl_warn('sbc_ice_alloc: allocation of arrays failed')
   END FUNCTION sbc_ice_alloc


   !!======================================================================
END MODULE sbc_ice
