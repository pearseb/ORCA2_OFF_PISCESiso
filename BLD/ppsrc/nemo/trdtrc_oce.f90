# 1 "trdtrc_oce.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "trdtrc_oce.F90"
MODULE trdtrc_oce
   !!======================================================================
   !!                   ***  MODULE trdtrc_oce  ***
   !! Ocean trends :   set tracer and momentum trend variables
   !!======================================================================

   !!----------------------------------------------------------------------
   !!   'key_top'                                                TOP models
   !!----------------------------------------------------------------------
   USE par_oce       ! ocean parameters
   USE par_trc       ! passive tracers parameters

   IMPLICIT NONE
   PUBLIC

   !                                         !!* Namelist namtoptrd:  diagnostics on passive tracers trends
   INTEGER  ::    nn_trd_trc                  !: time step frequency dynamics and tracers trends
   INTEGER  ::    nn_ctls_trc                 !: control surface type for trends vertical integration
   REAL(wp) ::    rn_ucf_trc                  !: unit conversion factor (for netCDF trends outputs)
   LOGICAL  ::    ln_trdmxl_trc_instant       !: flag to diagnose inst./mean ML trc trends
   LOGICAL  ::    ln_trdmxl_trc_restart       !: flag to restart mixed-layer trc diagnostics
   CHARACTER(len=50) ::  cn_trdrst_trc_in     !: suffix of pass. tracer restart name (input)
   CHARACTER(len=50) ::  cn_trdrst_trc_out    !: suffix of pass. tracer restart name (output)
   LOGICAL, DIMENSION(:), ALLOCATABLE ::   ln_trdtrc   !: large trends diagnostic to write or not (namelist)




   LOGICAL, PARAMETER ::   lk_trdtrc = .FALSE.   !: ML trend flag


# 105 "trdtrc_oce.F90"
   LOGICAL, PARAMETER ::   lk_trdmxl_trc = .FALSE.   !: ML trend flag


   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: trdtrc_oce.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION trd_trc_oce_alloc()
      !!----------------------------------------------------------------------
      !!         *** ROUTINE trd_trc_oce_alloc ***
      !!----------------------------------------------------------------------
      USE lib_mpp, ONLY: ctl_stop
      INTEGER :: ierr(2)
      !!----------------------------------------------------------------------
      ierr(:) = 0
      !
# 145 "trdtrc_oce.F90"
      !
      trd_trc_oce_alloc = MAXVAL(ierr)
      !
      IF( trd_trc_oce_alloc /= 0 )   CALL ctl_stop( 'STOP', 'trd_trc_oce_alloc: failed to allocate arrays' )
      !



      !
   END FUNCTION trd_trc_oce_alloc







   !!======================================================================
END MODULE trdtrc_oce
