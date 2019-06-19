# 1 "sbcflx.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "sbcflx.F90"
MODULE sbcflx
   !!======================================================================
   !!                       ***  MODULE  sbcflx  ***
   !! Ocean forcing:  momentum, heat and freshwater flux formulation
   !!=====================================================================
   !! History :  1.0  !  2006-06  (G. Madec)  Original code
   !!            3.3  !  2010-10  (S. Masson)  add diurnal cycle
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   namflx   : flux formulation namlist
   !!   sbc_flx  : flux formulation as ocean surface boundary condition (forced mode, fluxes read in NetCDF files)
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce         ! surface boundary condition: ocean fields
   USE sbcdcy          ! surface boundary condition: diurnal cycle on qsr
   USE phycst          ! physical constants
   !
   USE fldread         ! read input fields
   USE iom             ! IOM library
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC sbc_flx       ! routine called by step.F90

   INTEGER , PARAMETER ::   jpfld   = 5   ! maximum number of files to read 
   INTEGER , PARAMETER ::   jp_utau = 1   ! index of wind stress (i-component) file
   INTEGER , PARAMETER ::   jp_vtau = 2   ! index of wind stress (j-component) file
   INTEGER , PARAMETER ::   jp_qtot = 3   ! index of total (non solar+solar) heat file
   INTEGER , PARAMETER ::   jp_qsr  = 4   ! index of solar heat file
   INTEGER , PARAMETER ::   jp_emp  = 5   ! index of evaporation-precipation file
   TYPE(FLD), ALLOCATABLE, DIMENSION(:) ::   sf    ! structure of input fields (file informations, fields read)

   !! * Substitutions

# 1 "vectopt_loop_substitute.h90" 1
   !!----------------------------------------------------------------------
   !!                   ***  vectopt_loop_substitute  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute the inner loop start/end indices with CPP macro
   !!                allow unrolling of do-loop (useful with vector processors)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: vectopt_loop_substitute.h90 10068 2018-08-28 14:09:04Z nicolasmartin $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
# 41 "sbcflx.F90" 2
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: sbcflx.F90 10425 2018-12-19 21:54:16Z smasson $
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_flx( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_flx  ***
      !!                   
      !! ** Purpose :   provide at each time step the surface ocean fluxes
      !!                (momentum, heat, freshwater and runoff) 
      !!
      !! ** Method  : - READ each fluxes in NetCDF files:
      !!                   i-component of the stress              utau  (N/m2)
      !!                   j-component of the stress              vtau  (N/m2)
      !!                   net downward heat flux                 qtot  (watt/m2)
      !!                   net downward radiative flux            qsr   (watt/m2)
      !!                   net upward freshwater (evapo - precip) emp   (kg/m2/s)
      !!
      !!      CAUTION :  - never mask the surface stress fields
      !!                 - the stress is assumed to be in the (i,j) mesh referential
      !!
      !! ** Action  :   update at each time-step
      !!              - utau, vtau  i- and j-component of the wind stress
      !!              - taum        wind stress module at T-point
      !!              - wndm        10m wind module at T-point
      !!              - qns         non solar heat flux including heat flux due to emp
      !!              - qsr         solar heat flux
      !!              - emp         upward mass flux (evap. - precip.)
      !!              - sfx         salt flux; set to zero at nit000 but possibly non-zero
      !!                            if ice is present
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!
      INTEGER  ::   ji, jj, jf            ! dummy indices
      INTEGER  ::   ierror                ! return error code
      INTEGER  ::   ios                   ! Local integer output status for namelist read
      REAL(wp) ::   zfact                 ! temporary scalar
      REAL(wp) ::   zrhoa  = 1.22         ! Air density kg/m3
      REAL(wp) ::   zcdrag = 1.5e-3       ! drag coefficient
      REAL(wp) ::   ztx, zty, zmod, zcoef ! temporary variables
      !!
      CHARACTER(len=100) ::  cn_dir                               ! Root directory for location of flx files
      TYPE(FLD_N), DIMENSION(jpfld) ::   slf_i                    ! array of namelist information structures
      TYPE(FLD_N) ::   sn_utau, sn_vtau, sn_qtot, sn_qsr, sn_emp  ! informations about the fields to be read
      NAMELIST/namsbc_flx/ cn_dir, sn_utau, sn_vtau, sn_qtot, sn_qsr, sn_emp
      !!---------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN                ! First call kt=nit000  
         ! set file information
         REWIND( numnam_ref )              ! Namelist namsbc_flx in reference namelist : Files for fluxes
         READ  ( numnam_ref, namsbc_flx, IOSTAT = ios, ERR = 901)
901      IF( ios /= 0 )   CALL ctl_nam ( ios , 'namsbc_flx in reference namelist', lwp )

         REWIND( numnam_cfg )              ! Namelist namsbc_flx in configuration namelist : Files for fluxes
         READ  ( numnam_cfg, namsbc_flx, IOSTAT = ios, ERR = 902 )
902      IF( ios >  0 )   CALL ctl_nam ( ios , 'namsbc_flx in configuration namelist', lwp )
         IF(lwm) WRITE ( numond, namsbc_flx ) 
         !
         !                                         ! check: do we plan to use ln_dm2dc with non-daily forcing?
         IF( ln_dm2dc .AND. sn_qsr%nfreqh /= 24 )   &
            &   CALL ctl_stop( 'sbc_blk_core: ln_dm2dc can be activated only with daily short-wave forcing' ) 
         !
         !                                         ! store namelist information in an array
         slf_i(jp_utau) = sn_utau   ;   slf_i(jp_vtau) = sn_vtau
         slf_i(jp_qtot) = sn_qtot   ;   slf_i(jp_qsr ) = sn_qsr 
         slf_i(jp_emp ) = sn_emp
         !
         ALLOCATE( sf(jpfld), STAT=ierror )        ! set sf structure
         IF( ierror > 0 ) THEN   
            CALL ctl_stop( 'sbc_flx: unable to allocate sf structure' )   ;   RETURN  
         ENDIF
         DO ji= 1, jpfld
            ALLOCATE( sf(ji)%fnow(jpi,jpj,1) )
            IF( slf_i(ji)%ln_tint ) ALLOCATE( sf(ji)%fdta(jpi,jpj,1,2) )
         END DO
         !                                         ! fill sf with slf_i and control print
         CALL fld_fill( sf, slf_i, cn_dir, 'sbc_flx', 'flux formulation for ocean surface boundary condition', 'namsbc_flx' )
         !
         sfx(:,:) = 0.0_wp                         ! salt flux due to freezing/melting (non-zero only if ice is present)
         !
      ENDIF

      CALL fld_read( kt, nn_fsbc, sf )                            ! input fields provided at the current time-step
     
      IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN                        ! update ocean fluxes at each SBC frequency

         IF( ln_dm2dc ) THEN   ;   qsr(:,:) = sbc_dcy( sf(jp_qsr)%fnow(:,:,1) )   ! modify now Qsr to include the diurnal cycle
         ELSE                  ;   qsr(:,:) =          sf(jp_qsr)%fnow(:,:,1)
         ENDIF
         DO jj = 1, jpj                                           ! set the ocean fluxes from read fields
            DO ji = 1, jpi
               utau(ji,jj) = sf(jp_utau)%fnow(ji,jj,1)
               vtau(ji,jj) = sf(jp_vtau)%fnow(ji,jj,1)
               qns (ji,jj) = sf(jp_qtot)%fnow(ji,jj,1) - sf(jp_qsr)%fnow(ji,jj,1)
               emp (ji,jj) = sf(jp_emp )%fnow(ji,jj,1)
            END DO
         END DO
         !                                                        ! add to qns the heat due to e-p
         qns(:,:) = qns(:,:) - emp(:,:) * sst_m(:,:) * rcp        ! mass flux is at SST
         !
         qns(:,:) = qns(:,:) * tmask(:,:,1)
         emp(:,:) = emp(:,:) * tmask(:,:,1)
         !
         !                                                        ! module of wind stress and wind speed at T-point
         zcoef = 1. / ( zrhoa * zcdrag )
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vect. opt.
               ztx = utau(ji-1,jj  ) + utau(ji,jj) 
               zty = vtau(ji  ,jj-1) + vtau(ji,jj) 
               zmod = 0.5 * SQRT( ztx * ztx + zty * zty )
               taum(ji,jj) = zmod
               wndm(ji,jj) = SQRT( zmod * zcoef )
            END DO
         END DO
         taum(:,:) = taum(:,:) * tmask(:,:,1) ; wndm(:,:) = wndm(:,:) * tmask(:,:,1)
         CALL lbc_lnk( 'sbcflx', taum(:,:), 'T', 1. )   ;   CALL lbc_lnk( 'sbcflx', wndm(:,:), 'T', 1. )

         IF( nitend-nit000 <= 100 .AND. lwp ) THEN                ! control print (if less than 100 time-step asked)
            WRITE(numout,*) 
            WRITE(numout,*) '        read daily momentum, heat and freshwater fluxes OK'
            DO jf = 1, jpfld
               IF( jf == jp_utau .OR. jf == jp_vtau )   zfact =     1.
               IF( jf == jp_qtot .OR. jf == jp_qsr  )   zfact =     0.1
               IF( jf == jp_emp                     )   zfact = 86400.
               WRITE(numout,*) 
               WRITE(numout,*) ' day: ', ndastp , TRIM(sf(jf)%clvar), ' * ', zfact
            END DO
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE sbc_flx

   !!======================================================================
END MODULE sbcflx
