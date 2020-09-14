MODULE p4zsed
   !!======================================================================
   !!                         ***  MODULE p4sed  ***
   !! TOP :   PISCES Compute loss of organic matter in the sediments
   !!======================================================================
   !! History :   1.0  !  2004-03 (O. Aumont) Original code
   !!             2.0  !  2007-12 (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06 (C. Ethe) USE of fldread
   !!             3.5  !  2012-07 (O. Aumont) improvment of river input of nutrients 
   !!----------------------------------------------------------------------
   !!   p4z_sed        :  Compute loss of organic matter in the sediments
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zlim          !  Co-limitations of differents nutrients
   USE p4zsbc          !  External source of nutrients 
   USE p4zint          !  interpolation and computation of various fields
   USE sed             !  Sediment module
   USE iom             !  I/O manager
   USE prtctl_trc      !  print control for debugging

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_sed  
   PUBLIC   p4z_sed_alloc
 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: nitrpot    !: Nitrogen fixation 
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:  ) :: sdenit     !: Nitrate reduction in the sediments
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: zpdenit     !: Nitrate reduction in the sediments
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: zpdenit15   !: Nitrate reduction in the sediments
   REAL(wp) :: r1_rday                  !: inverse of rday
   LOGICAL, SAVE :: lk_sed

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zsed.F90 10780 2019-03-20 17:53:44Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_sed( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed  ***
      !!
      !! ** Purpose :   Compute loss of organic matter in the sediments. This
      !!              is by no way a sediment model. The loss is simply 
      !!              computed to balance the inout from rivers and dust
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      !
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      INTEGER  ::  ji, jj, jk, ikt
      REAL(wp) ::  zrivalk, zrivsil, zrivno3
      REAL(wp) ::  zwflux, zlim, zfact, zfactcal
      REAL(wp) ::  zo2, zno3, zflx, z1pdenit 
      REAL(wp) ::  zsiloss, zcaloss, zws3, zws4, zwsc, zdep
      REAL(wp) ::  zwstpoc, zwstpon, zwstpop
      REAL(wp) ::  zwstpoc15, zr15_no3, zr15_rain
      REAL(wp) ::  zwstpoc13, zr13_rain, zr13_cal
      REAL(wp) ::  ztrfer, ztrpo4s, ztrdp, zwdust, zmudia, ztemp
      REAL(wp) ::  xdiano3, xdianh4
      !
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(jpi,jpj    ) :: zdenit2d, zbureff, zwork
      REAL(wp), DIMENSION(jpi,jpj    ) :: zwsbio3, zwsbio4
      REAL(wp), DIMENSION(jpi,jpj    ) :: zsedcal, zsedsi, zsedc
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zsoufer, zlight
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zolimit, zolimit15
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zseddiss, zseddiss13
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zsedremin, zsedremin13
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: ztrpo4, ztrdop, zirondep, zpdep
      REAL(wp), ALLOCATABLE, DIMENSION(:,:  ) :: zsidep, zironice
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )  CALL timing_start('p4z_sed')
      !
      IF( kt == nittrc000 .AND. knt == 1 )   THEN
          r1_rday  = 1. / rday
          IF (ln_sediment .AND. ln_sed_2way) THEN
             lk_sed = .TRUE.
          ELSE
             lk_sed = .FALSE.
          ENDIF
      ENDIF
      !
      IF( kt == nittrc000 .AND. knt == 1 )   r1_rday  = 1. / rday
      !
      ! Allocate temporary workspace
      ALLOCATE( ztrpo4(jpi,jpj,jpk) )
      IF( ln_p5z )    ALLOCATE( ztrdop(jpi,jpj,jpk) )

      zdenit2d(:,:) = 0.e0
      zbureff (:,:) = 0.e0
      zwork   (:,:) = 0.e0
      zsedsi  (:,:) = 0.e0
      zsedcal (:,:) = 0.e0
      zsedc   (:,:) = 0.e0
      
      zseddiss(:,:,:) = 0.e0
      zseddiss13(:,:,:) = 0.e0
      zsedremin(:,:,:) = 0.e0
      zsedremin13(:,:,:) = 0.e0

      ! Iron input/uptake due to sea ice : Crude parameterization based on Lancelot et al.
      ! ----------------------------------------------------
      IF( ln_ironice ) THEN  
         !                                              
         ALLOCATE( zironice(jpi,jpj) )
         !                                              
         DO jj = 1, jpj
            DO ji = 1, jpi
               zdep    = rfact2 / e3t_n(ji,jj,1)
               zwflux  = fmmflx(ji,jj) / 1000._wp
               zironice(ji,jj) =  MAX( -0.99 * trb(ji,jj,1,jpfer), -zwflux * icefeinput * zdep )
            END DO
         END DO
         !
         tra(:,:,1,jpfer) = tra(:,:,1,jpfer) + zironice(:,:) 
         ! 
         IF( lk_iomput .AND. knt == nrdttrc .AND. iom_use( "Ironice" ) )   &
            &   CALL iom_put( "Ironice", zironice(:,:) * 1.e+3 * rfact2r * e3t_n(:,:,1) * tmask(:,:,1) ) ! iron flux from ice
         !
         DEALLOCATE( zironice )
         !                                              
      ENDIF

      ! Add the external input of nutrients from dust deposition
      ! ----------------------------------------------------------
      IF( ln_dust ) THEN
         !                                              
         ALLOCATE( zsidep(jpi,jpj), zpdep(jpi,jpj,jpk), zirondep(jpi,jpj,jpk) )
         !                                              ! Iron and Si deposition at the surface
         IF( ln_solub ) THEN
            zirondep(:,:,1) = solub(:,:) * dust(:,:) * mfrac * rfact2 / e3t_n(:,:,1) / 55.85 + 3.e-10 * r1_ryyss 
         ELSE
            zirondep(:,:,1) = dustsolub  * dust(:,:) * mfrac * rfact2 / e3t_n(:,:,1) / 55.85 + 3.e-10 * r1_ryyss 
         ENDIF
         zsidep(:,:)   = 8.8 * 0.075 * dust(:,:) * mfrac * rfact2 / e3t_n(:,:,1) / 28.1 
         zpdep (:,:,1) = 0.1 * 0.021 * dust(:,:) * mfrac * rfact2 / e3t_n(:,:,1) / 31. / po4r 
         !                                              ! Iron solubilization of particles in the water column
         !                                              ! dust in kg/m2/s ---> 1/55.85 to put in mol/Fe ;  wdust in m/j
         zwdust = 0.03 * rday / ( wdust * 55.85 ) / ( 270. * rday )
         DO jk = 2, jpkm1
            zirondep(:,:,jk) = dust(:,:) * mfrac * zwdust * rfact2 * EXP( -gdept_n(:,:,jk) / 540. )
            zpdep   (:,:,jk) = zirondep(:,:,jk) * 0.023
         END DO
         !                                              ! Iron solubilization of particles in the water column
         tra(:,:,1,jpsil) = tra(:,:,1,jpsil) + zsidep  (:,:)
         DO jk = 1, jpkm1
            tra(:,:,jk,jppo4) = tra(:,:,jk,jppo4) + zpdep   (:,:,jk)
            tra(:,:,jk,jpfer) = tra(:,:,jk,jpfer) + zirondep(:,:,jk) 
         ENDDO
         ! 
         IF( lk_iomput ) THEN
            IF( knt == nrdttrc ) THEN
                IF( iom_use( "Irondep" ) )   &
                &  CALL iom_put( "Irondep", zirondep(:,:,1) * 1.e+3 * rfact2r * e3t_n(:,:,1) * tmask(:,:,1) ) ! surface downward dust depo of iron
                IF( iom_use( "pdust" ) )   &
                &  CALL iom_put( "pdust"  , dust(:,:) / ( wdust * rday )  * tmask(:,:,1) ) ! dust concentration at surface
            ENDIF
         ENDIF
         DEALLOCATE( zsidep, zpdep, zirondep )
         !                                              
      ENDIF
     
      ! Add the external input of nutrients from river
      ! ----------------------------------------------------------
      IF( ln_river ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               DO jk = 1, nk_rnf(ji,jj)
                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) +  rivdip(ji,jj) * rfact2
                  tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) +  rivdin(ji,jj) * rfact2
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) +  rivdic(ji,jj) * 5.e-5 * rfact2
                  tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) +  rivdsi(ji,jj) * rfact2
                  tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) +  rivdic(ji,jj) * rfact2
                  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) +  ( rivalk(ji,jj) - rno3 * rivdin(ji,jj) ) * rfact2
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) +  rivdoc(ji,jj) * rfact2

                  IF ( ln_n15 ) THEN
                     tra(ji,jj,jk,jp15no3) = tra(ji,jj,jk,jp15no3) + (1. + d15n_riv*1e-3)*rivdin(ji,jj) * rfact2
                     tra(ji,jj,jk,jp15doc) = tra(ji,jj,jk,jp15doc) + (1. + d15n_riv*1e-3)*rivdoc(ji,jj) * rfact2
                  ENDIF 
                  IF ( ln_c13 ) THEN
                     tra(ji,jj,jk,jp13dic) = tra(ji,jj,jk,jp13dic) + (1. + d13c_rivdic*1e-3)*rivdic(ji,jj) * rfact2
                     tra(ji,jj,jk,jp13doc) = tra(ji,jj,jk,jp13doc) + (1. + d13c_rivdoc*1e-3)*rivdoc(ji,jj) * rfact2
                  ENDIF 

               ENDDO
            ENDDO
         ENDDO
         IF (ln_ligand) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi
                  DO jk = 1, nk_rnf(ji,jj)
                     tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) +  rivdic(ji,jj) * 5.e-5 * rfact2
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
         IF( ln_p5z ) THEN
            DO jj = 1, jpj
               DO ji = 1, jpi
                  DO jk = 1, nk_rnf(ji,jj)
                     tra(ji,jj,jk,jpdop) = tra(ji,jj,jk,jpdop) + rivdop(ji,jj) * rfact2
                     tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + rivdon(ji,jj) * rfact2
                  ENDDO
               ENDDO
            ENDDO
         ENDIF
      ENDIF
      
      ! Add the external input of nutrients from nitrogen deposition
      ! ----------------------------------------------------------
      IF( ln_ndepo ) THEN
         tra(:,:,1,jpno3) = tra(:,:,1,jpno3) + nitdep(:,:) * rfact2
         tra(:,:,1,jptal) = tra(:,:,1,jptal) - rno3 * nitdep(:,:) * rfact2

         IF ( ln_n15 ) THEN
            !! pjb
            !DO jj = 1, jpj
            !   DO ji = 1, jpi
            !      d15n_dep = ( (trb(ji,jj,1,jp15no3)+rtrn) / (trb(ji,jj,1,jpno3)+rtrn)-1.0 )*1e3
            !      tra(ji,jj,1,jp15no3) = tra(ji,jj,1,jp15no3) + (1. + d15n_dep*1e-3)*nitdep(ji,jj) * rfact2
            !   ENDDO
            !ENDDO
            !! pjb
            tra(:,:,1,jp15no3) = tra(:,:,1,jp15no3) + (1. + d15n_dep*1e-3)*nitdep(:,:) * rfact2
         ENDIF

      ENDIF

      ! Add the external input of iron from hydrothermal vents
      ! ------------------------------------------------------
      IF( ln_hydrofe ) THEN
            tra(:,:,:,jpfer) = tra(:,:,:,jpfer) + hydrofe(:,:,:) * rfact2
         IF( ln_ligand ) THEN
            tra(:,:,:,jplgw) = tra(:,:,:,jplgw) + ( hydrofe(:,:,:) * lgw_rath ) * rfact2
         ENDIF
         !
         IF( lk_iomput .AND. knt == nrdttrc .AND. iom_use( "HYDR" ) )   &
            &   CALL iom_put( "HYDR", hydrofe(:,:,:) * 1.e+3 * tmask(:,:,:) ) ! hydrothermal iron input
      ENDIF

      ! OA: Warning, the following part is necessary to avoid CFL problems above the sediments
      ! --------------------------------------------------------------------
      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt  = mbkt(ji,jj)
            zdep = e3t_n(ji,jj,ikt) / xstep
            zwsbio4(ji,jj) = MIN( 0.99 * zdep, wsbio4(ji,jj,ikt) ) ! final sinking rate 
            zwsbio3(ji,jj) = MIN( 0.99 * zdep, wsbio3(ji,jj,ikt) )
         END DO
      END DO
      !
      IF( .NOT.lk_sed ) THEN
!
         ! Add the external input of iron from sediment mobilization
         ! ------------------------------------------------------
         IF( ln_ironsed ) THEN
                            tra(:,:,:,jpfer) = tra(:,:,:,jpfer) + ironsed(:,:,:) * rfact2
            !
            IF( lk_iomput .AND. knt == nrdttrc .AND. iom_use( "Ironsed" ) )   &
               &   CALL iom_put( "Ironsed", ironsed(:,:,:) * 1.e+3 * tmask(:,:,:) ) ! iron inputs from sediments
         ENDIF

         ! Computation of the sediment denitrification proportion: The metamodel from midlleburg (2006) is being used
         ! Computation of the fraction of organic matter that is permanently buried from Dunne's model
         ! -------------------------------------------------------
         DO jj = 1, jpj
            DO ji = 1, jpi
              IF( tmask(ji,jj,1) == 1 ) THEN
                 ikt = mbkt(ji,jj) ! final k level of i,j grid
                 zflx = (  trb(ji,jj,ikt,jpgoc) * zwsbio4(ji,jj)   &
                   &     + trb(ji,jj,ikt,jppoc) * zwsbio3(ji,jj) )  * 1E3 * 1E6 / 1E4
                   ! zflx = incident POM to sediments (umol C / cm2 / s) 
                   !        (molC/L * m/s) * 1e3 * 1e6 = (umolC/m2/s)
                   !        (umolC/m2/s) / 1e4 = umolC/cm2/s
                 zflx  = LOG10( MAX( 1E-3, zflx ) )  ! minimum rate of 0.001 umol C / cm2 / s 
                 zo2   = LOG10( MAX( 10. , trb(ji,jj,ikt,jpoxy) * 1E6 ) )
                 zno3  = LOG10( MAX( 1.  , trb(ji,jj,ikt,jpno3) * 1E6 * rno3 ) )
                 zdep  = LOG10( gdepw_n(ji,jj,ikt+1) )
                 zdenit2d(ji,jj) = -2.2567 - 1.185 * zflx - 0.221 * zflx**2 - 0.3995 * zno3 * zo2 + 1.25 * zno3    &
                   &                + 0.4721 * zo2 - 0.0996 * zdep + 0.4256 * zflx * zo2
                 zdenit2d(ji,jj) = 10.0**( zdenit2d(ji,jj) ) ! proportion of POM denitrified in the sediments
                   !
                 zflx = (  trb(ji,jj,ikt,jpgoc) * zwsbio4(ji,jj)   &
                   &     + trb(ji,jj,ikt,jppoc) * zwsbio3(ji,jj) ) * 1E6
                   ! zflx = incident POM to sediments (mmol C / m2 / s) 
                   !        (molC/L * m/s) * 1e6 = (mmol C / m2 / s)

                 zbureff(ji,jj) = 0.013 + 0.53 * zflx**2 / ( 7.0 + zflx )**2
              ENDIF
            END DO
         END DO 
         !
      ENDIF

      ! This loss is scaled at each bottom grid cell for equilibrating the total budget of silica in the ocean.
      ! Thus, the amount of silica lost in the sediments equal the supply at the surface (dust+rivers)
      ! ------------------------------------------------------
      IF( .NOT.lk_sed )  zrivsil = 1._wp - sedsilfrac

      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt  = mbkt(ji,jj)
            zdep = xstep / e3t_n(ji,jj,ikt) 
            zwsc = zwsbio4(ji,jj) * zdep
            zsiloss = trb(ji,jj,ikt,jpgsi) * zwsc
            zcaloss = trb(ji,jj,ikt,jpcal) * zwsc
            !
            tra(ji,jj,ikt,jpgsi) = tra(ji,jj,ikt,jpgsi) - zsiloss
            tra(ji,jj,ikt,jpcal) = tra(ji,jj,ikt,jpcal) - zcaloss
            IF ( ln_c13 ) THEN
               zr13_cal = ( (trb(ji,jj,ikt,jp13cal)+rtrn) / (trb(ji,jj,ikt,jpcal)+rtrn) )
               tra(ji,jj,ikt,jp13cal) = tra(ji,jj,ikt,jp13cal) - zcaloss * zr13_cal
            ENDIF
         END DO
      END DO
      !
      IF( .NOT.lk_sed ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               ikt  = mbkt(ji,jj)
               zdep = xstep / e3t_n(ji,jj,ikt) 
               zwsc = zwsbio4(ji,jj) * zdep
               zsiloss = trb(ji,jj,ikt,jpgsi) * zwsc
               zcaloss = trb(ji,jj,ikt,jpcal) * zwsc
               tra(ji,jj,ikt,jpsil) = tra(ji,jj,ikt,jpsil) + zsiloss * zrivsil 
               !
               zfactcal = MIN( excess(ji,jj,ikt), 0.2 )
               zfactcal = MIN( 1., 1.3 * ( 0.2 - zfactcal ) / ( 0.4 - zfactcal ) )
               zrivalk  = sedcalfrac * zfactcal
               tra(ji,jj,ikt,jptal) =  tra(ji,jj,ikt,jptal) + zcaloss * zrivalk * 2.0
               tra(ji,jj,ikt,jpdic) =  tra(ji,jj,ikt,jpdic) + zcaloss * zrivalk
               zsedcal(ji,jj) = (1.0 - zrivalk) * zcaloss * e3t_n(ji,jj,ikt) 
               zsedsi (ji,jj) = (1.0 - zrivsil) * zsiloss * e3t_n(ji,jj,ikt) 
               zseddiss(ji,jj,ikt) = zrivalk * zcaloss 
               IF ( ln_c13 ) THEN
                  zr13_cal = ( (trb(ji,jj,ikt,jp13cal)+rtrn) / (trb(ji,jj,ikt,jpcal)+rtrn) )
                  tra(ji,jj,ikt,jp13dic) =  tra(ji,jj,ikt,jp13dic) + zcaloss * zrivalk * zr13_cal
                  zseddiss13(ji,jj,ikt) = zcaloss * zrivalk * zr13_cal
               ENDIF
            END DO
         END DO
      ENDIF
      ! remove the POC and GOC that hits sediments from the POC and GOC pools
      DO jj = 1, jpj
         DO ji = 1, jpi
            ikt  = mbkt(ji,jj)
            zdep = xstep / e3t_n(ji,jj,ikt) ! timesteps/thickness
            zws4 = zwsbio4(ji,jj) * zdep ! m/s * s/m = proportion of large particles in last depth level hitting sediment
            zws3 = zwsbio3(ji,jj) * zdep ! proportion of small particles in last depth level hitting sediment
            tra(ji,jj,ikt,jpgoc) = tra(ji,jj,ikt,jpgoc) - trb(ji,jj,ikt,jpgoc) * zws4 ! remove GOC from water column
            tra(ji,jj,ikt,jppoc) = tra(ji,jj,ikt,jppoc) - trb(ji,jj,ikt,jppoc) * zws3 ! remove POC from water column
            tra(ji,jj,ikt,jpbfe) = tra(ji,jj,ikt,jpbfe) - trb(ji,jj,ikt,jpbfe) * zws4 ! remove BFe from water column
            tra(ji,jj,ikt,jpsfe) = tra(ji,jj,ikt,jpsfe) - trb(ji,jj,ikt,jpsfe) * zws3 ! remove SFe from water column

            IF ( ln_n15 ) THEN
               tra(ji,jj,ikt,jp15goc) = tra(ji,jj,ikt,jp15goc) - trb(ji,jj,ikt,jp15goc) * zws4 ! remove previous GOC_15 
               tra(ji,jj,ikt,jp15poc) = tra(ji,jj,ikt,jp15poc) - trb(ji,jj,ikt,jp15poc) * zws3 ! remove previous POC_15
            ENDIF
            IF ( ln_c13 ) THEN
               tra(ji,jj,ikt,jp13goc) = tra(ji,jj,ikt,jp13goc) - trb(ji,jj,ikt,jp13goc) * zws4 ! remove previous GOC_13 
               tra(ji,jj,ikt,jp13poc) = tra(ji,jj,ikt,jp13poc) - trb(ji,jj,ikt,jp13poc) * zws3 ! remove previous POC_13
            ENDIF

         END DO
      END DO
      !
      IF( ln_p5z ) THEN
         DO jj = 1, jpj
            DO ji = 1, jpi
               ikt  = mbkt(ji,jj)
               zdep = xstep / e3t_n(ji,jj,ikt) 
               zws4 = zwsbio4(ji,jj) * zdep
               zws3 = zwsbio3(ji,jj) * zdep
               tra(ji,jj,ikt,jpgon) = tra(ji,jj,ikt,jpgon) - trb(ji,jj,ikt,jpgon) * zws4
               tra(ji,jj,ikt,jppon) = tra(ji,jj,ikt,jppon) - trb(ji,jj,ikt,jppon) * zws3
               tra(ji,jj,ikt,jpgop) = tra(ji,jj,ikt,jpgop) - trb(ji,jj,ikt,jpgop) * zws4
               tra(ji,jj,ikt,jppop) = tra(ji,jj,ikt,jppop) - trb(ji,jj,ikt,jppop) * zws3
            END DO
         END DO
      ENDIF

      IF( .NOT.lk_sed ) THEN
         ! The 0.5 factor in zpdenit is to avoid negative NO3 concentration after
         ! denitrification in the sediments. Not very clever, but simpliest option.
         DO jj = 1, jpj
            DO ji = 1, jpi
               zpdenit(ji,jj,:) = 0.0 ! set all values in water column to zero
               zpdenit15(ji,jj,:) = 0.0 ! set all values in water column to zero
               zolimit(ji,jj,:) = 0.0 ! set all values in water column to zero
               zolimit15(ji,jj,:) = 0.0 ! set all values in water column to zero
               ikt  = mbkt(ji,jj) ! we are in deepest grid cell of ocean
               zdep = xstep / e3t_n(ji,jj,ikt) ! s/m (timestep/thickness) 
               zws4 = zwsbio4(ji,jj) * zdep ! proportion of large particles in last depth level hitting sediment
               zws3 = zwsbio3(ji,jj) * zdep ! proportion of small particles in last depth level hitting sediment 
               zrivno3 = 1. - zbureff(ji,jj) ! proportion of sedimenting particles not buried
               zwstpoc = trb(ji,jj,ikt,jpgoc) * zws4 + trb(ji,jj,ikt,jppoc) * zws3 !molC/L 
                 ! zwstpoc = OM * proportion POC/GOC sedimented = molC/L of GOC+POC hitting sediment
               zpdenit(ji,jj,ikt)  = MIN( 0.5 * ( trb(ji,jj,ikt,jpno3) - rtrn ) / rdenit, zdenit2d(ji,jj) * zwstpoc * zrivno3 )
               z1pdenit = zwstpoc * zrivno3 - zpdenit(ji,jj,ikt)
                 ! zrivno3 = proportion of particulate carbon that isn't buried and must be remineralised 
                 ! zpdenit = amount of material (POC+GOC) that is denitrified (molC/L)
                 ! z1pdenit = amount of material (POC+GOC) that is remineralised without denitrification (molC/L)
               zolimit(ji,jj,ikt) = MIN( ( trb(ji,jj,ikt,jpoxy) - rtrn ) / o2ut, z1pdenit * ( 1.- nitrfac(ji,jj,ikt) ) )
                 ! zolimit = adjusted z1pdenit that occurs through oxic remineralisation

               ! zpdenit + zolimit + (z1pdenit-zolimit) = zwstpoc * zrivno3
               ! OM_denit + OM_oxy + OM_other_remin     = OM_total_hitting_sed_that_is_not_buried

               tra(ji,jj,ikt,jpdoc) = tra(ji,jj,ikt,jpdoc) + z1pdenit - zolimit(ji,jj,ikt)
                 ! z1pdenit (non-nitrogeneous remineralisation) - zolimit (oxic remineralisation)
                 !  thus... the component of remineralisation that cannot be accommodated by N or O2 goes to DOC
               tra(ji,jj,ikt,jppo4) = tra(ji,jj,ikt,jppo4) + zpdenit(ji,jj,ikt) + zolimit(ji,jj,ikt)
               tra(ji,jj,ikt,jpnh4) = tra(ji,jj,ikt,jpnh4) + zpdenit(ji,jj,ikt) + zolimit(ji,jj,ikt)
               tra(ji,jj,ikt,jpno3) = tra(ji,jj,ikt,jpno3) - rdenit * zpdenit(ji,jj,ikt)
          
               tra(ji,jj,ikt,jpoxy) = tra(ji,jj,ikt,jpoxy) - zolimit(ji,jj,ikt) * o2ut
               tra(ji,jj,ikt,jptal) = tra(ji,jj,ikt,jptal) + rno3 * (zolimit(ji,jj,ikt) + (1.+rdenit) * zpdenit(ji,jj,ikt) )
               tra(ji,jj,ikt,jpdic) = tra(ji,jj,ikt,jpdic) + zpdenit(ji,jj,ikt) + zolimit(ji,jj,ikt) 
               zsedremin(ji,jj,ikt) = zpdenit(ji,jj,ikt) + zolimit(ji,jj,ikt)
               sdenit(ji,jj) = rdenit * zpdenit(ji,jj,ikt) * e3t_n(ji,jj,ikt)
               zsedc(ji,jj)   = (1. - zrivno3) * zwstpoc * e3t_n(ji,jj,ikt)

               IF( ln_p5z ) THEN
                  zwstpop              = trb(ji,jj,ikt,jpgop) * zws4 + trb(ji,jj,ikt,jppop) * zws3
                  zwstpon              = trb(ji,jj,ikt,jpgon) * zws4 + trb(ji,jj,ikt,jppon) * zws3
                  tra(ji,jj,ikt,jpdon) = tra(ji,jj,ikt,jpdon) + ( z1pdenit - zolimit(ji,jj,ikt) ) * zwstpon / (zwstpoc + rtrn)
                  tra(ji,jj,ikt,jpdop) = tra(ji,jj,ikt,jpdop) + ( z1pdenit - zolimit(ji,jj,ikt) ) * zwstpop / (zwstpoc + rtrn)
               ENDIF

               IF ( ln_n15 ) THEN
                  zr15_no3 = ( (trb(ji,jj,ikt,jp15no3)+rtrn) / (trb(ji,jj,ikt,jpno3)+rtrn) )
                  zwstpoc15 = trb(ji,jj,ikt,jp15goc)*zws4 + trb(ji,jj,ikt,jp15poc)*zws3 !POC+GOC hitting sediment 
                  zr15_rain = ( (zwstpoc15+rtrn) / (zwstpoc+rtrn) )
                  tra(ji,jj,ikt,jp15doc) = tra(ji,jj,ikt,jp15doc) + zwstpoc15 * zrivno3           &
                  &                        - zpdenit(ji,jj,ikt) * zr15_rain - zolimit(ji,jj,ikt) * zr15_rain
                  tra(ji,jj,ikt,jp15nh4) = tra(ji,jj,ikt,jp15nh4) + zpdenit(ji,jj,ikt) * zr15_rain + zolimit(ji,jj,ikt) * zr15_rain
                  tra(ji,jj,ikt,jp15no3) = tra(ji,jj,ikt,jp15no3) - rdenit * zpdenit(ji,jj,ikt) * ( 1.0 - e15n_ben/1000.0 ) * zr15_no3 
                  zpdenit15(ji,jj,ikt) = zpdenit(ji,jj,ikt) * ( 1.0 - e15n_ben/1000.0 ) * zr15_no3 
                  zolimit15(ji,jj,ikt) = zolimit(ji,jj,ikt) * zr15_rain
               ENDIF 
               IF ( ln_c13 ) THEN
                  zwstpoc13 = trb(ji,jj,ikt,jp13goc)*zws4 + trb(ji,jj,ikt,jp13poc)*zws3 !POC+GOC hitting sediment 
                  zr13_rain = ( (zwstpoc13+rtrn) / (zwstpoc+rtrn) )
                  tra(ji,jj,ikt,jp13doc) = tra(ji,jj,ikt,jp13doc) + zwstpoc13 * zrivno3           &
                  &                        - zpdenit(ji,jj,ikt) * zr13_rain - zolimit(ji,jj,ikt) * zr13_rain
                  tra(ji,jj,ikt,jp13dic) = tra(ji,jj,ikt,jp13dic) + zpdenit(ji,jj,ikt) * zr13_rain + zolimit(ji,jj,ikt) * zr13_rain
                  zsedremin13(ji,jj,ikt) = zpdenit(ji,jj,ikt) * zr13_rain + zolimit(ji,jj,ikt) * zr13_rain
               ENDIF 

            END DO
         END DO
       ENDIF


      ! Nitrogen fixation process
      ! Small source iron from particulate inorganic iron
      !-----------------------------------
      DO jk = 1, jpkm1
         zlight (:,:,jk) =  ( 1.- EXP( -etot_ndcy(:,:,jk) / diazolight ) ) * ( 1. - fr_i(:,:) ) 
         zsoufer(:,:,jk) = zlight(:,:,jk) * 2E-11 / ( 2E-11 + biron(:,:,jk) )
      ENDDO
      IF( ln_p4z ) THEN
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! Potential nitrogen fixation dependant on temperature and iron
                  ztemp = tsn(ji,jj,jk,jp_tem) ! in situ temperature !pjb
                  zmudia = MAX( 0.,-0.001096*ztemp**2 + 0.057*ztemp -0.637 ) * 7.625 ! growth rate (T)
                  ! Potential nitrogen fixation dependant on temperature and iron
                  xdianh4 = trb(ji,jj,jk,jpnh4) / ( concnnh4 + trb(ji,jj,jk,jpnh4) )
                  xdiano3 = trb(ji,jj,jk,jpno3) / ( concnno3 + trb(ji,jj,jk,jpno3) ) * (1. - xdianh4)
                  zlim = ( 1.- xdiano3 - xdianh4 )
                  IF( zlim <= 0.1 )   zlim = 0.01   ! too much NH4 and NO3 and so N2fixation is minimal 
                  zfact = zlim * rfact2 
                    ! rfact2 = seconds per timestep (21600 secs)
                    ! zlim = limitation of N2 fixation [0,1] (zlim --> when NO3 and NH4 are very low )
                  ztrfer = biron(ji,jj,jk) / ( concfediaz + biron(ji,jj,jk) )
                    ! biron = bioavailble iron
                  ztrpo4(ji,jj,jk) = trb(ji,jj,jk,jppo4) / ( concphdiaz + trb(ji,jj,jk,jppo4) )
                  ztrdp = ztrpo4(ji,jj,jk)
                  nitrpot(ji,jj,jk) =  zmudia * r1_rday * zfact * MIN( ztrfer, ztrdp ) * zlight(ji,jj,jk)
                    ! r1_rday = seconds per day (1/86400)
                    ! zmudia = max growth rate set by temperature * r1_rday = growth rate per sec
                    ! zfact = nitrogen limitation multiplied by number of seconds per timestep
                    ! MIN( ztrfer, ztrdp ) = leibig law of minimum nutrient limitation
                    ! zlight(ji,jj,jk) = dependence on light [0,1] (low light --> 0)
               END DO
            END DO
         END DO
      ELSE       ! p5z
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  !                      ! Potential nitrogen fixation dependant on temperature and iron
                  ztemp = tsn(ji,jj,jk,jp_tem)
                  zmudia = MAX( 0.,-0.001096*ztemp**2 + 0.057*ztemp -0.637 ) * 7.625
                  !       Potential nitrogen fixation dependant on temperature and iron
                  xdianh4 = trb(ji,jj,jk,jpnh4) / ( concnnh4 + trb(ji,jj,jk,jpnh4) )
                  xdiano3 = trb(ji,jj,jk,jpno3) / ( concnno3 + trb(ji,jj,jk,jpno3) ) * (1. - xdianh4)
                  zlim = ( 1.- xdiano3 - xdianh4 )
                  IF( zlim <= 0.1 )   zlim = 0.01
                  zfact = zlim * rfact2
                  ztrfer = biron(ji,jj,jk) / ( concfediaz + biron(ji,jj,jk) )
                  ztrpo4(ji,jj,jk) = trb(ji,jj,jk,jppo4) / ( 1E-6 + trb(ji,jj,jk,jppo4) )
                  ztrdop(ji,jj,jk) = trb(ji,jj,jk,jpdop) / ( 1E-6 + trb(ji,jj,jk,jpdop) ) * (1. - ztrpo4(ji,jj,jk))
                  ztrdp = ztrpo4(ji,jj,jk) + ztrdop(ji,jj,jk)
                  nitrpot(ji,jj,jk) =  zmudia * r1_rday * zfact * MIN( ztrfer, ztrdp ) * zlight(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF

      ! Nitrogen change due to nitrogen fixation
      ! ----------------------------------------
      IF( ln_p4z ) THEN
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zfact = nitrpot(ji,jj,jk) * nitrfix
                  tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zfact / 3.0
                  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * zfact / 3.0
                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) - zfact * 2.0 / 3.0
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zfact * 1.0 / 3.0
                  tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zfact * 1.0 / 3.0 * 2.0 / 3.0
                  tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zfact * 1.0 / 3.0 * 1.0 / 3.0
                  tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + ( o2ut + o2nit ) * zfact * 2.0 / 3.0 + o2nit * zfact / 3.0
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - 30E-6 * zfact * 1.0 / 3.0
                  tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + 30E-6 * zfact * 1.0 / 3.0 * 2.0 / 3.0
                  tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + 30E-6 * zfact * 1.0 / 3.0 * 1.0 / 3.0
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + 0.002 * 4E-10 * zsoufer(ji,jj,jk) * rfact2 / rday
                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + concdnh4 / ( concdnh4 + trb(ji,jj,jk,jppo4) ) &
                  &                     * 0.001 * trb(ji,jj,jk,jpdoc) * xstep
                   
                    ! Divisions are routing nitrogen fixed matter directly to
                    ! different pools:
                    !   one third to NH4, 
                    !   one third to DOC, 
                    !   two ninths to POC,
                    !   one ninths to GOC 

                  IF ( ln_n15 ) THEN
                     !! pjb
                     !d15n_fix = ( (trb(ji,jj,jk,jp15no3)+rtrn) / (trb(ji,jj,jk,jpno3)+rtrn) - 1.0 )*1e3
                     !! pjb
                     tra(ji,jj,jk,jp15nh4) = tra(ji,jj,jk,jp15nh4) + zfact * (1. + d15n_fix*1e-3) * 1./3.  
                     tra(ji,jj,jk,jp15doc) = tra(ji,jj,jk,jp15doc) + zfact * (1. + d15n_fix*1e-3) * 1./3.
                     tra(ji,jj,jk,jp15poc) = tra(ji,jj,jk,jp15poc) + zfact * (1. + d15n_fix*1e-3) * 2./9.
                     tra(ji,jj,jk,jp15goc) = tra(ji,jj,jk,jp15goc) + zfact * (1. + d15n_fix*1e-3) * 1./9.
                  ENDIF
                  IF ( ln_c13 ) THEN
                     !! pjb
                     !d13c_fix = ( (trb(ji,jj,jk,jp13dic)+rtrn) / (trb(ji,jj,jk,jpdic)+rtrn) - 1.0 )*1e3
                     !! pjb
                     tra(ji,jj,jk,jp13doc) = tra(ji,jj,jk,jp13doc) + zfact * (1. + d13c_fix*1e-3) * 1./3.
                     tra(ji,jj,jk,jp13poc) = tra(ji,jj,jk,jp13poc) + zfact * (1. + d13c_fix*1e-3) * 2./9.
                     tra(ji,jj,jk,jp13goc) = tra(ji,jj,jk,jp13goc) + zfact * (1. + d13c_fix*1e-3) * 1./9.
                  ENDIF

              END DO
            END DO 
         END DO
      ELSE    ! p5z
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zfact = nitrpot(ji,jj,jk) * nitrfix
                  tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zfact / 3.0
                  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * zfact / 3.0
                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) - 16.0 / 46.0 * zfact * ( 1.0 - 1.0 / 3.0 ) &
                  &                     * ztrpo4(ji,jj,jk) / (ztrpo4(ji,jj,jk) + ztrdop(ji,jj,jk) + rtrn)
                  tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) + zfact * 1.0 / 3.0
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zfact * 1.0 / 3.0
                  tra(ji,jj,jk,jpdop) = tra(ji,jj,jk,jpdop) + 16.0 / 46.0 * zfact / 3.0  &
                  &                     - 16.0 / 46.0 * zfact * ztrdop(ji,jj,jk)   &
                  &                     / (ztrpo4(ji,jj,jk) + ztrdop(ji,jj,jk) + rtrn)
                  tra(ji,jj,jk,jppoc) = tra(ji,jj,jk,jppoc) + zfact * 1.0 / 3.0 * 2.0 / 3.0
                  tra(ji,jj,jk,jppon) = tra(ji,jj,jk,jppon) + zfact * 1.0 / 3.0 * 2.0 /3.0
                  tra(ji,jj,jk,jppop) = tra(ji,jj,jk,jppop) + 16.0 / 46.0 * zfact * 1.0 / 3.0 * 2.0 /3.0
                  tra(ji,jj,jk,jpgoc) = tra(ji,jj,jk,jpgoc) + zfact * 1.0 / 3.0 * 1.0 / 3.0
                  tra(ji,jj,jk,jpgon) = tra(ji,jj,jk,jpgon) + zfact * 1.0 / 3.0 * 1.0 /3.0
                  tra(ji,jj,jk,jpgop) = tra(ji,jj,jk,jpgop) + 16.0 / 46.0 * zfact * 1.0 / 3.0 * 1.0 /3.0
                  tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + ( o2ut + o2nit ) * zfact * 2.0 / 3.0 + o2nit * zfact / 3.0
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - 30E-6 * zfact * 1.0 / 3.0 
                  tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + 30E-6 * zfact * 1.0 / 3.0 * 2.0 / 3.0
                  tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + 30E-6 * zfact * 1.0 / 3.0 * 1.0 / 3.0
                  tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) + 0.002 * 4E-10 * zsoufer(ji,jj,jk) * rfact2 / rday
              END DO
            END DO 
         END DO
         !
      ENDIF

      IF( lk_iomput ) THEN
         IF( knt == nrdttrc ) THEN
            zfact = 1.e+3 * rfact2r !  conversion from molC/l/kt  to molC/m3/s
            IF( iom_use("Nfix"   ) ) CALL iom_put( "Nfix", nitrpot(:,:,:) * nitrfix * rno3 * zfact * tmask(:,:,:) )  ! nitrogen fixation 
            IF( iom_use("Nfix_15NH4"   ) ) CALL iom_put( "Nfix_15NH4", nitrpot(:,:,:) * nitrfix * rno3 * zfact * tmask(:,:,:) * (1. + d15n_fix*1e-3) * 1./3. )  ! nitrogen fixation to 15NH4
            IF( iom_use("INTNFIX") ) THEN   ! nitrogen fixation rate in ocean ( vertically integrated )
               zwork(:,:) = 0.
               DO jk = 1, jpkm1
                 zwork(:,:) = zwork(:,:) + nitrpot(:,:,jk) * nitrfix * rno3 * zfact * e3t_n(:,:,jk) * tmask(:,:,jk)
               ENDDO
               CALL iom_put( "INTNFIX" , zwork ) 
            ENDIF
            IF( iom_use("SedCal" ) ) CALL iom_put( "SedCal", zsedcal(:,:) * zfact )
            IF( iom_use("SedSi" ) )  CALL iom_put( "SedSi",  zsedsi (:,:) * zfact )
            IF( iom_use("SedC" ) )   CALL iom_put( "SedC",   zsedc  (:,:) * zfact )
            IF( iom_use("Sdenit" ) ) CALL iom_put( "Sdenit", sdenit (:,:) * zfact * rno3 )
            IF( iom_use("SedDiss" ) ) CALL iom_put( "SedDiss", zseddiss(:,:,:) * zfact * tmask(:,:,:) )
            IF( iom_use("SedDiss_13CAL" ) ) CALL iom_put( "SedDiss_13CAL", zseddiss13(:,:,:) * zfact * tmask(:,:,:) )
            IF( iom_use("SedReminC" ) ) CALL iom_put( "SedReminC", zsedremin(:,:,:) * zfact * tmask(:,:,:) )
            IF( iom_use("SedReminC_13POM" ) ) CALL iom_put( "SedReminC_13POM", zsedremin13(:,:,:) * zfact * tmask(:,:,:) )
            IF( iom_use("SDEN3D" ) ) CALL iom_put( "SDEN3D", zpdenit(:,:,:) * rdenit * zfact * rno3 * tmask(:,:,:) )
            IF( iom_use("SDEN3D_15NO3" ) ) CALL iom_put( "SDEN3D_15NO3", zpdenit15(:,:,:) * rdenit * zfact * rno3 * tmask(:,:,:) )
            IF( iom_use("SREM3D" ) ) CALL iom_put( "SREM3D", zolimit(:,:,:) * zfact * tmask(:,:,:) )
            IF( iom_use("SREM3D_15DOC" ) ) CALL iom_put( "SREM3D_15DOC", zolimit15(:,:,:) * zfact * tmask(:,:,:) )
            IF( iom_use("River_DIC" ) ) CALL iom_put( "River_DIC", rivdic(:,:) * zfact )
            IF( iom_use("River_13DIC" ) ) CALL iom_put( "River_13DIC", rivdic (:,:) * zfact * ( 1.0 + d13c_rivdic/1000.0 ) )
            IF( iom_use("River_NO3" ) ) CALL iom_put( "River_NO3", rivdin(:,:) * zfact )
            IF( iom_use("River_15NO3" ) ) CALL iom_put( "River_15NO3", rivdin (:,:) * zfact * ( 1.0 + d15n_riv/1000.0 ) )
            IF( iom_use("Ndep_NO3" ) ) CALL iom_put( "Ndep_NO3", nitdep(:,:) * zfact )
            IF( iom_use("Ndep_15NO3" ) ) CALL iom_put( "Ndep_15NO3", nitdep (:,:) * zfact * ( 1.0 + d15n_dep/1000.0 ) )
         ENDIF
      ENDIF
      !
      IF(ln_ctl) THEN  ! print mean trends (USEd for debugging)
         WRITE(charout, fmt="('sed ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_p5z )    DEALLOCATE( ztrpo4, ztrdop )
      !
      IF( ln_timing )  CALL timing_stop('p4z_sed')
      !
   END SUBROUTINE p4z_sed


   INTEGER FUNCTION p4z_sed_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_sed_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( nitrpot(jpi,jpj,jpk), sdenit(jpi,jpj), zpdenit(jpi,jpj,jpk), zpdenit15(jpi,jpj,jpk), STAT=p4z_sed_alloc )
      !
      IF( p4z_sed_alloc /= 0 )   CALL ctl_stop( 'STOP', 'p4z_sed_alloc: failed to allocate arrays' )
      !
   END FUNCTION p4z_sed_alloc

   !!======================================================================
END MODULE p4zsed
