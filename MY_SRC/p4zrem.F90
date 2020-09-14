MODULE p4zrem
   !!======================================================================
   !!                         ***  MODULE p4zrem  ***
   !! TOP :   PISCES Compute remineralization/dissolution of organic compounds
   !!=========================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-06  (O. Aumont, C. Ethe) Quota model for iron
   !!----------------------------------------------------------------------
   !!   p4z_rem       :  Compute remineralization/dissolution of organic compounds
   !!   p4z_rem_init  :  Initialisation of parameters for remineralisation
   !!   p4z_rem_alloc :  Allocate remineralisation variables
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zche          !  chemical model
   USE p4zprod         !  Growth rate of the 2 phyto groups
   USE p4zlim
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager


   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_rem         ! called in p4zbio.F90
   PUBLIC   p4z_rem_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_rem_alloc

   REAL(wp), PUBLIC ::   xremikc    !: remineralisation rate of DOC 
   REAL(wp), PUBLIC ::   xremikn    !: remineralisation rate of DON 
   REAL(wp), PUBLIC ::   xremikp    !: remineralisation rate of DOP 
   REAL(wp), PUBLIC ::   xremik     !: remineralisation rate of POC 
   REAL(wp), PUBLIC ::   nitrif     !: NH4 nitrification rate 
   REAL(wp), PUBLIC ::   xsirem     !: remineralisation rate of POC 
   REAL(wp), PUBLIC ::   xsiremlab  !: fast remineralisation rate of POC 
   REAL(wp), PUBLIC ::   xsilab     !: fraction of labile biogenic silica 
   REAL(wp), PUBLIC ::   feratb     !: Fe/C quota in bacteria
   REAL(wp), PUBLIC ::   xkferb     !: Half-saturation constant for bacteria Fe/C
   REAL(wp), PUBLIC ::   e15n_den   !: N15 denitrification fractionation
   REAL(wp), PUBLIC ::   e15n_nit   !: N15 nitrification fractionation

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr   !: denitrification array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zonitr   !: nitrification array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zdenitnh4!: anoxic nitrification array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zoxyremc !: anoxic remin (alternative) array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr15   !: denitrification array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   denitr15nh4!: denitrification (NH4 addition) array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zonitr15   !: nitrification array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zoxyremc15 !: anoxic remin (alternative) array
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zdenitnh415no3!: anoxic nitrification array (remove NO3)
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   zdenitnh415nh4!: anoxic nitrification array (remove NH4)

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zrem.F90 10425 2018-12-19 21:54:16Z smasson $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_rem( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem  ***
      !!
      !! ** Purpose :   Compute remineralization/scavenging of organic compounds
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt ! ocean time step
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zremik, zremikc, zremikn, zremikp, zsiremin, zfact 
      REAL(wp) ::   zsatur, zsatur2, znusil, znusil2, zdep, zdepmin, zfactdep
      REAL(wp) ::   zbactfer, zolimit, zrfact2
      REAL(wp) ::   zammonic, zoxyremn, zoxyremp
      REAL(wp) ::   zosil, ztem, zolimic, zolimin, zolimip, zdenitrn, zdenitrp
      REAL(wp) ::   zr15_doc, zr15_no3, zr15_nh4
      REAL(wp) ::   zr13_doc
      CHARACTER (len=25) :: charout
      REAL(wp), DIMENSION(jpi,jpj    ) :: ztempbac
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zdepbac, zolimi, zdepprod, zfacsi, zfacsib, zdepeff, zfebact
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zolimi15
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zremindic, zremindic13
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_rem')
      !
      ! Initialisation of arrys
      zdepprod(:,:,:) = 1._wp
      zdepeff (:,:,:) = 0.3_wp
      ztempbac(:,:)   = 0._wp
      zfacsib(:,:,:)  = xsilab / ( 1.0 - xsilab )
      zfebact(:,:,:)  = 0._wp
      zremindic(:,:,:)  = 0._wp
      zfacsi(:,:,:)   = xsilab

      IF ( ln_c13 ) THEN
         zremindic13(:,:,:)  = 0._wp
      ENDIF
      IF ( ln_n15 ) THEN
         zolimi15(:,:,:)  = 0._wp
      ENDIF

      ! Computation of the mean phytoplankton concentration as
      ! a crude estimate of the bacterial biomass
      ! this parameterization has been deduced from a model version
      ! that was modeling explicitely bacteria
      ! -------------------------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zdep = MAX( hmld(ji,jj), heup(ji,jj) ) 
               IF( gdept_n(ji,jj,jk) < zdep ) THEN  ! if depth is less than mixedlayer and euphotic zone
                  zdepbac(ji,jj,jk) = MIN( 0.7 * ( trb(ji,jj,jk,jpzoo) + 2.* trb(ji,jj,jk,jpmes) ), 4.e-6 )
                  ztempbac(ji,jj)   = zdepbac(ji,jj,jk)
               ELSE  ! if deeper than euphotic zone and mixed layer
                  zdepmin = MIN( 1., zdep / gdept_n(ji,jj,jk) ) ! fraction [0,1] (deeper --> 0)
                  zdepbac (ji,jj,jk) = zdepmin**0.683 * ztempbac(ji,jj) ! power law increase in bacteria with depth 
                  zdepprod(ji,jj,jk) = zdepmin**0.273
                  zdepeff (ji,jj,jk) = zdepeff(ji,jj,jk) * zdepmin**0.3 ! [0,0.3] (deeper --> 0)
               ENDIF
            END DO
         END DO
      END DO

      IF( ln_p4z ) THEN
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! DOC ammonification. Depends on depth, phytoplankton biomass
                  ! and a limitation term which is supposed to be a parameterization of the bacterial activity. 
                  zremik = xremik * xstep / 1.e-6 * xlimbac(ji,jj,jk) * zdepbac(ji,jj,jk)  
                  zremik = MAX( zremik, 2.74e-4 * xstep )  ! basic remineralisation rate
                  ! Ammonification in oxic waters with oxygen consumption
                  ! -----------------------------------------------------
                  zolimit = zremik * ( 1.- nitrfac(ji,jj,jk) ) * trb(ji,jj,jk,jpdoc) 
                    ! apply oxygen limitation and multiply rate by amount of DOC available
                  zolimi(ji,jj,jk) = MIN( ( trb(ji,jj,jk,jpoxy) - rtrn ) / o2ut, zolimit ) 
                    ! don't remove more O2 than is available
                  ! Ammonification in suboxic waters with denitrification
                  ! -------------------------------------------------------
                  zammonic = zremik * nitrfac(ji,jj,jk) * trb(ji,jj,jk,jpdoc) !DOC remineralised by denitrifiers
                  denitr(ji,jj,jk)  = zammonic * ( 1. - nitrfac2(ji,jj,jk) ) !apply NO3 limitation to denitrification too
                  denitr(ji,jj,jk)  = MIN( ( trb(ji,jj,jk,jpno3) - rtrn ) / rdenit, denitr(ji,jj,jk) ) 
                    ! don't remove more NO3 than is available
                  zoxyremc(ji,jj,jk) = zammonic - denitr(ji,jj,jk) 
                    ! reallocate material that can't be remineralised using NO3 to oxygen
                  !
                  zolimi(ji,jj,jk) = MAX( 0.e0, zolimi(ji,jj,jk) )  ! make sure all arrays are positive
                  denitr(ji,jj,jk) = MAX( 0.e0, denitr(ji,jj,jk) )
                  zoxyremc(ji,jj,jk) = MAX( 0.e0, zoxyremc(ji,jj,jk) )
                    ! zolimi = oxic remineralisation of DOC --> NH4 using oxygen
                    ! denitr = suboxic remineralisation of DOC --> NH4 using nitrate
                    ! zoxyremc = suboxic remineralisation of DOC --> NH4 not using nitrate 
 
                  !
                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zolimi(ji,jj,jk) + denitr(ji,jj,jk) + zoxyremc(ji,jj,jk)
                  tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zolimi(ji,jj,jk) + denitr(ji,jj,jk) + zoxyremc(ji,jj,jk)
                  tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - denitr(ji,jj,jk) * rdenit
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) - zolimi(ji,jj,jk) - denitr(ji,jj,jk) - zoxyremc(ji,jj,jk)
                  tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - zolimi(ji,jj,jk) * o2ut
                  tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zolimi(ji,jj,jk) + denitr(ji,jj,jk) + zoxyremc(ji,jj,jk)
                  zremindic(ji,jj,jk) = zolimi(ji,jj,jk) + denitr(ji,jj,jk) + zoxyremc(ji,jj,jk) 
                  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * ( zolimi(ji,jj,jk) + zoxyremc(ji,jj,jk)    &
                  &                     + ( rdenit + 1.) * denitr(ji,jj,jk) )

                  IF( ln_n15) THEN  ! should add isotope effect to ammonification ( Knapp et al 2011, Mobius 2013)
                     zr15_doc = ( (trb(ji,jj,jk,jp15doc)+rtrn) / (trb(ji,jj,jk,jpdoc)+rtrn) )
                     zr15_no3 = ( (trb(ji,jj,jk,jp15no3)+rtrn) / (trb(ji,jj,jk,jpno3)+rtrn) )
                     tra(ji,jj,jk,jp15nh4) = tra(ji,jj,jk,jp15nh4) +  & 
                     &                       ( zolimi(ji,jj,jk) + denitr(ji,jj,jk) + zoxyremc(ji,jj,jk) ) * zr15_doc
                     tra(ji,jj,jk,jp15no3) = tra(ji,jj,jk,jp15no3) - denitr(ji,jj,jk) * rdenit  &
                     &                       * ( 1.0 - e15n_den/1000.0 ) * zr15_no3
                     tra(ji,jj,jk,jp15doc) = tra(ji,jj,jk,jp15doc) -  &
                     &                       ( zolimi(ji,jj,jk) + denitr(ji,jj,jk) + zoxyremc(ji,jj,jk) ) * zr15_doc
 
                     zolimi15(ji,jj,jk) = zolimi(ji,jj,jk) * zr15_doc
                     denitr15(ji,jj,jk) = denitr(ji,jj,jk) * ( 1.0 - e15n_den/1000.0 ) * zr15_no3
                     denitr15nh4(ji,jj,jk) = denitr(ji,jj,jk) * zr15_doc
                     zoxyremc15(ji,jj,jk) = zoxyremc(ji,jj,jk) * zr15_doc
                  ENDIF
                  IF( ln_c13) THEN
                     zr13_doc = ( (trb(ji,jj,jk,jp13doc)+rtrn) / (trb(ji,jj,jk,jpdoc)+rtrn) )
                     tra(ji,jj,jk,jp13doc) = tra(ji,jj,jk,jp13doc) -  &
                     &                       ( zolimi(ji,jj,jk) + denitr(ji,jj,jk) + zoxyremc(ji,jj,jk) ) * zr13_doc 
                     tra(ji,jj,jk,jp13dic) = tra(ji,jj,jk,jp13dic) +  &
                     &                       ( zolimi(ji,jj,jk) + denitr(ji,jj,jk) + zoxyremc(ji,jj,jk) ) * zr13_doc 
                     zremindic13(ji,jj,jk) = ( zolimi(ji,jj,jk) + denitr(ji,jj,jk) + zoxyremc(ji,jj,jk) ) * zr13_doc
                  ENDIF

               END DO
            END DO
         END DO
      ELSE
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ! DOC ammonification. Depends on depth, phytoplankton biomass
                  ! and a limitation term which is supposed to be a parameterization of the bacterial activity. 
                  ! -----------------------------------------------------------------
                  zremik = xstep / 1.e-6 * MAX(0.01, xlimbac(ji,jj,jk)) * zdepbac(ji,jj,jk) 
                  zremik = MAX( zremik, 2.74e-4 * xstep / xremikc )

                  zremikc = xremikc * zremik
                  zremikn = xremikn / xremikc
                  zremikp = xremikp / xremikc

                  ! Ammonification in oxic waters with oxygen consumption
                  ! -----------------------------------------------------
                  zolimit = zremikc * ( 1.- nitrfac(ji,jj,jk) ) * trb(ji,jj,jk,jpdoc) 
                  zolimic = MAX( 0.e0, MIN( ( trb(ji,jj,jk,jpoxy) - rtrn ) / o2ut, zolimit ) ) 
                  zolimi(ji,jj,jk) = zolimic
                  zolimin = zremikn * zolimic * trb(ji,jj,jk,jpdon) / ( trb(ji,jj,jk,jpdoc) + rtrn )
                  zolimip = zremikp * zolimic * trb(ji,jj,jk,jpdop) / ( trb(ji,jj,jk,jpdoc) + rtrn ) 

                  ! Ammonification in suboxic waters with denitrification
                  ! -------------------------------------------------------
                  zammonic = zremikc * nitrfac(ji,jj,jk) * trb(ji,jj,jk,jpdoc)
                  denitr(ji,jj,jk)  = zammonic * ( 1. - nitrfac2(ji,jj,jk) )
                  denitr(ji,jj,jk)  = MAX(0., MIN(  ( trb(ji,jj,jk,jpno3) - rtrn ) / rdenit, denitr(ji,jj,jk) ) )
                  zoxyremc(ji,jj,jk) = MAX(0., zammonic - denitr(ji,jj,jk))
                  zdenitrn  = zremikn * denitr(ji,jj,jk) * trb(ji,jj,jk,jpdon) / ( trb(ji,jj,jk,jpdoc) + rtrn )
                  zdenitrp  = zremikp * denitr(ji,jj,jk) * trb(ji,jj,jk,jpdop) / ( trb(ji,jj,jk,jpdoc) + rtrn )
                  zoxyremn  = zremikn * zoxyremc(ji,jj,jk) * trb(ji,jj,jk,jpdon) / ( trb(ji,jj,jk,jpdoc) + rtrn )
                  zoxyremp  = zremikp * zoxyremc(ji,jj,jk) * trb(ji,jj,jk,jpdop) / ( trb(ji,jj,jk,jpdoc) + rtrn )

                  tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) + zolimip + zdenitrp + zoxyremp
                  tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) + zolimin + zdenitrn + zoxyremn
                  tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - denitr(ji,jj,jk) * rdenit
                  tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) - zolimic - denitr(ji,jj,jk) - zoxyremc(ji,jj,jk)
                  tra(ji,jj,jk,jpdon) = tra(ji,jj,jk,jpdon) - zolimin - zdenitrn - zoxyremn
                  tra(ji,jj,jk,jpdop) = tra(ji,jj,jk,jpdop) - zolimip - zdenitrp - zoxyremp
                  tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - zolimic * o2ut
                  tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) + zolimic + denitr(ji,jj,jk) + zoxyremc(ji,jj,jk)
                  tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * ( zolimin + zoxyremn + ( rdenit + 1.) * zdenitrn )
               END DO
            END DO
         END DO
         !
      ENDIF


      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               ! NH4 nitrification to NO3. Ceased for oxygen concentrations
               ! below 2 umol/L. Inhibited at strong light 
               ! ----------------------------------------------------------
               zonitr(ji,jj,jk)  = nitrif * xstep * trb(ji,jj,jk,jpnh4) * ( 1.- nitrfac(ji,jj,jk) )  &
               &         / ( 1.+ emoy(ji,jj,jk) ) * ( 1. + fr_i(ji,jj) * emoy(ji,jj,jk) )
               !! pjb
               !! Nitrification dependent on chemical equilibrium of NH3-NH4, which is dependent on pH.
               !zonitr(ji,jj,jk) = zonitr(ji,jj,jk) * MIN(1.2, MAX(0.0,                               &
               !&          ( 10**( (-1.0*LOG10( MAX(hi(ji,jj,jk),rtrn) ) ) - 9.3) / 10**(8.1-9.3) ) ))
               !! pjb
               zdenitnh4(ji,jj,jk) = nitrif * xstep * trb(ji,jj,jk,jpnh4) * nitrfac(ji,jj,jk)
               zdenitnh4(ji,jj,jk) = MIN(  ( trb(ji,jj,jk,jpno3) - rtrn ) / rdenita, zdenitnh4(ji,jj,jk) ) 

               ! zonitri = nitrification of NH4 --> NO3 under oxic conditions
               ! zdenitnh4 = nitrification of NH4 under anoxic conditions NH4 --> NO3
               !             ... except that this NH4 is removed from the system (NH4 --> N2)
               !             ... nitrifier-denitrification (NH4 --> NO2 --> NO  --> N2O --> N2),
               !             ... or coupled nitrification-denitrification
               !             here is where I should add anammox when I add NO2
               ! Update of the tracers trends
               ! ----------------------------
               tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - zonitr(ji,jj,jk) - zdenitnh4(ji,jj,jk)
               tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) + zonitr(ji,jj,jk) - rdenita * zdenitnh4(ji,jj,jk)
               tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) - o2nit * zonitr(ji,jj,jk)
               tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) - 2 * rno3 * zonitr(ji,jj,jk) + rno3 * ( rdenita - 1. ) * zdenitnh4(ji,jj,jk)

               IF( ln_n15 ) THEN
                  zr15_nh4 = ( (trb(ji,jj,jk,jp15nh4)+rtrn) / (trb(ji,jj,jk,jpnh4)+rtrn) )
                  zr15_no3 = ( (trb(ji,jj,jk,jp15no3)+rtrn) / (trb(ji,jj,jk,jpno3)+rtrn) )
                  tra(ji,jj,jk,jp15nh4) = tra(ji,jj,jk,jp15nh4)                                      & 
                  &                       - ( zonitr(ji,jj,jk) * ( 1.0 - e15n_nit/1000.0 ) + zdenitnh4(ji,jj,jk) ) * zr15_nh4   
                  tra(ji,jj,jk,jp15no3) = tra(ji,jj,jk,jp15no3)                                      &
                  &                       + zonitr(ji,jj,jk) * zr15_nh4 * ( 1.0 - e15n_nit/1000.0 )  &
                  &                       - rdenita * zdenitnh4(ji,jj,jk) * zr15_no3                                     
                  zonitr15(ji,jj,jk) = zonitr(ji,jj,jk) * zr15_nh4 * ( 1.0 - e15n_nit/1000.0 ) 
                  zdenitnh415no3(ji,jj,jk) = zdenitnh4(ji,jj,jk) * zr15_no3                                     
                  zdenitnh415nh4(ji,jj,jk) = zdenitnh4(ji,jj,jk) * zr15_nh4                                     
               ENDIF

            END DO
         END DO
      END DO

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem1')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

               ! Bacterial uptake of iron. No iron is available in DOC. So
               ! Bacteries are obliged to take up iron from the water. Some
               ! studies (especially at Papa) have shown this uptake to be significant
               ! ----------------------------------------------------------
               zbactfer = feratb *  rfact2 * 0.6_wp / rday * tgfunc(ji,jj,jk) * xlimbacl(ji,jj,jk)     &
                  &              * trb(ji,jj,jk,jpfer) / ( xkferb + trb(ji,jj,jk,jpfer) )    &
                  &              * zdepprod(ji,jj,jk) * zdepeff(ji,jj,jk) * zdepbac(ji,jj,jk)
               tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zbactfer*0.33
               tra(ji,jj,jk,jpsfe) = tra(ji,jj,jk,jpsfe) + zbactfer*0.25
               tra(ji,jj,jk,jpbfe) = tra(ji,jj,jk,jpbfe) + zbactfer*0.08
               zfebact(ji,jj,jk)   = zbactfer * 0.33
               blim(ji,jj,jk)      = xlimbacl(ji,jj,jk)  * zdepbac(ji,jj,jk) / 1.e-6 * zdepprod(ji,jj,jk)
            END DO
         END DO
      END DO

       IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem2')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      ! Initialization of the array which contains the labile fraction
      ! of bSi. Set to a constant in the upper ocean
      ! ---------------------------------------------------------------

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               zdep     = MAX( hmld(ji,jj), heup_01(ji,jj) )
               zsatur   = MAX( rtrn, ( sio3eq(ji,jj,jk) - trb(ji,jj,jk,jpsil) ) / ( sio3eq(ji,jj,jk) + rtrn ) )
               zsatur2  = ( 1. + tsn(ji,jj,jk,jp_tem) / 400.)**37
               znusil   = 0.225  * ( 1. + tsn(ji,jj,jk,jp_tem) / 15.) * zsatur + 0.775 * zsatur2 * zsatur**9.25
               ! Remineralization rate of BSi depedant on T and saturation
               ! ---------------------------------------------------------
               IF ( gdept_n(ji,jj,jk) > zdep ) THEN
                  zfacsib(ji,jj,jk) = zfacsib(ji,jj,jk-1) * EXP( -0.5 * ( xsiremlab - xsirem )  &
                  &                   * znusil * e3t_n(ji,jj,jk) / wsbio4(ji,jj,jk) )
                  zfacsi(ji,jj,jk)  = zfacsib(ji,jj,jk) / ( 1.0 + zfacsib(ji,jj,jk) )
                  zfacsib(ji,jj,jk) = zfacsib(ji,jj,jk) * EXP( -0.5 * ( xsiremlab - xsirem )    &
                  &                   * znusil * e3t_n(ji,jj,jk) / wsbio4(ji,jj,jk) )
               ENDIF
               zsiremin = ( xsiremlab * zfacsi(ji,jj,jk) + xsirem * ( 1. - zfacsi(ji,jj,jk) ) ) * xstep * znusil
               zosil    = zsiremin * trb(ji,jj,jk,jpgsi)
               !
               tra(ji,jj,jk,jpgsi) = tra(ji,jj,jk,jpgsi) - zosil
               tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) + zosil
            END DO
         END DO
      END DO

      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('rem3')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
       ENDIF

      IF( knt == nrdttrc ) THEN
          zrfact2 = 1.e3 * rfact2r
          ALLOCATE( zw3d(jpi,jpj,jpk) )
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "REMIN" ) )  THEN
              zw3d(:,:,:) = zolimi(:,:,:) * tmask(:,:,:) * zfact !  Remineralisation rate
              CALL iom_put( "REMIN"  , zw3d )
          ENDIF
          IF( iom_use( "REMIN_15DOC" ) )  THEN
              zw3d(:,:,:) = zolimi15(:,:,:) * tmask(:,:,:) * zfact !  Remineralisation rate
              CALL iom_put( "REMIN_15DOC"  , zw3d ) 
          ENDIF
          IF( iom_use( "REMINC" ) )  THEN
              zw3d(:,:,:) = zremindic(:,:,:) * tmask(:,:,:) * zfact !  Remineralisation rate
              CALL iom_put( "REMINC"  , zw3d )
          ENDIF
          IF( iom_use( "REMINC_13DOC" ) )  THEN
              zw3d(:,:,:) = zremindic13(:,:,:) * tmask(:,:,:) * zfact !  Remineralisation rate
              CALL iom_put( "REMINC_13DOC"  , zw3d )
          ENDIF
          IF( iom_use( "NITR" ) )  THEN
              zw3d(:,:,:) = zonitr(:,:,:) * rno3 * tmask(:,:,:) * zfact !
              CALL iom_put( "NITR"  , zw3d )
          ENDIF
          IF( iom_use( "NITR_15NH4" ) )  THEN
              zw3d(:,:,:) = zonitr15(:,:,:) * rno3 * tmask(:,:,:) * zfact !
              CALL iom_put( "NITR_15NH4"  , zw3d )
          ENDIF
          IF( iom_use( "DENIT" ) )  THEN
              zw3d(:,:,:) = denitr(:,:,:) * rdenit * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "DENIT"  , zw3d )
          ENDIF
          IF( iom_use( "DENIT_15NO3" ) )  THEN
              zw3d(:,:,:) = denitr15(:,:,:) * rdenit * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "DENIT_15NO3"  , zw3d )
          ENDIF
          IF( iom_use( "DENpNH4" ) )  THEN
              zw3d(:,:,:) = denitr(:,:,:) * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "DENpNH4"  , zw3d )
          ENDIF
          IF( iom_use( "DENpNH4_15DOC" ) )  THEN
              zw3d(:,:,:) = denitr15nh4(:,:,:) * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "DENpNH4_15DOC"  , zw3d )
          ENDIF
          IF( iom_use( "DENrNO3" ) )  THEN
              zw3d(:,:,:) = zdenitnh4(:,:,:) * rdenita * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "DENrNO3"  , zw3d )
          ENDIF
          IF( iom_use( "DENrNO3_15NO3" ) )  THEN
              zw3d(:,:,:) = zdenitnh415no3(:,:,:) * rdenita * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "DENrNO3_15NO3"  , zw3d )
          ENDIF
          IF( iom_use( "DENrNH4" ) )  THEN
              zw3d(:,:,:) = zdenitnh4(:,:,:) * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "DENrNH4"  , zw3d )
          ENDIF
          IF( iom_use( "DENrNH4_15NH4" ) )  THEN
              zw3d(:,:,:) = zdenitnh415nh4(:,:,:) * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "DENrNH4_15NH4"  , zw3d )
          ENDIF
          IF( iom_use( "ANOXREM" ) )  THEN
              zw3d(:,:,:) = zoxyremc(:,:,:) * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "ANOXREM"  , zw3d )
          ENDIF
          IF( iom_use( "ANOXREM_15DOC" ) )  THEN
              zw3d(:,:,:) = zoxyremc15(:,:,:) * rno3 * tmask(:,:,:) * zfact ! Denitrification
              CALL iom_put( "ANOXREM_15DOC"  , zw3d )
          ENDIF
          IF( iom_use( "BACT" ) )  THEN
               zw3d(:,:,:) = zdepbac(:,:,:) * 1.E6 * tmask(:,:,:)  ! Bacterial biomass
               CALL iom_put( "BACT", zw3d )
          ENDIF
          IF( iom_use( "FEBACT" ) )  THEN
               zw3d(:,:,:) = zfebact(:,:,:) * 1E9 * tmask(:,:,:) * zrfact2   ! Bacterial iron consumption
               CALL iom_put( "FEBACT" , zw3d )
          ENDIF
          !
          DEALLOCATE( zw3d )
       ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_rem')
      !
   END SUBROUTINE p4z_rem


   SUBROUTINE p4z_rem_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_rem_init  ***
      !!
      !! ** Purpose :   Initialization of remineralization parameters
      !!
      !! ** Method  :   Read the nampisrem namelist and check the parameters
      !!      called at the first timestep
      !!
      !! ** input   :   Namelist nampisrem
      !!
      !!----------------------------------------------------------------------
      NAMELIST/nampisrem/ xremik, nitrif, xsirem, xsiremlab, xsilab, feratb, xkferb, & 
         &                xremikc, xremikn, xremikp, e15n_den, e15n_nit
      INTEGER :: ios                 ! Local integer output status for namelist read
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_rem_init : Initialization of remineralization parameters'
         WRITE(numout,*) '~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampisrem in reference namelist : Pisces remineralization
      READ  ( numnatp_ref, nampisrem, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'nampisrem in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampisrem in configuration namelist : Pisces remineralization
      READ  ( numnatp_cfg, nampisrem, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'nampisrem in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, nampisrem )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist parameters for remineralization, nampisrem'
         IF( ln_p4z ) THEN
            WRITE(numout,*) '      remineralization rate of DOC              xremik    =', xremik
         ELSE
            WRITE(numout,*) '      remineralization rate of DOC              xremikc   =', xremikc
            WRITE(numout,*) '      remineralization rate of DON              xremikn   =', xremikn
            WRITE(numout,*) '      remineralization rate of DOP              xremikp   =', xremikp
         ENDIF
         WRITE(numout,*) '      remineralization rate of Si               xsirem    =', xsirem
         WRITE(numout,*) '      fast remineralization rate of Si          xsiremlab =', xsiremlab
         WRITE(numout,*) '      fraction of labile biogenic silica        xsilab    =', xsilab
         WRITE(numout,*) '      NH4 nitrification rate                    nitrif    =', nitrif
         WRITE(numout,*) '      Bacterial Fe/C ratio                      feratb    =', feratb
         WRITE(numout,*) '      N15 denitrification fractionation         e15n_den  =', e15n_den
         WRITE(numout,*) '      N15 nitrification fractionation           e15n_nit  =', e15n_nit
      ENDIF
      !
      denitr(:,:,:) = 0._wp
      zonitr(:,:,:) = 0._wp
      zdenitnh4(:,:,:) = 0._wp
      zoxyremc(:,:,:) = 0._wp
      denitr15(:,:,:) = 0._wp
      denitr15nh4(:,:,:) = 0._wp
      zonitr15(:,:,:) = 0._wp
      zdenitnh415no3(:,:,:) = 0._wp
      zdenitnh415nh4(:,:,:) = 0._wp
      zoxyremc15(:,:,:) = 0._wp
      !
   END SUBROUTINE p4z_rem_init


   INTEGER FUNCTION p4z_rem_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_rem_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( denitr(jpi,jpj,jpk), zonitr(jpi,jpj,jpk), zdenitnh4(jpi,jpj,jpk), zoxyremc(jpi,jpj,jpk),   &
      &         denitr15(jpi,jpj,jpk), zonitr15(jpi,jpj,jpk), zdenitnh415no3(jpi,jpj,jpk),                 &
      &         zdenitnh415nh4(jpi,jpj,jpk), zoxyremc15(jpi,jpj,jpk), denitr15nh4(jpi,jpj,jpk), STAT=p4z_rem_alloc )
      !
      IF( p4z_rem_alloc /= 0 )   CALL ctl_stop( 'STOP', 'p4z_rem_alloc: failed to allocate arrays' )
      !
   END FUNCTION p4z_rem_alloc

   !!======================================================================
END MODULE p4zrem
