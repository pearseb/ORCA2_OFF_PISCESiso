MODULE p4zprod
   !!======================================================================
   !!                         ***  MODULE p4zprod  ***
   !! TOP :  Growth Rate of the two phytoplanktons groups 
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!             3.4  !  2011-05  (O. Aumont, C. Ethe) New parameterization of light limitation
   !!----------------------------------------------------------------------
   !!   p4z_prod       : Compute the growth Rate of the two phytoplanktons groups
   !!   p4z_prod_init  : Initialization of the parameters for growth
   !!   p4z_prod_alloc : Allocate variables for growth
   !!----------------------------------------------------------------------
   USE oce_trc         ! shared variables between ocean and passive tracers
   USE trc             ! passive tracers common variables 
   USE sms_pisces      ! PISCES Source Minus Sink variables
   USE p4zlim          ! Co-limitations of differents nutrients
   USE p4zche
   USE prtctl_trc      ! print control for debugging
   USE iom             ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   p4z_prod         ! called in p4zbio.F90
   PUBLIC   p4z_prod_init    ! called in trcsms_pisces.F90
   PUBLIC   p4z_prod_alloc

   REAL(wp), PUBLIC ::   pislopen     !:
   REAL(wp), PUBLIC ::   pisloped     !:
   REAL(wp), PUBLIC ::   xadap        !:
   REAL(wp), PUBLIC ::   excretn      !:
   REAL(wp), PUBLIC ::   excretd      !:
   REAL(wp), PUBLIC ::   bresp        !:
   REAL(wp), PUBLIC ::   chlcnm       !:
   REAL(wp), PUBLIC ::   chlcdm       !:
   REAL(wp), PUBLIC ::   chlcmin      !:
   REAL(wp), PUBLIC ::   fecnm        !:
   REAL(wp), PUBLIC ::   fecdm        !:
   REAL(wp), PUBLIC ::   grosip       !:
   REAL(wp), PUBLIC ::   e15n_nprod   !:
   REAL(wp), PUBLIC ::   e15n_rprod   !:
   REAL(wp), PUBLIC ::   e13c_min     !:
   REAL(wp), PUBLIC ::   e13c_max     !:
   REAL(wp), PUBLIC ::   xxprod       !: pjb - factor increase in NPP

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   quotan   !: proxy of N quota in Nanophyto
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   quotad   !: proxy of N quota in diatomee
   
   REAL(wp) ::   r1_rday    ! 1 / rday
   REAL(wp) ::   texcretn   ! 1 - excretn 
   REAL(wp) ::   texcretd   ! 1 - excretd        

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zprod.F90 10401 2018-12-17 15:57:34Z cetlod $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_prod( kt , knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod  ***
      !!
      !! ** Purpose :   Compute the phytoplankton production depending on
      !!              light, temperature and nutrient availability
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt, knt   !
      !
      INTEGER  ::   ji, jj, jk
      REAL(wp) ::   zsilfac, znanotot, zdiattot, zconctemp, zconctemp2
      REAL(wp) ::   zratio, zmax, zsilim, ztn, zadap, zlim, zsilfac2, zsiborn
      REAL(wp) ::   zprod, zprochln, zprochld
      REAL(wp) ::   zmaxday, zdocprod, zpislopen, zpisloped
      REAL(wp) ::   zmxltst, zmxlday
      REAL(wp) ::   zrum, zcodel, zargu, zval, zfeup, chlcnm_n, chlcdm_n
      REAL(wp) ::   zfact
      REAL(wp) ::   zu_15, zun_15, zr15_new, zr15_reg
      REAL(wp) ::   zr13_dic, zr13, zr13_2
      REAL(wp) ::   ztc, zft, zrhop, zfco3, zbot, zdic, zph, zalka, zalk, zah2
      CHARACTER (len=25) :: charout
      REAL(wp), ALLOCATABLE, DIMENSION(:,:) :: zw2d
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) :: zw3d
      REAL(wp), DIMENSION(jpi,jpj    ) :: zstrn, zmixnano, zmixdiat
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprmaxn,zprmaxd
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpislopeadn, zpislopeadd, zysopt  
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprdia, zprbio, zprdch, zprnch   
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprorcan, zprorcad, zprofed, zprofen
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zprorcan13, zprorcad13
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpronewn, zpronewd
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zproregn, zproregd
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpronewn15, zpronewd15
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zproregn15, zproregd15
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zmxl_fac, zmxl_chl
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: zpligprod1, zpligprod2
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ze15nprod1, ze15nprod2
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: za_g, za_dic, zh2co3, z_co3
      REAL(wp), DIMENSION(jpi,jpj,jpk) :: ze13cprod1, ze13cprod2
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_prod')
      !
      !  Allocate temporary workspace
      !
      zprorcan(:,:,:) = 0._wp ; zprorcad(:,:,:) = 0._wp ; zprofed (:,:,:) = 0._wp
      zprofen (:,:,:) = 0._wp ; zysopt  (:,:,:) = 0._wp
      zpronewn(:,:,:) = 0._wp ; zpronewd(:,:,:) = 0._wp ; zprdia  (:,:,:) = 0._wp
      zproregn(:,:,:) = 0._wp ; zproregd(:,:,:) = 0._wp 
      zprbio  (:,:,:) = 0._wp ; zprdch  (:,:,:) = 0._wp ; zprnch  (:,:,:) = 0._wp 
      zmxl_fac(:,:,:) = 0._wp ; zmxl_chl(:,:,:) = 0._wp 
      IF( ln_ligand ) THEN
           zpligprod1(:,:,:)  = 0.     ;      zpligprod2(:,:,:) = 0.
      ENDIF
      IF( ln_n15 ) THEN
           ze15nprod1(:,:,:)  = 0.     ;      ze15nprod2(:,:,:) = 0.
           zpronewn15(:,:,:) = 0._wp ; zpronewd15(:,:,:) = 0._wp 
           zproregn15(:,:,:) = 0._wp ; zproregd15(:,:,:) = 0._wp 
      ENDIF
      IF( ln_c13 ) THEN
           ze13cprod1(:,:,:)  = 0.     ;      ze13cprod2(:,:,:) = 0.
           zprorcan13(:,:,:)  = 0.     ;      zprorcad13(:,:,:) = 0.
      ENDIF
      ! Computation of the optimal production (function of temperature - Eppley curve)
      zprmaxn(:,:,:) = 0.8_wp * r1_rday * tgfunc(:,:,:) * xxprod
      zprmaxd(:,:,:) = zprmaxn(:,:,:)

      ! compute the day length depending on latitude and the day
      zrum = REAL( nday_year - 80, wp ) / REAL( nyear_len(1), wp )
      zcodel = ASIN(  SIN( zrum * rpi * 2._wp ) * SIN( rad * 23.5_wp )  )

      ! day length in hours
      zstrn(:,:) = 0.
      DO jj = 1, jpj
         DO ji = 1, jpi
            zargu = TAN( zcodel ) * TAN( gphit(ji,jj) * rad )
            zargu = MAX( -1., MIN(  1., zargu ) )
            zstrn(ji,jj) = MAX( 0.0, 24. - 2. * ACOS( zargu ) / rad / 15. )
         END DO
      END DO

      ! Impact of the day duration and light intermittency on phytoplankton growth
      DO jk = 1, jpkm1
         DO jj = 1 ,jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN ! if PAR over 24h cycle > 0.001
                  zval = MAX( 1., zstrn(ji,jj) )      ! at least 1 hour of photosynthesis
                  IF( gdept_n(ji,jj,jk) <= hmld(ji,jj) ) THEN ! if depth <= mixed layer depth
                     zval = zval * MIN(1., heup_01(ji,jj) / ( hmld(ji,jj) + rtrn ))
                     ! zval (temporary variable) returns effective hours of day
                     ! by reducing the number of hours of light if the euphotic
                     ! zone depth is less than the mixed layer, as deeper mixed
                     ! layers will cause phytoplankton to see less light
                  ENDIF
                  zmxl_chl(ji,jj,jk) = zval / 24.
                  zmxl_fac(ji,jj,jk) = 1.5 * zval / ( 12. + zval )
                  ! zmxl_fac [0,1] is the limitation of growth by light, and
                  ! includes: day length (function of latitude and time of year) 
                  !           depth of euphotic zone / mixed layer depth
                  ! it solves for [0,1] range using Monod equation above
               ENDIF
            END DO
         END DO
      END DO

      ! nanophytoplankton
      zprbio(:,:,:) = zprmaxn(:,:,:) * zmxl_fac(:,:,:) ! max growth (/s) * light lim [0,1]
      ! diatoms
      zprdia(:,:,:) = zprmaxd(:,:,:) * zmxl_fac(:,:,:)

      ! Maximum light intensity
      WHERE( zstrn(:,:) < 1.e0 ) zstrn(:,:) = 24.

      ! Computation of the P-I slope for nanos and diatoms
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN  ! if PAR over 24h cycle > 0.001
                  ztn         = MAX( 0., tsn(ji,jj,jk,jp_tem) - 15. )
                  zadap       = xadap * ztn / ( 2.+ ztn )  ! this is zero because xadap is zero in namelist
                  zconctemp   = MAX( 0.e0 , trb(ji,jj,jk,jpdia) - xsizedia ) ! get concentration of diatoms above minimum amount
                  zconctemp2  = trb(ji,jj,jk,jpdia) - zconctemp ! returns minimum diatom concentration (xsizedia above) (?) 
                  ! nanophytoplankton
                  zpislopeadn(ji,jj,jk) = pislopen * ( 1.+ zadap  * EXP( -0.25 * enano(ji,jj,jk) ) )  &
                  &                   * trb(ji,jj,jk,jpnch) /( trb(ji,jj,jk,jpphy) * 12. + rtrn) 
                                                                         ! the 12 is to convert mol C to g C
                  ! diatoms
                  zpislopeadd(ji,jj,jk) = (pislopen * zconctemp2 + pisloped * zconctemp) / ( trb(ji,jj,jk,jpdia) + rtrn )   &
                  &                   * trb(ji,jj,jk,jpdch) /( trb(ji,jj,jk,jpdia) * 12. + rtrn)
                  ! In these equations, the initial slope of the PI curve (2.0)
                  ! is multiplied by the Chla:C ratio. I guess if there is more
                  ! chlorophyll-a relative to C in the cell, then the slope
                  ! should be steeper?
               ENDIF
            END DO
         END DO
      END DO

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                   ! Computation of production function for Carbon
                   !  ---------------------------------------------
                   ! adding a basic respiration rate to both groups of phyto
                   ! that is scaled by light limitation, as more respiration will accompany 
                   ! production under high light conditions
                   zpislopen = zpislopeadn(ji,jj,jk) / ( ( r1_rday + bresp * r1_rday ) &
                   &            * zmxl_fac(ji,jj,jk) * rday + rtrn)
                   zpisloped = zpislopeadd(ji,jj,jk) / ( ( r1_rday + bresp * r1_rday ) &
                   &            * zmxl_fac(ji,jj,jk) * rday + rtrn)
                   ! GROWTH RATE (/s) after light lim + respiration + PI slope (chl-a) + PAR
                   zprbio(ji,jj,jk) = zprbio(ji,jj,jk) * ( 1.- EXP( -zpislopen * enano(ji,jj,jk) )  ) ! eq 2a (GMD paper)
                   zprdia(ji,jj,jk) = zprdia(ji,jj,jk) * ( 1.- EXP( -zpisloped * ediat(ji,jj,jk) )  )
                   !  Computation of production function for Chlorophyll
                   !--------------------------------------------------
                   zpislopen = zpislopeadn(ji,jj,jk) / ( zprmaxn(ji,jj,jk) * zmxl_chl(ji,jj,jk) * rday + rtrn )
                   zpisloped = zpislopeadd(ji,jj,jk) / ( zprmaxd(ji,jj,jk) * zmxl_chl(ji,jj,jk) * rday + rtrn )
                   zprnch(ji,jj,jk) = zprmaxn(ji,jj,jk) * ( 1.- EXP( -zpislopen * enanom(ji,jj,jk) ) )
                   zprdch(ji,jj,jk) = zprmaxd(ji,jj,jk) * ( 1.- EXP( -zpisloped * ediatm(ji,jj,jk) ) )
               ENDIF
            END DO
         END DO
      END DO

      !  Computation of a proxy of the N/C ratio
      !  ---------------------------------------
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
                ! nutrient limited production / light limited production
                !  value > 1.0 means that light is more limiting
                !  value < 1.0 means that nutrients (P, N) are more limiting
                zval = MIN( xnanopo4(ji,jj,jk), ( xnanonh4(ji,jj,jk) + xnanono3(ji,jj,jk) ) )   &
                &      * zprmaxn(ji,jj,jk) / ( zprbio(ji,jj,jk) + rtrn )
                quotan(ji,jj,jk) = MIN( 1., 0.2 + 0.8 * zval ) ! < 1 when nutrient limitation is stronger than light
                zval = MIN( xdiatpo4(ji,jj,jk), ( xdiatnh4(ji,jj,jk) + xdiatno3(ji,jj,jk) ) )   &
                &      * zprmaxd(ji,jj,jk) / ( zprdia(ji,jj,jk) + rtrn )
                quotad(ji,jj,jk) = MIN( 1., 0.2 + 0.8 * zval )
                ! lower values of N:C when nutrients are more limiting than light
            END DO
         END DO
      END DO


      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi

                IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN ! when light available
                   !    Si/C of diatoms
                   !    ------------------------
                   !    Si/C increases with iron stress and silicate availability
                   !    Si/C is arbitrariliy increased for very high Si concentrations
                   !    to mimic the very high ratios observed in the Southern Ocean (silpot2)
                  zlim  = trb(ji,jj,jk,jpsil) / ( trb(ji,jj,jk,jpsil) + xksi1 ) ! Monod lim function for Si
                  zsilim = MIN( zprdia(ji,jj,jk) / ( zprmaxd(ji,jj,jk) + rtrn ), xlimsi(ji,jj,jk) )
                  zsilfac = 4.4 * EXP( -4.23 * zsilim ) * MAX( 0.e0, MIN( 1., 2.2 * ( zlim - 0.5 ) )  ) + 1.e0
                  ! this value (zsilfac) varies between 1 and 5.4.
                  ! If zlim (si limitation) is less than 0.5, meaning that Si is
                  ! relatively limiting to growth, then zsilfac==1.0
                  ! If zlim is nearer 1.0 and therefore not limiting to growth,
                  ! then the value is dependent on zsilim, which measures
                  ! the limitation of either light or other nutrients (N,P,Fe).
                  ! If light/nutrient limitation is strong, then the value of
                  ! zsilfac is closer to 5.4.
                  zsiborn = trb(ji,jj,jk,jpsil) * trb(ji,jj,jk,jpsil) * trb(ji,jj,jk,jpsil) ! si**3
                  IF (gphit(ji,jj) < -30 ) THEN ! if latitude less than 30S
                    zsilfac2 = 1. + 2. * zsiborn / ( zsiborn + xksi2**3 ) !varies from 1 (si=0) to ~3 (si>40)
                  ELSE
                    zsilfac2 = 1. +      zsiborn / ( zsiborn + xksi2**3 ) !varies from 1 (si=0) to ~2 (si>40)
                  ENDIF
                  ! Si:C ratio = mean Si:C ratio[0.159] * SiLim[0,1] * LightNutrientFactor[1,5.4] * SiExcessFactor[1,2/3]
                  zysopt(ji,jj,jk) = grosip * zlim * zsilfac * zsilfac2
                  ! higher Si:C ratios when light and nutrients (N,P,Fe) are
                  ! strongly limiting and when Si concentrations are very high
              ENDIF
            END DO
         END DO
      END DO

      !  Mixed-layer effect on production 
      !  Sea-ice effect on production

      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi  
               zprbio(ji,jj,jk) = zprbio(ji,jj,jk) * ( 1. - fr_i(ji,jj) )
               zprdia(ji,jj,jk) = zprdia(ji,jj,jk) * ( 1. - fr_i(ji,jj) )
            END DO
         END DO
      END DO


      IF ( ln_c13 ) THEN

!      DO jm = 1,10
         DO jk = 1,jpkm1
            DO jj = 1,jpj
               DO ji = 1,jpi

                  ! DUMMY VARIABLES FOR DIC, H+, AND BORATE
                  zbot  = borat(ji,jj,jk)
                  zrhop = rhop(ji,jj,jk) / 1000. + rtrn
                  zdic  = trb(ji,jj,jk,jpdic) / zrhop
                  zph   = MAX( hi(ji,jj,jk), 1.e-10 ) / zrhop
                  zalka = trb(ji,jj,jk,jptal) / zrhop

                  ! CALCULATE [ALK]([CO3--], [HCO3-])
                  zalk  = zalka - (  akw3(ji,jj,jk) / zph - zph + zbot / ( 1.+ zph / akb3(ji,jj,jk) )  )

                  ! CALCULATE [H+] AND [H2CO3]
                  zah2   = SQRT(  (zdic-zalk)**2 + 4.* ( zalk * ak23(ji,jj,jk)   &
                  &                               / ak13(ji,jj,jk) ) * ( 2.*zdic - zalk )  )
                  zah2   = 0.5 * ak13(ji,jj,jk) / zalk * ( ( zdic - zalk ) + zah2 )
                  zh2co3(ji,jj,jk) = ( 2.* zdic - zalk ) / ( 2.+ ak13(ji,jj,jk) / zah2) * zrhop
                  z_co3(ji,jj,jk) = zalk / ( 2. + zah2 / ak23(ji,jj,jk) ) * zrhop
                  hi(ji,jj,jk)   = zah2 * zrhop

               ENDDO
            ENDDO
         ENDDO
!      ENDDO
           
      DO jk = 1,jpkm1
         DO jj = 1,jpj
            DO ji = 1,jpi

               ! Compute fractionation factors for C13 from Zhang et al. 1995

               zfco3 = MAX(0.05,( (z_co3(ji,jj,jk)+rtrn) / (trb(ji,jj,jk,jpdic)+rtrn) ))
               zfco3 = MIN(0.2 , zfco3)
               ztc = MIN( 35., tsn(ji,jj,jk,jp_tem) )
               zft = MIN(25.,ztc)
               zft = MAX( 5.,zft)
               za_g(ji,jj,jk) = 1.   + ( -0.0049 * zft - 1.31 ) / 1000.
               za_dic(ji,jj,jk) = 1. + ( 0.014 * zft * zfco3 - 0.105 * zft + 10.53 ) / 1000.


            ENDDO
         ENDDO
      ENDDO

      ENDIF


      ! Computation of the various production terms 
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  !  production terms for nanophyto. (C)
                  ! zprbio = growth rate (max * daylength/ice/light/respiration/PIslope factors)
                  ! xlimphy = nutrient limitation factor (P,N,Fe) (Liebig law of the minimum)
                  ! trb(ji,jj,jk,jpphy) = concentration of phytoplankton
                  ! rfact = number of seconds per timestep (21600 secs, or 6 hours)
                  zprorcan(ji,jj,jk) = zprbio(ji,jj,jk)  * xlimphy(ji,jj,jk) * trb(ji,jj,jk,jpphy) * rfact2
                  ! zprorcan = production of C per timestep (growth rate * molC/L nanophytoplankton)
                  ! xnanono3 = nutrient limitation due to NO3
                  ! xnanonh4 = nutrient limitation due to NH4 
                  zpronewn(ji,jj,jk)  = zprorcan(ji,jj,jk)* xnanono3(ji,jj,jk) / ( xnanono3(ji,jj,jk) + xnanonh4(ji,jj,jk) + rtrn )
                  ! Ahh... so new production here is amount of total production
                  ! (zprorcan) that is due solely to NO3, not NH4
                  !
                  !
                  zratio = trb(ji,jj,jk,jpnfe) / ( trb(ji,jj,jk,jpphy) * fecnm + rtrn )
                  zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) ) 
                  zprofen(ji,jj,jk) = fecnm * zprmaxn(ji,jj,jk) * ( 1.0 - fr_i(ji,jj) )  &
                  &             * ( 4. - 4.5 * xlimnfe(ji,jj,jk) / ( xlimnfe(ji,jj,jk) + 0.5 ) )    &
                  &             * biron(ji,jj,jk) / ( biron(ji,jj,jk) + concnfe(ji,jj,jk) )  &
                  &             * zmax * trb(ji,jj,jk,jpphy) * rfact2
                  !  production terms for diatoms (C)
                  zprorcad(ji,jj,jk) = zprdia(ji,jj,jk) * xlimdia(ji,jj,jk) * trb(ji,jj,jk,jpdia) * rfact2
                  zpronewd(ji,jj,jk) = zprorcad(ji,jj,jk) * xdiatno3(ji,jj,jk) / ( xdiatno3(ji,jj,jk) + xdiatnh4(ji,jj,jk) + rtrn )
                  !
                  zratio = trb(ji,jj,jk,jpdfe) / ( trb(ji,jj,jk,jpdia) * fecdm + rtrn )
                  zmax   = MAX( 0., ( 1. - zratio ) / ABS( 1.05 - zratio ) ) 
                  zprofed(ji,jj,jk) = fecdm * zprmaxd(ji,jj,jk) * ( 1.0 - fr_i(ji,jj) )  &
                  &             * ( 4. - 4.5 * xlimdfe(ji,jj,jk) / ( xlimdfe(ji,jj,jk) + 0.5 ) )    &
                  &             * biron(ji,jj,jk) / ( biron(ji,jj,jk) + concdfe(ji,jj,jk) )  &
                  &             * zmax * trb(ji,jj,jk,jpdia) * rfact2
               ENDIF
            END DO
         END DO
      END DO

      ! Computation of the chlorophyll production terms
      DO jk = 1, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
               IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                  !  production terms for nanophyto. ( chlorophyll )
                  znanotot = enanom(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zprod    = rday * zprorcan(ji,jj,jk) * zprnch(ji,jj,jk) * xlimphy(ji,jj,jk)
                  zprochln = chlcmin * 12. * zprorcan (ji,jj,jk)
                  chlcnm_n   = MIN ( chlcnm, ( chlcnm / (1. - 1.14 / 43.4 *tsn(ji,jj,jk,jp_tem))) * (1. - 1.14 / 43.4 * 20.))
                  zprochln = zprochln + (chlcnm_n-chlcmin) * 12. * zprod / &
                                        & (  zpislopeadn(ji,jj,jk) * znanotot +rtrn)
                  !  production terms for diatoms ( chlorophyll )
                  zdiattot = ediatm(ji,jj,jk) / ( zmxl_chl(ji,jj,jk) + rtrn )
                  zprod    = rday * zprorcad(ji,jj,jk) * zprdch(ji,jj,jk) * xlimdia(ji,jj,jk)
                  zprochld = chlcmin * 12. * zprorcad(ji,jj,jk)
                  chlcdm_n   = MIN ( chlcdm, ( chlcdm / (1. - 1.14 / 43.4 * tsn(ji,jj,jk,jp_tem))) * (1. - 1.14 / 43.4 * 20.))
                  zprochld = zprochld + (chlcdm_n-chlcmin) * 12. * zprod / &
                                        & ( zpislopeadd(ji,jj,jk) * zdiattot +rtrn )
                  !   Update the arrays TRA which contain the Chla sources and sinks
                  tra(ji,jj,jk,jpnch) = tra(ji,jj,jk,jpnch) + zprochln * texcretn
                  tra(ji,jj,jk,jpdch) = tra(ji,jj,jk,jpdch) + zprochld * texcretd
               ENDIF
            END DO
         END DO
      END DO

      !   Update the arrays TRA which contain the biological sources and sinks
      DO jk = 1, jpkm1
         DO jj = 1, jpj
           DO ji =1 ,jpi
              IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                 zproregn(ji,jj,jk) = zprorcan(ji,jj,jk) - zpronewn(ji,jj,jk)
                 zproregd(ji,jj,jk) = zprorcad(ji,jj,jk) - zpronewd(ji,jj,jk)
                 zdocprod = excretd * zprorcad(ji,jj,jk) + excretn * zprorcan(ji,jj,jk)
                 ! zprorcan & zprorcad are total DIC taken up by nano & diat
                 ! zproregn & zproregd are total NH4 taken up by nano & diat
                 ! zpronewn & zpronewd are total NO3 taken up by nano & diat
                 ! zprorcan * texcretn = total carbon uptake after excretion [DIC] --> [PHY]
                 ! zprorcad * texcretd = total carbon uptake after excretion [DIC] --> [PHY2]
                 ! zdocprod = total carbon excreted as DOC [DIC] --> [PHY+PHY2] --> [DOC]
                 tra(ji,jj,jk,jppo4) = tra(ji,jj,jk,jppo4) - zprorcan(ji,jj,jk) - zprorcad(ji,jj,jk)
                 tra(ji,jj,jk,jpno3) = tra(ji,jj,jk,jpno3) - zpronewn(ji,jj,jk) - zpronewd(ji,jj,jk)
                 tra(ji,jj,jk,jpnh4) = tra(ji,jj,jk,jpnh4) - zproregn(ji,jj,jk) - zproregd(ji,jj,jk)
                 tra(ji,jj,jk,jpphy) = tra(ji,jj,jk,jpphy) + zprorcan(ji,jj,jk) * texcretn
                 tra(ji,jj,jk,jpnfe) = tra(ji,jj,jk,jpnfe) + zprofen(ji,jj,jk) * texcretn
                 tra(ji,jj,jk,jpdia) = tra(ji,jj,jk,jpdia) + zprorcad(ji,jj,jk) * texcretd
                 tra(ji,jj,jk,jpdfe) = tra(ji,jj,jk,jpdfe) + zprofed(ji,jj,jk) * texcretd
                 tra(ji,jj,jk,jpdsi) = tra(ji,jj,jk,jpdsi) + zprorcad(ji,jj,jk) * zysopt(ji,jj,jk) * texcretd
                 tra(ji,jj,jk,jpdoc) = tra(ji,jj,jk,jpdoc) + zdocprod
                 tra(ji,jj,jk,jpoxy) = tra(ji,jj,jk,jpoxy) + o2ut * ( zproregn(ji,jj,jk) + zproregd(ji,jj,jk)) &
                 &                   + ( o2ut + o2nit ) * ( zpronewn(ji,jj,jk) + zpronewd(ji,jj,jk) )
                 !
                 zfeup = texcretn * zprofen(ji,jj,jk) + texcretd * zprofed(ji,jj,jk)
                 tra(ji,jj,jk,jpfer) = tra(ji,jj,jk,jpfer) - zfeup
                 tra(ji,jj,jk,jpsil) = tra(ji,jj,jk,jpsil) - texcretd * zprorcad(ji,jj,jk) * zysopt(ji,jj,jk)
                 tra(ji,jj,jk,jpdic) = tra(ji,jj,jk,jpdic) - zprorcan(ji,jj,jk) - zprorcad(ji,jj,jk)
                 tra(ji,jj,jk,jptal) = tra(ji,jj,jk,jptal) + rno3 * ( zpronewn(ji,jj,jk) + zpronewd(ji,jj,jk) ) &
                 &                                         - rno3 * ( zproregn(ji,jj,jk) + zproregd(ji,jj,jk) )

                 IF ( ln_n15 ) THEN
                    ! First, calculate utilisation factors for new and regenerated production.
                    zu_15 = MIN(1.0, MAX(0.0, 1.0 - (zpronewn(ji,jj,jk) + zpronewd(ji,jj,jk) + rtrn)   &
                    &                             / (trb(ji,jj,jk,jpno3) + rtrn) ) )
                    zun_15 = MIN(1.0, MAX(0.0, 1.0 - (zproregn(ji,jj,jk) + zproregd(ji,jj,jk) + rtrn)   &
                    &                              / (trb(ji,jj,jk,jpnh4) + rtrn) ) )
                    ! save the fractionation factor due to new and regenerated production
                    ze15nprod1(ji,jj,jk) = e15n_nprod * zu_15
                    ze15nprod2(ji,jj,jk) = e15n_rprod * zun_15
                    
                    !! pjb
                    !IF ( gphit(ji,jj) > -20 .and. gphit(ji,jj) < 20 ) THEN ! select tropical latitudes 
                    !  ze15nprod1(ji,jj,jk) = max(e15n_nprod,5.0) * zu_15
                    !  ze15nprod2(ji,jj,jk) = max(e15n_rprod,5.0) * zun_15
                    !ENDIF
                    !! pjb

                    ! Second, apply utilisation factors to fractionation factors
                    zr15_new = ( 1.0 - ze15nprod1(ji,jj,jk) / 1000.0 )                                 &
                    &           * ( (trb(ji,jj,jk,jp15no3)+rtrn) / (trb(ji,jj,jk,jpno3)+rtrn) )
                    zr15_reg = ( 1.0 - ze15nprod2(ji,jj,jk) / 1000.0 )                                 &
                    &           * ( (trb(ji,jj,jk,jp15nh4)+rtrn) / (trb(ji,jj,jk,jpnh4)+rtrn) )
                    ! Third, apply the fractionation factors to the state variables
                    tra(ji,jj,jk,jp15phy) = tra(ji,jj,jk,jp15phy) + zr15_new * zpronewn(ji,jj,jk) * texcretn &
                    &                                             + zr15_reg * zproregn(ji,jj,jk) * texcretn
                    tra(ji,jj,jk,jp15dia) = tra(ji,jj,jk,jp15dia) + zr15_new * zpronewd(ji,jj,jk) * texcretd &
                    &                                             + zr15_reg * zproregd(ji,jj,jk) * texcretd
                    tra(ji,jj,jk,jp15doc) = tra(ji,jj,jk,jp15doc) + zr15_new * zpronewn(ji,jj,jk) * excretn &
                    &                                             + zr15_reg * zproregn(ji,jj,jk) * excretn &
                    &                                             + zr15_new * zpronewd(ji,jj,jk) * excretd &
                    &                                             + zr15_reg * zproregd(ji,jj,jk) * excretd 
                    tra(ji,jj,jk,jp15no3) = tra(ji,jj,jk,jp15no3) - zr15_new * ( zpronewn(ji,jj,jk) + zpronewd(ji,jj,jk) )
                    tra(ji,jj,jk,jp15nh4) = tra(ji,jj,jk,jp15nh4) - zr15_reg * ( zproregn(ji,jj,jk) + zproregd(ji,jj,jk) )
                    zpronewn15(ji,jj,jk) = zr15_new * zpronewn(ji,jj,jk)
                    zpronewd15(ji,jj,jk) = zr15_new * zpronewd(ji,jj,jk) 
                    zproregn15(ji,jj,jk) = zr15_reg * zproregn(ji,jj,jk)
                    zproregd15(ji,jj,jk) = zr15_reg * zproregd(ji,jj,jk) 
                 ENDIF
                 IF ( ln_c13 ) THEN
                    zr13_dic = ( (trb(ji,jj,jk,jp13dic)+rtrn) / (trb(ji,jj,jk,jpdic)+rtrn) )
                    
                    ze13cprod1(ji,jj,jk) = max(e13c_min, min(e13c_max, (                                    &
                    &             ( (86400 * zprorcan(ji,jj,jk)) / (trb(ji,jj,jk,jpphy) + rtrn) /           &
                    &             ( rfact2 * (zh2co3(ji,jj,jk)/1025*1e9 + rtrn)) - 0.371 ) / (-0.015) )))
                    ze13cprod2(ji,jj,jk) = max(e13c_min, min(e13c_max, (                                    &
                    &             ( (86400 * zprorcad(ji,jj,jk)) / (trb(ji,jj,jk,jpdia) + rtrn) /           &
                    &             ( rfact2 * (zh2co3(ji,jj,jk)/1025*1e9 + rtrn)) - 0.371 ) / (-0.015) )))
                    !!pjb
                    !IF ( gphit(ji,jj) > -20 .and. gphit(ji,jj) < 20 ) THEN ! select tropical latitudes 
                    !  ze13cprod1(ji,jj,jk) = 25.0 
                    !  ze13cprod2(ji,jj,jk) = 25.0 
                    !ENDIF
                    !!pjb
                
                    zr13   = ( 1.0 - ze13cprod1(ji,jj,jk)/1000.0 ) * zr13_dic 
                    zr13_2 = ( 1.0 - ze13cprod2(ji,jj,jk)/1000.0 ) * zr13_dic 

                    tra(ji,jj,jk,jp13phy) = tra(ji,jj,jk,jp13phy) + zprorcan(ji,jj,jk) * texcretn * zr13
                    tra(ji,jj,jk,jp13dia) = tra(ji,jj,jk,jp13dia) + zprorcad(ji,jj,jk) * texcretd * zr13_2
                    tra(ji,jj,jk,jp13doc) = tra(ji,jj,jk,jp13doc) + zprorcan(ji,jj,jk) * excretn * zr13      &
                    &                                             + zprorcad(ji,jj,jk) * excretd * zr13_2
                    tra(ji,jj,jk,jp13dic) = tra(ji,jj,jk,jp13dic) - zprorcan(ji,jj,jk) * zr13                &
                    &                                             - zprorcad(ji,jj,jk) * zr13_2 
                    zprorcan13(ji,jj,jk) = zprorcan(ji,jj,jk) * zr13
                    zprorcad13(ji,jj,jk) = zprorcad(ji,jj,jk) * zr13_2
                 ENDIF

              ENDIF
           END DO
        END DO
     END DO
     !
     IF( ln_ligand ) THEN
         DO jk = 1, jpkm1
            DO jj = 1, jpj
              DO ji =1 ,jpi
                 IF( etot_ndcy(ji,jj,jk) > 1.E-3 ) THEN
                    zdocprod = excretd * zprorcad(ji,jj,jk) + excretn * zprorcan(ji,jj,jk)
                    zfeup    = texcretn * zprofen(ji,jj,jk) + texcretd * zprofed(ji,jj,jk)
                    tra(ji,jj,jk,jplgw) = tra(ji,jj,jk,jplgw) + zdocprod * ldocp    &
                    &       - zfeup * plig(ji,jj,jk) / ( rtrn + plig(ji,jj,jk) + 2.E3 * (1.0 - plig(ji,jj,jk) ) ) 
                    zpligprod1(ji,jj,jk) = zdocprod * ldocp
                    zpligprod2(ji,jj,jk) = zfeup * plig(ji,jj,jk) / ( rtrn + plig(ji,jj,jk) + 2.E3 * (1.0 - plig(ji,jj,jk) ) )
                 ENDIF
              END DO
           END DO
        END DO
     ENDIF


    ! Total primary production per year
    IF( iom_use( "tintpp" ) .OR. ( ln_check_mass .AND. kt == nitend .AND. knt == nrdttrc )  )  &
         & tpp = glob_sum( 'p4zprod', ( zprorcan(:,:,:) + zprorcad(:,:,:) ) * cvol(:,:,:) )

    IF( lk_iomput ) THEN
       IF( knt == nrdttrc ) THEN
          ALLOCATE( zw2d(jpi,jpj), zw3d(jpi,jpj,jpk) )
          zfact = 1.e+3 * rfact2r  !  conversion from mol/l/kt to  mol/m3/s
          !
          IF( iom_use( "PPPHYN" ) .OR. iom_use( "PPPHYD" ) )  THEN
              zw3d(:,:,:) = zprorcan(:,:,:) * zfact * tmask(:,:,:)  ! primary production by nanophyto
              CALL iom_put( "PPPHYN"  , zw3d )
              !
              zw3d(:,:,:) = zprorcad(:,:,:) * zfact * tmask(:,:,:)  ! primary production by diatomes
              CALL iom_put( "PPPHYD"  , zw3d )
          ENDIF
          IF( iom_use( "PPPHYN_13DIC" ) .OR. iom_use( "PPPHYD_13DIC" ) )  THEN
              zw3d(:,:,:) = zprorcan13(:,:,:) * zfact * tmask(:,:,:)  ! primary production by nanophyto
              CALL iom_put( "PPPHYN_13DIC"  , zw3d )
              !
              zw3d(:,:,:) = zprorcad13(:,:,:) * zfact * tmask(:,:,:)  ! primary production by diatomes
              CALL iom_put( "PPPHYD_13DIC"  , zw3d )
          ENDIF
          IF( iom_use( "PPNEWN" ) .OR. iom_use( "PPNEWD" ) )  THEN
              zw3d(:,:,:) = zpronewn(:,:,:) * zfact * tmask(:,:,:)  ! new primary production by nanophyto
              CALL iom_put( "PPNEWN"  , zw3d )
              !
              zw3d(:,:,:) = zpronewd(:,:,:) * zfact * tmask(:,:,:)  ! new primary production by diatomes
              CALL iom_put( "PPNEWD"  , zw3d )
          ENDIF
          IF( iom_use( "PPNEWN_15NO3" ) .OR. iom_use( "PPNEWD_15NO3" ) )  THEN
              zw3d(:,:,:) = zpronewn15(:,:,:) * zfact * tmask(:,:,:)  ! new primary production by nanophyto
              CALL iom_put( "PPNEWN_15NO3"  , zw3d )
              !
              zw3d(:,:,:) = zpronewd15(:,:,:) * zfact * tmask(:,:,:)  ! new primary production by diatomes
              CALL iom_put( "PPNEWD_15NO3"  , zw3d )
          ENDIF
          IF( iom_use( "PPREGN" ) .OR. iom_use( "PPREGD" ) )  THEN
              zw3d(:,:,:) = zproregn(:,:,:) * zfact * tmask(:,:,:)  ! regenerated primary production by nanophyto
              CALL iom_put( "PPREGN"  , zw3d )
              !
              zw3d(:,:,:) = zproregd(:,:,:) * zfact * tmask(:,:,:)  ! regenerated primary production by diatomes
              CALL iom_put( "PPREGD"  , zw3d )
          ENDIF
          IF( iom_use( "PPREGN_15NH4" ) .OR. iom_use( "PPREGD_15NH4" ) )  THEN
              zw3d(:,:,:) = zproregn15(:,:,:) *zfact * tmask(:,:,:)  ! regenerated primary production by nanophyto
              CALL iom_put( "PPREGN_15NH4"  , zw3d )
              !
              zw3d(:,:,:) = zproregd15(:,:,:) * zfact * tmask(:,:,:)  ! regenerated primary production by diatomes
              CALL iom_put( "PPREGD_15NH4"  , zw3d )
          ENDIF
          IF( iom_use( "PBSi" ) )  THEN
              zw3d(:,:,:) = zprorcad(:,:,:) * zfact * tmask(:,:,:) * zysopt(:,:,:) ! biogenic silica production
              CALL iom_put( "PBSi"  , zw3d )
          ENDIF
          IF( iom_use( "PFeN" ) .OR. iom_use( "PFeD" ) )  THEN
              zw3d(:,:,:) = zprofen(:,:,:) * zfact * tmask(:,:,:)  ! biogenic iron production by nanophyto
              CALL iom_put( "PFeN"  , zw3d )
              !
              zw3d(:,:,:) = zprofed(:,:,:) * zfact * tmask(:,:,:)  ! biogenic iron production by  diatomes
              CALL iom_put( "PFeD"  , zw3d )
          ENDIF
          IF( iom_use( "LPRODP" ) )  THEN
              zw3d(:,:,:) = zpligprod1(:,:,:) * 1e9 * zfact * tmask(:,:,:)
              CALL iom_put( "LPRODP"  , zw3d )
          ENDIF
          IF( iom_use( "LDETP" ) )  THEN
              zw3d(:,:,:) = zpligprod2(:,:,:) * 1e9 * zfact * tmask(:,:,:)
              CALL iom_put( "LDETP"  , zw3d )
          ENDIF
          IF( iom_use( "E15Nnew" ) )  THEN
              zw3d(:,:,:) = ze15nprod1(:,:,:) * tmask(:,:,:)
              CALL iom_put( "E15Nnew"  , zw3d )
          ENDIF
          IF( iom_use( "E15Nreg" ) )  THEN
              zw3d(:,:,:) = ze15nprod2(:,:,:) * tmask(:,:,:)
              CALL iom_put( "E15Nreg"  , zw3d )
          ENDIF
          IF( iom_use( "E13Cphy" ) )  THEN
              zw3d(:,:,:) = ze13cprod1(:,:,:) * tmask(:,:,:)
              CALL iom_put( "E13Cphy"  , zw3d )
          ENDIF
          IF( iom_use( "E13Cdia" ) )  THEN
              zw3d(:,:,:) = ze13cprod2(:,:,:) * tmask(:,:,:)
              CALL iom_put( "E13Cdia"  , zw3d )
          ENDIF
          IF( iom_use( "Mumax" ) )  THEN
              zw3d(:,:,:) = zprmaxn(:,:,:) * tmask(:,:,:)   ! Maximum growth rate
              CALL iom_put( "Mumax"  , zw3d )
          ENDIF
          IF( iom_use( "MuN" ) .OR. iom_use( "MuD" ) )  THEN
              zw3d(:,:,:) = zprbio(:,:,:) * xlimphy(:,:,:) * tmask(:,:,:)  ! Realized growth rate for nanophyto
              CALL iom_put( "MuN"  , zw3d )
              !
              zw3d(:,:,:) =  zprdia(:,:,:) * xlimdia(:,:,:) * tmask(:,:,:)  ! Realized growth rate for diatoms
              CALL iom_put( "MuD"  , zw3d )
          ENDIF
          IF( iom_use( "LNlight" ) .OR. iom_use( "LDlight" ) )  THEN
              zw3d(:,:,:) = zprbio (:,:,:) / (zprmaxn(:,:,:) + rtrn) * tmask(:,:,:) ! light limitation term
              CALL iom_put( "LNlight"  , zw3d )
              !
              zw3d(:,:,:) = zprdia (:,:,:) / (zprmaxd(:,:,:) + rtrn) * tmask(:,:,:)  ! light limitation term
              CALL iom_put( "LDlight"  , zw3d )
          ENDIF
          IF( iom_use( "TPP" ) )  THEN
              zw3d(:,:,:) = ( zprorcan(:,:,:) + zprorcad(:,:,:) ) * zfact * tmask(:,:,:)  ! total primary production
              CALL iom_put( "TPP"  , zw3d )
          ENDIF
          IF( iom_use( "TPNEW" ) )  THEN
              zw3d(:,:,:) = ( zpronewn(:,:,:) + zpronewd(:,:,:) ) * zfact * tmask(:,:,:)  ! total new production
              CALL iom_put( "TPNEW"  , zw3d )
          ENDIF
          IF( iom_use( "TPBFE" ) )  THEN
              zw3d(:,:,:) = ( zprofen(:,:,:) + zprofed(:,:,:) ) * zfact * tmask(:,:,:)  ! total biogenic iron production
              CALL iom_put( "TPBFE"  , zw3d )
          ENDIF
          IF( iom_use( "INTPPPHYN" ) .OR. iom_use( "INTPPPHYD" ) ) THEN  
             zw2d(:,:) = 0.
             DO jk = 1, jpkm1
               zw2d(:,:) = zw2d(:,:) + zprorcan(:,:,jk) * e3t_n(:,:,jk) * zfact * tmask(:,:,jk)  ! vert. integrated  primary produc. by nano
             ENDDO
             CALL iom_put( "INTPPPHYN" , zw2d )
             !
             zw2d(:,:) = 0.
             DO jk = 1, jpkm1
                zw2d(:,:) = zw2d(:,:) + zprorcad(:,:,jk) * e3t_n(:,:,jk) * zfact * tmask(:,:,jk) ! vert. integrated  primary produc. by diatom
             ENDDO
             CALL iom_put( "INTPPPHYD" , zw2d )
          ENDIF
          IF( iom_use( "INTPP" ) ) THEN   
             zw2d(:,:) = 0.
             DO jk = 1, jpkm1
                zw2d(:,:) = zw2d(:,:) + ( zprorcan(:,:,jk) + zprorcad(:,:,jk) ) * e3t_n(:,:,jk) * zfact * tmask(:,:,jk) ! vert. integrated pp
             ENDDO
             CALL iom_put( "INTPP" , zw2d )
          ENDIF
          IF( iom_use( "INTPNEW" ) ) THEN    
             zw2d(:,:) = 0.
             DO jk = 1, jpkm1
                zw2d(:,:) = zw2d(:,:) + ( zpronewn(:,:,jk) + zpronewd(:,:,jk) ) * e3t_n(:,:,jk) * zfact * tmask(:,:,jk)  ! vert. integrated new prod
             ENDDO
             CALL iom_put( "INTPNEW" , zw2d )
          ENDIF
          IF( iom_use( "INTPBFE" ) ) THEN           !   total biogenic iron production  ( vertically integrated )
             zw2d(:,:) = 0.
             DO jk = 1, jpkm1
                zw2d(:,:) = zw2d(:,:) + ( zprofen(:,:,jk) + zprofed(:,:,jk) ) * e3t_n(:,:,jk) * zfact * tmask(:,:,jk) ! vert integr. bfe prod
             ENDDO
            CALL iom_put( "INTPBFE" , zw2d )
          ENDIF
          IF( iom_use( "INTPBSI" ) ) THEN           !   total biogenic silica production  ( vertically integrated )
             zw2d(:,:) = 0.
             DO jk = 1, jpkm1
                zw2d(:,:) = zw2d(:,:) + zprorcad(:,:,jk) * zysopt(:,:,jk) * e3t_n(:,:,jk) * zfact * tmask(:,:,jk)  ! vert integr. bsi prod
             ENDDO
             CALL iom_put( "INTPBSI" , zw2d )
          ENDIF
          IF( iom_use( "tintpp" ) )  CALL iom_put( "tintpp" , tpp * zfact )  !  global total integrated primary production molC/s
          !
          DEALLOCATE( zw2d, zw3d )
       ENDIF
     ENDIF

     IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('prod')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
     ENDIF
      !
      IF( ln_timing )  CALL timing_stop('p4z_prod')
      !
   END SUBROUTINE p4z_prod


   SUBROUTINE p4z_prod_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE p4z_prod_init  ***
      !!
      !! ** Purpose :   Initialization of phytoplankton production parameters
      !!
      !! ** Method  :   Read the nampisprod namelist and check the parameters
      !!      called at the first timestep (nittrc000)
      !!
      !! ** input   :   Namelist nampisprod
      !!----------------------------------------------------------------------
      INTEGER ::   ios   ! Local integer
      !
      NAMELIST/namp4zprod/ pislopen, pisloped, xadap, bresp, excretn, excretd,  &
         &                 chlcnm, chlcdm, chlcmin, fecnm, fecdm, grosip,       &
         &                 e15n_nprod, e15n_rprod, e13c_min, e13c_max, xxprod
      !!----------------------------------------------------------------------
      !
      IF(lwp) THEN                         ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'p4z_prod_init : phytoplankton growth'
         WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      !
      REWIND( numnatp_ref )              ! Namelist nampisprod in reference namelist : Pisces phytoplankton production
      READ  ( numnatp_ref, namp4zprod, IOSTAT = ios, ERR = 901)
901   IF( ios /= 0 )   CALL ctl_nam ( ios , 'namp4zprod in reference namelist', lwp )
      REWIND( numnatp_cfg )              ! Namelist nampisprod in configuration namelist : Pisces phytoplankton production
      READ  ( numnatp_cfg, namp4zprod, IOSTAT = ios, ERR = 902 )
902   IF( ios >  0 )   CALL ctl_nam ( ios , 'namp4zprod in configuration namelist', lwp )
      IF(lwm) WRITE( numonp, namp4zprod )

      IF(lwp) THEN                         ! control print
         WRITE(numout,*) '   Namelist : namp4zprod'
         WRITE(numout,*) '      mean Si/C ratio                           grosip       =', grosip
         WRITE(numout,*) '      P-I slope                                 pislopen     =', pislopen
         WRITE(numout,*) '      Acclimation factor to low light           xadap        =', xadap
         WRITE(numout,*) '      excretion ratio of nanophytoplankton      excretn      =', excretn
         WRITE(numout,*) '      excretion ratio of diatoms                excretd      =', excretd
         WRITE(numout,*) '      basal respiration in phytoplankton        bresp        =', bresp
         WRITE(numout,*) '      Maximum Chl/C in phytoplankton            chlcmin      =', chlcmin
         WRITE(numout,*) '      P-I slope  for diatoms                    pisloped     =', pisloped
         WRITE(numout,*) '      Minimum Chl/C in nanophytoplankton        chlcnm       =', chlcnm
         WRITE(numout,*) '      Minimum Chl/C in diatoms                  chlcdm       =', chlcdm
         WRITE(numout,*) '      Maximum Fe/C in nanophytoplankton         fecnm        =', fecnm
         WRITE(numout,*) '      Minimum Fe/C in diatoms                   fecdm        =', fecdm
         WRITE(numout,*) '      N15 fractionation by new production       e15n_nprod   =', e15n_nprod
         WRITE(numout,*) '      N15 fractionation by reg production       e15n_rprod   =', e15n_rprod
         WRITE(numout,*) '      C13 assimilation fractionation min        e13c_min     =', e13c_min
         WRITE(numout,*) '      C13 assimilation fractionation max        e13c_max     =', e13c_max
         WRITE(numout,*) '      pjb - factor increase in NPP              xxprod       =', xxprod
      ENDIF
      !
      r1_rday   = 1._wp / rday 
      texcretn  = 1._wp - excretn
      texcretd  = 1._wp - excretd
      tpp       = 0._wp
      !
   END SUBROUTINE p4z_prod_init


   INTEGER FUNCTION p4z_prod_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_prod_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( quotan(jpi,jpj,jpk), quotad(jpi,jpj,jpk), STAT = p4z_prod_alloc )
      !
      IF( p4z_prod_alloc /= 0 ) CALL ctl_warn('p4z_prod_alloc : failed to allocate arrays.')
      !
   END FUNCTION p4z_prod_alloc

   !!======================================================================
END MODULE p4zprod
