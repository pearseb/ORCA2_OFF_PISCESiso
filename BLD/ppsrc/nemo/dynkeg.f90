# 1 "dynkeg.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "dynkeg.F90"
MODULE dynkeg
   !!======================================================================
   !!                       ***  MODULE  dynkeg  ***
   !! Ocean dynamics:  kinetic energy gradient trend
   !!======================================================================
   !! History :  1.0  !  1987-09  (P. Andrich, M.-A. Foujols)  Original code
   !!            7.0  !  1997-05  (G. Madec)  Split dynber into dynkeg and dynhpg
   !!  NEMO      1.0  !  2002-07  (G. Madec)  F90: Free form and module
   !!            3.6  !  2015-05  (N. Ducousso, G. Madec)  add Hollingsworth scheme as an option 
   !!----------------------------------------------------------------------
   
   !!----------------------------------------------------------------------
   !!   dyn_keg      : update the momentum trend with the horizontal tke
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE trd_oce         ! trends: ocean variables
   USE trddyn          ! trend manager: dynamics
   !
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! MPP library
   USE prtctl          ! Print control
   USE timing          ! Timing
   USE bdy_oce         ! ocean open boundary conditions

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_keg    ! routine called by step module
   
   INTEGER, PARAMETER, PUBLIC  ::   nkeg_C2  = 0   !: 2nd order centered scheme (standard scheme)
   INTEGER, PARAMETER, PUBLIC  ::   nkeg_HW  = 1   !: Hollingsworth et al., QJRMS, 1983
   !
   REAL(wp) ::   r1_48 = 1._wp / 48._wp   !: =1/(4*2*6)
   
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
# 39 "dynkeg.F90" 2
   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: dynkeg.F90 10997 2019-05-17 16:16:13Z girrmann $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_keg( kt, kscheme )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_keg  ***
      !!
      !! ** Purpose :   Compute the now momentum trend due to the horizontal
      !!      gradient of the horizontal kinetic energy and add it to the 
      !!      general momentum trend.
      !!
      !! ** Method  : * kscheme = nkeg_C2 : 2nd order centered scheme that 
      !!      conserve kinetic energy. Compute the now horizontal kinetic energy 
      !!         zhke = 1/2 [ mi-1( un^2 ) + mj-1( vn^2 ) ]
      !!              * kscheme = nkeg_HW : Hollingsworth correction following
      !!      Arakawa (2001). The now horizontal kinetic energy is given by:
      !!         zhke = 1/6 [ mi-1(  2 * un^2 + ((un(j+1)+un(j-1))/2)^2  )
      !!                    + mj-1(  2 * vn^2 + ((vn(i+1)+vn(i-1))/2)^2  ) ]
      !!      
      !!      Take its horizontal gradient and add it to the general momentum
      !!      trend (ua,va).
      !!         ua = ua - 1/e1u di[ zhke ]
      !!         va = va - 1/e2v dj[ zhke ]
      !!
      !! ** Action : - Update the (ua, va) with the hor. ke gradient trend
      !!             - send this trends to trd_dyn (l_trddyn=T) for post-processing
      !!
      !! ** References : Arakawa, A., International Geophysics 2001.
      !!                 Hollingsworth et al., Quart. J. Roy. Meteor. Soc., 1983.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt        ! ocean time-step index
      INTEGER, INTENT( in ) ::   kscheme   ! =0/1   type of KEG scheme 
      !
      INTEGER  ::   ji, jj, jk, jb           ! dummy loop indices
      INTEGER  ::   ifu, ifv, igrd, ib_bdy   ! local integers
      REAL(wp) ::   zu, zv                   ! local scalars
      REAL(wp), DIMENSION(jpi,jpj,jpk)        ::   zhke
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:) ::   ztrdu, ztrdv 
      REAL(wp)  :: zweightu, zweightv
      !!----------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('dyn_keg')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_keg : kinetic energy gradient trend, scheme number=', kscheme
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF

      IF( l_trddyn ) THEN           ! Save the input trends
         ALLOCATE( ztrdu(jpi,jpj,jpk) , ztrdv(jpi,jpj,jpk) )
         ztrdu(:,:,:) = ua(:,:,:) 
         ztrdv(:,:,:) = va(:,:,:) 
      ENDIF
      
      zhke(:,:,jpk) = 0._wp

      SELECT CASE ( kscheme )             !== Horizontal kinetic energy at T-point  ==!
      !
      CASE ( nkeg_C2 )                          !--  Standard scheme  --!
         DO jk = 1, jpkm1
            DO jj = 2, jpj
               DO ji = 2, jpi   ! vector opt.
                  zu =    un(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)   &
                     &  + un(ji  ,jj  ,jk) * un(ji  ,jj  ,jk)
                  zv =    vn(ji  ,jj-1,jk) * vn(ji  ,jj-1,jk)   &
                     &  + vn(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk)
                  zhke(ji,jj,jk) = 0.25_wp * ( zv + zu )
               END DO  
            END DO
         END DO
         !
         IF (ln_bdy) THEN
            ! Maria Luneva & Fred Wobus: July-2016
            ! compensate for lack of turbulent kinetic energy on liquid bdy points
            DO ib_bdy = 1, nb_bdy
               IF( cn_dyn3d(ib_bdy) /= 'none' ) THEN
                  igrd = 1           ! compensating null velocity on the bdy
                  DO jb = 1, idx_bdy(ib_bdy)%nblenrim(igrd)
                     ji   = idx_bdy(ib_bdy)%nbi(jb,igrd)   ! maximum extent : from 2 to jpi-1
                     jj   = idx_bdy(ib_bdy)%nbj(jb,igrd)   ! maximum extent : from 2 to jpj-1
                     DO jk = 1, jpkm1
                        zhke(ji,jj,jk) = 0._wp
                        zweightu = umask(ji-1,jj  ,jk) + umask(ji,jj,jk)
                        zweightv = vmask(ji  ,jj-1,jk) + vmask(ji,jj,jk)
                        zu = un(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)  +  un(ji  ,jj  ,jk) * un(ji  ,jj  ,jk)
                        zv = vn(ji  ,jj-1,jk) * vn(ji  ,jj-1,jk)  +  vn(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk)
                        IF( zweightu > 0._wp )   zhke(ji,jj,jk) =  zhke(ji,jj,jk) + zu / (2._wp * zweightu) 
                        IF( zweightv > 0._wp )   zhke(ji,jj,jk) =  zhke(ji,jj,jk) + zv / (2._wp * zweightv) 
                     END DO
                  END DO
               END IF
               CALL lbc_bdy_lnk( 'dynkeg', zhke, 'T', 1., ib_bdy )   ! send 2 and recv jpi, jpj used in the computation of the speed tendencies
            END DO
         END IF
         !
      CASE ( nkeg_HW )                          !--  Hollingsworth scheme  --!
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1       
               DO ji = 2, jpim1   ! vector opt.
                  zu = 8._wp * ( un(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)    &
                     &         + un(ji  ,jj  ,jk) * un(ji  ,jj  ,jk) )  &
                     &   +     ( un(ji-1,jj-1,jk) + un(ji-1,jj+1,jk) ) * ( un(ji-1,jj-1,jk) + un(ji-1,jj+1,jk) )   &
                     &   +     ( un(ji  ,jj-1,jk) + un(ji  ,jj+1,jk) ) * ( un(ji  ,jj-1,jk) + un(ji  ,jj+1,jk) )
                     !
                  zv = 8._wp * ( vn(ji  ,jj-1,jk) * vn(ji  ,jj-1,jk)    &
                     &         + vn(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk) )  &
                     &  +      ( vn(ji-1,jj-1,jk) + vn(ji+1,jj-1,jk) ) * ( vn(ji-1,jj-1,jk) + vn(ji+1,jj-1,jk) )   &
                     &  +      ( vn(ji-1,jj  ,jk) + vn(ji+1,jj  ,jk) ) * ( vn(ji-1,jj  ,jk) + vn(ji+1,jj  ,jk) )
                  zhke(ji,jj,jk) = r1_48 * ( zv + zu )
               END DO  
            END DO
         END DO
         IF (ln_bdy) THEN
            ! Maria Luneva & Fred Wobus: July-2016
            ! compensate for lack of turbulent kinetic energy on liquid bdy points
            DO ib_bdy = 1, nb_bdy
               IF( cn_dyn3d(ib_bdy) /= 'none' ) THEN
                  igrd = 1           ! compensation null velocity on land at the bdy
                  DO jb = 1, idx_bdy(ib_bdy)%nblenrim(igrd)
                     ji   = idx_bdy(ib_bdy)%nbi(jb,igrd)   ! maximum extent : from 2 to jpi-1
                     jj   = idx_bdy(ib_bdy)%nbj(jb,igrd)   ! maximum extent : from 2 to jpj-1
                     DO jk = 1, jpkm1
                        zhke(ji,jj,jk) = 0._wp
                        zweightu = 8._wp * ( umask(ji-1,jj  ,jk) + umask(ji  ,jj  ,jk) ) &
                             &   + 2._wp * ( umask(ji-1,jj-1,jk) + umask(ji-1,jj+1,jk) + umask(ji  ,jj-1,jk) + umask(ji  ,jj+1,jk) ) 
                        zweightv = 8._wp * ( vmask(ji  ,jj-1,jk) + vmask(ji  ,jj-1,jk) ) &
                             &   + 2._wp * ( vmask(ji-1,jj-1,jk) + vmask(ji+1,jj-1,jk) + vmask(ji-1,jj  ,jk) + vmask(ji+1,jj  ,jk) )
                        zu = 8._wp * ( un(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)    &
                           &         + un(ji  ,jj  ,jk) * un(ji  ,jj  ,jk) )  &
                           &   +     ( un(ji-1,jj-1,jk) + un(ji-1,jj+1,jk) ) * ( un(ji-1,jj-1,jk) + un(ji-1,jj+1,jk) )   &
                           &   +     ( un(ji  ,jj-1,jk) + un(ji  ,jj+1,jk) ) * ( un(ji  ,jj-1,jk) + un(ji  ,jj+1,jk) )
                        zv = 8._wp * ( vn(ji  ,jj-1,jk) * vn(ji  ,jj-1,jk)    &
                           &         + vn(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk) )  &
                           &  +      ( vn(ji-1,jj-1,jk) + vn(ji+1,jj-1,jk) ) * ( vn(ji-1,jj-1,jk) + vn(ji+1,jj-1,jk) )   &
                           &  +      ( vn(ji-1,jj  ,jk) + vn(ji+1,jj  ,jk) ) * ( vn(ji-1,jj  ,jk) + vn(ji+1,jj  ,jk) )
                        IF( zweightu > 0._wp )   zhke(ji,jj,jk) =  zhke(ji,jj,jk) +  zu / ( 2._wp * zweightu )
                        IF( zweightv > 0._wp )   zhke(ji,jj,jk) =  zhke(ji,jj,jk) +  zv / ( 2._wp * zweightv )
                     END DO
                  END DO
               END IF
            END DO
         END IF
         CALL lbc_lnk( 'dynkeg', zhke, 'T', 1. )
         !
      END SELECT 
      !
      DO jk = 1, jpkm1                    !==  grad( KE ) added to the general momentum trends  ==!
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               ua(ji,jj,jk) = ua(ji,jj,jk) - ( zhke(ji+1,jj  ,jk) - zhke(ji,jj,jk) ) / e1u(ji,jj)
               va(ji,jj,jk) = va(ji,jj,jk) - ( zhke(ji  ,jj+1,jk) - zhke(ji,jj,jk) ) / e2v(ji,jj)
            END DO 
         END DO
      END DO
      !
      IF( l_trddyn ) THEN                 ! save the Kinetic Energy trends for diagnostic
         ztrdu(:,:,:) = ua(:,:,:) - ztrdu(:,:,:)
         ztrdv(:,:,:) = va(:,:,:) - ztrdv(:,:,:)
         CALL trd_dyn( ztrdu, ztrdv, jpdyn_keg, kt )
         DEALLOCATE( ztrdu , ztrdv )
      ENDIF
      !
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=ua, clinfo1=' keg  - Ua: ', mask1=umask,   &
         &                       tab3d_2=va, clinfo2=       ' Va: ', mask2=vmask, clinfo3='dyn' )
      !
      IF( ln_timing )   CALL timing_stop('dyn_keg')
      !
   END SUBROUTINE dyn_keg

   !!======================================================================
END MODULE dynkeg
