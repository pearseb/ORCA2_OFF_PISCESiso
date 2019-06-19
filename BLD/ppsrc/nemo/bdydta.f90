# 1 "bdydta.F90"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "bdydta.F90"
MODULE bdydta
   !!======================================================================
   !!                       ***  MODULE bdydta  ***
   !! Open boundary data : read the data for the unstructured open boundaries.
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!             -   !  2007-01  (D. Storkey) Update to use IOM module
   !!             -   !  2007-07  (D. Storkey) add bdy_dta_fla
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.3  !  2010-09  (E.O'Dea) modifications for Shelf configurations 
   !!            3.3  !  2010-09  (D.Storkey) add ice boundary conditions
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!            3.6  !  2012-01  (C. Rousset) add ice boundary conditions for sea ice
   !!            4.0  !  2018     (C. Rousset) SI3 compatibility
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!    bdy_dta      : read external data along open boundaries from file
   !!    bdy_dta_init : initialise arrays etc for reading of external data
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constants
   USE sbcapr         ! atmospheric pressure forcing
   USE sbctide        ! Tidal forcing or not
   USE bdy_oce        ! ocean open boundary conditions  
   USE bdytides       ! tidal forcing at boundaries




   !
   USE lib_mpp, ONLY: ctl_stop, ctl_nam
   USE fldread        ! read input fields
   USE iom            ! IOM library
   USE in_out_manager ! I/O logical units
   USE timing         ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   bdy_dta          ! routine called by step.F90 and dynspg_ts.F90
   PUBLIC   bdy_dta_init     ! routine called by nemogcm.F90

   INTEGER, ALLOCATABLE, DIMENSION(:)   ::   nb_bdy_fld        ! Number of fields to update for each boundary set.
   INTEGER                              ::   nb_bdy_fld_sum    ! Total number of fields to update for all boundary sets.
   LOGICAL,           DIMENSION(jp_bdy) ::   ln_full_vel_array ! =T => full velocities in 3D boundary conditions
                                                               ! =F => baroclinic velocities in 3D boundary conditions
!$AGRIF_DO_NOT_TREAT
   TYPE(FLD), PUBLIC, ALLOCATABLE, DIMENSION(:), TARGET ::   bf        ! structure of input fields (file informations, fields read)
!$AGRIF_END_DO_NOT_TREAT
   TYPE(MAP_POINTER), ALLOCATABLE, DIMENSION(:) :: nbmap_ptr   ! array of pointers to nbmap







   !!----------------------------------------------------------------------
   !! NEMO/OCE 4.0 , NEMO Consortium (2018)
   !! $Id: bdydta.F90 10952 2019-05-09 10:29:57Z jamesharle $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

      SUBROUTINE bdy_dta( kt, jit, time_offset )
      !!----------------------------------------------------------------------
      !!                   ***  SUBROUTINE bdy_dta  ***
      !!                    
      !! ** Purpose :   Update external data for open boundary conditions
      !!
      !! ** Method  :   Use fldread.F90
      !!                
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in)           ::   kt           ! ocean time-step index 
      INTEGER, INTENT(in), OPTIONAL ::   jit          ! subcycle time-step index (for timesplitting option)
      INTEGER, INTENT(in), OPTIONAL ::   time_offset  ! time offset in units of timesteps. NB. if jit
      !                                               ! is present then units = subcycle timesteps.
      !                                               ! time_offset = 0 => get data at "now" time level
      !                                               ! time_offset = -1 => get data at "before" time level
      !                                               ! time_offset = +1 => get data at "after" time level
      !                                               ! etc.
      !
      INTEGER ::  jbdy, jfld, jstart, jend, ib, jl  ! dummy loop indices
      INTEGER ::  ii, ij, ik, igrd                  ! local integers
      INTEGER,          DIMENSION(jpbgrd) ::   ilen1 
      INTEGER, POINTER, DIMENSION(:)      ::   nblen, nblenrim  ! short cuts
      TYPE(OBC_DATA), POINTER             ::   dta              ! short cut
      !!---------------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('bdy_dta')
      !
      ! Initialise data arrays once for all from initial conditions where required
      !---------------------------------------------------------------------------
      IF( kt == nit000 .AND. .NOT.PRESENT(jit) ) THEN

         ! Calculate depth-mean currents
         !-----------------------------
         
         DO jbdy = 1, nb_bdy
            !
            nblen    => idx_bdy(jbdy)%nblen
            nblenrim => idx_bdy(jbdy)%nblenrim
            dta      => dta_bdy(jbdy)
            !
            IF( nn_dyn2d_dta(jbdy) == 0 ) THEN 
               ilen1(:) = nblen(:)
               IF( dta%ll_ssh ) THEN 
                  igrd = 1
                  DO ib = 1, ilen1(igrd)
                     ii = idx_bdy(jbdy)%nbi(ib,igrd)
                     ij = idx_bdy(jbdy)%nbj(ib,igrd)
                     dta_bdy(jbdy)%ssh(ib) = sshn(ii,ij) * tmask(ii,ij,1)         
                  END DO 
               ENDIF
               IF( dta%ll_u2d ) THEN 
                  igrd = 2
                  DO ib = 1, ilen1(igrd)
                     ii = idx_bdy(jbdy)%nbi(ib,igrd)
                     ij = idx_bdy(jbdy)%nbj(ib,igrd)
                     dta_bdy(jbdy)%u2d(ib) = un_b(ii,ij) * umask(ii,ij,1)         
                  END DO 
               ENDIF
               IF( dta%ll_v2d ) THEN 
                  igrd = 3
                  DO ib = 1, ilen1(igrd)
                     ii = idx_bdy(jbdy)%nbi(ib,igrd)
                     ij = idx_bdy(jbdy)%nbj(ib,igrd)
                     dta_bdy(jbdy)%v2d(ib) = vn_b(ii,ij) * vmask(ii,ij,1)         
                  END DO 
               ENDIF
            ENDIF
            !
            IF( nn_dyn3d_dta(jbdy) == 0 ) THEN 
               ilen1(:) = nblen(:)
               IF( dta%ll_u3d ) THEN 
                  igrd = 2 
                  DO ib = 1, ilen1(igrd)
                     DO ik = 1, jpkm1
                        ii = idx_bdy(jbdy)%nbi(ib,igrd)
                        ij = idx_bdy(jbdy)%nbj(ib,igrd)
                        dta_bdy(jbdy)%u3d(ib,ik) =  ( un(ii,ij,ik) - un_b(ii,ij) ) * umask(ii,ij,ik)         
                     END DO
                  END DO 
               ENDIF
               IF( dta%ll_v3d ) THEN 
                  igrd = 3 
                  DO ib = 1, ilen1(igrd)
                     DO ik = 1, jpkm1
                        ii = idx_bdy(jbdy)%nbi(ib,igrd)
                        ij = idx_bdy(jbdy)%nbj(ib,igrd)
                        dta_bdy(jbdy)%v3d(ib,ik) =  ( vn(ii,ij,ik) - vn_b(ii,ij) ) * vmask(ii,ij,ik)         
                        END DO
                  END DO 
               ENDIF
            ENDIF

            IF( nn_tra_dta(jbdy) == 0 ) THEN 
               ilen1(:) = nblen(:)
               IF( dta%ll_tem ) THEN
                  igrd = 1 
                  DO ib = 1, ilen1(igrd)
                     DO ik = 1, jpkm1
                        ii = idx_bdy(jbdy)%nbi(ib,igrd)
                        ij = idx_bdy(jbdy)%nbj(ib,igrd)
                        dta_bdy(jbdy)%tem(ib,ik) = tsn(ii,ij,ik,jp_tem) * tmask(ii,ij,ik)         
                     END DO
                  END DO 
               ENDIF
               IF( dta%ll_sal ) THEN
                  igrd = 1 
                  DO ib = 1, ilen1(igrd)
                     DO ik = 1, jpkm1
                        ii = idx_bdy(jbdy)%nbi(ib,igrd)
                        ij = idx_bdy(jbdy)%nbj(ib,igrd)
                        dta_bdy(jbdy)%sal(ib,ik) = tsn(ii,ij,ik,jp_sal) * tmask(ii,ij,ik)         
                     END DO
                  END DO 
               ENDIF
            ENDIF

# 218 "bdydta.F90"
         END DO ! jbdy
         !
      ENDIF ! kt == nit000

      ! update external data from files
      !--------------------------------
     
      jstart = 1
      DO jbdy = 1, nb_bdy   
         dta => dta_bdy(jbdy)
         IF( nn_dta(jbdy) == 1 ) THEN ! skip this bit if no external data required
      
            IF( PRESENT(jit) ) THEN
               ! Update barotropic boundary conditions only
               ! jit is optional argument for fld_read and bdytide_update
               IF( cn_dyn2d(jbdy) /= 'none' ) THEN
                  IF( nn_dyn2d_dta(jbdy) == 2 ) THEN ! tidal harmonic forcing ONLY: initialise arrays
                     IF( dta%ll_ssh ) dta%ssh(:) = 0._wp
                     IF( dta%ll_u2d ) dta%u2d(:) = 0._wp
                     IF( dta%ll_u3d ) dta%v2d(:) = 0._wp
                  ENDIF
                  IF (cn_tra(jbdy) /= 'runoff') THEN
                     IF( nn_dyn2d_dta(jbdy) == 1 .OR. nn_dyn2d_dta(jbdy) == 3 ) THEN

                        jend = jstart + dta%nread(2) - 1
                        IF( ln_full_vel_array(jbdy) ) THEN
                           CALL fld_read( kt=kt, kn_fsbc=1, sd=bf(jstart:jend), map=nbmap_ptr(jstart:jend),  &
                                     & kit=jit, kt_offset=time_offset , jpk_bdy=nb_jpk_bdy(jbdy),   &
                                     & fvl=ln_full_vel_array(jbdy)  )
                        ELSE
                           CALL fld_read( kt=kt, kn_fsbc=1, sd=bf(jstart:jend), map=nbmap_ptr(jstart:jend),  &
                                     & kit=jit, kt_offset=time_offset  )
                        ENDIF

                        ! If full velocities in boundary data then extract barotropic velocities from 3D fields
                        IF( ln_full_vel_array(jbdy) .AND.                                             &
                          &    ( nn_dyn2d_dta(jbdy) == 1 .OR. nn_dyn2d_dta(jbdy) == 3 .OR.  &
                          &      nn_dyn3d_dta(jbdy) == 1 ) )THEN

                           igrd = 2                      ! zonal velocity
                           dta%u2d(:) = 0._wp
                           DO ib = 1, idx_bdy(jbdy)%nblen(igrd)
                              ii   = idx_bdy(jbdy)%nbi(ib,igrd)
                              ij   = idx_bdy(jbdy)%nbj(ib,igrd)
                              DO ik = 1, jpkm1
                                 dta%u2d(ib) = dta%u2d(ib) &
                       &                          + e3u_n(ii,ij,ik) * umask(ii,ij,ik) * dta%u3d(ib,ik)
                              END DO
                              dta%u2d(ib) =  dta%u2d(ib) * r1_hu_n(ii,ij)
                           END DO
                           igrd = 3                      ! meridional velocity
                           dta%v2d(:) = 0._wp
                           DO ib = 1, idx_bdy(jbdy)%nblen(igrd)
                              ii   = idx_bdy(jbdy)%nbi(ib,igrd)
                              ij   = idx_bdy(jbdy)%nbj(ib,igrd)
                              DO ik = 1, jpkm1
                                 dta%v2d(ib) = dta%v2d(ib) &
                       &                       + e3v_n(ii,ij,ik) * vmask(ii,ij,ik) * dta%v3d(ib,ik)
                              END DO
                              dta%v2d(ib) =  dta%v2d(ib) * r1_hv_n(ii,ij)
                           END DO
                        ENDIF                    
                     ENDIF
                     IF( nn_dyn2d_dta(jbdy) .ge. 2 ) THEN ! update tidal harmonic forcing
                        CALL bdytide_update( kt=kt, idx=idx_bdy(jbdy), dta=dta, td=tides(jbdy),   & 
                          &                 jit=jit, time_offset=time_offset )
                     ENDIF
                  ENDIF
               ENDIF
            ELSE
               IF (cn_tra(jbdy) == 'runoff') then      ! runoff condition
                  jend = nb_bdy_fld(jbdy)
                  CALL fld_read( kt=kt, kn_fsbc=1, sd=bf(jstart:jend),  &
                               & map=nbmap_ptr(jstart:jend), kt_offset=time_offset )
                  !
                  igrd = 2                      ! zonal velocity
                  DO ib = 1, idx_bdy(jbdy)%nblen(igrd)
                     ii   = idx_bdy(jbdy)%nbi(ib,igrd)
                     ij   = idx_bdy(jbdy)%nbj(ib,igrd)
                     dta%u2d(ib) = dta%u2d(ib) / ( e2u(ii,ij) * hu_0(ii,ij) )
                  END DO
                  !
                  igrd = 3                      ! meridional velocity
                  DO ib = 1, idx_bdy(jbdy)%nblen(igrd)
                     ii   = idx_bdy(jbdy)%nbi(ib,igrd)
                     ij   = idx_bdy(jbdy)%nbj(ib,igrd)
                     dta%v2d(ib) = dta%v2d(ib) / ( e1v(ii,ij) * hv_0(ii,ij) )
                  END DO
               ELSE
                  IF( nn_dyn2d_dta(jbdy) == 2 ) THEN ! tidal harmonic forcing ONLY: initialise arrays
                     IF( dta%ll_ssh ) dta%ssh(:) = 0._wp
                     IF( dta%ll_u2d ) dta%u2d(:) = 0._wp
                     IF( dta%ll_v2d ) dta%v2d(:) = 0._wp
                  ENDIF
                  IF( dta%nread(1) .gt. 0 ) THEN ! update external data
                     jend = jstart + dta%nread(1) - 1
                     CALL fld_read( kt=kt, kn_fsbc=1, sd=bf(jstart:jend), &
                                  & map=nbmap_ptr(jstart:jend), kt_offset=time_offset, jpk_bdy=nb_jpk_bdy(jbdy),   &
                                  & fvl=ln_full_vel_array(jbdy) )
                  ENDIF
                  ! If full velocities in boundary data then split into barotropic and baroclinic data
                  IF( ln_full_vel_array(jbdy) .and.                                             &
                    & ( nn_dyn2d_dta(jbdy) == 1 .OR. nn_dyn2d_dta(jbdy) == 3 .OR. &
                    &   nn_dyn3d_dta(jbdy) == 1 ) ) THEN
                     igrd = 2                      ! zonal velocity
                     dta%u2d(:) = 0._wp
                     DO ib = 1, idx_bdy(jbdy)%nblen(igrd)
                        ii   = idx_bdy(jbdy)%nbi(ib,igrd)
                        ij   = idx_bdy(jbdy)%nbj(ib,igrd)
                        DO ik = 1, jpkm1
                           dta%u2d(ib) = dta%u2d(ib) &
                 &                       + e3u_n(ii,ij,ik) * umask(ii,ij,ik) * dta%u3d(ib,ik)
                        END DO
                        dta%u2d(ib) =  dta%u2d(ib) * r1_hu_n(ii,ij)
                        DO ik = 1, jpkm1
                           dta%u3d(ib,ik) = dta%u3d(ib,ik) - dta%u2d(ib)
                        END DO
                     END DO
                     igrd = 3                      ! meridional velocity
                     dta%v2d(:) = 0._wp
                     DO ib = 1, idx_bdy(jbdy)%nblen(igrd)
                        ii   = idx_bdy(jbdy)%nbi(ib,igrd)
                        ij   = idx_bdy(jbdy)%nbj(ib,igrd)
                        DO ik = 1, jpkm1
                           dta%v2d(ib) = dta%v2d(ib) &
                 &                       + e3v_n(ii,ij,ik) * vmask(ii,ij,ik) * dta%v3d(ib,ik)
                        END DO
                        dta%v2d(ib) =  dta%v2d(ib) * r1_hv_n(ii,ij)
                        DO ik = 1, jpkm1
                           dta%v3d(ib,ik) = dta%v3d(ib,ik) - dta%v2d(ib)
                        END DO
                     END DO
                  ENDIF

               ENDIF
# 368 "bdydta.F90"
            ENDIF
            jstart = jstart + dta%nread(1)
         ENDIF    ! nn_dta(jbdy) = 1
      END DO  ! jbdy

      IF ( ln_apr_obc ) THEN
         DO jbdy = 1, nb_bdy
            IF (cn_tra(jbdy) /= 'runoff')THEN
               igrd = 1                      ! meridional velocity
               DO ib = 1, idx_bdy(jbdy)%nblenrim(igrd)
                  ii   = idx_bdy(jbdy)%nbi(ib,igrd)
                  ij   = idx_bdy(jbdy)%nbj(ib,igrd)
                  dta_bdy(jbdy)%ssh(ib) = dta_bdy(jbdy)%ssh(ib) + ssh_ib(ii,ij)
               END DO
            ENDIF
         END DO
      ENDIF

      IF ( ln_tide ) THEN
         IF (ln_dynspg_ts) THEN      ! Fill temporary arrays with slow-varying bdy data                           
            DO jbdy = 1, nb_bdy    ! Tidal component added in ts loop
               IF ( nn_dyn2d_dta(jbdy) .ge. 2 ) THEN
                  nblen => idx_bdy(jbdy)%nblen
                  nblenrim => idx_bdy(jbdy)%nblenrim
                  IF( cn_dyn2d(jbdy) == 'frs' ) THEN; ilen1(:)=nblen(:) ; ELSE ; ilen1(:)=nblenrim(:) ; ENDIF 
                  IF ( dta_bdy(jbdy)%ll_ssh ) dta_bdy_s(jbdy)%ssh(1:ilen1(1)) = dta_bdy(jbdy)%ssh(1:ilen1(1))
                  IF ( dta_bdy(jbdy)%ll_u2d ) dta_bdy_s(jbdy)%u2d(1:ilen1(2)) = dta_bdy(jbdy)%u2d(1:ilen1(2))
                  IF ( dta_bdy(jbdy)%ll_v2d ) dta_bdy_s(jbdy)%v2d(1:ilen1(3)) = dta_bdy(jbdy)%v2d(1:ilen1(3))
               ENDIF
            END DO
         ELSE ! Add tides if not split-explicit free surface else this is done in ts loop
            !
            CALL bdy_dta_tides( kt=kt, time_offset=time_offset )
         ENDIF
      ENDIF

      !
      IF( ln_timing )   CALL timing_stop('bdy_dta')
      !
   END SUBROUTINE bdy_dta


   SUBROUTINE bdy_dta_init
      !!----------------------------------------------------------------------
      !!                   ***  SUBROUTINE bdy_dta_init  ***
      !!                    
      !! ** Purpose :   Initialise arrays for reading of external data 
      !!                for open boundary conditions
      !!
      !! ** Method  :   
      !!                
      !!----------------------------------------------------------------------
      INTEGER ::   jbdy, jfld, jstart, jend, ierror, ios     ! Local integers
      !
      CHARACTER(len=100)                     ::   cn_dir        ! Root directory for location of data files
      CHARACTER(len=100), DIMENSION(nb_bdy)  ::   cn_dir_array  ! Root directory for location of data files
      CHARACTER(len = 256)::   clname                           ! temporary file name
      LOGICAL                                ::   ln_full_vel   ! =T => full velocities in 3D boundary data
      !                                                         ! =F => baroclinic velocities in 3D boundary data
      INTEGER                                ::   ilen_global   ! Max length required for global bdy dta arrays
      INTEGER, ALLOCATABLE, DIMENSION(:)     ::   ilen1, ilen3  ! size of 1st and 3rd dimensions of local arrays
      INTEGER, ALLOCATABLE, DIMENSION(:)     ::   ibdy           ! bdy set for a particular jfld
      INTEGER, ALLOCATABLE, DIMENSION(:)     ::   igrid         ! index for grid type (1,2,3 = T,U,V)
      INTEGER, POINTER, DIMENSION(:)         ::   nblen, nblenrim  ! short cuts
      TYPE(OBC_DATA), POINTER                ::   dta           ! short cut





      TYPE(FLD_N), ALLOCATABLE, DIMENSION(:) ::   blf_i         !  array of namelist information structures
      TYPE(FLD_N) ::   bn_tem, bn_sal, bn_u3d, bn_v3d   ! 
      TYPE(FLD_N) ::   bn_ssh, bn_u2d, bn_v2d           ! informations about the fields to be read



      NAMELIST/nambdy_dta/ cn_dir, bn_tem, bn_sal, bn_u3d, bn_v3d, bn_ssh, bn_u2d, bn_v2d 



      NAMELIST/nambdy_dta/ ln_full_vel
      !!---------------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'bdy_dta_ini : initialization of data at the open boundaries'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
      IF(lwp) WRITE(numout,*) ''

      ! Set nn_dta
      DO jbdy = 1, nb_bdy
         nn_dta(jbdy) = MAX(   nn_dyn2d_dta  (jbdy)    &
            &                , nn_dyn3d_dta  (jbdy)    &
            &                , nn_tra_dta    (jbdy)    &



                              )
         IF(nn_dta(jbdy) > 1)   nn_dta(jbdy) = 1
      END DO

      ! Work out upper bound of how many fields there are to read in and allocate arrays
      ! ---------------------------------------------------------------------------
      ALLOCATE( nb_bdy_fld(nb_bdy) )
      nb_bdy_fld(:) = 0
      DO jbdy = 1, nb_bdy         
         IF( cn_dyn2d(jbdy) /= 'none' .AND. ( nn_dyn2d_dta(jbdy) == 1 .or. nn_dyn2d_dta(jbdy) == 3 ) ) THEN
            nb_bdy_fld(jbdy) = nb_bdy_fld(jbdy) + 3
         ENDIF
         IF( cn_dyn3d(jbdy) /= 'none' .AND. nn_dyn3d_dta(jbdy) == 1 ) THEN
            nb_bdy_fld(jbdy) = nb_bdy_fld(jbdy) + 2
         ENDIF
         IF( cn_tra(jbdy) /= 'none' .AND. nn_tra_dta(jbdy) == 1  ) THEN
            nb_bdy_fld(jbdy) = nb_bdy_fld(jbdy) + 2
         ENDIF





         IF(lwp) WRITE(numout,*) 'Maximum number of files to open =', nb_bdy_fld(jbdy)
      END DO            

      nb_bdy_fld_sum = SUM( nb_bdy_fld )

      ALLOCATE( bf(nb_bdy_fld_sum), STAT=ierror )
      IF( ierror > 0 ) THEN   
         CALL ctl_stop( 'bdy_dta: unable to allocate bf structure' )   ;   RETURN  
      ENDIF
      ALLOCATE( blf_i(nb_bdy_fld_sum), STAT=ierror )
      IF( ierror > 0 ) THEN   
         CALL ctl_stop( 'bdy_dta: unable to allocate blf_i structure' )   ;   RETURN  
      ENDIF
      ALLOCATE( nbmap_ptr(nb_bdy_fld_sum), STAT=ierror )
      IF( ierror > 0 ) THEN   
         CALL ctl_stop( 'bdy_dta: unable to allocate nbmap_ptr structure' )   ;   RETURN  
      ENDIF
      ALLOCATE( ilen1(nb_bdy_fld_sum), ilen3(nb_bdy_fld_sum) ) 
      ALLOCATE( ibdy(nb_bdy_fld_sum) ) 
      ALLOCATE( igrid(nb_bdy_fld_sum) ) 

      ! Read namelists
      ! --------------
      REWIND(numnam_cfg)
      jfld = 0 
      DO jbdy = 1, nb_bdy         
         IF( nn_dta(jbdy) == 1 ) THEN
            REWIND(numnam_ref)
            READ  ( numnam_ref, nambdy_dta, IOSTAT = ios, ERR = 901)
901         IF( ios /= 0 )   CALL ctl_nam ( ios , 'nambdy_dta in reference namelist', lwp )
            READ  ( numnam_cfg, nambdy_dta, IOSTAT = ios, ERR = 902 )
902         IF( ios >  0 )   CALL ctl_nam ( ios , 'nambdy_dta in configuration namelist', lwp )
            IF(lwm) WRITE( numond, nambdy_dta )

            cn_dir_array(jbdy) = cn_dir
            ln_full_vel_array(jbdy) = ln_full_vel

            nblen => idx_bdy(jbdy)%nblen
            nblenrim => idx_bdy(jbdy)%nblenrim
            dta => dta_bdy(jbdy)
            dta%nread(2) = 0

            ! Only read in necessary fields for this set.
            ! Important that barotropic variables come first.
            IF( nn_dyn2d_dta(jbdy) == 1 .or. nn_dyn2d_dta(jbdy) == 3 ) THEN 

               IF( dta%ll_ssh ) THEN 
                  if(lwp) write(numout,*) '++++++ reading in ssh field'
                  jfld = jfld + 1
                  blf_i(jfld) = bn_ssh
                  ibdy(jfld) = jbdy
                  igrid(jfld) = 1
                  ilen1(jfld) = nblen(igrid(jfld))
                  ilen3(jfld) = 1
                  dta%nread(2) = dta%nread(2) + 1
               ENDIF

               IF( dta%ll_u2d .and. .not. ln_full_vel_array(jbdy) ) THEN
                  if(lwp) write(numout,*) '++++++ reading in u2d field'
                  jfld = jfld + 1
                  blf_i(jfld) = bn_u2d
                  ibdy(jfld) = jbdy
                  igrid(jfld) = 2
                  ilen1(jfld) = nblen(igrid(jfld))
                  ilen3(jfld) = 1
                  dta%nread(2) = dta%nread(2) + 1
               ENDIF

               IF( dta%ll_v2d .and. .not. ln_full_vel_array(jbdy) ) THEN
                  if(lwp) write(numout,*) '++++++ reading in v2d field'
                  jfld = jfld + 1
                  blf_i(jfld) = bn_v2d
                  ibdy(jfld) = jbdy
                  igrid(jfld) = 3
                  ilen1(jfld) = nblen(igrid(jfld))
                  ilen3(jfld) = 1
                  dta%nread(2) = dta%nread(2) + 1
               ENDIF

            ENDIF

            ! read 3D velocities if baroclinic velocities require OR if
            ! barotropic velocities required and ln_full_vel set to .true.
            IF( nn_dyn3d_dta(jbdy) == 1 .OR. &
           &  ( ln_full_vel_array(jbdy) .AND. ( nn_dyn2d_dta(jbdy) == 1 .OR. nn_dyn2d_dta(jbdy) == 3 ) ) ) THEN

               IF( dta%ll_u3d .OR. ( ln_full_vel_array(jbdy) .and. dta%ll_u2d ) ) THEN
                  if(lwp) write(numout,*) '++++++ reading in u3d field'
                  jfld = jfld + 1
                  blf_i(jfld) = bn_u3d
                  ibdy(jfld) = jbdy
                  igrid(jfld) = 2
                  ilen1(jfld) = nblen(igrid(jfld))
                  ilen3(jfld) = jpk
                  IF( ln_full_vel_array(jbdy) .and. dta%ll_u2d ) dta%nread(2) = dta%nread(2) + 1
               ENDIF

               IF( dta%ll_v3d .OR. ( ln_full_vel_array(jbdy) .and. dta%ll_v2d ) ) THEN
                  if(lwp) write(numout,*) '++++++ reading in v3d field'
                  jfld = jfld + 1
                  blf_i(jfld) = bn_v3d
                  ibdy(jfld) = jbdy
                  igrid(jfld) = 3
                  ilen1(jfld) = nblen(igrid(jfld))
                  ilen3(jfld) = jpk
                  IF( ln_full_vel_array(jbdy) .and. dta%ll_v2d ) dta%nread(2) = dta%nread(2) + 1
               ENDIF

            ENDIF

            ! temperature and salinity
            IF( nn_tra_dta(jbdy) == 1 ) THEN

               IF( dta%ll_tem ) THEN
                  if(lwp) write(numout,*) '++++++ reading in tem field'
                  jfld = jfld + 1
                  blf_i(jfld) = bn_tem
                  ibdy(jfld) = jbdy
                  igrid(jfld) = 1
                  ilen1(jfld) = nblen(igrid(jfld))
                  ilen3(jfld) = jpk
               ENDIF

               IF( dta%ll_sal ) THEN
                  if(lwp) write(numout,*) '++++++ reading in sal field'
                  jfld = jfld + 1
                  blf_i(jfld) = bn_sal
                  ibdy(jfld) = jbdy
                  igrid(jfld) = 1
                  ilen1(jfld) = nblen(igrid(jfld))
                  ilen3(jfld) = jpk
               ENDIF

            ENDIF

# 677 "bdydta.F90"
            ! Recalculate field counts
            !-------------------------
            IF( jbdy == 1 ) THEN 
               nb_bdy_fld_sum = 0
               nb_bdy_fld(jbdy) = jfld
               nb_bdy_fld_sum     = jfld              
            ELSE
               nb_bdy_fld(jbdy) = jfld - nb_bdy_fld_sum
               nb_bdy_fld_sum = nb_bdy_fld_sum + nb_bdy_fld(jbdy)
            ENDIF

            dta%nread(1) = nb_bdy_fld(jbdy)

         ENDIF ! nn_dta == 1
      ENDDO ! jbdy

      DO jfld = 1, nb_bdy_fld_sum
         ALLOCATE( bf(jfld)%fnow(ilen1(jfld),1,ilen3(jfld)) )
         IF( blf_i(jfld)%ln_tint ) ALLOCATE( bf(jfld)%fdta(ilen1(jfld),1,ilen3(jfld),2) )
         nbmap_ptr(jfld)%ptr => idx_bdy(ibdy(jfld))%nbmap(:,igrid(jfld))
         nbmap_ptr(jfld)%ll_unstruc = ln_coords_file(ibdy(jfld))
      ENDDO

      ! fill bf with blf_i and control print
      !-------------------------------------
      jstart = 1
      DO jbdy = 1, nb_bdy
         jend = jstart - 1 + nb_bdy_fld(jbdy) 
         CALL fld_fill( bf(jstart:jend), blf_i(jstart:jend), cn_dir_array(jbdy), 'bdy_dta',   &
         &              'open boundary conditions', 'nambdy_dta' )
         jstart = jend + 1
      ENDDO

      DO jfld = 1, nb_bdy_fld_sum
         bf(jfld)%igrd = igrid(jfld)
         bf(jfld)%ibdy = ibdy(jfld)
      ENDDO

      ! Initialise local boundary data arrays
      ! nn_xxx_dta=0 : allocate space - will be filled from initial conditions later
      ! nn_xxx_dta=1 : point to "fnow" arrays
      !-------------------------------------

      jfld = 0
      DO jbdy=1, nb_bdy

         nblen => idx_bdy(jbdy)%nblen
         dta => dta_bdy(jbdy)

         if(lwp) then
            write(numout,*) '++++++ dta%ll_ssh = ',dta%ll_ssh
            write(numout,*) '++++++ dta%ll_u2d = ',dta%ll_u2d
            write(numout,*) '++++++ dta%ll_v2d = ',dta%ll_v2d
            write(numout,*) '++++++ dta%ll_u3d = ',dta%ll_u3d
            write(numout,*) '++++++ dta%ll_v3d = ',dta%ll_v3d
            write(numout,*) '++++++ dta%ll_tem = ',dta%ll_tem
            write(numout,*) '++++++ dta%ll_sal = ',dta%ll_sal
         endif

         IF ( nn_dyn2d_dta(jbdy) == 0 .or. nn_dyn2d_dta(jbdy) == 2 ) THEN
            if(lwp) write(numout,*) '++++++ dta%ssh/u2d/u3d allocated space'
            IF( dta%ll_ssh ) ALLOCATE( dta%ssh(nblen(1)) )
            IF( dta%ll_u2d ) ALLOCATE( dta%u2d(nblen(2)) )
            IF( dta%ll_v2d ) ALLOCATE( dta%v2d(nblen(3)) )
         ENDIF
         IF ( nn_dyn2d_dta(jbdy) == 1 .or. nn_dyn2d_dta(jbdy) == 3 ) THEN
            IF( dta%ll_ssh ) THEN
               if(lwp) write(numout,*) '++++++ dta%ssh pointing to fnow'
               jfld = jfld + 1
               dta%ssh => bf(jfld)%fnow(:,1,1)
            ENDIF
            IF ( dta%ll_u2d ) THEN
               IF ( ln_full_vel_array(jbdy) ) THEN
                  if(lwp) write(numout,*) '++++++ dta%u2d allocated space'
                  ALLOCATE( dta%u2d(nblen(2)) )
               ELSE
                  if(lwp) write(numout,*) '++++++ dta%u2d pointing to fnow'
                  jfld = jfld + 1
                  dta%u2d => bf(jfld)%fnow(:,1,1)
               ENDIF
            ENDIF
            IF ( dta%ll_v2d ) THEN
               IF ( ln_full_vel_array(jbdy) ) THEN
                  if(lwp) write(numout,*) '++++++ dta%v2d allocated space'
                  ALLOCATE( dta%v2d(nblen(3)) )
               ELSE
                  if(lwp) write(numout,*) '++++++ dta%v2d pointing to fnow'
                  jfld = jfld + 1
                  dta%v2d => bf(jfld)%fnow(:,1,1)
               ENDIF
            ENDIF
         ENDIF

         IF ( nn_dyn3d_dta(jbdy) == 0 ) THEN
            if(lwp) write(numout,*) '++++++ dta%u3d/v3d allocated space'
            IF( dta%ll_u3d ) ALLOCATE( dta_bdy(jbdy)%u3d(nblen(2),jpk) )
            IF( dta%ll_v3d ) ALLOCATE( dta_bdy(jbdy)%v3d(nblen(3),jpk) )
         ENDIF
         IF ( nn_dyn3d_dta(jbdy) == 1 .or. &
           &  ( ln_full_vel_array(jbdy) .and. ( nn_dyn2d_dta(jbdy) == 1 .or. nn_dyn2d_dta(jbdy) == 3 ) ) ) THEN
            IF ( dta%ll_u3d .or. ( ln_full_vel_array(jbdy) .and. dta%ll_u2d ) ) THEN
               if(lwp) write(numout,*) '++++++ dta%u3d pointing to fnow'
               jfld = jfld + 1
               dta_bdy(jbdy)%u3d => bf(jfld)%fnow(:,1,:)
            ENDIF
            IF ( dta%ll_v3d .or. ( ln_full_vel_array(jbdy) .and. dta%ll_v2d ) ) THEN
               if(lwp) write(numout,*) '++++++ dta%v3d pointing to fnow'
               jfld = jfld + 1
               dta_bdy(jbdy)%v3d => bf(jfld)%fnow(:,1,:)
            ENDIF
         ENDIF

         IF( nn_tra_dta(jbdy) == 0 ) THEN
            if(lwp) write(numout,*) '++++++ dta%tem/sal allocated space'
            IF( dta%ll_tem ) ALLOCATE( dta_bdy(jbdy)%tem(nblen(1),jpk) )
            IF( dta%ll_sal ) ALLOCATE( dta_bdy(jbdy)%sal(nblen(1),jpk) )
         ELSE
            IF( dta%ll_tem ) THEN
               if(lwp) write(numout,*) '++++++ dta%tem pointing to fnow'
               jfld = jfld + 1
               dta_bdy(jbdy)%tem => bf(jfld)%fnow(:,1,:)
            ENDIF
            IF( dta%ll_sal ) THEN 
               if(lwp) write(numout,*) '++++++ dta%sal pointing to fnow'
               jfld = jfld + 1
               dta_bdy(jbdy)%sal => bf(jfld)%fnow(:,1,:)
            ENDIF
         ENDIF

# 836 "bdydta.F90"
         !
      END DO ! jbdy 
      !
   END SUBROUTINE bdy_dta_init

   !!==============================================================================
END MODULE bdydta
