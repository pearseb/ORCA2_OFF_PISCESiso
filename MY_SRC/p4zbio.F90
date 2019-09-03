MODULE p4zbio
   !!======================================================================
   !!                         ***  MODULE p4zbio  ***
   !! TOP :   PISCES bio-model
   !!======================================================================
   !! History :   1.0  !  2004     (O. Aumont) Original code
   !!             2.0  !  2007-12  (C. Ethe, G. Madec)  F90
   !!----------------------------------------------------------------------
   !!   p4z_bio        :   computes the interactions between the different
   !!                      compartments of PISCES
   !!----------------------------------------------------------------------
   USE oce_trc         !  shared variables between ocean and passive tracers
   USE trc             !  passive tracers common variables 
   USE sms_pisces      !  PISCES Source Minus Sink variables
   USE p4zsink         !  vertical flux of particulate matter due to sinking
   USE p4zopt          !  optical model
   USE p4zlim          !  Co-limitations of differents nutrients
   USE p4zprod         !  Growth rate of the 2 phyto groups
   USE p4zmort         !  Mortality terms for phytoplankton
   USE p4zmicro        !  Sources and sinks of microzooplankton
   USE p4zmeso         !  Sources and sinks of mesozooplankton
   USE p5zlim          !  Co-limitations of differents nutrients
   USE p5zprod         !  Growth rate of the 2 phyto groups
   USE p5zmort         !  Mortality terms for phytoplankton
   USE p5zmicro        !  Sources and sinks of microzooplankton
   USE p5zmeso         !  Sources and sinks of mesozooplankton
   USE p4zrem          !  Remineralisation of organic matter
   USE p4zpoc          !  Remineralization of organic particles
   USE p4zagg          !  Aggregation of particles
   USE p4zfechem
   USE p4zligand       !  Prognostic ligand model
   USE prtctl_trc      !  print control for debugging
   USE iom             !  I/O manager
  
   IMPLICIT NONE
   PRIVATE

   PUBLIC  p4z_bio    

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: p4zbio.F90 10227 2018-10-25 14:42:24Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE p4z_bio ( kt, knt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE p4z_bio  ***
      !!
      !! ** Purpose :   Ecosystem model in the whole ocean: computes the
      !!              different interactions between the different compartments
      !!              of PISCES
      !!
      !! ** Method  : - ???
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) :: kt, knt
      !
      INTEGER             :: ji, jj, jk, jn
      LOGICAL :: diag
      INTEGER :: iii,jjj,kkk
      CHARACTER (len=25) :: charout
      !!---------------------------------------------------------------------
      !
      IF( ln_timing )   CALL timing_start('p4z_bio')
      !
      !     ASSIGN THE SHEAR RATE THAT IS USED FOR AGGREGATION
      !     OF PHYTOPLANKTON AND DETRITUS

      xdiss(:,:,:) = 1.
!!gm the use of nmld should be better here?
      DO jk = 2, jpkm1
         DO jj = 1, jpj
            DO ji = 1, jpi
!!gm  :  use nmln  and test on jk ...  less memory acces
               IF( gdepw_n(ji,jj,jk+1) > hmld(ji,jj) )   xdiss(ji,jj,jk) = 0.01
            END DO 
         END DO
      END DO

      diag = .false.
      iii = 117 
      jjj = 75
      kkk = 1

      CALL p4z_opt     ( kt, knt )     ! Optic: PAR in the water column
      IF( ln_c13 .and. diag ) THEN
        print*, " " 
        print*, " before call to sinking routine " 
        print*, "NO3 " 
        print*, trb(iii,jjj,kkk,jpdic),trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jpdic)-trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jp13dic)/trb(iii,jjj,kkk,jpdic)
        print*, trn(iii,jjj,kkk,jpdic),trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jpdic)-trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jp13dic)/trn(iii,jjj,kkk,jpdic)
        print*, tra(iii,jjj,kkk,jpdic),tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jpdic)-tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jp13dic)/tra(iii,jjj,kkk,jpdic)
        print*, "NH4 "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpcal),trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jpcal)-trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jp13cal)/trb(iii,jjj,kkk,jpcal)
        print*, trn(iii,jjj,kkk,jpcal),trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jpcal)-trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jp13cal)/trn(iii,jjj,kkk,jpcal)
        print*, tra(iii,jjj,kkk,jpcal),tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jpcal)-tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jp13cal)/tra(iii,jjj,kkk,jpcal)
        print*, "DOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpdoc),trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jpdoc)-trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jp13doc)/trb(iii,jjj,kkk,jpdoc)
        print*, trn(iii,jjj,kkk,jpdoc),trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jpdoc)-trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jp13doc)/trn(iii,jjj,kkk,jpdoc)
        print*, tra(iii,jjj,kkk,jpdoc),tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jpdoc)-tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jp13doc)/tra(iii,jjj,kkk,jpdoc)
        print*, "POC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jppoc),trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jppoc)-trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jp13poc)/trb(iii,jjj,kkk,jppoc)
        print*, trn(iii,jjj,kkk,jppoc),trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jppoc)-trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jp13poc)/trn(iii,jjj,kkk,jppoc)
        print*, tra(iii,jjj,kkk,jppoc),tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jppoc)-tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jp13poc)/tra(iii,jjj,kkk,jppoc)
        print*, "GOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpgoc),trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jpgoc)-trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jp13goc)/trb(iii,jjj,kkk,jpgoc)
        print*, trn(iii,jjj,kkk,jpgoc),trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jpgoc)-trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jp13goc)/trn(iii,jjj,kkk,jpgoc)
        print*, tra(iii,jjj,kkk,jpgoc),tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jpgoc)-tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jp13goc)/tra(iii,jjj,kkk,jpgoc)
        print*, "PHY "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpphy),trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jpphy)-trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jp13phy)/trb(iii,jjj,kkk,jpphy)
        print*, trn(iii,jjj,kkk,jpphy),trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jpphy)-trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jp13phy)/trn(iii,jjj,kkk,jpphy)
        print*, tra(iii,jjj,kkk,jpphy),tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jpphy)-tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jp13phy)/tra(iii,jjj,kkk,jpphy)
        print*, "PHY2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpdia),trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jpdia)-trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jp13dia)/trb(iii,jjj,kkk,jpdia)
        print*, trn(iii,jjj,kkk,jpdia),trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jpdia)-trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jp13dia)/trn(iii,jjj,kkk,jpdia)
        print*, tra(iii,jjj,kkk,jpdia),tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jpdia)-tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jp13dia)/tra(iii,jjj,kkk,jpdia)
        print*, "ZOO "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpzoo),trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jpzoo)-trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jp13zoo)/trb(iii,jjj,kkk,jpzoo)
        print*, trn(iii,jjj,kkk,jpzoo),trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jpzoo)-trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jp13zoo)/trn(iii,jjj,kkk,jpzoo)
        print*, tra(iii,jjj,kkk,jpzoo),tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jpzoo)-tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jp13zoo)/tra(iii,jjj,kkk,jpzoo)
        print*, "ZOO2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpmes),trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jpmes)-trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jp13mes)/trb(iii,jjj,kkk,jpmes)
        print*, trn(iii,jjj,kkk,jpmes),trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jpmes)-trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jp13mes)/trn(iii,jjj,kkk,jpmes)
        print*, tra(iii,jjj,kkk,jpmes),tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jpmes)-tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jp13mes)/tra(iii,jjj,kkk,jpmes)
        print*, " "
      ENDIF
      CALL p4z_sink    ( kt, knt )     ! vertical flux of particulate organic matter
      IF( ln_c13 .and. diag ) THEN
        print*, " " 
        print*, rtrn
        print*, " after call to sinking routine " 
        print*, "NO3 " 
        print*, trb(iii,jjj,kkk,jpdic),trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jpdic)-trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jp13dic)/trb(iii,jjj,kkk,jpdic)
        print*, trn(iii,jjj,kkk,jpdic),trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jpdic)-trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jp13dic)/trn(iii,jjj,kkk,jpdic)
        print*, tra(iii,jjj,kkk,jpdic),tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jpdic)-tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jp13dic)/tra(iii,jjj,kkk,jpdic)
        print*, "NH4 "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpcal),trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jpcal)-trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jp13cal)/trb(iii,jjj,kkk,jpcal)
        print*, trn(iii,jjj,kkk,jpcal),trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jpcal)-trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jp13cal)/trn(iii,jjj,kkk,jpcal)
        print*, tra(iii,jjj,kkk,jpcal),tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jpcal)-tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jp13cal)/tra(iii,jjj,kkk,jpcal)
        print*, "DOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpdoc),trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jpdoc)-trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jp13doc)/trb(iii,jjj,kkk,jpdoc)
        print*, trn(iii,jjj,kkk,jpdoc),trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jpdoc)-trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jp13doc)/trn(iii,jjj,kkk,jpdoc)
        print*, tra(iii,jjj,kkk,jpdoc),tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jpdoc)-tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jp13doc)/tra(iii,jjj,kkk,jpdoc)
        print*, "POC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jppoc),trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jppoc)-trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jp13poc)/trb(iii,jjj,kkk,jppoc)
        print*, trn(iii,jjj,kkk,jppoc),trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jppoc)-trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jp13poc)/trn(iii,jjj,kkk,jppoc)
        print*, tra(iii,jjj,kkk,jppoc),tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jppoc)-tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jp13poc)/tra(iii,jjj,kkk,jppoc)
        print*, "GOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpgoc),trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jpgoc)-trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jp13goc)/trb(iii,jjj,kkk,jpgoc)
        print*, trn(iii,jjj,kkk,jpgoc),trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jpgoc)-trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jp13goc)/trn(iii,jjj,kkk,jpgoc)
        print*, tra(iii,jjj,kkk,jpgoc),tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jpgoc)-tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jp13goc)/tra(iii,jjj,kkk,jpgoc)
        print*, "PHY "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpphy),trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jpphy)-trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jp13phy)/trb(iii,jjj,kkk,jpphy)
        print*, trn(iii,jjj,kkk,jpphy),trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jpphy)-trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jp13phy)/trn(iii,jjj,kkk,jpphy)
        print*, tra(iii,jjj,kkk,jpphy),tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jpphy)-tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jp13phy)/tra(iii,jjj,kkk,jpphy)
        print*, "PHY2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpdia),trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jpdia)-trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jp13dia)/trb(iii,jjj,kkk,jpdia)
        print*, trn(iii,jjj,kkk,jpdia),trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jpdia)-trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jp13dia)/trn(iii,jjj,kkk,jpdia)
        print*, tra(iii,jjj,kkk,jpdia),tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jpdia)-tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jp13dia)/tra(iii,jjj,kkk,jpdia)
        print*, "ZOO "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpzoo),trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jpzoo)-trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jp13zoo)/trb(iii,jjj,kkk,jpzoo)
        print*, trn(iii,jjj,kkk,jpzoo),trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jpzoo)-trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jp13zoo)/trn(iii,jjj,kkk,jpzoo)
        print*, tra(iii,jjj,kkk,jpzoo),tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jpzoo)-tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jp13zoo)/tra(iii,jjj,kkk,jpzoo)
        print*, "ZOO2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpmes),trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jpmes)-trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jp13mes)/trb(iii,jjj,kkk,jpmes)
        print*, trn(iii,jjj,kkk,jpmes),trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jpmes)-trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jp13mes)/trn(iii,jjj,kkk,jpmes)
        print*, tra(iii,jjj,kkk,jpmes),tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jpmes)-tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jp13mes)/tra(iii,jjj,kkk,jpmes)
        print*, " "
      ENDIF
      CALL p4z_fechem  ( kt, knt )     ! Iron chemistry/scavenging
      !
      IF( ln_c13 .and. diag ) THEN
        print*, " " 
        print*, rtrn
        print*, " before calls to biological routines " 
        print*, "NO3 " 
        print*, trb(iii,jjj,kkk,jpdic),trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jpdic)-trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jp13dic)/trb(iii,jjj,kkk,jpdic)
        print*, trn(iii,jjj,kkk,jpdic),trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jpdic)-trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jp13dic)/trn(iii,jjj,kkk,jpdic)
        print*, tra(iii,jjj,kkk,jpdic),tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jpdic)-tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jp13dic)/tra(iii,jjj,kkk,jpdic)
        print*, "NH4 "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpcal),trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jpcal)-trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jp13cal)/trb(iii,jjj,kkk,jpcal)
        print*, trn(iii,jjj,kkk,jpcal),trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jpcal)-trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jp13cal)/trn(iii,jjj,kkk,jpcal)
        print*, tra(iii,jjj,kkk,jpcal),tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jpcal)-tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jp13cal)/tra(iii,jjj,kkk,jpcal)
        print*, "DOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpdoc),trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jpdoc)-trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jp13doc)/trb(iii,jjj,kkk,jpdoc)
        print*, trn(iii,jjj,kkk,jpdoc),trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jpdoc)-trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jp13doc)/trn(iii,jjj,kkk,jpdoc)
        print*, tra(iii,jjj,kkk,jpdoc),tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jpdoc)-tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jp13doc)/tra(iii,jjj,kkk,jpdoc)
        print*, "POC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jppoc),trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jppoc)-trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jp13poc)/trb(iii,jjj,kkk,jppoc)
        print*, trn(iii,jjj,kkk,jppoc),trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jppoc)-trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jp13poc)/trn(iii,jjj,kkk,jppoc)
        print*, tra(iii,jjj,kkk,jppoc),tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jppoc)-tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jp13poc)/tra(iii,jjj,kkk,jppoc)
        print*, "GOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpgoc),trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jpgoc)-trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jp13goc)/trb(iii,jjj,kkk,jpgoc)
        print*, trn(iii,jjj,kkk,jpgoc),trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jpgoc)-trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jp13goc)/trn(iii,jjj,kkk,jpgoc)
        print*, tra(iii,jjj,kkk,jpgoc),tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jpgoc)-tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jp13goc)/tra(iii,jjj,kkk,jpgoc)
        print*, "PHY "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpphy),trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jpphy)-trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jp13phy)/trb(iii,jjj,kkk,jpphy)
        print*, trn(iii,jjj,kkk,jpphy),trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jpphy)-trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jp13phy)/trn(iii,jjj,kkk,jpphy)
        print*, tra(iii,jjj,kkk,jpphy),tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jpphy)-tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jp13phy)/tra(iii,jjj,kkk,jpphy)
        print*, "PHY2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpdia),trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jpdia)-trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jp13dia)/trb(iii,jjj,kkk,jpdia)
        print*, trn(iii,jjj,kkk,jpdia),trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jpdia)-trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jp13dia)/trn(iii,jjj,kkk,jpdia)
        print*, tra(iii,jjj,kkk,jpdia),tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jpdia)-tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jp13dia)/tra(iii,jjj,kkk,jpdia)
        print*, "ZOO "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpzoo),trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jpzoo)-trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jp13zoo)/trb(iii,jjj,kkk,jpzoo)
        print*, trn(iii,jjj,kkk,jpzoo),trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jpzoo)-trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jp13zoo)/trn(iii,jjj,kkk,jpzoo)
        print*, tra(iii,jjj,kkk,jpzoo),tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jpzoo)-tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jp13zoo)/tra(iii,jjj,kkk,jpzoo)
        print*, "ZOO2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpmes),trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jpmes)-trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jp13mes)/trb(iii,jjj,kkk,jpmes)
        print*, trn(iii,jjj,kkk,jpmes),trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jpmes)-trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jp13mes)/trn(iii,jjj,kkk,jpmes)
        print*, tra(iii,jjj,kkk,jpmes),tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jpmes)-tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jp13mes)/tra(iii,jjj,kkk,jpmes)
        print*, " "
      ENDIF
      IF( ln_p4z ) THEN
         CALL p4z_lim  ( kt, knt )     ! co-limitations by the various nutrients
      IF( ln_c13 .and. diag ) THEN
        print*, " " 
        print*, " after call to p4zlim routines " 
        print*, "NO3 " 
        print*, trb(iii,jjj,kkk,jpdic),trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jpdic)-trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jp13dic)/trb(iii,jjj,kkk,jpdic)
        print*, trn(iii,jjj,kkk,jpdic),trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jpdic)-trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jp13dic)/trn(iii,jjj,kkk,jpdic)
        print*, tra(iii,jjj,kkk,jpdic),tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jpdic)-tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jp13dic)/tra(iii,jjj,kkk,jpdic)
        print*, "NH4 "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpcal),trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jpcal)-trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jp13cal)/trb(iii,jjj,kkk,jpcal)
        print*, trn(iii,jjj,kkk,jpcal),trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jpcal)-trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jp13cal)/trn(iii,jjj,kkk,jpcal)
        print*, tra(iii,jjj,kkk,jpcal),tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jpcal)-tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jp13cal)/tra(iii,jjj,kkk,jpcal)
        print*, "DOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpdoc),trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jpdoc)-trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jp13doc)/trb(iii,jjj,kkk,jpdoc)
        print*, trn(iii,jjj,kkk,jpdoc),trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jpdoc)-trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jp13doc)/trn(iii,jjj,kkk,jpdoc)
        print*, tra(iii,jjj,kkk,jpdoc),tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jpdoc)-tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jp13doc)/tra(iii,jjj,kkk,jpdoc)
        print*, "POC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jppoc),trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jppoc)-trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jp13poc)/trb(iii,jjj,kkk,jppoc)
        print*, trn(iii,jjj,kkk,jppoc),trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jppoc)-trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jp13poc)/trn(iii,jjj,kkk,jppoc)
        print*, tra(iii,jjj,kkk,jppoc),tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jppoc)-tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jp13poc)/tra(iii,jjj,kkk,jppoc)
        print*, "GOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpgoc),trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jpgoc)-trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jp13goc)/trb(iii,jjj,kkk,jpgoc)
        print*, trn(iii,jjj,kkk,jpgoc),trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jpgoc)-trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jp13goc)/trn(iii,jjj,kkk,jpgoc)
        print*, tra(iii,jjj,kkk,jpgoc),tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jpgoc)-tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jp13goc)/tra(iii,jjj,kkk,jpgoc)
        print*, "PHY "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpphy),trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jpphy)-trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jp13phy)/trb(iii,jjj,kkk,jpphy)
        print*, trn(iii,jjj,kkk,jpphy),trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jpphy)-trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jp13phy)/trn(iii,jjj,kkk,jpphy)
        print*, tra(iii,jjj,kkk,jpphy),tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jpphy)-tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jp13phy)/tra(iii,jjj,kkk,jpphy)
        print*, "PHY2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpdia),trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jpdia)-trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jp13dia)/trb(iii,jjj,kkk,jpdia)
        print*, trn(iii,jjj,kkk,jpdia),trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jpdia)-trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jp13dia)/trn(iii,jjj,kkk,jpdia)
        print*, tra(iii,jjj,kkk,jpdia),tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jpdia)-tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jp13dia)/tra(iii,jjj,kkk,jpdia)
        print*, "ZOO "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpzoo),trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jpzoo)-trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jp13zoo)/trb(iii,jjj,kkk,jpzoo)
        print*, trn(iii,jjj,kkk,jpzoo),trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jpzoo)-trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jp13zoo)/trn(iii,jjj,kkk,jpzoo)
        print*, tra(iii,jjj,kkk,jpzoo),tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jpzoo)-tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jp13zoo)/tra(iii,jjj,kkk,jpzoo)
        print*, "ZOO2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpmes),trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jpmes)-trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jp13mes)/trb(iii,jjj,kkk,jpmes)
        print*, trn(iii,jjj,kkk,jpmes),trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jpmes)-trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jp13mes)/trn(iii,jjj,kkk,jpmes)
        print*, tra(iii,jjj,kkk,jpmes),tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jpmes)-tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jp13mes)/tra(iii,jjj,kkk,jpmes)
        print*, " "
      ENDIF
         CALL p4z_prod ( kt, knt )     ! phytoplankton growth rate over the global ocean. 
      IF( ln_c13 .and. diag ) THEN
        print*, " " 
        print*, " after calls to p4zprod routines " 
        print*, "NO3 " 
        print*, trb(iii,jjj,kkk,jpdic),trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jpdic)-trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jp13dic)/trb(iii,jjj,kkk,jpdic)
        print*, trn(iii,jjj,kkk,jpdic),trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jpdic)-trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jp13dic)/trn(iii,jjj,kkk,jpdic)
        print*, tra(iii,jjj,kkk,jpdic),tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jpdic)-tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jp13dic)/tra(iii,jjj,kkk,jpdic)
        print*, "NH4 "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpcal),trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jpcal)-trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jp13cal)/trb(iii,jjj,kkk,jpcal)
        print*, trn(iii,jjj,kkk,jpcal),trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jpcal)-trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jp13cal)/trn(iii,jjj,kkk,jpcal)
        print*, tra(iii,jjj,kkk,jpcal),tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jpcal)-tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jp13cal)/tra(iii,jjj,kkk,jpcal)
        print*, "DOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpdoc),trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jpdoc)-trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jp13doc)/trb(iii,jjj,kkk,jpdoc)
        print*, trn(iii,jjj,kkk,jpdoc),trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jpdoc)-trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jp13doc)/trn(iii,jjj,kkk,jpdoc)
        print*, tra(iii,jjj,kkk,jpdoc),tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jpdoc)-tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jp13doc)/tra(iii,jjj,kkk,jpdoc)
        print*, "POC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jppoc),trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jppoc)-trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jp13poc)/trb(iii,jjj,kkk,jppoc)
        print*, trn(iii,jjj,kkk,jppoc),trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jppoc)-trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jp13poc)/trn(iii,jjj,kkk,jppoc)
        print*, tra(iii,jjj,kkk,jppoc),tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jppoc)-tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jp13poc)/tra(iii,jjj,kkk,jppoc)
        print*, "GOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpgoc),trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jpgoc)-trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jp13goc)/trb(iii,jjj,kkk,jpgoc)
        print*, trn(iii,jjj,kkk,jpgoc),trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jpgoc)-trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jp13goc)/trn(iii,jjj,kkk,jpgoc)
        print*, tra(iii,jjj,kkk,jpgoc),tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jpgoc)-tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jp13goc)/tra(iii,jjj,kkk,jpgoc)
        print*, "PHY "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpphy),trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jpphy)-trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jp13phy)/trb(iii,jjj,kkk,jpphy)
        print*, trn(iii,jjj,kkk,jpphy),trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jpphy)-trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jp13phy)/trn(iii,jjj,kkk,jpphy)
        print*, tra(iii,jjj,kkk,jpphy),tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jpphy)-tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jp13phy)/tra(iii,jjj,kkk,jpphy)
        print*, "PHY2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpdia),trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jpdia)-trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jp13dia)/trb(iii,jjj,kkk,jpdia)
        print*, trn(iii,jjj,kkk,jpdia),trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jpdia)-trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jp13dia)/trn(iii,jjj,kkk,jpdia)
        print*, tra(iii,jjj,kkk,jpdia),tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jpdia)-tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jp13dia)/tra(iii,jjj,kkk,jpdia)
        print*, "ZOO "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpzoo),trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jpzoo)-trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jp13zoo)/trb(iii,jjj,kkk,jpzoo)
        print*, trn(iii,jjj,kkk,jpzoo),trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jpzoo)-trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jp13zoo)/trn(iii,jjj,kkk,jpzoo)
        print*, tra(iii,jjj,kkk,jpzoo),tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jpzoo)-tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jp13zoo)/tra(iii,jjj,kkk,jpzoo)
        print*, "ZOO2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpmes),trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jpmes)-trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jp13mes)/trb(iii,jjj,kkk,jpmes)
        print*, trn(iii,jjj,kkk,jpmes),trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jpmes)-trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jp13mes)/trn(iii,jjj,kkk,jpmes)
        print*, tra(iii,jjj,kkk,jpmes),tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jpmes)-tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jp13mes)/tra(iii,jjj,kkk,jpmes)
        print*, " "
      ENDIF
         !                             ! (for each element : C, Si, Fe, Chl )
         CALL p4z_mort ( kt      )     ! phytoplankton mortality
      IF( ln_c13 .and. diag ) THEN
        print*, " " 
        print*, " after calls to p4zmort routines " 
        print*, "NO3 " 
        print*, trb(iii,jjj,kkk,jpdic),trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jpdic)-trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jp13dic)/trb(iii,jjj,kkk,jpdic)
        print*, trn(iii,jjj,kkk,jpdic),trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jpdic)-trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jp13dic)/trn(iii,jjj,kkk,jpdic)
        print*, tra(iii,jjj,kkk,jpdic),tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jpdic)-tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jp13dic)/tra(iii,jjj,kkk,jpdic)
        print*, "NH4 "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpcal),trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jpcal)-trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jp13cal)/trb(iii,jjj,kkk,jpcal)
        print*, trn(iii,jjj,kkk,jpcal),trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jpcal)-trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jp13cal)/trn(iii,jjj,kkk,jpcal)
        print*, tra(iii,jjj,kkk,jpcal),tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jpcal)-tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jp13cal)/tra(iii,jjj,kkk,jpcal)
        print*, "DOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpdoc),trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jpdoc)-trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jp13doc)/trb(iii,jjj,kkk,jpdoc)
        print*, trn(iii,jjj,kkk,jpdoc),trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jpdoc)-trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jp13doc)/trn(iii,jjj,kkk,jpdoc)
        print*, tra(iii,jjj,kkk,jpdoc),tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jpdoc)-tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jp13doc)/tra(iii,jjj,kkk,jpdoc)
        print*, "POC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jppoc),trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jppoc)-trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jp13poc)/trb(iii,jjj,kkk,jppoc)
        print*, trn(iii,jjj,kkk,jppoc),trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jppoc)-trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jp13poc)/trn(iii,jjj,kkk,jppoc)
        print*, tra(iii,jjj,kkk,jppoc),tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jppoc)-tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jp13poc)/tra(iii,jjj,kkk,jppoc)
        print*, "GOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpgoc),trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jpgoc)-trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jp13goc)/trb(iii,jjj,kkk,jpgoc)
        print*, trn(iii,jjj,kkk,jpgoc),trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jpgoc)-trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jp13goc)/trn(iii,jjj,kkk,jpgoc)
        print*, tra(iii,jjj,kkk,jpgoc),tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jpgoc)-tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jp13goc)/tra(iii,jjj,kkk,jpgoc)
        print*, "PHY "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpphy),trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jpphy)-trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jp13phy)/trb(iii,jjj,kkk,jpphy)
        print*, trn(iii,jjj,kkk,jpphy),trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jpphy)-trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jp13phy)/trn(iii,jjj,kkk,jpphy)
        print*, tra(iii,jjj,kkk,jpphy),tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jpphy)-tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jp13phy)/tra(iii,jjj,kkk,jpphy)
        print*, "PHY2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpdia),trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jpdia)-trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jp13dia)/trb(iii,jjj,kkk,jpdia)
        print*, trn(iii,jjj,kkk,jpdia),trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jpdia)-trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jp13dia)/trn(iii,jjj,kkk,jpdia)
        print*, tra(iii,jjj,kkk,jpdia),tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jpdia)-tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jp13dia)/tra(iii,jjj,kkk,jpdia)
        print*, "ZOO "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpzoo),trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jpzoo)-trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jp13zoo)/trb(iii,jjj,kkk,jpzoo)
        print*, trn(iii,jjj,kkk,jpzoo),trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jpzoo)-trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jp13zoo)/trn(iii,jjj,kkk,jpzoo)
        print*, tra(iii,jjj,kkk,jpzoo),tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jpzoo)-tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jp13zoo)/tra(iii,jjj,kkk,jpzoo)
        print*, "ZOO2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpmes),trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jpmes)-trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jp13mes)/trb(iii,jjj,kkk,jpmes)
        print*, trn(iii,jjj,kkk,jpmes),trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jpmes)-trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jp13mes)/trn(iii,jjj,kkk,jpmes)
        print*, tra(iii,jjj,kkk,jpmes),tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jpmes)-tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jp13mes)/tra(iii,jjj,kkk,jpmes)
        print*, " "
      ENDIF
         !                             ! zooplankton sources/sinks routines 
         CALL p4z_micro( kt, knt )           ! microzooplankton
      IF( ln_c13 .and. diag ) THEN
        print*, " " 
        print*, " after calls to p4zmicro routines " 
        print*, "NO3 " 
        print*, trb(iii,jjj,kkk,jpdic),trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jpdic)-trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jp13dic)/trb(iii,jjj,kkk,jpdic)
        print*, trn(iii,jjj,kkk,jpdic),trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jpdic)-trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jp13dic)/trn(iii,jjj,kkk,jpdic)
        print*, tra(iii,jjj,kkk,jpdic),tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jpdic)-tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jp13dic)/tra(iii,jjj,kkk,jpdic)
        print*, "NH4 "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpcal),trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jpcal)-trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jp13cal)/trb(iii,jjj,kkk,jpcal)
        print*, trn(iii,jjj,kkk,jpcal),trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jpcal)-trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jp13cal)/trn(iii,jjj,kkk,jpcal)
        print*, tra(iii,jjj,kkk,jpcal),tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jpcal)-tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jp13cal)/tra(iii,jjj,kkk,jpcal)
        print*, "DOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpdoc),trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jpdoc)-trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jp13doc)/trb(iii,jjj,kkk,jpdoc)
        print*, trn(iii,jjj,kkk,jpdoc),trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jpdoc)-trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jp13doc)/trn(iii,jjj,kkk,jpdoc)
        print*, tra(iii,jjj,kkk,jpdoc),tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jpdoc)-tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jp13doc)/tra(iii,jjj,kkk,jpdoc)
        print*, "POC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jppoc),trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jppoc)-trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jp13poc)/trb(iii,jjj,kkk,jppoc)
        print*, trn(iii,jjj,kkk,jppoc),trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jppoc)-trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jp13poc)/trn(iii,jjj,kkk,jppoc)
        print*, tra(iii,jjj,kkk,jppoc),tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jppoc)-tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jp13poc)/tra(iii,jjj,kkk,jppoc)
        print*, "GOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpgoc),trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jpgoc)-trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jp13goc)/trb(iii,jjj,kkk,jpgoc)
        print*, trn(iii,jjj,kkk,jpgoc),trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jpgoc)-trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jp13goc)/trn(iii,jjj,kkk,jpgoc)
        print*, tra(iii,jjj,kkk,jpgoc),tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jpgoc)-tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jp13goc)/tra(iii,jjj,kkk,jpgoc)
        print*, "PHY "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpphy),trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jpphy)-trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jp13phy)/trb(iii,jjj,kkk,jpphy)
        print*, trn(iii,jjj,kkk,jpphy),trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jpphy)-trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jp13phy)/trn(iii,jjj,kkk,jpphy)
        print*, tra(iii,jjj,kkk,jpphy),tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jpphy)-tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jp13phy)/tra(iii,jjj,kkk,jpphy)
        print*, "PHY2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpdia),trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jpdia)-trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jp13dia)/trb(iii,jjj,kkk,jpdia)
        print*, trn(iii,jjj,kkk,jpdia),trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jpdia)-trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jp13dia)/trn(iii,jjj,kkk,jpdia)
        print*, tra(iii,jjj,kkk,jpdia),tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jpdia)-tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jp13dia)/tra(iii,jjj,kkk,jpdia)
        print*, "ZOO "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpzoo),trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jpzoo)-trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jp13zoo)/trb(iii,jjj,kkk,jpzoo)
        print*, trn(iii,jjj,kkk,jpzoo),trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jpzoo)-trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jp13zoo)/trn(iii,jjj,kkk,jpzoo)
        print*, tra(iii,jjj,kkk,jpzoo),tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jpzoo)-tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jp13zoo)/tra(iii,jjj,kkk,jpzoo)
        print*, "ZOO2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpmes),trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jpmes)-trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jp13mes)/trb(iii,jjj,kkk,jpmes)
        print*, trn(iii,jjj,kkk,jpmes),trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jpmes)-trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jp13mes)/trn(iii,jjj,kkk,jpmes)
        print*, tra(iii,jjj,kkk,jpmes),tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jpmes)-tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jp13mes)/tra(iii,jjj,kkk,jpmes)
        print*, " "
      ENDIF
         CALL p4z_meso ( kt, knt )           ! mesozooplankton
      IF( ln_c13 .and. diag ) THEN
        print*, " " 
        print*, " after calls to p4zmeso routines " 
        print*, "NO3 " 
        print*, trb(iii,jjj,kkk,jpdic),trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jpdic)-trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jp13dic)/trb(iii,jjj,kkk,jpdic)
        print*, trn(iii,jjj,kkk,jpdic),trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jpdic)-trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jp13dic)/trn(iii,jjj,kkk,jpdic)
        print*, tra(iii,jjj,kkk,jpdic),tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jpdic)-tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jp13dic)/tra(iii,jjj,kkk,jpdic)
        print*, "NH4 "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpcal),trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jpcal)-trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jp13cal)/trb(iii,jjj,kkk,jpcal)
        print*, trn(iii,jjj,kkk,jpcal),trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jpcal)-trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jp13cal)/trn(iii,jjj,kkk,jpcal)
        print*, tra(iii,jjj,kkk,jpcal),tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jpcal)-tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jp13cal)/tra(iii,jjj,kkk,jpcal)
        print*, "DOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpdoc),trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jpdoc)-trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jp13doc)/trb(iii,jjj,kkk,jpdoc)
        print*, trn(iii,jjj,kkk,jpdoc),trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jpdoc)-trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jp13doc)/trn(iii,jjj,kkk,jpdoc)
        print*, tra(iii,jjj,kkk,jpdoc),tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jpdoc)-tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jp13doc)/tra(iii,jjj,kkk,jpdoc)
        print*, "POC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jppoc),trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jppoc)-trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jp13poc)/trb(iii,jjj,kkk,jppoc)
        print*, trn(iii,jjj,kkk,jppoc),trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jppoc)-trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jp13poc)/trn(iii,jjj,kkk,jppoc)
        print*, tra(iii,jjj,kkk,jppoc),tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jppoc)-tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jp13poc)/tra(iii,jjj,kkk,jppoc)
        print*, "GOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpgoc),trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jpgoc)-trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jp13goc)/trb(iii,jjj,kkk,jpgoc)
        print*, trn(iii,jjj,kkk,jpgoc),trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jpgoc)-trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jp13goc)/trn(iii,jjj,kkk,jpgoc)
        print*, tra(iii,jjj,kkk,jpgoc),tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jpgoc)-tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jp13goc)/tra(iii,jjj,kkk,jpgoc)
        print*, "PHY "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpphy),trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jpphy)-trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jp13phy)/trb(iii,jjj,kkk,jpphy)
        print*, trn(iii,jjj,kkk,jpphy),trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jpphy)-trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jp13phy)/trn(iii,jjj,kkk,jpphy)
        print*, tra(iii,jjj,kkk,jpphy),tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jpphy)-tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jp13phy)/tra(iii,jjj,kkk,jpphy)
        print*, "PHY2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpdia),trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jpdia)-trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jp13dia)/trb(iii,jjj,kkk,jpdia)
        print*, trn(iii,jjj,kkk,jpdia),trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jpdia)-trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jp13dia)/trn(iii,jjj,kkk,jpdia)
        print*, tra(iii,jjj,kkk,jpdia),tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jpdia)-tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jp13dia)/tra(iii,jjj,kkk,jpdia)
        print*, "ZOO "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpzoo),trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jpzoo)-trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jp13zoo)/trb(iii,jjj,kkk,jpzoo)
        print*, trn(iii,jjj,kkk,jpzoo),trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jpzoo)-trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jp13zoo)/trn(iii,jjj,kkk,jpzoo)
        print*, tra(iii,jjj,kkk,jpzoo),tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jpzoo)-tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jp13zoo)/tra(iii,jjj,kkk,jpzoo)
        print*, "ZOO2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpmes),trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jpmes)-trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jp13mes)/trb(iii,jjj,kkk,jpmes)
        print*, trn(iii,jjj,kkk,jpmes),trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jpmes)-trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jp13mes)/trn(iii,jjj,kkk,jpmes)
        print*, tra(iii,jjj,kkk,jpmes),tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jpmes)-tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jp13mes)/tra(iii,jjj,kkk,jpmes)
        print*, " "
      ENDIF
      ELSE
         CALL p5z_lim  ( kt, knt )     ! co-limitations by the various nutrients
         CALL p5z_prod ( kt, knt )     ! phytoplankton growth rate over the global ocean. 
         !                             ! (for each element : C, Si, Fe, Chl )
         CALL p5z_mort ( kt      )     ! phytoplankton mortality
         !                             ! zooplankton sources/sinks routines 
         CALL p5z_micro( kt, knt )           ! microzooplankton
         CALL p5z_meso ( kt, knt )           ! mesozooplankton
      ENDIF
      IF( ln_c13 .and. diag ) THEN
        print*, " " 
        print*, " after calls to biological routines " 
        print*, "NO3 " 
        print*, trb(iii,jjj,kkk,jpdic),trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jpdic)-trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jp13dic)/trb(iii,jjj,kkk,jpdic)
        print*, trn(iii,jjj,kkk,jpdic),trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jpdic)-trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jp13dic)/trn(iii,jjj,kkk,jpdic)
        print*, tra(iii,jjj,kkk,jpdic),tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jpdic)-tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jp13dic)/tra(iii,jjj,kkk,jpdic)
        print*, "NH4 "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpcal),trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jpcal)-trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jp13cal)/trb(iii,jjj,kkk,jpcal)
        print*, trn(iii,jjj,kkk,jpcal),trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jpcal)-trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jp13cal)/trn(iii,jjj,kkk,jpcal)
        print*, tra(iii,jjj,kkk,jpcal),tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jpcal)-tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jp13cal)/tra(iii,jjj,kkk,jpcal)
        print*, "DOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpdoc),trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jpdoc)-trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jp13doc)/trb(iii,jjj,kkk,jpdoc)
        print*, trn(iii,jjj,kkk,jpdoc),trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jpdoc)-trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jp13doc)/trn(iii,jjj,kkk,jpdoc)
        print*, tra(iii,jjj,kkk,jpdoc),tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jpdoc)-tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jp13doc)/tra(iii,jjj,kkk,jpdoc)
        print*, "POC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jppoc),trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jppoc)-trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jp13poc)/trb(iii,jjj,kkk,jppoc)
        print*, trn(iii,jjj,kkk,jppoc),trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jppoc)-trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jp13poc)/trn(iii,jjj,kkk,jppoc)
        print*, tra(iii,jjj,kkk,jppoc),tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jppoc)-tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jp13poc)/tra(iii,jjj,kkk,jppoc)
        print*, "GOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpgoc),trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jpgoc)-trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jp13goc)/trb(iii,jjj,kkk,jpgoc)
        print*, trn(iii,jjj,kkk,jpgoc),trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jpgoc)-trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jp13goc)/trn(iii,jjj,kkk,jpgoc)
        print*, tra(iii,jjj,kkk,jpgoc),tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jpgoc)-tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jp13goc)/tra(iii,jjj,kkk,jpgoc)
        print*, "PHY "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpphy),trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jpphy)-trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jp13phy)/trb(iii,jjj,kkk,jpphy)
        print*, trn(iii,jjj,kkk,jpphy),trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jpphy)-trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jp13phy)/trn(iii,jjj,kkk,jpphy)
        print*, tra(iii,jjj,kkk,jpphy),tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jpphy)-tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jp13phy)/tra(iii,jjj,kkk,jpphy)
        print*, "PHY2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpdia),trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jpdia)-trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jp13dia)/trb(iii,jjj,kkk,jpdia)
        print*, trn(iii,jjj,kkk,jpdia),trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jpdia)-trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jp13dia)/trn(iii,jjj,kkk,jpdia)
        print*, tra(iii,jjj,kkk,jpdia),tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jpdia)-tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jp13dia)/tra(iii,jjj,kkk,jpdia)
        print*, "ZOO "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpzoo),trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jpzoo)-trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jp13zoo)/trb(iii,jjj,kkk,jpzoo)
        print*, trn(iii,jjj,kkk,jpzoo),trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jpzoo)-trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jp13zoo)/trn(iii,jjj,kkk,jpzoo)
        print*, tra(iii,jjj,kkk,jpzoo),tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jpzoo)-tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jp13zoo)/tra(iii,jjj,kkk,jpzoo)
        print*, "ZOO2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpmes),trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jpmes)-trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jp13mes)/trb(iii,jjj,kkk,jpmes)
        print*, trn(iii,jjj,kkk,jpmes),trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jpmes)-trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jp13mes)/trn(iii,jjj,kkk,jpmes)
        print*, tra(iii,jjj,kkk,jpmes),tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jpmes)-tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jp13mes)/tra(iii,jjj,kkk,jpmes)
        print*, " "
      ENDIF
      !
      CALL p4z_agg     ( kt, knt )     ! Aggregation of particles
      IF( ln_c13 .and. diag ) THEN
        print*, " " 
        print*, " after call to p4zagg " 
        print*, "NO3 " 
        print*, trb(iii,jjj,kkk,jpdic),trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jpdic)-trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jp13dic)/trb(iii,jjj,kkk,jpdic)
        print*, trn(iii,jjj,kkk,jpdic),trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jpdic)-trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jp13dic)/trn(iii,jjj,kkk,jpdic)
        print*, tra(iii,jjj,kkk,jpdic),tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jpdic)-tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jp13dic)/tra(iii,jjj,kkk,jpdic)
        print*, "NH4 "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpcal),trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jpcal)-trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jp13cal)/trb(iii,jjj,kkk,jpcal)
        print*, trn(iii,jjj,kkk,jpcal),trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jpcal)-trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jp13cal)/trn(iii,jjj,kkk,jpcal)
        print*, tra(iii,jjj,kkk,jpcal),tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jpcal)-tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jp13cal)/tra(iii,jjj,kkk,jpcal)
        print*, "DOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpdoc),trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jpdoc)-trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jp13doc)/trb(iii,jjj,kkk,jpdoc)
        print*, trn(iii,jjj,kkk,jpdoc),trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jpdoc)-trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jp13doc)/trn(iii,jjj,kkk,jpdoc)
        print*, tra(iii,jjj,kkk,jpdoc),tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jpdoc)-tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jp13doc)/tra(iii,jjj,kkk,jpdoc)
        print*, "POC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jppoc),trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jppoc)-trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jp13poc)/trb(iii,jjj,kkk,jppoc)
        print*, trn(iii,jjj,kkk,jppoc),trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jppoc)-trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jp13poc)/trn(iii,jjj,kkk,jppoc)
        print*, tra(iii,jjj,kkk,jppoc),tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jppoc)-tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jp13poc)/tra(iii,jjj,kkk,jppoc)
        print*, "GOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpgoc),trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jpgoc)-trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jp13goc)/trb(iii,jjj,kkk,jpgoc)
        print*, trn(iii,jjj,kkk,jpgoc),trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jpgoc)-trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jp13goc)/trn(iii,jjj,kkk,jpgoc)
        print*, tra(iii,jjj,kkk,jpgoc),tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jpgoc)-tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jp13goc)/tra(iii,jjj,kkk,jpgoc)
        print*, "PHY "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpphy),trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jpphy)-trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jp13phy)/trb(iii,jjj,kkk,jpphy)
        print*, trn(iii,jjj,kkk,jpphy),trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jpphy)-trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jp13phy)/trn(iii,jjj,kkk,jpphy)
        print*, tra(iii,jjj,kkk,jpphy),tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jpphy)-tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jp13phy)/tra(iii,jjj,kkk,jpphy)
        print*, "PHY2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpdia),trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jpdia)-trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jp13dia)/trb(iii,jjj,kkk,jpdia)
        print*, trn(iii,jjj,kkk,jpdia),trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jpdia)-trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jp13dia)/trn(iii,jjj,kkk,jpdia)
        print*, tra(iii,jjj,kkk,jpdia),tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jpdia)-tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jp13dia)/tra(iii,jjj,kkk,jpdia)
        print*, "ZOO "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpzoo),trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jpzoo)-trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jp13zoo)/trb(iii,jjj,kkk,jpzoo)
        print*, trn(iii,jjj,kkk,jpzoo),trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jpzoo)-trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jp13zoo)/trn(iii,jjj,kkk,jpzoo)
        print*, tra(iii,jjj,kkk,jpzoo),tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jpzoo)-tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jp13zoo)/tra(iii,jjj,kkk,jpzoo)
        print*, "ZOO2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpmes),trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jpmes)-trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jp13mes)/trb(iii,jjj,kkk,jpmes)
        print*, trn(iii,jjj,kkk,jpmes),trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jpmes)-trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jp13mes)/trn(iii,jjj,kkk,jpmes)
        print*, tra(iii,jjj,kkk,jpmes),tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jpmes)-tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jp13mes)/tra(iii,jjj,kkk,jpmes)
        print*, " "
      ENDIF
      CALL p4z_rem     ( kt, knt )     ! remineralization terms of organic matter+scavenging of Fe
      IF( ln_c13 .and. diag ) THEN
        print*, " " 
        print*, " after call to p4zrem " 
        print*, "NO3 " 
        print*, trb(iii,jjj,kkk,jpdic),trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jpdic)-trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jp13dic)/trb(iii,jjj,kkk,jpdic)
        print*, trn(iii,jjj,kkk,jpdic),trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jpdic)-trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jp13dic)/trn(iii,jjj,kkk,jpdic)
        print*, tra(iii,jjj,kkk,jpdic),tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jpdic)-tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jp13dic)/tra(iii,jjj,kkk,jpdic)
        print*, "NH4 "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpcal),trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jpcal)-trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jp13cal)/trb(iii,jjj,kkk,jpcal)
        print*, trn(iii,jjj,kkk,jpcal),trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jpcal)-trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jp13cal)/trn(iii,jjj,kkk,jpcal)
        print*, tra(iii,jjj,kkk,jpcal),tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jpcal)-tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jp13cal)/tra(iii,jjj,kkk,jpcal)
        print*, "DOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpdoc),trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jpdoc)-trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jp13doc)/trb(iii,jjj,kkk,jpdoc)
        print*, trn(iii,jjj,kkk,jpdoc),trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jpdoc)-trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jp13doc)/trn(iii,jjj,kkk,jpdoc)
        print*, tra(iii,jjj,kkk,jpdoc),tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jpdoc)-tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jp13doc)/tra(iii,jjj,kkk,jpdoc)
        print*, "POC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jppoc),trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jppoc)-trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jp13poc)/trb(iii,jjj,kkk,jppoc)
        print*, trn(iii,jjj,kkk,jppoc),trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jppoc)-trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jp13poc)/trn(iii,jjj,kkk,jppoc)
        print*, tra(iii,jjj,kkk,jppoc),tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jppoc)-tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jp13poc)/tra(iii,jjj,kkk,jppoc)
        print*, "GOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpgoc),trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jpgoc)-trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jp13goc)/trb(iii,jjj,kkk,jpgoc)
        print*, trn(iii,jjj,kkk,jpgoc),trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jpgoc)-trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jp13goc)/trn(iii,jjj,kkk,jpgoc)
        print*, tra(iii,jjj,kkk,jpgoc),tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jpgoc)-tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jp13goc)/tra(iii,jjj,kkk,jpgoc)
        print*, "PHY "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpphy),trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jpphy)-trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jp13phy)/trb(iii,jjj,kkk,jpphy)
        print*, trn(iii,jjj,kkk,jpphy),trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jpphy)-trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jp13phy)/trn(iii,jjj,kkk,jpphy)
        print*, tra(iii,jjj,kkk,jpphy),tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jpphy)-tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jp13phy)/tra(iii,jjj,kkk,jpphy)
        print*, "PHY2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpdia),trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jpdia)-trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jp13dia)/trb(iii,jjj,kkk,jpdia)
        print*, trn(iii,jjj,kkk,jpdia),trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jpdia)-trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jp13dia)/trn(iii,jjj,kkk,jpdia)
        print*, tra(iii,jjj,kkk,jpdia),tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jpdia)-tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jp13dia)/tra(iii,jjj,kkk,jpdia)
        print*, "ZOO "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpzoo),trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jpzoo)-trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jp13zoo)/trb(iii,jjj,kkk,jpzoo)
        print*, trn(iii,jjj,kkk,jpzoo),trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jpzoo)-trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jp13zoo)/trn(iii,jjj,kkk,jpzoo)
        print*, tra(iii,jjj,kkk,jpzoo),tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jpzoo)-tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jp13zoo)/tra(iii,jjj,kkk,jpzoo)
        print*, "ZOO2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpmes),trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jpmes)-trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jp13mes)/trb(iii,jjj,kkk,jpmes)
        print*, trn(iii,jjj,kkk,jpmes),trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jpmes)-trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jp13mes)/trn(iii,jjj,kkk,jpmes)
        print*, tra(iii,jjj,kkk,jpmes),tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jpmes)-tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jp13mes)/tra(iii,jjj,kkk,jpmes)
        print*, " "
      ENDIF
      CALL p4z_poc     ( kt, knt )     ! Remineralization of organic particles
      IF( ln_c13 .and. diag ) THEN
        print*, " " 
        print*, " after call to p4zpoc " 
        print*, "NO3 " 
        print*, trb(iii,jjj,kkk,jpdic),trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jpdic)-trb(iii,jjj,kkk,jp13dic),trb(iii,jjj,kkk,jp13dic)/trb(iii,jjj,kkk,jpdic)
        print*, trn(iii,jjj,kkk,jpdic),trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jpdic)-trn(iii,jjj,kkk,jp13dic),trn(iii,jjj,kkk,jp13dic)/trn(iii,jjj,kkk,jpdic)
        print*, tra(iii,jjj,kkk,jpdic),tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jpdic)-tra(iii,jjj,kkk,jp13dic),tra(iii,jjj,kkk,jp13dic)/tra(iii,jjj,kkk,jpdic)
        print*, "NH4 "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpcal),trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jpcal)-trb(iii,jjj,kkk,jp13cal),trb(iii,jjj,kkk,jp13cal)/trb(iii,jjj,kkk,jpcal)
        print*, trn(iii,jjj,kkk,jpcal),trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jpcal)-trn(iii,jjj,kkk,jp13cal),trn(iii,jjj,kkk,jp13cal)/trn(iii,jjj,kkk,jpcal)
        print*, tra(iii,jjj,kkk,jpcal),tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jpcal)-tra(iii,jjj,kkk,jp13cal),tra(iii,jjj,kkk,jp13cal)/tra(iii,jjj,kkk,jpcal)
        print*, "DOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpdoc),trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jpdoc)-trb(iii,jjj,kkk,jp13doc),trb(iii,jjj,kkk,jp13doc)/trb(iii,jjj,kkk,jpdoc)
        print*, trn(iii,jjj,kkk,jpdoc),trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jpdoc)-trn(iii,jjj,kkk,jp13doc),trn(iii,jjj,kkk,jp13doc)/trn(iii,jjj,kkk,jpdoc)
        print*, tra(iii,jjj,kkk,jpdoc),tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jpdoc)-tra(iii,jjj,kkk,jp13doc),tra(iii,jjj,kkk,jp13doc)/tra(iii,jjj,kkk,jpdoc)
        print*, "POC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jppoc),trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jppoc)-trb(iii,jjj,kkk,jp13poc),trb(iii,jjj,kkk,jp13poc)/trb(iii,jjj,kkk,jppoc)
        print*, trn(iii,jjj,kkk,jppoc),trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jppoc)-trn(iii,jjj,kkk,jp13poc),trn(iii,jjj,kkk,jp13poc)/trn(iii,jjj,kkk,jppoc)
        print*, tra(iii,jjj,kkk,jppoc),tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jppoc)-tra(iii,jjj,kkk,jp13poc),tra(iii,jjj,kkk,jp13poc)/tra(iii,jjj,kkk,jppoc)
        print*, "GOC "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpgoc),trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jpgoc)-trb(iii,jjj,kkk,jp13goc),trb(iii,jjj,kkk,jp13goc)/trb(iii,jjj,kkk,jpgoc)
        print*, trn(iii,jjj,kkk,jpgoc),trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jpgoc)-trn(iii,jjj,kkk,jp13goc),trn(iii,jjj,kkk,jp13goc)/trn(iii,jjj,kkk,jpgoc)
        print*, tra(iii,jjj,kkk,jpgoc),tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jpgoc)-tra(iii,jjj,kkk,jp13goc),tra(iii,jjj,kkk,jp13goc)/tra(iii,jjj,kkk,jpgoc)
        print*, "PHY "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpphy),trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jpphy)-trb(iii,jjj,kkk,jp13phy),trb(iii,jjj,kkk,jp13phy)/trb(iii,jjj,kkk,jpphy)
        print*, trn(iii,jjj,kkk,jpphy),trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jpphy)-trn(iii,jjj,kkk,jp13phy),trn(iii,jjj,kkk,jp13phy)/trn(iii,jjj,kkk,jpphy)
        print*, tra(iii,jjj,kkk,jpphy),tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jpphy)-tra(iii,jjj,kkk,jp13phy),tra(iii,jjj,kkk,jp13phy)/tra(iii,jjj,kkk,jpphy)
        print*, "PHY2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpdia),trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jpdia)-trb(iii,jjj,kkk,jp13dia),trb(iii,jjj,kkk,jp13dia)/trb(iii,jjj,kkk,jpdia)
        print*, trn(iii,jjj,kkk,jpdia),trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jpdia)-trn(iii,jjj,kkk,jp13dia),trn(iii,jjj,kkk,jp13dia)/trn(iii,jjj,kkk,jpdia)
        print*, tra(iii,jjj,kkk,jpdia),tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jpdia)-tra(iii,jjj,kkk,jp13dia),tra(iii,jjj,kkk,jp13dia)/tra(iii,jjj,kkk,jpdia)
        print*, "ZOO "                                                                                                                                      
        print*, trb(iii,jjj,kkk,jpzoo),trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jpzoo)-trb(iii,jjj,kkk,jp13zoo),trb(iii,jjj,kkk,jp13zoo)/trb(iii,jjj,kkk,jpzoo)
        print*, trn(iii,jjj,kkk,jpzoo),trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jpzoo)-trn(iii,jjj,kkk,jp13zoo),trn(iii,jjj,kkk,jp13zoo)/trn(iii,jjj,kkk,jpzoo)
        print*, tra(iii,jjj,kkk,jpzoo),tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jpzoo)-tra(iii,jjj,kkk,jp13zoo),tra(iii,jjj,kkk,jp13zoo)/tra(iii,jjj,kkk,jpzoo)
        print*, "ZOO2 "                                                                                                                                     
        print*, trb(iii,jjj,kkk,jpmes),trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jpmes)-trb(iii,jjj,kkk,jp13mes),trb(iii,jjj,kkk,jp13mes)/trb(iii,jjj,kkk,jpmes)
        print*, trn(iii,jjj,kkk,jpmes),trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jpmes)-trn(iii,jjj,kkk,jp13mes),trn(iii,jjj,kkk,jp13mes)/trn(iii,jjj,kkk,jpmes)
        print*, tra(iii,jjj,kkk,jpmes),tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jpmes)-tra(iii,jjj,kkk,jp13mes),tra(iii,jjj,kkk,jp13mes)/tra(iii,jjj,kkk,jpmes)
        print*, " "
      ENDIF
      !
      IF( ln_ligand )  &
      & CALL p4z_ligand( kt, knt )
      !                                                             !
      IF(ln_ctl)   THEN  ! print mean trends (used for debugging)
         WRITE(charout, FMT="('bio ')")
         CALL prt_ctl_trc_info(charout)
         CALL prt_ctl_trc(tab4d=tra, mask=tmask, clinfo=ctrcnm)
      ENDIF
      !
      IF( ln_timing )   CALL timing_stop('p4z_bio')
      !
   END SUBROUTINE p4z_bio

   !!======================================================================
END MODULE p4zbio
