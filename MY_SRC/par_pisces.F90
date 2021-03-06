MODULE par_pisces
   !!======================================================================
   !!                        ***  par_pisces  ***
   !! TOP :   set the PISCES parameters
   !!======================================================================
   !! History :   2.0  !  2007-12  (C. Ethe, G. Madec)  revised architecture
   !!----------------------------------------------------------------------

   IMPLICIT NONE

   ! productive layer depth
   INTEGER, PUBLIC ::   jpkb       !: first vertical layers where biology is active
   INTEGER, PUBLIC ::   jpkbm1     !: first vertical layers where biology is active

   ! assign an index in trc arrays for each LOBSTER prognostic variables
   INTEGER, PUBLIC ::   jpdet     !: detritus                   
   INTEGER, PUBLIC ::   jpdom     !: dissolved organic matter 
   INTEGER, PUBLIC ::   jpdic     !: dissolved inoganic carbon concentration 
   INTEGER, PUBLIC ::   jptal     !: total alkalinity 
   INTEGER, PUBLIC ::   jpoxy     !: oxygen carbon concentration 
   INTEGER, PUBLIC ::   jpao2     !: oxygen carbon concentration 
   INTEGER, PUBLIC ::   jpcal     !: calcite  concentration 
   INTEGER, PUBLIC ::   jppo4     !: phosphate concentration 
   INTEGER, PUBLIC ::   jppoc     !: small particulate organic phosphate concentration
   INTEGER, PUBLIC ::   jpsil     !: silicate concentration
   INTEGER, PUBLIC ::   jpphy     !: phytoplancton concentration 
   INTEGER, PUBLIC ::   jpzoo     !: zooplancton concentration
   INTEGER, PUBLIC ::   jpdoc     !: dissolved organic carbon concentration 
   INTEGER, PUBLIC ::   jpdia     !: Diatoms Concentration
   INTEGER, PUBLIC ::   jpmes     !: Mesozooplankton Concentration
   INTEGER, PUBLIC ::   jpdsi     !: Diatoms Silicate Concentration
   INTEGER, PUBLIC ::   jpfer     !: Iron Concentration
   INTEGER, PUBLIC ::   jpbfe     !: Big iron particles Concentration
   INTEGER, PUBLIC ::   jpgoc     !: big particulate organic phosphate concentration
   INTEGER, PUBLIC ::   jpsfe     !: Small iron particles Concentration
   INTEGER, PUBLIC ::   jpdfe     !: Diatoms iron Concentration
   INTEGER, PUBLIC ::   jpgsi     !: (big) Silicate Concentration
   INTEGER, PUBLIC ::   jpnfe     !: Nano iron Concentration
   INTEGER, PUBLIC ::   jpnch     !: Nano Chlorophyll Concentration
   INTEGER, PUBLIC ::   jpdch     !: Diatoms Chlorophyll Concentration
   INTEGER, PUBLIC ::   jpno3     !: Nitrates Concentration
   INTEGER, PUBLIC ::   jpnh4     !: Ammonium Concentration
   INTEGER, PUBLIC ::   jpdon     !: dissolved organic nitrogen concentration
   INTEGER, PUBLIC ::   jpdop     !: dissolved organic phosphorus concentration
   INTEGER, PUBLIC ::   jppon     !: small particulate organic nitrogen concentration
   INTEGER, PUBLIC ::   jppop     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpnph     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppph     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpndi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppdi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppic     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpnpi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpppi     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppfe     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jppch     !: small particulate organic phosphorus concentration
   INTEGER, PUBLIC ::   jpgon     !: Big nitrogen particles Concentration
   INTEGER, PUBLIC ::   jpgop     !: Big phosphorus particles Concentration
   INTEGER, PUBLIC ::   jplgw     !: Weak Ligands

   INTEGER, PUBLIC ::   jp15poc   !: 15N small particulate organic concentration
   INTEGER, PUBLIC ::   jp15phy   !: 15N phytoplankton concentration
   INTEGER, PUBLIC ::   jp15zoo   !: 15N zooplankton concentration
   INTEGER, PUBLIC ::   jp15doc   !: 15N dissolved organic matter concentration
   INTEGER, PUBLIC ::   jp15dia   !: 15N diatoms concentration
   INTEGER, PUBLIC ::   jp15mes   !: 15N mesozooplankton concentration
   INTEGER, PUBLIC ::   jp15goc   !: 15N big particulate organic concentration
   INTEGER, PUBLIC ::   jp15no3   !: 15N Nitrates Concentration
   INTEGER, PUBLIC ::   jp15nh4   !: 15N Ammonium Concentration
   
   INTEGER, PUBLIC ::   jp13poc   !: 13C small particulate organic concentration
   INTEGER, PUBLIC ::   jp13phy   !: 13C phytoplankton concentration
   INTEGER, PUBLIC ::   jp13zoo   !: 13C zooplankton concentration
   INTEGER, PUBLIC ::   jp13doc   !: 13C dissolved organic matter concentration
   INTEGER, PUBLIC ::   jp13dia   !: 13C diatoms concentration
   INTEGER, PUBLIC ::   jp13mes   !: 13C mesozooplankton concentration
   INTEGER, PUBLIC ::   jp13goc   !: 13C big particulate organic concentration
   INTEGER, PUBLIC ::   jp13dic   !: 13C dissolved inorganic Concentration
   INTEGER, PUBLIC ::   jp13cal   !: 13C calcite Concentration
   

   !!---------------------------------------------------------------------
   !!   Default                                   No CFC geochemical model
   ! Starting/ending PISCES do-loop indices (N.B. no PISCES : jpl_pcs < jpf_pcs the do-loop are never done)
   INTEGER, PUBLIC  ::   jp_pcs0  !: First index of PISCES tracers
   INTEGER, PUBLIC  ::   jp_pcs1  !: Last  index of PISCES tracers

   !!----------------------------------------------------------------------
   !! NEMO/TOP 4.0 , NEMO Consortium (2018)
   !! $Id: par_pisces.F90 10416 2018-12-19 11:45:43Z aumont $ 
   !! Software governed by the CeCILL license (see ./LICENSE)
   !!======================================================================
END MODULE par_pisces
