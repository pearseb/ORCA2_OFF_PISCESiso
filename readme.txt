Files to be altered:

   routines within /users/pearseb/NEMOv4.0/src/TOP/PISCES/
~	par_pisces.F90  -->  sets global parameters for PISCES and defines tracers [module] 
~	trcini_pisces.F90  -->  initialisation of the PISCES biogeochemical model 
~	trcwri_pisces.F90  -->  defines outputs of tracers

   routines within /users/pearseb/NEMOv4.0/src/TOP/PISCES/P4Z/
	p4zsms.F90  -->  calls biological routines to calculate PISCES tracers
	    CALLS:
~		p4zsbc.F90  -->  boundary conditions (atmospheric and riverine fluxes)
~		p4zbio.F90  -->  controls the interactions between compartments of PISCES
	   	    CALLS:
~			p4zsink.F90  -->  vertical flux of particulate organic matter
                            CALLS:
~                               trcsink.F90 --> time splitting routine for sinking flux within/across depth levels
~			p4zprod.F90  -->  growth rate of the 2 phytoplankton
~			p4zmort.F90  -->  mortality of the 2 phytoplankton
~			p4zmicro.F90  -->  sources and sinks of microzooplankton
~			p4zmeso.F90  -->  sources and sinks of mesozooplankton
~                       p4zagg.F90  -->  Aggregation of particles
~			p4zrem.F90  -->  remineralisation of organic matter
~			p4zpoc.F90  -->  remineralisation of organic particles
~		p4zsed.F90  -->  sedimentary processes (NH4 release, denitrification, burial)


STEPS TAKEN...
	1. run model with additional dummy tracers of jp15no3 and jp15nh4
        2. added namelist (logical) control of n15 routines
        3. added full list of delta15N tracers without including fractionation
		- jp15no3	Nitrate
		- jp15nh4	Ammonium
		- jp15phy	Nanophytoplankton
		- jp15zoo	Microzooplankton
		- jp15dia	Diatoms
		- jp15mes	Mesozooplankton	
		- jp15poc	Small particulate organic carbon
		- jp15goc	Large particulate organic carbon
		- jp15doc	Dissolved organic carbon
        4. added d15n signatures of external sources (rivers + atmospheric deposition + n2fixation)
        5. added fractionation due to utilisation 
		- p4zprod.F90
        6. added fractionation due to water column denitrification and nitrification
		- p4zrem.F90
        7. added fractionation due to benthic denitrification and included fluxes at sediment
		- p4zsed.F90
        8. modified existing additions to avoid anomolous fractionation at small concentrations
		- BEFORE: trb(ji,jj,jk,jp15xxx) / ( trb(ji,jj,j,jpxxx) + rtrn )
		- AFTER: ( trb(ji,jj,jk,jp15xxx) + rtrn ) / ( trb(ji,jj,j,jpxxx) + rtrn )
        9. carried d15n signatures between DOC, POC and GOC during aggregation and disaggregation
		- p4zagg.F90
		- p4zpoc.F90
        10. carried d15n signatures between PHY, PHY2, POC and GOC during phytoplankton mortality
		- p4zmort.F90
        11. carried d15n signatures between NH4, DOC, PHY, PHY2, ZOO, POC and GOC in microzooplankton lifecycle
		- p4zmicro.F90
        12. carried d15n signatures between NH4, DOC, PHY, PHY2, ZOO, ZOO2, POC and GOC in mesozooplankton lifecycle
		- p4zmeso.F90
        13. added fractionation during zooplankton feeding
		- p4zmicro.F90
		- p4zmeso.F90
     
	


