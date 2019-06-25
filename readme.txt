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
~		p4zsed.F90  -->  sedimentary processes (NH4 release, denitrification, burial)
~		p4zpoc.F90  -->  remineralisation of organic particles


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
	


