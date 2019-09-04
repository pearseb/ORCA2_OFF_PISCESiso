Files to be altered for N15 routines:

   routines within /users/pearseb/NEMOv4.0/src/TOP/PISCES/
~	par_pisces.F90  -->  sets global parameters for PISCES and defines tracers [module] 
~	trcini_pisces.F90  -->  initialisation of the PISCES biogeochemical model 
~	trcwri_pisces.F90  -->  defines outputs of tracers
~	sms_pisces.F90  -->  defines flags (ln_n15) and min/max of tracer (nn_n15min/nn_n15max)

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
        14. placed constraints on min/max d15n of NH4 (set by nn_n15min & nn_n15max in namelist_pisces_ref)
		- p4zsms.F90
        15. removed contraints on min/max d15n for NH4, and
            added rtrn to numerator of calculations of variables "zu_15" and "zun_15" in p4zprod.F90
                - p4zsms.F90
                - p4zprod.F90 (lines 366-370)
        16. bug found and corrected (lines 255-257), ratio of 15N:N missing for mortality of mesozooplankton
                - p4zmeso.F90
        17. bug found and corrected (line 413), variable zr15_no3 used incorrect k-level index (jk --> ikt)
                - p4zsed.F90

     
	


Files to be altered for C13 routines:

   routines within /users/pearseb/NEMOv4.0/src/TOP/PISCES/
~	par_pisces.F90  -->  sets global parameters for PISCES and defines tracers [module] 
~	trcini_pisces.F90  -->  initialisation of the PISCES biogeochemical model 
~	sms_pisces.F90  -->  defines flags (ln_c13)

   routines within /users/pearseb/NEMOv4.0/src/TOP/PISCES/P4Z/
	p4zsms.F90  -->  calls biological routines to calculate PISCES tracers
	    CALLS:
		p4zsbc.F90  -->  boundary conditions (atmospheric and riverine fluxes)
		p4zbio.F90  -->  controls the interactions between compartments of PISCES
	   	    CALLS:
~			p4zsink.F90  -->  vertical flux of particulate organic matter
                            CALLS:
                                trcsink.F90 --> time splitting routine for sinking flux within/across depth levels
~			p4zprod.F90  -->  growth rate of the 2 phytoplankton
~			p4zmort.F90  -->  mortality of the 2 phytoplankton
~			p4zmicro.F90  -->  sources and sinks of microzooplankton
~			p4zmeso.F90  -->  sources and sinks of mesozooplankton
~                       p4zagg.F90  -->  Aggregation of particles
~			p4zrem.F90  -->  remineralisation of organic matter
~			p4zpoc.F90  -->  remineralisation of organic particles
~		p4zsed.F90  -->  sedimentary processes (NH4 release, denitrification, burial)
~		p4zlys.F90  -->  Compute CaCO3 saturation
~		p4zflx.F90  -->  Compute surface fluxes (gas exchange)


STEPS TAKEN...
~	1. run model with dummy tracers and added namelist (logical) control of c13 routines:
		- jp13dic	Dissolved inorganic carbon
		- jp13cal	Calcite
		- jp13phy	Nanophytoplankton
		- jp13zoo	Microzooplankton
		- jp13dia	Diatoms
		- jp13mes	Mesozooplankton	
		- jp13poc	Small particulate organic carbon
		- jp13goc	Large particulate organic carbon
		- jp13doc	Dissolved organic carbon
~	2. added namelist variables for controlling d13c and e13c values
		- p4zflx.F90 	(d13c_co2)
		- p4zsbc.F90 	(d13c_rivdoc, d13c_rivdic, d13c_fix)
		- p4zsed.F90	(d13c_rivdoc, d13c_rivdic, d13c_fix)
		- p4zprod.F90	(e13c_min, e13c_max)
		- p4zmort.F90	(e13c_cal)
		- p4zmicro.F90	(e13c_calz)
		- p4zmeso.F90	(e13c_cal2)
~	3. ensured that riverine fluxes do not bias d13c variables
	4. added fractionation factors to inner carbon cycling
~		- p4zprod.F90
~		- p4zmort.F90
~		- p4zmicro.F90
~		- p4zmeso.F90
	5. added fractionation during air-sea gas exchange (Zhang et al. 1995)
		- p4zflx.F90



