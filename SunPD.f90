MODULE SunPD
	
	DOUBLE PRECISION, PARAMETER :: TeffSun=5778. !K
	DOUBLE PRECISION, PARAMETER :: loggSun=4.43191399643523D0 !cgs
	DOUBLE PRECISION, PARAMETER :: RSun=6.9598D10 !cm
	DOUBLE PRECISION, PARAMETER :: RSunm=6.9598D8 !m
	DOUBLE PRECISION, PARAMETER :: RSunKm=6.9598D5 !km
	DOUBLE PRECISION, PARAMETER :: Grav=6.67384D-8 !cgs
	DOUBLE PRECISION, PARAMETER :: MSun=10.D0**loggSun*RSun**2/Grav !g !!! MSun=1.96214783383218D+33
	DOUBLE PRECISION, PARAMETER :: MSunkg=MSun/1000.D0 !1.96214783383218D+30 !kg
	DOUBLE PRECISION, PARAMETER, PRIVATE :: pi=3.14159265358979 !Greek pi
	DOUBLE PRECISION, PARAMETER :: rhoSun=MSun/(4./3.*pi*RSun**3) !g/cm3   1.38948171766240D0 !g/cm3   !!
	!!Absolute magnitude
	!(1) Johnson V-band interpolated in the grid as from Bonfanti+ (2015)
	!(2) Gaia G-band as taken from Casagrande&VandenBerg (2018)
	DOUBLE PRECISION, DIMENSION(2), PARAMETER :: MabsSun=(/4.833,4.67/)
	
	DOUBLE PRECISION, PARAMETER :: ZSun=0.01524D0 !present day solar metallicity
                           			!(Caffau+, 2011). The initial one is 0.01774
    DOUBLE PRECISION, PARAMETER  :: costZ=-log10(ZSun) !1.81701503299642D0 !!!Z=10.**(FeH-1.817)!!!
    !!
    DOUBLE PRECISION, PARAMETER  :: dYdZ=1.78,He0=0.2485
    DOUBLE PRECISION, PARAMETER  :: numaxSun=3090 !microHz Chaplin 14. In Kjeldsen&Bedding 3050
	DOUBLE PRECISION, PARAMETER  :: DnuSun=135.1 !microHz Chaplin 14. In Kjeldsen&Bedding 134.9
                        	
END MODULE SunPD
