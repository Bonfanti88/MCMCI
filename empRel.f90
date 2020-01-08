MODULE empRel
	
	!1. Johnson (1966): B-V logTeff relation
	!2. Jordi+ (2010): G_BP-G_RP logTeff relation
	DOUBLE PRECISION, DIMENSION(4,2), PARAMETER :: cJ=reshape((/ 0.D0, 0.D0, -0.234D0, 3.908D0, &
													&	-0.316D0, 0.709D0, -0.654D0, 3.999D0 /), &
													& (/4,2/))
	DOUBLE PRECISION :: DTeJ=0.02 !error bar amplitude in both relations
	
	!! Mamayek (2008) !!
	DOUBLE PRECISION, PARAMETER :: logtSunPD=9.65
	DOUBLE PRECISION, PARAMETER :: logRHKSun=-4.9126
	DOUBLE PRECISION, PARAMETER :: logtSunM=-38.053-17.912*logRHKSun-1.6675*logRHKSun**2 !Sun age according to Mamayek relation
	DOUBLE PRECISION, PARAMETER :: c0=-(logtSunM-logtSunPD)-38.053 !correction | Mamayek Relation finds Padova Sun_Age
	DOUBLE PRECISION, PARAMETER :: DlR=0.17 !Conservative level on magnetic activity
	DOUBLE PRECISION, PARAMETER :: logtCutoff5=8.75 !logt=8.7=>t=0.5Gyr (so to cut up to 8.7 included).
					 								!Max cut that may derive from magnetic activity
	!! Denissenkov (2010) !!
	DOUBLE PRECISION, PARAMETER :: OmSun=2.86e-6 !rad/s. Assume PSun=25.4d
	DOUBLE PRECISION, PARAMETER :: OmeZAMS=4.7*OmSun
	DOUBLE PRECISION, PARAMETER :: A=(OmeZAMS/OmSun)**2-1
	DOUBLE PRECISION, PARAMETER :: B=0 !!!!In fact OmeZAMS<=OmSat
	DOUBLE PRECISION, PARAMETER :: logtZAMS=7.6
	DOUBLE PRECISION, PARAMETER :: tZAMS=10.**logtZAMS
	
	!! Barnes (2010) !!
	DOUBLE PRECISION, PARAMETER :: kC=0.646 ![d/Myr]
	DOUBLE PRECISION, PARAMETER :: kI=452 ![Myr/d]
	DOUBLE PRECISION, PARAMETER :: P0=1.1 ![d]
	DOUBLE PRECISION, PARAMETER :: logtCutoff25=9.45 !logt=9.4=>t=2.5Gyr (so to cut up to 9.4 included). Max cut that may derive
													 !from rotation, as also shown on Nature (29/01/15vol517,589)
	DOUBLE PRECISION, PARAMETER :: DtBarnes=0.125 !125Myr age uncertainty given by Barnes
	
	!! Nissen (2016) [Y/Mg]=a+b*t !!
	DOUBLE PRECISION, PARAMETER :: aYMg=0.170;
	DOUBLE PRECISION, PARAMETER :: bYMg=-0.0371;
	DOUBLE PRECISION, PARAMETER :: DtNissen=0.6e9; ![yrs] Nissen (2016) uncertainty on age
	
	!! Torres (2009) !!
	DOUBLE PRECISION, DIMENSION(7), PARAMETER :: bT=(/ 2.4427, 0.6679, 0.1771, 0.705, -0.21415, 0.02306, 0.04173 /)
	!! My coefficients as inferred from least squared interpolation from Torres sample
	!R is given in solar radii
	!input logg or logrho must be in cgs 
	!R=R(logTeff,logg,[Fe/H]) 
	DOUBLE PRECISION, DIMENSION(10), PARAMETER :: blogg=(/ 4.588468, 0.698017, 0.290507, 0.250746, -2.037976, &
												& 0.400437, -0.036721, 0.031149, -0.080481, 0.013993 /)
	!R=R(logTeff,logrho,[Fe/H]) 
	DOUBLE PRECISION, DIMENSION(10), PARAMETER :: brho=(/ 0.1922952, 0.4429884, 0.1975146, 0.0951324, -0.3927193, &
												& -0.0190245, -0.0052578, 0.0206659, -0.0544062, 0.0014269 /)
	
	!! Eker+ (2015) !! logM=aE*logL+bE
	DOUBLE PRECISION, DIMENSION(3), PARAMETER :: logLsup=(/ 7.6577D-2, 1.6436D0, 3.4683D0 /)
	DOUBLE PRECISION, DIMENSION(4), PARAMETER :: aE=(/ 2.0657D-1, 2.3105D-1, 2.5240D-1, 3.6684D-1 /)
	DOUBLE PRECISION, DIMENSION(4), PARAMETER :: bE=(/ 5.3708D-3, 4.6211D-4, -3.0288D-2, -4.5378D-1 /)

END MODULE empRel
