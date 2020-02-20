MODULE IsoPD
	
	INTEGER, PARAMETER :: ct=2,cM=4,clogL=5,clogTe=6,clogg=7
	INTEGER, PARAMETER :: cmag1=26,cmag2=25,cmbol=23,cmagx=24 !V->26; B->25
	INTEGER, PARAMETER :: logtVal=0
	INTEGER, PARAMETER :: cgZAMS=5,crhoZAMS=6,cLZAMS=3
	INTEGER, PARAMETER :: ct_T=3,cM_T=2,clogL_T=4 !ct_T yrs, not log
	INTEGER, PARAMETER :: clogTe_T=5,clogR_T=6
	INTEGER, PARAMETER :: cVelL=2,cVelg=3,cVelrho=4
	INTEGER, PARAMETER :: cBTeff=1,cBBmV=2,ctau=4
	INTEGER, PARAMETER :: clab=8
		
	INTEGER, PRIVATE :: idd
    DOUBLE PRECISION, DIMENSION(38) :: MLeo=(/ (dble(idd)/100, idd=55,240,5) /)
    DOUBLE PRECISION, PARAMETER :: MinfED=0.55,MsupED=2.4 !element diffusion
    DOUBLE PRECISION, PARAMETER :: MminT=0.09,MmaxT=350.,MlowMS=0.5 !for lower mass value R=R(t) is considered constant
    DOUBLE PRECISION, PARAMETER :: logt_halfMS=9.7 !t=5e9 yrs
    DOUBLE PRECISION, PARAMETER :: logtstep=0.05
    
    DOUBLE PRECISION, PARAMETER :: BmVinf=-0.3,BmVsup=2.5
    DOUBLE PRECISION, PARAMETER :: Zi_inf=5.E-4,Zi_sup=0.07
        
    ! Low limits to avoid numerical errors during interpolation !
    DOUBLE PRECISION, PARAMETER :: dlim=5.D-4 !minimum horizontal and vertical distance that is allowed between reference isoc.
									!If lower than this limit, other ref isoc are considered
	DOUBLE PRECISION, PARAMETER :: dCVvlim=1.D-4,dCLlim=1.D-4,dCglim=1.D-4
	DOUBLE PRECISION, PARAMETER :: varTrlim=1.D-8 !minimum variation requested in Track parameters
	
	DOUBLE PRECISION, PARAMETER :: logtlimCal=7.,DlogL=0.001,DMTr=0.05
	DOUBLE PRECISION, PARAMETER :: vlab=9.,tlab=1.D-4
	
	LOGICAL, PARAMETER :: photIsocAvail=.true.,elDiff=.true.
	
END MODULE IsoPD
