SUBROUTINE SCPmcmcPD12S(SCP,Intestaz,IsocTab,Zvec,Zndxi,Zndxf,TrackTab,nM,Mav,indxZt,Tndxi,Tndxf, &
			& velTrackTab,Vndxi,Vndxf,GyroTab,ZAMStab,ZAndxi,ZAndxf,ZtevoTab,Endxi,Endxf,MeAv,M_star,I_M_star, &
			& R_star,I_R_star,Teff,I_Teff,t_star,I_t_star,L,I_L,link,row,acc)
	
	USE SunPD
	USE IsoPD
	USE empRel
	
	IMPLICIT NONE
	CHARACTER(LEN=171) :: Intestaz
	
	INTEGER :: link,row,acc
	INTEGER, DIMENSION(:), INTENT(in) :: Zndxi,Zndxf,nM,indxZt,Vndxi,Vndxf
	INTEGER, DIMENSION(:), INTENT(in) :: ZAndxi,ZAndxf,Endxi,Endxf
	INTEGER, DIMENSION(:,:), INTENT(in) :: Tndxi,Tndxf
	
	DOUBLE PRECISION M_star,R_star,I_M_star,I_R_star,Teff,I_Teff,t_star,I_t_star,L,I_L
	DOUBLE PRECISION, DIMENSION(29) :: SCP
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: Zvec,MeAv
	DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: IsocTab,TrackTab,Mav,velTrackTab,GyroTab
	DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: ZAMStab,ZtevoTab
	!!!
	CHARACTER(LEN=8)  :: Isoch="./Isoch/"
	CHARACTER(LEN=4)  :: FileTxt="Star"
	
	CHARACTER(LEN=4)  :: NamesFile="Name"
	CHARACTER(LEN=20) :: sdum !equal to length of fName
	CHARACTER(LEN=1000) :: titolo,titoloN,titoloAge,titoloAgeN
			
	CHARACTER(LEN=:), ALLOCATABLE :: model,percorso,Est
	CHARACTER(LEN=:), ALLOCATABLE :: outNames
	CHARACTER(LEN=:), ALLOCATABLE :: TLRM,MetCut
	CHARACTER(LEN=:), ALLOCATABLE :: outAge,outAgeN,outI,outIN
	
	INTEGER :: IncZ=0,indxZ,indxM,xMTi,xMTf
	INTEGER :: Bin,cycleW
	INTEGER :: cRif,cBrif
	INTEGER :: ios,lines
	INTEGER :: rowN,row_i,i
	INTEGER :: sismo,idCol,cmag0
	INTEGER :: rowAge,Nage,ndxDU,age_s,age_s1,jd,kD
	INTEGER :: ti0,tf0,statev,state,rowDv,rowD
	INTEGER :: check
	INTEGER :: mid_n,ndxTot,jrho,ndxrho,klim,klim0,resto,c_i
	INTEGER :: logtCutoff0,iIso,iIso8
	INTEGER :: cZAMS,cyT
	INTEGER :: ageIsoc,kTr,Mndx,iZg,indxZgini,kwhile
	INTEGER :: gFlag,checkM,Zf,Zi,xVTi,xVTf,xZAi,xZAf
	INTEGER :: jj,jk,kwl
	
	INTEGER, DIMENSION(2) :: ti,tf,id
	INTEGER, DIMENSION(:), ALLOCATABLE :: dist_ndx,ix,dist_ndxU,ix_sorted
	INTEGER, DIMENSION(:), ALLOCATABLE :: ndxrhoi,ndxrhof
	INTEGER, DIMENSION(:), ALLOCATABLE :: ndxDv,ndxD,Iso_indx
	INTEGER, DIMENSION(:), ALLOCATABLE :: indxZg,cumInt
		
	DOUBLE PRECISION, PARAMETER :: e=2.71828182845905 !Nepero number
	DOUBLE PRECISION, PARAMETER :: pi=3.14159265358979 !Greek pi
			
	DOUBLE PRECISION :: Z_isoLOW,Z_isoUP
	DOUBLE PRECISION :: YMg,t_Nissen,logt_Nissen
	DOUBLE PRECISION :: logL,R,I_R,M,I_M
	DOUBLE PRECISION :: V,BmV,Teffinput,Teffinput0,Hipf,d,I_d,DM0,Vass,VMag
	DOUBLE PRECISION :: numax,I_numax,Dnu,I_Dnu,Rf,R1,I_R1,R2,I_R2,Ri,I_Ri,Rc,I_Rc
	DOUBLE PRECISION :: rhof,rho,I_rho,logrho,I_logrho,loggf,logg,I_logg,g,I_g
	DOUBLE PRECISION :: rhoS,I_rhoS,logrhoS,I_logrhoS,loggS,I_loggS,gS,I_gS
	DOUBLE PRECISION :: rhoAS,I_rhoAS,logrhoAS,I_logrhoAS,loggAS,I_loggAS,gAS,I_gAS
	DOUBLE PRECISION :: logTeff,I_logTeff,I_logL
	DOUBLE PRECISION :: Z,Zact,Yact,Yini
	DOUBLE PRECISION :: first_logt,last_logt
	DOUBLE PRECISION :: Vf,BmVf,rhoI,loggI,hstar,hstarlim,logTeJ
	DOUBLE PRECISION :: logRHK,vsini,Uvsini,P,UP
	DOUBLE PRECISION :: x,x_i1,x_i2,x_i,dist0,xtau
	DOUBLE PRECISION :: BmV_i1,V_i1,logTeff_i1,logL_i1,M_i1,logg_i1,logrho_i1
	DOUBLE PRECISION :: BmV_i2,V_i2,logTeff_i2,logL_i2,M_i2,logg_i2,logrho_i2
	DOUBLE PRECISION :: mBV,mT,mg,mL,mM,mrh,qBV,qg,qL,qM,qrho,qT
	DOUBLE PRECISION :: age_iso1,age_iso2,tmp_age1,agePivot
	DOUBLE PRECISION :: BC
	DOUBLE PRECISION :: dCV_isv,dCL_isv,dCL_is,dHR_is
	DOUBLE PRECISION :: loggf0,logg0,I_logg0,g0,I_g0,rhof0,rho0,I_rho0,logrho0,I_logrho0
	DOUBLE PRECISION :: x_is2,y1_is2,y2_is2
	DOUBLE PRECISION :: Mlogg,I_Mlogg,Mrho,I_Mrho
	DOUBLE PRECISION :: y,yl,y_lim,dClim
	DOUBLE PRECISION :: logtlim,logtMmod,logtlim1,logtlim2,logtlim3,logtlim4
	DOUBLE PRECISION :: logY,logY_soglia,logY2,logY2_soglia
	DOUBLE PRECISION :: Rv,MilowMS,logY_sogliaVsini,xlow,xup
	DOUBLE PRECISION :: R_Tr,rho_Tr,XT
	DOUBLE PRECISION :: Omega,tau,logt_Barnes,rhot,logtCutoff,logt8,xx,yy,Zlow,Zup
	DOUBLE PRECISION :: logY0,logY0_soglia
	DOUBLE PRECISION :: vRif
	DOUBLE PRECISION :: Spesi,M_starTr,I_y,M_Tracks,M_starTr2
	DOUBLE PRECISION :: Teff_star,L_star,g_star,rho_star,BmV_star,BC_star,logg_star
	DOUBLE PRECISION :: I_Teff_star,I_L_star,I_g_star,I_rho_star,I_BmV_star,I_BC_star
	DOUBLE PRECISION :: M_starNC,I_M_starNC,R_starNC,I_R_starNC
	DOUBLE PRECISION :: I_logg_star,Zg1,Zg2,rZ,Zgini1,Zgini2,Zgini,Z2,FeH
	DOUBLE PRECISION :: t_star2,Teff_star2,L_star2,M_star2,g_star2,rho_star2,BmV_star2,BC_star2,logg_star2
	DOUBLE PRECISION :: R_star2,I_t_star2,I_Teff_star2,I_L_star2,I_M_star2,I_g_star2,I_rho_star2,I_BmV_star2
	DOUBLE PRECISION :: I_BC_star2,I_logg_star2,I_R_star2,tol
	DOUBLE PRECISION :: FeH_,I_FeH_,Z_,I_Z_,useColor,InclogTeff
	
	DOUBLE PRECISION, DIMENSION(10) :: vecT
	DOUBLE PRECISION, DIMENSION(4) :: limAges
	DOUBLE PRECISION, DIMENSION(size(SCP)+12) :: SCPInputtmp,SCPInputNtmp
	DOUBLE PRECISION, DIMENSION(25) :: Stelle,StelleN
	DOUBLE PRECISION, DIMENSION(2) :: BmV_is1v,VMag_is1v,logL_is1v,VMag_isv,logL_isv
	DOUBLE PRECISION, DIMENSION(2) :: BmV_is2v,VMag_is2v,logL_is2v
	DOUBLE PRECISION, DIMENSION(2) :: logL_is1,BmV_is1,logT_is1,BmV_is,logT_is,logy_is1
	DOUBLE PRECISION, DIMENSION(2) :: logL_is2,BmV_is2,logT_is2,logy_is2
	DOUBLE PRECISION, DIMENSION(1) :: logR_Tr,logY_sogliaVsini1,ylow,yup
	DOUBLE PRECISION, DIMENSION(1) :: tau1,taulow,tauup
	DOUBLE PRECISION, DIMENSION(4) :: t_starvec
	
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: logTeTr,logTeTrtmp
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: t_iso,logt_iso
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: M_iso,logL_iso,L_iso
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: logTeff_iso,Teff_iso
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: logg_iso,g_iso,R_iso
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rho_iso,logrho_iso
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: BmV_iso,V_iso,BC_iso
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VMag_iso,ageMinDist
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: MTr,RTr
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: loggTr,logRhoTr
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: t_iTmp,BC_iTmp,logTeff_iTmp,BmV_iTmp
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: logg_iTmp,logL_iTmp,M_iTmp,logrho_iTmp
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: t_i,BC_i,BmV_i,logL_i
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: logTeff_i,M_i,logg_i
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: logrho_i,Zguess
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: x_iso,xx_iso
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Teff_i,L_i,g_i,rho_i,y_i
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dsorted
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Diff_BmV,Diff_logL
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: BmVD,VMagD
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: VMag_isor,BmV_isor,logL_isor,logT_isor
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: logy_iso,logy_isor,Diff_logy,y_iso
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: logY_isot,rho_isot
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: yT_iso,Gauss,w
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: t_Tracks,logTeff_Tracks,logL_Tracks,vEvo
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: R_Tracks,logY_Tracks
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Zini,tsp,tmax,MTrAv
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: SCPInput,SCPInputN,Mvec
	
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Isoc,dist,distTmp,TabVel,ZAMSz
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Track05,Iso_i,Iso_i2,Tracks,Ztevo
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: TrRhoG,Isoc_i,a1,a2
	
	LOGICAL :: loggAvailCal,rhoAvailCal,loggAvail0,rhoAvail0,chooseLogg
	LOGICAL :: calibHRD,calibNoD,calibSPEC,gProxyAvail,rhoAvail,loggAvail
	LOGICAL :: isEq,BmVuguali,logLuguali,logYuguali,caliblogL
	LOGICAL :: vcheck,Mavail,rhoAvailAS,loggAvailAS,bestLogg,agreeM
	
	LOGICAL, DIMENSION(:), ALLOCATABLE :: TeB
	
	interface
	
		subroutine loadMatrix(fileName,matrix,head)
			IMPLICIT NONE
			CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: fileName
			DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: matrix
			CHARACTER(LEN=1000) :: head
		end subroutine loadMatrix
		
		subroutine trovaIndici(v,ndxi,ndxf)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: v
			INTEGER, DIMENSION(:), ALLOCATABLE :: ndxi,ndxf
		end subroutine trovaIndici
		
		subroutine num2str(x,d,str)
			IMPLICIT NONE
			DOUBLE PRECISION x
			CHARACTER(LEN=:), ALLOCATABLE :: str !was 40
			INTEGER d
		end subroutine num2str
		
		subroutine openTracksPD(Z,M,PD,Tracks,M_in,Is)
			IMPLICIT NONE
			DOUBLE PRECISION Z,M,M_in
			INTEGER PD
			DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Tracks
			CHARACTER(LEN=:), ALLOCATABLE :: Is
		end subroutine
		subroutine indexx(arr,indx)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: arr
			INTEGER indx(size(arr))
 		end subroutine indexx
 		subroutine sort1(arr)
 			IMPLICIT NONE
 			DOUBLE PRECISION, DIMENSION(:), INTENT(inout) :: arr
 		end subroutine sort1
 		subroutine sort0(arr)
			IMPLICIT NONE
			INTEGER, DIMENSION(:), INTENT(inout) :: arr
		end subroutine sort0
 		
 		subroutine findYgivenX_i(M,x,xx,ny,y)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: M
			DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: y
			DOUBLE PRECISION x
			INTEGER xx,ny 
		end subroutine findYgivenX_i
		subroutine findYgivenX_v(M,x,xx,ny,y)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: M
			DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: y
			DOUBLE PRECISION, DIMENSION(size(M,1)) :: xx
			DOUBLE PRECISION x
			INTEGER ny
		end subroutine findYgivenX_v
		
		subroutine find(b,ix)
			IMPLICIT NONE
			LOGICAL, DIMENSION(:), INTENT(in) :: b
			INTEGER, DIMENSION(:), ALLOCATABLE :: ix
		end subroutine find
		
		subroutine uniqueFast(list,d,indices,first)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: list
			INTEGER d
			INTEGER, DIMENSION(:), ALLOCATABLE :: indices
			LOGICAL first
		end subroutine uniqueFast
		
		subroutine choose_i2Col(ndxDU,x,x_iso,xx_iso,calibSPEC,hstar,hstarlim,V_iso,logL_iso, &
			& M_iso,logg_iso,logrho_iso,logt_iso,BmV_i2,V_i2,logTeff_i2,logL_i2,M_i2,logg_i2,logrho_i2)
			IMPLICIT NONE
			INTEGER ndxDU
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: x_iso
			DOUBLE PRECISION, DIMENSION(size(x_iso)) :: xx_iso
			DOUBLE PRECISION, DIMENSION(size(x_iso)) :: V_iso,logL_iso
			DOUBLE PRECISION, DIMENSION(size(x_iso)) :: M_iso,logg_iso,logrho_iso,logt_iso
			DOUBLE PRECISION x,BmV_i2,V_i2,logTeff_i2,logL_i2,M_i2,logg_i2,logrho_i2
			DOUBLE PRECISION hstar,hstarlim
			LOGICAL calibSPEC
		end subroutine choose_i2Col
		subroutine choose_i2xy1y2(x,x_isor,y1_isor,y2_isor,i_jd,hstar,hstarlim, & !y1,
			& x_is1,x_is2,y1_is1,y1_is2,y2_is1,y2_is2) !_is1 _is2 two points on the same isoc
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: x_isor,y1_isor,y2_isor
			DOUBLE PRECISION x,hstar,hstarlim,x_is1,x_is2,y1_is1,y1_is2,y2_is1,y2_is2 !y1,
			INTEGER i_jd
		end subroutine choose_i2xy1y2
		
		subroutine selectIsoc(logt_iso,age,logtstep,last_logt,ti0,tf0)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: logt_iso
			DOUBLE PRECISION age,logtstep,last_logt
			INTEGER ti0,tf0
		end subroutine selectIsoc
		
		subroutine searchTOold(y_iso,logTeff_iso,caliblogL,logt_iso,last_logt,logY_soglia)
			IMPLICIT NONE
			DOUBLE PRECISION last_logt,logY_soglia
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: y_iso,logTeff_iso,logt_iso
			LOGICAL caliblogL
		end subroutine searchTOold
		
		subroutine InterpLin_M(fileM,x,colx,OutRange,colY,y,xlow,ylow,xup,yup)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: fileM
			DOUBLE PRECISION x,xlow,xup
			INTEGER, DIMENSION(:), INTENT(in) :: colY
			DOUBLE PRECISION, DIMENSION(size(colY)) :: ylow,yup,y
			INTEGER colx,OutRange
		end subroutine InterpLin_M
		
		subroutine M2R(logt_iso,mid_n,Isoc,Mi,caliblogL,loggAvail,loggAvailCal,logY,Rv)
			IMPLICIT NONE
			DOUBLE PRECISION Mi,logY,Rv
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: logt_iso
			DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: Isoc
			INTEGER mid_n
			LOGICAL caliblogL,loggAvail,loggAvailCal
		end subroutine M2R
		
		subroutine selectMfromTracks(MTrAv,M,I_M,Mvec)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: MTrAv
			DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Mvec
			DOUBLE PRECISION M,I_M
		end subroutine selectMfromTracks
		
		subroutine multipleZisocPhSCP(FeH,I_FeH,Z_iso,model,x,y,logt8,percorso,k,symmZ, &
			& calibHRD,calibNoD,calibSPEC,loggAvail,caliblogL,hstar,hstarlim,Isoc_i,Z_low,Z_up, &
			& useColor,idCol)
			IMPLICIT NONE
			CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: model,percorso
			INTEGER symmZ,idCol
			DOUBLE PRECISION FeH,I_FeH,x,y,k,hstar,hstarlim,Z_low,Z_up,logt8
			DOUBLE PRECISION useColor
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: Z_iso
			DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Isoc_i
			LOGICAL calibHRD,calibNoD,calibSPEC,loggAvail,caliblogL
		end subroutine multipleZisocPhSCP
		
		subroutine append(A,B)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: A
			DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: B
		end subroutine append
		
		subroutine append1D(A,B)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(inout) :: A
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: B
		end subroutine append1D
		
		subroutine ifComputeIsocAge(ZAMSz,cZAMS,logY0,logt_iso,logTeff_iso,yT_iso,I_logTeff,ageIsoc)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: ZAMSz
			DOUBLE PRECISION logY0,I_logTeff
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: logt_iso,logTeff_iso,yT_iso
			INTEGER cZAMS,ageIsoc
		end subroutine ifComputeIsocAge
		
		subroutine stepEnoughVar(tab,ix,varLim,colD,stepMax,backDir)
			IMPLICIT NONE
			DOUBLE PRECISION varLim
			DOUBLE PRECISION, DIMENSION(:,:) :: tab
			INTEGER ix,colD,stepMax
			LOGICAL backDir
		end subroutine stepEnoughVar
		
		subroutine computeVevo(t_Tracks,y_Tracks,logTeff_Tracks,t_i,varTrlim,vEvo)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: t_Tracks,y_Tracks,logTeff_Tracks,t_i
			DOUBLE PRECISION, DIMENSION(size(t_i)) :: vEvo
			DOUBLE PRECISION varTrlim
		end subroutine computeVevo
		subroutine computeVrif(TabVel,cRif,M,vRif)
			IMPLICIT NONE
			DOUBLE PRECISION M,vRif
			DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: TabVel
			INTEGER cRif
		end subroutine computeVrif
		
		subroutine removeColumn(matrix,titolo,field,newMatrix)
			IMPLICIT NONE
			DOUBLE PRECISION field
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: matrix
			DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: newMatrix
			CHARACTER(LEN=1000) :: titolo
		end subroutine removeColumn
		
		subroutine consistentM(TeB,cumInt,Mvec,M,agreeM)
			IMPLICIT NONE
			LOGICAL, DIMENSION(:), INTENT(in) :: TeB
			LOGICAL agreeM
			INTEGER, DIMENSION(:), INTENT(in) :: cumInt
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: Mvec
			DOUBLE PRECISION M
		end subroutine consistentM
		
	end interface
	
	interface
		function isEq_v(x1,x2,n)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: x1
			DOUBLE PRECISION, INTENT(in) :: x2
			INTEGER, INTENT(in) :: n
			LOGICAL, DIMENSION(size(x1)) :: isEq_v
		end function isEq_v
		function strcmp(str1,str2)
			IMPLICIT NONE
			CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: str1,str2
			LOGICAL strcmp
		end function strcmp
		function ones(n)
			IMPLICIT NONE
			INTEGER n
			DOUBLE PRECISION, DIMENSION(n) :: ones
		end function ones
		function polyval(v,x)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: v
			DOUBLE PRECISION x,polyval
		end function polyval
	end interface
	
!!!	TYPE item
!!!		CHARACTER(LEN=8) colLab
!!!		INTEGER colMag(8)
!!!	END TYPE
!!!	TYPE(item) Smag
	
	TYPE itName
		CHARACTER(LEN=20) fName
	END TYPE
	TYPE(itName), ALLOCATABLE :: StarName(:)
	
	model="PD1.2S"
	percorso="/home/bonfanti/Documents/PostDocLiegi/PARSEC1.2S/IsocroneRed/Z" !reduced tables
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	!!!!!!!!!!! Input/Output directories !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	!!!DIRECTORIES FOR ISOCHRONES AND STELLAR AGES 
	
	Est=trim(model)//'.txt' !!!approach
	
!!!	!!!!Defining plot labels in CMD 
!!!	Smag%colLab="UBVRIJHK"
!!!	Smag%colMag=(/(idum, idum=24,31)/) !column numbers where to find magnitudes from 24 to 31 in PD1.2S
!!!	cmag=(/ cmag1,cmag2 /) 
!!!	do i=1,2 
!!!		jS=1
!!!		do while (Smag%colMag(jS).ne.cmag(i) .and. jS<size(Smag%colMag)) 
!!!			jS=jS+1
!!!		end do 
!!!		if (Smag%colMag(jS).ne.cmag(i)) then 
!!!			write(*,*) 'Unable to find the selected magnitude in the isochrone grid'
!!!			stop
!!!		end if 
!!!		if (i==1) then 
!!!			mag1Label=Smag%colLab(jS:jS)
!!!			MAGLabel="M_{"//mag1Label//"}" !{} added if the subscript contains more than one character
!!!		else !i==2
!!!			mag2Label=Smag%colLab(jS:jS)
!!!		end if 
!!!	end do 
!!!	colLabel=mag2Label//"-"//mag1Label
!!!	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

	Z_isoLOW=minval(Zvec)
	Z_isoUP=maxval(Zvec)
	
	FeH_=SCP(2)
	if (IncZ==1) then 
		I_FeH_=SCP(3) !0.05;!
	end if
!	costZ=-log10(ZSun) !1.817; !!!Z=10.**(FeH-1.817)!!!
	Z_=10.**(FeH_-costZ)
	if (IncZ==1) then
		I_Z_=Z_*log(10.)*I_FeH_
	end if 
	
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	!!!!!!!!!!!!!!!! Name stars !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
	if (link.eq.1) then !just the 1st of the chain
	
		allocate(StarName(1)) !size(SCP,1)
		OPEN(UNIT=3, FILE=Isoch//NamesFile//'.txt', STATUS='OLD', IOSTAT=ios)
		if (ios.eq.0) then !file reporting star names exists
			do lines=1,1 !number of lines equals the size(SCP,1)
				READ(3,'(A)') StarName(lines)%fName
			end do 
		else
			do lines=1,1 !size(SCP,1)
				write(sdum,'(A3,I0)') "row",lines
				StarName(lines)%fName=sdum
			end do 
		end if
		CLOSE(3) 
		!Output stellar names
		outNames="Output"//trim(Est)
		OPEN(UNIT=3, FILE=Isoch//NamesFile//trim(outNames))
		OPEN(UNIT=4, FILE=Isoch//NamesFile//"DISC"//trim(outNames))
	
	end if
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	row=0
	rowN=0
	acc=1
	
	idCol=IDINT(SCP(29))
	select case (idCol)
	case (1)
		cmag0=cmag1
	case (2)
		cmag0=cmagx
	case default
		cmag0=cmag1
	end select
	
	do i=1,1
		indxZ=minloc(abs(Z_-Zvec),1)
		Zi=Zndxi(indxZ)
		Zf=Zndxf(indxZ)
		allocate(Isoc(Zf-Zi+1,size(IsocTab,2)))
		Isoc=IsocTab(Zi:Zf,:)
		Z=Zvec(indxZ)
		Zact=Z !Z of the very first isochrone grid to be opened. If element diffusion plays a role,
		!it could differ from the last and definitive metallic grid of isoc that will be used to
		!compute stellar parameters. In that case (Zact!=Zini) 
		Yact=dnint((He0+dYdZ*Zact)*100)/100
		Yini=Yact
		
		allocate(MTrAv(nM(indxZt(indxZ))))
		MTrAv=Mav(indxZt(indxZ),1:nM(indxZt(indxZ)))
		
		xVTi=Vndxi(indxZt(indxZ))
		xVTf=Vndxf(indxZt(indxZ))
		allocate(TabVel(xVTf-xVTi+1,size(velTrackTab,2)))
		TabVel=velTrackTab(xVTi:xVTf,:)
		
		xZAi=ZAndxi(indxZt(indxZ))
		xZAf=ZAndxf(indxZt(indxZ))
		
		allocate(ZAMSz(xZAf-xZAi+1,size(ZAMStab,2)))
		ZAMSz=ZAMStab(xZAi:xZAf,:)
		
		row_i=size(Isoc,1)
		allocate(t_iso(row_i)); allocate(logt_iso(row_i))
				
		t_iso=Isoc(:,ct)
		logt_iso=dnint(log10(t_iso)*100)/100
		t_iso=10.**logt_iso !to avoid numerical errors: in this way t_iso SHOULD come from logt
		! whose 2nd decimal digit is exactly 0 or 5, but Check!!!
		first_logt=logt_iso(1)
		last_logt=logt_iso(row_i)
		
		allocate(M_iso(row_i)); allocate(logL_iso(row_i))
		allocate(L_iso(row_i)); allocate(logTeff_iso(row_i))
		allocate(Teff_iso(row_i)); allocate(logg_iso(row_i))
		allocate(g_iso(row_i)); allocate(R_iso(row_i))
		allocate(rho_iso(row_i)); allocate(logrho_iso(row_i))
		
		M_iso=Isoc(:,cM)
		logL_iso=Isoc(:,clogL)
		L_iso=10.**(logL_iso)
		logTeff_iso=Isoc(:,clogTe)
		Teff_iso=10.**Isoc(:,clogTe)
		logg_iso=Isoc(:,clogg)
		g_iso=10.**Isoc(:,clogg)
		R_iso=sqrt(L_iso/(Teff_iso/TeffSun)**4) !RSun
		rho_iso=M_iso/(R_iso**3)*rhoSun !g/cm3
		logrho_iso=log10(rho_iso)
		if (photIsocAvail) then
			useColor=SCP(28)
			allocate(V_iso(row_i)); allocate(BmV_iso(row_i))
			allocate(BC_iso(row_i))
			V_iso=Isoc(:,cmag0)
			if (isEq(useColor,1.D0,2)) then
				BmV_iso=Isoc(:,cmag2)-Isoc(:,cmag1) !Isoc(:,9)-Isoc(:,11);! !10,11:B-V ->PDisoc 9,11:J-K ->2MASSisoc
			else
				BmV_iso=logTeff_iso
				InclogTeff=SCP(5)/SCP(4)*log10(e) !uncertainty in the input logTeff
				DTeJ=InclogTeff !in checking logg-logrho compatibility use the real known Teff uncertainty
			end if
			!!! 
			BC_iso=Isoc(:,cmbol)-Isoc(:,cmag0)
			!!! 
		end if 

		!calibHRD: input V, B-V, d. (Mv, B-V)-->(L, Teff). => R. gProxy preferable, though not compulsory
		!calibNoD: input V, B-V. gProxy compulsory. (gProxy, B-V)-->(gProxy, Teff) 
		!calibSPEC: input Teff, gProxy 
		calibHRD=.false.
		calibNoD=.false.
		calibSPEC=.false.
		gProxyAvail=.true. !.true. as default. It could become .false. later on if neither logg nor rho are avail
		!this could happen if  (calibHRD), where I will perform just a geometrical calibration in the (Teff,L)
		!plane. The estimate of age is not so reliable.
		rhoAvail=.false.
		loggAvail=.false.
		rhoAvailAS=.false.
		loggAvailAS=.false.
		loggAvailCal=.false. !it becomes .true. if it's inferred from calibration
							 !(e.g. after calibHRD, if rho is known, M and logg are inferred)
		rhoAvailCal=.false.
		bestLogg=.false.
		
		!If asteroseismology is available, I recover logg from Dnu, NuMax, Teffinput 
		!Otherwise I'll use spectroscopic logg
		!!!!!!!!!!!!!!!!!!!!!!!!!!! 
		!!!!CALIBRATION CHOICE!!!!	29/4/15-->Added the possibility to enter just spectroscopic parameters
		!!!!!!!!!!!!!!!!!!!!!!!!!!!
		sismo=SCP(20) !single(SCP(i,)); 0=NOT available; 1=available
		Vf=SCP(23)
		if (isEq(useColor,1.D0,2)) then
			BmVf=SCP(24)
		else
			BmVf=log10(SCP(4)) !actually logTeff
		end if
		!!!Inizialized to -1. Assume the respective values rho or logg if they are actually available
		rhoI=-1.
		loggI=-1.
		if ((.not.isEq(Vf,-100.D0,2)).and.(.not.isEq(BmVf,-100.D0,2))) then !Photometry is available
			V=Vf
			BmV=BmVf
			Teffinput=-1 !SCP(88);! It should be available if I use asteroseismo
			!To be compared to the first calibrated Teff from BmV in case of doubtful photometry
			Hipf=SCP(27)  !Hipparcos Flag. -1=no plx distance available
			if (.not.isEq(Hipf,-1.D0,2)) then !plx distance available
				calibHRD=.true.
				
				Rf=SCP(21) !differ from -1 if available from input (measured through interferometry)
							
				d=SCP(25)
				I_d=SCP(26)
				DM0=5*log10(d)-5
				Vass=V-DM0
				VMag=Vass
				rhof=SCP(6) !-1;!		 -1=rho NOT available
				if (.not.isEq(rhof,-1.D0,2)) then
					rhoS=rhof
					rhoI=rhoS
					I_rhoS=SCP(7) !-1;!
					logrhoS=log10(rhoS)
					I_logrhoS=log10(e)*I_rhoS/rhoS
					rhoAvail=.true.
				end if 
				loggf=SCP(8) !-1;!
				if (.not.isEq(loggf,-1.D0,2)) then
					loggS=loggf
					loggI=loggS
					I_loggS=SCP(9) !0.15;!0.06;!
					gS=10.**loggS
					I_gS=gS*I_loggS*log(10.)
					loggAvail=.true. !logg surely available. rho is optional
				end if
				if (sismo.ne.0) then 
					!!!!ASTEROSEISMOLOGY available!!! 
					numax=SCP(16)
					I_numax=SCP(17)
					Dnu=SCP(18)
					I_Dnu=SCP(19)
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					if (.not.isEq(Dnu,-1.D0,2)) then !Dnu available
						rhoAS=(Dnu/DnuSun)**2*rhoSun
						rhoI=rhoAS
						I_rhoAS=rhoAS*2.*I_Dnu/Dnu
						logrhoAS=log10(rhoAS)
						I_logrhoAS=log10(e)*I_rhoAS/rhoAS
						rhoAvailAS=.true.
					end if
					if (.not.isEq(numax,-1D0,2)) then !numax available
						if (isEq(Teffinput,-1.D0,2)) then !Teff input not available
							if (isEq(useColor,1.D0,2)) then
								logTeJ=polyval(cJ(:,idCol),BmV)
								!3.908-0.234*BmV !logTe-BmV relation according to Johnson (1966)
								Teffinput0=10.**logTeJ
							else
								Teffinput0=10.**BmV !actually the temperature
							end if
						else
							Teffinput0=Teffinput
						end if
						gAS=numax/numaxSun*sqrt(Teffinput0/TeffSun)*10.**loggSun
						I_gAS=gAS*(I_numax/numax+.5*.02) !Inc rel Teff ~2%
						loggAS=log10(gAS)
						loggI=loggAS
						I_loggAS=I_gAS/gAS*log10(e)
						loggAvailAS=.true.
					end if
				end if
				if (loggAvail) then
					if (loggAvailAS) then !both available: weighted mean between two logg determinations
						call wMean(gS,I_gS,gAS,I_gAS,g,I_g)
!						pS=I_gS**(-2)
!						pAS=I_gAS**(-2)
!						g=(gS*pS+gAS*pAS)/(pS+pAS)
!						I_g=1./sqrt(pS+pAS)
						logg=log10(g)
						I_logg=I_g/g*log10(e)
					else !only logg from spectroscopy
						g=gS
						I_g=I_gS
						logg=loggS
						I_logg=I_loggS
					end if
				else
					if (loggAvailAS) then !only logg from asteroseismology
						g=gAS
						I_g=I_gAS
						logg=loggAS
						I_logg=I_loggAS
					end if
				end if
				if (rhoAvail) then
					if (rhoAvailAS) then !both available => weighted mean
						call wMean(rhoS,I_rhoS,rhoAS,I_rhoAS,rho,I_rho)
!						pS=I_rhoS**(-2)
!						pAS=I_rhoAS**(-2)
!						rho=(pS*rhoS+pAS*rhoAS)/(pS+pAS)
!						I_rho=1./sqrt(pS+pAS)
						logrho=log10(rho)
						I_logrho=I_rho/rho*log10(e)
					else !only rho from transit
						rho=rhoS
						I_rho=I_rhoS
						logrho=logrhoS
						I_logrho=I_logrhoS
					end if
				else !only rho from asteroseismology
					if (rhoAvailAS) then
						rho=rhoAS
						I_rho=I_rhoAS
						logrho=logrhoAS
						I_logrho=I_logrhoAS
					end if
				end if
				if (.not.rhoAvail.and..not.loggAvail.and..not.rhoAvailAS.and..not.loggAvailAS) then !neither logg, nor rho
					gProxyAvail=.false.
				end if 
				!!!	(rho,logg) available; I could determine input values for M, R. But it's better to 
				!!!	find them after the calibration M_V-->L, BmV-->Teff 
				!!	R=g/10^loggSun*rhoSun/rho; 
				!!	I_R=R*(I_g/g+I_rho/rho); 
				!!	M=(g/10^loggSun)^3*(rhoSun/rho)^2; 
				!!	I_M=M*(3*I_g/g+2*I_rho/rho); 

				hstar=VMag !Threshold to go forward/backward with the indices during the choice of ref isoc
				hstarlim=7 !In this case the vertical coord of the CMD is the absolute mag
			else !I have photometry + gProxy, but NOT the plx distance
				calibNoD=.true. !so the calibration is just related to BmV-->Teff from the gProxy(or logL if Rinp)-BmV plane
				
				rhof=SCP(6) !-1;!		 -1=rho NOT available
				if (.not.isEq(rhof,-1.D0,2)) then
					rhoS=rhof
					rhoI=rhoS
					I_rhoS=SCP(7) !-1;!
					logrhoS=log10(rhoS)
					I_logrhoS=log10(e)*I_rhoS/rhoS
					rhoAvail=.true.
				end if 
				loggf=SCP(8) !-1;!
				if (.not.isEq(loggf,-1.D0,2)) then
					loggS=loggf
					loggI=loggS
					I_loggS=SCP(9) !0.15;!0.06;!
					gS=10.**loggS
					I_gS=gS*I_loggS*log(10.)
					loggAvail=.true. !logg surely available. rho is optional
				end if
				if (sismo.ne.0) then 
					!!!!ASTEROSEISMOLOGY available!!! 
					numax=SCP(16)
					I_numax=SCP(17)
					Dnu=SCP(18)
					I_Dnu=SCP(19)
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
					if (.not.isEq(Dnu,-1.D0,2)) then !Dnu available
						rhoAS=(Dnu/DnuSun)**2*rhoSun
						rhoI=rhoAS
						I_rhoAS=rhoAS*2.*I_Dnu/Dnu
						logrhoAS=log10(rhoAS)
						I_logrhoAS=log10(e)*I_rhoAS/rhoAS
						rhoAvailAS=.true.
					end if
					if (.not.isEq(numax,-1D0,2)) then !numax available
						if (isEq(Teffinput,-1.D0,2)) then !Teff input not available
							logTeJ=polyval(cJ(:,idCol),BmV)
							!3.908-0.234*BmV !logTe-BmV relation according to Johnson (1966)
							Teffinput0=10.**logTeJ
						else
							Teffinput0=Teffinput
						end if
						gAS=numax/numaxSun*sqrt(Teffinput0/TeffSun)*10.**loggSun
						I_gAS=gAS*(I_numax/numax+.5*.02) !Inc rel Teff ~2%
						loggAS=log10(gAS)
						loggI=loggAS
						I_loggAS=I_gAS/gAS*log10(e)
						loggAvailAS=.true.
					end if 
				end if
				if (loggAvail) then
					if (loggAvailAS) then !both available: weighted mean between two logg determinations
						call wMean(gS,I_gS,gAS,I_gAS,g,I_g)
						logg=log10(g)
						I_logg=I_g/g*log10(e)
					else !only logg from spectroscopy
						g=gS
						I_g=I_gS
						logg=loggS
						I_logg=I_loggS
					end if
				else
					if (loggAvailAS) then !only logg from asteroseismology
						g=gAS
						I_g=I_gAS
						logg=loggAS
						I_logg=I_loggAS
					end if
				end if
				if (rhoAvail) then
					if (rhoAvailAS) then !both available => weighted mean
						call wMean(rhoS,I_rhoS,rhoAS,I_rhoAS,rho,I_rho)
						logrho=log10(rho)
						I_logrho=I_rho/rho*log10(e)
					else !only rho from transit
						rho=rhoS
						I_rho=I_rhoS
						logrho=logrhoS
						I_logrho=I_logrhoS
					end if
				else !only rho from asteroseismology
					if (rhoAvailAS) then
						rho=rhoAS
						I_rho=I_rhoAS
						logrho=logrhoAS
						I_logrho=I_logrhoAS
					end if
				end if
				
				Rf=SCP(21) !differ from -1 if available from input (measured through interferometry)
				if (.not.(isEq(Rf,-1.D0,2))) then
					R1=Rf
					I_R1=SCP(22)
					logTeJ=polyval(cJ(:,idCol),BmV)
					!3.908-0.234*BmV !logTe-BmV relation according to Johnson (1966)
					logL=2.*log10(R1)+4.*(logTeJ-log10(TeffSun)) !rough estimate so far	
				end if
				if (.not.loggAvail.and..not.rhoAvail.and..not.loggAvailAS.and..not.rhoAvailAS) then !neither logg, nor rho
					gProxyAvail=.false.
					if (isEq(Rf,-1.D0,2)) then !neither R to then compute L once Teff has been calibrated
						print*,'No Luminosity, nor gProxy are available'; stop
					end if
				end if
				
				if (.not.isEq(Rf,-1.D0,2)) then !necessarily R available => logL avail
					hstar=-logL !minus to let hstar>hstarlim refer to low MS
					hstarlim=-(-0.7D0)
				else
					if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. .not.bestLogg) &
						& .or. ((rhoAvail.or.rhoAvailAS).and.(.not.(loggAvail.or.loggAvailAS)))) then
						hstar=logrho
						if (Z<0.010) then 
							hstarlim=-0.25D0
						else
							hstarlim=-1.D0
						end if 
					else
						hstar=logg
						if (Z<0.010) then 
							hstarlim=4.4D0
						else
							hstarlim=3.9D0
						end if 
					end if
				end if
			end if 
		else !No photometry => just spectroscopy
			calibSPEC=.true.
			Teff=SCP(4)
			I_Teff=SCP(5) !0.01*Teff;!
			logTeff=log10(Teff)
			I_logTeff=I_Teff/Teff*log10(e)

			rhof=SCP(6) !-1;!		 -1=rho NOT available
			if (.not.isEq(rhof,-1.D0,2)) then
				rhoS=rhof
				rhoI=rhoS
				I_rhoS=SCP(7) !-1;!
				logrhoS=log10(rhoS)
				I_logrhoS=log10(e)*I_rhoS/rhoS
				rhoAvail=.true.
			end if 
			loggf=SCP(8) !-1;!
			if (.not.isEq(loggf,-1.D0,2)) then
				loggS=loggf
				loggI=loggS
				I_loggS=SCP(9) !0.15;!0.06;!
				gS=10.**loggS
				I_gS=gS*I_loggS*log(10.)
				loggAvail=.true. !logg surely available. rho is optional
			end if
			if (sismo.ne.0) then 
				!!!!ASTEROSEISMOLOGY available!!! 
				numax=SCP(16)
				I_numax=SCP(17)
				Dnu=SCP(18)
				I_Dnu=SCP(19)
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
				if (.not.isEq(Dnu,-1.D0,2)) then !Dnu available
					rhoAS=(Dnu/DnuSun)**2*rhoSun
					rhoI=rhoAS
					I_rhoAS=rhoAS*2.*I_Dnu/Dnu
					logrhoAS=log10(rhoAS)
					I_logrhoAS=log10(e)*I_rhoAS/rhoAS
					rhoAvailAS=.true.
				end if
				if (.not.isEq(numax,-1D0,2)) then !numax available
					gAS=numax/numaxSun*sqrt(Teff/TeffSun)*10.**loggSun
					I_gAS=gAS*(I_numax/numax+.5*I_Teff/Teff) !Inc rel Teff ~2%
					loggAS=log10(gAS)
					loggI=loggAS
					I_loggAS=I_gAS/gAS*log10(e)
					loggAvailAS=.true.
				end if 
			end if
			if (loggAvail) then
				if (loggAvailAS) then !both available: weighted mean between two logg determinations
					call wMean(gS,I_gS,gAS,I_gAS,g,I_g)
					logg=log10(g)
					I_logg=I_g/g*log10(e)
				else !only logg from spectroscopy
					g=gS
					I_g=I_gS
					logg=loggS
					I_logg=I_loggS
				end if
			else
				if (loggAvailAS) then !only logg from asteroseismology
					g=gAS
					I_g=I_gAS
					logg=loggAS
					I_logg=I_loggAS
				end if
			end if
			if (rhoAvail) then
				if (rhoAvailAS) then !both available => weighted mean
					call wMean(rhoS,I_rhoS,rhoAS,I_rhoAS,rho,I_rho)
					logrho=log10(rho)
					I_logrho=I_rho/rho*log10(e)
				else !only rho from transit
					rho=rhoS
					I_rho=I_rhoS
					logrho=logrhoS
					I_logrho=I_logrhoS
				end if
			else !only rho from asteroseismology
				if (rhoAvailAS) then
					rho=rhoAS
					I_rho=I_rhoAS
					logrho=logrhoAS
					I_logrho=I_logrhoAS
				end if
			end if
			Rf=SCP(21)
			if (.not.(isEq(Rf,-1.D0,2))) then
				R1=Rf
				I_R1=SCP(22)	
			end if
			if (.not.rhoAvail.and..not.loggAvail.and..not.rhoAvailAS.and..not.loggAvailAS) then !logg, nor rho available
				gProxyAvail=.false.
				if (isEq(Rf,-1.D0,2)) then
					print*,'No Luminosity, nor gProxy are available'; stop
				end if
			end if
			if ((loggAvail.or.loggAvailAS).and.(rhoAvail.or.rhoAvailAS)) then 
				!!rho, logg both available. I determine input values for M, R 
				R2=g/10.**loggSun*rhoSun/rho
				I_R2=R2*(I_g/g+I_rho/rho)
				if (isEq(Rf,-1.D0,2)) then
					R=R2
					I_R=I_R2
				else
					call wMean(R1,I_R1,R2,I_R2,R,I_R)
				end if
				M=(g/10.**loggSun)**3*(rhoSun/rho)**2
				I_M=M*(3.*I_g/g+2.*I_rho/rho)
				L=R**2*(Teff/TeffSun)**4
				I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
				logL=log10(L)
				I_logL=I_L/L*log10(e)
				!!!!!!Check compatibility between rho and logg 16/10/17 
				!!Modified to account for the mass uncertainty 26/10/2018 
				if (.not.(M+I_M.lt.minval(MTrAv) .or. M-I_M.gt.maxval(MTrAv))) then !open tracks only if inside MassRange 
					call selectMfromTracks(MTrAv,M,I_M,Mvec)
					allocate(cumInt(size(Mvec)))
					jk=0
					do jj=1,size(Mvec)
						indxM=minloc(abs(Mvec(jj)-MTrAv),1)
						xMTi=Tndxi(indxZt(indxZ),indxM)
						xMTf=Tndxf(indxZt(indxZ),indxM)
						allocate(TrRhoG(xMTf-xMTi+1,size(TrackTab,2)))
						TrRhoG=TrackTab(xMTi:xMTf,:)
	
						allocate(MTr(size(TrRhoG,1))); allocate(RTr(size(TrRhoG,1)))
						allocate(logRhoTr(size(TrRhoG,1))); allocate(loggTr(size(TrRhoG,1)))
						
						MTr=TrRhoG(:,cM_T) !Msun
						RTr=10.**TrRhoG(:,clogR_T) !cm
						logRhoTr=log10(MTr/(RTr/RSun)**3*rhoSun)
						loggTr=log10(MTr/(RTr/RSun)**2)+loggSun
						call findYgivenX_v(TrRhoG,logg,loggTr,clogTe_T,logTeTrtmp)
						if (allocated(logTeTrtmp)) then
							jk=jk+1
							if (jk.eq.1) then
								allocate(logTeTr(size(logTeTrtmp)))
								logTeTr=logTeTrtmp
							else
								call append1D(logTeTr,logTeTrtmp)
							end if
							cumInt(jj)=size(logTeTr)
							deallocate(logTeTrtmp)
						else
							cumInt(jj)=0
						end if
						
						deallocate(TrRhoG)
						deallocate(MTr); deallocate(RTr)
						deallocate(logRhoTr); deallocate(loggTr)
					end do
					logTeJ=logTeff !I don't use Johnson (1966) relation because I don't have colour index
					DTeJ=I_logTeff !directly set to the true uncertainty
					if (allocated(logTeTr)) then 
						allocate(TeB(size(logTeTr)))
						TeB=((logTeTr>logTeJ-DTeJ).and.(logTeTr<logTeJ+DTeJ))
						if (.not.any(TeB)) then !all elements of TeB are 0 => logg and rho inconsistent
							if (I_logg>I_logrho) then !logrho is best determined
								if (isEq(Rf,-1.D0,2)) then
									call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
									call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
									loggf=-1.
									call setNaN(logg); call setNaN(I_logg)
									call setNaN(g); call setNaN(I_g)
									loggAvail=.false.
									loggAvailAS=.false.
								else !input R is available
									!Only using rho: determine M; re-determine logg (discard the input value) 
									loggf=-1.
									loggAvail=.false.
									loggAvailAS=.false.
									R=R1 !use just the input reliable value
									I_R=I_R1
									call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
									loggAvailCal=.true. !!logg now available thanks to calibration
								end if
							else !logg is best determined
								if (isEq(Rf,-1.D0,2)) then
									call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
									call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
									rhof=-1.
									call setNaN(rho); call setNaN(I_rho)
									call setNaN(logrho); call setNaN(I_logrho)
									rhoAvail=.false.
									rhoAvailAS=.false.
								else !input R is available
									!Only using logg: determine M; re-determine rho (discard the input value) 
									rhof=-1.
									rhoAvail=.false.
									rhoAvailAS=.false.
									R=R1 !use just the input reliable value
									I_R=I_R1
									call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
									rhoAvailCal=.true.
								end if
							end if 
						else
							!!!!!!!!!
							call consistentM(TeB,cumInt,Mvec,M,agreeM)
							if (.not.agreeM) then !new M has been recomputed inside the subroutine in case agreeM is false
											      !This M will be used in case R is not available directly from input
								if (I_logg.gt.I_logrho) then !logrho is best determined
									if (isEq(Rf,-1.D0,2)) then
										loggAvail=.false.
										loggAvailAS=.false.
										!g: mantain original input uncertainty
										call computeLRgfromMrho(M,I_M,rho,I_rho,Teff,I_Teff,L,logL,I_L,I_logL,R,I_R,g,logg)
										loggAvailCal=.true.
									else !input R is available
										!Only using rho: determine M; re-determine logg (discard the input value) 
										loggf=-1.
										loggAvail=.false.
										loggAvailAS=.false.
										R=R1
										I_R=I_R1
										call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
										loggAvailCal=.true. !!logg now available thanks to calibration
									end if
								else
									if (isEq(Rf,-1.D0,2)) then
										rhoAvail=.false.
										rhoAvailAS=.false.
										!rho in g/cm3. Mantain original input uncertainty on rho
										call computeLRrhofromMg(M,I_M,g,I_g,Teff,I_Teff,L,logL,I_L,I_logL,R,I_R,rho,logrho)
										rhoAvailCal=.true.
									else !input R is available
										!Only using logg: determine M; re-determine rho (discard the input value) 
										rhof=-1.
										rhoAvail=.false.
										rhoAvailAS=.false.
										R=R1
										I_R=I_R1
										call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
										rhoAvailCal=.true.
									end if
								end if
								if (isEq(Rf,-1.D0,2)) then
									R=g/10**loggSun*rhoSun/rho
									I_R=R*(I_g/g+I_rho/rho)
									L=R**2*(Teff/TeffSun)**4
									I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
									logL=log10(L)
									I_logL=I_L/L*log10(e)
								end if
							else
								if (I_logg.lt.I_logrho) then
									bestLogg=.true.
								else
									bestLogg=.false.
								end if
								if (.not.isEq(Rf,-1.D0,2)) then
								!input R is available: weighted mean to infer M
									Mlogg=g/10.**loggSun*R**2
									I_Mlogg=Mlogg*(I_g/g+2.*I_R/R)
									Mrho=rho/rhoSun*R**3 !Mo
									I_Mrho=Mrho*(I_rho/rho+3.*I_R/R)
									call wMean(Mlogg,I_Mlogg,Mrho,I_Mrho,M,I_M)
									!L already defined since the beginning from R
								end if
							end if
							!!!!!!!!!!!!!
						end if
						deallocate(logTeTr); deallocate(TeB)
					else !Track has been opened, but no logTe compatible with logg were found
						if (I_logg>I_logrho) then !logrho is best determined
							if (isEq(Rf,-1.D0,2)) then
								call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
								call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
								loggf=-1.
								call setNaN(logg); call setNaN(I_logg)
								call setNaN(g); call setNaN(I_g)
								loggAvail=.false.
								loggAvailAS=.false.
							else !input R is available
								!Only using rho: determine M; re-determine logg (discard the input value) 
								loggf=-1.
								loggAvail=.false.
								loggAvailAS=.false.
								R=R1 !use just the input reliable value
								I_R=I_R1
								call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
								loggAvailCal=.true. !!logg now available thanks to calibration
							end if
						else !logg is best determined
							if (isEq(Rf,-1.D0,2)) then
								call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
								call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
								rhof=-1.
								call setNaN(rho); call setNaN(I_rho)
								call setNaN(logrho); call setNaN(I_logrho)
								rhoAvail=.false.
								rhoAvailAS=.false.
							else !input R is available
								!Only using logg: determine M; re-determine rho (discard the input value) 
								rhof=-1.
								rhoAvail=.false.
								rhoAvailAS=.false.
								R=R1 !use just the input reliable value
								I_R=I_R1
								call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
								rhoAvailCal=.true.
							end if
						end if 
					end if
					deallocate(cumInt) 
				else
					if (I_logg>I_logrho) then !logrho is best determined
						if (isEq(Rf,-1.D0,2)) then
							call setNan(R); call setNan(I_R); call setNan(M); call setNan(I_M)
							call setNan(L); call setNan(I_L); call setNan(logL); call setNan(I_logL)
							loggf=-1.
							call setNan(logg); call setNan(I_logg)
							call setNan(g); call setNan(I_g)
							loggAvail=.false.
							loggAvailAS=.false.
						else !input R is available
							!Only using rho: determine M; re-determine logg (discard the input value) 
							loggf=-1.
							loggAvail=.false.
							loggAvailAS=.false.
							R=R1 !use just the input reliable value
							I_R=I_R1
							call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
							loggAvailCal=.true. !!logg now available thanks to calibration
						end if
					else !logg is best determined
						if (isEq(Rf,-1.D0,2)) then
							call setNan(R); call setNan(I_R); call setNan(M); call setNan(I_M)
							call setNan(L); call setNan(I_L); call setNan(logL); call setNan(I_logL)
							rhof=-1.
							call setNan(rho); call setNan(I_rho)
							call setNan(logrho); call setNan(I_logrho)
							rhoAvail=.false.
							rhoAvailAS=.false.
						else !input R is available
							!Only using logg: determine M; re-determine rho (discard the input value) 
							rhof=-1.
							rhoAvail=.false.
							rhoAvailAS=.false.
							R=R1 !use just the input reliable value
							I_R=I_R1
							call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
							rhoAvailCal=.true.
						end if
					end if 
				end if 
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			else
				if (.not.isEq(Rf,-1.D0,2)) then
					R=R1
					I_R=I_R1
					if (loggAvail.or.loggAvailAS) then
						call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
						rhoAvailCal=.true.
					else if (rhoAvail.or.rhoAvailAS) then
						call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
						loggAvailCal=.true. !!logg available thanks to calibration
					else !noGproxy
						L=R**2*(Teff/TeffSun)**4
						I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
						logL=log10(L)
						I_logL=I_L/L*log10(e)
					end if
				end if
			end if 

			if (.not.isEq(Rf,-1.D0,2)) then !necessarily R available => logL avail
				hstar=-logL !minus to let hstar>hstarlim refer to low MS
				hstarlim=-(-0.7D0)
			else
				if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. .not.bestLogg) &
					& .or. ((rhoAvail.or.rhoAvailAS).and.(.not.(loggAvail.or.loggAvailAS)))) then 
					hstar=logrho
					if (Z<0.010) then 
						hstarlim=-0.25D0
					else
						hstarlim=-1.D0
					end if 
				else
					hstar=logg
					if (Z<0.010) then 
						hstarlim=4.4D0
					else
						hstarlim=3.9D0
					end if 
				end if
			end if
		end if 
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
		!!!!!! END CALIBRATION CHOICE !!!!!!!!! 
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		logRHK=SCP(14) !0;!
		vsini=SCP(10) !-1;!
		Uvsini=SCP(11) !-1;!
		P=SCP(12) !-1. ! !!!Stellar rotational period [days]
		UP=SCP(13) !-1. !
		YMg=SCP(15)
		
		Bin=-1 !SCP(99);! !!!flag changed to -1 if not binary. Before it was 0 if not binary
		!!!Selection of reference isochrones!!! 
		!!First evaluate which isochrones are the closest to the star to be analysed
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
		!!!!!!!NEW CODE 12/2/14!!!!!!!!!!!!!!! 
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
		allocate(dist_ndx(size(logTeff_iso))) !equal to any other size of _iso
		allocate(x_iso(size(logTeff_iso))); allocate(xx_iso(size(logTeff_iso)))
		if (calibHRD) then
			call indexx(sqrt((BmV-BmV_iso)**2+(Vass-V_iso)**2), dist_ndx) 
			x=BmV
			x_iso=BmV_iso
			xx_iso=logTeff_iso
			if (isEq(useColor,1.D0,2)) then
				if (idCol.eq.1) then
					cBrif=cBBmV
				else
					cBrif=cBTeff
				end if
			else
				cBrif=cBTeff
			end if
		end if 
		if (calibNoD) then
			if (.not.isEq(Rf,-1.D0,2)) then
				call indexx(sqrt((BmV-BmV_iso)**2+(logL-logL_iso)**2),dist_ndx) !logL so far just estimated
			else if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
				& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))) then 
				call indexx(sqrt((BmV-BmV_iso)**2+(logg-logg_iso)**2),dist_ndx) 
			else !rhoAvail
				call indexx(sqrt((BmV-BmV_iso)**2+(logrho-logrho_iso)**2),dist_ndx) 
			end if 
			x=BmV
			x_iso=BmV_iso
			xx_iso=logTeff_iso
			cBrif=cBBmV
		end if 
		if (calibSPEC) then
			if (.not.isEq(Rf,-1.D0,2)) then
			 	call indexx(sqrt((logTeff-logTeff_iso)**2+(logL-logL_iso)**2),dist_ndx)
			else if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
				& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))) then 
				call indexx(sqrt((logTeff-logTeff_iso)**2+(logg-logg_iso)**2), dist_ndx) 
			else !rhoAvail
				call indexx(sqrt((logTeff-logTeff_iso)**2+(logrho-logrho_iso)**2), dist_ndx) 
			end if 
			x=logTeff
			xtau=x
			x_iso=logTeff_iso
			xx_iso=BmV_iso
			cBrif=cBTeff
		end if
	
		allocate(ageMinDist(size(logt_iso)))
		ageMinDist=logt_iso(dist_ndx)
		
		call uniqueFast(ageMinDist,2,ix,.true.)
		call sort0(ix) !so to display elements of ageMinDist in the oroginal order
		allocate(dist_ndxU(size(ix))) 
		dist_ndxU=dist_ndx(ix)
		
		deallocate(dist_ndx); deallocate(ageMinDist)
		deallocate(ix)
		rowAge=1
		! allocation of already defined variables
		allocate(t_iTmp(size(dist_ndxU)));allocate(distTmp(size(dist_ndxU),2))
		allocate(logTeff_iTmp(size(dist_ndxU)));allocate(logL_iTmp(size(dist_ndxU)))
		allocate(M_iTmp(size(dist_ndxU)));allocate(logg_iTmp(size(dist_ndxU)))
		allocate(logrho_iTmp(size(dist_ndxU)))
		if (photIsocAvail) then 
			allocate(BmV_iTmp(size(dist_ndxU)));allocate(BC_iTmp(size(dist_ndxU)))
		end if
		!
		do Nage=1,size(dist_ndxU) 
			ndxDU=dist_ndxU(Nage)
			if (photIsocAvail) then 
				BmV_i1=BmV_iso(ndxDU)
				V_i1=V_iso(ndxDU)
			end if 
			logTeff_i1=logTeff_iso(ndxDU)
			logL_i1=logL_iso(ndxDU)
			M_i1=M_iso(ndxDU)
			logg_i1=logg_iso(ndxDU)
			logrho_i1=logrho_iso(ndxDU)
			
			call choose_i2Col(ndxDU,x,x_iso,xx_iso,calibSPEC,hstar,hstarlim,V_iso,logL_iso, &
				& M_iso,logg_iso,logrho_iso,logt_iso,BmV_i2,V_i2,logTeff_i2,logL_i2,M_i2,logg_i2,logrho_i2)
			if (calibSPEC) then
				x_i1=logTeff_i1
				x_i2=logTeff_i2
			else
				x_i1=BmV_i1
				x_i2=BmV_i2
			end if
			if (isEq(x_i1,x_i2,4)) then 
				cycle
			end if 
			t_iTmp(rowAge)=t_iso(ndxDU)
			if (photIsocAvail) then 
				BC_iTmp(rowAge)=Isoc(ndxDU,cmbol)-Isoc(ndxDU,cmag0) !just rough.
				!When I recover BC in the logg-BmV plane, I'll make the interpolation using BC_i1 e BC_i2
			end if 
			if (calibHRD) then !!Condiz sulla calibrazione da fare 15/12/14
				call findX_i(BmV,Vass,BmV_i1,BmV_i2,V_i1,V_i2,x_i,dist0)
				BmV_iTmp(rowAge)=x_i
			end if 
			if (calibNoD) then
				if (.not.(isEq(Rf,-1.D0,2))) then
					call findX_i(BmV,logL,BmV_i1,BmV_i2,logL_i1,logL_i2,x_i,dist0) 
				else if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
				& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))) then 
					call findX_i(BmV,logg,BmV_i1,BmV_i2,logg_i1,logg_i2,x_i,dist0)
				else !rhoAvail
					call findX_i(BmV,logrho,BmV_i1,BmV_i2,logrho_i1,logrho_i2,x_i,dist0)
				end if 
				BmV_iTmp(rowAge)=x_i
			end if 
			if (calibSPEC) then
				if (.not.(isEq(Rf,-1.D0,2))) then
					call findX_i(logTeff,logL,logTeff_i1,logTeff_i2,logL_i1,logL_i2,x_i,dist0) 
				else if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
				& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))) then 
					call findX_i(logTeff,logg,logTeff_i1,logTeff_i2,logg_i1,logg_i2,x_i,dist0)
				else !rhoAvail
					call findX_i(logTeff,logrho,logTeff_i1,logTeff_i2,logrho_i1,logrho_i2,x_i,dist0)
				end if 
				logTeff_iTmp(rowAge)=x_i 
			end if
			distTmp(rowAge,:)=(/ dist0, logt_iso(ndxDU) /)
						
			!!18/3/15 The following if cancel the criterion of perpendicular distance
			!!between a star and the isochrone if the theoretical point to be interpolated
			!!doesn't fall inside the isochrone (i.e. between x_i1 and x_i2), but on its
			!!extension. This possibility, in fact, would build fictitious isochrones that
			!!could be erroneously close to the star
			if (.not.((x_i>=x_i1.and.x_i<=x_i2).or.(x_i>=x_i2.and.x_i<=x_i1))) then 
				if (photIsocAvail) then 
					BmV_iTmp(rowAge)=BmV_i1
				end if 
				logg_iTmp(rowAge)=logg_i1
				logTeff_iTmp(rowAge)=logTeff_i1
				logL_iTmp(rowAge)=logL_i1
				M_iTmp(rowAge)=M_i1
				logrho_iTmp(rowAge)=logrho_i1
			else
				if (calibSPEC) then 
					if (photIsocAvail) then 
						mBV=(BmV_i2-BmV_i1)/(x_i2-x_i1)
						qBV=-mBV*x_i1+BmV_i1
						BmV_iTmp(rowAge)=mBV*x_i+qBV
					end if 
				else
					mT=(logTeff_i2-logTeff_i1)/(x_i2-x_i1)
					qT=-mT*x_i1+logTeff_i1
					logTeff_iTmp(rowAge)=mT*x_i+qT
				end if
				mg=(logg_i2-logg_i1)/(x_i2-x_i1)
				qg=-mg*x_i1+logg_i1
				logg_iTmp(rowAge)=mg*x_i+qg
				mL=(logL_i2-logL_i1)/(x_i2-x_i1)
				qL=-mL*x_i1+logL_i1
				logL_iTmp(rowAge)=mL*x_i+qL
				mM=(M_i2-M_i1)/(x_i2-x_i1)
				qM=-mM*x_i1+M_i1
				M_iTmp(rowAge)=mM*x_i+qM
				mrh=(logrho_i2-logrho_i1)/(x_i2-x_i1)
				qrho=-mrh*x_i1+logrho_i1
				logrho_iTmp(rowAge)=mrh*x_i+qrho
			end if 

			rowAge=rowAge+1
		end do
		
		deallocate(x_iso); deallocate(xx_iso)
				
		allocate(t_i(rowAge-1));allocate(dist(rowAge-1,2))
		allocate(logTeff_i(rowAge-1));allocate(logL_i(rowAge-1))
		allocate(M_i(rowAge-1));allocate(logg_i(rowAge-1))
		allocate(logrho_i(rowAge-1))
		allocate(BmV_i(rowAge-1));allocate(BC_i(rowAge-1))
		t_i=t_iTmp(1:rowAge-1)
		dist=distTmp(1:rowAge-1,:)
		logTeff_i=logTeff_iTmp(1:rowAge-1)
		logL_i=logL_iTmp(1:rowAge-1)
		M_i=M_iTmp(1:rowAge-1)
		logg_i=logg_iTmp(1:rowAge-1)
		logrho_i=logrho_iTmp(1:rowAge-1)
		BmV_i=BmV_iTmp(1:rowAge-1)
		BC_i=BC_iTmp(1:rowAge-1)
		deallocate(t_iTmp);deallocate(distTmp);deallocate(logTeff_iTmp)
		deallocate(logL_iTmp);deallocate(M_iTmp);deallocate(logg_iTmp)
		deallocate(logrho_iTmp)
		if (photIsocAvail) then
			deallocate(BmV_iTmp);deallocate(BC_iTmp)
		end if
		
		allocate(Teff_i(rowAge-1));allocate(L_i(rowAge-1))
		allocate(g_i(rowAge-1));allocate(rho_i(rowAge-1))
		
		Teff_i=10.**logTeff_i
		L_i=10.**logL_i
		g_i=10.**logg_i
		rho_i=10.**logrho_i
		
		deallocate(dist_ndxU)
		
		allocate(dsorted(size(dist)))
		dsorted=dist(:,1) !then it will be overwritten by sort so that it will be actually sorted
		call sort1(dsorted)
		allocate(ix_sorted(size(dist,1)))
		call indexx(dist(:,1), ix_sorted)
		!!Avoid extremely young isochrones (t<10 Myr) in the following calibration process
		age_s1=1 
		age_iso1=dist(ix_sorted(age_s1),2)
		do while (age_iso1.lt.logtlimCal)
			age_s1=age_s1+1
			age_iso1=dist(ix_sorted(age_s1),2)
		end do
		age_s=age_s1+1
		if (dsorted(age_s1)<1.e-6) then 
			do while ((log10(dsorted(age_s))-log10(dsorted(age_s1))<4.5.or.dist(ix_sorted(age_s),2).lt.logtlimCal) &
					& .and. age_s<size(ix_sorted))
				age_s=age_s+1
			end do 
		else
			do while (dist(ix_sorted(age_s),2).lt.logtlimCal .and. age_s<size(ix_sorted))
				age_s=age_s+1
			end do
		end if
		deallocate(dsorted)
		age_iso2=dist(ix_sorted(age_s),2)
		if (age_iso1>age_iso2) then 
			tmp_age1=age_iso1
			age_iso1=age_iso2
			age_iso2=tmp_age1 !age_iso2 BECOMES the isoc with the minimum distance from the star
			agePivot=age_iso2
		else
			agePivot=age_iso1 !age_iso1 REMAINS the isoc with the minimum distance from the star
		end if
		
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
		!!!!!!!!!!!END NEW CODE 
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
		
		call selectIsoc(logt_iso,age_iso1,logtstep,last_logt,ti0,tf0)
		ti(1)=ti0
		tf(1)=tf0
		call selectIsoc(logt_iso,age_iso2,logtstep,last_logt,ti0,tf0)
		ti(2)=ti0
		tf(2)=tf0
				
		if (calibHRD) then !!15/12/14
			VMag=V-DM0
			allocate(VMag_iso(size(Isoc,1)))
			VMag_iso=Isoc(:,cmag0)
			BmVuguali=.false.
			do jd=1,2 
				allocate(VMag_isor(tf(jd)-ti(jd)+1)); allocate(logL_isor(tf(jd)-ti(jd)+1))
				allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(Diff_BmV(tf(jd)-ti(jd)+1))
				VMag_isor=VMag_iso(ti(jd):tf(jd))
				logL_isor=logL_iso(ti(jd):tf(jd))
				BmV_isor=BmV_iso(ti(jd):tf(jd))
				Diff_BmV=sqrt((BmV-BmV_isor)**2+(VMag-VMag_isor)**2) !!!!!!!!!!!!18/2
				if (jd==1) then 
					id(jd)=minloc(Diff_BmV,1) 
					BmV_is1v(jd)=BmV_isor(id(jd))
				else
					allocate(ndxDv(size(Diff_BmV)))
					call indexx(Diff_BmV,ndxDv) 
					rowDv=size(ndxDv)
					allocate(VMagD(size(VMag_isor)))
					VMagD=VMag_isor(ndxDv)
					kD=1
					id(jd)=ndxDv(kD)
					BmV_is1v(jd)=BmV_isor(id(jd))
					do while ((abs(VMag_is1v(1)-VMagD(kD))<=0.01 .or. abs(VMag_is2v(1)-VMagD(kD))<=0.01) .and. &
							& kD<rowDv .and. VMag<2)!!!0.01 instead of 0.1!18/2 
						kD=kD+1
						id(jd)=ndxDv(kD)
						BmV_is1v(jd)=BmV_isor(id(jd))
					end do
					deallocate(ndxDv); deallocate(VMagD)
				end if 
				VMag_is1v(jd)=VMag_isor(id(jd))
				logL_is1v(jd)=logL_isor(id(jd))
				!!!
				call choose_i2xy1y2(BmV,BmV_isor,VMag_isor,logL_isor,id(jd),hstar,hstarlim, & !VMag,
					& BmV_is1v(jd),x_is2,VMag_is1v(jd),y1_is2,logL_is1v(jd),y2_is2)
				BmV_is2v(jd)=x_is2
				VMag_is2v(jd)=y1_is2
				logL_is2v(jd)=y2_is2
				!!!
												
				if (isEq(BmV_is1v(jd),BmV_is2v(jd),4)) then 	!!!if 10/2/14 -modified 7/11/2018
					BmVuguali=.true.
					VMag_isv(jd)=VMag
					call InterpLin_M(reshape((/VMag_is1v(jd),VMag_is2v(jd),logL_is1v(jd),logL_is2v(jd)/), &
						& (/2,2/)),VMag,1,1,(/2/),logL_isv(jd),xlow,ylow,xup,yup)
				else
					VMag_isv(jd)=(VMag_is2v(jd)-VMag_is1v(jd))/(BmV_is2v(jd)-BmV_is1v(jd))* &
							& (BmV-BmV_is1v(jd))+VMag_is1v(jd)
					logL_isv(jd)=(logL_is2v(jd)-logL_is1v(jd))/(BmV_is2v(jd)-BmV_is1v(jd))* &
							& (BmV-BmV_is1v(jd))+logL_is1v(jd)	
				end if 
								
				deallocate(VMag_isor); deallocate(logL_isor)
				deallocate(BmV_isor); deallocate(Diff_BmV)
			end do
						
			dCV_isv=abs(VMag_isv(1)-VMag_isv(2))  !vertical dist between isoc in the Color-Mag diagr
			dCL_isv=abs(logL_isv(1)-logL_isv(2))  !vertical dist between isoc in the Color-logL diagr
			if ((dCV_isv>dlim .and. dCL_isv>dlim)) then  !condition added 17/2/14
	!				dCVv=abs(VMag-VMag_isv)  		  !vertical dist between the star and the 2 isoc
				if (VMag<VMag_isv(1) .and. VMag<VMag_isv(2)) then 
					statev=1
				else if (VMag>VMag_isv(1) .and. VMag>VMag_isv(2)) then 
					statev=2
				else
					statev=3
				end if
				
				call calibrateDiag(statev,dCVvlim,VMag_isv,logL_isv,VMag,logL) 
				 
				L=10.**logL
				BC=-2.5*log10(L/d**2)-V-0.23
				I_L=L*(0.4*log(10.)*0.03+2.*I_d/d)
				I_logL=I_L/L*log10(e)
				logLuguali=.false.
				if (isEq(useColor,1.D0,2)) then
					do jd=1,2
						allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(logT_isor(tf(jd)-ti(jd)+1))
						allocate(logL_isor(tf(jd)-ti(jd)+1)); allocate(Diff_logL(tf(jd)-ti(jd)+1))
						
						BmV_isor=BmV_iso(ti(jd):tf(jd))
						logT_isor=logTeff_iso(ti(jd):tf(jd))
						logL_isor=logL_iso(ti(jd):tf(jd))
						Diff_logL=sqrt((logL-logL_isor)**2+(BmV-BmV_isor)**2) !!!!!!!!!!!!!!!!18/2/14
						if (jd==1) then
							id(jd)=minloc(Diff_logL,1) 
							logL_is1(jd)=logL_isor(id(jd))
						else
							allocate(ndxD(size(Diff_logL)))
							call indexx(Diff_logL,ndxD) 
							rowD=size(ndxD)
							allocate(BmVD(size(BmV_isor)))
							BmVD=BmV_isor(ndxD)
							kD=1
							id(jd)=ndxD(kD)
							logL_is1(jd)=logL_isor(id(jd))
							do while ((abs(BmV_is1(1)-BmVD(kD))<=0.04 .or. abs(BmV_is2(1)-BmVD(kD))<=0.04) .and. &
									& kD<rowD .and. logL>1.5) 
								kD=kD+1
								id(jd)=ndxD(kD)
								logL_is1(jd)=logL_isor(id(jd))
							end do
							deallocate(ndxD); deallocate(BmVD)
						end if 
						BmV_is1(jd)=BmV_isor(id(jd))
						logT_is1(jd)=logT_isor(id(jd))
						!!!
						call choose_i2xy1y2(logL,logL_isor,BmV_isor,logT_isor,id(jd),hstar,hstarlim, & !BmV,
							& logL_is1(jd),x_is2,BmV_is1(jd),y1_is2,logT_is1(jd),y2_is2)
						logL_is2(jd)=x_is2
						BmV_is2(jd)=y1_is2
						logT_is2(jd)=y2_is2
						!!!
													
						if (isEq(logL_is1(jd),logL_is2(jd),4)) then 
							logLuguali=.true.
							BmV_is(jd)=BmV
							call InterpLin_M(reshape((/BmV_is1(jd),BmV_is2(jd),logT_is1(jd),logT_is2(jd)/), &
								& (/2,2/)),BmV,1,1,(/2/),logT_is(jd),xlow,ylow,xup,yup)
						else
							BmV_is(jd)=(logL-logL_is1(jd))*(BmV_is2(jd)-BmV_is1(jd))/ &
										& (logL_is2(jd)-logL_is1(jd))+BmV_is1(jd)
							logT_is(jd)=(logL-logL_is1(jd))*(logT_is2(jd)-logT_is1(jd))/ &
										& (logL_is2(jd)-logL_is1(jd))+logT_is1(jd)
						end if 
											
						deallocate(BmV_isor); deallocate(logT_isor)
						deallocate(logL_isor); deallocate(Diff_logL)
					end do
										
					dCL_is=abs(BmV_is(1)-BmV_is(2))		!dist between isoc in the 'CMD'
					dHR_is=abs(logT_is(1)-logT_is(2))	!dist between isoc in the HRD
					if (dCL_is>dlim .and. dHR_is>dlim) then !condition added 17/2/14
		!					dCL=abs(BmV-BmV_is)				!vector of distances between the star and the 2 isoc
						if (BmV<BmV_is(1) .and. BmV<BmV_is(2)) then 
							state=1
						else if (BmV>BmV_is(1) .and. BmV>BmV_is(2)) then 
							state=2
						else
							state=3
						end if
						call calibrateDiag(state,dCLlim,BmV_is,logT_is,BmV,logTeff)
						
					else if (age_s==size(ix_sorted)) then 
						cycle
					end if
				else
					logTeff=BmV
					dCL_is=2.*dlim
					dHR_is=2.*dlim
				end if
				if (idCol.eq.1) then
					xtau=BmV
				else
					xtau=logTeff
				end if
			else if (age_s==size(ix_sorted)) then 
				cycle
			else
				dCL_is=0.
				dHR_is=0.
			end if
			cycleW=0
			do while ((dCV_isv<=dlim .or. dCL_isv<=dlim .or. dCL_is<=dlim .or. dHR_is<=dlim) .and. &
						& age_s<size(ix_sorted)) 
				age_s=age_s+1
				do while (dist(ix_sorted(age_s),2).lt.logtlimCal .and. age_s<size(ix_sorted))
					age_s=age_s+1
				end do
				if (isEq(age_iso1,agePivot,2)) then 
					age_iso2=dist(ix_sorted(age_s),2)
					if (age_iso1>age_iso2) then 
						tmp_age1=age_iso1
						age_iso1=age_iso2
						age_iso2=tmp_age1
						agePivot=age_iso2
					end if 
				else
					age_iso1=dist(ix_sorted(age_s),2)
					if (age_iso1>age_iso2) then 
						tmp_age1=age_iso1
						age_iso1=age_iso2
						age_iso2=tmp_age1
						agePivot=age_iso1
					end if 
				end if 
				
				call selectIsoc(logt_iso,age_iso1,logtstep,last_logt,ti0,tf0)
				ti(1)=ti0
				tf(1)=tf0
				call selectIsoc(logt_iso,age_iso2,logtstep,last_logt,ti0,tf0)
				ti(2)=ti0
				tf(2)=tf0
				
				VMag=V-DM0
				VMag_iso=Isoc(:,cmag0)
				BmVuguali=.false.
				do jd=1,2
					allocate(VMag_isor(tf(jd)-ti(jd)+1)); allocate(logL_isor(tf(jd)-ti(jd)+1))
					allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(Diff_BmV(tf(jd)-ti(jd)+1))					
					VMag_isor=VMag_iso(ti(jd):tf(jd))
					logL_isor=logL_iso(ti(jd):tf(jd))
					BmV_isor=BmV_iso(ti(jd):tf(jd))
					Diff_BmV=sqrt((BmV-BmV_isor)**2+(VMag-VMag_isor)**2) !!!!!!!!!!!!!!!!!18/2/14
					if (jd==1) then 
						id(jd)=minloc(Diff_BmV,1) 
						BmV_is1v(jd)=BmV_isor(id(jd))
					else
						allocate(ndxDv(size(Diff_BmV)))
						call indexx(Diff_BmV,ndxDv) 
						rowDv=size(ndxDv)
						allocate(VMagD(size(VMag_isor)))
						VMagD=VMag_isor(ndxDv)
						kD=1
						id(jd)=ndxDv(kD)
						BmV_is1v(jd)=BmV_isor(id(jd))
						do while ((abs(VMag_is1v(1)-VMagD(kD))<=0.01 .or. abs(VMag_is2v(1)-VMagD(kD))<=0.01) .and. &
								& kD<rowDv .and. VMag<2)!0.01 instead of 0.1 18/2 
							kD=kD+1
							id(jd)=ndxDv(kD)
							BmV_is1v(jd)=BmV_isor(id(jd))
						end do
						deallocate(ndxDv); deallocate(VMagD)
					end if 
					VMag_is1v(jd)=VMag_isor(id(jd))
					logL_is1v(jd)=logL_isor(id(jd))
					!!!
					call choose_i2xy1y2(BmV,BmV_isor,VMag_isor,logL_isor,id(jd),hstar,hstarlim, & !VMag,
						& BmV_is1v(jd),x_is2,VMag_is1v(jd),y1_is2,logL_is1v(jd),y2_is2) !!!
					BmV_is2v(jd)=x_is2
					VMag_is2v(jd)=y1_is2
					logL_is2v(jd)=y2_is2
					!!!
					if (isEq(BmV_is1v(jd),BmV_is2v(jd),4)) then 	!!!if 10/2/14 -modified 7/11/2018
						BmVuguali=.true.
						VMag_isv(jd)=VMag
						call InterpLin_M(reshape((/VMag_is1v(jd),VMag_is2v(jd),logL_is1v(jd),logL_is2v(jd)/), &
							& (/2,2/)),VMag,1,1,(/2/),logL_isv(jd),xlow,ylow,xup,yup)
					else
						VMag_isv(jd)=(VMag_is2v(jd)-VMag_is1v(jd))/(BmV_is2v(jd)-BmV_is1v(jd))* &
								& (BmV-BmV_is1v(jd))+VMag_is1v(jd)
						logL_isv(jd)=(logL_is2v(jd)-logL_is1v(jd))/(BmV_is2v(jd)-BmV_is1v(jd))* &
								& (BmV-BmV_is1v(jd))+logL_is1v(jd)	
					end if
					
					deallocate(VMag_isor); deallocate(logL_isor)
					deallocate(BmV_isor); deallocate(Diff_BmV)
				end do
								
				dCV_isv=abs(VMag_isv(1)-VMag_isv(2))	!vertical dist between isochrones in the CMD
				dCL_isv=abs(logL_isv(1)-logL_isv(2))	!vertical dist between isochrones in the Color-logL diagr
				if ((dCV_isv>dlim .and. dCL_isv>dlim)) then !condition added 17/2/14
	!					dCVv=abs(VMag-VMag_isv)					!vertical dist between the star and the 2 isoc
					if (VMag<VMag_isv(1) .and. VMag<VMag_isv(2)) then 
						statev=1
					else if (VMag>VMag_isv(1) .and. VMag>VMag_isv(2)) then 
						statev=2
					else
						statev=3
					end if
!					print*,'statev',statev
!					print*,'VMag_is1v 1st and 2nd isoc',VMag_is1v
!					print*,'VMag_is2v 1st and 2nd isoc',VMag_is2v
!					print*,'logL_is1v 1st and 2nd isoc',logL_is1v
!					print*,'logL_is2v 1st and 2nd isoc',logL_is2v
!					print*,'VMag_isv',VMag_isv
!					print*,'logL_isv',logL_isv
!					print*,'dCVvlim',dCVvlim
					call calibrateDiag(statev,dCVvlim,VMag_isv,logL_isv,VMag,logL)
					 
					L=10.**logL
					BC=-2.5*log10(L/d**2)-V-0.23
					I_L=L*(0.4*log(10.)*0.03+2.*I_d/d)
					I_logL=I_L/L*log10(e)
					logLuguali=.false.
					if (isEq(useColor,1.D0,2)) then
						do jd=1,2
							allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(logT_isor(tf(jd)-ti(jd)+1))
							allocate(logL_isor(tf(jd)-ti(jd)+1)); allocate(Diff_logL(tf(jd)-ti(jd)+1))
							logL_isor=logL_iso(ti(jd):tf(jd))
							logT_isor=logTeff_iso(ti(jd):tf(jd))
							BmV_isor=BmV_iso(ti(jd):tf(jd))
							Diff_logL=sqrt((logL-logL_isor)**2+(BmV-BmV_isor)**2) !!!!!!!!!!!!!!!18/2
							if (jd==1) then 
								id(jd)=minloc(Diff_logL,1) 
								logL_is1(jd)=logL_isor(id(jd))
							else
								allocate(ndxD(size(Diff_logL)))
								call indexx(Diff_logL,ndxD) 
								rowD=size(ndxD)
								allocate(BmVD(size(BmV_isor)))
								BmVD=BmV_isor(ndxD)
								kD=1
								id(jd)=ndxD(kD)
								logL_is1(jd)=logL_isor(id(jd))
								do while ((abs(BmV_is1(1)-BmVD(kD))<=0.04 .or. abs(BmV_is2(1)-BmVD(kD))<=0.04) .and. &
										& kD<rowD .and. logL>1.5) 
									kD=kD+1
									id(jd)=ndxD(kD)
									logL_is1(jd)=logL_isor(id(jd))
								end do
								deallocate(ndxD); deallocate(BmVD)
							end if 
							BmV_is1(jd)=BmV_isor(id(jd))
							logT_is1(jd)=logT_isor(id(jd))
							!!!
							call choose_i2xy1y2(logL,logL_isor,BmV_isor,logT_isor,id(jd),hstar,hstarlim, & !BmV,
								& logL_is1(jd),x_is2,BmV_is1(jd),y1_is2,logT_is1(jd),y2_is2)
							logL_is2(jd)=x_is2
							BmV_is2(jd)=y1_is2
							logT_is2(jd)=y2_is2
							!!!
							if (isEq(logL_is1(jd),logL_is2(jd),4)) then 
								logLuguali=.true.
								BmV_is(jd)=BmV
								call InterpLin_M(reshape((/BmV_is1(jd),BmV_is2(jd),logT_is1(jd),logT_is2(jd)/), &
									& (/2,2/)),BmV,1,1,(/2/),logT_is(jd),xlow,ylow,xup,yup)
							else
								BmV_is(jd)=(logL-logL_is1(jd))*(BmV_is2(jd)-BmV_is1(jd))/ &
											& (logL_is2(jd)-logL_is1(jd))+BmV_is1(jd)
								logT_is(jd)=(logL-logL_is1(jd))*(logT_is2(jd)-logT_is1(jd))/ &
											& (logL_is2(jd)-logL_is1(jd))+logT_is1(jd)
							end if 
												
							deallocate(BmV_isor); deallocate(logT_isor)
							deallocate(logL_isor); deallocate(Diff_logL)
						end do
												
						dCL_is=abs(BmV_is(1)-BmV_is(2))		!dist between isoc in the 'CMD'
						dHR_is=abs(logT_is(1)-logT_is(2))	!dist between isoc in the HRD
						if ((dCL_is>dlim .and. dHR_is>dlim)) then !condition added 17/2/14
		!						dCL=abs(BmV-BmV_is)					!vector of the distances between star-->2 isoc
							if (BmV<BmV_is(1) .and. BmV<BmV_is(2)) then 
								state=1
							else if (BmV>BmV_is(1) .and. BmV>BmV_is(2)) then 
								state=2
							else
								state=3
							end if
!							print*,'state',state
!							print*,'BmV_is1 1st and 2nd isoc',BmV_is1
!							print*,'BmV_is2 1st and 2nd isoc',BmV_is2
!							print*,'logT_is1 1st and 2nd isoc',logT_is1
!							print*,'logT_is2 1st and 2nd isoc',logT_is2
!							print*,'BmV_is',BmV_is
!							print*,'logT_is',logT_is
!							print*,'dCLlim',dCLlim
							call calibrateDiag(state,dCLlim,BmV_is,logT_is,BmV,logTeff)
														
						else if (age_s==size(ix_sorted)) then 
							cycleW=1
							cycle
						end if
					else
						logTeff=BmV
						dCL_is=2.*dlim
						dHR_is=2.*dlim
					end if
					if (idCol.eq.1) then
						xtau=BmV
					else
						xtau=logTeff
					end if
				else if (age_s==size(ix_sorted)) then 
					cycleW=1
					cycle
				end if 
			end do
			if (cycleW.eq.1) then
				cycle
			end if
			
			Teff=10.**logTeff
			if (.not.calibSPEC .and. Teffinput.ne.-1) then !Teffinput condtion added 14/7/2015
				if (abs(Teff-Teffinput)>300) then !Input B-V totally inconsistent with spectroscopic Teff
					cycle
				end if 
			end if
			if (isEq(useColor,1.D0,2)) then
				I_Teff=0.01*Teff
				I_logTeff=0.01*log10(e)
			else
				I_logTeff=InclogTeff
				I_Teff=Teff*log(10.)*I_logTeff
			end if
			Rc=sqrt(L/(Teff/TeffSun)**4)
			I_Rc=Rc*(0.5*I_L/L+2.*I_Teff/Teff)
			if (isEq(Rf,-1.D0,2)) then
				R1=Rc
				I_R1=I_Rc	
			else
				Ri=Rf
				I_Ri=SCP(22)
				call wMean(Ri,I_Ri,Rc,I_Rc,R1,I_R1)
			end if
			
			if (gProxyAvail) then  !!!14/7/2015; modified 18/9/2017
				loggAvail0=loggAvail.or.loggAvailAS
				rhoAvail0=rhoAvail.or.rhoAvailAS
				if ((loggAvail.or.loggAvailAS).and.(rhoAvail.or.rhoAvailAS)) then 
					!!!!!!Check compatibility between rho and logg 17/10/17 
					loggf0=loggf
					rhof0=rhof
					logg0=logg
					I_logg0=I_logg
					g0=g
					I_g0=I_g
					rho0=rho
					I_rho0=I_rho
					logrho0=logrho
					I_logrho0=I_logrho
					!!rho, logg both available. I determine input values for M, R 
					R2=g0/10.**loggSun*rhoSun/rho0
					I_R2=R2*(I_g0/g0+I_rho0/rho0)
					call wMean(R1,I_R1,R2,I_R2,R,I_R)
					
					M=(g0/10.**loggSun)**3*(rhoSun/rho0)**2
					I_M=M*(3.*I_g0/g0+2.*I_rho0/rho0)
!!					L=R**2*(Teff/TeffSun)**4
!!					I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
!!					logL=log10(L)
!!					I_logL=I_L/L*log10(e)
					!!!!!!Check compatibility between rho and logg 16/10/17 
					
					!!Modified to account for the mass uncertainty 26/10/2018 
					if (.not.(M+I_M.lt.minval(MTrAv) .or. M-I_M.gt.maxval(MTrAv))) then !open tracks only if inside MassRange 
						call selectMfromTracks(MTrAv,M,I_M,Mvec)
						allocate(cumInt(size(Mvec)))
						jk=0
						do jj=1,size(Mvec)
							indxM=minloc(abs(Mvec(jj)-MTrAv),1)
							xMTi=Tndxi(indxZt(indxZ),indxM)
							xMTf=Tndxf(indxZt(indxZ),indxM)
							allocate(TrRhoG(xMTf-xMTi+1,size(TrackTab,2)))
							TrRhoG=TrackTab(xMTi:xMTf,:)
		
							allocate(MTr(size(TrRhoG,1))); allocate(RTr(size(TrRhoG,1)))
							allocate(logRhoTr(size(TrRhoG,1))); allocate(loggTr(size(TrRhoG,1)))
							
							MTr=TrRhoG(:,cM_T) !Msun
							RTr=10.**TrRhoG(:,clogR_T) !cm
							logRhoTr=log10(MTr/(RTr/RSun)**3*rhoSun)
							loggTr=log10(MTr/(RTr/RSun)**2)+loggSun
							call findYgivenX_v(TrRhoG,logg,loggTr,clogTe_T,logTeTrtmp)
							if (allocated(logTeTrtmp)) then
								jk=jk+1
								if (jk.eq.1) then
									allocate(logTeTr(size(logTeTrtmp)))
									logTeTr=logTeTrtmp
								else
									call append1D(logTeTr,logTeTrtmp)
								end if
								cumInt(jj)=size(logTeTr)
								deallocate(logTeTrtmp)
							else
								cumInt(jj)=0
							end if
							
							deallocate(TrRhoG)
							deallocate(MTr); deallocate(RTr)
							deallocate(logRhoTr); deallocate(loggTr)
						end do
						!
						logTeJ=logTeff !I don't use Johnson (1966) relation because
						! I've just calibrate the correct isochronal Teff
						if (allocated(logTeTr)) then 
							allocate(TeB(size(logTeTr)))
							TeB=((logTeTr>logTeJ-DTeJ).and.(logTeTr<logTeJ+DTeJ))
							call consistentM(TeB,cumInt,Mvec,M,agreeM)
							if (.not.any(TeB) .or. .not.agreeM) then !all elements of TeB are 0 => logg and rho inconsistent
								!or even if there's some TeB, however (g,rho) no fully consistent
								R=R1 !ignore the R2
								I_R=I_R1
								if (I_logg>I_logrho) then !logrho is best determined
									!Only using rho: determine M; re-determine logg (discard the input value) 
									loggf=-1.
									loggAvail=.false.
									loggAvailAS=.false.
									if (.not.isEq(Rf,-1.D0,2)) then
										call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
									else
										call computeMg(R,I_R,Teff,I_Teff,rho,I_rho,M,I_M,g,I_g,logg,I_logg)
									end if
									loggAvailCal=.true. !!logg now available thanks to calibration
								else !logg is best determined
									!Only using logg: determine M; re-determine rho (discard the input value) 
									rhof=-1.
									rhoAvail=.false.
									rhoAvailAS=.false.
									if (.not.isEq(Rf,-1.D0,2)) then
										call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
									else
										call computeMrho(R,I_R,Teff,I_Teff,g,I_g,M,I_M,rho,I_rho,logrho,I_logrho)
									end if
									rhoAvailCal=.true.
								end if 
							else !input logg and rho are consistent. Compute M through a weighted mean between
								 ! M=M(logg) and M=M(rho)
								Mlogg=g/10.**loggSun*R**2
								I_Mlogg=Mlogg*(I_g/g+2.*I_R/R)
								Mrho=rho/rhoSun*R**3 !Mo
								I_Mrho=Mrho*(I_rho/rho+3.*I_R/R)
								call wMean(Mlogg,I_Mlogg,Mrho,I_Mrho,M,I_M)
							end if
							deallocate(TeB); deallocate(logTeTr)
						else !Track has been opened, but no logTe compatible with logg were found
							R=R1
							I_R=I_R1
							if (I_logg>I_logrho) then !logrho is best determined
								!Only using rho: determine M; re-determine logg (discard the input value) 
								loggf=-1.
								loggAvail=.false.
								loggAvailAS=.false.
								if (.not.isEq(Rf,-1.D0,2)) then
									call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
								else
									call computeMg(R,I_R,Teff,I_Teff,rho,I_rho,M,I_M,g,I_g,logg,I_logg)
								end if
								loggAvailCal=.true. !!logg now available thanks to calibration
							else !logg is best determined
								!Only using logg: determine M; re-determine rho (discard the input value) 
								rhof=-1.
								rhoAvail=.false.
								rhoAvailAS=.false.
								if (.not.isEq(Rf,-1.D0,2)) then
									call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
								else
									call computeMrho(R,I_R,Teff,I_Teff,g,I_g,M,I_M,rho,I_rho,logrho,I_logrho)
								end if
								rhoAvailCal=.true.
							end if 
						end if
						deallocate(cumInt); deallocate(Mvec)
					else !M is out of mass range of tracks
						R=R1
						I_R=I_R1
						if (I_logg>I_logrho) then !logrho is best determined
							!Only using rho: determine M; re-determine logg (discard the input value) 
							loggf=-1.
							loggAvail=.false.
							loggAvailAS=.false.
							if (.not.isEq(Rf,-1.D0,2)) then
								call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
							else
								call computeMg(R,I_R,Teff,I_Teff,rho,I_rho,M,I_M,g,I_g,logg,I_logg)
							end if
							loggAvailCal=.true. !!logg now available thanks to calibration
						else !logg is best determined
							!Only using logg: determine M; re-determine rho (discard the input value) 
							rhof=-1.
							rhoAvail=.false.
							rhoAvailAS=.false.
							if (.not.isEq(Rf,-1.D0,2)) then
								call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
							else
								call computeMrho(R,I_R,Teff,I_Teff,g,I_g,M,I_M,rho,I_rho,logrho,I_logrho)
							end if
							rhoAvailCal=.true.
						end if 
					end if 
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
				else !only one among logg and rho is present
					R=R1
					I_R=I_R1
					if (loggAvail.or.loggAvailAS) then 
						!	loggAvail0=1; 
						!	rhoAvail0=0;
						if (.not.isEq(Rf,-1.D0,2)) then
							call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
						else
							call computeMrho(R,I_R,Teff,I_Teff,g,I_g,M,I_M,rho,I_rho,logrho,I_logrho)
						end if
						rhoAvailCal=.true.
					else !rhoAvail since gProxyAvail
						!	loggAvail0=0; 
						!	rhoAvail0=1;
						if (.not.isEq(Rf,-1.D0,2)) then
							call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
						else
							call computeMg(R,I_R,Teff,I_Teff,rho,I_rho,M,I_M,g,I_g,logg,I_logg)
						end if
						loggAvailCal=.true. !!logg available thanks to calibration
					end if 
				end if !endif (loggAvail .and. rhoAvail)
			else
				R=R1
				I_R=I_R1
				if (.not.isEq(Rf,-1.D0,2)) then
					L=R**2*(Teff/TeffSun)**4
					I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
					logL=log10(L)
					I_logL=I_L/L*log10(e)
				end if
			end if !endif gProxy
		else if (calibNoD) then 
			allocate(logy_iso(size(Isoc,1)))
			if (.not.isEq(Rf,-1.D0,2)) then
				y=logL
				yl=-y
				logy_iso=logL_iso
				y_lim=-1.75
				dClim=dCLlim
			else if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
				& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))) then  !18/9/2017 - 30/10/18
				y=logg
				yl=y
				logy_iso=logg_iso
				y_lim=3.
				dClim=dCglim
			else !rhoAvail
				y=logrho
				yl=y
				logy_iso=logrho_iso
				y_lim=-2.2
				dClim=dCglim
			end if 
			logYuguali=.false.
			do jd=1,2
				allocate(logy_isor(tf(jd)-ti(jd)+1)); allocate(logT_isor(tf(jd)-ti(jd)+1))
				allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(Diff_logy(tf(jd)-ti(jd)+1))
				
				logy_isor=logy_iso(ti(jd):tf(jd))
				logT_isor=logTeff_iso(ti(jd):tf(jd))
				BmV_isor=BmV_iso(ti(jd):tf(jd))
				Diff_logy=sqrt((y-logy_isor)**2+(BmV-BmV_isor)**2) !!!!!!!!!!!!!!!!18/2/14
				if (jd==1) then 
					id(jd)=minloc(Diff_logy,1) 
					logy_is1(jd)=logy_isor(id(jd))
				else
					allocate(ndxD(size(Diff_logy)))
					call indexx(Diff_logy,ndxD) 
					rowD=size(ndxD)
					allocate(BmVD(size(BmV_isor)))
					BmVD=BmV_isor(ndxD)
					kD=1
					id(jd)=ndxD(kD)
					logy_is1(jd)=logy_isor(id(jd))
					do while ((abs(BmV_is1(1)-BmVD(kD))<=0.04 .or. abs(BmV_is2(1)-BmVD(kD))<=0.04) &
							& .and. kD<rowD .and. yl<y_lim) 
						kD=kD+1
						id(jd)=ndxD(kD)
						logy_is1(jd)=logy_isor(id(jd))
					end do
					deallocate(ndxD); deallocate(BmVD)
				end if 
				BmV_is1(jd)=BmV_isor(id(jd))
				logT_is1(jd)=logT_isor(id(jd))
				!!!
				call choose_i2xy1y2(y,logy_isor,BmV_isor,logT_isor,id(jd),hstar,hstarlim, & !BmV,
					& logy_is1(jd),x_is2,BmV_is1(jd),y1_is2,logT_is1(jd),y2_is2) !!!
				logy_is2(jd)=x_is2
				BmV_is2(jd)=y1_is2
				logT_is2(jd)=y2_is2
				!!!
				if (isEq(logy_is1(jd),logy_is2(jd),4)) then 
					logYuguali=.true.
					BmV_is(jd)=BmV
					call InterpLin_M(reshape((/BmV_is1(jd),BmV_is2(jd),logT_is1(jd),logT_is2(jd)/), &
						& (/2,2/)),BmV,1,1,(/2/),logT_is(jd),xlow,ylow,xup,yup)
				else
					BmV_is(jd)=(y-logy_is1(jd))*(BmV_is2(jd)-BmV_is1(jd))/(logy_is2(jd)-logy_is1(jd)) &
								& +BmV_is1(jd)
					logT_is(jd)=(y-logy_is1(jd))*(logT_is2(jd)-logT_is1(jd))/(logy_is2(jd)-logy_is1(jd)) &
								& +logT_is1(jd)
				end if
				
				deallocate(BmV_isor); deallocate(logy_isor);
				deallocate(logT_isor); deallocate(Diff_logy)
			end do 
								 
			dCL_is=abs(BmV_is(1)-BmV_is(2))		!dist between isoc in the BmV-logg Diagram
			dHR_is=abs(logT_is(1)-logT_is(2))	!dist between isoc in the logTeff-logg Diagram
			if ((dCL_is>dlim .and. dHR_is>dlim)) then !condition added 17/2/14
	!				dCL=abs(BmV-BmV_is)				!vector of distances between the star and the two isoc
				if (BmV<BmV_is(1) .and. BmV<BmV_is(2)) then 
					state=1
				else if (BmV>BmV_is(1) .and. BmV>BmV_is(2)) then 
					state=2
				else
					state=3
				end if
				!!
				call calibrateDiag(state,dClim,BmV_is,logT_is,BmV,logTeff)
				
				if (idCol.eq.1) then
					xtau=BmV
				else
					xtau=logTeff
				end if
				!! 
			else if (age_s==size(ix_sorted,1)) then 
				cycle
			else
				dCL_is=0.
				dHR_is=0.
			end if
			cycleW=0
			do while ((dCL_is<=dlim .or. dHR_is<=dlim) .and. age_s<size(ix_sorted,1)) 
				age_s=age_s+1
				do while (dist(ix_sorted(age_s),2).lt.logtlimCal .and. age_s<size(ix_sorted))
					age_s=age_s+1
				end do
				if (isEq(age_iso1,agePivot,2)) then 
					age_iso2=dist(ix_sorted(age_s),2)
					if (age_iso1>age_iso2) then 
						tmp_age1=age_iso1
						age_iso1=age_iso2
						age_iso2=tmp_age1
						agePivot=age_iso2
					end if 
				else
					age_iso1=dist(ix_sorted(age_s),2)
					if (age_iso1>age_iso2) then 
						tmp_age1=age_iso1
						age_iso1=age_iso2
						age_iso2=tmp_age1
						agePivot=age_iso1
					end if 
				end if 
				
				call selectIsoc(logt_iso,age_iso1,logtstep,last_logt,ti0,tf0)
				ti(1)=ti0
				tf(1)=tf0
				call selectIsoc(logt_iso,age_iso2,logtstep,last_logt,ti0,tf0)
				ti(2)=ti0
				tf(2)=tf0
				 
				logYuguali=.false.
				do jd=1,2
					allocate(logy_isor(tf(jd)-ti(jd)+1)); allocate(logT_isor(tf(jd)-ti(jd)+1))
					allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(Diff_logy(tf(jd)-ti(jd)+1))
				
					logy_isor=logy_iso(ti(jd):tf(jd))
					logT_isor=logTeff_iso(ti(jd):tf(jd))
					BmV_isor=BmV_iso(ti(jd):tf(jd))
					Diff_logy=sqrt((y-logy_isor)**2+(BmV-BmV_isor)**2) !!!!!!!!!!!!!!!!18/2/14
					if (jd==1) then 
						id(jd)=minloc(Diff_logy,1) 
						logy_is1(jd)=logy_isor(id(jd))
					else
						allocate(ndxD(size(Diff_logy)))
						call indexx(Diff_logy,ndxD) 
						rowD=size(ndxD)
						allocate(BmVD(size(BmV_isor)))
						BmVD=BmV_isor(ndxD)
						kD=1
						id(jd)=ndxD(kD)
						logy_is1(jd)=logy_isor(id(jd))
						do while ((abs(BmV_is1(1)-BmVD(kD))<=0.04 .or. abs(BmV_is2(1)-BmVD(kD))<=0.04) &
								& .and. kD<rowD .and. yl<y_lim) 
							kD=kD+1
							id(jd)=ndxD(kD)
							logy_is1(jd)=logy_isor(id(jd))
						end do
						deallocate(ndxD); deallocate(BmVD)
					end if 
					BmV_is1(jd)=BmV_isor(id(jd))
					logT_is1(jd)=logT_isor(id(jd))
					!!!
					call choose_i2xy1y2(y,logy_isor,BmV_isor,logT_isor,id(jd),hstar,hstarlim, & !BmV,
						& logy_is1(jd),x_is2,BmV_is1(jd),y1_is2,logT_is1(jd),y2_is2) !!!
					logy_is2(jd)=x_is2
					BmV_is2(jd)=y1_is2
					logT_is2(jd)=y2_is2
					!!!
					if (isEq(logy_is1(jd),logy_is2(jd),4)) then 
						logYuguali=.true.
						BmV_is(jd)=BmV
						call InterpLin_M(reshape((/BmV_is1(jd),BmV_is2(jd),logT_is1(jd),logT_is2(jd)/), &
							& (/2,2/)),BmV,1,1,(/2/),logT_is(jd),xlow,ylow,xup,yup)
					else
						BmV_is(jd)=(y-logy_is1(jd))*(BmV_is2(jd)-BmV_is1(jd))/(logy_is2(jd)-logy_is1(jd)) &
									& +BmV_is1(jd)
						logT_is(jd)=(y-logy_is1(jd))*(logT_is2(jd)-logT_is1(jd))/(logy_is2(jd)-logy_is1(jd)) &
									& +logT_is1(jd)
					end if 
				
					deallocate(BmV_isor); deallocate(logy_isor);
					deallocate(logT_isor); deallocate(Diff_logy)
				end do 
				
				dCL_is=abs(BmV_is(1)-BmV_is(2))		!dist between isoc in the BmV-logg Diagram
				dHR_is=abs(logT_is(1)-logT_is(2))	!dist between isoc in the logTeff-logg Diagram
				if ((dCL_is>dlim .and. dHR_is>dlim)) then !condition added 17/2/14
	!				dCL=abs(BmV-BmV_is)				!vector of distances between the star and the two isoc
					if (BmV<BmV_is(1) .and. BmV<BmV_is(2)) then 
						state=1
					else if (BmV>BmV_is(1) .and. BmV>BmV_is(2)) then 
						state=2
					else
						state=3
					end if
					!!
					call calibrateDiag(state,dClim,BmV_is,logT_is,BmV,logTeff)
					if (idCol.eq.1) then
						xtau=BmV
					else
						xtau=logTeff
					end if
					!! 
				else if (age_s==size(ix_sorted,1)) then 
					cycleW=1
					cycle
				end if 
			end do
			if (cycleW.eq.1) then
				cycle
			end if
			Teff=10.**logTeff
			if ((.not.calibSPEC .and. Teffinput.ne.-1)) then !Inside calibNoD => .not.calibSPEC is always true
				if (abs(Teff-Teffinput)>300) then !Inconsistency between BmV and spectrscopic Teff
					cycle
				end if 
			end if 
			I_Teff=0.01*Teff
			I_logTeff=0.01*log10(e)
			!!Compatibility
			loggAvail0=loggAvail.or.loggAvailAS
			rhoAvail0=rhoAvail.or.rhoAvailAS
			if ((loggAvail.or.loggAvailAS).and.(rhoAvail.or.rhoAvailAS)) then 
				!!rho, logg both available. I determine input values for M, R 
				loggf0=loggf
				rhof0=rhof
				logg0=logg
				I_logg0=I_logg
				g0=g
				I_g0=I_g
				rho0=rho
				I_rho0=I_rho
				logrho0=logrho
				I_logrho0=I_logrho
				R2=g/10.**loggSun*rhoSun/rho
				I_R2=R2*(I_g/g+I_rho/rho)
				if (isEq(Rf,-1.D0,2)) then
					R=R2
					I_R=I_R2
				else
					call wMean(R1,I_R1,R2,I_R2,R,I_R)
				end if
				M=(g/10.**loggSun)**3*(rhoSun/rho)**2
				I_M=M*(3.*I_g/g+2.*I_rho/rho)
				L=R**2*(Teff/TeffSun)**4
				I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
				logL=log10(L)
				I_logL=I_L/L*log10(e)
				!!!!!!Check compatibility between rho and logg 16/10/17 
				!!Modified to account for the mass uncertainty 26/10/2018 
				if (.not.(M+I_M.lt.minval(MTrAv) .or. M-I_M.gt.maxval(MTrAv))) then !open tracks only if inside MassRange 
					call selectMfromTracks(MTrAv,M,I_M,Mvec)
					allocate(cumInt(size(Mvec)))
					jk=0
					do jj=1,size(Mvec)
						indxM=minloc(abs(Mvec(jj)-MTrAv),1)
						xMTi=Tndxi(indxZt(indxZ),indxM)
						xMTf=Tndxf(indxZt(indxZ),indxM)
						allocate(TrRhoG(xMTf-xMTi+1,size(TrackTab,2)))
						TrRhoG=TrackTab(xMTi:xMTf,:)
	
						allocate(MTr(size(TrRhoG,1))); allocate(RTr(size(TrRhoG,1)))
						allocate(logRhoTr(size(TrRhoG,1))); allocate(loggTr(size(TrRhoG,1)))
						
						MTr=TrRhoG(:,cM_T) !Msun
						RTr=10.**TrRhoG(:,clogR_T) !cm
						logRhoTr=log10(MTr/(RTr/RSun)**3*rhoSun)
						loggTr=log10(MTr/(RTr/RSun)**2)+loggSun
						call findYgivenX_v(TrRhoG,logg,loggTr,clogTe_T,logTeTrtmp)
						if (allocated(logTeTrtmp)) then
							jk=jk+1
							if (jk.eq.1) then
								allocate(logTeTr(size(logTeTrtmp)))
								logTeTr=logTeTrtmp
							else
								call append1D(logTeTr,logTeTrtmp)
							end if
							cumInt(jj)=size(logTeTr)
							deallocate(logTeTrtmp)
						else
							cumInt(jj)=0
						end if
						
						deallocate(TrRhoG)
						deallocate(MTr); deallocate(RTr)
						deallocate(logRhoTr); deallocate(loggTr)
					end do
					logTeJ=logTeff !I don't use Johnson (1966) relation because I've just calibrated Teff
					DTeJ=I_logTeff !directly set to the true uncertainty
					if (allocated(logTeTr)) then 
						allocate(TeB(size(logTeTr)))
						TeB=((logTeTr>logTeJ-DTeJ).and.(logTeTr<logTeJ+DTeJ))
						if (.not.any(TeB)) then !all elements of TeB are 0 => logg and rho inconsistent
							if (I_logg>I_logrho) then !logrho is best determined
								if (isEq(Rf,-1.D0,2)) then
									call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
									call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
									loggf=-1.
									call setNaN(logg); call setNaN(I_logg)
									call setNaN(g); call setNaN(I_g)
									loggAvail=.false.
									loggAvailAS=.false.
								else !input R is available
									!Only using rho: determine M; re-determine logg (discard the input value) 
									loggf=-1.
									loggAvail=.false.
									loggAvailAS=.false.
									R=R1 !use just the input reliable value
									I_R=I_R1
									call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
									loggAvailCal=.true. !!logg now available thanks to calibration
								end if
							else !logg is best determined
								if (isEq(Rf,-1.D0,2)) then
									call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
									call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
									rhof=-1.
									call setNaN(rho); call setNaN(I_rho)
									call setNaN(logrho); call setNaN(I_logrho)
									rhoAvail=.false.
									rhoAvailAS=.false.
								else !input R is available
									!Only using logg: determine M; re-determine rho (discard the input value) 
									rhof=-1.
									rhoAvail=.false.
									rhoAvailAS=.false.
									R=R1 !use just the input reliable value
									I_R=I_R1
									call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
									rhoAvailCal=.true.
								end if
							end if 
						else
							!!!!!!!!!
							call consistentM(TeB,cumInt,Mvec,M,agreeM)
							if (.not.agreeM) then !new M has been recomputed inside the subroutine in case agreeM is false
											      !This M will be used in case R is not available directly from input
								if (I_logg.gt.I_logrho) then !logrho is best determined
									if (isEq(Rf,-1.D0,2)) then
										loggAvail=.false.
										loggAvailAS=.false.
										!g: mantain original input uncertainty
										call computeLRgfromMrho(M,I_M,rho,I_rho,Teff,I_Teff,L,logL,I_L,I_logL,R,I_R,g,logg)
										loggAvailCal=.true.
									else !input R is available
										!Only using rho: determine M; re-determine logg (discard the input value) 
										loggf=-1.
										loggAvail=.false.
										loggAvailAS=.false.
										R=R1
										I_R=I_R1
										call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
										loggAvailCal=.true. !!logg now available thanks to calibration
									end if
								else
									if (isEq(Rf,-1.D0,2)) then
										rhoAvail=.false.
										rhoAvailAS=.false.
										!rho in g/cm3. Mantain original input uncertainty on rho
										call computeLRrhofromMg(M,I_M,g,I_g,Teff,I_Teff,L,logL,I_L,I_logL,R,I_R,rho,logrho)
										rhoAvailCal=.true.
									else !input R is available
										!Only using logg: determine M; re-determine rho (discard the input value) 
										rhof=-1.
										rhoAvail=.false.
										rhoAvailAS=.false.
										R=R1
										I_R=I_R1
										call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
										rhoAvailCal=.true.
									end if
								end if
								if (isEq(Rf,-1.D0,2)) then
									R=g/10**loggSun*rhoSun/rho
									I_R=R*(I_g/g+I_rho/rho)
									L=R**2*(Teff/TeffSun)**4
									I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
									logL=log10(L)
									I_logL=I_L/L*log10(e)
								end if
							else
								if (I_logg.lt.I_logrho) then
									bestLogg=.true.
								else
									bestLogg=.false.
								end if
								if (.not.isEq(Rf,-1.D0,2)) then
								!input R is available: weighted mean to infer M
									Mlogg=g/10.**loggSun*R**2
									I_Mlogg=Mlogg*(I_g/g+2.*I_R/R)
									Mrho=rho/rhoSun*R**3 !Mo
									I_Mrho=Mrho*(I_rho/rho+3.*I_R/R)
									call wMean(Mlogg,I_Mlogg,Mrho,I_Mrho,M,I_M)
									!L already defined since the beginning from R
								end if
							end if
							!!!!!!!!!!!!!
						end if
						deallocate(logTeTr); deallocate(TeB)
					else !Track has been opened, but no logTe compatible with logg were found
						if (I_logg>I_logrho) then !logrho is best determined
							if (isEq(Rf,-1.D0,2)) then
								call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
								call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
								loggf=-1.
								call setNaN(logg); call setNaN(I_logg)
								call setNaN(g); call setNaN(I_g)
								loggAvail=.false.
								loggAvailAS=.false.
							else !input R is available
								!Only using rho: determine M; re-determine logg (discard the input value) 
								loggf=-1.
								loggAvail=.false.
								loggAvailAS=.false.
								R=R1 !use just the input reliable value
								I_R=I_R1
								call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
								loggAvailCal=.true. !!logg now available thanks to calibration
							end if
						else !logg is best determined
							if (isEq(Rf,-1.D0,2)) then
								call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
								call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
								rhof=-1.
								call setNaN(rho); call setNaN(I_rho)
								call setNaN(logrho); call setNaN(I_logrho)
								rhoAvail=.false.
								rhoAvailAS=.false.
							else !input R is available
								!Only using logg: determine M; re-determine rho (discard the input value) 
								rhof=-1.
								rhoAvail=.false.
								rhoAvailAS=.false.
								R=R1 !use just the input reliable value
								I_R=I_R1
								call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
								rhoAvailCal=.true.
							end if
						end if 
					end if
					deallocate(cumInt) 
				else
					if (I_logg>I_logrho) then !logrho is best determined
						if (isEq(Rf,-1.D0,2)) then
							call setNan(R); call setNan(I_R); call setNan(M); call setNan(I_M)
							call setNan(L); call setNan(I_L); call setNan(logL); call setNan(I_logL)
							loggf=-1.
							call setNan(logg); call setNan(I_logg)
							call setNan(g); call setNan(I_g)
							loggAvail=.false.
							loggAvailAS=.false.
						else !input R is available
							!Only using rho: determine M; re-determine logg (discard the input value) 
							loggf=-1.
							loggAvail=.false.
							loggAvailAS=.false.
							R=R1 !use just the input reliable value
							I_R=I_R1
							call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
							loggAvailCal=.true. !!logg now available thanks to calibration
						end if
					else !logg is best determined
						if (isEq(Rf,-1.D0,2)) then
							call setNan(R); call setNan(I_R); call setNan(M); call setNan(I_M)
							call setNan(L); call setNan(I_L); call setNan(logL); call setNan(I_logL)
							rhof=-1.
							call setNan(rho); call setNan(I_rho)
							call setNan(logrho); call setNan(I_logrho)
							rhoAvail=.false.
							rhoAvailAS=.false.
						else !input R is available
							!Only using logg: determine M; re-determine rho (discard the input value) 
							rhof=-1.
							rhoAvail=.false.
							rhoAvailAS=.false.
							R=R1 !use just the input reliable value
							I_R=I_R1
							call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
							rhoAvailCal=.true.
						end if
					end if 
				end if 
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			else
				if (.not.isEq(Rf,-1.D0,2)) then
					R=R1
					I_R=I_R1
					if (loggAvail.or.loggAvailAS) then
						call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
						rhoAvailCal=.true.
					else if (rhoAvail.or.rhoAvailAS) then
						call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
						loggAvailCal=.true. !!logg available thanks to calibration
					else !noGproxy
						L=R**2*(Teff/TeffSun)**4
						I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
						logL=log10(L)
						I_logL=I_L/L*log10(e)
					end if
				end if
			end if
			!!End compatibility
			if (.not.isEq(Rf,-1.D0,2)) then
				kwl=0
				cycleW=0
				age_s=1
				do while (abs(logL-y).gt.DlogL.and.kwl.lt.10.and.age_s<size(ix_sorted))
					y=logL
					kwl=kwl+1
					deallocate(t_i); deallocate(rho_i); deallocate(dist);deallocate(ix_sorted)
					deallocate(logTeff_i); deallocate(logL_i); deallocate(M_i); deallocate(logg_i)
					deallocate(logrho_i); deallocate(Teff_i); deallocate(g_i); deallocate(L_i)
					deallocate(logy_iso)
					if (photIsocAvail) then
						deallocate(BC_i); deallocate(BmV_i)
					end if
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
					!!!!!!!NEW CODE 12/2/14!!!!!!!!!!!!!!! 
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
					allocate(dist_ndx(size(logTeff_iso))) !equal to any other size of _iso
					allocate(x_iso(size(logTeff_iso))); allocate(xx_iso(size(logTeff_iso)))
					call indexx(sqrt((BmV-BmV_iso)**2+(logL-logL_iso)**2),dist_ndx)
					x=BmV
					x_iso=BmV_iso
					xx_iso=logTeff_iso
					cBrif=cBBmV
					allocate(ageMinDist(size(logt_iso)))
					ageMinDist=logt_iso(dist_ndx)
					
					call uniqueFast(ageMinDist,2,ix,.true.)
					call sort0(ix) !so to display elements of ageMinDist in the oroginal order
					allocate(dist_ndxU(size(ix))) 
					dist_ndxU=dist_ndx(ix)
					
					deallocate(dist_ndx); deallocate(ageMinDist)
					deallocate(ix)
					rowAge=1
					! allocation of already defined variables
					allocate(t_iTmp(size(dist_ndxU)));allocate(distTmp(size(dist_ndxU),2))
					allocate(logTeff_iTmp(size(dist_ndxU)));allocate(logL_iTmp(size(dist_ndxU)))
					allocate(M_iTmp(size(dist_ndxU)));allocate(logg_iTmp(size(dist_ndxU)))
					allocate(logrho_iTmp(size(dist_ndxU)))
					if (photIsocAvail) then 
						allocate(BmV_iTmp(size(dist_ndxU)));allocate(BC_iTmp(size(dist_ndxU)))
					end if
					!
					do Nage=1,size(dist_ndxU) 
						ndxDU=dist_ndxU(Nage)
						if (photIsocAvail) then 
							BmV_i1=BmV_iso(ndxDU)
							V_i1=V_iso(ndxDU)
						end if 
						logTeff_i1=logTeff_iso(ndxDU)
						logL_i1=logL_iso(ndxDU)
						M_i1=M_iso(ndxDU)
						logg_i1=logg_iso(ndxDU)
						logrho_i1=logrho_iso(ndxDU)
						
						call choose_i2Col(ndxDU,x,x_iso,xx_iso,calibSPEC,hstar,hstarlim,V_iso,logL_iso, &
							& M_iso,logg_iso,logrho_iso,logt_iso,BmV_i2,V_i2,logTeff_i2,logL_i2,M_i2,logg_i2,logrho_i2)
						if (calibSPEC) then
							x_i1=logTeff_i1
							x_i2=logTeff_i2
						else
							x_i1=BmV_i1
							x_i2=BmV_i2
						end if
						if (isEq(x_i1,x_i2,4)) then 
							cycle
						end if 
						t_iTmp(rowAge)=t_iso(ndxDU)
						if (photIsocAvail) then 
							BC_iTmp(rowAge)=Isoc(ndxDU,cmbol)-Isoc(ndxDU,cmag0) !just rough.
							!When I recover BC in the logg-BmV plane, I'll make the interpolation using BC_i1 e BC_i2
						end if 
						
						call findX_i(BmV,logL,BmV_i1,BmV_i2,logL_i1,logL_i2,x_i,dist0)
						BmV_iTmp(rowAge)=x_i 
						
						distTmp(rowAge,:)=(/ dist0, logt_iso(ndxDU) /)
									
						!!18/3/15 The following if cancel the criterion of perpendicular distance
						!!between a star and the isochrone if the theoretical point to be interpolated
						!!doesn't fall inside the isochrone (i.e. between x_i1 and x_i2), but on its
						!!extension. This possibility, in fact, would build fictitious isochrones that
						!!could be erroneously close to the star
						if (.not.((x_i>=x_i1.and.x_i<=x_i2).or.(x_i>=x_i2.and.x_i<=x_i1))) then 
							if (photIsocAvail) then 
								BmV_iTmp(rowAge)=BmV_i1
							end if 
							logg_iTmp(rowAge)=logg_i1
							logTeff_iTmp(rowAge)=logTeff_i1
							logL_iTmp(rowAge)=logL_i1
							M_iTmp(rowAge)=M_i1
							logrho_iTmp(rowAge)=logrho_i1
						else
							mT=(logTeff_i2-logTeff_i1)/(x_i2-x_i1)
							qT=-mT*x_i1+logTeff_i1
							logTeff_iTmp(rowAge)=mT*x_i+qT
							
							mg=(logg_i2-logg_i1)/(x_i2-x_i1)
							qg=-mg*x_i1+logg_i1
							logg_iTmp(rowAge)=mg*x_i+qg
							mL=(logL_i2-logL_i1)/(x_i2-x_i1)
							qL=-mL*x_i1+logL_i1
							logL_iTmp(rowAge)=mL*x_i+qL
							mM=(M_i2-M_i1)/(x_i2-x_i1)
							qM=-mM*x_i1+M_i1
							M_iTmp(rowAge)=mM*x_i+qM
							mrh=(logrho_i2-logrho_i1)/(x_i2-x_i1)
							qrho=-mrh*x_i1+logrho_i1
							logrho_iTmp(rowAge)=mrh*x_i+qrho
						end if 

						rowAge=rowAge+1
					end do
					
					deallocate(x_iso); deallocate(xx_iso)
							
					allocate(t_i(rowAge-1));allocate(dist(rowAge-1,2))
					allocate(logTeff_i(rowAge-1));allocate(logL_i(rowAge-1))
					allocate(M_i(rowAge-1));allocate(logg_i(rowAge-1))
					allocate(logrho_i(rowAge-1))
					allocate(BmV_i(rowAge-1));allocate(BC_i(rowAge-1))
					t_i=t_iTmp(1:rowAge-1)
					dist=distTmp(1:rowAge-1,:)
					logTeff_i=logTeff_iTmp(1:rowAge-1)
					logL_i=logL_iTmp(1:rowAge-1)
					M_i=M_iTmp(1:rowAge-1)
					logg_i=logg_iTmp(1:rowAge-1)
					logrho_i=logrho_iTmp(1:rowAge-1)
					BmV_i=BmV_iTmp(1:rowAge-1)
					BC_i=BC_iTmp(1:rowAge-1)
					deallocate(t_iTmp);deallocate(distTmp);deallocate(logTeff_iTmp)
					deallocate(logL_iTmp);deallocate(M_iTmp);deallocate(logg_iTmp)
					deallocate(logrho_iTmp)
					if (photIsocAvail) then
						deallocate(BmV_iTmp);deallocate(BC_iTmp)
					end if
					
					allocate(Teff_i(rowAge-1));allocate(L_i(rowAge-1))
					allocate(g_i(rowAge-1));allocate(rho_i(rowAge-1))
					
					Teff_i=10.**logTeff_i
					L_i=10.**logL_i
					g_i=10.**logg_i
					rho_i=10.**logrho_i
					
					deallocate(dist_ndxU)
					
					allocate(dsorted(size(dist)))
					dsorted=dist(:,1) !then it will be overwritten by sort so that it will be actually sorted
					call sort1(dsorted)
					allocate(ix_sorted(size(dist,1)))
					call indexx(dist(:,1), ix_sorted) 
					!!Avoid extremely young isochrones (t<10 Myr) in the following calibration process
					age_s1=1 
					age_iso1=dist(ix_sorted(age_s1),2)
					do while (age_iso1.lt.logtlimCal)
						age_s1=age_s1+1
						age_iso1=dist(ix_sorted(age_s1),2)
					end do
					age_s=age_s1+1
					if (dsorted(age_s1)<1.e-6) then 
						do while ((log10(dsorted(age_s))-log10(dsorted(age_s1))<4.5.or.dist(ix_sorted(age_s),2).lt.logtlimCal) &
								& .and. age_s<size(ix_sorted))
							age_s=age_s+1
						end do 
					else
						do while (dist(ix_sorted(age_s),2).lt.logtlimCal .and. age_s<size(ix_sorted))
							age_s=age_s+1
						end do
					end if
					deallocate(dsorted)
					age_iso2=dist(ix_sorted(age_s),2)
					if (age_iso1>age_iso2) then 
						tmp_age1=age_iso1
						age_iso1=age_iso2
						age_iso2=tmp_age1 !age_iso2 BECOMES the isoc with the minimum distance from the star
						agePivot=age_iso2
					else
						agePivot=age_iso1 !age_iso1 REMAINS the isoc with the minimum distance from the star
					end if
					
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
					!!!!!!!!!!!END NEW CODE 
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
					
					call selectIsoc(logt_iso,age_iso1,logtstep,last_logt,ti0,tf0)
					ti(1)=ti0
					tf(1)=tf0
					call selectIsoc(logt_iso,age_iso2,logtstep,last_logt,ti0,tf0)
					ti(2)=ti0
					tf(2)=tf0
										
					allocate(logy_iso(size(Isoc,1)))
					
					y=logL !already stated before
					yl=-y
					logy_iso=logL_iso
					y_lim=-1.75
					
					logYuguali=.false.
					do jd=1,2
						allocate(logy_isor(tf(jd)-ti(jd)+1)); allocate(logT_isor(tf(jd)-ti(jd)+1))
						allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(Diff_logy(tf(jd)-ti(jd)+1))
						
						logy_isor=logy_iso(ti(jd):tf(jd))
						logT_isor=logTeff_iso(ti(jd):tf(jd))
						BmV_isor=BmV_iso(ti(jd):tf(jd))
						Diff_logy=sqrt((y-logy_isor)**2+(BmV-BmV_isor)**2) !!!!!!!!!!!!!!!!18/2/14
						if (jd==1) then 
							id(jd)=minloc(Diff_logy,1) 
							logy_is1(jd)=logy_isor(id(jd))
						else
							allocate(ndxD(size(Diff_logy)))
							call indexx(Diff_logy,ndxD) 
							rowD=size(ndxD)
							allocate(BmVD(size(BmV_isor)))
							BmVD=BmV_isor(ndxD)
							kD=1
							id(jd)=ndxD(kD)
							logy_is1(jd)=logy_isor(id(jd))
							do while ((abs(BmV_is1(1)-BmVD(kD))<=0.04 .or. abs(BmV_is2(1)-BmVD(kD))<=0.04) &
									& .and. kD<rowD .and. yl<y_lim) 
								kD=kD+1
								id(jd)=ndxD(kD)
								logy_is1(jd)=logy_isor(id(jd))
							end do
							deallocate(ndxD); deallocate(BmVD)
						end if 
						BmV_is1(jd)=BmV_isor(id(jd))
						logT_is1(jd)=logT_isor(id(jd))
						!!!
						call choose_i2xy1y2(y,logy_isor,BmV_isor,logT_isor,id(jd),hstar,hstarlim, & !BmV,
							& logy_is1(jd),x_is2,BmV_is1(jd),y1_is2,logT_is1(jd),y2_is2) !!!
						logy_is2(jd)=x_is2
						BmV_is2(jd)=y1_is2
						logT_is2(jd)=y2_is2
						!!!
						if (isEq(logy_is1(jd),logy_is2(jd),4)) then 
							logYuguali=.true.
							BmV_is(jd)=BmV
							call InterpLin_M(reshape((/BmV_is1(jd),BmV_is2(jd),logT_is1(jd),logT_is2(jd)/), &
								& (/2,2/)),BmV,1,1,(/2/),logT_is(jd),xlow,ylow,xup,yup)
						else
							BmV_is(jd)=(y-logy_is1(jd))*(BmV_is2(jd)-BmV_is1(jd))/(logy_is2(jd)-logy_is1(jd)) &
										& +BmV_is1(jd)
							logT_is(jd)=(y-logy_is1(jd))*(logT_is2(jd)-logT_is1(jd))/(logy_is2(jd)-logy_is1(jd)) &
										& +logT_is1(jd)
						end if
						
						deallocate(BmV_isor); deallocate(logy_isor);
						deallocate(logT_isor); deallocate(Diff_logy)
					end do 
												 
					dCL_is=abs(BmV_is(1)-BmV_is(2))		!dist between isoc in the BmV-logg Diagram
					dHR_is=abs(logT_is(1)-logT_is(2))	!dist between isoc in the logTeff-logg Diagram
					if ((dCL_is>dlim .and. dHR_is>dlim)) then !condition added 17/2/14
			!				dCL=abs(BmV-BmV_is)				!vector of distances between the star and the two isoc
						if (BmV<BmV_is(1) .and. BmV<BmV_is(2)) then 
							state=1
						else if (BmV>BmV_is(1) .and. BmV>BmV_is(2)) then 
							state=2
						else
							state=3
						end if
						!!
						call calibrateDiag(state,dClim,BmV_is,logT_is,BmV,logTeff)
						if (idCol.eq.1) then
							xtau=BmV
						else
							xtau=logTeff
						end if
						!! 
					else if (age_s==size(ix_sorted,1)) then 
						cycleW=1
						cycle
					else
						dCL_is=0.
						dHR_is=0.
					end if 
					do while ((dCL_is<=dlim .or. dHR_is<=dlim) .and. age_s<size(ix_sorted,1)) 
						age_s=age_s+1
						do while (dist(ix_sorted(age_s),2).lt.logtlimCal .and. age_s<size(ix_sorted))
							age_s=age_s+1
						end do
						if (isEq(age_iso1,agePivot,2)) then 
							age_iso2=dist(ix_sorted(age_s),2)
							if (age_iso1>age_iso2) then 
								tmp_age1=age_iso1
								age_iso1=age_iso2
								age_iso2=tmp_age1
								agePivot=age_iso2
							end if 
						else
							age_iso1=dist(ix_sorted(age_s),2)
							if (age_iso1>age_iso2) then 
								tmp_age1=age_iso1
								age_iso1=age_iso2
								age_iso2=tmp_age1
								agePivot=age_iso1
							end if 
						end if 
						
						call selectIsoc(logt_iso,age_iso1,logtstep,last_logt,ti0,tf0)
						ti(1)=ti0
						tf(1)=tf0
						call selectIsoc(logt_iso,age_iso2,logtstep,last_logt,ti0,tf0)
						ti(2)=ti0
						tf(2)=tf0
						 
						logYuguali=.false.
						do jd=1,2
							allocate(logy_isor(tf(jd)-ti(jd)+1)); allocate(logT_isor(tf(jd)-ti(jd)+1))
							allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(Diff_logy(tf(jd)-ti(jd)+1))
						
							logy_isor=logy_iso(ti(jd):tf(jd))
							logT_isor=logTeff_iso(ti(jd):tf(jd))
							BmV_isor=BmV_iso(ti(jd):tf(jd))
							Diff_logy=sqrt((y-logy_isor)**2+(BmV-BmV_isor)**2) !!!!!!!!!!!!!!!!18/2/14
							if (jd==1) then 
								id(jd)=minloc(Diff_logy,1) 
								logy_is1(jd)=logy_isor(id(jd))
							else
								allocate(ndxD(size(Diff_logy)))
								call indexx(Diff_logy,ndxD) 
								rowD=size(ndxD)
								allocate(BmVD(size(BmV_isor)))
								BmVD=BmV_isor(ndxD)
								kD=1
								id(jd)=ndxD(kD)
								logy_is1(jd)=logy_isor(id(jd))
								do while ((abs(BmV_is1(1)-BmVD(kD))<=0.04 .or. abs(BmV_is2(1)-BmVD(kD))<=0.04) &
										& .and. kD<rowD .and. yl<y_lim) 
									kD=kD+1
									id(jd)=ndxD(kD)
									logy_is1(jd)=logy_isor(id(jd))
								end do
								deallocate(ndxD); deallocate(BmVD)
							end if 
							BmV_is1(jd)=BmV_isor(id(jd))
							logT_is1(jd)=logT_isor(id(jd))
							!!!
							call choose_i2xy1y2(y,logy_isor,BmV_isor,logT_isor,id(jd),hstar,hstarlim, & !BmV,
								& logy_is1(jd),x_is2,BmV_is1(jd),y1_is2,logT_is1(jd),y2_is2) !!!
							logy_is2(jd)=x_is2
							BmV_is2(jd)=y1_is2
							logT_is2(jd)=y2_is2
							!!!
							if (isEq(logy_is1(jd),logy_is2(jd),4)) then 
								logYuguali=.true.
								BmV_is(jd)=BmV
								call InterpLin_M(reshape((/BmV_is1(jd),BmV_is2(jd),logT_is1(jd),logT_is2(jd)/), &
									& (/2,2/)),BmV,1,1,(/2/),logT_is(jd),xlow,ylow,xup,yup)
							else
								BmV_is(jd)=(y-logy_is1(jd))*(BmV_is2(jd)-BmV_is1(jd))/(logy_is2(jd)-logy_is1(jd)) &
											& +BmV_is1(jd)
								logT_is(jd)=(y-logy_is1(jd))*(logT_is2(jd)-logT_is1(jd))/(logy_is2(jd)-logy_is1(jd)) &
											& +logT_is1(jd)
							end if 
						
							deallocate(BmV_isor); deallocate(logy_isor);
							deallocate(logT_isor); deallocate(Diff_logy)
						end do 
						
						dCL_is=abs(BmV_is(1)-BmV_is(2))		!dist between isoc in the BmV-logg Diagram
						dHR_is=abs(logT_is(1)-logT_is(2))	!dist between isoc in the logTeff-logg Diagram
						if ((dCL_is>dlim .and. dHR_is>dlim)) then !condition added 17/2/14
			!				dCL=abs(BmV-BmV_is)				!vector of distances between the star and the two isoc
							if (BmV<BmV_is(1) .and. BmV<BmV_is(2)) then 
								state=1
							else if (BmV>BmV_is(1) .and. BmV>BmV_is(2)) then 
								state=2
							else
								state=3
							end if
							!!
							call calibrateDiag(state,dClim,BmV_is,logT_is,BmV,logTeff)
							if (idCol.eq.1) then
								xtau=BmV
							else
								xtau=logTeff
							end if
							!! 
							if (isnan(logTeff)) then 
								cycle
							end if 
						else if (age_s==size(ix_sorted,1)) then 
							cycleW=1
							cycle
						end if 
					end do
					if (cycleW.eq.1) then
						cycle
					end if
					Teff=10.**logTeff
					if ((.not.calibSPEC .and. Teffinput.ne.-1)) then !Inside calibNoD => .not.calibSPEC is always true
						if (abs(Teff-Teffinput)>300) then !Inconsistency between BmV and spectrscopic Teff
							cycle
						end if 
					end if 
					I_Teff=0.01*Teff
					I_logTeff=0.01*log10(e)
					!!Compatibility
					if (loggAvail0.and.rhoAvail0) then 
						!!rho, logg both available. I determine input values for M, R 
						loggf=loggf0
						rhof=rhof0
						R2=g0/10.**loggSun*rhoSun/rho0
						I_R2=R2*(I_g0/g0+I_rho0/rho0)
						if (isEq(Rf,-1.D0,2)) then
							R=R2
							I_R=I_R2
						else
							call wMean(R1,I_R1,R2,I_R2,R,I_R)
						end if
						M=(g0/10.**loggSun)**3*(rhoSun/rho0)**2
						I_M=M*(3.*I_g0/g0+2.*I_rho0/rho0)
						L=R**2*(Teff/TeffSun)**4
						I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
						logL=log10(L)
						I_logL=I_L/L*log10(e)
						!!!!!!Check compatibility between rho and logg 16/10/17 
						!!Modified to account for the mass uncertainty 26/10/2018 
						if (.not.(M+I_M.lt.minval(MTrAv) .or. M-I_M.gt.maxval(MTrAv))) then !open tracks only if inside MassRange 
							call selectMfromTracks(MTrAv,M,I_M,Mvec)
							allocate(cumInt(size(Mvec)))
							jk=0
							do jj=1,size(Mvec)
								indxM=minloc(abs(Mvec(jj)-MTrAv),1)
								xMTi=Tndxi(indxZt(indxZ),indxM)
								xMTf=Tndxf(indxZt(indxZ),indxM)
								allocate(TrRhoG(xMTf-xMTi+1,size(TrackTab,2)))
								TrRhoG=TrackTab(xMTi:xMTf,:)
			
								allocate(MTr(size(TrRhoG,1))); allocate(RTr(size(TrRhoG,1)))
								allocate(logRhoTr(size(TrRhoG,1))); allocate(loggTr(size(TrRhoG,1)))
								
								MTr=TrRhoG(:,cM_T) !Msun
								RTr=10.**TrRhoG(:,clogR_T) !cm
								logRhoTr=log10(MTr/(RTr/RSun)**3*rhoSun)
								loggTr=log10(MTr/(RTr/RSun)**2)+loggSun
								call findYgivenX_v(TrRhoG,logg,loggTr,clogTe_T,logTeTrtmp)
								if (allocated(logTeTrtmp)) then
									jk=jk+1
									if (jk.eq.1) then
										allocate(logTeTr(size(logTeTrtmp)))
										logTeTr=logTeTrtmp
									else
										call append1D(logTeTr,logTeTrtmp)
									end if
									cumInt(jj)=size(logTeTr)
									deallocate(logTeTrtmp)
								else
									cumInt(jj)=0
								end if
								
								deallocate(TrRhoG)
								deallocate(MTr); deallocate(RTr)
								deallocate(logRhoTr); deallocate(loggTr)
							end do
							logTeJ=logTeff !I don't use Johnson (1966) relation because I've just calibrated Teff
							DTeJ=I_logTeff !directly set to the true uncertainty
							if (allocated(logTeTr)) then 
								allocate(TeB(size(logTeTr)))
								TeB=((logTeTr>logTeJ-DTeJ).and.(logTeTr<logTeJ+DTeJ))
								if (.not.any(TeB)) then !all elements of TeB are 0 => logg and rho inconsistent
									if (I_logg0>I_logrho0) then !logrho is best determined
										if (isEq(Rf,-1.D0,2)) then
											call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
											call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
											loggf=-1.
											call setNaN(logg); call setNaN(I_logg)
											call setNaN(g); call setNaN(I_g)
											loggAvail=.false.
											loggAvailAS=.false.
										else !input R is available
											!Only using rho: determine M; re-determine logg (discard the input value) 
											loggf=-1.
											loggAvail=.false.
											loggAvailAS=.false.
											R=R1 !use just the input reliable value
											I_R=I_R1
											call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
											loggAvailCal=.true. !!logg now available thanks to calibration
										end if
									else !logg is best determined
										if (isEq(Rf,-1.D0,2)) then
											call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
											call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
											rhof=-1.
											call setNaN(rho); call setNaN(I_rho)
											call setNaN(logrho); call setNaN(I_logrho)
											rhoAvail=.false.
											rhoAvailAS=.false.
										else !input R is available
											!Only using logg: determine M; re-determine rho (discard the input value) 
											rhof=-1.
											rhoAvail=.false.
											rhoAvailAS=.false.
											R=R1 !use just the input reliable value
											I_R=I_R1
											call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
											rhoAvailCal=.true.
										end if
									end if 
								else
									!!!!!!!!!
									call consistentM(TeB,cumInt,Mvec,M,agreeM)
									if (.not.agreeM) then !new M has been recomputed inside the subroutine in case agreeM is false
														  !This M will be used in case R is not available directly from input
										if (I_logg0.gt.I_logrho0) then !logrho is best determined
											if (isEq(Rf,-1.D0,2)) then
												loggAvail=.false.
												loggAvailAS=.false.
												!g: mantain original input uncertainty
												call computeLRgfromMrho(M,I_M,rho,I_rho,Teff,I_Teff,L,logL,I_L,I_logL,R,I_R,g,logg)
												loggAvailCal=.true.
											else !input R is available
												!Only using rho: determine M; re-determine logg (discard the input value) 
												loggf=-1.
												loggAvail=.false.
												loggAvailAS=.false.
												R=R1
												I_R=I_R1
												call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
												loggAvailCal=.true. !!logg now available thanks to calibration
											end if
										else
											if (isEq(Rf,-1.D0,2)) then
												rhoAvail=.false.
												rhoAvailAS=.false.
												!rho in g/cm3. Mantain original input uncertainty on rho
												call computeLRrhofromMg(M,I_M,g,I_g,Teff,I_Teff,L,logL,I_L,I_logL,R,I_R,rho,logrho)
												rhoAvailCal=.true.
											else !input R is available
												!Only using logg: determine M; re-determine rho (discard the input value) 
												rhof=-1.
												rhoAvail=.false.
												rhoAvailAS=.false.
												R=R1
												I_R=I_R1
												call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
												rhoAvailCal=.true.
											end if
										end if
										if (isEq(Rf,-1.D0,2)) then
											R=g/10**loggSun*rhoSun/rho
											I_R=R*(I_g/g+I_rho/rho)
											L=R**2*(Teff/TeffSun)**4
											I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
											logL=log10(L)
											I_logL=I_L/L*log10(e)
										end if
									else
										if (I_logg0.lt.I_logrho0) then
											bestLogg=.true.
										else
											bestLogg=.false.
										end if
										if (.not.isEq(Rf,-1.D0,2)) then
										!input R is available: weighted mean to infer M
											Mlogg=g0/10.**loggSun*R**2
											I_Mlogg=Mlogg*(I_g0/g0+2.*I_R/R)
											Mrho=rho0/rhoSun*R**3 !Mo
											I_Mrho=Mrho*(I_rho0/rho0+3.*I_R/R)
											call wMean(Mlogg,I_Mlogg,Mrho,I_Mrho,M,I_M)
											!L already defined since the beginning from R
										end if
									end if
									!!!!!!!!!!!!!
								end if
								deallocate(logTeTr); deallocate(TeB)
							else !Track has been opened, but no logTe compatible with logg were found
								if (I_logg0>I_logrho0) then !logrho is best determined
									if (isEq(Rf,-1.D0,2)) then
										call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
										call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
										loggf=-1.
										call setNaN(logg); call setNaN(I_logg)
										call setNaN(g); call setNaN(I_g)
										loggAvail=.false.
										loggAvailAS=.false.
									else !input R is available
										!Only using rho: determine M; re-determine logg (discard the input value) 
										loggf=-1.
										loggAvail=.false.
										loggAvailAS=.false.
										R=R1 !use just the input reliable value
										I_R=I_R1
										call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
										loggAvailCal=.true. !!logg now available thanks to calibration
									end if
								else !logg is best determined
									if (isEq(Rf,-1.D0,2)) then
										call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
										call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
										rhof=-1.
										call setNaN(rho); call setNaN(I_rho)
										call setNaN(logrho); call setNaN(I_logrho)
										rhoAvail=.false.
										rhoAvailAS=.false.
									else !input R is available
										!Only using logg: determine M; re-determine rho (discard the input value) 
										rhof=-1.
										rhoAvail=.false.
										rhoAvailAS=.false.
										R=R1 !use just the input reliable value
										I_R=I_R1
										call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
										rhoAvailCal=.true.
									end if
								end if 
							end if
							deallocate(cumInt) 
						else
							if (I_logg0>I_logrho0) then !logrho is best determined
								if (isEq(Rf,-1.D0,2)) then
									call setNan(R); call setNan(I_R); call setNan(M); call setNan(I_M)
									call setNan(L); call setNan(I_L); call setNan(logL); call setNan(I_logL)
									loggf=-1.
									call setNan(logg); call setNan(I_logg)
									call setNan(g); call setNan(I_g)
									loggAvail=.false.
									loggAvailAS=.false.
								else !input R is available
									!Only using rho: determine M; re-determine logg (discard the input value) 
									loggf=-1.
									loggAvail=.false.
									loggAvailAS=.false.
									R=R1 !use just the input reliable value
									I_R=I_R1
									call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
									loggAvailCal=.true. !!logg now available thanks to calibration
								end if
							else !logg is best determined
								if (isEq(Rf,-1.D0,2)) then
									call setNan(R); call setNan(I_R); call setNan(M); call setNan(I_M)
									call setNan(L); call setNan(I_L); call setNan(logL); call setNan(I_logL)
									rhof=-1.
									call setNan(rho); call setNan(I_rho)
									call setNan(logrho); call setNan(I_logrho)
									rhoAvail=.false.
									rhoAvailAS=.false.
								else !input R is available
									!Only using logg: determine M; re-determine rho (discard the input value) 
									rhof=-1.
									rhoAvail=.false.
									rhoAvailAS=.false.
									R=R1 !use just the input reliable value
									I_R=I_R1
									call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
									rhoAvailCal=.true.
								end if
							end if 
						end if 
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
					else
						if (.not.isEq(Rf,-1.D0,2)) then
							R=R1
							I_R=I_R1
							if (loggAvail.or.loggAvailAS) then
								call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
								rhoAvailCal=.true.
							else if (rhoAvail.or.rhoAvailAS) then
								call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
								loggAvailCal=.true. !!logg available thanks to calibration
							else !noGproxy
								L=R**2*(Teff/TeffSun)**4
								I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
								logL=log10(L)
								I_logL=I_L/L*log10(e)
							end if
						end if
					end if
					!!End compatibility
				end do
				if (cycleW.eq.1) then
					cycle
				end if
			end if
			deallocate(logy_iso)
		end if
		deallocate(ix_sorted)
		
		!!!Check Activity/vsini. If inactive stars or slow rotators, consider only isochrones with logt>=logtCutoff 
		check=0
		logtlim=-1
		if (.not.isEq(logRHK,0.D0,2)) then 
			logtMmod=c0-17.912*(logRHK+DlR)-1.6675*(logRHK+DlR)**2 !modified Mamayek
			if (logtMmod<logtCutoff5) then 
				logtlim1=logtMmod
			else
				logtlim1=logtCutoff5
			end if 
		else
			logtlim1=-1
		end if 
		
		!
		Mavail=(calibHRD.and.(loggAvail.or.loggAvailAS.or.rhoAvail.or.rhoAvailAS)).or. &
				& ((loggAvail.or.loggAvailAS.or.loggAvailCal).and.(rhoAvail.or.rhoAvailAS.or.rhoAvailCal)) &
				& .or.(calibSPEC.and.gProxyAvail.and..not.isEq(Rf,-1.D0,2)) & 
				& .or.(calibNoD.and.gProxyAvail.and..not.isEq(Rf,-1.D0,2))
		!
		!Searching T.O. of the oldest isochrone!!!!!!!!!!!!!!!!
		allocate(y_iso(size(Isoc,1)))
		caliblogL=calibHRD.or.((loggAvail.or.loggAvailAS.or.loggAvailCal).and.(rhoAvail.or.rhoAvailAS.or.rhoAvailCal)) &
					& .or.((calibSPEC.or.calibNoD).and..not.isEq(Rf,-1.D0,2))
		if (caliblogL) then !19/9/2017 !logL avail
			logY=logL
			y_iso=logL_iso
		else if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
				& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))) then 
			logY=-logg !to let following inequality logY< refer to lowMS stars
			y_iso=logg_iso
		else
			logY=-logrho !to let following inequality logY< refer to lowMS stars
			y_iso=logrho_iso
		end if 
		call searchTOold(y_iso,logTeff_iso,caliblogL,logt_iso,last_logt,logY_soglia)
		
		logY2=logY
		logY2_soglia=logY_soglia !consider these values later for low MS stars

		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

		if (.not.isEq(vsini,-1.D0,2) .or. .not.isEq(P,-1.D0,2)) then 
			vcheck=.true.
			!Rv: radius to be used in vsini relation
			if (caliblogL) then 
				Rv=R
			else	!infer a rough estimate for stellar radius to be used for gyrochronology
				if (MminT<MlowMS) then !tracks of lower mass stars are available
					!define logY_sogliaVsini: open M=0.5Mo track; Z=Zstar and search logL/logg(last_logt)
					! that will be logY_sogliaVsini
					indxM=minloc(abs(MlowMS-MTrAv),1)
					xMTi=Tndxi(indxZt(indxZ),indxM)
					xMTf=Tndxf(indxZt(indxZ),indxM)
					allocate(Track05(xMTf-xMTi+1,size(TrackTab,2)))
					Track05=TrackTab(xMTi:xMTf,:)
					MilowMS=MTrAv(indxM)
						
					!Establish maximim luminosity for a ~0.5solarMass star 
					if (caliblogL) then !Never enter this condition. Could only in the following two 
						call InterpLin_M( reshape((/Track05(:,ct_T),Track05(:,clogL_T)/),(/size(Track05,1),2/)), &
										& 10.**last_logt,1,1,(/ 2/),logY_sogliaVsini1,xlow,ylow,xup,yup)
						logY_sogliaVsini=logY_sogliaVsini1(1)
					else if (loggAvail.or.loggAvailAS.or.loggAvailCal) then 
						call InterpLin_M( reshape((/Track05(:,ct_T),Track05(:,clogR_T)/),(/size(Track05,1),2/)), &
										& 10.**last_logt,1,1,(/ 2/),logR_Tr,xlow,ylow,xup,yup)!cm
						R_Tr=10.**logR_Tr(1)
						logY_sogliaVsini=-(log10(MilowMS/(R_Tr/RSun)**2)+loggSun) !minus to let hold inequality
																				  ! logY>logY_sogliaVsini
					else
						call InterpLin_M( reshape((/Track05(:,ct_T),Track05(:,clogR_T)/),(/size(Track05,1),2/)), &
										& 10.**last_logt,1,1,(/ 2/),logR_Tr,xlow,ylow,xup,yup)!cm
						R_Tr=10.**logR_Tr(1)
						rho_Tr=MilowMS/(R_Tr/RSun)**3*rhoSun
						logY_sogliaVsini=-log10(rho_Tr)
					end if
					deallocate(Track05)
									
					if (logY<logY_sogliaVsini) then !very low MS star. R=R(t) almost constant
						!Select the middle-grid isochrone containing only the lines where M<MilowMS
						!and choose R in correspondence of logY 
						mid_n=nint((logt_halfMS-first_logt)/logtstep+1)
						
						call M2R(logt_iso,mid_n,Isoc,MilowMS,caliblogL,loggAvail.or.loggAvailAS,loggAvailCal,logY,Rv)
					else
						XT=logTeff-4.1
						if (rhoAvail.or.rhoAvailAS.or.rhoAvailCal) then !implementation my relation R=R(logTeff,logrho,[Fe/H])
							vecT=(/ 1.D0,XT,XT**2,XT**3,logrho,logrho**2,logrho**3,FeH_,FeH_**2,FeH_**3 /)
							Rv=10.**(dot_product(brho,vecT))
						else if (loggAvail.or.loggAvailAS.or.loggAvailCal) then !implementation my relation R=R(logTeff,logg,[Fe/H])
							vecT=(/ 1.D0,XT,XT**2,XT**3,logg,logg**2,logg**3,FeH_,FeH_**2,FeH_**3 /)
							Rv=10.**(dot_product(blogg,vecT))
						else !Should never enter
							vcheck=.false.
						end if 
					end if
				else
					XT=logTeff-4.1
					if (rhoAvail.or.rhoAvailAS.or.rhoAvailCal) then !implementation of my relation R=R(logTeff,logrho,[Fe/H])
						vecT=(/ 1.D0,XT,XT**2,XT**3,logrho,logrho**2,logrho**3,FeH_,FeH_**2,FeH_**3 /)
						Rv=10.**(dot_product(brho,vecT))
					else if (loggAvail.or.loggAvailAS.or.loggAvailCal) then !implementation of my relation R=R(logTeff,logg,[Fe/H])
						vecT=(/ 1.D0,XT,XT**2,XT**3,logg,logg**2,logg**3,FeH_,FeH_**2,FeH_**3 /)
						Rv=10.**(dot_product(blogg,vecT))
					else !Should never enter
						vcheck=.false.
					end if 
				end if 
			end if 
			if (vcheck.eqv..true.) then !if vsini or P available should be always true
				if (isEq(vsini,0.D0,2)) then !slow rotating stars for which vsini is reported as zero
					logtlim2=logtCutoff25 !No +0.05. 
										  !It's already been taken into account in the definition of logtCutoff25
				else
					if (.not.isEq(vsini,-1.D0,2).and.(isEq(P,-1.D0,2))) then 
						Omega=4./pi*vsini/(Rv*RSunKm)
						!	logt_Deniss=log10(((OmSun/Omega)^2-1)*(1+A)/(2*(A/2+B))*(tSunLG-tZAMS)+tSunLG); 
						!	logtlim=logt_Deniss; 
						P=2.*pi/Omega/86400 !days
					end if
					call InterpLin_M(GyroTab,xtau,cBrif,1,(/ ctau /),tau1,xlow,taulow,xup,tauup)
					tau=tau1(1)
					if (.not.isEq(tau,0.D0,2) .and. x>=xlow .and. x<=xup) then 
						logt_Barnes=log10(tau/kC*log(P/P0)+kI/(2.*tau)*(P**2-P0**2))-3 ![Gyr]
						if (10.**logt_Barnes>DtBarnes) then 
							logt_Barnes=log10(10.**logt_Barnes-DtBarnes)+9 ![yrs]
						else
							logt_Barnes=first_logt
						end if 
						if (logt_Barnes<logtCutoff25) then 
							logtlim2=logt_Barnes
						else
							logtlim2=logtCutoff25
						end if 
					else
						logtlim2=-1 !tau=0 => Barnes relation not defined. I don't set any gyro age
					end if
				end if 
			else
				logtlim2=-1
			end if 
		else
			logtlim2=-1
		end if 
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
		if ((rhoAvail .or. rhoAvailCal)) then !rhof.ne.-1
			if (logY<logY_soglia) then
				if (.not.caliblogL) then
					!re-set correct values of logg (or logrho) and of logg_soglia (or logrho_soglia)
					logY=-logY
					logY_soglia=-logY_soglia
				end if 
				call trovaIndici(logt_iso,ndxrhoi,ndxrhof)
				ndxTot=size(ndxrhoi) !=size(ndxrhof)
				jrho=1
				rhot=-100 !initialize like this so that it enters the cycle at least once
				do while (logt_iso(ndxrhoi(jrho))<8 .and. rho-I_rho>rhot .and. jrho<=ndxTot) 
					!Consider only isochrones with t<100Myr. jrho<=ndxTot will always hold
					allocate(logY_isot(ndxrhof(jrho)-ndxrhoi(jrho)+1))
					allocate(rho_isot(ndxrhof(jrho)-ndxrhoi(jrho)+1))
					logY_isot=y_iso(ndxrhoi(jrho):ndxrhof(jrho))
					rho_isot=rho_iso(ndxrhoi(jrho):ndxrhof(jrho))
					ndxrho=minloc(abs(logY-logY_isot),1) 
					!logYt (useless) is the value reported by the isochrone t, that is closer to stellar logY
					!In correspondence of this value, I'll select the corresponding isochronal stellar density
	!					logYt=logY_isot(ndxrho)
					!rhotmp(jrho)=rho_isot(ndxrho) -->useless to consider it as a vector
					
					rhot=rho_isot(ndxrho)
					jrho=jrho+1
					deallocate(logY_isot); deallocate(rho_isot)
				end do
												
				if (jrho==2) then !While cycle done only once. Don't discard any isochrone
					logtlim3=logt_iso(1) !first age value reported by the isochrones=>don't discard any isochrone!!
				else
					logtlim3=logt_iso(ndxrhoi(jrho-1))+logtstep !added logtstep so that all isochrones with
																!logt<=logt_iso(ndxrhoi(ii)) are discarded. In fact
																!isochrones mantained start from logtlim included
					if (logtlim3>logtCutoff5) then 
						logtlim3=logtCutoff5
					end if 
				end if
				deallocate(ndxrhoi); deallocate(ndxrhof)
			else
				logtlim3=-1
			end if
		else
			logtlim3=-1
		end if
		deallocate(y_iso)
		if (.not.isEq(YMg,-100.D0,2).and.FeH_>=-0.2 .and. FeH_<=0.2) then !metall condition
			t_Nissen=(YMg-aYMg)/bYMg*1.e9 ![yrs]
			if (t_Nissen>log10(DtNissen)) then
				logt_Nissen=log10(t_Nissen-DtNissen) !at least this age
			else
				logt_Nissen=first_logt
			end if
			if (logt_Nissen<logtCutoff25) then
				logtlim4=logt_Nissen
			else
				logtlim4=logtCutoff25
			end if
		else
			logtlim4=-1
		end if
						 
		limAges=(/ logtlim1,logtlim2,logtlim3,logtlim4 /)
		klim0=0
		do klim=1,size(limAges) 
			if (.not.isEq(limAges(klim),-1.D0,2)) then 
				klim0=1
				exit
			end if 
		end do
		
		if (klim0==1) then 
			logtlim=maxval(limAges)
			check=maxloc(limAges,1)
		else
			logtlim=-1
		end if
		
		if (.not.isEq(logtlim,-1.D0,2)) then !activity check done
			logtCutoff0=idnint(logtlim*100)
			resto=mod(logtCutoff0,idnint(logtstep*100))
			if (resto>floor(anint(logtstep*100)/2)) then 
				logtCutoff=(logtCutoff0+(dnint(logtstep*100)-resto))/100.
			else
				logtCutoff=(logtCutoff0-resto)/100.
			end if 
			if (photIsocAvail) then !anint should be useless
				c_i=13
				allocate(Iso_i(size(t_i),c_i))
				Iso_i=reshape((/ dnint(log10(t_i)*100)/100,t_i,logTeff_i,Teff_i,logL_i,L_i,logg_i,g_i,M_i, &
							& logrho_i,rho_i,BmV_i,BC_i /),(/size(Iso_i,1),c_i/))
				!added rho_i 19/9/2017
			else
				c_i=11
				allocate(Iso_i(size(t_i),c_i))
				Iso_i=reshape((/ dnint(log10(t_i)*100)/100,t_i,logTeff_i,Teff_i,logL_i,L_i,logg_i,g_i,M_i, &
							& logrho_i,rho_i /),(/size(Iso_i,1),c_i/))
				!added rho_i 19/9/2017
			end if
			allocate(Iso_indx(size(Iso_i,1))) 
			call indexx(Iso_i(:,1),Iso_indx)
			!Check if Iso_i it's overwritten correctly
			Iso_i=Iso_i(Iso_indx,:)
			deallocate(Iso_indx)
			do iIso=1,size(Iso_i,1) 
				if (isEq(Iso_i(iIso,1),logtCutoff,2)) then 
					iIso8=iIso
					exit
				end if 
			end do
			allocate(Iso_i2(size(Iso_i,1)-iIso8+1,c_i))
			Iso_i2=Iso_i(iIso8:size(Iso_i,1),:)
			deallocate(Iso_i)
			!dealloc _i and then re-alloc with proper dimension (lower than before in case of act check)
			deallocate(t_i); deallocate(logTeff_i); deallocate(Teff_i); deallocate(logL_i)
			deallocate(L_i); deallocate(logg_i); deallocate(g_i); deallocate(M_i)
			deallocate(logrho_i); deallocate(rho_i)
			allocate(t_i(size(Iso_i2,1))); allocate(logTeff_i(size(Iso_i2,1))); allocate(Teff_i(size(Iso_i2,1)))
			allocate(logL_i(size(Iso_i2,1))); allocate(L_i(size(Iso_i2,1))); allocate(logg_i(size(Iso_i2,1)))
			allocate(g_i(size(Iso_i2,1))); allocate(M_i(size(Iso_i2,1))); allocate(logrho_i(size(Iso_i2,1)))
			allocate(rho_i(size(Iso_i2,1)))
			t_i=Iso_i2(:,2)
			logTeff_i=Iso_i2(:,3)
			Teff_i=Iso_i2(:,4)
			logL_i=Iso_i2(:,5)
			L_i=Iso_i2(:,6)
			logg_i=Iso_i2(:,7)
			g_i=Iso_i2(:,8)
			M_i=Iso_i2(:,9)
			logrho_i=Iso_i2(:,10) !19/9/2017
			rho_i=Iso_i2(:,11) !19/9/2017
			if (photIsocAvail) then 
				deallocate(BmV_i); deallocate(BC_i)
				allocate(BmV_i(size(Iso_i2,1))); allocate(BC_i(size(Iso_i2,1)))
				BmV_i=Iso_i2(:,12) !19/9/2017
				BC_i=Iso_i2(:,13) !19/9/2017
			end if 
			logt8=logtCutoff
			deallocate(Iso_i2)
		else
			logt8=0
		end if
		!!!!End check activity/vsini
				
		if (IncZ==1) then !load further metallic grids, accounting for I_FeH.
						  !Error bar extension=1 sigma and are "symmetric" in Z plane
			!both flag in multipleZisocPhSCP input parameters are set to 1
			chooseLogg=((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
						& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))
			if (I_FeH_>0) then !avoid that flags are considered as actual uncertainties
				if (calibHRD) then 
					xx=BmV
					yy=Vass
				end if 
				if (calibNoD) then 
					xx=BmV
					if (caliblogL) then
						yy=logL
					else if (chooseLogg) then 
						yy=logg
					else
						yy=logrho
					end if 
				end if 
				if (calibSPEC) then 
					xx=logTeff
					if (caliblogL) then
						yy=logL
					else if (chooseLogg) then 
						yy=logg
					else
						yy=logrho
					end if 
				end if
				call multipleZisocPhSCP( FeH_, I_FeH_, Zvec, model, xx, yy, &
					& logt8, percorso, 1.D0, 1, calibHRD, calibNoD, calibSPEC, chooseLogg, &
					& caliblogL,hstar, hstarlim, Isoc_i, Zlow, Zup,useColor,idCol)
				if (allocated(Isoc_i)) then !isochrone grids (WITH the reference one INCLUDED) are loaded
					deallocate(t_i); deallocate(logTeff_i); deallocate(Teff_i)
					deallocate(logL_i); deallocate(L_i); deallocate(logg_i)
					deallocate(g_i); deallocate(M_i); deallocate(logrho_i)
					deallocate(rho_i)
					allocate(t_i(size(Isoc_i,1))); allocate(logTeff_i(size(Isoc_i,1)))
					allocate(Teff_i(size(Isoc_i,1))); allocate(logL_i(size(Isoc_i,1)))
					allocate(L_i(size(Isoc_i,1))); allocate(logg_i(size(Isoc_i,1)))
					allocate(g_i(size(Isoc_i,1))); allocate(M_i(size(Isoc_i,1)))
					allocate(logrho_i(size(Isoc_i,1))); allocate(rho_i(size(Isoc_i,1)))
					t_i=Isoc_i(:,2)
					logTeff_i=Isoc_i(:,3)
					Teff_i=Isoc_i(:,4)
					logL_i=Isoc_i(:,5)
					L_i=Isoc_i(:,6)
					logg_i=Isoc_i(:,7)
					g_i=Isoc_i(:,8)
					M_i=Isoc_i(:,9)
					logrho_i=Isoc_i(:,10) !19/9/2017
					rho_i=Isoc_i(:,11) !19/9/2017
					if (photIsocAvail) then
						deallocate(BmV_i); deallocate(BC_i)
						allocate(BmV_i(size(Isoc_i,1))); allocate(BC_i(size(Isoc_i,1)))
						BmV_i=Isoc_i(:,12) !19/9/2017
						BC_i=Isoc_i(:,13) !19/9/2017
					end if
					deallocate(Isoc_i) 
				end if 
			end if 
		end if 
		
		!!!only theoretical data _i
		allocate(Gauss(size(logTeff_i)))
		allocate(w(size(logTeff_i)))
		if (gProxyAvail) then 
			if (calibHRD .or. ((loggAvail.or.loggAvailAS.or.loggAvailCal) .and. &
				& (rhoAvail.or.rhoAvailAS.or.rhoAvailCal)).or.((calibSPEC.or.calibNoD).and..not.isEq(Rf,-1.D0,2))) then
				!I can recover Teff, L, M, logg
				if (logY2<logY2_soglia) then 
					!!!Find ZAMS in the Track: logL starts increasing after the decreasing along the Hayashi line 
					call setThreshold(caliblogL,(((loggAvail.or.loggAvailAS).and.(rhoAvail.or.rhoAvailAS).and.bestLogg) &
									& .or.((loggAvail.or.loggAvailAS).and..not.(rhoAvail.or.rhoAvailAS))),logY2, &
									& logY2_soglia,logY0,logY0_soglia,cZAMS,cyT)
					
					allocate(yT_iso(size(Isoc,1)))
					if (cyT.eq.0) then !happens for logrho that hasn't a direct column in the isoch grid
						yT_iso=logrho_iso
					else
						yT_iso=Isoc(:,cyT)
					end if
					
					call ifComputeIsocAge(ZAMSz,cZAMS,logY0,logt_iso,logTeff_iso,yT_iso,I_logTeff,ageIsoc)
					
					deallocate(yT_iso)
				else
					ageIsoc=1
				end if 
				indxM=minloc(abs(M-MTrAv),1)
				xMTi=Tndxi(indxZt(indxZ),indxM)
				xMTf=Tndxf(indxZt(indxZ),indxM)
				allocate(Tracks(xMTf-xMTi+1,size(TrackTab,2)))
				Tracks=TrackTab(xMTi:xMTf,:)
								
				allocate(t_Tracks(size(Tracks,1)));allocate(logL_Tracks(size(Tracks,1)))
				allocate(logTeff_Tracks(size(Tracks,1)));allocate(vEvo(size(t_i)))
				t_Tracks=Tracks(:,ct_T)
				logL_Tracks=Tracks(:,clogL_T)
				logTeff_Tracks=Tracks(:,clogTe_T)
				deallocate(Tracks)
				
				call computeVevo(t_Tracks,logL_Tracks,logTeff_Tracks,t_i,varTrlim,vEvo)
				
				deallocate(t_Tracks);deallocate(logL_Tracks);deallocate(logTeff_Tracks)
				
				cRif=cVelL
				
				call computeVrif(TabVel,cRif,M,vRif)
								
				Gauss=1./(2.*pi*I_logTeff*I_logL)*e**(-0.5*((logTeff-logTeff_i)/I_logTeff)**2)* &
						& e**(-0.5*((logL-logL_i)/I_logL)**2)
				w=1./(((logL-logL_i)/I_logL)**2+((logTeff-logTeff_i)/I_logTeff)**2+((M-M_i)/I_M)**2+ &
						& ((logg-logg_i)/I_logg)**2+(log10(vRif/vEvo))**2)
				deallocate(vEvo)
			else
				if (logY2<logY2_soglia) then 
					call setThreshold(caliblogL,loggAvail.or.loggAvailAS,logY2,logY2_soglia,logY0,logY0_soglia,cZAMS,cyT)
					allocate(yT_iso(size(Isoc,1)))
					if (cyT.eq.0) then !happens for logrho that hasn't a direct column in the isoch grid
						yT_iso=logrho_iso
					else
						yT_iso=Isoc(:,cyT)
					end if
					
					call ifComputeIsocAge(ZAMSz,cZAMS,logY0,logt_iso,logTeff_iso,yT_iso,I_logTeff,ageIsoc)
															
					deallocate(yT_iso)
				else
					ageIsoc=1
				end if 
				
				allocate(y_i(size(logg_i))) !equal to the size of logrho_i
				if (loggAvail.or.loggAvailAS) then !only logg available
					y=logg
					I_y=I_logg
					y_i=logg_i
					cRif=cVelg
				else !only logrho available
					y=logrho
					I_y=I_logrho
					y_i=logrho_i
					cRif=cVelrho
				end if 
				!!!!!!!!!!!!!!!! 
				!!!!!!!!!!!!!!!! 
				Gauss=1./(2.*pi*I_logTeff*I_y)*e**(-0.5*((logTeff-logTeff_i)/I_logTeff)**2)*e**(-0.5*((y-y_i)/I_y)**2)
				w=1./(((logTeff-logTeff_i)/I_logTeff)**2+((y-y_i)/I_y)**2)
				Spesi=sum(w*Gauss)
				if (Spesi<=1.D-200) then 
					cycle
				end if
				
				M_starTr=dot_product(w,M_i*Gauss)/Spesi
				
				indxM=minloc(abs(M_starTr-MTrAv),1)
				xMTi=Tndxi(indxZt(indxZ),indxM)
				xMTf=Tndxf(indxZt(indxZ),indxM)
				allocate(Tracks(xMTf-xMTi+1,size(TrackTab,2)))
				Tracks=TrackTab(xMTi:xMTf,:)
!				M_in=MTrAv(indxM)
	
				allocate(t_Tracks(size(Tracks,1))); allocate(logTeff_Tracks(size(Tracks,1)))
				allocate(R_Tracks(size(Tracks,1))); allocate(logY_Tracks(size(Tracks,1)))
				allocate(vEvo(size(t_i)))
				t_Tracks=Tracks(:,ct_T)
				logTeff_Tracks=Tracks(:,clogTe_T)
				R_Tracks=10.**Tracks(:,clogR_T) !cm
				M_Tracks=Tracks(1,cM_T)
				deallocate(Tracks)
				if (loggAvail.or.loggAvailAS) then !only logg available
					logY_Tracks=log10(M_Tracks/(R_Tracks/RSun)**2)+loggSun
				else !only logrho available
					logY_Tracks=log10(M_Tracks/(R_Tracks/RSun)**3)+log10(rhoSun)
				end if
				call computeVevo(t_Tracks,logY_Tracks,logTeff_Tracks,t_i,varTrlim,vEvo)
				deallocate(t_Tracks); deallocate(logTeff_Tracks); deallocate(R_Tracks)
				deallocate(logY_Tracks)
				
				call computeVrif(TabVel,cRif,M_starTr,vRif)
				
				w=1./(((logTeff-logTeff_i)/I_logTeff)**2+((y-y_i)/I_y)**2+(log10(vRif/vEvo))**2)
				deallocate(vEvo)
				Spesi=sum(w*Gauss)
				if (Spesi<=1.D-200) then 
					cycle
				end if
				M_starTr2=dot_product(w,M_i*Gauss)/Spesi
				kTr=1
				do while (abs((M_starTr2-M_starTr)/M_starTr)>DMTr .and. kTr<10)
					if (logY2<logY2_soglia) then 
						call setThreshold(caliblogL,loggAvail.or.loggAvailAS,logY2,logY2_soglia,logY0,logY0_soglia,cZAMS,cyT)
						allocate(yT_iso(size(Isoc,1)))
						if (cyT.eq.0) then !happens for logrho that hasn't a direct column in the isoch grid
							yT_iso=logrho_iso
						else
							yT_iso=Isoc(:,cyT)
						end if
						
						call ifComputeIsocAge(ZAMSz,cZAMS,logY0,logt_iso,logTeff_iso,yT_iso,I_logTeff,ageIsoc)
										
						deallocate(yT_iso)
					else
						ageIsoc=1
					end if 
										
					M_starTr=M_starTr2
					indxM=minloc(abs(M_starTr-MTrAv),1)
					xMTi=Tndxi(indxZt(indxZ),indxM)
					xMTf=Tndxf(indxZt(indxZ),indxM)
					allocate(Tracks(xMTf-xMTi+1,size(TrackTab,2)))
					Tracks=TrackTab(xMTi:xMTf,:)
	
					allocate(t_Tracks(size(Tracks,1))); allocate(logTeff_Tracks(size(Tracks,1)))
					allocate(R_Tracks(size(Tracks,1))); allocate(logY_Tracks(size(Tracks,1)));
					allocate(vEvo(size(t_i)))
					t_Tracks=Tracks(:,ct_T)
					logTeff_Tracks=Tracks(:,clogTe_T)
					R_Tracks=10.**Tracks(:,clogR_T) !cm
					M_Tracks=Tracks(1,cM_T)
					deallocate(Tracks)
					if (loggAvail.or.loggAvailAS) then !only logg available
						logY_Tracks=log10(M_Tracks/(R_Tracks/RSun)**2)+loggSun
					else !only logrho available
						logY_Tracks=log10(M_Tracks/(R_Tracks/RSun)**3)+log10(rhoSun)
					end if
					call computeVevo(t_tracks,logY_Tracks,logTeff_Tracks,t_i,varTrlim,vEvo)
					deallocate(t_Tracks); deallocate(logTeff_Tracks); deallocate(R_Tracks)
					deallocate(logY_Tracks)
					
					call computeVrif(TabVel,cRif,M_starTr,vRif)
					
					w=1./(((logTeff-logTeff_i)/I_logTeff)**2+((y-y_i)/I_y)**2+(log10(vRif/vEvo))**2)
					deallocate(vEvo)
					Spesi=sum(w*Gauss)
					if (Spesi<=1.D-200) then 
						cycle
					end if 
					M_starTr2=dot_product(w,M_i*Gauss)/Spesi
					kTr=kTr+1
				end do
				!!!!!!!!!!!!!!!!!!! 
				!!!!!!!!!!!!!!!!!!! 
			end if 
		else
			if (logY2<logY2_soglia) then 
				call setThreshold(caliblogL,loggAvail.or.loggAvailAS,logY2,logY2_soglia,logY0,logY0_soglia,cZAMS,cyT)
				allocate(yT_iso(size(Isoc,1)))
				if (cyT.eq.0) then !happens for logrho that hasn't a direct column in the isoch grid
					yT_iso=logrho_iso
				else
					yT_iso=Isoc(:,cyT)
				end if
				 
				call ifComputeIsocAge(ZAMSz,cZAMS,logY0,logt_iso,logTeff_iso,yT_iso,I_logTeff,ageIsoc)
									
				deallocate(yT_iso)
			else
				ageIsoc=1
			end if 
			
			Gauss=1./(2.*pi*I_logTeff*I_logL)*e**(-0.5*((logTeff-logTeff_i)/I_logTeff)**2)* &
					& e**(-0.5*((logL-logL_i)/I_logL)**2)
			w=1./(((logL-logL_i)/I_logL)**2+((logTeff-logTeff_i)/I_logTeff)**2)
			Spesi=sum(w*Gauss)
			if (Spesi<=1.D-200) then 
				cycle
			end if
			M_starTr=dot_product(w,M_i*Gauss)/Spesi
			indxM=minloc(abs(M_starTr-MTrAv),1)
			xMTi=Tndxi(indxZt(indxZ),indxM)
			xMTf=Tndxf(indxZt(indxZ),indxM)
			allocate(Tracks(xMTf-xMTi+1,size(TrackTab,2)))
			Tracks=TrackTab(xMTi:xMTf,:)
			allocate(t_Tracks(size(Tracks,1)));allocate(logL_Tracks(size(Tracks,1)))
			allocate(logTeff_Tracks(size(Tracks,1)));allocate(vEvo(size(t_i)))
			
			t_Tracks=Tracks(:,ct_T)
			logL_Tracks=Tracks(:,clogL_T)
			logTeff_Tracks=Tracks(:,clogTe_T)
			deallocate(Tracks)
			
			call computeVevo(t_Tracks,logL_Tracks,logTeff_Tracks,t_i,varTrlim,vEvo)
			
			deallocate(t_Tracks);deallocate(logL_Tracks);deallocate(logTeff_Tracks) 
			
			cRif=cVelL
			call computeVrif(TabVel,cRif,M_starTr,vRif)
			
			w=1./(((logL-logL_i)/I_logL)**2+((logTeff-logTeff_i)/I_logTeff)**2+(log10(vRif/vEvo))**2)
			deallocate(vEvo)
			Spesi=sum(w*Gauss)
			if (Spesi<=1.D-200) then 
				cycle
			end if 
			M_starTr2=dot_product(w,M_i*Gauss)/Spesi
			kTr=1
			do while (abs((M_starTr2-M_starTr)/M_starTr)>DMTr .and. kTr<10) 
				if (logY2<logY2_soglia) then 
					call setThreshold(caliblogL,loggAvail.or.loggAvailAS,logY2,logY2_soglia,logY0,logY0_soglia,cZAMS,cyT)
					allocate(yT_iso(size(Isoc,1)))
					if (cyT.eq.0) then !happens for logrho that hasn't a direct column in the isoch grid
						yT_iso=logrho_iso
					else
						yT_iso=Isoc(:,cyT)
					end if
					 
					call ifComputeIsocAge(ZAMSz,cZAMS,logY0,logt_iso,logTeff_iso,yT_iso,I_logTeff,ageIsoc)
										
					deallocate(yT_iso)
				else
					ageIsoc=1
				end if 
				
				M_starTr=M_starTr2
				indxM=minloc(abs(M_starTr-MTrAv),1)
				xMTi=Tndxi(indxZt(indxZ),indxM)
				xMTf=Tndxf(indxZt(indxZ),indxM)
				allocate(Tracks(xMTf-xMTi+1,size(TrackTab,2)))
				Tracks=TrackTab(xMTi:xMTf,:)
				allocate(t_Tracks(size(Tracks,1)));allocate(logL_Tracks(size(Tracks,1)))
				allocate(logTeff_Tracks(size(Tracks,1)));allocate(vEvo(size(t_i)))
	
				t_Tracks=Tracks(:,ct_T)
				logL_Tracks=Tracks(:,clogL_T)
				logTeff_Tracks=Tracks(:,clogTe_T)
				deallocate(Tracks)
				call computeVevo(t_Tracks,logL_Tracks,logTeff_Tracks,t_i,varTrlim,vEvo)
				deallocate(t_Tracks);deallocate(logL_Tracks);deallocate(logTeff_Tracks) 
			
				cRif=cVelL
				call computeVrif(TabVel,cRif,M_starTr,vRif)
				
				w=1./(((logL-logL_i)/I_logL)**2+((logTeff-logTeff_i)/I_logTeff)**2+(log10(vRif/vEvo))**2)
				deallocate(vEvo)
				Spesi=sum(w*Gauss)
				if (Spesi<=1.D-200) then 
					cycle
				end if 
				M_starTr2=dot_product(w,M_i*Gauss)/Spesi
				kTr=kTr+1
			end do
		end if
		
		Spesi=sum(w*Gauss)
		if (Spesi<=1.D-200) then 
			cycle
		end if
		t_star=dot_product(w,t_i*Gauss)/Spesi !filter the selected isochrone through a Gaussian
		!and then weight it through w. It's like as attributing it the weight w*Gauss 
		Teff_star=dot_product(w,Teff_i*Gauss)/Spesi
		L_star=dot_product(w,L_i*Gauss)/Spesi
		M_star=dot_product(w,M_i*Gauss)/Spesi
		g_star=dot_product(w,g_i*Gauss)/Spesi
		rho_star=dot_product(w,rho_i*Gauss)/Spesi
		!!! 
		if (photIsocAvail) then 
			BmV_star=dot_product(w,BmV_i*Gauss)/Spesi
			BC_star=dot_product(w,BC_i*Gauss)/Spesi
		end if 
		!!! 
		I_t_star=sqrt(dot_product(w*Gauss,(t_i-t_star)**2)/Spesi)
		I_Teff_star=sqrt(dot_product(w*Gauss,(Teff_i-Teff_star)**2)/Spesi)
		I_L_star=sqrt(dot_product(w*Gauss,(L_i-L_star)**2)/Spesi)
		I_M_star=sqrt(dot_product(w*Gauss,(M_i-M_star)**2)/Spesi)
		I_g_star=sqrt(dot_product(w*Gauss,(g_i-g_star)**2)/Spesi)
		I_rho_star=sqrt(dot_product(w*Gauss,(rho_i-rho_star)**2)/Spesi)
		!!! 
		if (photIsocAvail) then 
			I_BmV_star=sqrt(dot_product(w*Gauss,(BmV_i-BmV_star)**2)/Spesi)
			I_BC_star=sqrt(dot_product(w*Gauss,(BC_i-BC_star)**2)/Spesi)
		end if 

		logg_star=log10(g_star)
		R_star=sqrt(L_star/(Teff_star/TeffSun)**4)
		I_logg_star=log10(e)*I_g_star/g_star
		I_R_star=R_star*(0.5*I_L_star/L_star+2.*I_Teff_star/Teff_star)
				
		if (M_star>=MinfED .and. M_star<=MsupED .and.elDiff) then !elDiff=.false. to switch-off element diffusion
			Mndx=minloc(abs(M_star-MeAv),1)
			allocate(Ztevo(Endxf(Mndx)-Endxi(Mndx)+1,size(ZtevoTab,2)))
			Ztevo=ZtevoTab(Endxi(Mndx):Endxf(Mndx),:)
			!!Polynomials interpolating the Z=Z(t) have degree=3
			allocate(Zini(size(Ztevo,1))); allocate(tsp(size(Ztevo,1))); allocate(tmax(size(Ztevo,1)))
			allocate(a1(size(Ztevo),4)); allocate(a2(size(Ztevo),4))
			Zini=Ztevo(:,1)
			tsp=Ztevo(:,2)
			tmax=Ztevo(:,3)
			a1=Ztevo(:,4:7)
			a2=Ztevo(:,8:11)
			t_starvec=(/1.D0,t_star,t_star**2,t_star**3/)
			!
			allocate(Zguess(size(Ztevo,1)))
			!
			do iZg=1,size(Ztevo,1) 
				if (t_star>tmax(iZg)) then 
					Zguess(iZg)=Zini(iZg)
				else if (tsp(iZg)==-1) then 
					Zguess(iZg)=dot_product(a1(iZg,:),t_starvec)
				else
					if (t_star<=tsp(iZg)) then 
						Zguess(iZg)=dot_product(a1(iZg,:),t_starvec)
					else
						Zguess(iZg)=dot_product(a2(iZg,:),t_starvec)
					end if 
				end if 
			end do
			deallocate(Ztevo)
			allocate(indxZg(size(Zguess))) 
			call indexx(abs(Zguess-Z_),indxZg) 
			Zg1=Zguess(indxZg(1))
			Zg2=Zguess(indxZg(2))
			rZ=(Z_-Zg1)/(Zg2-Zg1)
			Zgini1=Zini(indxZg(1))
			Zgini2=Zini(indxZg(2))
			Zgini=Zgini1+rZ*(Zgini2-Zgini1)
			indxZgini=minloc(abs(Zgini-Zvec),1) 
			Z2=Zvec(indxZgini)
			
			deallocate(Zini); deallocate(tsp); deallocate(tmax)
			deallocate(a1); deallocate(a2); deallocate(indxZg)
			if (Z2.ne.Z) then
				deallocate(Isoc)
				deallocate(t_iso); deallocate(logt_iso)
				deallocate(M_iso); deallocate(logL_iso)
				deallocate(L_iso); deallocate(logTeff_iso)
				deallocate(Teff_iso); deallocate(logg_iso)
				deallocate(g_iso); deallocate(R_iso)
				deallocate(rho_iso); deallocate(logrho_iso)
				if (allocated(BmV_iso)) then
					deallocate(BmV_iso); deallocate(V_iso)
					deallocate(BC_iso)
				end if
				deallocate(t_i); deallocate(dist)
				deallocate(logTeff_i); deallocate(logL_i); deallocate(M_i); deallocate(logg_i)
				deallocate(logrho_i); deallocate(Teff_i); deallocate(g_i); deallocate(L_i)
				deallocate(rho_i); deallocate(Zguess)
				if (allocated(y_i)) then
					deallocate(y_i)
				end if
				deallocate(Gauss); deallocate(w)
				if (photIsocAvail) then
					deallocate(BC_i); deallocate(BmV_i)
				end if
				
				deallocate(MTrAv)
				allocate(MTrAv(nM(indxZt(indxZgini))))
				MTrAv=Mav(indxZt(indxZgini),1:nM(indxZt(indxZgini)))
				Z=Z2
				FeH=log10(Z)+costZ
				Yini=dnint((He0+dYdZ*Z)*100)/100
				
				deallocate(TabVel)
				xVTi=Vndxi(indxZt(indxZgini))
				xVTf=Vndxf(indxZt(indxZgini))
				allocate(TabVel(xVTf-xVTi+1,size(velTrackTab,2)))
				TabVel=velTrackTab(xVTi:xVTf,:)
				
				deallocate(ZAMSz)
				xZAi=ZAndxi(indxZt(indxZgini))
				xZAf=ZAndxf(indxZt(indxZgini))
				allocate(ZAMSz(xZAf-xZAi+1,size(ZAMStab,2)))
				ZAMSz=ZAMStab(xZAi:xZAf,:)
				
				Zi=Zndxi(indxZgini)
				Zf=Zndxf(indxZgini)
				allocate(Isoc(Zf-Zi+1,size(IsocTab,2)))
				Isoc=IsocTab(Zi:Zf,:)
	
				row_i=size(Isoc,1)
				allocate(t_iso(row_i)); allocate(logt_iso(row_i))
				
				t_iso=Isoc(:,ct)
				logt_iso=dnint(log10(t_iso)*100)/100
				t_iso=10.**logt_iso
				first_logt=logt_iso(1)
				last_logt=logt_iso(row_i)

				allocate(M_iso(row_i)); allocate(logL_iso(row_i))
				allocate(L_iso(row_i)); allocate(logTeff_iso(row_i))
				allocate(Teff_iso(row_i)); allocate(logg_iso(row_i))
				allocate(g_iso(row_i)); allocate(R_iso(row_i))
				allocate(rho_iso(row_i)); allocate(logrho_iso(row_i))
		
				M_iso=Isoc(:,cM)
				logL_iso=Isoc(:,clogL)
				L_iso=10.**(logL_iso)
				logTeff_iso=Isoc(:,clogTe)
				Teff_iso=10.**Isoc(:,clogTe)
				logg_iso=Isoc(:,clogg)
				g_iso=10.**Isoc(:,clogg)
				R_iso=sqrt(L_iso/(Teff_iso/TeffSun)**4) !raggi solari
				rho_iso=M_iso/(R_iso**3)*rhoSun !g/cm3
				logrho_iso=log10(rho_iso)
				if (photIsocAvail) then
					allocate(V_iso(row_i)); allocate(BmV_iso(row_i))
					allocate(BC_iso(row_i))
					V_iso=Isoc(:,cmag0)
					if (isEq(useColor,1.D0,2)) then
						BmV_iso=Isoc(:,cmag2)-Isoc(:,cmag1)
					else
						BmV_iso=logTeff_iso
					end if
					!!! 
					BC_iso=Isoc(:,cmbol)-Isoc(:,cmag0)
					!!! 
				end if 

				!!!Selection of reference isochrones!!! 
				!!First evaluate which isochrones are the closest to the star to be analysed
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
				!!!!!!!NEW CODE 12/2/14!!!!!!!!!!!!!!! 
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
				allocate(dist_ndx(size(logTeff_iso))) !equal to any other size of _iso
				allocate(x_iso(size(logTeff_iso))); allocate(xx_iso(size(logTeff_iso)))
				if (calibHRD) then
					call indexx(sqrt((BmV-BmV_iso)**2+(Vass-V_iso)**2), dist_ndx) 
					x=BmV
					x_iso=BmV_iso
					xx_iso=logTeff_iso
					if (isEq(useColor,1.D0,2)) then
						if (idCol.eq.1) then
							cBrif=cBBmV
						else
							cBrif=cBTeff
						end if
					else
						cBrif=cBTeff
					end if
				end if 
				if (calibNoD) then
					if (.not.isEq(Rf,-1.D0,2)) then
						call indexx(sqrt((BmV-BmV_iso)**2+(logL-logL_iso)**2),dist_ndx)
					else if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
						& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))) then 
						call indexx(sqrt((BmV-BmV_iso)**2+(logg-logg_iso)**2),dist_ndx) 
					else !rhoAvail
						call indexx(sqrt((BmV-BmV_iso)**2+(logrho-logrho_iso)**2),dist_ndx) 
					end if 
					x=BmV
					x_iso=BmV_iso
					xx_iso=logTeff_iso
					cBrif=cBBmV
				end if 
				if (calibSPEC) then
					if (.not.isEq(Rf,-1.D0,2)) then
						call indexx(sqrt((logTeff-logTeff_iso)**2+(logL-logL_iso)**2),dist_ndx)
					else if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
						& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))) then 
						call indexx(sqrt((logTeff-logTeff_iso)**2+(logg-logg_iso)**2), dist_ndx) 
					else !rhoAvail
						call indexx(sqrt((logTeff-logTeff_iso)**2+(logrho-logrho_iso)**2), dist_ndx) 
					end if 
					x=logTeff
					xtau=x
					x_iso=logTeff_iso
					xx_iso=BmV_iso
					cBrif=cBTeff
				end if
				allocate(ageMinDist(size(logt_iso)))
				ageMinDist=logt_iso(dist_ndx)
				call uniqueFast(ageMinDist,2,ix,.true.)
				allocate(dist_ndxU(size(ix))) 
				dist_ndxU=dist_ndx(ix)
		
				deallocate(dist_ndx); deallocate(ageMinDist)
				deallocate(ix)
				rowAge=1
				!
				allocate(t_iTmp(size(dist_ndxU)));allocate(distTmp(size(dist_ndxU),2))
				allocate(logTeff_iTmp(size(dist_ndxU)));allocate(logL_iTmp(size(dist_ndxU)))
				allocate(M_iTmp(size(dist_ndxU)));allocate(logg_iTmp(size(dist_ndxU)))
				allocate(logrho_iTmp(size(dist_ndxU)))
				if (photIsocAvail) then 
					allocate(BmV_iTmp(size(dist_ndxU)));allocate(BC_iTmp(size(dist_ndxU)))
				end if
				!
				do Nage=1,size(dist_ndxU) 
					ndxDU=dist_ndxU(Nage)
					if (photIsocAvail) then 
						BmV_i1=BmV_iso(ndxDU)
						V_i1=V_iso(ndxDU)
					end if 
					logTeff_i1=logTeff_iso(ndxDU)
					logL_i1=logL_iso(ndxDU)
					M_i1=M_iso(ndxDU)
					logg_i1=logg_iso(ndxDU)
					logrho_i1=logrho_iso(ndxDU)
						
					call choose_i2Col(ndxDU,x,x_iso,xx_iso,calibSPEC,hstar,hstarlim,V_iso,logL_iso, &
						& M_iso,logg_iso,logrho_iso,logt_iso,BmV_i2,V_i2,logTeff_i2,logL_i2,M_i2,logg_i2,logrho_i2)
					if (calibSPEC) then
						x_i1=logTeff_i1
						x_i2=logTeff_i2
					else
						x_i1=BmV_i1
						x_i2=BmV_i2
					end if
					if (isEq(x_i1,x_i2,4)) then 
						cycle
					end if 
					t_iTmp(rowAge)=t_iso(ndxDU)
					if (photIsocAvail) then 
						BC_iTmp(rowAge)=Isoc(ndxDU,cmbol)-Isoc(ndxDU,cmag0) !just rough. 
						!When I recover BC in the logg-BmV plane, I'll make the interpolation using BC_i1 e BC_i2
					end if 
					if (calibHRD) then !!Condiz sulla calibrazione da fare 15/12/14
						call findX_i(BmV,Vass,BmV_i1,BmV_i2,V_i1,V_i2,x_i,dist0)
						BmV_iTmp(rowAge)=x_i
					end if 
					if (calibNoD) then
						if (.not.(isEq(Rf,-1.D0,2))) then
							call findX_i(BmV,logL,BmV_i1,BmV_i2,logL_i1,logL_i2,x_i,dist0)
						else if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
							& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))) then 
							call findX_i(BmV,logg,BmV_i1,BmV_i2,logg_i1,logg_i2,x_i,dist0)
						else !rhoAvail
							call findX_i(BmV,logrho,BmV_i1,BmV_i2,logrho_i1,logrho_i2,x_i,dist0)
						end if 
						BmV_iTmp(rowAge)=x_i
					end if 
					if (calibSPEC) then
						if (.not.(isEq(Rf,-1.D0,2))) then
							call findX_i(logTeff,logL,logTeff_i1,logTeff_i2,logL_i1,logL_i2,x_i,dist0)
						else if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
							& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))) then 
							call findX_i(logTeff,logg,logTeff_i1,logTeff_i2,logg_i1,logg_i2,x_i,dist0)
						else !rhoAvail
							call findX_i(logTeff,logrho,logTeff_i1,logTeff_i2,logrho_i1,logrho_i2,x_i,dist0)
						end if 
						logTeff_iTmp(rowAge)=x_i 
					end if
			
					distTmp(rowAge,:)=(/ dist0, logt_iso(ndxDU) /)
		
					!!18/3/15 The following if cancel the criterion of perpendicular distance
					!!between a star and the isochrone if the theoretical point to be interpolated
					!!doesn't fall inside the isochrone (i.e. between x_i1 and x_i2), but on its
					!!extension. This possibility, in fact, would build fictitious isochrones that
					!!could be erroneously close to the star
					if (.not.((x_i>=x_i1.and.x_i<=x_i2).or.(x_i>=x_i2.and.x_i<=x_i1))) then 
						if (photIsocAvail) then 
							BmV_iTmp(rowAge)=BmV_i1
						end if 
						logg_iTmp(rowAge)=logg_i1
						logTeff_iTmp(rowAge)=logTeff_i1
						logL_iTmp(rowAge)=logL_i1
						M_iTmp(rowAge)=M_i1
						logrho_iTmp(rowAge)=logrho_i1
					else
						if (calibSPEC) then 
							if (photIsocAvail) then 
								mBV=(BmV_i2-BmV_i1)/(x_i2-x_i1)
								qBV=-mBV*x_i1+BmV_i1
								BmV_iTmp(rowAge)=mBV*x_i+qBV
							end if 
						else
							mT=(logTeff_i2-logTeff_i1)/(x_i2-x_i1)
							qT=-mT*x_i1+logTeff_i1
							logTeff_iTmp(rowAge)=mT*x_i+qT
						end if 
						mg=(logg_i2-logg_i1)/(x_i2-x_i1)
						qg=-mg*x_i1+logg_i1
						logg_iTmp(rowAge)=mg*x_i+qg
						mL=(logL_i2-logL_i1)/(x_i2-x_i1)
						qL=-mL*x_i1+logL_i1
						logL_iTmp(rowAge)=mL*x_i+qL
						mM=(M_i2-M_i1)/(x_i2-x_i1)
						qM=-mM*x_i1+M_i1
						M_iTmp(rowAge)=mM*x_i+qM
						mrh=(logrho_i2-logrho_i1)/(x_i2-x_i1)
						qrho=-mrh*x_i1+logrho_i1
						logrho_iTmp(rowAge)=mrh*x_i+qrho
					end if 

					rowAge=rowAge+1
				end do
				deallocate(x_iso); deallocate(xx_iso)
		
				allocate(t_i(rowAge-1));allocate(dist(rowAge-1,2))
				allocate(logTeff_i(rowAge-1));allocate(logL_i(rowAge-1))
				allocate(M_i(rowAge-1));allocate(logg_i(rowAge-1))
				allocate(logrho_i(rowAge-1))
				allocate(BmV_i(rowAge-1));allocate(BC_i(rowAge-1))
				t_i=t_iTmp(1:rowAge-1)
				dist=distTmp(1:rowAge-1,:)
				logTeff_i=logTeff_iTmp(1:rowAge-1)
				logL_i=logL_iTmp(1:rowAge-1)
				M_i=M_iTmp(1:rowAge-1)
				logg_i=logg_iTmp(1:rowAge-1)
				logrho_i=logrho_iTmp(1:rowAge-1)
				BmV_i=BmV_iTmp(1:rowAge-1)
				BC_i=BC_iTmp(1:rowAge-1)
				deallocate(t_iTmp);deallocate(distTmp);deallocate(logTeff_iTmp)
				deallocate(logL_iTmp);deallocate(M_iTmp);deallocate(logg_iTmp)
				deallocate(logrho_iTmp)
				if (photIsocAvail) then
					deallocate(BmV_iTmp);deallocate(BC_iTmp)
				end if
		
				allocate(Teff_i(rowAge-1));allocate(L_i(rowAge-1))
				allocate(g_i(rowAge-1));allocate(rho_i(rowAge-1))
		
				Teff_i=10.**logTeff_i
				L_i=10.**logL_i
				g_i=10.**logg_i
				rho_i=10.**logrho_i
		
				deallocate(dist_ndxU)
		
				allocate(dsorted(size(dist)))
				dsorted=dist(:,1) !then it will be overwritten by sort so that it will be actually sorted
				call sort1(dsorted)
				allocate(ix_sorted(size(dist,1)))
				call indexx(dist(:,1), ix_sorted) 
				!!Avoid extremely young isochrones (t<10 Myr) in the following calibration process
				age_s1=1 
				age_iso1=dist(ix_sorted(age_s1),2)
				do while (age_iso1.lt.logtlimCal)
					age_s1=age_s1+1
					age_iso1=dist(ix_sorted(age_s1),2)
				end do
				age_s=age_s1+1
				if (dsorted(age_s1)<1.e-6) then 
					do while ((log10(dsorted(age_s))-log10(dsorted(age_s1))<4.5.or.dist(ix_sorted(age_s),2).lt.logtlimCal) &
							& .and. age_s<size(ix_sorted))
						age_s=age_s+1
					end do 
				else
					do while (dist(ix_sorted(age_s),2).lt.logtlimCal .and. age_s<size(ix_sorted))
						age_s=age_s+1
					end do
				end if
				deallocate(dsorted)
				age_iso2=dist(ix_sorted(age_s),2)
				if (age_iso1>age_iso2) then 
					tmp_age1=age_iso1
					age_iso1=age_iso2
					age_iso2=tmp_age1 !age_iso2 BECOMES the isoc with the minimum distance from the star
					agePivot=age_iso2
				else
					agePivot=age_iso1 !age_iso1 REMAINS the isoc with the minimum distance from the star
				end if 

				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
				!!!!!!!!!!!END NEW CODE 
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
		
				call selectIsoc(logt_iso,age_iso1,logtstep,last_logt,ti0,tf0)
				ti(1)=ti0
				tf(1)=tf0
				call selectIsoc(logt_iso,age_iso2,logtstep,last_logt,ti0,tf0)
				ti(2)=ti0
				tf(2)=tf0
						
				if (calibHRD) then !!15/12/14
					VMag=V-DM0
					deallocate(VMag_iso)
					allocate(VMag_iso(size(Isoc,1)))
					VMag_iso=Isoc(:,cmag0)
					BmVuguali=.false.
					do jd=1,2 
						allocate(VMag_isor(tf(jd)-ti(jd)+1)); allocate(logL_isor(tf(jd)-ti(jd)+1))
						allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(Diff_BmV(tf(jd)-ti(jd)+1))
						VMag_isor=VMag_iso(ti(jd):tf(jd))
						logL_isor=logL_iso(ti(jd):tf(jd))
						BmV_isor=BmV_iso(ti(jd):tf(jd))
						Diff_BmV=sqrt((BmV-BmV_isor)**2+(VMag-VMag_isor)**2) !!!!!!!!!!!!18/2
						if (jd==1) then
							id(jd)=minloc(Diff_BmV,1) 
							BmV_is1v(jd)=BmV_isor(id(jd))
						else
							allocate(ndxDv(size(Diff_BmV)))
							call indexx(Diff_BmV,ndxDv)
							rowDv=size(ndxDv)
							allocate(VMagD(size(VMag_isor)))
							VMagD=VMag_isor(ndxDv)
							kD=1
							id(jd)=ndxDv(kD)
							BmV_is1v(jd)=BmV_isor(id(jd))
							do while ((abs(VMag_is1v(1)-VMagD(kD))<=0.01 .or. abs(VMag_is2v(1)-VMagD(kD))<=0.01) &
									& .and. kD<rowDv .and. VMag<2)!!!0.01 instead of 0.1!18/2 
								kD=kD+1
								id(jd)=ndxDv(kD)
								BmV_is1v(jd)=BmV_isor(id(jd))
							end do
							deallocate(ndxDv); deallocate(VMagD)
						end if
						VMag_is1v(jd)=VMag_isor(id(jd))
						logL_is1v(jd)=logL_isor(id(jd))
						!!!
						call choose_i2xy1y2(BmV,BmV_isor,VMag_isor,logL_isor,id(jd),hstar,hstarlim, & !VMag,
							& BmV_is1v(jd),x_is2,VMag_is1v(jd),y1_is2,logL_is1v(jd),y2_is2)
						
						BmV_is2v(jd)=x_is2
						VMag_is2v(jd)=y1_is2
						logL_is2v(jd)=y2_is2
						!!!
						if (isEq(BmV_is1v(jd),BmV_is2v(jd),4)) then 	!!!if 10/2/14 -modified 7/11/2018
							BmVuguali=.true.
							VMag_isv(jd)=VMag
							call InterpLin_M(reshape((/VMag_is1v(jd),VMag_is2v(jd),logL_is1v(jd),logL_is2v(jd)/), &
								& (/2,2/)),VMag,1,1,(/2/),logL_isv(jd),xlow,ylow,xup,yup)
						else
							VMag_isv(jd)=(VMag_is2v(jd)-VMag_is1v(jd))/(BmV_is2v(jd)-BmV_is1v(jd))* &
									& (BmV-BmV_is1v(jd))+VMag_is1v(jd)
							logL_isv(jd)=(logL_is2v(jd)-logL_is1v(jd))/(BmV_is2v(jd)-BmV_is1v(jd))* &
									& (BmV-BmV_is1v(jd))+logL_is1v(jd)	
						end if
						deallocate(VMag_isor); deallocate(logL_isor)
						deallocate(BmV_isor); deallocate(Diff_BmV)
					end do
										
					dCV_isv=abs(VMag_isv(1)-VMag_isv(2))  !vertical dist between isoc in the Color-Mag diagr
					dCL_isv=abs(logL_isv(1)-logL_isv(2))  !vertical dist between isoc in the Color-logL diagr
					if ((dCV_isv>dlim .and. dCL_isv>dlim)) then  !condition added 17/2/14
		!				dCVv=abs(VMag-VMag_isv)  		  !vertical dist between the star and the 2 isoc
						if (VMag<VMag_isv(1) .and. VMag<VMag_isv(2)) then 
							statev=1
						else if (VMag>VMag_isv(1) .and. VMag>VMag_isv(2)) then 
							statev=2
						else
							statev=3
						end if
						call calibrateDiag(statev,dCVvlim,VMag_isv,logL_isv,VMag,logL) 
						 
						L=10.**logL
						BC=-2.5*log10(L/d**2)-V-0.23
						I_L=L*(0.4*log(10.)*0.03+2.*I_d/d)
						I_logL=I_L/L*log10(e)
						logLuguali=.false.
						if (isEq(useColor,1.D0,2)) then
							do jd=1,2
								allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(logT_isor(tf(jd)-ti(jd)+1))
								allocate(logL_isor(tf(jd)-ti(jd)+1)); allocate(Diff_logL(tf(jd)-ti(jd)+1))
						
								BmV_isor=BmV_iso(ti(jd):tf(jd))
								logT_isor=logTeff_iso(ti(jd):tf(jd))
								logL_isor=logL_iso(ti(jd):tf(jd))
								Diff_logL=sqrt((logL-logL_isor)**2+(BmV-BmV_isor)**2) !!!!!!!!!!!!!!!!18/2/14
								if (jd==1) then
									id(jd)=minloc(Diff_logL,1) 
									logL_is1(jd)=logL_isor(id(jd))
								else
									allocate(ndxD(size(Diff_logL)))
									call indexx(Diff_logL,ndxD) 
									rowD=size(ndxD)
									allocate(BmVD(size(BmV_isor)))
									BmVD=BmV_isor(ndxD)
									kD=1
									id(jd)=ndxD(kD)
									logL_is1(jd)=logL_isor(id(jd))
									do while ((abs(BmV_is1(1)-BmVD(kD))<=0.04 .or. abs(BmV_is2(1)-BmVD(kD))<=0.04) .and. &
											& kD<rowD .and. logL>1.5) 
										kD=kD+1
										id(jd)=ndxD(kD)
										logL_is1(jd)=logL_isor(id(jd))
									end do
									deallocate(ndxD); deallocate(BmVD)
								end if 
								BmV_is1(jd)=BmV_isor(id(jd))
								logT_is1(jd)=logT_isor(id(jd))
								!!!
								call choose_i2xy1y2(logL,logL_isor,BmV_isor,logT_isor,id(jd),hstar,hstarlim, & !BmV,
									& logL_is1(jd),x_is2,BmV_is1(jd),y1_is2,logT_is1(jd),y2_is2)
								logL_is2(jd)=x_is2
								BmV_is2(jd)=y1_is2
								logT_is2(jd)=y2_is2
								!!!
								if (isEq(logL_is1(jd),logL_is2(jd),4)) then 
									logLuguali=.true.
									BmV_is(jd)=BmV
									call InterpLin_M(reshape((/BmV_is1(jd),BmV_is2(jd),logT_is1(jd),logT_is2(jd)/), &
										& (/2,2/)),BmV,1,1,(/2/),logT_is(jd),xlow,ylow,xup,yup)
								else
									BmV_is(jd)=(logL-logL_is1(jd))*(BmV_is2(jd)-BmV_is1(jd))/ &
												& (logL_is2(jd)-logL_is1(jd))+BmV_is1(jd)
									logT_is(jd)=(logL-logL_is1(jd))*(logT_is2(jd)-logT_is1(jd))/ &
												& (logL_is2(jd)-logL_is1(jd))+logT_is1(jd)
								end if 
													
								deallocate(BmV_isor); deallocate(logT_isor)
								deallocate(logL_isor); deallocate(Diff_logL)
							end do
														
							dCL_is=abs(BmV_is(1)-BmV_is(2))		!dist between isoc in the 'CMD'
							dHR_is=abs(logT_is(1)-logT_is(2))	!dist between isoc in the HRD
							if (dCL_is>dlim .and. dHR_is>dlim) then !condition added 17/2/14
			!					dCL=abs(BmV-BmV_is)				!vector of distances between the star and the 2 isoc
								if (BmV<BmV_is(1) .and. BmV<BmV_is(2)) then 
									state=1
								else if (BmV>BmV_is(1) .and. BmV>BmV_is(2)) then 
									state=2
								else
									state=3
								end if
								call calibrateDiag(state,dCLlim,BmV_is,logT_is,BmV,logTeff)
								 
!!!								if (isnan(logTeff)) then 
!!!									cycle
!!!								end if 
							else if (age_s==size(ix_sorted)) then 
								cycle
							end if
						else
							logTeff=BmV
							dCL_is=2.*dlim
							dHR_is=2.*dlim
						end if
						if (idCol.eq.1) then
							xtau=BmV
						else
							xtau=logTeff
						end if
					else if (age_s==size(ix_sorted)) then 
						cycle
					else
						dCL_is=0.
						dHR_is=0.
					end if
					cycleW=0
					do while ((dCV_isv<=dlim .or. dCL_isv<=dlim .or. dCL_is<=dlim .or. dHR_is<=dlim) .and. &
							& age_s<size(ix_sorted))
						age_s=age_s+1
						do while (dist(ix_sorted(age_s),2).lt.logtlimCal .and. age_s<size(ix_sorted))
							age_s=age_s+1
						end do
						if (isEq(age_iso1,agePivot,2)) then 
							age_iso2=dist(ix_sorted(age_s),2)
							if (age_iso1>age_iso2) then 
								tmp_age1=age_iso1
								age_iso1=age_iso2
								age_iso2=tmp_age1
								agePivot=age_iso2
							end if 
						else
							age_iso1=dist(ix_sorted(age_s),2)
							if (age_iso1>age_iso2) then 
								tmp_age1=age_iso1
								age_iso1=age_iso2
								age_iso2=tmp_age1
								agePivot=age_iso1
							end if 
						end if 
				
						call selectIsoc(logt_iso,age_iso1,logtstep,last_logt,ti0,tf0)
						ti(1)=ti0
						tf(1)=tf0
						call selectIsoc(logt_iso,age_iso2,logtstep,last_logt,ti0,tf0)
						ti(2)=ti0
						tf(2)=tf0
				
						VMag=V-DM0
						VMag_iso=Isoc(:,cmag0)
						BmVuguali=.false.
						do jd=1,2
							allocate(VMag_isor(tf(jd)-ti(jd)+1)); allocate(logL_isor(tf(jd)-ti(jd)+1))
							allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(Diff_BmV(tf(jd)-ti(jd)+1))					
							VMag_isor=VMag_iso(ti(jd):tf(jd))
							logL_isor=logL_iso(ti(jd):tf(jd))
							BmV_isor=BmV_iso(ti(jd):tf(jd))
							Diff_BmV=sqrt((BmV-BmV_isor)**2+(VMag-VMag_isor)**2) !!!!!!!!!!!!!!!!!18/2/14
							if (jd==1) then 
								id(jd)=minloc(Diff_BmV,1) 
								BmV_is1v(jd)=BmV_isor(id(jd))
							else
								allocate(ndxDv(size(Diff_BmV)))
								call indexx(Diff_BmV,ndxDv) 
								rowDv=size(ndxDv)
								allocate(VMagD(size(VMag_isor)))
								VMagD=VMag_isor(ndxDv)
								kD=1
								id(jd)=ndxDv(kD)
								BmV_is1v(jd)=BmV_isor(id(jd))
								do while ((abs(VMag_is1v(1)-VMagD(kD))<=0.01 .or. abs(VMag_is2v(1)-VMagD(kD))<=0.01) &
										& .and. kD<rowDv .and. VMag<2)!0.01 instead of 0.1 18/2 
									kD=kD+1
									id(jd)=ndxDv(kD)
									BmV_is1v(jd)=BmV_isor(id(jd))
								end do
								deallocate(ndxDv); deallocate(VMagD)
							end if 
							VMag_is1v(jd)=VMag_isor(id(jd))
							logL_is1v(jd)=logL_isor(id(jd))
							!!!
							call choose_i2xy1y2(BmV,BmV_isor,VMag_isor,logL_isor,id(jd),hstar,hstarlim, & !VMag,
								& BmV_is1v(jd),x_is2,VMag_is1v(jd),y1_is2,logL_is1v(jd),y2_is2) !!!
							BmV_is2v(jd)=x_is2
							VMag_is2v(jd)=y1_is2
							logL_is2v(jd)=y2_is2
							!!!
							if (isEq(BmV_is1v(jd),BmV_is2v(jd),4)) then 	!!!if 10/2/14 -modified 7/11/2018
								BmVuguali=.true.
								VMag_isv(jd)=VMag
								call InterpLin_M(reshape((/VMag_is1v(jd),VMag_is2v(jd),logL_is1v(jd),logL_is2v(jd)/), &
									& (/2,2/)),VMag,1,1,(/2/),logL_isv(jd),xlow,ylow,xup,yup)
							else
								VMag_isv(jd)=(VMag_is2v(jd)-VMag_is1v(jd))/(BmV_is2v(jd)-BmV_is1v(jd))* &
										& (BmV-BmV_is1v(jd))+VMag_is1v(jd)
								logL_isv(jd)=(logL_is2v(jd)-logL_is1v(jd))/(BmV_is2v(jd)-BmV_is1v(jd))* &
										& (BmV-BmV_is1v(jd))+logL_is1v(jd)	
							end if
					
							deallocate(VMag_isor); deallocate(logL_isor)
							deallocate(BmV_isor); deallocate(Diff_BmV)
						end do 
						
						dCV_isv=abs(VMag_isv(1)-VMag_isv(2))	!vertical dist between isochrones in the CMD
						dCL_isv=abs(logL_isv(1)-logL_isv(2))	!vert dist between isoch in the Color-logL diagr
						if ((dCV_isv>dlim .and. dCL_isv>dlim)) then !condition added 17/2/14
		!					dCVv=abs(VMag-VMag_isv)					!vertical dist between the star and the 2 isoc
							if (VMag<VMag_isv(1) .and. VMag<VMag_isv(2)) then 
								statev=1
							else if (VMag>VMag_isv(1) .and. VMag>VMag_isv(2)) then 
								statev=2
							else
								statev=3
							end if
							call calibrateDiag(statev,dCVvlim,VMag_isv,logL_isv,VMag,logL)
							 
!!!							if (isnan(logL)) then 
!!!								cycle
!!!							end if 
							L=10.**logL
							BC=-2.5*log10(L/d**2)-V-0.23
							I_L=L*(0.4*log(10.)*0.03+2.*I_d/d)
							I_logL=I_L/L*log10(e)
							logLuguali=.false.
							if (isEq(useColor,1.D0,2)) then
								do jd=1,2
									allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(logT_isor(tf(jd)-ti(jd)+1))
									allocate(logL_isor(tf(jd)-ti(jd)+1)); allocate(Diff_logL(tf(jd)-ti(jd)+1))
									logL_isor=logL_iso(ti(jd):tf(jd))
									logT_isor=logTeff_iso(ti(jd):tf(jd))
									BmV_isor=BmV_iso(ti(jd):tf(jd))
									Diff_logL=sqrt((logL-logL_isor)**2+(BmV-BmV_isor)**2) !!!!!!!!!!!!!!!18/2
									if (jd==1) then 
										id(jd)=minloc(Diff_logL,1) 
										logL_is1(jd)=logL_isor(id(jd))
									else
										allocate(ndxD(size(Diff_logL)))
										call indexx(Diff_logL,ndxD) 
										rowD=size(ndxD)
										allocate(BmVD(size(BmV_isor)))
										BmVD=BmV_isor(ndxD)
										kD=1
										id(jd)=ndxD(kD)
										logL_is1(jd)=logL_isor(id(jd))
										do while ((abs(BmV_is1(1)-BmVD(kD))<=0.04 .or. abs(BmV_is2(1)-BmVD(kD))<=0.04) &
												& .and. kD<rowD .and. logL>1.5) 
											kD=kD+1
											id(jd)=ndxD(kD)
											logL_is1(jd)=logL_isor(id(jd))
										end do
										deallocate(ndxD); deallocate(BmVD)
									end if 
									BmV_is1(jd)=BmV_isor(id(jd))
									logT_is1(jd)=logT_isor(id(jd))
									!!!
									call choose_i2xy1y2(logL,logL_isor,BmV_isor,logT_isor,id(jd),hstar,hstarlim, & !BmV,
										& logL_is1(jd),x_is2,BmV_is1(jd),y1_is2,logT_is1(jd),y2_is2)
									logL_is2(jd)=x_is2
									BmV_is2(jd)=y1_is2
									logT_is2(jd)=y2_is2
									!!!
									if (isEq(logL_is1(jd),logL_is2(jd),4)) then 
										logLuguali=.true.
										BmV_is(jd)=BmV
										call InterpLin_M(reshape((/BmV_is1(jd),BmV_is2(jd),logT_is1(jd),logT_is2(jd)/), &
											& (/2,2/)),BmV,1,1,(/2/),logT_is(jd),xlow,ylow,xup,yup)
									else
										BmV_is(jd)=(logL-logL_is1(jd))*(BmV_is2(jd)-BmV_is1(jd))/ &
													& (logL_is2(jd)-logL_is1(jd))+BmV_is1(jd)
										logT_is(jd)=(logL-logL_is1(jd))*(logT_is2(jd)-logT_is1(jd))/ &
													& (logL_is2(jd)-logL_is1(jd))+logT_is1(jd)
									end if 
														
									deallocate(BmV_isor); deallocate(logT_isor)
									deallocate(logL_isor); deallocate(Diff_logL)
								end do 
								
								dCL_is=abs(BmV_is(1)-BmV_is(2))		!dist between isoc in the 'CMD'
								dHR_is=abs(logT_is(1)-logT_is(2))	!dist between isoc in the HRD
								if ((dCL_is>dlim .and. dHR_is>dlim)) then !condition added 17/2/14
			!						dCL=abs(BmV-BmV_is)					!vector of the distances between star-->2 isoc
									if (BmV<BmV_is(1) .and. BmV<BmV_is(2)) then 
										state=1
									else if (BmV>BmV_is(1) .and. BmV>BmV_is(2)) then 
										state=2
									else
										state=3
									end if
							
									call calibrateDiag(state,dCLlim,BmV_is,logT_is,BmV,logTeff)
							
!!!									if (isnan(logTeff)) then 
!!!										cycle
!!!									end if 
								else if (age_s==size(ix_sorted)) then 
									cycleW=1
									cycle
								end if
							else
								logTeff=BmV
								dCL_is=2.*dlim
								dHR_is=2.*dlim
							end if
							if (idCol.eq.1) then
								xtau=BmV
							else
								xtau=logTeff
							end if
						else if (age_s==size(ix_sorted)) then 
							cycleW=1
							cycle
						end if 
					end do
					if (cycleW.eq.1) then
						cycle
					end if
					
					Teff=10.**logTeff
					if (.not.calibSPEC .and. Teffinput.ne.-1) then !Teffinput condtion added 14/7/2015
						if (abs(Teff-Teffinput)>300) then !Input B-V totally inconsistent with spectroscopic Teff
							cycle
						end if 
					end if
					if (isEq(useColor,1.D0,2)) then
						I_Teff=0.01*Teff
						I_logTeff=0.01*log10(e)
					else
						I_logTeff=InclogTeff
						I_Teff=Teff*log(10.)*I_logTeff
					end if
					Rc=sqrt(L/(Teff/TeffSun)**4)
					I_Rc=Rc*(0.5*I_L/L+2.*I_Teff/Teff)
					if (isEq(Rf,-1.D0,2)) then
						R1=Rc
						I_R1=I_Rc	
					else
						Ri=Rf
						I_Ri=SCP(22)
						call wMean(Ri,I_Ri,Rc,I_Rc,R1,I_R1)
					end if
					if (gProxyAvail) then !modified 18/9/2017 !Re-check compatibility between rho and logg 17/10/17
						if (loggAvail0 .and. rhoAvail0) then 
							loggf=loggf0
							rhof=rhof0
							R2=g0/10.**loggSun*rhoSun/rho0
							I_R2=R2*(I_g0/g0+I_rho0/rho0)
							call wMean(R1,I_R1,R2,I_R2,R,I_R)
							
							M=(g0/10.**loggSun)**3*(rhoSun/rho0)**2
							I_M=M*(3.*I_g0/g0+2.*I_rho0/rho0)
!!							L=R**2*(Teff/TeffSun)**4
!!							I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
!!							logL=log10(L)
!!							I_logL=I_L/L*log10(e)
							!!Modified to account for the mass uncertainty 26/10/2018 
							if (.not.(M+I_M.lt.minval(MTrAv) .or. M-I_M.gt.maxval(MTrAv))) then !open tracks only if inside MassRange 
								call selectMfromTracks(MTrAv,M,I_M,Mvec)
								allocate(cumInt(size(Mvec)))
								jk=0
								do jj=1,size(Mvec)
									indxM=minloc(abs(Mvec(jj)-MTrAv),1)
									xMTi=Tndxi(indxZt(indxZ),indxM)
									xMTf=Tndxf(indxZt(indxZ),indxM)
									allocate(TrRhoG(xMTf-xMTi+1,size(TrackTab,2)))
									TrRhoG=TrackTab(xMTi:xMTf,:)
				
									allocate(MTr(size(TrRhoG,1))); allocate(RTr(size(TrRhoG,1)))
									allocate(logRhoTr(size(TrRhoG,1))); allocate(loggTr(size(TrRhoG,1)))
									
									MTr=TrRhoG(:,cM_T) !Msun
									RTr=10.**TrRhoG(:,clogR_T) !cm
									logRhoTr=log10(MTr/(RTr/RSun)**3*rhoSun)
									loggTr=log10(MTr/(RTr/RSun)**2)+loggSun
									call findYgivenX_v(TrRhoG,logg,loggTr,clogTe_T,logTeTrtmp)
									if (allocated(logTeTrtmp)) then
										jk=jk+1
										if (jk.eq.1) then
											allocate(logTeTr(size(logTeTrtmp)))
											logTeTr=logTeTrtmp
										else
											call append1D(logTeTr,logTeTrtmp)
										end if
										cumInt(jj)=size(logTeTr)
										deallocate(logTeTrtmp)
									else
										cumInt(jj)=0
									end if
									
									deallocate(TrRhoG)
									deallocate(MTr); deallocate(RTr)
									deallocate(logRhoTr); deallocate(loggTr)
								end do
								!
								logTeJ=logTeff !I don't use Johnson (1966) relation because
								! I've just calibrate the correct isochronal Teff
								if (allocated(logTeTr)) then 
									allocate(TeB(size(logTeTr)))
									TeB=((logTeTr>logTeJ-DTeJ).and.(logTeTr<logTeJ+DTeJ))
									call consistentM(TeB,cumInt,Mvec,M,agreeM)
									if (.not.any(TeB) .or. .not.agreeM) then !all elements of TeB are 0=>logg and rho inconsistent
										!or even if there's some TeB, however (g,rho) no fully consistent
										R=R1
										I_R=I_R1
										if (I_logg0>I_logrho0) then !logrho is best determined
											!Only using rho: determine M; re-determine logg (discard the input value) 
											loggf=-1.
											loggAvail=.false.
											loggAvailAS=.false.
											if (.not.isEq(Rf,-1.D0,2)) then
												call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
											else
												call computeMg(R,I_R,Teff,I_Teff,rho,I_rho,M,I_M,g,I_g,logg,I_logg)
											end if
											loggAvailCal=.true. !!logg has become available thanks to calibration
										else !logg is best determined
											!Only using logg: determine M; re-determine rho (discard the input value)
											rhof=-1.
											rhoAvail=.false.
											rhoAvailAS=.false.
											if (.not.isEq(Rf,-1.D0,2)) then
												call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
											else
												call computeMrho(R,I_R,Teff,I_Teff,g,I_g,M,I_M,rho,I_rho,logrho,I_logrho)
											end if
											rhoAvailCal=.true.
										end if 
									else !input logg and rho are consistent.
										 ! Compute M through a weighted mean between M=M(logg) and M=M(rho)
										Mlogg=g0/10.**loggSun*R**2
										I_Mlogg=Mlogg*(I_g0/g0+2.*I_R/R)
										Mrho=rho0/rhoSun*R**3 !Mo
										I_Mrho=Mrho*(I_rho0/rho0+3.*I_R/R)
										call wMean(Mlogg,I_Mlogg,Mrho,I_Mrho,M,I_M)
									end if
									deallocate(TeB); deallocate(logTeTr)
								else !Track has been opened, but no logTe compatible with logg were found
									R=R1
									I_R=I_R1
									if (I_logg0>I_logrho0) then !logrho is best determined
										!Only using rho: determine M; re-determine logg (discard the input value) 
										loggf=-1.
										loggAvail=.false.
										loggAvailAS=.false.
										if (.not.isEq(Rf,-1.D0,2)) then
											call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
										else
											call computeMg(R,I_R,Teff,I_Teff,rho,I_rho,M,I_M,g,I_g,logg,I_logg)
										end if
										loggAvailCal=.true. !!logg has become available thanks to calibration
									else !logg is best determined
										!Only using logg: determine M; re-determine rho (discard the input value)
										rhof=-1.
										rhoAvail=.false.
										rhoAvailAS=.false.
										if (.not.isEq(Rf,-1.D0,2)) then
											call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho, &
												& logrho,I_logrho)
										else
											call computeMrho(R,I_R,Teff,I_Teff,g,I_g,M,I_M,rho,I_rho,logrho,I_logrho)
										end if
										rhoAvailCal=.true.
									end if 
								end if
								deallocate(cumInt); deallocate(Mvec) 
							else !M is out of mass range of tracks
								R=R1
								I_R=I_R1
								if (I_logg0>I_logrho0) then !logrho is best determined
									!Only using rho: determine M; re-determine logg (discard the input value) 
									loggf=-1.
									loggAvail=.false.
									loggAvailAS=.false.
									if (.not.isEq(Rf,-1.D0,2)) then
										call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
									else
										call computeMg(R,I_R,Teff,I_Teff,rho,I_rho,M,I_M,g,I_g,logg,I_logg)
									end if
									loggAvailCal=.true. !!logg has become available thanks to calibration
								else !logg is best determined
									!Only using logg: determine M; re-determine rho (discard the input value) 
									rhof=-1.
									rhoAvail=.false.
									rhoAvailAS=.false.
									if (.not.isEq(Rf,-1.D0,2)) then
										call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho, &
											& I_logrho)
									else
										call computeMrho(R,I_R,Teff,I_Teff,g,I_g,M,I_M,rho,I_rho,logrho,I_logrho)
									end if
									rhoAvailCal=.true.
								end if 
							end if 
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
						else !only one among logg and rho is available
							R=R1
							I_R=I_R1
							if (loggAvail0) then !g is still the original input one because it's not been overwritten
								if (.not.isEq(Rf,-1.D0,2)) then
									call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
								else
									call computeMrho(R,I_R,Teff,I_Teff,g,I_g,M,I_M,rho,I_rho,logrho,I_logrho)
								end if
								rhoAvailCal=.true.
							else !rhoAvail0 visto che vale gProxyAvail
								if (.not.isEq(Rf,-1.D0,2)) then
									call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
								else
									call computeMg(R,I_R,Teff,I_Teff,rho,I_rho,M,I_M,g,I_g,logg,I_logg)
								end if
								loggAvailCal=.true. !!logg has become available thanks to calibration
							end if 
						end if !endif (loggAvail .and. rhoAvail)
					else
						R=R1
						I_R=I_R1
						if (.not.isEq(Rf,-1.D0,2)) then
							L=R**2*(Teff/TeffSun)**4
							I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
							logL=log10(L)
							I_logL=I_L/L*log10(e)
						end if
					end if !endif gProxy
				else if (calibNoD) then 
					allocate(logy_iso(size(Isoc,1)))
					if (.not.isEq(Rf,-1.D0,2)) then
						y=logL
						yl=-y
						logy_iso=logL_iso
						y_lim=-1.75
					else if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
						& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))) then  !18/9/2017
						y=logg
						yl=y
						logy_iso=logg_iso
						y_lim=3.
					else !rhoAvail
						y=logrho
						yl=y
						logy_iso=logrho_iso
						y_lim=-2.2
					end if 
					logYuguali=.false.
					do jd=1,2
						allocate(logy_isor(tf(jd)-ti(jd)+1)); allocate(logT_isor(tf(jd)-ti(jd)+1))
						allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(Diff_logy(tf(jd)-ti(jd)+1))
				
						logy_isor=logy_iso(ti(jd):tf(jd))
						logT_isor=logTeff_iso(ti(jd):tf(jd))
						BmV_isor=BmV_iso(ti(jd):tf(jd))
						Diff_logy=sqrt((y-logy_isor)**2+(BmV-BmV_isor)**2) !!!!!!!!!!!!!!!!18/2/14
						if (jd==1) then 
							id(jd)=minloc(Diff_logy,1) 
							logy_is1(jd)=logy_isor(id(jd))
						else
							allocate(ndxD(size(Diff_logy)))
							call indexx(Diff_logy,ndxD) 
							rowD=size(ndxD)
							allocate(BmVD(size(BmV_isor)))
							BmVD=BmV_isor(ndxD)
							kD=1
							id(jd)=ndxD(kD)
							logy_is1(jd)=logy_isor(id(jd))
							do while ((abs(BmV_is1(1)-BmVD(kD))<=0.04 .or. abs(BmV_is2(1)-BmVD(kD))<=0.04) &
									& .and. kD<rowD .and. yl<y_lim) 
								kD=kD+1
								id(jd)=ndxD(kD)
								logy_is1(jd)=logy_isor(id(jd))
							end do
							deallocate(ndxD); deallocate(BmVD)
						end if 
						BmV_is1(jd)=BmV_isor(id(jd))
						logT_is1(jd)=logT_isor(id(jd))
						!!!
						call choose_i2xy1y2(y,logy_isor,BmV_isor,logT_isor,id(jd),hstar,hstarlim, & !BmV,
							& logy_is1(jd),x_is2,BmV_is1(jd),y1_is2,logT_is1(jd),y2_is2) !!!
						logy_is2(jd)=x_is2
						BmV_is2(jd)=y1_is2
						logT_is2(jd)=y2_is2
						!!!
						if (isEq(logy_is1(jd),logy_is2(jd),4)) then 
							logYuguali=.true.
							BmV_is(jd)=BmV
							call InterpLin_M(reshape((/BmV_is1(jd),BmV_is2(jd),logT_is1(jd),logT_is2(jd)/), &
								& (/2,2/)),BmV,1,1,(/2/),logT_is(jd),xlow,ylow,xup,yup)
						else
							BmV_is(jd)=(y-logy_is1(jd))*(BmV_is2(jd)-BmV_is1(jd))/ &
										& (logy_is2(jd)-logy_is1(jd))+BmV_is1(jd)
							logT_is(jd)=(y-logy_is1(jd))*(logT_is2(jd)-logT_is1(jd))/ &
										& (logy_is2(jd)-logy_is1(jd))+logT_is1(jd)
						end if 
				
						deallocate(BmV_isor); deallocate(logy_isor);
						deallocate(logT_isor); deallocate(Diff_logy)
					end do 
												 
					dCL_is=abs(BmV_is(1)-BmV_is(2))		!dist between isoc in the BmV-logg Diagram
					dHR_is=abs(logT_is(1)-logT_is(2))	!dist between isoc in the logTeff-logg Diagram
					if ((dCL_is>dlim .and. dHR_is>dlim)) then !condition added 17/2/14
		!				dCL=abs(BmV-BmV_is)				!vector of distances between the star and the two isoc
						if (BmV<BmV_is(1) .and. BmV<BmV_is(2)) then 
							state=1
						else if (BmV>BmV_is(1) .and. BmV>BmV_is(2)) then 
							state=2
						else
							state=3
						end if
						!!
						call calibrateDiag(state,dClim,BmV_is,logT_is,BmV,logTeff)
						
						if (idCol.eq.1) then
							xtau=BmV
						else
							xtau=logTeff
						end if
						!! 
					else if (age_s==size(ix_sorted,1)) then 
						cycle
					else
						dCL_is=0.
						dHR_is=0.
					end if
					cycleW=0
					do while ((dCL_is<=dlim .or. dHR_is<=dlim) .and. age_s<size(ix_sorted,1)) 
						age_s=age_s+1
						do while (dist(ix_sorted(age_s),2).lt.logtlimCal .and. age_s<size(ix_sorted))
							age_s=age_s+1
						end do
						if (isEq(age_iso1,agePivot,2)) then 
							age_iso2=dist(ix_sorted(age_s),2)
							if (age_iso1>age_iso2) then 
								tmp_age1=age_iso1
								age_iso1=age_iso2
								age_iso2=tmp_age1
								agePivot=age_iso2
							end if 
						else
							age_iso1=dist(ix_sorted(age_s),2)
							if (age_iso1>age_iso2) then 
								tmp_age1=age_iso1
								age_iso1=age_iso2
								age_iso2=tmp_age1
								agePivot=age_iso1
							end if 
						end if 
				
						call selectIsoc(logt_iso,age_iso1,logtstep,last_logt,ti0,tf0)
						ti(1)=ti0
						tf(1)=tf0
						call selectIsoc(logt_iso,age_iso2,logtstep,last_logt,ti0,tf0)
						ti(2)=ti0
						tf(2)=tf0
						 
						logYuguali=.false.
						do jd=1,2
							allocate(logy_isor(tf(jd)-ti(jd)+1)); allocate(logT_isor(tf(jd)-ti(jd)+1))
							allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(Diff_logy(tf(jd)-ti(jd)+1))
				
							logy_isor=logy_iso(ti(jd):tf(jd))
							logT_isor=logTeff_iso(ti(jd):tf(jd))
							BmV_isor=BmV_iso(ti(jd):tf(jd))
							Diff_logy=sqrt((y-logy_isor)**2+(BmV-BmV_isor)**2) !!!!!!!!!!!!!!!!18/2/14
							if (jd==1) then 
								id(jd)=minloc(Diff_logy,1) 
								logy_is1(jd)=logy_isor(id(jd))
							else
								allocate(ndxD(size(Diff_logy)))
								call indexx(Diff_logy,ndxD) 
								rowD=size(ndxD)
								allocate(BmVD(size(BmV_isor)))
								BmVD=BmV_isor(ndxD)
								kD=1
								id(jd)=ndxD(kD)
								logy_is1(jd)=logy_isor(id(jd))
								do while ((abs(BmV_is1(1)-BmVD(kD))<=0.04 .or. abs(BmV_is2(1)-BmVD(kD))<=0.04) &
										& .and. kD<rowD .and. yl<y_lim) 
									kD=kD+1
									id(jd)=ndxD(kD)
									logy_is1(jd)=logy_isor(id(jd))
								end do
								deallocate(ndxD); deallocate(BmVD)
							end if 
							BmV_is1(jd)=BmV_isor(id(jd))
							logT_is1(jd)=logT_isor(id(jd))
							!!!
							call choose_i2xy1y2(y,logy_isor,BmV_isor,logT_isor,id(jd),hstar,hstarlim, & !BmV,
								& logy_is1(jd),x_is2,BmV_is1(jd),y1_is2,logT_is1(jd),y2_is2) !!!
							logy_is2(jd)=x_is2
							BmV_is2(jd)=y1_is2
							logT_is2(jd)=y2_is2
							!!!
							if (isEq(logy_is1(jd),logy_is2(jd),4)) then 
								logYuguali=.true.
								BmV_is(jd)=BmV
								call InterpLin_M(reshape((/BmV_is1(jd),BmV_is2(jd),logT_is1(jd),logT_is2(jd)/), &
									& (/2,2/)),BmV,1,1,(/2/),logT_is(jd),xlow,ylow,xup,yup)
							else
								BmV_is(jd)=(y-logy_is1(jd))*(BmV_is2(jd)-BmV_is1(jd))/ &
											& (logy_is2(jd)-logy_is1(jd))+BmV_is1(jd)
								logT_is(jd)=(y-logy_is1(jd))*(logT_is2(jd)-logT_is1(jd))/ &
											& (logy_is2(jd)-logy_is1(jd))+logT_is1(jd)
							end if 
				
							deallocate(BmV_isor); deallocate(logy_isor);
							deallocate(logT_isor); deallocate(Diff_logy)
						end do 
						
						dCL_is=abs(BmV_is(1)-BmV_is(2))		!dist between isoc in the BmV-logg Diagram
						dHR_is=abs(logT_is(1)-logT_is(2))	!dist between isoc in the logTeff-logg Diagram
						if ((dCL_is>dlim .and. dHR_is>dlim)) then !condition added 17/2/14
			!				dCL=abs(BmV-BmV_is)				!vector of distances between the star and the two isoc
							if (BmV<BmV_is(1) .and. BmV<BmV_is(2)) then 
								state=1
							else if (BmV>BmV_is(1) .and. BmV>BmV_is(2)) then 
								state=2
							else
								state=3
							end if
							!!
							call calibrateDiag(state,dClim,BmV_is,logT_is,BmV,logTeff)
							if (idCol.eq.1) then
								xtau=BmV
							else
								xtau=logTeff
							end if
							!! 
						else if (age_s==size(ix_sorted,1)) then 
							cycleW=1
							cycle
						end if 
					end do
					if (cycleW.eq.1) then
						cycle
					end if
					Teff=10.**logTeff
					if ((.not.calibSPEC .and. Teffinput.ne.-1)) then !Inside calibNoD => .not.calibSPEC is always true
						if (abs(Teff-Teffinput)>300) then !Inconsistency between BmV and spectrscopic Teff
							cycle
						end if 
					end if 
					I_Teff=0.01*Teff
					I_logTeff=0.01*log10(e)
					!!Compatibility
					if (loggAvail0.and.rhoAvail0) then 
						!!rho, logg both available. I determine input values for M, R 
						loggf=loggf0
						rhof=rhof0
						R2=g0/10.**loggSun*rhoSun/rho0
						I_R2=R2*(I_g0/g0+I_rho0/rho0)
						if (isEq(Rf,-1.D0,2)) then
							R=R2
							I_R=I_R2
						else
							call wMean(R1,I_R1,R2,I_R2,R,I_R)
						end if
						M=(g0/10.**loggSun)**3*(rhoSun/rho0)**2
						I_M=M*(3.*I_g0/g0+2.*I_rho0/rho0)
						L=R**2*(Teff/TeffSun)**4
						I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
						logL=log10(L)
						I_logL=I_L/L*log10(e)
						!!!!!!Check compatibility between rho and logg 16/10/17 
						!!Modified to account for the mass uncertainty 26/10/2018 
						if (.not.(M+I_M.lt.minval(MTrAv) .or. M-I_M.gt.maxval(MTrAv))) then !open tracks only if inside MassRange 
							call selectMfromTracks(MTrAv,M,I_M,Mvec)
							allocate(cumInt(size(Mvec)))
							jk=0
							do jj=1,size(Mvec)
								indxM=minloc(abs(Mvec(jj)-MTrAv),1)
								xMTi=Tndxi(indxZt(indxZ),indxM)
								xMTf=Tndxf(indxZt(indxZ),indxM)
								allocate(TrRhoG(xMTf-xMTi+1,size(TrackTab,2)))
								TrRhoG=TrackTab(xMTi:xMTf,:)
			
								allocate(MTr(size(TrRhoG,1))); allocate(RTr(size(TrRhoG,1)))
								allocate(logRhoTr(size(TrRhoG,1))); allocate(loggTr(size(TrRhoG,1)))
								
								MTr=TrRhoG(:,cM_T) !Msun
								RTr=10.**TrRhoG(:,clogR_T) !cm
								logRhoTr=log10(MTr/(RTr/RSun)**3*rhoSun)
								loggTr=log10(MTr/(RTr/RSun)**2)+loggSun
								call findYgivenX_v(TrRhoG,logg,loggTr,clogTe_T,logTeTrtmp)
								if (allocated(logTeTrtmp)) then
									jk=jk+1
									if (jk.eq.1) then
										allocate(logTeTr(size(logTeTrtmp)))
										logTeTr=logTeTrtmp
									else
										call append1D(logTeTr,logTeTrtmp)
									end if
									cumInt(jj)=size(logTeTr)
									deallocate(logTeTrtmp)
								else
									cumInt(jj)=0
								end if
								
								deallocate(TrRhoG)
								deallocate(MTr); deallocate(RTr)
								deallocate(logRhoTr); deallocate(loggTr)
							end do
							logTeJ=logTeff !I don't use Johnson (1966) relation because I've just calibrated Teff
							DTeJ=I_logTeff !directly set to the true uncertainty
							if (allocated(logTeTr)) then 
								allocate(TeB(size(logTeTr)))
								TeB=((logTeTr>logTeJ-DTeJ).and.(logTeTr<logTeJ+DTeJ))
								if (.not.any(TeB)) then !all elements of TeB are 0 => logg and rho inconsistent
									if (I_logg0>I_logrho0) then !logrho is best determined
										if (isEq(Rf,-1.D0,2)) then
											call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
											call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
											loggf=-1.
											call setNaN(logg); call setNaN(I_logg)
											call setNaN(g); call setNaN(I_g)
											loggAvail=.false.
											loggAvailAS=.false.
										else !input R is available
											!Only using rho: determine M; re-determine logg (discard the input value) 
											loggf=-1.
											loggAvail=.false.
											loggAvailAS=.false.
											R=R1 !use just the input reliable value
											I_R=I_R1
											call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
											loggAvailCal=.true. !!logg now available thanks to calibration
										end if
									else !logg is best determined
										if (isEq(Rf,-1.D0,2)) then
											call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
											call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
											rhof=-1.
											call setNaN(rho); call setNaN(I_rho)
											call setNaN(logrho); call setNaN(I_logrho)
											rhoAvail=.false.
											rhoAvailAS=.false.
										else !input R is available
											!Only using logg: determine M; re-determine rho (discard the input value) 
											rhof=-1.
											rhoAvail=.false.
											rhoAvailAS=.false.
											R=R1 !use just the input reliable value
											I_R=I_R1
											call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
											rhoAvailCal=.true.
										end if
									end if 
								else
									!!!!!!!!!
									call consistentM(TeB,cumInt,Mvec,M,agreeM)
									if (.not.agreeM) then !new M has been recomputed inside the subroutine in case agreeM is false
														  !This M will be used in case R is not available directly from input
										if (I_logg0.gt.I_logrho0) then !logrho is best determined
											if (isEq(Rf,-1.D0,2)) then
												loggAvail=.false.
												loggAvailAS=.false.
												!g: mantain original input uncertainty
												call computeLRgfromMrho(M,I_M,rho,I_rho,Teff,I_Teff,L,logL,I_L,I_logL,R,I_R,g,logg)
												loggAvailCal=.true.
											else !input R is available
												!Only using rho: determine M; re-determine logg (discard the input value) 
												loggf=-1.
												loggAvail=.false.
												loggAvailAS=.false.
												R=R1
												I_R=I_R1
												call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
												loggAvailCal=.true. !!logg now available thanks to calibration
											end if
										else
											if (isEq(Rf,-1.D0,2)) then
												rhoAvail=.false.
												rhoAvailAS=.false.
												!rho in g/cm3. Mantain original input uncertainty on rho
												call computeLRrhofromMg(M,I_M,g,I_g,Teff,I_Teff,L,logL,I_L,I_logL,R,I_R,rho,logrho)
												rhoAvailCal=.true.
											else !input R is available
												!Only using logg: determine M; re-determine rho (discard the input value) 
												rhof=-1.
												rhoAvail=.false.
												rhoAvailAS=.false.
												R=R1
												I_R=I_R1
												call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
												rhoAvailCal=.true.
											end if
										end if
										if (isEq(Rf,-1.D0,2)) then
											R=g/10**loggSun*rhoSun/rho
											I_R=R*(I_g/g+I_rho/rho)
											L=R**2*(Teff/TeffSun)**4
											I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
											logL=log10(L)
											I_logL=I_L/L*log10(e)
										end if
									else
										if (I_logg0.lt.I_logrho0) then
											bestLogg=.true.
										else
											bestLogg=.false.
										end if
										if (.not.isEq(Rf,-1.D0,2)) then
										!input R is available: weighted mean to infer M
											Mlogg=g0/10.**loggSun*R**2
											I_Mlogg=Mlogg*(I_g0/g0+2.*I_R/R)
											Mrho=rho0/rhoSun*R**3 !Mo
											I_Mrho=Mrho*(I_rho0/rho0+3.*I_R/R)
											call wMean(Mlogg,I_Mlogg,Mrho,I_Mrho,M,I_M)
											!L already defined since the beginning from R
										end if
									end if
									!!!!!!!!!!!!!
								end if
								deallocate(logTeTr); deallocate(TeB)
							else !Track has been opened, but no logTe compatible with logg were found
								if (I_logg0>I_logrho0) then !logrho is best determined
									if (isEq(Rf,-1.D0,2)) then
										call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
										call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
										loggf=-1.
										call setNaN(logg); call setNaN(I_logg)
										call setNaN(g); call setNaN(I_g)
										loggAvail=.false.
										loggAvailAS=.false.
									else !input R is available
										!Only using rho: determine M; re-determine logg (discard the input value) 
										loggf=-1.
										loggAvail=.false.
										loggAvailAS=.false.
										R=R1 !use just the input reliable value
										I_R=I_R1
										call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
										loggAvailCal=.true. !!logg now available thanks to calibration
									end if
								else !logg is best determined
									if (isEq(Rf,-1.D0,2)) then
										call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
										call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
										rhof=-1.
										call setNaN(rho); call setNaN(I_rho)
										call setNaN(logrho); call setNaN(I_logrho)
										rhoAvail=.false.
										rhoAvailAS=.false.
									else !input R is available
										!Only using logg: determine M; re-determine rho (discard the input value) 
										rhof=-1.
										rhoAvail=.false.
										rhoAvailAS=.false.
										R=R1 !use just the input reliable value
										I_R=I_R1
										call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
										rhoAvailCal=.true.
									end if
								end if 
							end if
							deallocate(cumInt) 
						else
							if (I_logg0>I_logrho0) then !logrho is best determined
								if (isEq(Rf,-1.D0,2)) then
									call setNan(R); call setNan(I_R); call setNan(M); call setNan(I_M)
									call setNan(L); call setNan(I_L); call setNan(logL); call setNan(I_logL)
									loggf=-1.
									call setNan(logg); call setNan(I_logg)
									call setNan(g); call setNan(I_g)
									loggAvail=.false.
									loggAvailAS=.false.
								else !input R is available
									!Only using rho: determine M; re-determine logg (discard the input value) 
									loggf=-1.
									loggAvail=.false.
									loggAvailAS=.false.
									R=R1 !use just the input reliable value
									I_R=I_R1
									call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
									loggAvailCal=.true. !!logg now available thanks to calibration
								end if
							else !logg is best determined
								if (isEq(Rf,-1.D0,2)) then
									call setNan(R); call setNan(I_R); call setNan(M); call setNan(I_M)
									call setNan(L); call setNan(I_L); call setNan(logL); call setNan(I_logL)
									rhof=-1.
									call setNan(rho); call setNan(I_rho)
									call setNan(logrho); call setNan(I_logrho)
									rhoAvail=.false.
									rhoAvailAS=.false.
								else !input R is available
									!Only using logg: determine M; re-determine rho (discard the input value) 
									rhof=-1.
									rhoAvail=.false.
									rhoAvailAS=.false.
									R=R1 !use just the input reliable value
									I_R=I_R1
									call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
									rhoAvailCal=.true.
								end if
							end if 
						end if 
						!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
					else
						if (.not.isEq(Rf,-1.D0,2)) then
							R=R1
							I_R=I_R1
							if (loggAvail.or.loggAvailAS) then
								call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
								rhoAvailCal=.true.
							else if (rhoAvail.or.rhoAvailAS) then
								call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
								loggAvailCal=.true. !!logg available thanks to calibration
							else !noGproxy
								L=R**2*(Teff/TeffSun)**4
								I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
								logL=log10(L)
								I_logL=I_L/L*log10(e)
							end if
						end if
					end if
					!!End compatibility
					if (.not.isEq(Rf,-1.D0,2)) then
						kwl=0
						cycleW=0
						age_s=1
						do while (abs(logL-y).gt.DlogL.and.kwl.lt.10.and.age_s<size(ix_sorted))
							y=logL
							kwl=kwl+1
							deallocate(t_i); deallocate(rho_i); deallocate(dist);deallocate(ix_sorted)
							deallocate(logTeff_i); deallocate(logL_i); deallocate(M_i); deallocate(logg_i)
							deallocate(logrho_i); deallocate(Teff_i); deallocate(g_i); deallocate(L_i)
							deallocate(logy_iso)
							if (photIsocAvail) then
								deallocate(BC_i); deallocate(BmV_i)
							end if
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
							!!!!!!!NEW CODE 12/2/14!!!!!!!!!!!!!!! 
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
							allocate(dist_ndx(size(logTeff_iso))) !equal to any other size of _iso
							allocate(x_iso(size(logTeff_iso))); allocate(xx_iso(size(logTeff_iso)))
							call indexx(sqrt((BmV-BmV_iso)**2+(logL-logL_iso)**2),dist_ndx)
							x=BmV
							x_iso=BmV_iso
							xx_iso=logTeff_iso
							cBrif=cBBmV
							allocate(ageMinDist(size(logt_iso)))
							ageMinDist=logt_iso(dist_ndx)
							
							call uniqueFast(ageMinDist,2,ix,.true.)
							call sort0(ix) !so to display elements of ageMinDist in the oroginal order
							allocate(dist_ndxU(size(ix))) 
							dist_ndxU=dist_ndx(ix)
							
							deallocate(dist_ndx); deallocate(ageMinDist)
							deallocate(ix)
							rowAge=1
							! allocation of already defined variables
							allocate(t_iTmp(size(dist_ndxU)));allocate(distTmp(size(dist_ndxU),2))
							allocate(logTeff_iTmp(size(dist_ndxU)));allocate(logL_iTmp(size(dist_ndxU)))
							allocate(M_iTmp(size(dist_ndxU)));allocate(logg_iTmp(size(dist_ndxU)))
							allocate(logrho_iTmp(size(dist_ndxU)))
							if (photIsocAvail) then 
								allocate(BmV_iTmp(size(dist_ndxU)));allocate(BC_iTmp(size(dist_ndxU)))
							end if
							!
							do Nage=1,size(dist_ndxU) 
								ndxDU=dist_ndxU(Nage)
								if (photIsocAvail) then 
									BmV_i1=BmV_iso(ndxDU)
									V_i1=V_iso(ndxDU)
								end if 
								logTeff_i1=logTeff_iso(ndxDU)
								logL_i1=logL_iso(ndxDU)
								M_i1=M_iso(ndxDU)
								logg_i1=logg_iso(ndxDU)
								logrho_i1=logrho_iso(ndxDU)
								
								call choose_i2Col(ndxDU,x,x_iso,xx_iso,calibSPEC,hstar,hstarlim,V_iso,logL_iso, &
									& M_iso,logg_iso,logrho_iso,logt_iso,BmV_i2,V_i2,logTeff_i2,logL_i2,M_i2,logg_i2,logrho_i2)
								if (calibSPEC) then
									x_i1=logTeff_i1
									x_i2=logTeff_i2
								else
									x_i1=BmV_i1
									x_i2=BmV_i2
								end if
								if (isEq(x_i1,x_i2,4)) then 
									cycle
								end if 
								t_iTmp(rowAge)=t_iso(ndxDU)
								if (photIsocAvail) then 
									BC_iTmp(rowAge)=Isoc(ndxDU,cmbol)-Isoc(ndxDU,cmag0) !just rough.
									!When I recover BC in the logg-BmV plane, I'll make the interpolation using BC_i1 e BC_i2
								end if 
								
								call findX_i(BmV,logL,BmV_i1,BmV_i2,logL_i1,logL_i2,x_i,dist0)
								BmV_iTmp(rowAge)=x_i 
								
								distTmp(rowAge,:)=(/ dist0, logt_iso(ndxDU) /)
											
								!!18/3/15 The following if cancel the criterion of perpendicular distance
								!!between a star and the isochrone if the theoretical point to be interpolated
								!!doesn't fall inside the isochrone (i.e. between x_i1 and x_i2), but on its
								!!extension. This possibility, in fact, would build fictitious isochrones that
								!!could be erroneously close to the star
								if (.not.((x_i>=x_i1.and.x_i<=x_i2).or.(x_i>=x_i2.and.x_i<=x_i1))) then 
									if (photIsocAvail) then 
										BmV_iTmp(rowAge)=BmV_i1
									end if 
									logg_iTmp(rowAge)=logg_i1
									logTeff_iTmp(rowAge)=logTeff_i1
									logL_iTmp(rowAge)=logL_i1
									M_iTmp(rowAge)=M_i1
									logrho_iTmp(rowAge)=logrho_i1
								else
									mT=(logTeff_i2-logTeff_i1)/(x_i2-x_i1)
									qT=-mT*x_i1+logTeff_i1
									logTeff_iTmp(rowAge)=mT*x_i+qT
									
									mg=(logg_i2-logg_i1)/(x_i2-x_i1)
									qg=-mg*x_i1+logg_i1
									logg_iTmp(rowAge)=mg*x_i+qg
									mL=(logL_i2-logL_i1)/(x_i2-x_i1)
									qL=-mL*x_i1+logL_i1
									logL_iTmp(rowAge)=mL*x_i+qL
									mM=(M_i2-M_i1)/(x_i2-x_i1)
									qM=-mM*x_i1+M_i1
									M_iTmp(rowAge)=mM*x_i+qM
									mrh=(logrho_i2-logrho_i1)/(x_i2-x_i1)
									qrho=-mrh*x_i1+logrho_i1
									logrho_iTmp(rowAge)=mrh*x_i+qrho
								end if 

								rowAge=rowAge+1
							end do
							
							deallocate(x_iso); deallocate(xx_iso)
									
							allocate(t_i(rowAge-1));allocate(dist(rowAge-1,2))
							allocate(logTeff_i(rowAge-1));allocate(logL_i(rowAge-1))
							allocate(M_i(rowAge-1));allocate(logg_i(rowAge-1))
							allocate(logrho_i(rowAge-1))
							allocate(BmV_i(rowAge-1));allocate(BC_i(rowAge-1))
							t_i=t_iTmp(1:rowAge-1)
							dist=distTmp(1:rowAge-1,:)
							logTeff_i=logTeff_iTmp(1:rowAge-1)
							logL_i=logL_iTmp(1:rowAge-1)
							M_i=M_iTmp(1:rowAge-1)
							logg_i=logg_iTmp(1:rowAge-1)
							logrho_i=logrho_iTmp(1:rowAge-1)
							BmV_i=BmV_iTmp(1:rowAge-1)
							BC_i=BC_iTmp(1:rowAge-1)
							deallocate(t_iTmp);deallocate(distTmp);deallocate(logTeff_iTmp)
							deallocate(logL_iTmp);deallocate(M_iTmp);deallocate(logg_iTmp)
							deallocate(logrho_iTmp)
							if (photIsocAvail) then
								deallocate(BmV_iTmp);deallocate(BC_iTmp)
							end if
							
							allocate(Teff_i(rowAge-1));allocate(L_i(rowAge-1))
							allocate(g_i(rowAge-1));allocate(rho_i(rowAge-1))
							
							Teff_i=10.**logTeff_i
							L_i=10.**logL_i
							g_i=10.**logg_i
							rho_i=10.**logrho_i
							
							deallocate(dist_ndxU)
							
							allocate(dsorted(size(dist)))
							dsorted=dist(:,1) !then it will be overwritten by sort so that it will be actually sorted
							call sort1(dsorted)
							allocate(ix_sorted(size(dist,1)))
							call indexx(dist(:,1), ix_sorted) 
							!!Avoid extremely young isochrones (t<10 Myr) in the following calibration process
							age_s1=1 
							age_iso1=dist(ix_sorted(age_s1),2)
							do while (age_iso1.lt.logtlimCal)
								age_s1=age_s1+1
								age_iso1=dist(ix_sorted(age_s1),2)
							end do
							age_s=age_s1+1
							if (dsorted(age_s1)<1.e-6) then 
								do while ((log10(dsorted(age_s))-log10(dsorted(age_s1))<4.5.or.dist(ix_sorted(age_s),2).lt.logtlimCal) &
										& .and. age_s<size(ix_sorted))
									age_s=age_s+1
								end do 
							else
								do while (dist(ix_sorted(age_s),2).lt.logtlimCal .and. age_s<size(ix_sorted))
									age_s=age_s+1
								end do
							end if
							deallocate(dsorted)
							age_iso2=dist(ix_sorted(age_s),2)
							if (age_iso1>age_iso2) then 
								tmp_age1=age_iso1
								age_iso1=age_iso2
								age_iso2=tmp_age1 !age_iso2 BECOMES the isoc with the minimum distance from the star
								agePivot=age_iso2
							else
								agePivot=age_iso1 !age_iso1 REMAINS the isoc with the minimum distance from the star
							end if
							
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
							!!!!!!!!!!!END NEW CODE 
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
							
							call selectIsoc(logt_iso,age_iso1,logtstep,last_logt,ti0,tf0)
							ti(1)=ti0
							tf(1)=tf0
							call selectIsoc(logt_iso,age_iso2,logtstep,last_logt,ti0,tf0)
							ti(2)=ti0
							tf(2)=tf0
														
							allocate(logy_iso(size(Isoc,1)))
							
							y=logL !already stated before
							yl=-y
							logy_iso=logL_iso
							y_lim=-1.75
							
							logYuguali=.false.
							do jd=1,2
								allocate(logy_isor(tf(jd)-ti(jd)+1)); allocate(logT_isor(tf(jd)-ti(jd)+1))
								allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(Diff_logy(tf(jd)-ti(jd)+1))
								
								logy_isor=logy_iso(ti(jd):tf(jd))
								logT_isor=logTeff_iso(ti(jd):tf(jd))
								BmV_isor=BmV_iso(ti(jd):tf(jd))
								Diff_logy=sqrt((y-logy_isor)**2+(BmV-BmV_isor)**2) !!!!!!!!!!!!!!!!18/2/14
								if (jd==1) then 
									id(jd)=minloc(Diff_logy,1) 
									logy_is1(jd)=logy_isor(id(jd))
								else
									allocate(ndxD(size(Diff_logy)))
									call indexx(Diff_logy,ndxD) 
									rowD=size(ndxD)
									allocate(BmVD(size(BmV_isor)))
									BmVD=BmV_isor(ndxD)
									kD=1
									id(jd)=ndxD(kD)
									logy_is1(jd)=logy_isor(id(jd))
									do while ((abs(BmV_is1(1)-BmVD(kD))<=0.04 .or. abs(BmV_is2(1)-BmVD(kD))<=0.04) &
											& .and. kD<rowD .and. yl<y_lim) 
										kD=kD+1
										id(jd)=ndxD(kD)
										logy_is1(jd)=logy_isor(id(jd))
									end do
									deallocate(ndxD); deallocate(BmVD)
								end if 
								BmV_is1(jd)=BmV_isor(id(jd))
								logT_is1(jd)=logT_isor(id(jd))
								!!!
								call choose_i2xy1y2(y,logy_isor,BmV_isor,logT_isor,id(jd),hstar,hstarlim, & !BmV,
									& logy_is1(jd),x_is2,BmV_is1(jd),y1_is2,logT_is1(jd),y2_is2) !!!
								logy_is2(jd)=x_is2
								BmV_is2(jd)=y1_is2
								logT_is2(jd)=y2_is2
								!!!
								if (isEq(logy_is1(jd),logy_is2(jd),4)) then 
									logYuguali=.true.
									BmV_is(jd)=BmV
									call InterpLin_M(reshape((/BmV_is1(jd),BmV_is2(jd),logT_is1(jd),logT_is2(jd)/), &
										& (/2,2/)),BmV,1,1,(/2/),logT_is(jd),xlow,ylow,xup,yup)
								else
									BmV_is(jd)=(y-logy_is1(jd))*(BmV_is2(jd)-BmV_is1(jd))/(logy_is2(jd)-logy_is1(jd)) &
												& +BmV_is1(jd)
									logT_is(jd)=(y-logy_is1(jd))*(logT_is2(jd)-logT_is1(jd))/(logy_is2(jd)-logy_is1(jd)) &
												& +logT_is1(jd)
								end if
								
								deallocate(BmV_isor); deallocate(logy_isor);
								deallocate(logT_isor); deallocate(Diff_logy)
							end do 
									 
							dCL_is=abs(BmV_is(1)-BmV_is(2))		!dist between isoc in the BmV-logg Diagram
							dHR_is=abs(logT_is(1)-logT_is(2))	!dist between isoc in the logTeff-logg Diagram
							if ((dCL_is>dlim .and. dHR_is>dlim)) then !condition added 17/2/14
					!				dCL=abs(BmV-BmV_is)				!vector of distances between the star and the two isoc
								if (BmV<BmV_is(1) .and. BmV<BmV_is(2)) then 
									state=1
								else if (BmV>BmV_is(1) .and. BmV>BmV_is(2)) then 
									state=2
								else
									state=3
								end if
								!!
								call calibrateDiag(state,dClim,BmV_is,logT_is,BmV,logTeff)
								if (idCol.eq.1) then
									xtau=BmV
								else
									xtau=logTeff
								end if
								!! 
							else if (age_s==size(ix_sorted,1)) then 
								cycleW=1
								cycle
							else
								dCL_is=0.
								dHR_is=0.
							end if
							do while ((dCL_is<=dlim .or. dHR_is<=dlim) .and. age_s<size(ix_sorted,1)) 
								age_s=age_s+1
								do while (dist(ix_sorted(age_s),2).lt.logtlimCal .and. age_s<size(ix_sorted))
									age_s=age_s+1
								end do
								if (isEq(age_iso1,agePivot,2)) then 
									age_iso2=dist(ix_sorted(age_s),2)
									if (age_iso1>age_iso2) then 
										tmp_age1=age_iso1
										age_iso1=age_iso2
										age_iso2=tmp_age1
										agePivot=age_iso2
									end if 
								else
									age_iso1=dist(ix_sorted(age_s),2)
									if (age_iso1>age_iso2) then 
										tmp_age1=age_iso1
										age_iso1=age_iso2
										age_iso2=tmp_age1
										agePivot=age_iso1
									end if 
								end if 
								
								call selectIsoc(logt_iso,age_iso1,logtstep,last_logt,ti0,tf0)
								ti(1)=ti0
								tf(1)=tf0
								call selectIsoc(logt_iso,age_iso2,logtstep,last_logt,ti0,tf0)
								ti(2)=ti0
								tf(2)=tf0
								 
								logYuguali=.false.
								do jd=1,2
									allocate(logy_isor(tf(jd)-ti(jd)+1)); allocate(logT_isor(tf(jd)-ti(jd)+1))
									allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(Diff_logy(tf(jd)-ti(jd)+1))
								
									logy_isor=logy_iso(ti(jd):tf(jd))
									logT_isor=logTeff_iso(ti(jd):tf(jd))
									BmV_isor=BmV_iso(ti(jd):tf(jd))
									Diff_logy=sqrt((y-logy_isor)**2+(BmV-BmV_isor)**2) !!!!!!!!!!!!!!!!18/2/14
									if (jd==1) then 
										id(jd)=minloc(Diff_logy,1) 
										logy_is1(jd)=logy_isor(id(jd))
									else
										allocate(ndxD(size(Diff_logy)))
										call indexx(Diff_logy,ndxD) 
										rowD=size(ndxD)
										allocate(BmVD(size(BmV_isor)))
										BmVD=BmV_isor(ndxD)
										kD=1
										id(jd)=ndxD(kD)
										logy_is1(jd)=logy_isor(id(jd))
										do while ((abs(BmV_is1(1)-BmVD(kD))<=0.04 .or. abs(BmV_is2(1)-BmVD(kD))<=0.04) &
												& .and. kD<rowD .and. yl<y_lim) 
											kD=kD+1
											id(jd)=ndxD(kD)
											logy_is1(jd)=logy_isor(id(jd))
										end do
										deallocate(ndxD); deallocate(BmVD)
									end if 
									BmV_is1(jd)=BmV_isor(id(jd))
									logT_is1(jd)=logT_isor(id(jd))
									!!!
									call choose_i2xy1y2(y,logy_isor,BmV_isor,logT_isor,id(jd),hstar,hstarlim, & !BmV,
										& logy_is1(jd),x_is2,BmV_is1(jd),y1_is2,logT_is1(jd),y2_is2) !!!
									logy_is2(jd)=x_is2
									BmV_is2(jd)=y1_is2
									logT_is2(jd)=y2_is2
									!!!
									if (isEq(logy_is1(jd),logy_is2(jd),4)) then 
										logYuguali=.true.
										BmV_is(jd)=BmV
										call InterpLin_M(reshape((/BmV_is1(jd),BmV_is2(jd),logT_is1(jd),logT_is2(jd)/), &
											& (/2,2/)),BmV,1,1,(/2/),logT_is(jd),xlow,ylow,xup,yup)
									else
										BmV_is(jd)=(y-logy_is1(jd))*(BmV_is2(jd)-BmV_is1(jd))/(logy_is2(jd)-logy_is1(jd)) &
													& +BmV_is1(jd)
										logT_is(jd)=(y-logy_is1(jd))*(logT_is2(jd)-logT_is1(jd))/(logy_is2(jd)-logy_is1(jd)) &
													& +logT_is1(jd)
									end if 
								
									deallocate(BmV_isor); deallocate(logy_isor);
									deallocate(logT_isor); deallocate(Diff_logy)
								end do 
								
								dCL_is=abs(BmV_is(1)-BmV_is(2))		!dist between isoc in the BmV-logg Diagram
								dHR_is=abs(logT_is(1)-logT_is(2))	!dist between isoc in the logTeff-logg Diagram
								if ((dCL_is>dlim .and. dHR_is>dlim)) then !condition added 17/2/14
					!				dCL=abs(BmV-BmV_is)				!vector of distances between the star and the two isoc
									if (BmV<BmV_is(1) .and. BmV<BmV_is(2)) then 
										state=1
									else if (BmV>BmV_is(1) .and. BmV>BmV_is(2)) then 
										state=2
									else
										state=3
									end if
									!!
									call calibrateDiag(state,dClim,BmV_is,logT_is,BmV,logTeff)
									if (idCol.eq.1) then
										xtau=BmV
									else
										xtau=logTeff
									end if
									!! 
								else if (age_s==size(ix_sorted,1)) then 
									cycleW=1
									cycle
								end if 
							end do
							if (cycleW.eq.1) then
								cycle
							end if
							Teff=10.**logTeff
							if ((.not.calibSPEC .and. Teffinput.ne.-1)) then !Inside calibNoD => .not.calibSPEC is always true
								if (abs(Teff-Teffinput)>300) then !Inconsistency between BmV and spectrscopic Teff
									cycle
								end if 
							end if 
							I_Teff=0.01*Teff
							I_logTeff=0.01*log10(e)
							!!Compatibility
							if (loggAvail0.and.rhoAvail0) then 
								!!rho, logg both available. I determine input values for M, R 
								loggf=loggf0
								rhof=rhof0
								R2=g0/10.**loggSun*rhoSun/rho0
								I_R2=R2*(I_g0/g0+I_rho0/rho0)
								if (isEq(Rf,-1.D0,2)) then
									R=R2
									I_R=I_R2
								else
									call wMean(R1,I_R1,R2,I_R2,R,I_R)
								end if
								M=(g0/10.**loggSun)**3*(rhoSun/rho0)**2
								I_M=M*(3.*I_g0/g0+2.*I_rho0/rho0)
								L=R**2*(Teff/TeffSun)**4
								I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
								logL=log10(L)
								I_logL=I_L/L*log10(e)
								!!!!!!Check compatibility between rho and logg 16/10/17 
								!!Modified to account for the mass uncertainty 26/10/2018 
								if (.not.(M+I_M.lt.minval(MTrAv) .or. M-I_M.gt.maxval(MTrAv))) then !open tracks only if inside MassRange 
									call selectMfromTracks(MTrAv,M,I_M,Mvec)
									allocate(cumInt(size(Mvec)))
									jk=0
									do jj=1,size(Mvec)
										indxM=minloc(abs(Mvec(jj)-MTrAv),1)
										xMTi=Tndxi(indxZt(indxZ),indxM)
										xMTf=Tndxf(indxZt(indxZ),indxM)
										allocate(TrRhoG(xMTf-xMTi+1,size(TrackTab,2)))
										TrRhoG=TrackTab(xMTi:xMTf,:)
					
										allocate(MTr(size(TrRhoG,1))); allocate(RTr(size(TrRhoG,1)))
										allocate(logRhoTr(size(TrRhoG,1))); allocate(loggTr(size(TrRhoG,1)))
										
										MTr=TrRhoG(:,cM_T) !Msun
										RTr=10.**TrRhoG(:,clogR_T) !cm
										logRhoTr=log10(MTr/(RTr/RSun)**3*rhoSun)
										loggTr=log10(MTr/(RTr/RSun)**2)+loggSun
										call findYgivenX_v(TrRhoG,logg,loggTr,clogTe_T,logTeTrtmp)
										if (allocated(logTeTrtmp)) then
											jk=jk+1
											if (jk.eq.1) then
												allocate(logTeTr(size(logTeTrtmp)))
												logTeTr=logTeTrtmp
											else
												call append1D(logTeTr,logTeTrtmp)
											end if
											cumInt(jj)=size(logTeTr)
											deallocate(logTeTrtmp)
										else
											cumInt(jj)=0
										end if
										
										deallocate(TrRhoG)
										deallocate(MTr); deallocate(RTr)
										deallocate(logRhoTr); deallocate(loggTr)
									end do
									logTeJ=logTeff !I don't use Johnson (1966) relation because I've just calibrated Teff
									DTeJ=I_logTeff !directly set to the true uncertainty
									if (allocated(logTeTr)) then 
										allocate(TeB(size(logTeTr)))
										TeB=((logTeTr>logTeJ-DTeJ).and.(logTeTr<logTeJ+DTeJ))
										if (.not.any(TeB)) then !all elements of TeB are 0 => logg and rho inconsistent
											if (I_logg0>I_logrho0) then !logrho is best determined
												if (isEq(Rf,-1.D0,2)) then
													call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
													call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
													loggf=-1.
													call setNaN(logg); call setNaN(I_logg)
													call setNaN(g); call setNaN(I_g)
													loggAvail=.false.
													loggAvailAS=.false.
												else !input R is available
													!Only using rho: determine M; re-determine logg (discard the input value) 
													loggf=-1.
													loggAvail=.false.
													loggAvailAS=.false.
													R=R1 !use just the input reliable value
													I_R=I_R1
													call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
													loggAvailCal=.true. !!logg now available thanks to calibration
												end if
											else !logg is best determined
												if (isEq(Rf,-1.D0,2)) then
													call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
													call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
													rhof=-1.
													call setNaN(rho); call setNaN(I_rho)
													call setNaN(logrho); call setNaN(I_logrho)
													rhoAvail=.false.
													rhoAvailAS=.false.
												else !input R is available
													!Only using logg: determine M; re-determine rho (discard the input value) 
													rhof=-1.
													rhoAvail=.false.
													rhoAvailAS=.false.
													R=R1 !use just the input reliable value
													I_R=I_R1
													call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
													rhoAvailCal=.true.
												end if
											end if 
										else
											!!!!!!!!!
											call consistentM(TeB,cumInt,Mvec,M,agreeM)
											if (.not.agreeM) then !new M has been recomputed inside the subroutine in case agreeM is false
																  !This M will be used in case R is not available directly from input
												if (I_logg0.gt.I_logrho0) then !logrho is best determined
													if (isEq(Rf,-1.D0,2)) then
														loggAvail=.false.
														loggAvailAS=.false.
														!g: mantain original input uncertainty
														call computeLRgfromMrho(M,I_M,rho,I_rho,Teff,I_Teff,L,logL,I_L,I_logL,R,I_R,g,logg)
														loggAvailCal=.true.
													else !input R is available
														!Only using rho: determine M; re-determine logg (discard the input value) 
														loggf=-1.
														loggAvail=.false.
														loggAvailAS=.false.
														R=R1
														I_R=I_R1
														call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
														loggAvailCal=.true. !!logg now available thanks to calibration
													end if
												else
													if (isEq(Rf,-1.D0,2)) then
														rhoAvail=.false.
														rhoAvailAS=.false.
														!rho in g/cm3. Mantain original input uncertainty on rho
														call computeLRrhofromMg(M,I_M,g,I_g,Teff,I_Teff,L,logL,I_L,I_logL,R,I_R,rho,logrho)
														rhoAvailCal=.true.
													else !input R is available
														!Only using logg: determine M; re-determine rho (discard the input value) 
														rhof=-1.
														rhoAvail=.false.
														rhoAvailAS=.false.
														R=R1
														I_R=I_R1
														call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
														rhoAvailCal=.true.
													end if
												end if
												if (isEq(Rf,-1.D0,2)) then
													R=g/10**loggSun*rhoSun/rho
													I_R=R*(I_g/g+I_rho/rho)
													L=R**2*(Teff/TeffSun)**4
													I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
													logL=log10(L)
													I_logL=I_L/L*log10(e)
												end if
											else
												if (I_logg0.lt.I_logrho0) then
													bestLogg=.true.
												else
													bestLogg=.false.
												end if
												if (.not.isEq(Rf,-1.D0,2)) then
												!input R is available: weighted mean to infer M
													Mlogg=g0/10.**loggSun*R**2
													I_Mlogg=Mlogg*(I_g0/g0+2.*I_R/R)
													Mrho=rho0/rhoSun*R**3 !Mo
													I_Mrho=Mrho*(I_rho0/rho0+3.*I_R/R)
													call wMean(Mlogg,I_Mlogg,Mrho,I_Mrho,M,I_M)
													!L already defined since the beginning from R
												end if
											end if
											!!!!!!!!!!!!!
										end if
										deallocate(logTeTr); deallocate(TeB)
									else !Track has been opened, but no logTe compatible with logg were found
										if (I_logg0>I_logrho0) then !logrho is best determined
											if (isEq(Rf,-1.D0,2)) then
												call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
												call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
												loggf=-1.
												call setNaN(logg); call setNaN(I_logg)
												call setNaN(g); call setNaN(I_g)
												loggAvail=.false.
												loggAvailAS=.false.
											else !input R is available
												!Only using rho: determine M; re-determine logg (discard the input value) 
												loggf=-1.
												loggAvail=.false.
												loggAvailAS=.false.
												R=R1 !use just the input reliable value
												I_R=I_R1
												call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
												loggAvailCal=.true. !!logg now available thanks to calibration
											end if
										else !logg is best determined
											if (isEq(Rf,-1.D0,2)) then
												call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
												call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
												rhof=-1.
												call setNaN(rho); call setNaN(I_rho)
												call setNaN(logrho); call setNaN(I_logrho)
												rhoAvail=.false.
												rhoAvailAS=.false.
											else !input R is available
												!Only using logg: determine M; re-determine rho (discard the input value) 
												rhof=-1.
												rhoAvail=.false.
												rhoAvailAS=.false.
												R=R1 !use just the input reliable value
												I_R=I_R1
												call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
												rhoAvailCal=.true.
											end if
										end if 
									end if
									deallocate(cumInt) 
								else
									if (I_logg0>I_logrho0) then !logrho is best determined
										if (isEq(Rf,-1.D0,2)) then
											call setNan(R); call setNan(I_R); call setNan(M); call setNan(I_M)
											call setNan(L); call setNan(I_L); call setNan(logL); call setNan(I_logL)
											loggf=-1.
											call setNan(logg); call setNan(I_logg)
											call setNan(g); call setNan(I_g)
											loggAvail=.false.
											loggAvailAS=.false.
										else !input R is available
											!Only using rho: determine M; re-determine logg (discard the input value) 
											loggf=-1.
											loggAvail=.false.
											loggAvailAS=.false.
											R=R1 !use just the input reliable value
											I_R=I_R1
											call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
											loggAvailCal=.true. !!logg now available thanks to calibration
										end if
									else !logg is best determined
										if (isEq(Rf,-1.D0,2)) then
											call setNan(R); call setNan(I_R); call setNan(M); call setNan(I_M)
											call setNan(L); call setNan(I_L); call setNan(logL); call setNan(I_logL)
											rhof=-1.
											call setNan(rho); call setNan(I_rho)
											call setNan(logrho); call setNan(I_logrho)
											rhoAvail=.false.
											rhoAvailAS=.false.
										else !input R is available
											!Only using logg: determine M; re-determine rho (discard the input value) 
											rhof=-1.
											rhoAvail=.false.
											rhoAvailAS=.false.
											R=R1 !use just the input reliable value
											I_R=I_R1
											call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
											rhoAvailCal=.true.
										end if
									end if 
								end if 
								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
							else
								if (.not.isEq(Rf,-1.D0,2)) then
									R=R1
									I_R=I_R1
									if (loggAvail.or.loggAvailAS) then
										call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
										rhoAvailCal=.true.
									else if (rhoAvail.or.rhoAvailAS) then
										call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
										loggAvailCal=.true. !!logg available thanks to calibration
									else !noGproxy
										L=R**2*(Teff/TeffSun)**4
										I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
										logL=log10(L)
										I_logL=I_L/L*log10(e)	
									end if
								end if
							end if
							!!End compatibility
						end do
						if (cycleW.eq.1) then
							cycle
						end if
					end if
					deallocate(logy_iso)
				end if
				deallocate(ix_sorted)

				!!!Check Activity/vsini.
				!If inactive stars or slow rotators, consider only isochrones with logt>=logtCutoff 
				check=0
				logtlim=-1
				if (.not.isEq(logRHK,0.D0,2)) then 
					logtMmod=c0-17.912*(logRHK+DlR)-1.6675*(logRHK+DlR)**2 !modified Mamayek
					if (logtMmod<logtCutoff5) then 
						logtlim1=logtMmod
					else
						logtlim1=logtCutoff5
					end if 
				else
					logtlim1=-1
				end if 

				!Searching T.O. of the oldest isochrone!!!!!!!!!!!!!!!!
				allocate(y_iso(size(Isoc,1)))
				caliblogL=calibHRD.or.((loggAvail.or.loggAvailAS.or.loggAvailCal).and.(rhoAvail.or.rhoAvailAS.or. & 
						& rhoAvailCal)).or.((calibSPEC.or.calibNoD).and..not.isEq(Rf,-1.D0,2))
				if (caliblogL) then !19/9/2017 !logL avail
					logY=logL
					y_iso=logL_iso
				else if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
					& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))) then 
					logY=-logg !to let following inequality logY< refer to lowMS stars
					y_iso=logg_iso
				else
					logY=-logrho !to let following inequality logY< refer to lowMS stars
					y_iso=logrho_iso
				end if 
				call searchTOold(y_iso,logTeff_iso,caliblogL,logt_iso,last_logt,logY_soglia)
		
				logY2=logY
				logY2_soglia=logY_soglia !consider these values later for low MS stars

				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

				if (.not.isEq(vsini,-1.D0,2) .or. .not.isEq(P,-1.D0,2)) then 
					vcheck=.true.
					!Rv: radius to be used in vsini relation 
					if (caliblogL) then 
						Rv=R
					else	!infer a rough estimate for stellar radius to be used for gyrochronology
						if (MminT<MlowMS) then !tracks of lower mass stars are available
							!define logY_sogliaVsini: open M=0.5Mo track; Z=Zstar and search logL/logg(last_logt)
							! that will be logY_sogliaVsini 
							indxM=minloc(abs(MlowMS-MTrAv),1)
							xMTi=Tndxi(indxZt(indxZgini),indxM)
							xMTf=Tndxf(indxZt(indxZgini),indxM)
							allocate(Track05(xMTf-xMTi+1,size(TrackTab,2)))
							Track05=TrackTab(xMTi:xMTf,:)
							MilowMS=MTrAv(indxM)
						
							!Establish maximim luminosity for a ~0.5solarMass star 
							if (caliblogL) then !Never enter this condition. Could only in the following two 
								call InterpLin_M( reshape((/Track05(:,ct_T),Track05(:,clogL_T)/), &
												& (/size(Track05,1),2/)),10.**last_logt,1,1,(/ 2/), &
												& logY_sogliaVsini1,xlow,ylow,xup,yup)
								logY_sogliaVsini=logY_sogliaVsini1(1)
							else if (loggAvail.or.loggAvailAS.or.loggAvailCal) then 
								call InterpLin_M( reshape((/Track05(:,ct_T),Track05(:,clogR_T)/), &
												& (/size(Track05,1),2/)),10.**last_logt,1,1,(/ 2/), &
												& logR_Tr,xlow,ylow,xup,yup)!cm
								R_Tr=10.**logR_Tr(1)
								logY_sogliaVsini=-(log10(MilowMS/(R_Tr/RSun)**2)+loggSun) !minus to let hold
												! inequality logY>logY_sogliaVsini
							else
								call InterpLin_M( reshape((/Track05(:,ct_T),Track05(:,clogR_T)/), &
												& (/size(Track05,1),2/)),10.**last_logt,1,1,(/ 2/), &
												& logR_Tr,xlow,ylow,xup,yup) !cm
								R_Tr=10.**logR_Tr(1)
								rho_Tr=MilowMS/(R_Tr/RSun)**3*rhoSun
								logY_sogliaVsini=-log10(rho_Tr)
							end if
							deallocate(Track05)
					
							if (logY<logY_sogliaVsini) then !very low MS star. R=R(t) almost constant
								!Select the middle-grid isochrone containing only the lines where M<MilowMS
								!and choose R in correspondence of logY 
								mid_n=nint((logt_halfMS-first_logt)/logtstep+1)
						
								call M2R(logt_iso,mid_n,Isoc,MilowMS,caliblogL,loggAvail.or.loggAvailAS,loggAvailCal,logY,Rv)
							else
								XT=logTeff-4.1
								if (rhoAvail.or.rhoAvailAS.or.rhoAvailCal) then
									!implementation of my relation R=R(logTeff,logrho,[Fe/H])
									vecT=(/ 1.D0,XT,XT**2,XT**3,logrho,logrho**2,logrho**3,FeH,FeH**2,FeH**3 /)
									Rv=10.**(dot_product(brho,vecT))
								else if (loggAvail.or.loggAvailAS.or.loggAvailCal) then 
									!implementation of my relation R=R(logTeff,logg,[Fe/H])
									vecT=(/ 1.D0,XT,XT**2,XT**3,logg,logg**2,logg**3,FeH,FeH**2,FeH**3 /)
									Rv=10.**(dot_product(blogg,vecT))
								else !Should never enter
									vcheck=.false.
								end if 
							end if 
						else
							XT=logTeff-4.1
							if (rhoAvail.or.rhoAvailAS.or.rhoAvailCal) then
								!implementation of my relation R=R(logTeff,logrho,[Fe/H])
								vecT=(/ 1.D0,XT,XT**2,XT**3,logrho,logrho**2,logrho**3,FeH,FeH**2,FeH**3 /)
								Rv=10.**(dot_product(brho,vecT))
							else if (loggAvail.or.loggAvailAS.or.loggAvailCal) then
								!implementation of my relation R=R(logTeff,logg,[Fe/H])
								vecT=(/ 1.D0,XT,XT**2,XT**3,logg,logg**2,logg**3,FeH,FeH**2,FeH**3 /)
								Rv=10.**(dot_product(blogg,vecT))
							else !Should never enter
								vcheck=.false.
							end if 
						end if 
					end if 
					if (vcheck.eqv..true.) then !if vsini or P available should be always true
						if (isEq(vsini,0.D0,2)) then !slow rotating stars for which vsini is reported as zero
							logtlim2=logtCutoff25 !no +0.05.
							!It's already been taken into account in defining logtCutoff25
						else
							if (.not.isEq(vsini,-1.D0,2).and.(isEq(P,-1.D0,2))) then 
								Omega=4./pi*vsini/(Rv*RSunKm)
								!logt_Deniss=log10(((OmSun/Omega)^2-1)*(1+A)/(2*(A/2+B))*(tSunLG-tZAMS)+tSunLG);
								!logtlim=logt_Deniss; 
								P=2.*pi/Omega/86400 !days
							end if
							call InterpLin_M(GyroTab,xtau,cBrif,1,(/ ctau /),tau1,xlow,taulow,xup,tauup)
							tau=tau1(1)
							if (.not.isEq(tau,0.D0,2) .and. x>=xlow .and. x<=xup) then 
								logt_Barnes=log10(tau/kC*log(P/P0)+kI/(2.*tau)*(P**2-P0**2))-3 ![Gyr]
								if (10.**logt_Barnes>DtBarnes) then 
									logt_Barnes=log10(10.**logt_Barnes-DtBarnes)+9 ![yrs]
								else
									logt_Barnes=first_logt
								end if 
								if (logt_Barnes<logtCutoff25) then 
									logtlim2=logt_Barnes
								else
									logtlim2=logtCutoff25
								end if 
							else
								logtlim2=-1 !tau=0 => Barnes relation not defined. I don't set any gyro age
							end if 
						end if 
					else
						logtlim2=-1
					end if 
				else
					logtlim2=-1
				end if 

				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
				if ((rhoAvail .or. rhoAvailAS .or. rhoAvailCal)) then !rhof.ne.-1
					if (logY<logY_soglia) then 
						if (.not.caliblogL) then
							!re-set correct values of logg (or logrho) and of logg_soglia (or logrho_soglia)
							logY=-logY
							logY_soglia=-logY_soglia
						end if 
						call trovaIndici(logt_iso,ndxrhoi,ndxrhof) 
						ndxTot=size(ndxrhoi) !=size(ndxrhof)
						jrho=1
						rhot=-100 !initialize like this so that it enters the cycle at least once
						do while (logt_iso(ndxrhoi(jrho))<8 .and. rho-I_rho>rhot .and. jrho<=ndxTot)
							!Consider only isochrones with t<100Myr
							!jrho<=ndxTot will always hold
							allocate(logY_isot(ndxrhof(jrho)-ndxrhoi(jrho)+1))
							allocate(rho_isot(ndxrhof(jrho)-ndxrhoi(jrho)+1))
							logY_isot=y_iso(ndxrhoi(jrho):ndxrhof(jrho))
							rho_isot=rho_iso(ndxrhoi(jrho):ndxrhof(jrho))
							ndxrho=minloc(abs(logY-logY_isot),1) 
							!logYt (useless) is the value reported by the isochrone t, that is closer to stellar logY
							!In correspondence of this value, I'll select the corresponding isochronal stellar density
		!					logYt=logY_isot(ndxrho)
							!rhotmp(jrho)=rho_isot(ndxrho) -->useless to consider it as a vector
					
							rhot=rho_isot(ndxrho)
							jrho=jrho+1
							deallocate(logY_isot); deallocate(rho_isot)
						end do
				
						if (jrho==2) then !While cycle done only once. Don't discard any isochrone
							logtlim3=logt_iso(1) !first age value reported by the isoc=>don't discard any isoc
						else
							logtlim3=logt_iso(ndxrhoi(jrho-1))+logtstep !added logtstep so that all isochrones with
																!logt<=logt_iso(ndxrhoi(ii)) are discarded. In fact
																!isochrones mantained start from logtlim included
							if (logtlim3>logtCutoff5) then 
								logtlim3=logtCutoff5
							end if 
						end if
						deallocate(ndxrhoi); deallocate(ndxrhof) 
					else
						logtlim3=-1
					end if
				else
					logtlim3=-1
				end if
				deallocate(y_iso)
				
				if (.not.isEq(YMg,-100.D0,2).and.FeH_>=-0.2 .and. FeH_<=0.2) then !metall condition
					t_Nissen=(YMg-aYMg)/bYMg*1.e9 ![yrs]
					if (t_Nissen>log10(DtNissen)) then
						logt_Nissen=log10(t_Nissen-DtNissen) !at least this age
					else
						logt_Nissen=first_logt
					end if
					if (logt_Nissen<logtCutoff25) then
						logtlim4=logt_Nissen
					else
						logtlim4=logtCutoff25
					end if
				else
					logtlim4=-1
				end if
								 
				limAges=(/ logtlim1,logtlim2,logtlim3,logtlim4 /)
				klim0=0
				do klim=1,size(limAges) 
					if (.not.isEq(limAges(klim),-1.D0,2)) then 
						klim0=1
						exit
					end if 
				end do
		
				if (klim0==1) then 
					logtlim=maxval(limAges)
					check=maxloc(limAges,1)
				else
					logtlim=-1
				end if 
				 
				if (.not.isEq(logtlim,-1.D0,2)) then !activity check done
					logtCutoff0=idnint(logtlim*100)
					resto=mod(logtCutoff0,idnint(logtstep*100))
					if (resto>floor(anint(logtstep*100)/2)) then 
						logtCutoff=(logtCutoff0+(dnint(logtstep*100)-resto))/100.
					else
						logtCutoff=(logtCutoff0-resto)/100.
					end if
					if (photIsocAvail) then !anint should be useless
						c_i=13
						allocate(Iso_i(size(t_i),c_i))
						Iso_i=reshape((/ dnint(log10(t_i)*100)/100,t_i,logTeff_i,Teff_i,logL_i,L_i,logg_i,g_i,M_i, &
									& logrho_i,rho_i,BmV_i,BC_i /),(/size(Iso_i,1),c_i/))
						!added rho_i 19/9/2017
					else
						c_i=11
						allocate(Iso_i(size(t_i),c_i))
						Iso_i=reshape((/ dnint(log10(t_i)*100)/100,t_i,logTeff_i,Teff_i,logL_i,L_i,logg_i,g_i,M_i, &
									& logrho_i,rho_i /),(/size(Iso_i,1),c_i/))
						!added rho_i 19/9/2017
					end if
					allocate(Iso_indx(size(Iso_i,1))) 
					call indexx(Iso_i(:,1),Iso_indx) 
					Iso_i=Iso_i(Iso_indx,:)
					deallocate(Iso_indx)
					do iIso=1,size(Iso_i,1) 
						if (isEq(Iso_i(iIso,1),logtCutoff,2)) then 
							iIso8=iIso
							exit
						end if 
					end do
					allocate(Iso_i2(size(Iso_i,1)-iIso8+1,c_i))
					Iso_i2=Iso_i(iIso8:size(Iso_i,1),:)
					deallocate(Iso_i)
					!dealloc _i and then re-alloc with proper dimension (lower than before in case of act check)
					deallocate(t_i); deallocate(logTeff_i); deallocate(Teff_i); deallocate(logL_i)
					deallocate(L_i); deallocate(logg_i); deallocate(g_i); deallocate(M_i)
					deallocate(logrho_i); deallocate(rho_i)
					allocate(t_i(size(Iso_i2,1))); allocate(logTeff_i(size(Iso_i2,1)))
					allocate(logL_i(size(Iso_i2,1))); allocate(L_i(size(Iso_i2,1))); allocate(logg_i(size(Iso_i2,1)))
					allocate(g_i(size(Iso_i2,1))); allocate(M_i(size(Iso_i2,1))); allocate(logrho_i(size(Iso_i2,1)))
					allocate(rho_i(size(Iso_i2,1))); allocate(Teff_i(size(Iso_i2,1)))
					t_i=Iso_i2(:,2)
					logTeff_i=Iso_i2(:,3)
					Teff_i=Iso_i2(:,4)
					logL_i=Iso_i2(:,5)
					L_i=Iso_i2(:,6)
					logg_i=Iso_i2(:,7)
					g_i=Iso_i2(:,8)
					M_i=Iso_i2(:,9)
					logrho_i=Iso_i2(:,10) !19/9/2017
					rho_i=Iso_i2(:,11) !19/9/2017
					if (photIsocAvail) then 
						deallocate(BmV_i); deallocate(BC_i)
						allocate(BmV_i(size(Iso_i2,1))); allocate(BC_i(size(Iso_i2,1)))
						BmV_i=Iso_i2(:,12) !19/9/2017
						BC_i=Iso_i2(:,13) !19/9/2017
					end if 
					logt8=logtCutoff
					deallocate(Iso_i2)
				else
					logt8=0
				end if
				!!!!End check activity/vsini 

				if (IncZ==1) then !load further metallic grids, accounting for I_FeH.
								  !Error bar extension=1 sigma and are "symmetric" in Z plane
					!both flag in multipleZisocPhSCP input parameters are set to 1
					chooseLogg=((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
							& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))
					if (I_FeH_>0) then !avoid that flags are considered as actual uncertainties
						if (calibHRD) then 
							xx=BmV
							yy=Vass
						end if 
						if (calibNoD) then 
							xx=BmV
							if (caliblogL) then
								yy=logL
							else if (chooseLogg) then 
								yy=logg
							else
								yy=logrho
							end if 
						end if 
						if (calibSPEC) then 
							xx=logTeff
							if (caliblogL) then
								yy=logL
							else if (chooseLogg) then 
								yy=logg
							else
								yy=logrho
							end if 
						end if 
						call multipleZisocPhSCP( FeH, I_FeH_, Zvec, model, xx, yy, &
							& logt8, percorso, 1.D0, 1, calibHRD, calibNoD, calibSPEC, chooseLogg, &
							& caliblogL,hstar, hstarlim, Isoc_i, Zlow, Zup,useColor,idCol)
						if (allocated(Isoc_i)) then !isochrone grids (WITH the reference one INCLUDED) are loaded
							deallocate(t_i); deallocate(logTeff_i); deallocate(Teff_i)
							deallocate(logL_i); deallocate(L_i); deallocate(logg_i)
							deallocate(g_i); deallocate(M_i); deallocate(logrho_i)
							deallocate(rho_i)
							allocate(t_i(size(Isoc_i,1))); allocate(logTeff_i(size(Isoc_i,1)))
							allocate(Teff_i(size(Isoc_i,1))); allocate(logL_i(size(Isoc_i,1)))
							allocate(L_i(size(Isoc_i,1))); allocate(logg_i(size(Isoc_i,1)))
							allocate(g_i(size(Isoc_i,1))); allocate(M_i(size(Isoc_i,1)))
							allocate(logrho_i(size(Isoc_i,1))); allocate(rho_i(size(Isoc_i,1)))
							t_i=Isoc_i(:,2)
							logTeff_i=Isoc_i(:,3)
							Teff_i=Isoc_i(:,4)
							logL_i=Isoc_i(:,5)
							L_i=Isoc_i(:,6)
							logg_i=Isoc_i(:,7)
							g_i=Isoc_i(:,8)
							M_i=Isoc_i(:,9)
							logrho_i=Isoc_i(:,10) !19/9/2017
							rho_i=Isoc_i(:,11) !19/9/2017
							if (photIsocAvail) then
								deallocate(BmV_i); deallocate(BC_i)
								allocate(BmV_i(size(Isoc_i,1))); allocate(BC_i(size(Isoc_i,1)))
								BmV_i=Isoc_i(:,12) !19/9/2017
								BC_i=Isoc_i(:,13) !19/9/2017
							end if
							deallocate(Isoc_i)
						end if 
					end if 
				end if 

				!!!only theoretical data _i
				allocate(Gauss(size(logTeff_i)))
				allocate(w(size(logTeff_i)))
				if (gProxyAvail) then 
					if (calibHRD .or. ((loggAvail.or.loggAvailAS.or.loggAvailCal) .and. &
						& (rhoAvail.or.rhoAvailAS.or.rhoAvailCal)).or.((calibSPEC.or.calibNoD).and..not.isEq(Rf,-1.D0,2))) then
						!I can recover Teff, L, M, logg
						if (logY2<logY2_soglia) then 
							!!!Find ZAMS in the Track:
							!!! logL starts increasing after the decreasing along the Hayashi line 
							call setThreshold(caliblogL,(((loggAvail.or.loggAvailAS).and.(rhoAvail.or.rhoAvailAS).and.bestLogg) &
									& .or.((loggAvail.or.loggAvailAS).and..not.(rhoAvail.or.rhoAvailAS))),logY2,logY2_soglia, &
									& logY0,logY0_soglia,cZAMS,cyT)
							allocate(yT_iso(size(Isoc,1)))
							if (cyT.eq.0) then !happens for logrho that hasn't a direct column in the isoch grid
								yT_iso=logrho_iso
							else
								yT_iso=Isoc(:,cyT)
							end if

							call ifComputeIsocAge(ZAMSz,cZAMS,logY0,logt_iso,logTeff_iso,yT_iso,I_logTeff,ageIsoc)
										
							deallocate(yT_iso)
						else
							ageIsoc=1
						end if 
				
						indxM=minloc(abs(M-MTrAv),1)
						xMTi=Tndxi(indxZt(indxZgini),indxM)
						xMTf=Tndxf(indxZt(indxZgini),indxM)
						allocate(Tracks(xMTf-xMTi+1,size(TrackTab,2)))
						Tracks=TrackTab(xMTi:xMTf,:)
	
						allocate(t_Tracks(size(Tracks,1)));allocate(logL_Tracks(size(Tracks,1)))
						allocate(logTeff_Tracks(size(Tracks,1)));allocate(vEvo(size(t_i)))
						t_Tracks=Tracks(:,ct_T)
						logL_Tracks=Tracks(:,clogL_T)
						logTeff_Tracks=Tracks(:,clogTe_T)
						deallocate(Tracks)
						call computeVevo(t_Tracks,logL_Tracks,logTeff_Tracks,t_i,varTrlim,vEvo)
						deallocate(t_Tracks);deallocate(logL_Tracks);deallocate(logTeff_Tracks)
				
						cRif=cVelL
						call computeVrif(TabVel,cRif,M,vRif)
						
						Gauss=1./(2.*pi*I_logTeff*I_logL)*e**(-0.5*((logTeff-logTeff_i)/I_logTeff)**2)* &
								& e**(-0.5*((logL-logL_i)/I_logL)**2)
						w=1./(((logL-logL_i)/I_logL)**2+((logTeff-logTeff_i)/I_logTeff)**2+ &
								& ((M-M_i)/I_M)**2+((logg-logg_i)/I_logg)**2+(log10(vRif/vEvo))**2)
						deallocate(vEvo)
					else
						if (logY2<logY2_soglia) then 
							call setThreshold(caliblogL,loggAvail.or.loggAvailAS,logY2,logY2_soglia,logY0,logY0_soglia,cZAMS,cyT)
							allocate(yT_iso(size(Isoc,1)))
							if (cyT.eq.0) then !happens for logrho that hasn't a direct column in the isoch grid
								yT_iso=logrho_iso
							else
								yT_iso=Isoc(:,cyT)
							end if
					
							call ifComputeIsocAge(ZAMSz,cZAMS,logY0,logt_iso,logTeff_iso,yT_iso,I_logTeff,ageIsoc)
										
							deallocate(yT_iso)
						else
							ageIsoc=1
						end if 
						
						allocate(y_i(size(logg_i)))
						if (loggAvail.or.loggAvailAS) then !only logg available
							y=logg
							I_y=I_logg
							y_i=logg_i
							cRif=cVelg
						else !only logrho available
							y=logrho
							I_y=I_logrho
							y_i=logrho_i
							cRif=cVelrho
						end if 
						!!!!!!!!!!!!!!!! 
						!!!!!!!!!!!!!!!!
						Gauss=1./(2.*pi*I_logTeff*I_y)*e**(-0.5*((logTeff-logTeff_i)/I_logTeff)**2)* &
								& e**(-0.5*((y-y_i)/I_y)**2)
						w=1./(((logTeff-logTeff_i)/I_logTeff)**2+((y-y_i)/I_y)**2)
						Spesi=sum(w*Gauss)
						if (Spesi<=1.D-200) then 
							cycle
						end if
				
						M_starTr=dot_product(w,M_i*Gauss)/Spesi
				
						indxM=minloc(abs(M_starTr-MTrAv),1)
						xMTi=Tndxi(indxZt(indxZgini),indxM)
						xMTf=Tndxf(indxZt(indxZgini),indxM)
						allocate(Tracks(xMTf-xMTi+1,size(TrackTab,2)))
						Tracks=TrackTab(xMTi:xMTf,:)
	
						allocate(t_Tracks(size(Tracks,1))); allocate(logTeff_Tracks(size(Tracks,1)))
						allocate(R_Tracks(size(Tracks,1))); allocate(logY_Tracks(size(Tracks,1)))
						allocate(vEvo(size(t_i)))
						t_Tracks=Tracks(:,ct_T)
						logTeff_Tracks=Tracks(:,clogTe_T)
						R_Tracks=10.**Tracks(:,clogR_T) !cm
						M_Tracks=Tracks(1,cM_T)
						deallocate(Tracks)
						if (loggAvail.or.loggAvailAS) then !only logg available
							logY_Tracks=log10(M_Tracks/(R_Tracks/RSun)**2)+loggSun
						else !only logrho available
							logY_Tracks=log10(M_Tracks/(R_Tracks/RSun)**3)+log10(rhoSun)
						end if
						call computeVevo(t_Tracks,logY_Tracks,logTeff_Tracks,t_i,varTrlim,vEvo)
						deallocate(t_Tracks); deallocate(logTeff_Tracks); deallocate(R_Tracks)
						deallocate(logY_Tracks)
				
						call computeVrif(TabVel,cRif,M_starTr,vRif)
				
						w=1./(((logTeff-logTeff_i)/I_logTeff)**2+((y-y_i)/I_y)**2+(log10(vRif/vEvo))**2)
						deallocate(vEvo)
						Spesi=sum(w*Gauss)
						if (Spesi<=1.D-200) then 
							cycle
						end if 
						M_starTr2=dot_product(w,M_i*Gauss)/Spesi
						kTr=1
						do while (abs((M_starTr2-M_starTr)/M_starTr)>DMTr .and. kTr<10) 
							if (logY2<logY2_soglia) then 
								call setThreshold(caliblogL,loggAvail.or.loggAvailAS,logY2,logY2_soglia,logY0,logY0_soglia,cZAMS,cyT)
								allocate(yT_iso(size(Isoc,1)))
								if (cyT.eq.0) then !happens for logrho that hasn't a direct column in the isoch grid
									yT_iso=logrho_iso
								else
									yT_iso=Isoc(:,cyT)
								end if
						
								call ifComputeIsocAge(ZAMSz,cZAMS,logY0,logt_iso,logTeff_iso,yT_iso, & 
													& I_logTeff,ageIsoc)
										
								deallocate(yT_iso)
							else
								ageIsoc=1
							end if 
					
							M_starTr=M_starTr2
							indxM=minloc(abs(M_starTr-MTrAv),1)
							xMTi=Tndxi(indxZt(indxZgini),indxM)
							xMTf=Tndxf(indxZt(indxZgini),indxM)
							allocate(Tracks(xMTf-xMTi+1,size(TrackTab,2)))
							Tracks=TrackTab(xMTi:xMTf,:)
	
							allocate(t_Tracks(size(Tracks,1))); allocate(logTeff_Tracks(size(Tracks,1)))
							allocate(R_Tracks(size(Tracks,1))); allocate(logY_Tracks(size(Tracks,1)));
							allocate(vEvo(size(t_i)))
							t_Tracks=Tracks(:,ct_T)
							logTeff_Tracks=Tracks(:,clogTe_T)
							R_Tracks=10.**Tracks(:,clogR_T) !cm
							M_Tracks=Tracks(1,cM_T)
							deallocate(Tracks)
							if (loggAvail.or.loggAvailAS) then !only logg available
								logY_Tracks=log10(M_Tracks/(R_Tracks/RSun)**2)+loggSun
							else !only logrho available
								logY_Tracks=log10(M_Tracks/(R_Tracks/RSun)**3)+log10(rhoSun)
							end if
							call computeVevo(t_tracks,logY_Tracks,logTeff_Tracks,t_i,varTrlim,Vevo)
							deallocate(t_Tracks); deallocate(logTeff_Tracks); deallocate(R_Tracks)
							deallocate(logY_Tracks)
					
							call computeVrif(TabVel,cRif,M_starTr,vRif)
					
							w=1./(((logTeff-logTeff_i)/I_logTeff)**2+((y-y_i)/I_y)**2+(log10(vRif/vEvo))**2)
							deallocate(vEvo)
							Spesi=sum(w*Gauss)
							if (Spesi<=1.D-200) then 
								cycle
							end if 
							M_starTr2=dot_product(w,M_i*Gauss)/Spesi
							kTr=kTr+1
						end do 
						!!!!!!!!!!!!!!!!!!! 
						!!!!!!!!!!!!!!!!!!! 
					end if 
				else
					if (logY2<logY2_soglia) then 
						call setThreshold(caliblogL,loggAvail.or.loggAvailAS,logY2,logY2_soglia,logY0,logY0_soglia,cZAMS,cyT)
						allocate(yT_iso(size(Isoc,1)))
						if (cyT.eq.0) then !happens for logrho that hasn't a direct column in the isoch grid
							yT_iso=logrho_iso
						else
							yT_iso=Isoc(:,cyT)
						end if
						 
						call ifComputeIsocAge(ZAMSz,cZAMS,logY0,logt_iso,logTeff_iso,yT_iso,I_logTeff,ageIsoc)
									
						deallocate(yT_iso)
					else
						ageIsoc=1
					end if 
					
					Gauss=1./(2.*pi*I_logTeff*I_logL)*e**(-0.5*((logTeff-logTeff_i)/I_logTeff)**2)* & 
							& e**(-0.5*((logL-logL_i)/I_logL)**2)
					w=1./(((logL-logL_i)/I_logL)**2+((logTeff-logTeff_i)/I_logTeff)**2)
					Spesi=sum(w*Gauss)
					if (Spesi<=1.D-200) then 
						cycle
					end if 
					M_starTr=dot_product(w,M_i*Gauss)/Spesi

					indxM=minloc(abs(M_starTr-MTrAv),1)
					xMTi=Tndxi(indxZt(indxZgini),indxM)
					xMTf=Tndxf(indxZt(indxZgini),indxM)
					allocate(Tracks(xMTf-xMTi+1,size(TrackTab,2)))
					Tracks=TrackTab(xMTi:xMTf,:)
					allocate(t_Tracks(size(Tracks,1)));allocate(logL_Tracks(size(Tracks,1)))
					allocate(logTeff_Tracks(size(Tracks,1)));allocate(vEvo(size(t_i)))
	
					t_Tracks=Tracks(:,ct_T)
					logL_Tracks=Tracks(:,clogL_T)
					logTeff_Tracks=Tracks(:,clogTe_T)
					deallocate(Tracks)
					call computeVevo(t_Tracks,logL_Tracks,logTeff_Tracks,t_i,varTrlim,vEvo)
					deallocate(t_Tracks);deallocate(logL_Tracks);deallocate(logTeff_Tracks) 
			
					cRif=cVelL
					call computeVrif(TabVel,cRif,M_starTr,vRif)
					w=1./(((logL-logL_i)/I_logL)**2+((logTeff-logTeff_i)/I_logTeff)**2+(log10(vRif/vEvo))**2)
					deallocate(vEvo)
					Spesi=sum(w*Gauss)
					if (Spesi<=1.D-200) then 
						cycle
					end if 
					M_starTr2=dot_product(w,M_i*Gauss)/Spesi
					kTr=1
					do while (abs((M_starTr2-M_starTr)/M_starTr)>DMTr .and. kTr<10) 
						if (logY2<logY2_soglia) then 
							call setThreshold(caliblogL,loggAvail.or.loggAvailAS,logY2,logY2_soglia,logY0,logY0_soglia,cZAMS,cyT)
							allocate(yT_iso(size(Isoc,1)))
							if (cyT.eq.0) then !happens for logrho that hasn't a direct column in the isoch grid
								yT_iso=logrho_iso
							else
								yT_iso=Isoc(:,cyT)
							end if
							 
							call ifComputeIsocAge(ZAMSz,cZAMS,logY0,logt_iso,logTeff_iso,yT_iso,I_logTeff,ageIsoc)
										
							deallocate(yT_iso)
						else
							ageIsoc=1
						end if 
				
						M_starTr=M_starTr2
						indxM=minloc(abs(M_starTr-MTrAv),1)
						xMTi=Tndxi(indxZt(indxZgini),indxM)
						xMTf=Tndxf(indxZt(indxZgini),indxM)
						allocate(Tracks(xMTf-xMTi+1,size(TrackTab,2)))
						Tracks=TrackTab(xMTi:xMTf,:)
						allocate(t_Tracks(size(Tracks,1)));allocate(logL_Tracks(size(Tracks,1)))
						allocate(logTeff_Tracks(size(Tracks,1)));allocate(vEvo(size(t_i)))
	
						t_Tracks=Tracks(:,ct_T)
						logL_Tracks=Tracks(:,clogL_T)
						logTeff_Tracks=Tracks(:,clogTe_T)
						deallocate(Tracks)
						call computeVevo(t_Tracks,logL_Tracks,logTeff_Tracks,t_i,varTrlim,vEvo)
						deallocate(t_Tracks);deallocate(logL_Tracks);deallocate(logTeff_Tracks) 
			
						cRif=cVelL
						call computeVrif(TabVel,cRif,M_starTr,vRif)
						w=1./(((logL-logL_i)/I_logL)**2+((logTeff-logTeff_i)/I_logTeff)**2+(log10(vRif/vEvo))**2)
						deallocate(vEvo)
						Spesi=sum(w*Gauss)
						if (Spesi<=1.D-200) then 
							cycle
						end if 
						M_starTr2=dot_product(w,M_i*Gauss)/Spesi
						kTr=kTr+1
					end do 
				end if
				Spesi=sum(w*Gauss)
				if (Spesi<=1.D-200) then 
					cycle
				end if 
				t_star2=dot_product(w,t_i*Gauss)/Spesi !filter the selected isochrone through a Gaussian
				!and then weight it through w. It's like as attributing it the weight w*Gauss
				Teff_star2=dot_product(w,Teff_i*Gauss)/Spesi
				L_star2=dot_product(w,L_i*Gauss)/Spesi
				M_star2=dot_product(w,M_i*Gauss)/Spesi
				g_star2=dot_product(w,g_i*Gauss)/Spesi
				rho_star2=dot_product(w,rho_i*Gauss)/Spesi
				!!! 
				if (photIsocAvail) then 
					BmV_star2=dot_product(w,BmV_i*Gauss)/Spesi
					BC_star2=dot_product(w,BC_i*Gauss)/Spesi
				end if 
				!!! 
				I_t_star2=sqrt(dot_product(w*Gauss,(t_i-t_star2)**2)/Spesi)
				I_Teff_star2=sqrt(dot_product(w*Gauss,(Teff_i-Teff_star2)**2)/Spesi)
				I_L_star2=sqrt(dot_product(w*Gauss,(L_i-L_star2)**2)/Spesi)
				I_M_star2=sqrt(dot_product(w*Gauss,(M_i-M_star2)**2)/Spesi)
				I_g_star2=sqrt(dot_product(w*Gauss,(g_i-g_star2)**2)/Spesi)
				I_rho_star2=sqrt(dot_product(w*Gauss,(rho_i-rho_star2)**2)/Spesi)
				!!! 
				if (photIsocAvail) then 
					I_BmV_star2=sqrt(dot_product(w*Gauss,(BmV_i-BmV_star2)**2)/Spesi)
					I_BC_star2=sqrt(dot_product(w*Gauss,(BC_i-BC_star2)**2)/Spesi)
				end if 

				logg_star2=log10(g_star2)
				R_star2=sqrt(L_star2/(Teff_star2/TeffSun)**4)
				I_logg_star2=log10(e)*I_g_star2/g_star2
				I_R_star2=R_star2*(0.5*I_L_star2/L_star2+2.*I_Teff_star2/Teff_star2)

				tol=0.02D0
				kwhile=1
				cycleW=0
				do while (abs(M_star2-M_star)>I_M_star .and. abs(R_star2-R_star)>I_R_star .and. kwhile<10 &
						& .and. M_star2>=MinfED .and. M_star2<=MsupED) 
					t_star=t_star2
					Teff_star=Teff_star2
					L_star=L_star2
					M_star=M_star2
					g_star=g_star2
					rho_star=rho_star2
					if (photIsocAvail) then 
						BmV_star=BmV_star2
						BC_star=BC_star2
					end if 
					I_t_star=I_t_star2
					I_Teff_star=I_Teff_star2
					I_L_star=I_L_star2
					I_M_star=I_M_star2
					I_g_star=I_g_star2
					I_rho_star=I_rho_star2
					if (photIsocAvail) then 
						I_BmV_star=I_BmV_star2
						I_BC_star=I_BC_star2
					end if 
					logg_star=logg_star2
					R_star=R_star2
					I_logg_star=I_logg_star2
					I_R_star=I_R_star2
					if (allocated(Zguess)) then
						deallocate(Zguess)
					end if
					Mndx=minloc(abs(M_star-MeAv),1)
					allocate(Ztevo(Endxf(Mndx)-Endxi(Mndx)+1,size(ZtevoTab,2)))
					Ztevo=ZtevoTab(Endxi(Mndx):Endxf(Mndx),:)
					!!Polynomials interpolating the Z=Z(t) have degree=3
					allocate(Zini(size(Ztevo,1))); allocate(tsp(size(Ztevo,1))); allocate(tmax(size(Ztevo,1)))
					allocate(a1(size(Ztevo),4)); allocate(a2(size(Ztevo),4))
					Zini=Ztevo(:,1)
					tsp=Ztevo(:,2)
					tmax=Ztevo(:,3)
					a1=Ztevo(:,4:7)
					a2=Ztevo(:,8:11)
					t_starvec=(/1.D0,t_star,t_star**2,t_star**3/)
					!
					allocate(Zguess(size(Ztevo,1)))
					!
					do iZg=1,size(Ztevo,1) 
						if (t_star>tmax(iZg)) then 
							Zguess(iZg)=Zini(iZg)
						else if (tsp(iZg)==-1) then 
							Zguess(iZg)=dot_product(a1(iZg,:),t_starvec)
						else
							if (t_star<=tsp(iZg)) then 
								Zguess(iZg)=dot_product(a1(iZg,:),t_starvec)
							else
								Zguess(iZg)=dot_product(a2(iZg,:),t_starvec)
							end if 
						end if 
					end do
					deallocate(Ztevo)
					allocate(indxZg(size(Zguess))) 
					call indexx(abs(Zguess-Z_),indxZg) 
					Zg1=Zguess(indxZg(1))
					Zg2=Zguess(indxZg(2))
					rZ=(Z_-Zg1)/(Zg2-Zg1)
					Zgini1=Zini(indxZg(1))
					Zgini2=Zini(indxZg(2))
					Zgini=Zgini1+rZ*(Zgini2-Zgini1)
					indxZgini=minloc(abs(Zgini-Zvec),1) 
					Z2=Zvec(indxZgini)
					
					deallocate(Zini); deallocate(tsp); deallocate(tmax)
					deallocate(a1); deallocate(a2); deallocate(indxZg)
					
					deallocate(Isoc)
					deallocate(t_iso); deallocate(logt_iso)
					deallocate(M_iso); deallocate(logL_iso)
					deallocate(L_iso); deallocate(logTeff_iso)
					deallocate(Teff_iso); deallocate(logg_iso)
					deallocate(g_iso); deallocate(R_iso)
					deallocate(rho_iso); deallocate(logrho_iso)
					if (allocated(BmV_iso)) then
						deallocate(BmV_iso); deallocate(V_iso)
						deallocate(BC_iso)
					end if
					deallocate(t_i); deallocate(dist)
					deallocate(logTeff_i); deallocate(logL_i); deallocate(M_i); deallocate(logg_i)
					deallocate(logrho_i); deallocate(Teff_i); deallocate(g_i); deallocate(L_i)
					deallocate(rho_i); deallocate(Zguess)
					if (allocated(y_i)) then
						deallocate(y_i)
					end if
					deallocate(Gauss); deallocate(w)
					if (photIsocAvail) then
						deallocate(BC_i); deallocate(BmV_i)
					end if
					
					deallocate(MTrAv)
					allocate(MTrAv(nM(indxZt(indxZgini))))
					MTrAv=Mav(indxZt(indxZgini),1:nM(indxZt(indxZgini)))
					Z=Z2
					FeH=log10(Z)+costZ
					Yini=dnint((He0+dYdZ*Z)*100)/100
					
					deallocate(TabVel)
					xVTi=Vndxi(indxZt(indxZgini))
					xVTf=Vndxf(indxZt(indxZgini))
					allocate(TabVel(xVTf-xVTi+1,size(velTrackTab,2)))
					TabVel=velTrackTab(xVTi:xVTf,:)
					
					deallocate(ZAMSz)
					xZAi=ZAndxi(indxZt(indxZgini))
					xZAf=ZAndxf(indxZt(indxZgini))
					allocate(ZAMSz(xZAf-xZAi+1,size(ZAMStab,2)))
					ZAMSz=ZAMStab(xZAi:xZAf,:)
					
					Zi=Zndxi(indxZgini)
					Zf=Zndxf(indxZgini)
					allocate(Isoc(Zf-Zi+1,size(IsocTab,2)))
					Isoc=IsocTab(Zi:Zf,:)
					row_i=size(Isoc,1)
					allocate(t_iso(row_i)); allocate(logt_iso(row_i))
				
					t_iso=Isoc(:,ct)
					logt_iso=dnint(log10(t_iso)*100)/100
					t_iso=10.**logt_iso
					first_logt=logt_iso(1)
					last_logt=logt_iso(row_i)

					allocate(M_iso(row_i)); allocate(logL_iso(row_i))
					allocate(L_iso(row_i)); allocate(logTeff_iso(row_i))
					allocate(Teff_iso(row_i)); allocate(logg_iso(row_i))
					allocate(g_iso(row_i)); allocate(R_iso(row_i))
					allocate(rho_iso(row_i)); allocate(logrho_iso(row_i))
		
					M_iso=Isoc(:,cM)
					logL_iso=Isoc(:,clogL)
					L_iso=10.**(logL_iso)
					logTeff_iso=Isoc(:,clogTe)
					Teff_iso=10.**Isoc(:,clogTe)
					logg_iso=Isoc(:,clogg)
					g_iso=10.**Isoc(:,clogg)
					R_iso=sqrt(L_iso/(Teff_iso/TeffSun)**4) !raggi solari
					rho_iso=M_iso/(R_iso**3)*rhoSun !g/cm3
					logrho_iso=log10(rho_iso)
					if (photIsocAvail) then
						allocate(V_iso(row_i)); allocate(BmV_iso(row_i))
						allocate(BC_iso(row_i))
						V_iso=Isoc(:,cmag0)
						if (isEq(useColor,1.D0,2)) then
							BmV_iso=Isoc(:,cmag2)-Isoc(:,cmag1)
						else
							BmV_iso=logTeff_iso
						end if
						!!! 
						BC_iso=Isoc(:,cmbol)-Isoc(:,cmag0)
						!!! 
					end if 

					!!!Selection of reference isochrones!!! 
					!!First evaluate which isochrones are the closest to the star to be analysed
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
					!!!!!!!NEW CODE 12/2/14!!!!!!!!!!!!!!! 
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
					allocate(dist_ndx(size(logTeff_iso))) !equal to any other size of _iso
					allocate(x_iso(size(logTeff_iso))); allocate(xx_iso(size(logTeff_iso)))
					if (calibHRD) then
						call indexx(sqrt((BmV-BmV_iso)**2+(Vass-V_iso)**2), dist_ndx) 
						x=BmV
						x_iso=BmV_iso
						xx_iso=logTeff_iso
						if (isEq(useColor,1.D0,2)) then
							if (idCol.eq.1) then
								cBrif=cBBmV
							else
								cBrif=cBTeff
							end if
						else
							cBrif=cBTeff
						end if
					end if 
					if (calibNoD) then
						if (.not.isEq(Rf,-1.D0,2)) then
							call indexx(sqrt((BmV-BmV_iso)**2+(logL-logL_iso)**2),dist_ndx)
						else if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
							& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))) then 
							call indexx(sqrt((BmV-BmV_iso)**2+(logg-logg_iso)**2),dist_ndx) 
						else !rhoAvail
							call indexx(sqrt((BmV-BmV_iso)**2+(logrho-logrho_iso)**2),dist_ndx) 
						end if 
						x=BmV
						x_iso=BmV_iso
						xx_iso=logTeff_iso
						cBrif=cBBmV
					end if 
					if (calibSPEC) then
						if (.not.isEq(Rf,-1.D0,2)) then
							call indexx(sqrt((logTeff-logTeff_iso)**2+(logL-logL_iso)**2),dist_ndx)
						else if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
							& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))) then 
							call indexx(sqrt((logTeff-logTeff_iso)**2+(logg-logg_iso)**2), dist_ndx) 
						else !rhoAvail
							call indexx(sqrt((logTeff-logTeff_iso)**2+(logrho-logrho_iso)**2), dist_ndx) 
						end if 
						x=logTeff
						xtau=x
						x_iso=logTeff_iso
						xx_iso=BmV_iso
						cBrif=cBTeff
					end if
					allocate(ageMinDist(size(logt_iso)))
					ageMinDist=logt_iso(dist_ndx)
					call uniqueFast(ageMinDist,2,ix,.true.)
					allocate(dist_ndxU(size(ix))) 
					dist_ndxU=dist_ndx(ix)
		
					deallocate(dist_ndx); deallocate(ageMinDist)
					deallocate(ix)
					rowAge=1
					! allocazione delle variabili già incontrate
					allocate(t_iTmp(size(dist_ndxU)));allocate(distTmp(size(dist_ndxU),2))
					allocate(logTeff_iTmp(size(dist_ndxU)));allocate(logL_iTmp(size(dist_ndxU)))
					allocate(M_iTmp(size(dist_ndxU)));allocate(logg_iTmp(size(dist_ndxU)))
					allocate(logrho_iTmp(size(dist_ndxU)))
					if (photIsocAvail) then 
						allocate(BmV_iTmp(size(dist_ndxU)));allocate(BC_iTmp(size(dist_ndxU)))
					end if
					!
					do Nage=1,size(dist_ndxU) 
						ndxDU=dist_ndxU(Nage)
						if (photIsocAvail) then 
							BmV_i1=BmV_iso(ndxDU)
							V_i1=V_iso(ndxDU)
						end if 
						logTeff_i1=logTeff_iso(ndxDU)
						logL_i1=logL_iso(ndxDU)
						M_i1=M_iso(ndxDU)
						logg_i1=logg_iso(ndxDU)
						logrho_i1=logrho_iso(ndxDU)
						
						call choose_i2Col(ndxDU,x,x_iso,xx_iso,calibSPEC,hstar,hstarlim,V_iso,logL_iso, &
							& M_iso,logg_iso,logrho_iso,logt_iso,BmV_i2,V_i2,logTeff_i2,logL_i2,M_i2, &
							& logg_i2,logrho_i2)
						if (calibSPEC) then
							x_i1=logTeff_i1
							x_i2=logTeff_i2
						else
							x_i1=BmV_i1
							x_i2=BmV_i2
						end if
						if (isEq(x_i1,x_i2,4)) then 
							cycle
						end if 
						t_iTmp(rowAge)=t_iso(ndxDU)
						if (photIsocAvail) then 
							BC_iTmp(rowAge)=Isoc(ndxDU,cmbol)-Isoc(ndxDU,cmag0) !just rough.  When I recover 
							! BC in the logg-BmV plane, I'll make the interpolation using BC_i1 e BC_i2
						end if 
						if (calibHRD) then !!Condiz sulla calibrazione da fare 15/12/14
							call findX_i(BmV,Vass,BmV_i1,BmV_i2,V_i1,V_i2,x_i,dist0)
							BmV_iTmp(rowAge)=x_i
						end if 
						if (calibNoD) then
							if (.not.(isEq(Rf,-1.D0,2))) then
								call findX_i(BmV,logL,BmV_i1,BmV_i2,logL_i1,logL_i2,x_i,dist0)
							else if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
								& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))) then 
								call findX_i(BmV,logg,BmV_i1,BmV_i2,logg_i1,logg_i2,x_i,dist0)
							else !rhoAvail
								call findX_i(BmV,logrho,BmV_i1,BmV_i2,logrho_i1,logrho_i2,x_i,dist0)
							end if 
							BmV_iTmp(rowAge)=x_i
						end if 
						if (calibSPEC) then
							if (.not.(isEq(Rf,-1.D0,2))) then
								call findX_i(logTeff,logL,logTeff_i1,logTeff_i2,logL_i1,logL_i2,x_i,dist0)
							else if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
								& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))) then 
								call findX_i(logTeff,logg,logTeff_i1,logTeff_i2,logg_i1,logg_i2,x_i,dist0)
							else !rhoAvail
								call findX_i(logTeff,logrho,logTeff_i1,logTeff_i2,logrho_i1,logrho_i2,x_i,dist0)
							end if 
							logTeff_iTmp(rowAge)=x_i 
						end if
			
						distTmp(rowAge,:)=(/ dist0, logt_iso(ndxDU) /)
		
						!!18/3/15 The following if cancel the criterion of perpendicular distance
						!!between a star and the isochrone if the theoretical point to be interpolated
						!!doesn't fall inside the isochrone (i.e. between x_i1 and x_i2), but on its
						!!extension. This possibility, in fact, would build fictitious isochrones that
						!!could be erroneously close to the star
						if (.not.((x_i>=x_i1.and.x_i<=x_i2).or.(x_i>=x_i2.and.x_i<=x_i1))) then 
							if (photIsocAvail) then 
								BmV_iTmp(rowAge)=BmV_i1
							end if 
							logg_iTmp(rowAge)=logg_i1
							logTeff_iTmp(rowAge)=logTeff_i1
							logL_iTmp(rowAge)=logL_i1
							M_iTmp(rowAge)=M_i1
							logrho_iTmp(rowAge)=logrho_i1
						else
							if (calibSPEC) then 
								if (photIsocAvail) then 
									mBV=(BmV_i2-BmV_i1)/(x_i2-x_i1)
									qBV=-mBV*x_i1+BmV_i1
									BmV_iTmp(rowAge)=mBV*x_i+qBV
								end if 
							else
								mT=(logTeff_i2-logTeff_i1)/(x_i2-x_i1)
								qT=-mT*x_i1+logTeff_i1
								logTeff_iTmp(rowAge)=mT*x_i+qT
							end if 
							mg=(logg_i2-logg_i1)/(x_i2-x_i1)
							qg=-mg*x_i1+logg_i1
							logg_iTmp(rowAge)=mg*x_i+qg
							mL=(logL_i2-logL_i1)/(x_i2-x_i1)
							qL=-mL*x_i1+logL_i1
							logL_iTmp(rowAge)=mL*x_i+qL
							mM=(M_i2-M_i1)/(x_i2-x_i1)
							qM=-mM*x_i1+M_i1
							M_iTmp(rowAge)=mM*x_i+qM
							mrh=(logrho_i2-logrho_i1)/(x_i2-x_i1)
							qrho=-mrh*x_i1+logrho_i1
							logrho_iTmp(rowAge)=mrh*x_i+qrho
						end if 

						rowAge=rowAge+1
					end do
					deallocate(x_iso); deallocate(xx_iso)
		
					allocate(t_i(rowAge-1));allocate(dist(rowAge-1,2))
					allocate(logTeff_i(rowAge-1));allocate(logL_i(rowAge-1))
					allocate(M_i(rowAge-1));allocate(logg_i(rowAge-1))
					allocate(logrho_i(rowAge-1))
					allocate(BmV_i(rowAge-1));allocate(BC_i(rowAge-1))
					t_i=t_iTmp(1:rowAge-1)
					dist=distTmp(1:rowAge-1,:)
					logTeff_i=logTeff_iTmp(1:rowAge-1)
					logL_i=logL_iTmp(1:rowAge-1)
					M_i=M_iTmp(1:rowAge-1)
					logg_i=logg_iTmp(1:rowAge-1)
					logrho_i=logrho_iTmp(1:rowAge-1)
					BmV_i=BmV_iTmp(1:rowAge-1)
					BC_i=BC_iTmp(1:rowAge-1)
					deallocate(t_iTmp);deallocate(distTmp);deallocate(logTeff_iTmp)
					deallocate(logL_iTmp);deallocate(M_iTmp);deallocate(logg_iTmp)
					deallocate(logrho_iTmp)
					if (photIsocAvail) then
						deallocate(BmV_iTmp);deallocate(BC_iTmp)
					end if
		
					allocate(Teff_i(rowAge-1));allocate(L_i(rowAge-1))
					allocate(g_i(rowAge-1));allocate(rho_i(rowAge-1))
		
					Teff_i=10.**logTeff_i
					L_i=10.**logL_i
					g_i=10.**logg_i
					rho_i=10.**logrho_i
		
					deallocate(dist_ndxU)
		
					allocate(dsorted(size(dist)))
					dsorted=dist(:,1) !then it will be overwritten by sort so that it will be actually sorted
					call sort1(dsorted)
					allocate(ix_sorted(size(dist,1)))
					call indexx(dist(:,1), ix_sorted) 
					!!Avoid extremely young isochrones (t<10 Myr) in the following calibration process
					age_s1=1 
					age_iso1=dist(ix_sorted(age_s1),2)
					do while (age_iso1.lt.logtlimCal)
						age_s1=age_s1+1
						age_iso1=dist(ix_sorted(age_s1),2)
					end do
					age_s=age_s1+1
					if (dsorted(age_s1)<1.e-6) then 
						do while ((log10(dsorted(age_s))-log10(dsorted(age_s1))<4.5.or.dist(ix_sorted(age_s),2).lt.logtlimCal) &
								& .and. age_s<size(ix_sorted))
							age_s=age_s+1
						end do 
					else
						do while (dist(ix_sorted(age_s),2).lt.logtlimCal .and. age_s<size(ix_sorted))
							age_s=age_s+1
						end do
					end if
					deallocate(dsorted)
					age_iso2=dist(ix_sorted(age_s),2)
					if (age_iso1>age_iso2) then 
						tmp_age1=age_iso1
						age_iso1=age_iso2
						age_iso2=tmp_age1 !age_iso2 BECOMES the isoc with the minimum distance from the star
						agePivot=age_iso2
					else
						agePivot=age_iso1 !age_iso1 REMAINS the isoc with the minimum distance from the star
					end if 

					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
					!!!!!!!!!!!END NEW CODE 
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
		
					call selectIsoc(logt_iso,age_iso1,logtstep,last_logt,ti0,tf0)
					ti(1)=ti0
					tf(1)=tf0
					call selectIsoc(logt_iso,age_iso2,logtstep,last_logt,ti0,tf0)
					ti(2)=ti0
					tf(2)=tf0
							
					if (calibHRD) then !!15/12/14
						VMag=V-DM0
						deallocate(VMag_iso)
						allocate(VMag_iso(size(Isoc,1)))
						VMag_iso=Isoc(:,cmag0)
						BmVuguali=.false.
						do jd=1,2 
							allocate(VMag_isor(tf(jd)-ti(jd)+1)); allocate(logL_isor(tf(jd)-ti(jd)+1))
							allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(Diff_BmV(tf(jd)-ti(jd)+1))
							VMag_isor=VMag_iso(ti(jd):tf(jd))
							logL_isor=logL_iso(ti(jd):tf(jd))
							BmV_isor=BmV_iso(ti(jd):tf(jd))
							Diff_BmV=sqrt((BmV-BmV_isor)**2+(VMag-VMag_isor)**2) !!!!!!!!!!!!18/2
							if (jd==1) then 
								id(jd)=minloc(Diff_BmV,1) 
								BmV_is1v(jd)=BmV_isor(id(jd))
							else
								allocate(ndxDv(size(Diff_BmV)))
								call indexx(Diff_BmV,ndxDv) 
								rowDv=size(ndxDv)
								allocate(VMagD(size(VMag_isor)))
								VMagD=VMag_isor(ndxDv)
								kD=1
								id(jd)=ndxDv(kD)
								BmV_is1v(jd)=BmV_isor(id(jd))
								do while ((abs(VMag_is1v(1)-VMagD(kD))<=0.01 .or. abs(VMag_is2v(1)-VMagD(kD))<=0.01) &
										& .and. kD<rowDv .and. VMag<2)!!!0.01 invece di 0.1!18/2 
									kD=kD+1
									id(jd)=ndxDv(kD)
									BmV_is1v(jd)=BmV_isor(id(jd))
								end do
								deallocate(ndxDv); deallocate(VMagD)
							end if 
							VMag_is1v(jd)=VMag_isor(id(jd))
							logL_is1v(jd)=logL_isor(id(jd))
							!!!
							call choose_i2xy1y2(BmV,BmV_isor,VMag_isor,logL_isor,id(jd),hstar,hstarlim, & !VMag,
											& BmV_is1v(jd),x_is2,VMag_is1v(jd),y1_is2,logL_is1v(jd),y2_is2)
							BmV_is2v(jd)=x_is2
							VMag_is2v(jd)=y1_is2
							logL_is2v(jd)=y2_is2
							!!!
				
							if (isEq(BmV_is1v(jd),BmV_is2v(jd),4)) then 	!!!if 10/2/14 -modified 7/11/2018
								BmVuguali=.true.
								VMag_isv(jd)=VMag
								call InterpLin_M(reshape((/VMag_is1v(jd),VMag_is2v(jd),logL_is1v(jd),logL_is2v(jd)/), &
									& (/2,2/)),VMag,1,1,(/2/),logL_isv(jd),xlow,ylow,xup,yup)
							else
								VMag_isv(jd)=(VMag_is2v(jd)-VMag_is1v(jd))/(BmV_is2v(jd)-BmV_is1v(jd))* &
										& (BmV-BmV_is1v(jd))+VMag_is1v(jd)
								logL_isv(jd)=(logL_is2v(jd)-logL_is1v(jd))/(BmV_is2v(jd)-BmV_is1v(jd))* &
										& (BmV-BmV_is1v(jd))+logL_is1v(jd)	
							end if
				
							deallocate(VMag_isor); deallocate(logL_isor)
							deallocate(BmV_isor); deallocate(Diff_BmV)
						end do
						 
						dCV_isv=abs(VMag_isv(1)-VMag_isv(2))  !vertical dist between isoc in the Color-Mag diagr
						dCL_isv=abs(logL_isv(1)-logL_isv(2))  !vertical dist between isoc in the Color-logL diagr
						if ((dCV_isv>dlim .and. dCL_isv>dlim)) then  !condition added 17/2/14
			!				dCVv=abs(VMag-VMag_isv)  		  !vertical dist between the star and the 2 isoc
							if (VMag<VMag_isv(1) .and. VMag<VMag_isv(2)) then 
								statev=1
							else if (VMag>VMag_isv(1) .and. VMag>VMag_isv(2)) then 
								statev=2
							else
								statev=3
							end if
							call calibrateDiag(statev,dCVvlim,VMag_isv,logL_isv,VMag,logL) 
							 
!!!							if (isnan(logL)) then 
!!!								cycle
!!!							end if
				
							L=10.**logL
							BC=-2.5*log10(L/d**2)-V-0.23
							I_L=L*(0.4*log(10.)*0.03+2.*I_d/d)
							I_logL=I_L/L*log10(e)
							logLuguali=.false.
							if (isEq(useColor,1.D0,2)) then
								do jd=1,2
									allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(logT_isor(tf(jd)-ti(jd)+1))
									allocate(logL_isor(tf(jd)-ti(jd)+1)); allocate(Diff_logL(tf(jd)-ti(jd)+1))
						
									BmV_isor=BmV_iso(ti(jd):tf(jd))
									logT_isor=logTeff_iso(ti(jd):tf(jd))
									logL_isor=logL_iso(ti(jd):tf(jd))
									Diff_logL=sqrt((logL-logL_isor)**2+(BmV-BmV_isor)**2) !!!!!!!!!!!!!!!!18/2/14
									if (jd==1) then
										id(jd)=minloc(Diff_logL,1) 
										logL_is1(jd)=logL_isor(id(jd))
									else
										allocate(ndxD(size(Diff_logL)))
										call indexx(Diff_logL,ndxD) 
										rowD=size(ndxD)
										allocate(BmVD(size(BmV_isor)))
										BmVD=BmV_isor(ndxD)
										kD=1
										id(jd)=ndxD(kD)
										logL_is1(jd)=logL_isor(id(jd))
										do while ((abs(BmV_is1(1)-BmVD(kD))<=0.04 .or. abs(BmV_is2(1)-BmVD(kD))<=0.04) &
												& .and. kD<rowD .and. logL>1.5) 
											kD=kD+1
											id(jd)=ndxD(kD)
											logL_is1(jd)=logL_isor(id(jd))
										end do
										deallocate(ndxD); deallocate(BmVD)
									end if 
									BmV_is1(jd)=BmV_isor(id(jd))
									logT_is1(jd)=logT_isor(id(jd))
									!!!
									call choose_i2xy1y2(logL,logL_isor,BmV_isor,logT_isor,id(jd),hstar,hstarlim, & !BmV,
										& logL_is1(jd),x_is2,BmV_is1(jd),y1_is2,logT_is1(jd),y2_is2)
									logL_is2(jd)=x_is2
									BmV_is2(jd)=y1_is2
									logT_is2(jd)=y2_is2
									!!!
									if (isEq(logL_is1(jd),logL_is2(jd),4)) then 
										logLuguali=.true.
										BmV_is(jd)=BmV
										call InterpLin_M(reshape((/BmV_is1(jd),BmV_is2(jd),logT_is1(jd),logT_is2(jd)/), &
											& (/2,2/)),BmV,1,1,(/2/),logT_is(jd),xlow,ylow,xup,yup)
									else
										BmV_is(jd)=(logL-logL_is1(jd))*(BmV_is2(jd)-BmV_is1(jd))/ &
													& (logL_is2(jd)-logL_is1(jd))+BmV_is1(jd)
										logT_is(jd)=(logL-logL_is1(jd))*(logT_is2(jd)-logT_is1(jd))/ &
													& (logL_is2(jd)-logL_is1(jd))+logT_is1(jd)
									end if 
														
									deallocate(BmV_isor); deallocate(logT_isor)
									deallocate(logL_isor); deallocate(Diff_logL)
								end do 
								
								dCL_is=abs(BmV_is(1)-BmV_is(2))		!dist between isoc in the 'CMD'
								dHR_is=abs(logT_is(1)-logT_is(2))	!dist between isoc in the HRD
								if (dCL_is>dlim .and. dHR_is>dlim) then !condition added 17/2/14
				!					dCL=abs(BmV-BmV_is)				!vector of distances between the star and the 2 isoc
									if (BmV<BmV_is(1) .and. BmV<BmV_is(2)) then 
										state=1
									else if (BmV>BmV_is(1) .and. BmV>BmV_is(2)) then 
										state=2
									else
										state=3
									end if
									call calibrateDiag(state,dCLlim,BmV_is,logT_is,BmV,logTeff)
									 
!!!									if (isnan(logTeff)) then 
!!!										cycle
!!!									end if 
								else if (age_s==size(ix_sorted)) then 
									cycleW=1
									cycle
								end if
							else
								logTeff=BmV
								dCL_is=2.*dlim
								dHR_is=2.*dlim
							end if
							if (idCol.eq.1) then
								xtau=BmV
							else
								xtau=logTeff
							end if
						else if (age_s==size(ix_sorted)) then 
							cycleW=1
							cycle
						else
							dCL_is=0.
							dHR_is=0.
						end if 
						do while ((dCV_isv<=dlim .or. dCL_isv<=dlim .or. dCL_is<=dlim .or. dHR_is<=dlim) .and. &
								& age_s<size(ix_sorted)) 
							age_s=age_s+1
							do while (dist(ix_sorted(age_s),2).lt.logtlimCal .and. age_s<size(ix_sorted))
								age_s=age_s+1
							end do
							if (isEq(age_iso1,agePivot,2)) then 
								age_iso2=dist(ix_sorted(age_s),2)
								if (age_iso1>age_iso2) then 
									tmp_age1=age_iso1
									age_iso1=age_iso2
									age_iso2=tmp_age1
									agePivot=age_iso2
								end if 
							else
								age_iso1=dist(ix_sorted(age_s),2)
								if (age_iso1>age_iso2) then 
									tmp_age1=age_iso1
									age_iso1=age_iso2
									age_iso2=tmp_age1
									agePivot=age_iso1
								end if 
							end if 
				
							call selectIsoc(logt_iso,age_iso1,logtstep,last_logt,ti0,tf0)
							ti(1)=ti0
							tf(1)=tf0
							call selectIsoc(logt_iso,age_iso2,logtstep,last_logt,ti0,tf0)
							ti(2)=ti0
							tf(2)=tf0
				
							VMag=V-DM0
							VMag_iso=Isoc(:,cmag0)
							BmVuguali=.false.
							do jd=1,2
								allocate(VMag_isor(tf(jd)-ti(jd)+1)); allocate(logL_isor(tf(jd)-ti(jd)+1))
								allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(Diff_BmV(tf(jd)-ti(jd)+1))					
								VMag_isor=VMag_iso(ti(jd):tf(jd))
								logL_isor=logL_iso(ti(jd):tf(jd))
								BmV_isor=BmV_iso(ti(jd):tf(jd))
								Diff_BmV=sqrt((BmV-BmV_isor)**2+(VMag-VMag_isor)**2) !!!!!!!!!!!!!!!!!18/2/14
								if (jd==1) then 
									id(jd)=minloc(Diff_BmV,1) 
									BmV_is1v(jd)=BmV_isor(id(jd))
								else
									allocate(ndxDv(size(Diff_BmV)))
									call indexx(Diff_BmV,ndxDv) 
									rowDv=size(ndxDv)
									allocate(VMagD(size(VMag_isor)))
									VMagD=VMag_isor(ndxDv)
									kD=1
									id(jd)=ndxDv(kD)
									BmV_is1v(jd)=BmV_isor(id(jd))
									do while ((abs(VMag_is1v(1)-VMagD(kD))<=0.01 .or. abs(VMag_is2v(1)-VMagD(kD)) &
											& <=0.01) .and. kD<rowDv .and. VMag<2)!0.01 invece di 0.1 18/2 
										kD=kD+1
										id(jd)=ndxDv(kD)
										BmV_is1v(jd)=BmV_isor(id(jd))
									end do
									deallocate(ndxDv); deallocate(VMagD)
								end if 
								VMag_is1v(jd)=VMag_isor(id(jd))
								logL_is1v(jd)=logL_isor(id(jd))
								!!!
								call choose_i2xy1y2(BmV,BmV_isor,VMag_isor,logL_isor,id(jd),hstar,hstarlim, & !VMag,
									& BmV_is1v(jd),x_is2,VMag_is1v(jd),y1_is2,logL_is1v(jd),y2_is2) !!!
								BmV_is2v(jd)=x_is2
								VMag_is2v(jd)=y1_is2
								logL_is2v(jd)=y2_is2
								!!!
								if (isEq(BmV_is1v(jd),BmV_is2v(jd),4)) then 	!!!if 10/2/14 -modified 7/11/2018
									BmVuguali=.true.
									VMag_isv(jd)=VMag
									call InterpLin_M(reshape((/VMag_is1v(jd),VMag_is2v(jd),logL_is1v(jd),logL_is2v(jd)/), &
										& (/2,2/)),VMag,1,1,(/2/),logL_isv(jd),xlow,ylow,xup,yup)
								else
									VMag_isv(jd)=(VMag_is2v(jd)-VMag_is1v(jd))/(BmV_is2v(jd)-BmV_is1v(jd))* &
											& (BmV-BmV_is1v(jd))+VMag_is1v(jd)
									logL_isv(jd)=(logL_is2v(jd)-logL_is1v(jd))/(BmV_is2v(jd)-BmV_is1v(jd))* &
											& (BmV-BmV_is1v(jd))+logL_is1v(jd)	
								end if
					
								deallocate(VMag_isor); deallocate(logL_isor)
								deallocate(BmV_isor); deallocate(Diff_BmV)
							end do 

							dCV_isv=abs(VMag_isv(1)-VMag_isv(2))	!vertical dist between isochrones in the CMD
							dCL_isv=abs(logL_isv(1)-logL_isv(2))	!vert dist between isoc in the Color-logL diagr
							if ((dCV_isv>dlim .and. dCL_isv>dlim)) then !condition added 17/2/14
			!					dCVv=abs(VMag-VMag_isv)				!vertical dist between the star and the 2 isoc
								if (VMag<VMag_isv(1) .and. VMag<VMag_isv(2)) then 
									statev=1
								else if (VMag>VMag_isv(1) .and. VMag>VMag_isv(2)) then 
									statev=2
								else
									statev=3
								end if
								call calibrateDiag(statev,dCVvlim,VMag_isv,logL_isv,VMag,logL)
								 
!!!								if (isnan(logL)) then 
!!!									cycle
!!!								end if 
								L=10.**logL
								BC=-2.5*log10(L/d**2)-V-0.23
								I_L=L*(0.4*log(10.)*0.03+2.*I_d/d)
								I_logL=I_L/L*log10(e)
								logLuguali=.false.
								if (isEq(useColor,1.D0,2)) then
									do jd=1,2
										allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(logT_isor(tf(jd)-ti(jd)+1))
										allocate(logL_isor(tf(jd)-ti(jd)+1)); allocate(Diff_logL(tf(jd)-ti(jd)+1))
										logL_isor=logL_iso(ti(jd):tf(jd))
										logT_isor=logTeff_iso(ti(jd):tf(jd))
										BmV_isor=BmV_iso(ti(jd):tf(jd))
										Diff_logL=sqrt((logL-logL_isor)**2+(BmV-BmV_isor)**2) !!!!!!!!!!!!!!!18/2
										if (jd==1) then 
											id(jd)=minloc(Diff_logL,1) 
											logL_is1(jd)=logL_isor(id(jd))
										else
											allocate(ndxD(size(Diff_logL)))
											call indexx(Diff_logL,ndxD) 
											rowD=size(ndxD)
											allocate(BmVD(size(BmV_isor)))
											BmVD=BmV_isor(ndxD)
											kD=1
											id(jd)=ndxD(kD)
											logL_is1(jd)=logL_isor(id(jd))
											do while ((abs(BmV_is1(1)-BmVD(kD))<=0.04 .or. abs(BmV_is2(1)-BmVD(kD)) &
													& <=0.04) .and. kD<rowD .and. logL>1.5) 
												kD=kD+1
												id(jd)=ndxD(kD)
												logL_is1(jd)=logL_isor(id(jd))
											end do
											deallocate(ndxD); deallocate(BmVD)
										end if 
										BmV_is1(jd)=BmV_isor(id(jd))
										logT_is1(jd)=logT_isor(id(jd))
										!!!
										call choose_i2xy1y2(logL,logL_isor,BmV_isor,logT_isor,id(jd),hstar,hstarlim, & !BmV,
											& logL_is1(jd),x_is2,BmV_is1(jd),y1_is2,logT_is1(jd),y2_is2)
										logL_is2(jd)=x_is2
										BmV_is2(jd)=y1_is2
										logT_is2(jd)=y2_is2
										!!!
										if (isEq(logL_is1(jd),logL_is2(jd),4)) then 
											logLuguali=.true.
											BmV_is(jd)=BmV
											call InterpLin_M(reshape((/BmV_is1(jd),BmV_is2(jd),logT_is1(jd),logT_is2(jd)/), &
												& (/2,2/)),BmV,1,1,(/2/),logT_is(jd),xlow,ylow,xup,yup)
										else
											BmV_is(jd)=(logL-logL_is1(jd))*(BmV_is2(jd)-BmV_is1(jd))/ &
														& (logL_is2(jd)-logL_is1(jd))+BmV_is1(jd)
											logT_is(jd)=(logL-logL_is1(jd))*(logT_is2(jd)-logT_is1(jd))/ &
														& (logL_is2(jd)-logL_is1(jd))+logT_is1(jd)
										end if 
															
										deallocate(BmV_isor); deallocate(logT_isor)
										deallocate(logL_isor); deallocate(Diff_logL)
									end do 
									
									dCL_is=abs(BmV_is(1)-BmV_is(2))		!dist between isoc in the 'CMD'
									dHR_is=abs(logT_is(1)-logT_is(2))	!dist between isoc in the HRD
									if ((dCL_is>dlim .and. dHR_is>dlim)) then !condition added 17/2/14
				!						dCL=abs(BmV-BmV_is)				!vector of the distances between star-->2 isoc
										if (BmV<BmV_is(1) .and. BmV<BmV_is(2)) then 
											state=1
										else if (BmV>BmV_is(1) .and. BmV>BmV_is(2)) then 
											state=2
										else
											state=3
										end if
							
										call calibrateDiag(state,dCLlim,BmV_is,logT_is,BmV,logTeff)
							
!!!										if (isnan(logTeff)) then 
!!!											cycle
!!!										end if 
									else if (age_s==size(ix_sorted)) then 
										cycleW=1
										cycle
									end if
								else
									logTeff=BmV
									dCL_is=2.*dlim
									dHR_is=2.*dlim
								end if
								if (idCol.eq.1) then
									xtau=BmV
								else
									xtau=logTeff
								end if
							else if (age_s==size(ix_sorted)) then 
								cycleW=1
								cycle
							end if 
						end do
						if (cycleW.eq.1) then
							cycle
						end if

						Teff=10.**logTeff
						if (.not.calibSPEC .and. Teffinput.ne.-1) then !Teffinput condition added 14/7/2015
							if (abs(Teff-Teffinput)>300) then !Input B-V totally inconsistent with spectroscopic Teff
								cycle
							end if 
						end if
						if (isEq(useColor,1.D0,2)) then 
							I_Teff=0.01*Teff
							I_logTeff=0.01*log10(e)
						else
							I_logTeff=InclogTeff
							I_Teff=Teff*log(10.)*I_logTeff
						end if
						Rc=sqrt(L/(Teff/TeffSun)**4)
						I_Rc=Rc*(0.5*I_L/L+2.*I_Teff/Teff)
						if (isEq(Rf,-1.D0,2)) then
							R1=Rc
							I_R1=I_Rc	
						else
							Ri=Rf
							I_Ri=SCP(22)
							call wMean(Ri,I_Ri,Rc,I_Rc,R1,I_R1)
						end if
						if (gProxyAvail) then !modified 18/9/2017 !Re-check compatibility between (rho,logg) 17/10/17
							if (loggAvail0 .and. rhoAvail0) then 
								loggf=loggf0
								rhof=rhof0
								R2=g0/10.**loggSun*rhoSun/rho0
								I_R2=R2*(I_g0/g0+I_rho0/rho0)
								call wMean(R1,I_R1,R2,I_R2,R,I_R)
								
								M=(g0/10.**loggSun)**3*(rhoSun/rho0)**2
								I_M=M*(3.*I_g0/g0+2.*I_rho0/rho0)
!!								L=R**2*(Teff/TeffSun)**4
!!								I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
!!								logL=log10(L)
!!								I_logL=I_L/L*log10(e)
								!!Modified to account for the mass uncertainty 26/10/2018 
								if (.not.(M+I_M.lt.minval(MTrAv) .or. M-I_M.gt.maxval(MTrAv))) then !inside MassRange
									call selectMfromTracks(MTrAv,M,I_M,Mvec)
									allocate(cumInt(size(Mvec)))
									jk=0
									do jj=1,size(Mvec)
										indxM=minloc(abs(Mvec(jj)-MTrAv),1)
										xMTi=Tndxi(indxZt(indxZ),indxM)
										xMTf=Tndxf(indxZt(indxZ),indxM)
										allocate(TrRhoG(xMTf-xMTi+1,size(TrackTab,2)))
										TrRhoG=TrackTab(xMTi:xMTf,:)
					
										allocate(MTr(size(TrRhoG,1))); allocate(RTr(size(TrRhoG,1)))
										allocate(logRhoTr(size(TrRhoG,1))); allocate(loggTr(size(TrRhoG,1)))
										
										MTr=TrRhoG(:,cM_T) !Msun
										RTr=10.**TrRhoG(:,clogR_T) !cm
										logRhoTr=log10(MTr/(RTr/RSun)**3*rhoSun)
										loggTr=log10(MTr/(RTr/RSun)**2)+loggSun
										call findYgivenX_v(TrRhoG,logg,loggTr,clogTe_T,logTeTrtmp)
										if (allocated(logTeTrtmp)) then
											jk=jk+1
											if (jk.eq.1) then
												allocate(logTeTr(size(logTeTrtmp)))
												logTeTr=logTeTrtmp
											else
												call append1D(logTeTr,logTeTrtmp)
											end if
											cumInt(jj)=size(logTeTr)
											deallocate(logTeTrtmp)
										else
											cumInt(jj)=0
										end if
										
										deallocate(TrRhoG)
										deallocate(MTr); deallocate(RTr)
										deallocate(logRhoTr); deallocate(loggTr)
									end do
									!
									logTeJ=logTeff !I don't use Johnson (1966) relation because
									! I've just calibrate the correct isochronal Teff
									if (allocated(logTeTr)) then 
										allocate(TeB(size(logTeTr)))
										TeB=((logTeTr>logTeJ-DTeJ).and.(logTeTr<logTeJ+DTeJ))
										call consistentM(TeB,cumInt,Mvec,M,agreeM)
										if (.not.any(TeB) .or. .not.agreeM) then !all elements of TeB are 0 => logg and rho inconsistent
											R=R1
											I_R=I_R1
											if (I_logg0>I_logrho0) then !logrho is best determined
											!Only using rho: determine M; re-determine logg (discard the input value) 
												loggf=-1.
												loggAvail=.false.
												loggAvailAS=.false.
												if (.not.isEq(Rf,-1.D0,2)) then
													call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M, &
														& g,I_g,logg,I_logg)
												else
													call computeMg(R,I_R,Teff,I_Teff,rho,I_rho,M,I_M,g,I_g,logg,I_logg)
												end if
												loggAvailCal=.true. !!logg has become available thanks to calibration
											else !logg is best determined
											!Only using logg: determine M; re-determine rho (discard the input value) 
												rhof=-1.
												rhoAvail=.false.
												rhoAvailAS=.false.
												if (.not.isEq(Rf,-1.D0,2)) then
													call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho, &
														& I_rho,logrho,I_logrho)
												else
													call computeMrho(R,I_R,Teff,I_Teff,g,I_g,M,I_M,rho,I_rho,logrho,I_logrho)
												end if
												rhoAvailCal=.true.
											end if 
										else !input logg and rho are consistent. Compute M through a weighted mean
											! between M=M(logg) and M=M(rho)
											Mlogg=g0/10.**loggSun*R**2
											I_Mlogg=Mlogg*(I_g0/g0+2.*I_R/R)
											Mrho=rho0/rhoSun*R**3 !Mo
											I_Mrho=Mrho*(I_rho0/rho0+3.*I_R/R)
											call wMean(Mlogg,I_Mlogg,Mrho,I_Mrho,M,I_M)
										end if
										deallocate(TeB); deallocate(logTeTr)
									else !Track has been opened, but no logTe compatible with logg were found
										R=R1
										I_R=I_R1
										if (I_logg0>I_logrho0) then !logrho is best determined
											!Only using rho: determine M; re-determine logg (discard the input value) 
											loggf=-1.
											loggAvail=.false.
											loggAvailAS=.false.
											if (.not.isEq(Rf,-1.D0,2)) then
												call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g &
													& ,logg,I_logg)
											else
												call computeMg(R,I_R,Teff,I_Teff,rho,I_rho,M,I_M,g,I_g,logg,I_logg)
											end if
											loggAvailCal=.true. !!logg now available thanks to calibration
										else !logg is best determined
											!Only using logg: determine M; re-determine rho (discard the input value) 
											rhof=-1.
											rhoAvail=.false.
											rhoAvailAS=.false.
											if (.not.isEq(Rf,-1.D0,2)) then
												call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M, &
													& rho,I_rho,logrho,I_logrho)
											else
												call computeMrho(R,I_R,Teff,I_Teff,g,I_g,M,I_M,rho,I_rho,logrho,I_logrho)
											end if
											rhoAvailCal=.true.
										end if 
									end if
									deallocate(cumInt); deallocate(Mvec)
								else !M is out of mass range of tracks
									R=R1
									I_R=I_R1
									if (I_logg0>I_logrho0) then !logrho is best determined
										!Only using rho: determine M; re-determine logg (discard the input value) 
										loggf=-1.
										loggAvail=.false.
										loggAvailAS=.false.
										if (.not.isEq(Rf,-1.D0,2)) then
											call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg, &
												& I_logg)
										else
											call computeMg(R,I_R,Teff,I_Teff,rho,I_rho,M,I_M,g,I_g,logg,I_logg)
										end if
										loggAvailCal=.true. !!logg now available thanks to calibration
									else !logg is best determined
										!Only using logg: determine M; re-determine rho (discard the input value) 
										rhof=-1.
										rhoAvail=.false.
										rhoAvailAS=.false.
										if (.not.isEq(Rf,-1.D0,2)) then
											call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho, &
												& I_rho,logrho,I_logrho)
										else
											call computeMrho(R,I_R,Teff,I_Teff,g,I_g,M,I_M,rho,I_rho,logrho,I_logrho)
										end if
										rhoAvailCal=.true.
									end if 
								end if 
								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
							else !only one among logg and rho is available
								R=R1
								I_R=I_R1
								if (loggAvail0) then !g still the original input one because it's not been overwritten
									if (.not.isEq(Rf,-1.D0,2)) then
										call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho, &
											& logrho,I_logrho)
									else
										call computeMrho(R,I_R,Teff,I_Teff,g,I_g,M,I_M,rho,I_rho,logrho,I_logrho)
									end if
									rhoAvailCal=.true.
								else !rhoAvail0 visto che vale gProxyAvail
									if (.not.isEq(Rf,-1.D0,2)) then
										call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
									else
										call computeMg(R,I_R,Teff,I_Teff,rho,I_rho,M,I_M,g,I_g,logg,I_logg)
									end if
									loggAvailCal=.true. !!logg now available thanks to calibration
								end if 
							end if !endif (loggAvail .and. rhoAvail)
						else
							R=R1
							I_R=I_R1
							if (.not.isEq(Rf,-1.D0,2)) then
								L=R**2*(Teff/TeffSun)**4
								I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
								logL=log10(L)
								I_logL=I_L/L*log10(e)
							end if
						end if !endif gProxy 
					else if (calibNoD) then 
						allocate(logy_iso(size(Isoc,1)))
						if (.not.isEq(Rf,-1.D0,2)) then
							y=logL
							yl=-y
							logy_iso=logL_iso
							y_lim=-1.75
						else if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
						& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))) then  !18/9/2017
							y=logg
							yl=y
							logy_iso=logg_iso
							y_lim=3.
						else !rhoAvail
							y=logrho
							yl=y
							logy_iso=logrho_iso
							y_lim=-2.2
						end if 
						logYuguali=.false.
						do jd=1,2
							allocate(logy_isor(tf(jd)-ti(jd)+1)); allocate(logT_isor(tf(jd)-ti(jd)+1))
							allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(Diff_logy(tf(jd)-ti(jd)+1))
				
							logy_isor=logy_iso(ti(jd):tf(jd))
							logT_isor=logTeff_iso(ti(jd):tf(jd))
							BmV_isor=BmV_iso(ti(jd):tf(jd))
							Diff_logy=sqrt((y-logy_isor)**2+(BmV-BmV_isor)**2) !!!!!!!!!!!!!!!!18/2/14
							if (jd==1) then 
								id(jd)=minloc(Diff_logy,1) 
								logy_is1(jd)=logy_isor(id(jd))
							else
								allocate(ndxD(size(Diff_logy)))
								call indexx(Diff_logy,ndxD) 
								rowD=size(ndxD)
								allocate(BmVD(size(BmV_isor)))
								BmVD=BmV_isor(ndxD)
								kD=1
								id(jd)=ndxD(kD)
								logy_is1(jd)=logy_isor(id(jd))
								do while ((abs(BmV_is1(1)-BmVD(kD))<=0.04 .or. abs(BmV_is2(1)-BmVD(kD))<=0.04) &
										& .and. kD<rowD .and. yl<y_lim) 
									kD=kD+1
									id(jd)=ndxD(kD)
									logy_is1(jd)=logy_isor(id(jd))
								end do
								deallocate(ndxD); deallocate(BmVD)
							end if 
							BmV_is1(jd)=BmV_isor(id(jd))
							logT_is1(jd)=logT_isor(id(jd))
							!!!
							call choose_i2xy1y2(y,logy_isor,BmV_isor,logT_isor,id(jd),hstar,hstarlim, & !BmV,
								& logy_is1(jd),x_is2,BmV_is1(jd),y1_is2,logT_is1(jd),y2_is2) !!!
							logy_is2(jd)=x_is2
							BmV_is2(jd)=y1_is2
							logT_is2(jd)=y2_is2
							!!!
							if (isEq(logy_is1(jd),logy_is2(jd),4)) then 
								logYuguali=.true.
								BmV_is(jd)=BmV
								call InterpLin_M(reshape((/BmV_is1(jd),BmV_is2(jd),logT_is1(jd),logT_is2(jd)/), &
									& (/2,2/)),BmV,1,1,(/2/),logT_is(jd),xlow,ylow,xup,yup)
							else
								BmV_is(jd)=(y-logy_is1(jd))*(BmV_is2(jd)-BmV_is1(jd))/ &
											& (logy_is2(jd)-logy_is1(jd))+BmV_is1(jd)
								logT_is(jd)=(y-logy_is1(jd))*(logT_is2(jd)-logT_is1(jd))/ &
											& (logy_is2(jd)-logy_is1(jd))+logT_is1(jd)
							end if 
				
							deallocate(BmV_isor); deallocate(logy_isor);
							deallocate(logT_isor); deallocate(Diff_logy)
						end do 
														 
						dCL_is=abs(BmV_is(1)-BmV_is(2))		!dist between isoc in the BmV-logg Diagram
						dHR_is=abs(logT_is(1)-logT_is(2))	!dist between isoc in the logTeff-logg Diagram
						if ((dCL_is>dlim .and. dHR_is>dlim)) then !condition added 17/2/14
			!				dCL=abs(BmV-BmV_is)				!vector of distances between the star and the two isoc
							if (BmV<BmV_is(1) .and. BmV<BmV_is(2)) then 
								state=1
							else if (BmV>BmV_is(1) .and. BmV>BmV_is(2)) then 
								state=2
							else
								state=3
							end if
							!!
							call calibrateDiag(state,dClim,BmV_is,logT_is,BmV,logTeff)
							if (idCol.eq.1) then
								xtau=BmV
							else
								xtau=logTeff
							end if
							!! 
!!!							if (isnan(logTeff)) then 
!!!								cycle
!!!							end if 
						else if (age_s==size(ix_sorted,1)) then 
							cycleW=1
							cycle
						else
							dCL_is=0.
							dHR_is=0.
						end if 
						do while ((dCL_is<=dlim .or. dHR_is<=dlim) .and. age_s<size(ix_sorted,1)) 
							age_s=age_s+1
							do while (dist(ix_sorted(age_s),2).lt.logtlimCal .and. age_s<size(ix_sorted))
								age_s=age_s+1
							end do
							if (isEq(age_iso1,agePivot,2)) then 
								age_iso2=dist(ix_sorted(age_s),2)
								if (age_iso1>age_iso2) then 
									tmp_age1=age_iso1
									age_iso1=age_iso2
									age_iso2=tmp_age1
									agePivot=age_iso2
								end if 
							else
								age_iso1=dist(ix_sorted(age_s),2)
								if (age_iso1>age_iso2) then 
									tmp_age1=age_iso1
									age_iso1=age_iso2
									age_iso2=tmp_age1
									agePivot=age_iso1
								end if 
							end if 
				
							call selectIsoc(logt_iso,age_iso1,logtstep,last_logt,ti0,tf0)
							ti(1)=ti0
							tf(1)=tf0
							call selectIsoc(logt_iso,age_iso2,logtstep,last_logt,ti0,tf0)
							ti(2)=ti0
							tf(2)=tf0
							 
							logYuguali=.false.
							do jd=1,2
								allocate(logy_isor(tf(jd)-ti(jd)+1)); allocate(logT_isor(tf(jd)-ti(jd)+1))
								allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(Diff_logy(tf(jd)-ti(jd)+1))
				
								logy_isor=logy_iso(ti(jd):tf(jd))
								logT_isor=logTeff_iso(ti(jd):tf(jd))
								BmV_isor=BmV_iso(ti(jd):tf(jd))
								Diff_logy=sqrt((y-logy_isor)**2+(BmV-BmV_isor)**2) !!!!!!!!!!!!!!!!18/2/14
								if (jd==1) then 
									id(jd)=minloc(Diff_logy,1) 
									logy_is1(jd)=logy_isor(id(jd))
								else
									allocate(ndxD(size(Diff_logy)))
									call indexx(Diff_logy,ndxD) 
									rowD=size(ndxD)
									allocate(BmVD(size(BmV_isor)))
									BmVD=BmV_isor(ndxD)
									kD=1
									id(jd)=ndxD(kD)
									logy_is1(jd)=logy_isor(id(jd))
									do while ((abs(BmV_is1(1)-BmVD(kD))<=0.04 .or. abs(BmV_is2(1)-BmVD(kD))<=0.04) &
											& .and. kD<rowD .and. yl<y_lim) 
										kD=kD+1
										id(jd)=ndxD(kD)
										logy_is1(jd)=logy_isor(id(jd))
									end do
									deallocate(ndxD); deallocate(BmVD)
								end if 
								BmV_is1(jd)=BmV_isor(id(jd))
								logT_is1(jd)=logT_isor(id(jd))
								!!!
								call choose_i2xy1y2(y,logy_isor,BmV_isor,logT_isor,id(jd),hstar,hstarlim, & !BmV,
									& logy_is1(jd),x_is2,BmV_is1(jd),y1_is2,logT_is1(jd),y2_is2) !!!
								logy_is2(jd)=x_is2
								BmV_is2(jd)=y1_is2
								logT_is2(jd)=y2_is2
								!!!
								if (isEq(logy_is1(jd),logy_is2(jd),4)) then 
									logYuguali=.true.
									BmV_is(jd)=BmV
									call InterpLin_M(reshape((/BmV_is1(jd),BmV_is2(jd),logT_is1(jd),logT_is2(jd)/), &
										& (/2,2/)),BmV,1,1,(/2/),logT_is(jd),xlow,ylow,xup,yup)
								else
									BmV_is(jd)=(y-logy_is1(jd))*(BmV_is2(jd)-BmV_is1(jd))/ &
												& (logy_is2(jd)-logy_is1(jd))+BmV_is1(jd)
									logT_is(jd)=(y-logy_is1(jd))*(logT_is2(jd)-logT_is1(jd))/ &
												& (logy_is2(jd)-logy_is1(jd))+logT_is1(jd)
								end if 
				
								deallocate(BmV_isor); deallocate(logy_isor);
								deallocate(logT_isor); deallocate(Diff_logy)
							end do 
							
							dCL_is=abs(BmV_is(1)-BmV_is(2))		!dist between isoc in the BmV-logg Diagram
							dHR_is=abs(logT_is(1)-logT_is(2))	!dist between isoc in the logTeff-logg Diagram
							if ((dCL_is>dlim .and. dHR_is>dlim)) then !condition added 17/2/14
				!				dCL=abs(BmV-BmV_is)				!vector of distances between the star and the two isoc
								if (BmV<BmV_is(1) .and. BmV<BmV_is(2)) then 
									state=1
								else if (BmV>BmV_is(1) .and. BmV>BmV_is(2)) then 
									state=2
								else
									state=3
								end if
								!!
								call calibrateDiag(state,dClim,BmV_is,logT_is,BmV,logTeff)
								if (idCol.eq.1) then
									xtau=BmV
								else
									xtau=logTeff
								end if
								!! 
!!!								if (isnan(logTeff)) then 
!!!									cycle
!!!								end if 
							else if (age_s==size(ix_sorted,1)) then 
								cycleW=1
								cycle
							end if 
						end do
						if (cycleW.eq.1) then
							cycle
						end if
						Teff=10.**logTeff
						if ((.not.calibSPEC .and. Teffinput.ne.-1)) then 
							!Inside calibNoD => .not.calibSPEC is always true
							if (abs(Teff-Teffinput)>300) then !Inconsistency between BmV and spectroscopic Teff
								cycle
							end if 
						end if 
						I_Teff=0.01*Teff
						I_logTeff=0.01*log10(e)
						!!Compatibility
						if (loggAvail0.and.rhoAvail0) then 
							!!rho, logg both available. I determine input values for M, R 
							loggf=loggf0
							rhof=rhof0
							R2=g0/10.**loggSun*rhoSun/rho0
							I_R2=R2*(I_g0/g0+I_rho0/rho0)
							if (isEq(Rf,-1.D0,2)) then
								R=R2
								I_R=I_R2
							else
								call wMean(R1,I_R1,R2,I_R2,R,I_R)
							end if
							M=(g0/10.**loggSun)**3*(rhoSun/rho0)**2
							I_M=M*(3.*I_g0/g0+2.*I_rho0/rho0)
							L=R**2*(Teff/TeffSun)**4
							I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
							logL=log10(L)
							I_logL=I_L/L*log10(e)
							!!!!!!Check compatibility between rho and logg 16/10/17 
							!!Modified to account for the mass uncertainty 26/10/2018 
							if (.not.(M+I_M.lt.minval(MTrAv) .or. M-I_M.gt.maxval(MTrAv))) then !open tracks only if inside MassRange 
								call selectMfromTracks(MTrAv,M,I_M,Mvec)
								allocate(cumInt(size(Mvec)))
								jk=0
								do jj=1,size(Mvec)
									indxM=minloc(abs(Mvec(jj)-MTrAv),1)
									xMTi=Tndxi(indxZt(indxZ),indxM)
									xMTf=Tndxf(indxZt(indxZ),indxM)
									allocate(TrRhoG(xMTf-xMTi+1,size(TrackTab,2)))
									TrRhoG=TrackTab(xMTi:xMTf,:)
				
									allocate(MTr(size(TrRhoG,1))); allocate(RTr(size(TrRhoG,1)))
									allocate(logRhoTr(size(TrRhoG,1))); allocate(loggTr(size(TrRhoG,1)))
									
									MTr=TrRhoG(:,cM_T) !Msun
									RTr=10.**TrRhoG(:,clogR_T) !cm
									logRhoTr=log10(MTr/(RTr/RSun)**3*rhoSun)
									loggTr=log10(MTr/(RTr/RSun)**2)+loggSun
									call findYgivenX_v(TrRhoG,logg,loggTr,clogTe_T,logTeTrtmp)
									if (allocated(logTeTrtmp)) then
										jk=jk+1
										if (jk.eq.1) then
											allocate(logTeTr(size(logTeTrtmp)))
											logTeTr=logTeTrtmp
										else
											call append1D(logTeTr,logTeTrtmp)
										end if
										cumInt(jj)=size(logTeTr)
										deallocate(logTeTrtmp)
									else
										cumInt(jj)=0
									end if
									
									deallocate(TrRhoG)
									deallocate(MTr); deallocate(RTr)
									deallocate(logRhoTr); deallocate(loggTr)
								end do
								logTeJ=logTeff !I don't use Johnson (1966) relation because I've just calibrated Teff
								DTeJ=I_logTeff !directly set to the true uncertainty
								if (allocated(logTeTr)) then 
									allocate(TeB(size(logTeTr)))
									TeB=((logTeTr>logTeJ-DTeJ).and.(logTeTr<logTeJ+DTeJ))
									if (.not.any(TeB)) then !all elements of TeB are 0 => logg and rho inconsistent
										if (I_logg0>I_logrho0) then !logrho is best determined
											if (isEq(Rf,-1.D0,2)) then
												call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
												call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
												loggf=-1.
												call setNaN(logg); call setNaN(I_logg)
												call setNaN(g); call setNaN(I_g)
												loggAvail=.false.
												loggAvailAS=.false.
											else !input R is available
												!Only using rho: determine M; re-determine logg (discard the input value) 
												loggf=-1.
												loggAvail=.false.
												loggAvailAS=.false.
												R=R1 !use just the input reliable value
												I_R=I_R1
												call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
												loggAvailCal=.true. !!logg now available thanks to calibration
											end if
										else !logg is best determined
											if (isEq(Rf,-1.D0,2)) then
												call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
												call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
												rhof=-1.
												call setNaN(rho); call setNaN(I_rho)
												call setNaN(logrho); call setNaN(I_logrho)
												rhoAvail=.false.
												rhoAvailAS=.false.
											else !input R is available
												!Only using logg: determine M; re-determine rho (discard the input value) 
												rhof=-1.
												rhoAvail=.false.
												rhoAvailAS=.false.
												R=R1 !use just the input reliable value
												I_R=I_R1
												call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
												rhoAvailCal=.true.
											end if
										end if 
									else
										!!!!!!!!!
										call consistentM(TeB,cumInt,Mvec,M,agreeM)
										if (.not.agreeM) then !new M has been recomputed inside the subroutine in case agreeM is false
															  !This M will be used in case R is not available directly from input
											if (I_logg0.gt.I_logrho0) then !logrho is best determined
												if (isEq(Rf,-1.D0,2)) then
													loggAvail=.false.
													loggAvailAS=.false.
													!g: mantain original input uncertainty
													call computeLRgfromMrho(M,I_M,rho,I_rho,Teff,I_Teff,L,logL,I_L,I_logL,R,I_R,g,logg)
													loggAvailCal=.true.
												else !input R is available
													!Only using rho: determine M; re-determine logg (discard the input value) 
													loggf=-1.
													loggAvail=.false.
													loggAvailAS=.false.
													R=R1
													I_R=I_R1
													call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
													loggAvailCal=.true. !!logg now available thanks to calibration
												end if
											else
												if (isEq(Rf,-1.D0,2)) then
													rhoAvail=.false.
													rhoAvailAS=.false.
													!rho in g/cm3. Mantain original input uncertainty on rho
													call computeLRrhofromMg(M,I_M,g,I_g,Teff,I_Teff,L,logL,I_L,I_logL,R,I_R,rho,logrho)
													rhoAvailCal=.true.
												else !input R is available
													!Only using logg: determine M; re-determine rho (discard the input value) 
													rhof=-1.
													rhoAvail=.false.
													rhoAvailAS=.false.
													R=R1
													I_R=I_R1
													call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
													rhoAvailCal=.true.
												end if
											end if
											if (isEq(Rf,-1.D0,2)) then
												R=g/10**loggSun*rhoSun/rho
												I_R=R*(I_g/g+I_rho/rho)
												L=R**2*(Teff/TeffSun)**4
												I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
												logL=log10(L)
												I_logL=I_L/L*log10(e)
											end if
										else
											if (I_logg0.lt.I_logrho0) then
												bestLogg=.true.
											else
												bestLogg=.false.
											end if
											if (.not.isEq(Rf,-1.D0,2)) then
											!input R is available: weighted mean to infer M
												Mlogg=g0/10.**loggSun*R**2
												I_Mlogg=Mlogg*(I_g0/g0+2.*I_R/R)
												Mrho=rho0/rhoSun*R**3 !Mo
												I_Mrho=Mrho*(I_rho0/rho0+3.*I_R/R)
												call wMean(Mlogg,I_Mlogg,Mrho,I_Mrho,M,I_M)
												!L already defined since the beginning from R
											end if
										end if
										!!!!!!!!!!!!!
									end if
									deallocate(logTeTr); deallocate(TeB)
								else !Track has been opened, but no logTe compatible with logg were found
									if (I_logg0>I_logrho0) then !logrho is best determined
										if (isEq(Rf,-1.D0,2)) then
											call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
											call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
											loggf=-1.
											call setNaN(logg); call setNaN(I_logg)
											call setNaN(g); call setNaN(I_g)
											loggAvail=.false.
											loggAvailAS=.false.
										else !input R is available
											!Only using rho: determine M; re-determine logg (discard the input value) 
											loggf=-1.
											loggAvail=.false.
											loggAvailAS=.false.
											R=R1 !use just the input reliable value
											I_R=I_R1
											call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
											loggAvailCal=.true. !!logg now available thanks to calibration
										end if
									else !logg is best determined
										if (isEq(Rf,-1.D0,2)) then
											call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
											call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
											rhof=-1.
											call setNaN(rho); call setNaN(I_rho)
											call setNaN(logrho); call setNaN(I_logrho)
											rhoAvail=.false.
											rhoAvailAS=.false.
										else !input R is available
											!Only using logg: determine M; re-determine rho (discard the input value) 
											rhof=-1.
											rhoAvail=.false.
											rhoAvailAS=.false.
											R=R1 !use just the input reliable value
											I_R=I_R1
											call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
											rhoAvailCal=.true.
										end if
									end if 
								end if
								deallocate(cumInt) 
							else
								if (I_logg0>I_logrho0) then !logrho is best determined
									if (isEq(Rf,-1.D0,2)) then
										call setNan(R); call setNan(I_R); call setNan(M); call setNan(I_M)
										call setNan(L); call setNan(I_L); call setNan(logL); call setNan(I_logL)
										loggf=-1.
										call setNan(logg); call setNan(I_logg)
										call setNan(g); call setNan(I_g)
										loggAvail=.false.
										loggAvailAS=.false.
									else !input R is available
										!Only using rho: determine M; re-determine logg (discard the input value) 
										loggf=-1.
										loggAvail=.false.
										loggAvailAS=.false.
										R=R1 !use just the input reliable value
										I_R=I_R1
										call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
										loggAvailCal=.true. !!logg now available thanks to calibration
									end if
								else !logg is best determined
									if (isEq(Rf,-1.D0,2)) then
										call setNan(R); call setNan(I_R); call setNan(M); call setNan(I_M)
										call setNan(L); call setNan(I_L); call setNan(logL); call setNan(I_logL)
										rhof=-1.
										call setNan(rho); call setNan(I_rho)
										call setNan(logrho); call setNan(I_logrho)
										rhoAvail=.false.
										rhoAvailAS=.false.
									else !input R is available
										!Only using logg: determine M; re-determine rho (discard the input value) 
										rhof=-1.
										rhoAvail=.false.
										rhoAvailAS=.false.
										R=R1 !use just the input reliable value
										I_R=I_R1
										call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
										rhoAvailCal=.true.
									end if
								end if 
							end if 
							!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
						else
							if (.not.isEq(Rf,-1.D0,2)) then
								R=R1
								I_R=I_R1
								if (loggAvail.or.loggAvailAS) then
									call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
									rhoAvailCal=.true.
								else if (rhoAvail.or.rhoAvailAS) then
									call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
									loggAvailCal=.true. !!logg available thanks to calibration
								else !noGproxy
									L=R**2*(Teff/TeffSun)**4
									I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
									logL=log10(L)
									I_logL=I_L/L*log10(e)
								end if
							end if
						end if						
						!!End compatibility
						if (.not.isEq(Rf,-1.D0,2)) then
							kwl=0
							cycleW=0
							age_s=1
							do while (abs(logL-y).gt.DlogL.and.kwl.lt.10.and.age_s<size(ix_sorted))
								y=logL
								kwl=kwl+1
								deallocate(t_i); deallocate(rho_i); deallocate(dist);deallocate(ix_sorted)
								deallocate(logTeff_i); deallocate(logL_i); deallocate(M_i); deallocate(logg_i)
								deallocate(logrho_i); deallocate(Teff_i); deallocate(g_i); deallocate(L_i)
								deallocate(logy_iso)
								if (photIsocAvail) then
									deallocate(BC_i); deallocate(BmV_i)
								end if
								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
								!!!!!!!NEW CODE 12/2/14!!!!!!!!!!!!!!! 
								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
								allocate(dist_ndx(size(logTeff_iso))) !equal to any other size of _iso
								allocate(x_iso(size(logTeff_iso))); allocate(xx_iso(size(logTeff_iso)))
								call indexx(sqrt((BmV-BmV_iso)**2+(logL-logL_iso)**2),dist_ndx)
								x=BmV
								x_iso=BmV_iso
								xx_iso=logTeff_iso
								cBrif=cBBmV
								allocate(ageMinDist(size(logt_iso)))
								ageMinDist=logt_iso(dist_ndx)
								
								call uniqueFast(ageMinDist,2,ix,.true.)
								call sort0(ix) !so to display elements of ageMinDist in the oroginal order
								allocate(dist_ndxU(size(ix))) 
								dist_ndxU=dist_ndx(ix)
								
								deallocate(dist_ndx); deallocate(ageMinDist)
								deallocate(ix)
								rowAge=1
								! allocation of already defined variables
								allocate(t_iTmp(size(dist_ndxU)));allocate(distTmp(size(dist_ndxU),2))
								allocate(logTeff_iTmp(size(dist_ndxU)));allocate(logL_iTmp(size(dist_ndxU)))
								allocate(M_iTmp(size(dist_ndxU)));allocate(logg_iTmp(size(dist_ndxU)))
								allocate(logrho_iTmp(size(dist_ndxU)))
								if (photIsocAvail) then 
									allocate(BmV_iTmp(size(dist_ndxU)));allocate(BC_iTmp(size(dist_ndxU)))
								end if
								!
								do Nage=1,size(dist_ndxU) 
									ndxDU=dist_ndxU(Nage)
									if (photIsocAvail) then 
										BmV_i1=BmV_iso(ndxDU)
										V_i1=V_iso(ndxDU)
									end if 
									logTeff_i1=logTeff_iso(ndxDU)
									logL_i1=logL_iso(ndxDU)
									M_i1=M_iso(ndxDU)
									logg_i1=logg_iso(ndxDU)
									logrho_i1=logrho_iso(ndxDU)
									
									call choose_i2Col(ndxDU,x,x_iso,xx_iso,calibSPEC,hstar,hstarlim,V_iso,logL_iso, &
										& M_iso,logg_iso,logrho_iso,logt_iso,BmV_i2,V_i2,logTeff_i2,logL_i2,M_i2,logg_i2,logrho_i2)
									if (calibSPEC) then
										x_i1=logTeff_i1
										x_i2=logTeff_i2
									else
										x_i1=BmV_i1
										x_i2=BmV_i2
									end if
									if (isEq(x_i1,x_i2,4)) then 
										cycle
									end if 
									t_iTmp(rowAge)=t_iso(ndxDU)
									if (photIsocAvail) then 
										BC_iTmp(rowAge)=Isoc(ndxDU,cmbol)-Isoc(ndxDU,cmag0) !just rough.
										!When I recover BC in the logg-BmV plane, I'll make the interpolation using BC_i1 e BC_i2
									end if 
									
									call findX_i(BmV,logL,BmV_i1,BmV_i2,logL_i1,logL_i2,x_i,dist0)
									BmV_iTmp(rowAge)=x_i
									
									distTmp(rowAge,:)=(/ dist0, logt_iso(ndxDU) /)
												
									!!18/3/15 The following if cancel the criterion of perpendicular distance
									!!between a star and the isochrone if the theoretical point to be interpolated
									!!doesn't fall inside the isochrone (i.e. between x_i1 and x_i2), but on its
									!!extension. This possibility, in fact, would build fictitious isochrones that
									!!could be erroneously close to the star
									if (.not.((x_i>=x_i1.and.x_i<=x_i2).or.(x_i>=x_i2.and.x_i<=x_i1))) then 
										if (photIsocAvail) then 
											BmV_iTmp(rowAge)=BmV_i1
										end if 
										logg_iTmp(rowAge)=logg_i1
										logTeff_iTmp(rowAge)=logTeff_i1
										logL_iTmp(rowAge)=logL_i1
										M_iTmp(rowAge)=M_i1
										logrho_iTmp(rowAge)=logrho_i1
									else
										mT=(logTeff_i2-logTeff_i1)/(x_i2-x_i1)
										qT=-mT*x_i1+logTeff_i1
										logTeff_iTmp(rowAge)=mT*x_i+qT
										
										mg=(logg_i2-logg_i1)/(x_i2-x_i1)
										qg=-mg*x_i1+logg_i1
										logg_iTmp(rowAge)=mg*x_i+qg
										mL=(logL_i2-logL_i1)/(x_i2-x_i1)
										qL=-mL*x_i1+logL_i1
										logL_iTmp(rowAge)=mL*x_i+qL
										mM=(M_i2-M_i1)/(x_i2-x_i1)
										qM=-mM*x_i1+M_i1
										M_iTmp(rowAge)=mM*x_i+qM
										mrh=(logrho_i2-logrho_i1)/(x_i2-x_i1)
										qrho=-mrh*x_i1+logrho_i1
										logrho_iTmp(rowAge)=mrh*x_i+qrho
									end if 

									rowAge=rowAge+1
								end do
								
								deallocate(x_iso); deallocate(xx_iso)
										
								allocate(t_i(rowAge-1));allocate(dist(rowAge-1,2))
								allocate(logTeff_i(rowAge-1));allocate(logL_i(rowAge-1))
								allocate(M_i(rowAge-1));allocate(logg_i(rowAge-1))
								allocate(logrho_i(rowAge-1))
								allocate(BmV_i(rowAge-1));allocate(BC_i(rowAge-1))
								t_i=t_iTmp(1:rowAge-1)
								dist=distTmp(1:rowAge-1,:)
								logTeff_i=logTeff_iTmp(1:rowAge-1)
								logL_i=logL_iTmp(1:rowAge-1)
								M_i=M_iTmp(1:rowAge-1)
								logg_i=logg_iTmp(1:rowAge-1)
								logrho_i=logrho_iTmp(1:rowAge-1)
								BmV_i=BmV_iTmp(1:rowAge-1)
								BC_i=BC_iTmp(1:rowAge-1)
								deallocate(t_iTmp);deallocate(distTmp);deallocate(logTeff_iTmp)
								deallocate(logL_iTmp);deallocate(M_iTmp);deallocate(logg_iTmp)
								deallocate(logrho_iTmp)
								if (photIsocAvail) then
									deallocate(BmV_iTmp);deallocate(BC_iTmp)
								end if
								
								allocate(Teff_i(rowAge-1));allocate(L_i(rowAge-1))
								allocate(g_i(rowAge-1));allocate(rho_i(rowAge-1))
								
								Teff_i=10.**logTeff_i
								L_i=10.**logL_i
								g_i=10.**logg_i
								rho_i=10.**logrho_i
								
								deallocate(dist_ndxU)
								
								allocate(dsorted(size(dist)))
								dsorted=dist(:,1) !then it will be overwritten by sort so that it will be actually sorted
								call sort1(dsorted)
								allocate(ix_sorted(size(dist,1)))
								call indexx(dist(:,1), ix_sorted) 
								!!Avoid extremely young isochrones (t<10 Myr) in the following calibration process
								age_s1=1 
								age_iso1=dist(ix_sorted(age_s1),2)
								do while (age_iso1.lt.logtlimCal)
									age_s1=age_s1+1
									age_iso1=dist(ix_sorted(age_s1),2)
								end do
								age_s=age_s1+1
								if (dsorted(age_s1)<1.e-6) then 
									do while ((log10(dsorted(age_s))-log10(dsorted(age_s1))<4.5.or.dist(ix_sorted(age_s),2).lt.logtlimCal) &
											& .and. age_s<size(ix_sorted))
										age_s=age_s+1
									end do 
								else
									do while (dist(ix_sorted(age_s),2).lt.logtlimCal .and. age_s<size(ix_sorted))
										age_s=age_s+1
									end do
								end if
								deallocate(dsorted)
								age_iso2=dist(ix_sorted(age_s),2)
								if (age_iso1>age_iso2) then 
									tmp_age1=age_iso1
									age_iso1=age_iso2
									age_iso2=tmp_age1 !age_iso2 BECOMES the isoc with the minimum distance from the star
									agePivot=age_iso2
								else
									agePivot=age_iso1 !age_iso1 REMAINS the isoc with the minimum distance from the star
								end if
								
								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
								!!!!!!!!!!!END NEW CODE 
								!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
								
								call selectIsoc(logt_iso,age_iso1,logtstep,last_logt,ti0,tf0)
								ti(1)=ti0
								tf(1)=tf0
								call selectIsoc(logt_iso,age_iso2,logtstep,last_logt,ti0,tf0)
								ti(2)=ti0
								tf(2)=tf0
																
								allocate(logy_iso(size(Isoc,1)))
								
								y=logL !already stated before
								yl=-y
								logy_iso=logL_iso
								y_lim=-1.75
								
								logYuguali=.false.
								do jd=1,2
									allocate(logy_isor(tf(jd)-ti(jd)+1)); allocate(logT_isor(tf(jd)-ti(jd)+1))
									allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(Diff_logy(tf(jd)-ti(jd)+1))
									
									logy_isor=logy_iso(ti(jd):tf(jd))
									logT_isor=logTeff_iso(ti(jd):tf(jd))
									BmV_isor=BmV_iso(ti(jd):tf(jd))
									Diff_logy=sqrt((y-logy_isor)**2+(BmV-BmV_isor)**2) !!!!!!!!!!!!!!!!18/2/14
									if (jd==1) then 
										id(jd)=minloc(Diff_logy,1) 
										logy_is1(jd)=logy_isor(id(jd))
									else
										allocate(ndxD(size(Diff_logy)))
										call indexx(Diff_logy,ndxD) 
										rowD=size(ndxD)
										allocate(BmVD(size(BmV_isor)))
										BmVD=BmV_isor(ndxD)
										kD=1
										id(jd)=ndxD(kD)
										logy_is1(jd)=logy_isor(id(jd))
										do while ((abs(BmV_is1(1)-BmVD(kD))<=0.04 .or. abs(BmV_is2(1)-BmVD(kD))<=0.04) &
												& .and. kD<rowD .and. yl<y_lim) 
											kD=kD+1
											id(jd)=ndxD(kD)
											logy_is1(jd)=logy_isor(id(jd))
										end do
										deallocate(ndxD); deallocate(BmVD)
									end if 
									BmV_is1(jd)=BmV_isor(id(jd))
									logT_is1(jd)=logT_isor(id(jd))
									!!!
									call choose_i2xy1y2(y,logy_isor,BmV_isor,logT_isor,id(jd),hstar,hstarlim, & !BmV,
										& logy_is1(jd),x_is2,BmV_is1(jd),y1_is2,logT_is1(jd),y2_is2) !!!
									logy_is2(jd)=x_is2
									BmV_is2(jd)=y1_is2
									logT_is2(jd)=y2_is2
									!!!
									if (isEq(logy_is1(jd),logy_is2(jd),4)) then 
										logYuguali=.true.
										BmV_is(jd)=BmV
										call InterpLin_M(reshape((/BmV_is1(jd),BmV_is2(jd),logT_is1(jd),logT_is2(jd)/), &
											& (/2,2/)),BmV,1,1,(/2/),logT_is(jd),xlow,ylow,xup,yup)
									else
										BmV_is(jd)=(y-logy_is1(jd))*(BmV_is2(jd)-BmV_is1(jd))/(logy_is2(jd)-logy_is1(jd)) &
													& +BmV_is1(jd)
										logT_is(jd)=(y-logy_is1(jd))*(logT_is2(jd)-logT_is1(jd))/(logy_is2(jd)-logy_is1(jd)) &
													& +logT_is1(jd)
									end if
									
									deallocate(BmV_isor); deallocate(logy_isor);
									deallocate(logT_isor); deallocate(Diff_logy)
								end do 
																		 
								dCL_is=abs(BmV_is(1)-BmV_is(2))		!dist between isoc in the BmV-logg Diagram
								dHR_is=abs(logT_is(1)-logT_is(2))	!dist between isoc in the logTeff-logg Diagram
								if ((dCL_is>dlim .and. dHR_is>dlim)) then !condition added 17/2/14
						!				dCL=abs(BmV-BmV_is)				!vector of distances between the star and the two isoc
									if (BmV<BmV_is(1) .and. BmV<BmV_is(2)) then 
										state=1
									else if (BmV>BmV_is(1) .and. BmV>BmV_is(2)) then 
										state=2
									else
										state=3
									end if
									!!
									call calibrateDiag(state,dClim,BmV_is,logT_is,BmV,logTeff)
									if (idCol.eq.1) then
										xtau=BmV
									else
										xtau=logTeff
									end if
									!! 
!!!									if (isnan(logTeff)) then 
!!!										cycle
!!!									end if 
								else if (age_s==size(ix_sorted,1)) then 
									cycleW=1
									cycle
								else
									dCL_is=0.
									dHR_is=0.
								end if 
								do while ((dCL_is<=dlim .or. dHR_is<=dlim) .and. age_s<size(ix_sorted,1)) 
									age_s=age_s+1
									do while (dist(ix_sorted(age_s),2).lt.logtlimCal .and. age_s<size(ix_sorted))
										age_s=age_s+1
									end do
									if (isEq(age_iso1,agePivot,2)) then 
										age_iso2=dist(ix_sorted(age_s),2)
										if (age_iso1>age_iso2) then 
											tmp_age1=age_iso1
											age_iso1=age_iso2
											age_iso2=tmp_age1
											agePivot=age_iso2
										end if 
									else
										age_iso1=dist(ix_sorted(age_s),2)
										if (age_iso1>age_iso2) then 
											tmp_age1=age_iso1
											age_iso1=age_iso2
											age_iso2=tmp_age1
											agePivot=age_iso1
										end if 
									end if 
									
									call selectIsoc(logt_iso,age_iso1,logtstep,last_logt,ti0,tf0)
									ti(1)=ti0
									tf(1)=tf0
									call selectIsoc(logt_iso,age_iso2,logtstep,last_logt,ti0,tf0)
									ti(2)=ti0
									tf(2)=tf0
									 
									logYuguali=.false.
									do jd=1,2
										allocate(logy_isor(tf(jd)-ti(jd)+1)); allocate(logT_isor(tf(jd)-ti(jd)+1))
										allocate(BmV_isor(tf(jd)-ti(jd)+1)); allocate(Diff_logy(tf(jd)-ti(jd)+1))
									
										logy_isor=logy_iso(ti(jd):tf(jd))
										logT_isor=logTeff_iso(ti(jd):tf(jd))
										BmV_isor=BmV_iso(ti(jd):tf(jd))
										Diff_logy=sqrt((y-logy_isor)**2+(BmV-BmV_isor)**2) !!!!!!!!!!!!!!!!18/2/14
										if (jd==1) then 
											id(jd)=minloc(Diff_logy,1) 
											logy_is1(jd)=logy_isor(id(jd))
										else
											allocate(ndxD(size(Diff_logy)))
											call indexx(Diff_logy,ndxD) 
											rowD=size(ndxD)
											allocate(BmVD(size(BmV_isor)))
											BmVD=BmV_isor(ndxD)
											kD=1
											id(jd)=ndxD(kD)
											logy_is1(jd)=logy_isor(id(jd))
											do while ((abs(BmV_is1(1)-BmVD(kD))<=0.04 .or. abs(BmV_is2(1)-BmVD(kD))<=0.04) &
													& .and. kD<rowD .and. yl<y_lim) 
												kD=kD+1
												id(jd)=ndxD(kD)
												logy_is1(jd)=logy_isor(id(jd))
											end do
											deallocate(ndxD); deallocate(BmVD)
										end if 
										BmV_is1(jd)=BmV_isor(id(jd))
										logT_is1(jd)=logT_isor(id(jd))
										!!!
										call choose_i2xy1y2(y,logy_isor,BmV_isor,logT_isor,id(jd),hstar,hstarlim, & !BmV,
											& logy_is1(jd),x_is2,BmV_is1(jd),y1_is2,logT_is1(jd),y2_is2) !!!
										logy_is2(jd)=x_is2
										BmV_is2(jd)=y1_is2
										logT_is2(jd)=y2_is2
										!!!
										if (isEq(logy_is1(jd),logy_is2(jd),4)) then 
											logYuguali=.true.
											BmV_is(jd)=BmV
											call InterpLin_M(reshape((/BmV_is1(jd),BmV_is2(jd),logT_is1(jd),logT_is2(jd)/), &
												& (/2,2/)),BmV,1,1,(/2/),logT_is(jd),xlow,ylow,xup,yup)
										else
											BmV_is(jd)=(y-logy_is1(jd))*(BmV_is2(jd)-BmV_is1(jd))/(logy_is2(jd)-logy_is1(jd)) &
														& +BmV_is1(jd)
											logT_is(jd)=(y-logy_is1(jd))*(logT_is2(jd)-logT_is1(jd))/(logy_is2(jd)-logy_is1(jd)) &
														& +logT_is1(jd)
										end if 
									
										deallocate(BmV_isor); deallocate(logy_isor);
										deallocate(logT_isor); deallocate(Diff_logy)
									end do 
									
									dCL_is=abs(BmV_is(1)-BmV_is(2))		!dist between isoc in the BmV-logg Diagram
									dHR_is=abs(logT_is(1)-logT_is(2))	!dist between isoc in the logTeff-logg Diagram
									if ((dCL_is>dlim .and. dHR_is>dlim)) then !condition added 17/2/14
						!				dCL=abs(BmV-BmV_is)				!vector of distances between the star and the two isoc
										if (BmV<BmV_is(1) .and. BmV<BmV_is(2)) then 
											state=1
										else if (BmV>BmV_is(1) .and. BmV>BmV_is(2)) then 
											state=2
										else
											state=3
										end if
										!!
										call calibrateDiag(state,dClim,BmV_is,logT_is,BmV,logTeff)
										if (idCol.eq.1) then
											xtau=BmV
										else
											xtau=logTeff
										end if
										!! 
!!!										if (isnan(logTeff)) then 
!!!											cycle
!!!										end if 
									else if (age_s==size(ix_sorted,1)) then 
										cycleW=1
										cycle
									end if 
								end do
								if (cycleW.eq.1) then
									cycle
								end if 
								Teff=10.**logTeff
								if ((.not.calibSPEC .and. Teffinput.ne.-1)) then !Inside calibNoD => .not.calibSPEC is always true
									if (abs(Teff-Teffinput)>300) then !Inconsistency between BmV and spectrscopic Teff
										cycle
									end if 
								end if 
								I_Teff=0.01*Teff
								I_logTeff=0.01*log10(e)
								!!Compatibility
								if (loggAvail0.and.rhoAvail0) then 
									!!rho, logg both available. I determine input values for M, R 
									loggf=loggf0
									rhof=rhof0
									R2=g0/10.**loggSun*rhoSun/rho0
									I_R2=R2*(I_g0/g0+I_rho0/rho0)
									if (isEq(Rf,-1.D0,2)) then
										R=R2
										I_R=I_R2
									else
										call wMean(R1,I_R1,R2,I_R2,R,I_R)
									end if
									M=(g0/10.**loggSun)**3*(rhoSun/rho0)**2
									I_M=M*(3.*I_g0/g0+2.*I_rho0/rho0)
									L=R**2*(Teff/TeffSun)**4
									I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
									logL=log10(L)
									I_logL=I_L/L*log10(e)
									!!!!!!Check compatibility between rho and logg 16/10/17 
									!!Modified to account for the mass uncertainty 26/10/2018 
									if (.not.(M+I_M.lt.minval(MTrAv) .or. M-I_M.gt.maxval(MTrAv))) then !open tracks only if inside MassRange 
										call selectMfromTracks(MTrAv,M,I_M,Mvec)
										allocate(cumInt(size(Mvec)))
										jk=0
										do jj=1,size(Mvec)
											indxM=minloc(abs(Mvec(jj)-MTrAv),1)
											xMTi=Tndxi(indxZt(indxZ),indxM)
											xMTf=Tndxf(indxZt(indxZ),indxM)
											allocate(TrRhoG(xMTf-xMTi+1,size(TrackTab,2)))
											TrRhoG=TrackTab(xMTi:xMTf,:)
						
											allocate(MTr(size(TrRhoG,1))); allocate(RTr(size(TrRhoG,1)))
											allocate(logRhoTr(size(TrRhoG,1))); allocate(loggTr(size(TrRhoG,1)))
											
											MTr=TrRhoG(:,cM_T) !Msun
											RTr=10.**TrRhoG(:,clogR_T) !cm
											logRhoTr=log10(MTr/(RTr/RSun)**3*rhoSun)
											loggTr=log10(MTr/(RTr/RSun)**2)+loggSun
											call findYgivenX_v(TrRhoG,logg,loggTr,clogTe_T,logTeTrtmp)
											if (allocated(logTeTrtmp)) then
												jk=jk+1
												if (jk.eq.1) then
													allocate(logTeTr(size(logTeTrtmp)))
													logTeTr=logTeTrtmp
												else
													call append1D(logTeTr,logTeTrtmp)
												end if
												cumInt(jj)=size(logTeTr)
												deallocate(logTeTrtmp)
											else
												cumInt(jj)=0
											end if
											
											deallocate(TrRhoG)
											deallocate(MTr); deallocate(RTr)
											deallocate(logRhoTr); deallocate(loggTr)
										end do
										logTeJ=logTeff !I don't use Johnson (1966) relation because I've just calibrated Teff
										DTeJ=I_logTeff !directly set to the true uncertainty
										if (allocated(logTeTr)) then 
											allocate(TeB(size(logTeTr)))
											TeB=((logTeTr>logTeJ-DTeJ).and.(logTeTr<logTeJ+DTeJ))
											if (.not.any(TeB)) then !all elements of TeB are 0 => logg and rho inconsistent
												if (I_logg0>I_logrho0) then !logrho is best determined
													if (isEq(Rf,-1.D0,2)) then
														call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
														call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
														loggf=-1.
														call setNaN(logg); call setNaN(I_logg)
														call setNaN(g); call setNaN(I_g)
														loggAvail=.false.
														loggAvailAS=.false.
													else !input R is available
														!Only using rho: determine M; re-determine logg (discard the input value) 
														loggf=-1.
														loggAvail=.false.
														loggAvailAS=.false.
														R=R1 !use just the input reliable value
														I_R=I_R1
														call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
														loggAvailCal=.true. !!logg now available thanks to calibration
													end if
												else !logg is best determined
													if (isEq(Rf,-1.D0,2)) then
														call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
														call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
														rhof=-1.
														call setNaN(rho); call setNaN(I_rho)
														call setNaN(logrho); call setNaN(I_logrho)
														rhoAvail=.false.
														rhoAvailAS=.false.
													else !input R is available
														!Only using logg: determine M; re-determine rho (discard the input value) 
														rhof=-1.
														rhoAvail=.false.
														rhoAvailAS=.false.
														R=R1 !use just the input reliable value
														I_R=I_R1
														call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
														rhoAvailCal=.true.
													end if
												end if 
											else
												!!!!!!!!!
												call consistentM(TeB,cumInt,Mvec,M,agreeM)
												if (.not.agreeM) then !new M has been recomputed inside the subroutine in case agreeM is false
																	  !This M will be used in case R is not available directly from input
													if (I_logg0.gt.I_logrho0) then !logrho is best determined
														if (isEq(Rf,-1.D0,2)) then
															loggAvail=.false.
															loggAvailAS=.false.
															!g: mantain original input uncertainty
															call computeLRgfromMrho(M,I_M,rho,I_rho,Teff,I_Teff,L,logL,I_L,I_logL,R,I_R,g,logg)
															loggAvailCal=.true.
														else !input R is available
															!Only using rho: determine M; re-determine logg (discard the input value) 
															loggf=-1.
															loggAvail=.false.
															loggAvailAS=.false.
															R=R1
															I_R=I_R1
															call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
															loggAvailCal=.true. !!logg now available thanks to calibration
														end if
													else
														if (isEq(Rf,-1.D0,2)) then
															rhoAvail=.false.
															rhoAvailAS=.false.
															!rho in g/cm3. Mantain original input uncertainty on rho
															call computeLRrhofromMg(M,I_M,g,I_g,Teff,I_Teff,L,logL,I_L,I_logL,R,I_R,rho,logrho)
															rhoAvailCal=.true.
														else !input R is available
															!Only using logg: determine M; re-determine rho (discard the input value) 
															rhof=-1.
															rhoAvail=.false.
															rhoAvailAS=.false.
															R=R1
															I_R=I_R1
															call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
															rhoAvailCal=.true.
														end if
													end if
													if (isEq(Rf,-1.D0,2)) then
														R=g/10**loggSun*rhoSun/rho
														I_R=R*(I_g/g+I_rho/rho)
														L=R**2*(Teff/TeffSun)**4
														I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
														logL=log10(L)
														I_logL=I_L/L*log10(e)
													end if
												else
													if (I_logg0.lt.I_logrho0) then
														bestLogg=.true.
													else
														bestLogg=.false.
													end if
													if (.not.isEq(Rf,-1.D0,2)) then
													!input R is available: weighted mean to infer M
														Mlogg=g0/10.**loggSun*R**2
														I_Mlogg=Mlogg*(I_g0/g0+2.*I_R/R)
														Mrho=rho0/rhoSun*R**3 !Mo
														I_Mrho=Mrho*(I_rho0/rho0+3.*I_R/R)
														call wMean(Mlogg,I_Mlogg,Mrho,I_Mrho,M,I_M)
														!L already defined since the beginning from R
													end if
												end if
												!!!!!!!!!!!!!
											end if
											deallocate(logTeTr); deallocate(TeB)
										else !Track has been opened, but no logTe compatible with logg were found
											if (I_logg0>I_logrho0) then !logrho is best determined
												if (isEq(Rf,-1.D0,2)) then
													call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
													call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
													loggf=-1.
													call setNaN(logg); call setNaN(I_logg)
													call setNaN(g); call setNaN(I_g)
													loggAvail=.false.
													loggAvailAS=.false.
												else !input R is available
													!Only using rho: determine M; re-determine logg (discard the input value) 
													loggf=-1.
													loggAvail=.false.
													loggAvailAS=.false.
													R=R1 !use just the input reliable value
													I_R=I_R1
													call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
													loggAvailCal=.true. !!logg now available thanks to calibration
												end if
											else !logg is best determined
												if (isEq(Rf,-1.D0,2)) then
													call setNaN(R); call setNaN(I_R); call setNaN(M); call setNaN(I_M)
													call setNaN(L); call setNaN(I_L); call setNaN(logL); call setNaN(I_logL)
													rhof=-1.
													call setNaN(rho); call setNaN(I_rho)
													call setNaN(logrho); call setNaN(I_logrho)
													rhoAvail=.false.
													rhoAvailAS=.false.
												else !input R is available
													!Only using logg: determine M; re-determine rho (discard the input value) 
													rhof=-1.
													rhoAvail=.false.
													rhoAvailAS=.false.
													R=R1 !use just the input reliable value
													I_R=I_R1
													call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
													rhoAvailCal=.true.
												end if
											end if 
										end if
										deallocate(cumInt) 
									else
										if (I_logg0>I_logrho0) then !logrho is best determined
											if (isEq(Rf,-1.D0,2)) then
												call setNan(R); call setNan(I_R); call setNan(M); call setNan(I_M)
												call setNan(L); call setNan(I_L); call setNan(logL); call setNan(I_logL)
												loggf=-1.
												call setNan(logg); call setNan(I_logg)
												call setNan(g); call setNan(I_g)
												loggAvail=.false.
												loggAvailAS=.false.
											else !input R is available
												!Only using rho: determine M; re-determine logg (discard the input value) 
												loggf=-1.
												loggAvail=.false.
												loggAvailAS=.false.
												R=R1 !use just the input reliable value
												I_R=I_R1
												call computeLMg(R,I_R,Teff,I_Teff,rho0,I_rho0,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
												loggAvailCal=.true. !!logg now available thanks to calibration
											end if
										else !logg is best determined
											if (isEq(Rf,-1.D0,2)) then
												call setNan(R); call setNan(I_R); call setNan(M); call setNan(I_M)
												call setNan(L); call setNan(I_L); call setNan(logL); call setNan(I_logL)
												rhof=-1.
												call setNan(rho); call setNan(I_rho)
												call setNan(logrho); call setNan(I_logrho)
												rhoAvail=.false.
												rhoAvailAS=.false.
											else !input R is available
												!Only using logg: determine M; re-determine rho (discard the input value) 
												rhof=-1.
												rhoAvail=.false.
												rhoAvailAS=.false.
												R=R1 !use just the input reliable value
												I_R=I_R1
												call computeLMrho(R,I_R,Teff,I_Teff,g0,I_g0,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
												rhoAvailCal=.true.
											end if
										end if 
									end if 
									!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
								else
									if (.not.isEq(Rf,-1.D0,2)) then
										R=R1
										I_R=I_R1
										if (loggAvail.or.loggAvailAS) then
											call computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
											rhoAvailCal=.true.
										else if (rhoAvail.or.rhoAvailAS) then
											call computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
											loggAvailCal=.true. !!logg available thanks to calibration
										else !noGproxy
											L=R**2*(Teff/TeffSun)**4
											I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
											logL=log10(L)
											I_logL=I_L/L*log10(e)
										end if
									end if
								end if
								!!End compatibility
							end do
							if (cycleW.eq.1) then
								cycle
							end if
						end if
						deallocate(logy_iso)
					end if
					deallocate(ix_sorted)

					!!!Check Activity/vsini. 
					!If inactive stars or slow rotators, consider only isochrones with logt>=logtCutoff 
					check=0
					logtlim=-1
					if (.not.isEq(logRHK,0.D0,2)) then 
						logtMmod=c0-17.912*(logRHK+DlR)-1.6675*(logRHK+DlR)**2 !modified Mamayek
						if (logtMmod<logtCutoff5) then 
							logtlim1=logtMmod
						else
							logtlim1=logtCutoff5
						end if 
					else
						logtlim1=-1
					end if 

					!Searching T.O. of the oldest isochrone!!!!!!!!!!!!!!!!
					allocate(y_iso(size(Isoc,1)))
					caliblogL=calibHRD.or.((loggAvail.or.loggAvailAS.or.loggAvailCal).and.(rhoAvail.or.rhoAvailAS.or. &
							& rhoAvailCal)).or.((calibSPEC.or.calibNoD).and..not.isEq(Rf,-1.D0,2))
					if (caliblogL) then !19/9/2017 !logL avail
						logY=logL
						y_iso=logL_iso
					else if (((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
						& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))) then 
						logY=-logg !to let following inequality logY< refer to lowMS stars
						y_iso=logg_iso
					else
						logY=-logrho !to let following inequality logY< refer to lowMS stars
						y_iso=logrho_iso
					end if 
					call searchTOold(y_iso,logTeff_iso,caliblogL,logt_iso,last_logt,logY_soglia)
		
					logY2=logY
					logY2_soglia=logY_soglia !consider these values later for low MS stars

					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

					if (.not.isEq(vsini,-1.D0,2) .or. .not.isEq(P,-1.D0,2)) then 
						vcheck=.true.
						!Rv: radius to be used in vsini relation 
						if (caliblogL) then 
							Rv=R
						else	!infer a rough estimate for stellar radius to be used for gyrochronology
							if (MminT<MlowMS) then !tracks of lower mass stars are available
								!define logY_sogliaVsini: open M=0.5Mo track; Z=Zstar and search logL/logg(last_logt)
								! that will be logY_sogliaVsini 
								indxM=minloc(abs(MlowMS-MTrAv),1)
								xMTi=Tndxi(indxZt(indxZgini),indxM)
								xMTf=Tndxf(indxZt(indxZgini),indxM)
								allocate(Track05(xMTf-xMTi+1,size(TrackTab,2)))
								Track05=TrackTab(xMTi:xMTf,:)
								MilowMS=MTrAv(indxM)
						
								!Establish maximim luminosity for a ~0.5solarMass star 
								if (caliblogL) then !Never enter this condition. Could only in the following two 
									call InterpLin_M( reshape((/Track05(:,ct_T),Track05(:,clogL_T)/), &
													& (/size(Track05,1),2/)),10.**last_logt,1,1,(/ 2/), &
													& logY_sogliaVsini1,xlow,ylow,xup,yup)
									logY_sogliaVsini=logY_sogliaVsini1(1)
								else if (loggAvail.or.loggAvailAS.or.loggAvailCal) then 
									call InterpLin_M( reshape((/Track05(:,ct_T),Track05(:,clogR_T)/), &
													& (/size(Track05,1),2/)),10.**last_logt,1,1,(/ 2/), &
													& logR_Tr,xlow,ylow,xup,yup) !cm
									R_Tr=10.**logR_Tr(1)
									logY_sogliaVsini=-(log10(MilowMS/(R_Tr/RSun)**2)+loggSun) !minus to let hold
																				! inequality logY>logY_sogliaVsini
								else
									call InterpLin_M( reshape((/Track05(:,ct_T),Track05(:,clogR_T)/), &
													& (/size(Track05,1),2/)),10.**last_logt,1,1,(/ 2/), &
													& logR_Tr,xlow,ylow,xup,yup) !cm
									R_Tr=10.**logR_Tr(1)
									rho_Tr=MilowMS/(R_Tr/RSun)**3*rhoSun
									logY_sogliaVsini=-log10(rho_Tr)
								end if
								deallocate(Track05)
					
								if (logY<logY_sogliaVsini) then !very low MS star. R=R(t) almost constant
									!Select the middle-grid isochrone containing only the lines where M<MilowMS
									!and choose R in correspondence of logY 
									mid_n=nint((logt_halfMS-first_logt)/logtstep+1)
						
									call M2R(logt_iso,mid_n,Isoc,MilowMS,caliblogL,loggAvail.or.loggAvailAS,loggAvailCal,logY,Rv)
								else
									XT=logTeff-4.1
									if (rhoAvail.or.rhoAvailAS.or.rhoAvailCal) then 
										!implementation of my relation R=R(logTeff,logrho,[Fe/H])
										vecT=(/ 1.D0,XT,XT**2,XT**3,logrho,logrho**2,logrho**3,FeH,FeH**2,FeH**3 /)
										Rv=10.**(dot_product(brho,vecT))
									else if (loggAvail.or.loggAvailAS.or.loggAvailCal) then 
										!implementation of my relation R=R(logTeff,logg,[Fe/H])
										vecT=(/ 1.D0,XT,XT**2,XT**3,logg,logg**2,logg**3,FeH,FeH**2,FeH**3 /)
										Rv=10.**(dot_product(blogg,vecT))
									else !Should never enter
										vcheck=.false.
									end if 
								end if 
							else
								XT=logTeff-4.1
								if (rhoAvail.or.rhoAvailAS.or.rhoAvailCal) then 
									!implementation of my relation R=R(logTeff,logrho,[Fe/H])
									vecT=(/ 1.D0,XT,XT**2,XT**3,logrho,logrho**2,logrho**3,FeH,FeH**2,FeH**3 /)
									Rv=10.**(dot_product(brho,vecT))
								else if (loggAvail.or.loggAvailAS.or.loggAvailCal) then 
									!implementation of my relation R=R(logTeff,logg,[Fe/H])
									vecT=(/ 1.D0,XT,XT**2,XT**3,logg,logg**2,logg**3,FeH,FeH**2,FeH**3 /)
									Rv=10.**(dot_product(blogg,vecT))
								else !Should never enter
									vcheck=.false.
								end if 
							end if 
						end if 
						if (vcheck.eqv..true.) then !if vsini or P available should be always true
							if (isEq(vsini,0.D0,2)) then !slow rotating stars for which vsini is reported as zero
								logtlim2=logtCutoff25 !No +0.05. 
								!It's already been taken into account in defining logtCutoff25
							else
								if (.not.isEq(vsini,-1.D0,2).and.(isEq(P,-1.D0,2))) then 
									Omega=4./pi*vsini/(Rv*RSunKm)
									!	logt_Deniss=log10(((OmSun/Omega)^2-1)*(1+A)/(2*(A/2+B))*(tSunLG-tZAMS)+tSunLG) 
									!	logtlim=logt_Deniss; 
									P=2.*pi/Omega/86400 !days
								end if 
								call InterpLin_M(GyroTab,xtau,cBrif,1,(/ ctau /),tau1,xlow,taulow,xup,tauup)
								tau=tau1(1)
								if (.not.isEq(tau,0.D0,2) .and. x>=xlow .and. x<=xup) then 
									logt_Barnes=log10(tau/kC*log(P/P0)+kI/(2.*tau)*(P**2-P0**2))-3 ![Gyr]
									if (10.**logt_Barnes>DtBarnes) then 
										logt_Barnes=log10(10.**logt_Barnes-DtBarnes)+9 ![yrs]
									else
										logt_Barnes=first_logt
									end if 
									if (logt_Barnes<logtCutoff25) then 
										logtlim2=logt_Barnes
									else
										logtlim2=logtCutoff25
									end if 
								else
									logtlim2=-1 !tau=0 => Barnes relation not defined. I don't set any gyro age
								end if 
							end if 
						else
							logtlim2=-1
						end if 
					else
						logtlim2=-1
					end if 
										
					!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
					if ((rhoAvail.or.rhoAvailAS .or. rhoAvailCal)) then !rhof.ne.-1
						if (logY<logY_soglia) then 
							if (.not.caliblogL) then
								!re-set correct values of logg (or logrho) and of logg_soglia (or logrho_soglia)
								logY=-logY
								logY_soglia=-logY_soglia
							end if 
							call trovaIndici(logt_iso,ndxrhoi,ndxrhof) 
							ndxTot=size(ndxrhoi) !=size(ndxrhof)
							jrho=1
							rhot=-100 !initialize like this so that it enters the cycle at least once
							do while (logt_iso(ndxrhoi(jrho))<8 .and. rho-I_rho>rhot .and. jrho<=ndxTot) 
							!Consider only isochrones with t<100Myr. jrho<=ndxTot will always hold
								allocate(logY_isot(ndxrhof(jrho)-ndxrhoi(jrho)+1))
								allocate(rho_isot(ndxrhof(jrho)-ndxrhoi(jrho)+1))
								logY_isot=y_iso(ndxrhoi(jrho):ndxrhof(jrho))
								rho_isot=rho_iso(ndxrhoi(jrho):ndxrhof(jrho))
								ndxrho=minloc(abs(logY-logY_isot),1) 
								!logYt (useless) is the value reported by the isochrone t, that is closer
								! to stellar logY. In correspondence of this value, I'll select the corresponding
								! isochronal stellar density
			!					logYt=logY_isot(ndxrho)
								!rhotmp(jrho)=rho_isot(ndxrho) -->useless to consider it as a vector
					
								rhot=rho_isot(ndxrho)
								jrho=jrho+1
								deallocate(logY_isot); deallocate(rho_isot)
							end do
				
							if (jrho==2) then !While cycle done only once. Don't discard any isochrone
								logtlim3=logt_iso(1) 
								!first age value reported by the isochrones=>don't discard any isochrone!!
							else
								logtlim3=logt_iso(ndxrhoi(jrho-1))+logtstep !added logtstep so that all isochrones with
																	!logt<=logt_iso(ndxrhoi(ii)) are discarded. In fact
																	!isochrones mantained start from logtlim included
								if (logtlim3>logtCutoff5) then 
									logtlim3=logtCutoff5
								end if 
							end if
							deallocate(ndxrhoi); deallocate(ndxrhof) 
						else
							logtlim3=-1
						end if
					else
						logtlim3=-1
					end if
					deallocate(y_iso)
					
					if (.not.isEq(YMg,-100.D0,2).and.FeH_>=-0.2 .and. FeH_<=0.2) then !metall condition
						t_Nissen=(YMg-aYMg)/bYMg*1.e9 ![yrs]
						if (t_Nissen>log10(DtNissen)) then
							logt_Nissen=log10(t_Nissen-DtNissen) !at least this age
						else
							logt_Nissen=first_logt
						end if
						if (logt_Nissen<logtCutoff25) then
							logtlim4=logt_Nissen
						else
							logtlim4=logtCutoff25
						end if
					else
						logtlim4=-1
					end if
										 
					limAges=(/ logtlim1,logtlim2,logtlim3,logtlim4 /)
					klim0=0
					do klim=1,size(limAges) 
						if (.not.isEq(limAges(klim),-1.D0,2)) then 
							klim0=1
							exit
						end if 
					end do
		
					if (klim0==1) then 
						logtlim=maxval(limAges)
						check=maxloc(limAges,1)
					else
						logtlim=-1
					end if
					
					if (.not.isEq(logtlim,-1.D0,2)) then !activity check done
						logtCutoff0=idnint(logtlim*100)
						resto=mod(logtCutoff0,idnint(logtstep*100))
						if (resto>floor(anint(logtstep*100)/2)) then 
							logtCutoff=(logtCutoff0+(dnint(logtstep*100)-resto))/100.
						else
							logtCutoff=(logtCutoff0-resto)/100.
						end if
						if (photIsocAvail) then !anint should be useless
							c_i=13
							allocate(Iso_i(size(t_i),c_i))
							Iso_i=reshape((/ dnint(log10(t_i)*100)/100,t_i,logTeff_i,Teff_i,logL_i,L_i,logg_i,g_i, &
										& M_i,logrho_i,rho_i,BmV_i,BC_i /),(/size(Iso_i,1),c_i/))
							!added rho_i 19/9/2017
						else
							c_i=11
							allocate(Iso_i(size(t_i),c_i))
							Iso_i=reshape((/ dnint(log10(t_i)*100)/100,t_i,logTeff_i,Teff_i,logL_i,L_i,logg_i,g_i, &
										& M_i,logrho_i,rho_i /),(/size(Iso_i,1),c_i/))
							!added rho_i 19/9/2017
						end if
						allocate(Iso_indx(size(Iso_i,1))) 
						call indexx(Iso_i(:,1),Iso_indx) 
						Iso_i=Iso_i(Iso_indx,:)
						deallocate(Iso_indx)
						do iIso=1,size(Iso_i,1) 
							if (isEq(Iso_i(iIso,1),logtCutoff,2)) then 
								iIso8=iIso
								exit
							end if 
						end do
						allocate(Iso_i2(size(Iso_i,1)-iIso8+1,c_i))
						Iso_i2=Iso_i(iIso8:size(Iso_i,1),:)
						deallocate(Iso_i)
						!dealloc _i and then re-alloc with proper dimension (lower than before in case of act check)
						deallocate(t_i); deallocate(logTeff_i); deallocate(Teff_i); deallocate(logL_i)
						deallocate(L_i); deallocate(logg_i); deallocate(g_i); deallocate(M_i)
						deallocate(logrho_i); deallocate(rho_i)
						allocate(t_i(size(Iso_i2,1))); allocate(logTeff_i(size(Iso_i2,1)))
						allocate(Teff_i(size(Iso_i2,1))); allocate(logL_i(size(Iso_i2,1)))
						allocate(L_i(size(Iso_i2,1))); allocate(logg_i(size(Iso_i2,1)))
						allocate(g_i(size(Iso_i2,1))); allocate(M_i(size(Iso_i2,1)))
						allocate(logrho_i(size(Iso_i2,1))); allocate(rho_i(size(Iso_i2,1)))
						t_i=Iso_i2(:,2)
						logTeff_i=Iso_i2(:,3)
						Teff_i=Iso_i2(:,4)
						logL_i=Iso_i2(:,5)
						L_i=Iso_i2(:,6)
						logg_i=Iso_i2(:,7)
						g_i=Iso_i2(:,8)
						M_i=Iso_i2(:,9)
						logrho_i=Iso_i2(:,10) !19/9/2017
						rho_i=Iso_i2(:,11) !19/9/2017
						if (photIsocAvail) then 
							deallocate(BmV_i); deallocate(BC_i)
							allocate(BmV_i(size(Iso_i2,1))); allocate(BC_i(size(Iso_i2,1)))
							BmV_i=Iso_i2(:,12) !19/9/2017
							BC_i=Iso_i2(:,13) !19/9/2017
						end if 
						logt8=logtCutoff
						deallocate(Iso_i2)
					else
						logt8=0
					end if
					!!!!End check activity/vsini
					
					if (IncZ==1) then !load further metallic grids, accounting for I_FeH.
									  !Error bar extension=1 sigma and are "symmetric" in Z plane
						!both flag in multipleZisocPhSCP input parameters are set to 1
						chooseLogg=((loggAvail .or. loggAvailAS) .and. (rhoAvail .or. rhoAvailAS) .and. bestLogg) &
									& .or. ((loggAvail.or.loggAvailAS) .and. .not.(rhoAvail.or.rhoAvailAS))
						if (I_FeH_>0) then !avoid that flags are considered as actual uncertainties
							if (calibHRD) then 
								xx=BmV
								yy=Vass
							end if 
							if (calibNoD) then 
								xx=BmV
								if (caliblogL) then
									yy=logL
								else if (chooseLogg) then 
									yy=logg
								else
									yy=logrho
								end if 
							end if 
							if (calibSPEC) then 
								xx=logTeff
								if (caliblogL) then
									yy=logL
								else if (chooseLogg) then 
									yy=logg
								else
									yy=logrho
								end if 
							end if
							call multipleZisocPhSCP( FeH, I_FeH_, Zvec, model, xx, yy, &
								& logt8, percorso, 1.D0, 1, calibHRD, calibNoD, calibSPEC, chooseLogg, &
								& caliblogL,hstar, hstarlim, Isoc_i, Zlow, Zup,useColor,idCol)
							if (allocated(Isoc_i)) then !isochrone grids (WITH the reference one INCLUDED) are loaded
								deallocate(t_i); deallocate(logTeff_i); deallocate(Teff_i)
								deallocate(logL_i); deallocate(L_i); deallocate(logg_i)
								deallocate(g_i); deallocate(M_i); deallocate(logrho_i)
								deallocate(rho_i)
								allocate(t_i(size(Isoc_i,1))); allocate(logTeff_i(size(Isoc_i,1)))
								allocate(Teff_i(size(Isoc_i,1))); allocate(logL_i(size(Isoc_i,1)))
								allocate(L_i(size(Isoc_i,1))); allocate(logg_i(size(Isoc_i,1)))
								allocate(g_i(size(Isoc_i,1))); allocate(M_i(size(Isoc_i,1)))
								allocate(logrho_i(size(Isoc_i,1))); allocate(rho_i(size(Isoc_i,1)))
								t_i=Isoc_i(:,2)
								logTeff_i=Isoc_i(:,3)
								Teff_i=Isoc_i(:,4)
								logL_i=Isoc_i(:,5)
								L_i=Isoc_i(:,6)
								logg_i=Isoc_i(:,7)
								g_i=Isoc_i(:,8)
								M_i=Isoc_i(:,9)
								logrho_i=Isoc_i(:,10) !19/9/2017
								rho_i=Isoc_i(:,11) !19/9/2017
								if (photIsocAvail) then
									deallocate(BmV_i); deallocate(BC_i)
									allocate(BmV_i(size(Isoc_i,1))); allocate(BC_i(size(Isoc_i,1)))
									BmV_i=Isoc_i(:,12) !19/9/2017
									BC_i=Isoc_i(:,13) !19/9/2017
								end if
								deallocate(Isoc_i)
							end if 
						end if 
					end if 

					!!!only theoretical data _i
					allocate(Gauss(size(logTeff_i)))
					allocate(w(size(logTeff_i)))
					if (gProxyAvail) then 
						if (calibHRD .or. ((loggAvail.or.loggAvailAS.or.loggAvailCal) .and. &
							& (rhoAvail.or.rhoAvailAS.or.rhoAvailCal)).or.((calibSPEC.or.calibNoD).and..not.isEq(Rf,-1.D0,2))) then
							!I can recover Teff, L, M, logg
							if (logY2<logY2_soglia) then 
								!!!Find ZAMS in the Track:
								! logL starts increasing after the decreasing along the Hayashi line 
								call setThreshold(caliblogL,(((loggAvail.or.loggAvailAS).and.(rhoAvail.or.rhoAvailAS).and.bestLogg) &
									& .or.((loggAvail.or.loggAvailAS).and..not.(rhoAvail.or.rhoAvailAS))),logY2,logY2_soglia, &
									& logY0,logY0_soglia,cZAMS,cyT)
								allocate(yT_iso(size(Isoc,1)))
								if (cyT.eq.0) then !happens for logrho that hasn't a direct column in the isoch grid
									yT_iso=logrho_iso
								else
									yT_iso=Isoc(:,cyT)
								end if

								call ifComputeIsocAge(ZAMSz,cZAMS,logY0,logt_iso,logTeff_iso,yT_iso, &
													& I_logTeff,ageIsoc)
										
								deallocate(yT_iso)
							else
								ageIsoc=1
							end if 
				
							indxM=minloc(abs(M-MTrAv),1)
							xMTi=Tndxi(indxZt(indxZgini),indxM)
							xMTf=Tndxf(indxZt(indxZgini),indxM)
							allocate(Tracks(xMTf-xMTi+1,size(TrackTab,2)))
							Tracks=TrackTab(xMTi:xMTf,:)
	
							allocate(t_Tracks(size(Tracks,1)));allocate(logL_Tracks(size(Tracks,1)))
							allocate(logTeff_Tracks(size(Tracks,1)));allocate(vEvo(size(t_i)))
							t_Tracks=Tracks(:,ct_T)
							logL_Tracks=Tracks(:,clogL_T)
							logTeff_Tracks=Tracks(:,clogTe_T)
							deallocate(Tracks)
							call computeVevo(t_Tracks,logL_Tracks,logTeff_Tracks,t_i,varTrlim,vEvo)
							deallocate(t_Tracks);deallocate(logL_Tracks);deallocate(logTeff_Tracks)
				
							cRif=cVelL
							call computeVrif(TabVel,cRif,M,vRif)
								
							Gauss=1./(2.*pi*I_logTeff*I_logL)*e**(-0.5*((logTeff-logTeff_i)/I_logTeff)**2)* &
									& e**(-0.5*((logL-logL_i)/I_logL)**2)
							w=1./(((logL-logL_i)/I_logL)**2+((logTeff-logTeff_i)/I_logTeff)**2+ &
								& ((M-M_i)/I_M)**2+((logg-logg_i)/I_logg)**2+(log10(vRif/vEvo))**2)
							deallocate(vEvo)
						else
							if (logY2<logY2_soglia) then 
								call setThreshold(caliblogL,loggAvail.or.loggAvailAS,logY2,logY2_soglia, &
									& logY0,logY0_soglia,cZAMS,cyT)
								allocate(yT_iso(size(Isoc,1)))
								if (cyT.eq.0) then !happens for logrho that hasn't a direct column in the isoch grid
									yT_iso=logrho_iso
								else
									yT_iso=Isoc(:,cyT)
								end if
					
								call ifComputeIsocAge(ZAMSz,cZAMS,logY0,logt_iso,logTeff_iso,yT_iso,I_logTeff, &
													& ageIsoc)
								deallocate(yT_iso)
							else
								ageIsoc=1
							end if 
							
							allocate(y_i(size(logg_i)))
							if (loggAvail.or.loggAvailAS) then !only logg available
								y=logg
								I_y=I_logg
								y_i=logg_i
								cRif=cVelg
							else !only logrho available
								y=logrho
								I_y=I_logrho
								y_i=logrho_i
								cRif=cVelrho
							end if 
							!!!!!!!!!!!!!!!! 
							!!!!!!!!!!!!!!!!
							Gauss=1./(2.*pi*I_logTeff*I_y)*e**(-0.5*((logTeff-logTeff_i)/I_logTeff)**2)* &
							& e**(-0.5*((y-y_i)/I_y)**2)
							w=1./(((logTeff-logTeff_i)/I_logTeff)**2+((y-y_i)/I_y)**2)
							Spesi=sum(w*Gauss)
							if (Spesi<=1.D-200) then 
								cycle
							end if
				
							M_starTr=dot_product(w,M_i*Gauss)/Spesi
				
							indxM=minloc(abs(M_starTr-MTrAv),1)
							xMTi=Tndxi(indxZt(indxZgini),indxM)
							xMTf=Tndxf(indxZt(indxZgini),indxM)
							allocate(Tracks(xMTf-xMTi+1,size(TrackTab,2)))
							Tracks=TrackTab(xMTi:xMTf,:)
	
							allocate(t_Tracks(size(Tracks,1))); allocate(logTeff_Tracks(size(Tracks,1)))
							allocate(R_Tracks(size(Tracks,1))); allocate(logY_Tracks(size(Tracks,1)))
							allocate(vEvo(size(t_i)))
							t_Tracks=Tracks(:,ct_T)
							logTeff_Tracks=Tracks(:,clogTe_T)
							R_Tracks=10.**Tracks(:,clogR_T) !cm
							M_Tracks=Tracks(1,cM_T)
							deallocate(Tracks)
							if (loggAvail.or.loggAvailAS) then !only logg available
								logY_Tracks=log10(M_Tracks/(R_Tracks/RSun)**2)+loggSun
							else !only logrho available
								logY_Tracks=log10(M_Tracks/(R_Tracks/RSun)**3)+log10(rhoSun)
							end if
							call computeVevo(t_Tracks,logY_Tracks,logTeff_Tracks,t_i,varTrlim,vEvo)
							deallocate(t_Tracks); deallocate(logTeff_Tracks); deallocate(R_Tracks)
							deallocate(logY_Tracks)
				
							call computeVrif(TabVel,cRif,M_starTr,vRif)
				
							w=1./(((logTeff-logTeff_i)/I_logTeff)**2+((y-y_i)/I_y)**2+(log10(vRif/vEvo))**2)
							deallocate(vEvo)
							Spesi=sum(w*Gauss)
							if (Spesi<=1.D-200) then 
								cycle
							end if 
							M_starTr2=dot_product(w,M_i*Gauss)/Spesi
							kTr=1
							do while (abs((M_starTr2-M_starTr)/M_starTr)>DMTr .and. kTr<10) 
								if (logY2<logY2_soglia) then 
									call setThreshold(caliblogL,loggAvail.or.loggAvailAS,logY2,logY2_soglia, &
										& logY0,logY0_soglia,cZAMS,cyT)
									allocate(yT_iso(size(Isoc,1)))
									if (cyT.eq.0) then 
									!happens for logrho that hasn't a direct column in the isoch grid
										yT_iso=logrho_iso
									else
										yT_iso=Isoc(:,cyT)
									end if
						
									call ifComputeIsocAge(ZAMSz,cZAMS,logY0,logt_iso,logTeff_iso,yT_iso, &
														& I_logTeff,ageIsoc)
										
									deallocate(yT_iso)
								else
									ageIsoc=1
								end if 
					
								M_starTr=M_starTr2
								indxM=minloc(abs(M_starTr-MTrAv),1)
								xMTi=Tndxi(indxZt(indxZgini),indxM)
								xMTf=Tndxf(indxZt(indxZgini),indxM)
								allocate(Tracks(xMTf-xMTi+1,size(TrackTab,2)))
								Tracks=TrackTab(xMTi:xMTf,:)
	
								allocate(t_Tracks(size(Tracks,1))); allocate(logTeff_Tracks(size(Tracks,1)))
								allocate(R_Tracks(size(Tracks,1))); allocate(logY_Tracks(size(Tracks,1)));
								allocate(vEvo(size(t_i)))
								t_Tracks=Tracks(:,ct_T)
								logTeff_Tracks=Tracks(:,clogTe_T)
								R_Tracks=10.**Tracks(:,clogR_T) !cm
								M_Tracks=Tracks(1,cM_T)
								deallocate(Tracks)
								if (loggAvail.or.loggAvailAS) then !only logg available
									logY_Tracks=log10(M_Tracks/(R_Tracks/RSun)**2)+loggSun
								else !only logrho available
									logY_Tracks=log10(M_Tracks/(R_Tracks/RSun)**3)+log10(rhoSun)
								end if
								call computeVevo(t_tracks,logY_Tracks,logTeff_Tracks,t_i,varTrlim,Vevo)
								deallocate(t_Tracks); deallocate(logTeff_Tracks); deallocate(R_Tracks)
								deallocate(logY_Tracks)
					
								call computeVrif(TabVel,cRif,M_starTr,vRif)
					
								w=1./(((logTeff-logTeff_i)/I_logTeff)**2+((y-y_i)/I_y)**2+(log10(vRif/vEvo))**2)
								deallocate(vEvo)
								Spesi=sum(w*Gauss)
								if (Spesi<=1.D-200) then 
									cycle
								end if 
								M_starTr2=dot_product(w,M_i*Gauss)/Spesi
								kTr=kTr+1
							end do 
							!!!!!!!!!!!!!!!!!!! 
							!!!!!!!!!!!!!!!!!!! 
						end if 
					else
						if (logY2<logY2_soglia) then 
							call setThreshold(caliblogL,loggAvail.or.loggAvailAS,logY2,logY2_soglia, &
									& logY0,logY0_soglia,cZAMS,cyT)
							allocate(yT_iso(size(Isoc,1)))
							if (cyT.eq.0) then !happens for logrho that hasn't a direct column in the isoch grid
								yT_iso=logrho_iso
							else
								yT_iso=Isoc(:,cyT)
							end if
							 
							call ifComputeIsocAge(ZAMSz,cZAMS,logY0,logt_iso,logTeff_iso,yT_iso,I_logTeff,ageIsoc)
									
							deallocate(yT_iso)
						else
							ageIsoc=1
						end if 
						
						Gauss=1./(2.*pi*I_logTeff*I_logL)*e**(-0.5*((logTeff-logTeff_i)/I_logTeff)**2)* &
								& e**(-0.5*((logL-logL_i)/I_logL)**2)
						w=1./(((logL-logL_i)/I_logL)**2+((logTeff-logTeff_i)/I_logTeff)**2)
						Spesi=sum(w*Gauss)
						if (Spesi<=1.D-200) then 
							cycle
						end if 
						M_starTr=dot_product(w,M_i*Gauss)/Spesi

						indxM=minloc(abs(M_starTr-MTrAv),1)
						xMTi=Tndxi(indxZt(indxZgini),indxM)
						xMTf=Tndxf(indxZt(indxZgini),indxM)
						allocate(Tracks(xMTf-xMTi+1,size(TrackTab,2)))
						Tracks=TrackTab(xMTi:xMTf,:)
						allocate(t_Tracks(size(Tracks,1)));allocate(logL_Tracks(size(Tracks,1)))
						allocate(logTeff_Tracks(size(Tracks,1)));allocate(vEvo(size(t_i)))
	
						t_Tracks=Tracks(:,ct_T)
						logL_Tracks=Tracks(:,clogL_T)
						logTeff_Tracks=Tracks(:,clogTe_T)
						deallocate(Tracks)
						call computeVevo(t_Tracks,logL_Tracks,logTeff_Tracks,t_i,varTrlim,vEvo)
						deallocate(t_Tracks);deallocate(logL_Tracks);deallocate(logTeff_Tracks) 
			
						cRif=cVelL
						call computeVrif(TabVel,cRif,M_starTr,vRif)
						
						w=1./(((logL-logL_i)/I_logL)**2+((logTeff-logTeff_i)/I_logTeff)**2+(log10(vRif/vEvo))**2)
						deallocate(vEvo)
						Spesi=sum(w*Gauss)
						if (Spesi<=1.D-200) then 
							cycle
						end if 
						M_starTr2=dot_product(w,M_i*Gauss)/Spesi
						kTr=1
						do while (abs((M_starTr2-M_starTr)/M_starTr)>DMTr .and. kTr<10) 
							if (logY2<logY2_soglia) then 
								call setThreshold(caliblogL,loggAvail.or.loggAvailAS,logY2,logY2_soglia, &
									& logY0,logY0_soglia,cZAMS,cyT)
								allocate(yT_iso(size(Isoc,1)))
								if (cyT.eq.0) then !happens for logrho that hasn't a direct column in the isoch grid
									yT_iso=logrho_iso
								else
									yT_iso=Isoc(:,cyT)
								end if
								 
								call ifComputeIsocAge(ZAMSz,cZAMS,logY0,logt_iso,logTeff_iso,yT_iso,I_logTeff, &
													& ageIsoc)
										
								deallocate(yT_iso)
							else
								ageIsoc=1
							end if 
				
							M_starTr=M_starTr2
							indxM=minloc(abs(M_starTr-MTrAv),1)
							xMTi=Tndxi(indxZt(indxZgini),indxM)
							xMTf=Tndxf(indxZt(indxZgini),indxM)
							allocate(Tracks(xMTf-xMTi+1,size(TrackTab,2)))
							Tracks=TrackTab(xMTi:xMTf,:)
							allocate(t_Tracks(size(Tracks,1)));allocate(logL_Tracks(size(Tracks,1)))
							allocate(logTeff_Tracks(size(Tracks,1)));allocate(vEvo(size(t_i)))
	
							t_Tracks=Tracks(:,ct_T)
							logL_Tracks=Tracks(:,clogL_T)
							logTeff_Tracks=Tracks(:,clogTe_T)
							deallocate(Tracks)
							call computeVevo(t_Tracks,logL_Tracks,logTeff_Tracks,t_i,varTrlim,vEvo)
							deallocate(t_Tracks);deallocate(logL_Tracks);deallocate(logTeff_Tracks) 
			
							cRif=cVelL
							call computeVrif(TabVel,cRif,M_starTr,vRif)
							w=1./(((logL-logL_i)/I_logL)**2+((logTeff-logTeff_i)/I_logTeff)**2+(log10(vRif/vEvo))**2)
							deallocate(vEvo)
							Spesi=sum(w*Gauss)
							if (Spesi<=1.D-200) then 
								cycle
							end if 
							M_starTr2=dot_product(w,M_i*Gauss)/Spesi
							kTr=kTr+1
						end do 
					end if 
					Spesi=sum(w*Gauss)
					if (Spesi<=1.D-200) then 
						cycle
					end if
					if (cycleW.eq.1) then
						cycle
					end if 
					t_star2=dot_product(w,t_i*Gauss)/Spesi !filter the selected isochrone through a Gaussian
					!and then weight it through w. It's like as attributing it the weight w*Gauss
					Teff_star2=dot_product(w,Teff_i*Gauss)/Spesi
					L_star2=dot_product(w,L_i*Gauss)/Spesi
					M_star2=dot_product(w,M_i*Gauss)/Spesi
					g_star2=dot_product(w,g_i*Gauss)/Spesi
					rho_star2=dot_product(w,rho_i*Gauss)/Spesi
					!!! 
					if (photIsocAvail) then 
						BmV_star2=dot_product(w,BmV_i*Gauss)/Spesi
						BC_star2=dot_product(w,BC_i*Gauss)/Spesi
					end if 
					!!! 
					I_t_star2=sqrt(dot_product(w*Gauss,(t_i-t_star2)**2)/Spesi)
					I_Teff_star2=sqrt(dot_product(w*Gauss,(Teff_i-Teff_star2)**2)/Spesi)
					I_L_star2=sqrt(dot_product(w*Gauss,(L_i-L_star2)**2)/Spesi)
					I_M_star2=sqrt(dot_product(w*Gauss,(M_i-M_star2)**2)/Spesi)
					I_g_star2=sqrt(dot_product(w*Gauss,(g_i-g_star2)**2)/Spesi)
					I_rho_star2=sqrt(dot_product(w*Gauss,(rho_i-rho_star2)**2)/Spesi)
					!!! 
					if (photIsocAvail) then 
						I_BmV_star2=sqrt(dot_product(w*Gauss,(BmV_i-BmV_star2)**2)/Spesi)
						I_BC_star2=sqrt(dot_product(w*Gauss,(BC_i-BC_star2)**2)/Spesi)
					end if 

					logg_star2=log10(g_star2)
					R_star2=sqrt(L_star2/(Teff_star2/TeffSun)**4)
					I_logg_star2=log10(e)*I_g_star2/g_star2
					I_R_star2=R_star2*(0.5*I_L_star2/L_star2+2.*I_Teff_star2/Teff_star2)
					kwhile=kwhile+1
				end do
				t_star=t_star2
				Teff_star=Teff_star2
				L_star=L_star2
				M_star=M_star2
				g_star=g_star2
				rho_star=rho_star2
				if (photIsocAvail) then 
					BmV_star=BmV_star2
					BC_star=BC_star2
				end if 
				I_t_star=I_t_star2
				I_Teff_star=I_Teff_star2
				I_L_star=I_L_star2
				I_M_star=I_M_star2
				I_g_star=I_g_star2
				I_rho_star=I_rho_star2
				if (photIsocAvail) then 
					I_BmV_star=I_BmV_star2
					I_BC_star=I_BC_star2
				end if 
				logg_star=logg_star2
				R_star=R_star2
				I_logg_star=I_logg_star2
				I_R_star=I_R_star2
			end if !endif element diffusion: Z2!=Z. if true, element diffusion is performed
		end if !endif: test if Z changes because of element diffusion: M_star>=MinfED & M_star<=MsupED
		if (Spesi<=1.D-200) then 
			cycle
		end if
		if (cycleW.eq.1) then
			cycle
		end if
		if (ageIsoc==0) then !do not consider stellar age from isochrone
			if (logtlim.ne.-1) then  !activity check done
				t_star=10.**logt8
				I_t_star=0 !flag to say that the stellar age is a lower limit
			else
				t_star=logt8 !=0 since activity checks hasn't been performed
				I_t_star=logt8 !=0
			end if 
		end if
						
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
		!!!    Posterior check if output mass and radius fall within the input interval of mass and radius    !!!
		!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
		Mavail=(calibHRD.and.(loggAvail.or.loggAvailAS.or.rhoAvail.or.rhoAvailAS)).or. &
				& ((loggAvail.or.loggAvailAS.or.loggAvailCal).and.(rhoAvail.or.rhoAvailAS.or.rhoAvailCal)) &
				& .or.(calibSPEC.and.gProxyAvail.and..not.isEq(Rf,-1.D0,2)) & 
				& .or.(calibNoD.and.gProxyAvail.and..not.isEq(Rf,-1.D0,2))
		if (caliblogL) then !R is defined
			if (I_R/R.lt.0.05) then
				checkM=1
				row=row+1
			else
				if (Mavail) then
					if ((abs(M_star-M)>I_M .and. abs(R_star-R)>I_R)) then 
					!in both cases no consistency between input/output
						checkM=0
						rowN=rowN+1
						M_starNC=M_star
						I_M_starNC=I_M_star
						R_starNC=R_star
						I_R_starNC=I_R_star
					else !output mass or radius inside the error box of the input value:
						 !consistency between input and output at least for one parameter
						checkM=1
						row=row+1
					end if 
				else !I have only R and not M
					if (abs(R_star-R)>I_R) then 
						checkM=0
						rowN=rowN+1
						M_starNC=M_star
						I_M_starNC=I_M_star
						R_starNC=R_star
						I_R_starNC=I_R_star
					else
						checkM=1
						row=row+1
					end if 
				end if
			end if
		else !if I don't have R, I have neither M
			checkM=1 !accept the star because there's no way to check it
			row=row+1
		end if
		
		if (link.eq.1) then !just the 1st of the chain
			TLRM=" Teff I_Teff L I_L R I_R M I_M"
			MetCut=" Zact Zini Yact Yini"
			titolo=trim(Intestaz)//trim(TLRM)//" logg I_logg rho I_rho"
			titoloN=titolo
			if (gProxyAvail) then 
				if (.not.calibSPEC) then 
					if ((loggAvail.or.loggAvailAS) .and. (rhoAvail.or.rhoAvailAS)) then
					!necesarily consistent, otherwise I'd have set one of them to false
					!Thus Avail means that the original values have been used
						if (checkM==1) then 
							SCPInputtmp=(/SCP, Teff, I_Teff, L, I_L, R, I_R, M, I_M, -1.*ones(4)/)
						else
							SCPInputNtmp=(/SCP, Teff, I_Teff, L, I_L, R, I_R, M, I_M, -1.*ones(4)/)
						end if 
					else if ((loggAvail.or.loggAvailAS) .and. calibHRD) then
					!only logg used as original value in calibHRD 
						if (checkM==1) then 
							SCPInputtmp=(/SCP, Teff, I_Teff, L, I_L, R, I_R, M, I_M, -1.*ones(2), &
												& rho, I_rho/)
						else
							SCPInputNtmp=(/SCP, Teff, I_Teff, L, I_L, R, I_R, M, I_M, -1.*ones(2), &
												& rho, I_rho/)
						end if 
					else if (loggAvail.or.loggAvailAS) then !only logg in calibNoD
						if (checkM==1) then
							if (rhoAvailCal) then
							!there was also the input rho, but it's been considered not fully consistent (and then recomputed)
								SCPInputtmp=(/SCP, Teff, I_Teff, L, I_L, R, I_R, M, I_M, -1.*ones(2), rho, I_rho/)
							else 
								SCPInputtmp=(/SCP, Teff, I_Teff, -1.*ones(10)/)
							end if
						else
							if (rhoAvailCal) then
								SCPInputNtmp=(/SCP, Teff, I_Teff, L, I_L, R, I_R, M, I_M, -1.*ones(2), rho, I_rho/)
							else
								SCPInputNtmp=(/SCP, Teff, I_Teff, -1.*ones(10)/)
							end if
						end if 
					else !only rhoAvail --> only original rho used
						if (calibHRD) then 
							if (checkM==1) then 
								SCPInputtmp=(/SCP, Teff, I_Teff, L, I_L, R, I_R, M, I_M, logg, I_logg, &
												& -1.*ones(2)/)
							else
								SCPInputNtmp=(/SCP, Teff, I_Teff, L, I_L, R, I_R, M, I_M, logg, I_logg, &
												& -1.*ones(2)/)
							end if 
						else !only rho in calibNoD
							if (checkM==1) then
								if (loggAvailCal) then
								!there was also the input logg, but it's been considered not fully consistent (and then recomputed)
									SCPInputtmp=(/SCP, Teff, I_Teff, L, I_L, R, I_R, M, I_M, logg, I_logg, -1.*ones(2)/)
								else !logg never been available
									SCPInputtmp=(/SCP, Teff, I_Teff, -1.*ones(10)/)
								end if
							else
								if (loggAvailCal) then
									SCPInputNtmp=(/SCP, Teff, I_Teff, L, I_L, R, I_R, M, I_M, logg, I_logg, -1.*ones(2)/)
								else
									SCPInputNtmp=(/SCP, Teff, I_Teff, -1.*ones(10)/)
								end if
							end if 
						end if 
					end if 
				else
					if ((loggAvail.or.loggAvailAS) .and. (rhoAvail.or.rhoAvailAS)) then
					!both used as original values 
						if (checkM==1) then 
							SCPInputtmp=(/SCP, -1.*ones(2), L, I_L, R, I_R, M, I_M, -1.*ones(4)/)
						else
							SCPInputNtmp=(/SCP, -1.*ones(2), L, I_L, R, I_R, M, I_M, -1.*ones(4)/)
						end if
					!In the other else: only rho xor logg Avail, i.e. just one of them used as original value.
					!The other could be either present but not consistent, or not present
					else if (loggAvail.or.loggAvailAS) then
						if (checkM==1) then
							if (rhoAvailCal) then
								SCPInputtmp=(/SCP, -1.*ones(2), L, I_L, R, I_R, M, I_M, -1.*ones(2), rho, I_rho/)
							else
								SCPInputtmp=(/SCP, -1.*ones(12)/)
							end if
						else
							if (rhoAvailCal) then
								SCPInputNtmp=(/SCP, -1.*ones(2), L, I_L, R, I_R, M, I_M, -1.*ones(2), rho, I_rho/)
							else
								SCPInputNtmp=(/SCP, -1.*ones(12)/)
							end if
						end if 
					else !only rho used as original input
						if (checkM==1) then
							if (loggAvailCal) then
								SCPInputtmp=(/SCP, -1.*ones(2), L, I_L, R, I_R, M, I_M, logg, I_logg, -1.*ones(2)/)
							else 
								SCPInputtmp=(/SCP, -1.*ones(12)/)
							end if
						else
							if (loggAvailCal) then
								SCPInputNtmp=(/SCP, -1.*ones(2), L, I_L, R, I_R, M, I_M, logg, I_logg, -1.*ones(2)/)
							else
								SCPInputNtmp=(/SCP, -1.*ones(12)/)
							end if
						end if 
					end if 
				end if 
			else
				if (checkM==1) then 
					SCPInputtmp=(/SCP, Teff, I_Teff, L, I_L, R, I_R, -1.*ones(6)/)
				else
					SCPInputNtmp=(/SCP, Teff, I_Teff, L, I_L, R, I_R, -1.*ones(6)/)
				end if 
			end if
		

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			!!!!!!!!!!!!		Writing gProxyFlag		!!!!!!!!!!!!!!!!!! 
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			if ((loggAvail.or.loggAvailAS) .and. (rhoAvail.or.rhoAvailAS)) then 
				gFlag=1
			else if ((loggAvail.or.loggAvailAS) .and. rhoAvailCal) then 
				if (isEq(rhoI,-1.D0,2)) then !rho was not available as input parameter. It's been computed
					gFlag=102
				else !rho was available, but was not judged consistent with logg. It's been recomputed then
					gFlag=122
				end if 
			else if ((rhoAvail.or.rhoAvailAS) .and. loggAvailCal) then 
				if (isEq(loggI,-1.D0,2)) then !logg was not available as input parameter. It's been computed
					gFlag=201
				else !logg was available, but was not judged consistent with rho. It's been recomputed then
					gFlag=211
				end if 
			else if ((loggAvail.or.loggAvailAS) .and. .not.rhoAvail .and. .not.rhoAvailAS .and. .not.rhoAvailCal) then !only logg used
				if (isEq(rhoI,-1.D0,2)) then !rho was not available as input parameter, nor it was possible to compute it
					gFlag=100
				else
					gFlag=120 !rho was available as input parameter, but judged inconsistent with logg
							  ! and was not possible to recompute it
				end if 
			else if ((rhoAvail.or.rhoAvailAS) .and. .not.loggAvail .and. .not.loggAvailAS .and. .not.loggAvailCal) then !only rho used
				if (isEq(loggI,-1.D0,2)) then !logg was not available as input parameter, nor it was possible to compute it
					gFlag=200
				else
					gFlag=210 !logg was available as input parameter, but judged inconsistent with rho
							  ! and was not possible to recompute it
				end if 
			else !no gProxy available
				gFlag=0
			end if 

			if (photIsocAvail) then 
				if (checkM==1) then 
					Stelle=(/ SCP(1),t_star,I_t_star,Teff_star,I_Teff_star,L_star,I_L_star,M_star,I_M_star, &
									& logg_star,I_logg_star,R_star,I_R_star,rho_star,I_rho_star,BmV_star,I_BmV_star, &
									& dble(check),logt8,dble(Bin),dble(gFlag),Zact,Z,Yact,Yini /) !, BC_star, I_BC_star
				else
					StelleN=(/ SCP(1),t_star,I_t_star,Teff_star,I_Teff_star,L_star,I_L_star,M_starNC,I_M_starNC, &
									& logg_star,I_logg_star,R_starNC,I_R_starNC,rho_star,I_rho_star,BmV_star,I_BmV_star, &
									& dble(check),logt8,dble(Bin),dble(gFlag),Zact,Z,Yact,Yini /)
				end if 
			else
				if (checkM==1) then 
					Stelle=(/ SCP(1),t_star,I_t_star,Teff_star,I_Teff_star,L_star,I_L_star,M_star,I_M_star, &
									& logg_star,I_logg_star,R_star,I_R_star,rho_star,I_rho_star,dble(-100),dble(-100), &
									& dble(check),logt8,dble(Bin),dble(gFlag),Zact,Z,Yact,Yini /)
									!!BmV_star and I_BmV_star not recoverable if not present in the isochrones
				else
					StelleN=(/ SCP(1),t_star,I_t_star,Teff_star,I_Teff_star,L_star,I_L_star,M_starNC,I_M_starNC, &
									& logg_star,I_logg_star,R_starNC,I_R_starNC,rho_star,I_rho_star,dble(-100),dble(-100), &
									& dble(check),logt8,dble(Bin),dble(gFlag),Zact,Z,Yact,Yini /)
				end if 
			end if 
	
			if (checkM==1) then 
				write(3,*) StarName(1)%fName
			else
				write(4,*) StarName(1)%fName
			end if
		end if
	end do
	if (row.eq.0.and.rowN.eq.0) then !Process not completed
		acc=0
	end if
	MetCut=" Zact Zini Yact Yini"
		
	if (mod(link,10000).eq.0) then
		print*,'link',link
	endif
	if (link.eq.1) then !just the 1st of the chain to check what's happening
		CLOSE(3)
		CLOSE(4)
		
		if (photIsocAvail) then 
		titoloAge="#ROW Age Inc_Age Teff Inc_Teff L Inc_L M Inc_M logg Inc_logg R Inc_R rho Inc_rho BmV Inc_BmV "// &
					& "Act_check AgeCutoff Binary gFlag"//trim(MetCut)
		else
			titoloAge="#ROW Age Inc_Age Teff Inc_Teff L Inc_L M Inc_M logg Inc_logg R Inc_R rho Inc_rho "// &
					& "Act_check AgeCutoff Binary gFlag"//trim(MetCut)
		end if 
		titoloAgeN=titoloAge

		if (row>0) then
			outI="Input"//FileTxt//trim(Est)
			outAge="Age"//FileTxt//trim(Est)

			!Check whether SCPInput has an all column made of -1 (then cutoff it with the title as well) 
			call removeColumn(SCPInputtmp,titolo,-1.D0,SCPInput)
		
			OPEN(UNIT=3, FILE=Isoch//trim(outI))
			write(3,*) trim(titolo)
			write(3,*) SCPInput
			CLOSE(3)
		
			OPEN(UNIT=3, FILE=Isoch//trim(outAge))
			write(3,*) trim(titoloAge)
			write(3,*) Stelle
			CLOSE(3)
		end if 
		if (rowN>0) then
			outIN="InputDISC"//FileTxt//trim(Est)
			outAgeN="AgeDISC"//FileTxt//trim(Est)

			!!!discarded ones 
			call removeColumn(SCPInputNtmp,titoloN,-1.D0,SCPInputN)
			OPEN(UNIT=3, FILE=Isoch//trim(outIN))
			write(3,*) trim(titoloN)
			write(3,*) SCPInputN
			CLOSE(3)

			!!!discarded ones 
			OPEN(UNIT=3, FILE=Isoch//trim(outAgeN))
			write(3,*) trim(titoloAgeN)
			write(3,*) StelleN
			CLOSE(3)
		end if
	end if
	
END SUBROUTINE SCPmcmcPD12S

function ones(n)
	IMPLICIT NONE
	INTEGER n
	DOUBLE PRECISION, DIMENSION(n) :: ones
	
	INTEGER i
	
	ones=(/ (1.D0, i=1,n) /)
end function ones

function strncmp(str1,str2,n)
	IMPLICIT NONE
	CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: str1,str2
	INTEGER n
	LOGICAL strncmp
	
	INTEGER l1,l2
	CHARACTER(LEN=:), ALLOCATABLE :: str1n,str2n
		
	interface
		function strcmp(str1,str2)
			IMPLICIT NONE
			CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: str1,str2
			LOGICAL strcmp
		end function strcmp
	end interface
	
	l1=LEN_TRIM(str1)
	l2=LEN_TRIM(str2)
	if ((n.gt.l1).or.(n.gt.l2)) then
		print*,'from strncmp: Specified n exceeds string length'
		strncmp=.false.
	else
		str1n=str1(1:n)
		str2n=str2(1:n)
		strncmp=strcmp(str1n,str2n)
	end if

	return
end function strncmp

function strcmp(str1,str2)
	!Case sensitive string comparison
	!True if str1 and str2 are equal (blanks at the beginning or at the end are ignored)
	IMPLICIT NONE
	CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: str1,str2
	LOGICAL strcmp
	
	INTEGER i,l,ug
	LOGICAL chCmp
	
	if (LEN_TRIM(str1).eq.LEN_TRIM(str2)) then
		ug=1
		l=LEN_TRIM(str1)
		do i=1,l
			chCmp=(str1(i:i).eq.str2(i:i))
			if (chCmp.eqv..false.) then
				ug=0
				exit
			end if
		end do
		if (ug.eq.1) then
			strcmp=.true.
		else
			strcmp=.false.
		end if
	else
		strcmp=.false.
	end if
	
	return
end function strcmp

LOGICAL function isEq(x1,x2,n) !comparison between two scalars
	IMPLICIT NONE
	DOUBLE PRECISION x1,x2
	INTEGER n !number of decimal digits up to which check the equality
	
	if (abs(x1-x2).lt.0.5D0*10.**(-n)) then
		isEq=.true.
	else
		isEq=.false.
	end if
	
	return
end function isEq

function isEq_v(x1,x2,n) !comparison between a vector and a scalar
	IMPLICIT NONE
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: x1
	DOUBLE PRECISION, INTENT(in) :: x2
	INTEGER, INTENT(in) :: n	
	INTEGER i
	LOGICAL, DIMENSION(size(x1)) :: isEq_v
	
	do i=1,size(x1)
		if (abs(x1(i)-x2).lt.0.5D0*10.**(-n)) then
			isEq_v(i)=.true.
		else
			isEq_v(i)=.false.
		end if
	end do
	
	return
end function isEq_v

function isEq_vv(x1,x2,n) !comparison between two vectors
	IMPLICIT NONE
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: x1,x2
	INTEGER, INTENT(in) :: n	
	INTEGER i
	LOGICAL, DIMENSION(size(x1)) :: isEq_vv
	
	do i=1,size(x1)
		if (abs(x1(i)-x2(i)).lt.0.5D0*10.**(-n)) then
			isEq_vv(i)=.true.
		else
			isEq_vv(i)=.false.
		end if
	end do
	
	return
end function isEq_vv

function strfind(str,ch)
	IMPLICIT NONE
	CHARACTER(LEN=1000) :: str
	CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: ch
	INTEGER, DIMENSION(:), ALLOCATABLE :: strfind
	
	INTEGER, DIMENSION(len_trim(str)) :: ix
	INTEGER is,is1,j
	
	interface
		subroutine allocVecI(vTmp,n,v)
			IMPLICIT NONE
			INTEGER, DIMENSION(:), INTENT(in) :: vTmp
			INTEGER, DIMENSION(:), ALLOCATABLE :: v
			INTEGER n
		end subroutine allocVecI
	end interface
	
	j=0
	is=scan(trim(str),ch)
	if (is.ne.0) then !the result of scan is zero if no 'ch' is found in str
		j=j+1
		ix(j)=is
	end if
	is1=1 !initialzed so to enter the while loop
	do while ((is1.ne.0).and.(is.lt.len_trim(str)))
		is1=scan(str(is+1:len_trim(str)),ch) !trim()
		if (is1.ne.0) then
			is=is+is1
			j=j+1
			ix(j)=is
		end if
	end do
	
	call allocVecI(ix,j,strfind)
	
end function strfind

function polyval(v,x)
	!v contains coefficients of a polynomial in descendent degree order
	!polyval evaluates the polynomial at the value of x
	IMPLICIT NONE
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: v
	DOUBLE PRECISION x,polyval
	
	INTEGER n,i
	DOUBLE PRECISION add
	
	n=size(v)
	polyval=0.D0
	do i=1,n
		add=v(i)*x**(n-i)
		polyval=polyval+add
	end do
	
	return
end function polyval

subroutine removeColumn(matrix,titolo,field,newMatrix)
	!Check if there are columns reporting field in each of their rows. In this case,
	! cut these columns e produce the newMatrix without these columns 
	!Adapt the column titles as well 

	IMPLICIT NONE
	! Arguments declarations
	DOUBLE PRECISION field
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: matrix
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: newMatrix
	CHARACTER(LEN=1000) :: titolo
	! Variable declarations
	INTEGER ncol,jc,jr,lt,i,kr,j,bn
	INTEGER, DIMENSION(size(matrix)) :: ictmp,irtmp
	INTEGER, DIMENSION(:), ALLOCATABLE :: ic,ir,ndxC,XndxR1
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: ndxR,ndxRR
	CHARACTER(LEN=:), ALLOCATABLE :: ch
	CHARACTER(LEN=1000) :: titAdd,newTitolo
	LOGICAL isEq
	
	interface
		subroutine allocVecI(vTmp,n,v)
			IMPLICIT NONE
			INTEGER, DIMENSION(:), INTENT(in) :: vTmp
			INTEGER, DIMENSION(:), ALLOCATABLE :: v
			INTEGER n
		end subroutine allocVecI
		subroutine indexI(arr,indx)
			IMPLICIT NONE
			INTEGER, DIMENSION(:), INTENT(in) :: arr
			INTEGER indx(size(arr))
 		end subroutine indexI
		
		function strfind(str,ch)
			IMPLICIT NONE
			CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: ch
			CHARACTER(LEN=1000) :: str
			INTEGER, DIMENSION(:), ALLOCATABLE :: strfind
		end function strfind
	end interface

	ncol=size(matrix)
	jc=0
	jr=0
	do i=1,ncol 
		if (.not.isEq(matrix(i),field,2)) then
			jc=jc+1
			ictmp(jc)=i
		else
			jr=jr+1
			irtmp(jr)=i !index of removed columns
		end if 
	end do
	if (jc.gt.0) then
		call allocVecI(ictmp,jc,ic)
		allocate(newMatrix(jc))
		newMatrix=matrix(ic)
	end if
	if (jr.gt.0) then
		call allocVecI(irtmp,jr,ir)
	end if
	
	lt=len_trim(titolo)
	
	!If every matrix element==field, newMatrix remains deallocated
	if (jr==0) then !no column rejection
		newTitolo=titolo
	else
		ch=' '
		ndxC=strfind(titolo,ch)
		kr=1
		allocate(ndxR(size(ir),2))
		do j=1,size(ir) !Remove a field of the type ___' ' except for the last that is removed without blank ___ 
			if (ir(j)==1) then
				ndxR(kr,:)=(/ 1,ndxC(ir(j)) /) !ir(j)=1
				kr=kr+1
			else if (ir(j)==size(matrix)) then !,2
				ndxR(kr,:)=(/ ndxC(ir(j)-1)+1,lt /) !if columns are n, blanks are n-1 => ir(j)-1=size(ndxC)
				kr=kr+1
			else
				ndxR(kr,:)=(/ ndxC(ir(j)-1)+1,ndxC(ir(j)) /) 
				kr=kr+1
			end if 
		end do
		allocate(XndxR1(size(ndxR,1)))
		call indexI(ndxR(:,1),XndxR1)
		allocate(ndxRR(size(ndxR,1),size(ndxR,2)))
		ndxRR=ndxR(XndxR1,:)
		if (size(ndxRR,1)==1) then 
			if (ndxRR(1,1)==1) then 
				newTitolo=titolo(ndxRR(1,2)+1:lt)
			else if (ndxRR(1,2)==lt) then
				newTitolo=titolo(1:ndxRR(1,1)-1)
			else
				newTitolo=titolo(1:ndxRR(1,1)-1)
				titAdd=titolo(ndxRR(1,2)+1:lt)
				newTitolo=trim(newTitolo)//' '//titAdd !' ' necessary because external blanks in strings are trimmed
			end if 
		else
			bn=0 !check if newTitolo already defined
			if (ndxRR(1,1).ne.1) then 
				newTitolo=titolo(1:ndxRR(1,1)-1)
				bn=1
			end if 
			j=1
			do while (j<=size(ndxRR,1)-1) 
				if (ndxRR(j,2)+1.ne.ndxRR(j+1,1)) then
						titAdd=titolo(ndxRR(j,2)+1:ndxRR(j+1,1)-1)
						if (bn.eq.1) then
							newTitolo=trim(newTitolo)//' '//titAdd
						else
							newTitolo=trim(titAdd)
							bn=1
						end if
				end if 
				j=j+1
			end do 
			if (ndxRR(size(ndxRR,1),2).ne.lt) then
				titAdd=titolo(ndxRR(size(ndxRR,1),2)+1:lt)
				if (bn.eq.1) then
					newTitolo=trim(newTitolo)//' '//titAdd
				else
					newTitolo=trim(titAdd)
				end if
			end if 
		end if
	end if
	
	titolo=newTitolo

	return
end subroutine removeColumn

subroutine avoidLastFlag(Isoc)
	USE IsoPD, only: clab,vlab,tlab	
	
	IMPLICIT NONE
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: Isoc
	
	DOUBLE PRECISION, DIMENSION(size(Isoc,1)) :: lab
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Isoc1
	LOGICAL, DIMENSION(size(Isoc,1)) :: b
	INTEGER, DIMENSION(:), ALLOCATABLE :: ndx
	
	interface
		subroutine find(b,ix)
			IMPLICIT NONE
			LOGICAL, DIMENSION(:), INTENT(in) :: b
			INTEGER, DIMENSION(:), ALLOCATABLE :: ix
		end subroutine find
	end interface

	lab=Isoc(:,clab)
	b=(abs(lab-vlab).gt.tlab)
	call find(b,ndx)
	allocate(Isoc1(size(ndx),size(Isoc,2)))
	Isoc1=Isoc(ndx,:)
	deallocate(Isoc)
	allocate(Isoc(size(Isoc1,1),size(Isoc1,2)))
	Isoc=Isoc1

end subroutine avoidLastFlag

subroutine wMean(x1,I_x1,x2,I_x2,x,I_x)
	!Compute the weighted mean between x1 and x2 considering
	!I_x1^-2 and I_x2^-2 as weights
	IMPLICIT NONE
	DOUBLE PRECISION x1,I_x1,x2,I_x2,x,I_x
	DOUBLE PRECISION w1,w2
	
	w1=I_x1**(-2)
	w2=I_x2**(-2)
	x=(x1*w1+x2*w2)/(w1+w2)
	I_x=1./sqrt(w1+w2)
end subroutine wMean

subroutine selectMfromTracks(MTrAv,M,I_M,Mvec)
	IMPLICIT NONE
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: MTrAv
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Mvec
	DOUBLE PRECISION M,I_M
	
	LOGICAL, DIMENSION(size(MTrAv)) :: Mb
	INTEGER, DIMENSION(:), ALLOCATABLE :: ix
		
	interface
		subroutine find(b,ix)
			IMPLICIT NONE
			LOGICAL, DIMENSION(:), INTENT(in) :: b
			INTEGER, DIMENSION(:), ALLOCATABLE :: ix
		end subroutine find
	end interface
	
	Mb=(MTrAv.ge.M-I_M .and. MTrAv.le.M+I_M)
	
	call find(Mb,ix)
	if (allocated(ix)) then
		allocate(Mvec(size(ix)))
		Mvec=MTrAv(ix)
	else !Mvec remains deallocated
		print*,'from selectMfromTracks: no mass available'
	end if
end subroutine selectMfromTracks

subroutine consistentM(TeB,cumInt,Mvec,M,agreeM)
	IMPLICIT NONE
	LOGICAL, DIMENSION(:), INTENT(in) :: TeB
	LOGICAL agreeM
	INTEGER, DIMENSION(:), INTENT(in) :: cumInt
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: Mvec
	DOUBLE PRECISION M
	
	LOGICAL, DIMENSION(size(cumInt)) :: cumb
	INTEGER, DIMENSION(:), ALLOCATABLE :: TeB1,ix,xM
	INTEGER jj
	DOUBLE PRECISION Mm
	interface
		subroutine find(b,ix)
			IMPLICIT NONE
			LOGICAL, DIMENSION(:), INTENT(in) :: b
			INTEGER, DIMENSION(:), ALLOCATABLE :: ix
		end subroutine find
	end interface
		
	call find(TeB,TeB1)
	if (allocated(TeB1)) then
		allocate(xM(size(TeB1)))
		do jj=1,size(TeB1)
			cumb=(cumInt.ge.TeB1(jj))
			call find(cumb,ix)
			xM(jj)=ix(1)
			deallocate(ix)
		end do
		Mm=sum(Mvec(xM))/size(Mvec(xM))

		if (abs((Mm-M)/M).lt.0.1) then
			agreeM=.true.
		else
			M=Mm
			agreeM=.false.
		end if
	else
		agreeM=.false.
	end if
end subroutine consistentM

subroutine find(b,ix)
	!Given an array of logicals, find all the indices where .true. is placed
	!and store all of them in the vector ix. If b contains only .false.
	!then ix remains deallocated
	IMPLICIT NONE
	LOGICAL, DIMENSION(:), INTENT(in) :: b
	INTEGER, DIMENSION(:), ALLOCATABLE :: ix
	
	LOGICAL, DIMENSION(size(b)) :: mask
	INTEGER sB,sI,i
		
	sB=size(b)
	mask=b.eqv..true.
	sI=count(mask)
	if (sI.gt.0) then !
		allocate(ix(sI))
		ix=pack((/(i,i=1,sB)/),mask)
	end if
	
end subroutine find

subroutine computeLRgfromMrho(M,I_M,rho,I_rho,Teff,I_Teff,L,logL,I_L,I_logL,R,I_R,g,logg)
	USE SunPD, only: TeffSun,rhoSun,loggSun
	
	IMPLICIT NONE
	DOUBLE PRECISION M,I_M,rho,I_rho,Teff,I_Teff
	DOUBLE PRECISION L,logL,I_L,I_logL,R,I_R,g,logg
	
	DOUBLE PRECISION, PARAMETER :: e=2.71828182845905 !Nepero number
	
	g=(M*(rho/rhoSun)**2)**(1./3.)*10**loggSun !mantain original input uncertainty
	logg=log10(g)
	R=(M/(rho/rhoSun))**(1./3.)
	I_R=1./3.*(I_M/M+I_rho/rho)*R
	L=R**2*(Teff/TeffSun)**4
	logL=log10(L)
	I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
	I_logL=I_L/L*log10(e)

end subroutine computeLRgfromMrho

subroutine computeLRrhofromMg(M,I_M,g,I_g,Teff,I_Teff,L,logL,I_L,I_logL,R,I_R,rho,logrho)
	USE SunPD, only: TeffSun,rhoSun,loggSun
	
	IMPLICIT NONE
	DOUBLE PRECISION M,I_M,g,I_g,Teff,I_Teff
	DOUBLE PRECISION L,logL,I_L,I_logL,R,I_R,rho,logrho
	
	DOUBLE PRECISION, PARAMETER :: e=2.71828182845905 !Nepero number
	
	rho=sqrt((g/10**loggSun)**3/M)*rhoSun !g/cm3. Mantain original input uncertainty
	logrho=log10(rho)
	R=sqrt(M/(g/10**loggSun))
	I_R=0.5*(I_M/M+I_g/g)*R
	L=R**2*(Teff/TeffSun)**4
	logL=log10(L)
	I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
	I_logL=I_L/L*log10(e)

end subroutine computeLRrhofromMg

subroutine computeLMg(R,I_R,Teff,I_Teff,rho,I_rho,L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg)
	USE SunPD, only: TeffSun,rhoSun,loggSun
	
	IMPLICIT NONE
	DOUBLE PRECISION R,I_R,Teff,I_Teff,rho,I_rho
	DOUBLE PRECISION L,logL,I_L,I_logL,M,I_M,g,I_g,logg,I_logg
	
	DOUBLE PRECISION, PARAMETER :: e=2.71828182845905 !Nepero number
	
	L=R**2*(Teff/TeffSun)**4
	logL=log10(L)
	I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
	I_logL=I_L/L*log10(e)
	M=rho/rhoSun*R**3 !Mo
	I_M=M*(I_rho/rho+3.*I_R/R)
	g=M/R**2*10.**loggSun !cgs
	I_g=g*(I_M/M+2.*I_R/R)
	logg=log10(g)
	I_logg=log10(e)*I_g/g
	
end subroutine computeLMg

subroutine computeMg(R,I_R,Teff,I_Teff,rho,I_rho,M,I_M,g,I_g,logg,I_logg)
	USE SunPD, only: TeffSun,rhoSun,loggSun
	
	IMPLICIT NONE
	DOUBLE PRECISION R,I_R,Teff,I_Teff,rho,I_rho
	DOUBLE PRECISION M,I_M,g,I_g,logg,I_logg
	
	DOUBLE PRECISION, PARAMETER :: e=2.71828182845905 !Nepero number
	
	M=rho/rhoSun*R**3 !Mo
	I_M=M*(I_rho/rho+3.*I_R/R)
	g=M/R**2*10.**loggSun !cgs
	I_g=g*(I_M/M+2.*I_R/R)
	logg=log10(g)
	I_logg=log10(e)*I_g/g

end subroutine computeMg

subroutine computeLMrho(R,I_R,Teff,I_Teff,g,I_g,L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho)
	USE SunPD, only: TeffSun,rhoSun,loggSun
	
	IMPLICIT NONE
	DOUBLE PRECISION R,I_R,Teff,I_Teff,g,I_g
	DOUBLE PRECISION L,logL,I_L,I_logL,M,I_M,rho,I_rho,logrho,I_logrho
	
	DOUBLE PRECISION, PARAMETER :: e=2.71828182845905 !Nepero number

	L=R**2*(Teff/TeffSun)**4
	logL=log10(L)
	I_L=L*(2.*I_R/R+4.*I_Teff/Teff)
	I_logL=I_L/L*log10(e)
	M=g/10.**loggSun*R**2
	I_M=M*(I_g/g+2.*I_R/R)
	rho=M/R**3*rhoSun
	I_rho=rho*(I_M/M+3.*I_R/R)
	logrho=log10(rho)
	I_logrho=log10(e)*I_rho/rho
	
end subroutine computeLMrho

subroutine computeMrho(R,I_R,Teff,I_Teff,g,I_g,M,I_M,rho,I_rho,logrho,I_logrho)
	USE SunPD, only: TeffSun,rhoSun,loggSun
	
	IMPLICIT NONE
	DOUBLE PRECISION R,I_R,Teff,I_Teff,g,I_g
	DOUBLE PRECISION M,I_M,rho,I_rho,logrho,I_logrho
	
	DOUBLE PRECISION, PARAMETER :: e=2.71828182845905 !Nepero number

	M=g/10.**loggSun*R**2
	I_M=M*(I_g/g+2.*I_R/R)
	rho=M/R**3*rhoSun
	I_rho=rho*(I_M/M+3.*I_R/R)
	logrho=log10(rho)
	I_logrho=log10(e)*I_rho/rho
	
end subroutine computeMrho

subroutine computeVrif(TabRif,cRif,M,vRif)
	IMPLICIT NONE
	DOUBLE PRECISION M,vRif
	INTEGER cRif
	
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Mvel,vvel
	DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: TabRif
	INTEGER kMRif
		
	allocate(Mvel(size(TabRif,1)));allocate(vvel(size(TabRif,1)))
	Mvel=TabRif(:,1)
	vvel=TabRif(:,cRif)
	kMRif=minloc(abs(M-Mvel),1) 
	vRif=vvel(kMRif)

end subroutine computeVrif

subroutine computeVevo(t_Tracks,y_Tracks,logTeff_Tracks,t_i,varTrlim,vEvo) !,back
	IMPLICIT NONE
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: t_Tracks,y_Tracks,logTeff_Tracks,t_i
	DOUBLE PRECISION, DIMENSION(size(t_i)) :: vEvo
	DOUBLE PRECISION varTrlim
		
	DOUBLE PRECISION vY,vT,dminTr,dminTr2 !!
	DOUBLE PRECISION, DIMENSION(size(t_Tracks)) :: d_Tracks
	INTEGER i,kminTr,stepTr
	LOGICAL back
	
	interface
		subroutine stepEnoughVar(tab,ix,varLim,colD,stepMax,backDir)
			IMPLICIT NONE
			DOUBLE PRECISION varLim
			DOUBLE PRECISION, DIMENSION(:,:) :: tab
			INTEGER ix,colD,stepMax
			LOGICAL backDir
		end subroutine stepEnoughVar
	end interface
	
	do i=1,size(t_i)
		d_Tracks=abs(t_i(i)-t_Tracks)
		dminTr=minval(d_Tracks)
		kminTr=minloc(d_Tracks,1)
		dminTr2=d_Tracks(kminTr) 
		call stepEnoughVar(reshape((/t_Tracks,y_Tracks,logTeff_Tracks/),(/size(t_Tracks,1),3/)), &
			& kminTr,varTrlim,1,stepTr,back) 
		if (.not.back) then !back is false -> go forward
			vY=(y_Tracks(kminTr+stepTr)-y_Tracks(kminTr))/(t_Tracks(kminTr+stepTr)-t_Tracks(kminTr))
			vT=(logTeff_Tracks(kminTr+stepTr)-logTeff_Tracks(kminTr))/(t_Tracks(kminTr+stepTr)-t_Tracks(kminTr))
		else
			vY=(y_Tracks(kminTr)-y_Tracks(kminTr-stepTr))/(t_Tracks(kminTr)-t_Tracks(kminTr-stepTr))
			vT=(logTeff_Tracks(kminTr)-logTeff_Tracks(kminTr-stepTr))/(t_Tracks(kminTr)-t_Tracks(kminTr-stepTr))
		end if 
		vEvo(i)=sqrt(vT**2+vY**2)
	end do

end subroutine computeVevo

subroutine stepEnoughVar(tab,ix_,varLim,colD,stepMax,backDir)
	IMPLICIT NONE
	! Arguments declarations
	DOUBLE PRECISION varLim
	DOUBLE PRECISION, DIMENSION(:,:) :: tab
	INTEGER ix_,colD,stepMax
	LOGICAL backDir
	! Variable declarations
	DOUBLE PRECISION x1,x2,vRel
	DOUBLE PRECISION, DIMENSION(size(tab,1)) :: x
	INTEGER ix,ix0,i,step,ndxStep,bt,j
	INTEGER, DIMENSION(size(tab,2)) :: steps
	INTEGER, DIMENSION(:), ALLOCATABLE :: btx,stepsBck
	LOGICAL back
	LOGICAL, DIMENSION(size(tab,2)) :: backs
	
	if (colD.ne.0) then !Put colD column as 1st column (if it isn't yet)
		if (colD.ne.1) then 
			if (colD==size(tab,2)) then 
				tab=reshape((/tab(:,colD),tab(:,1:colD-1)/),(/size(tab,1),size(tab,2)/))
			else
				tab=reshape((/tab(:,colD),tab(:,1:colD-1),tab(:,colD+1:size(tab,2))/),(/size(tab,1),size(tab,2)/))
			end if 
		end if 
	end if 
	
	ix=ix_ !to avoid that the input ix_ will be changed in this routine and used in the main program as changed value
	ix0=ix
	back=.true.
	do i=1,size(tab,2)
		x=tab(:,i)
		x1=x(ix0)
		step=1
		if ((colD.ne.0).and.(i>1).and.(back)) then !variable at denominator requires
		! to go backward so that differences between its values are non zeros
			ix=ix0-step
			x2=x(ix)
			vRel=abs((x2-x1)/x1)
			do while (vRel<varLim .and. ix>1) !not possible that ix becomes 1 
				x1=x2
				step=step+1
				ix=ix0-step
				x2=x(ix)
				vRel=abs((x2-x1)/x1)
			end do 
		else
			if (ix0==size(x)-1) then 
				ix=ix0+step
				x2=x(ix)
				vRel=abs((x2-x1)/x1)
				if (vRel<varLim) then 
				!reset to initial values and consider variations going backwards in the following
					back=.true.
					step=1
					ix=ix0
					x1=x(ix)
				else
					back=.false.
				end if 
			else if (ix0<size(x)-1) then !at least two steps forward are possible
				back=.false.
				ix=ix0+step
				x2=x(ix)
				vRel=abs((x2-x1)/x1)
				do while (vRel<varLim .and. ix<size(x)) 
					x1=x2
					step=step+1
					ix=ix0+step
					x2=x(ix)
					vRel=abs((x2-x1)/x1)
					if (ix==size(x)) then 
						if (vRel<varLim) then 
							back=.true. !reset to initial values
							step=1
							ix=ix0
							x1=x(ix)
						else
							back=.false.
						end if 
						exit
					end if 
				end do 
			end if 
			if (ix0==size(x) .or. back) then ! (only going backwards is possible)
				ix=ix0-step
				x2=x(ix)
				vRel=abs((x2-x1)/x1)
				do while (vRel<varLim .and. ix>1) !not possible that ix becomes 1 
					x1=x2
					step=step+1
					ix=ix0-step
					x2=x(ix)
					vRel=abs((x2-x1)/x1)
				end do 
			end if 
		end if 
		steps(i)=step
		backs(i)=back
	end do
	if (colD.ne.0 .and. backs(1).eqv..true.) then  
	!there's denominator and you have to go backwards so that Denominator!=0
		stepMax=maxval(steps)
		ndxStep=maxloc(steps,1)
		backDir=backs(ndxStep) !should be 1
	else
		bt=count(backs)
		if (bt==0) then !all steps are forward
			stepMax=maxval(steps)
			ndxStep=maxloc(steps,1)
			backDir=backs(ndxStep) !should be 0
		else
			allocate(btx(bt))
			btx=pack((/(j,j=1,size(backs))/),backs) !indices of backs for which backs elements are true
			allocate(stepsBck(bt))
			stepsBck=steps(btx)
			stepMax=maxval(stepsBck)
			ndxStep=maxloc(stepsBck,1)
			backDir=backs(ndxStep) !should be 1
		end if 
	end if 

	return
end subroutine stepEnoughVar

subroutine ifComputeIsocAge(ZAMStab,cZAMS,logY0,logt_iso,logTeff_iso,yT_iso,I_logTeff,ageIsoc)
	!Decide whether to compute isochronal age, according to the quickness of the star evolution
	!The speed is evaluated considering the difference in stellar Teff at the ZAMS phase
	!and at the age of the universe
	IMPLICIT NONE
	DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: ZAMStab
	DOUBLE PRECISION logY0,I_logTeff
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: logt_iso,logTeff_iso,yT_iso
	INTEGER cZAMS,ageIsoc
	
	DOUBLE PRECISION logt_ZAMS_
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: logt_ZAMSv,logTeMS1,logTeMS2,logTeMS !logt_ZAMS,logY_ZAMS
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: yTTeff1,yTTeff2 !ZAMStab,
	INTEGER jiMS,niMS1,nfMS1,niMS2,nfMS2
	INTEGER, DIMENSION(:), ALLOCATABLE :: niMSv,nfMSv
		
	interface
		subroutine findYgivenX_i(M,x,xx,ny,y)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: M
			DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: y
			DOUBLE PRECISION x
			INTEGER xx,ny 
		end subroutine findYgivenX_i
		subroutine trovaIndici(v,ndxi,ndxf)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: v
			INTEGER, DIMENSION(:), ALLOCATABLE :: ndxi,ndxf
		end subroutine trovaIndici
	end interface

	call findYgivenX_i(ZAMStab,logY0,cZAMS,2,logt_ZAMSv) !,1,2)
	if (.not.allocated(logt_ZAMSv)) then !logY0 out of available range. Since the
		!if statement logY0<logy0_soglia where this subroutine works implies low MS, 
		!the star is even dimmer than the values reported by the tracks 
		logt_ZAMS_=9 !very low MS star => very slow evolution. 0.1Mo-stars have logtZAMS.not.(9)
	else
		logt_ZAMS_=minval(logt_ZAMSv)
	end if
			
	call trovaIndici(logt_iso,niMSv,nfMSv) 
	jiMS=minloc(abs(logt_ZAMS_-logt_iso(niMSv)),1) 
	niMS1=niMSv(jiMS)
	nfMS1=nfMSv(jiMS)
	jiMS=size(logt_iso(niMSv)) !easier than the two previous instructions
	niMS2=niMSv(jiMS)
	nfMS2=nfMSv(jiMS)
	allocate(yTTeff1(nfMS1-niMS1+1,2)); allocate(yTTeff2(nfMS2-niMS2+1,2))
	yTTeff1=reshape((/ yT_iso(niMS1:nfMS1),logTeff_iso(niMS1:nfMS1) /),(/nfMS1-niMS1+1,2 /))
	yTTeff2=reshape((/ yT_iso(niMS2:nfMS2),logTeff_iso(niMS2:nfMS2) /),(/nfMS2-niMS2+1,2 /))
		
	call findYgivenX_i(yTTeff1,logY0,1,2,logTeMS1)
	call findYgivenX_i(yTTeff2,logY0,1,2,logTeMS2)
	deallocate(yTTeff1); deallocate(yTTeff2)
	
	allocate(logTeMS(size(logTeMS1)+size(logTeMS2)))
					
	logTeMS=(/ logTeMS1,logTeMS2 /) !reshape( ,(/ /))
	if (maxval(logTeMS)-minval(logTeMS)<I_logTeff) then 
	!!!(more restrictive--> .or. max(logLSub_T)-min(logLSub_T)<I_logL)
		ageIsoc=0 !do not compute isochronal age in the case of slowly-evolving stars
	else
		ageIsoc=1
	end if

end subroutine ifComputeIsocAge

subroutine setThreshold(caliblogL,loggAvail,logY2,logY2_soglia,logY0,logY0_soglia,cZAMS,cyT)
!Re-change sign to logg or rho in case they are used as threshold value
!Set the column index for ZAMS tab
!Set the column index for _iso
	USE IsoPD
	
	IMPLICIT NONE
	LOGICAL caliblogL,loggAvail
	DOUBLE PRECISION logY2,logY2_soglia,logY0,logY0_soglia
	INTEGER cZAMS,cyT
	
	if (.not.caliblogL) then
		logY0=-logY2
		logY0_soglia=-logY2_soglia
		if (loggAvail) then 
			cZAMS=cgZAMS !column index of logg
			cyT=clogg
		else
			cZAMS=crhoZAMS !column index of logrho
			cyT=0 !there isn't a direct column index for rho in the isochrone grid
		end if 
	else !logL does exist!
		logY0=logY2
		logY0_soglia=logY2_soglia
		cZAMS=cLZAMS !column index of logL
		cyT=clogL
	end if

end subroutine setThreshold

subroutine multipleZisocPhSCP(FeH,I_FeH,Z_iso,model,x,y,logt8,percorso,k,symmZ, &
		& calibHRD,calibNoD,calibSPEC,loggAvail,caliblogL,hstar,hstarlim,Isoc_i, &
		& Z_low,Z_up,useColor,idCol)
	!symmZ=0. Load isochrones considering symmetric error bars in FeH (and not in Z) 
	!symmZ=1. Load isochrones considering "symmetric" error bars in Z (and not in FeH)
	!actually the most metallic is missing to avoid biases 
	!x abscissa variable according to the calibration (e.g. BmV)
	!y ordinata variable according to the calibration
	!x_iso abscissa isochrone variable which is the same of x (e.g. BmV)
	!xx_iso abscissa isochrone variable which is the other one (e.g. logTeff)
	!So x <--> x_iso; instead xx_iso is the other variable
	USE SunPD, ONLY: TeffSun,rhoSun,costZ,YHe=>He0
	USE IsoPD, clogt=>ct
	
	IMPLICIT NONE
	! Arguments declarations
	CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: model,percorso !
		
	INTEGER symmZ,idCol
	
	DOUBLE PRECISION FeH,I_FeH,x,y,k,hstar,hstarlim,Z_low,Z_up,logt8
	DOUBLE PRECISION useColor
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: Z_iso
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Isoc_i
	
	LOGICAL calibHRD,calibNoD,calibSPEC,loggAvail,caliblogL
	
	! Variable declarations
	CHARACTER(LEN=:), ALLOCATABLE :: Ystr,Z_lstr,Z_ustr,Z_inputS
	CHARACTER(LEN=:), ALLOCATABLE :: strMod1,strMod2,strMod3,strMod0
	CHARACTER(LEN=:), ALLOCATABLE :: fn1,fn2,fn3
	CHARACTER(LEN=5) Z_inputSLG
	CHARACTER(LEN=1000) head
	
	INTEGER indxZlow,indxZup,indxZinput
	INTEGER n_is,n,ilogt8,rowAge,Nage,ndxDU,j,cmag0
	INTEGER, DIMENSION(:), ALLOCATABLE :: dist_ndx,dist_ndxU,ix !ndxti,ndxtf,
	
	DOUBLE PRECISION Z_isoLOW,Z_isoUP
	DOUBLE PRECISION qgp,logt8c,Inct !,logtstep
	DOUBLE PRECISION Z_low_,Z_up_,Z_input,FeH_step,dist0
	DOUBLE PRECISION x_i,x_i1,x_i2,BmV_i1,BmV_i2,logg_i1,logg_i2
	DOUBLE PRECISION logL_i1,logL_i2,logrho_i1,logrho_i2,logTeff_i1,logTeff_i2
	DOUBLE PRECISION M_i1,M_i2,V_i1,V_i2,mBV,mg,mL,mM,mrho,mT,qBV,qg,qL,qM,qrho,qT
	DOUBLE PRECISION Z_input2n,Z_input2
	
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: logtiso,FeHcol,ageMinDist !,Vsteps
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: t_iso,logt_iso,M_iso,L_iso,logL_iso
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Teff_iso,logTeff_iso,g_iso,logg_iso
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: R_iso,rho_iso,logrho_iso,V_iso
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: BmV_iso,BC_iso,x_iso,xx_iso,y_iso
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: t_iTmp,BC_iTmp,logTeff_iTmp,BmV_iTmp
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: logg_iTmp,logL_iTmp,M_iTmp,logrho_iTmp
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: t_i,BC_i,logTeff_i,BmV_i,logg_i,logL_i
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: M_i,logrho_i,Teff_i,L_i,g_i,rho_i
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Isoc,IsocTmp,dist,distTmp
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Isoc1,Isoc_i1
	
	LOGICAL isEq

	interface
		subroutine num2str(x,d,str)
			IMPLICIT NONE
			DOUBLE PRECISION x
			CHARACTER(LEN=:), ALLOCATABLE :: str !was 40
			INTEGER d
		end subroutine num2str
		subroutine loadMatrix(fileName,matrix,head)
			IMPLICIT NONE
			CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: fileName
			DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: matrix
			CHARACTER(LEN=1000) :: head
		end subroutine loadMatrix
		subroutine trovaIndici(v,ndxi,ndxf)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: v
			INTEGER, DIMENSION(:), ALLOCATABLE :: ndxi,ndxf
		end subroutine trovaIndici
		subroutine indexx(arr,indx)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: arr !(n)
			INTEGER indx(size(arr))
		end subroutine indexx
		subroutine uniqueFast(list,d,indices,first)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: list
			INTEGER d
			INTEGER, DIMENSION(:), ALLOCATABLE :: indices
			LOGICAL first
		end subroutine uniqueFast
		subroutine choose_i2Col(ndxDU,x,x_iso,xx_iso,calibSPEC,hstar,hstarlim,V_iso,logL_iso, &
			& M_iso,logg_iso,logrho_iso,logt_iso,BmV_i2,V_i2,logTeff_i2,logL_i2,M_i2,logg_i2,logrho_i2)
			IMPLICIT NONE
			INTEGER ndxDU
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: x_iso
			DOUBLE PRECISION, DIMENSION(size(x_iso)) :: xx_iso
			DOUBLE PRECISION, DIMENSION(size(x_iso)) :: V_iso,logL_iso
			DOUBLE PRECISION, DIMENSION(size(x_iso)) :: M_iso,logg_iso,logrho_iso,logt_iso
			DOUBLE PRECISION x,BmV_i2,V_i2,logTeff_i2,logL_i2,M_i2,logg_i2,logrho_i2
			DOUBLE PRECISION hstar,hstarlim
			LOGICAL calibSPEC
		end subroutine choose_i2Col
		subroutine append(A,B)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: A
			DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: B
		end subroutine append
	end interface
	interface
		function strncmp(str1,str2,n)
			IMPLICIT NONE
			CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: str1,str2
			INTEGER n
			LOGICAL strncmp
		end function strncmp
		function strcmp(str1,str2)
			IMPLICIT NONE
			CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: str1,str2
			LOGICAL strcmp
		end function strcmp
		
		function isEq_v(x1,x2,n)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: x1
			DOUBLE PRECISION, INTENT(in) :: x2
			INTEGER, INTENT(in) :: n
			LOGICAL, DIMENSION(size(x1)) :: isEq_v
		end function isEq_v
		
		function ones(n)
			IMPLICIT NONE
			INTEGER n
			DOUBLE PRECISION, DIMENSION(n) :: ones
		end function
	end interface
	
	strMod1="PD"
	strMod2="LG"
	strMod3="PD1.2S"
	strMod0="PD"
	
	select case (idCol)
	case (1)
		cmag0=cmag1
	case (2)
		cmag0=cmagx
	case default
		cmag0=cmag1
	end select

	Z_isoLOW=minval(Z_iso)
	Z_isoUP=maxval(Z_iso)
	logt8c=logt8 !if the activity check has been done, logt8c=0. 
		!!logt8 will change in the following as 10**logt8 in the case of PD1.2S model
	if (strcmp(model,strMod3)) then 
		logt8=10.**logt8 !PD1.2S reports t and not logt
	end if
	
	if (strcmp(model,strMod2)) then 
		call num2str(YHe,3,Ystr) !3 is max number of decimals. If less, right zeros are trimmed
	end if 

	if (FeH-k*I_FeH<log10(0.010)+costZ .and. strncmp(model,strMod0,2)) then 
		Z_low_=10**(FeH-k*I_FeH-costZ)
		indxZlow=minloc(abs(Z_low_-Z_iso),1) 
		Z_low=Z_iso(indxZlow)
		Z_lstr=repeat(achar(48),7)
		write(Z_lstr,'(F7.5)') Z_low
		if (symmZ==1) then 
			Z_up_=10**(FeH-costZ)+(10**(FeH-costZ)-Z_low_)
			indxZup=minloc(abs(Z_up_-Z_iso),1) 
			Z_up=Z_iso(indxZup)
			if (Z_up<0.010) then
				Z_ustr=repeat(achar(48),7)
				write(Z_ustr,'(F7.5)') Z_up
			else
				Z_ustr=repeat(achar(48),5)
				write(Z_ustr,'(F5.3)') Z_up
			end if 
		end if
	else
		Z_low_=10**(FeH-k*I_FeH-costZ)
		indxZlow=minloc(abs(Z_low_-Z_iso),1)
		Z_low=Z_iso(indxZlow)
		Z_lstr=repeat(achar(48),5)
		write(Z_lstr,'(F5.3)') Z_low
		if (symmZ==1) then 
			Z_up_=10**(FeH-costZ)+(10**(FeH-costZ)-Z_low_)
			indxZup=minloc(abs(Z_up_-Z_iso),1) 
			Z_up=Z_iso(indxZup)
			Z_ustr=repeat(achar(48),5)
			write(Z_ustr,'(F5.3)') Z_up 
		end if
	end if 

	if (symmZ==0) then
		if (FeH+k*I_FeH<log10(0.010)+costZ .and. strncmp(model,strMod0,2)) then 
			Z_up_=10**(FeH+k*I_FeH-costZ)
			indxZup=minloc(abs(Z_up_-Z_iso),1) 
			Z_up=Z_iso(indxZup)
			Z_ustr=repeat(achar(48),7)
			write(Z_ustr,'(F7.5)') Z_up
		else
			Z_up_=10**(FeH+k*I_FeH-costZ)
			indxZup=minloc(abs(Z_up_-Z_iso),1) 
			Z_up=Z_iso(indxZup)
			Z_ustr=repeat(achar(48),5)
			write(Z_ustr,'(F5.3)') Z_up
		end if 
	end if 

	! dZ=max(round((0.04*log(10.)*Z)*1000)/1000,0.001); 
	! n_is=floor((Z_up-Z_low)/dZ); !numbers of isochrones of different metallicity
	! to be considered besides the one with Z=Z_low 
	! Considering the steps in Z, isochrone separated by [Fe/H]_step=0.05 are considered
	! Jorgensen&Lindegren2005 considered a step of 0.04
	!1.115=1+FeH_step*ln(10) is the reason of geometric progression in Z (FeH_step=0.05)
	!!!!!!!!!!
	FeH_step=0.05
	!!!!!!!!!!

	qgp=1.+FeH_step*log(10.)
	n_is=floor(log(Z_up/Z_low)/log(qgp))
	
	if (n_is.ne.0) then 
		Z_input=Z_low
		n=0
		do while (Z_input<Z_isoLOW .and. n<n_is) 
			n=n+1
			Z_input=Z_low*qgp**n
		end do 
		indxZinput=minloc(abs(Z_input-Z_iso),1) 
		Z_input=Z_iso(indxZinput)
		if (n_is==n .or. Z_low>=Z_isoUP) then !I just deal with ONLY ONE grid of isochrones
			!=> _i have already been produced in the main
			!Isoc_i remains de-allocated
			return
		else
!			n_is=n_is+1 !in this case load the necessary grids of isochrones so that the
!			!metallicity interval is exactly symmetric in Z. The effect is to load an
!			!additional metallic grid of isochrone and this bias the age. It's better not
!			!to consider this option so that the most metallic grid to be load Zlim<ZupZS<SupFehS
!			!where ZupZS (ZupFehS) is the right limit that guarantee the symmetry in Z (in FeH)
!			!So the interval of metallic grids to be loaded is rather asymmetric giving
!			!preferennces to not very metallic grids around the stellar metalicity, but this choice
!			!is the best one in terms of rms and of avoiding biases (8/3/18)
			if (strcmp(model,strMod1)) then
				call num2str(Z_input,5,Z_inputS)
				fn1=trim(percorso)//trim(Z_inputS)//"M.txt"
				call loadMatrix(fn1,Isoc,head)
			else if (strcmp(model,strMod2)) then
				write(Z_inputSLG,'(F5.3)') Z_input 
				fn2=trim(percorso)//Ystr//"Z"//Z_inputSLG//".txt"
				call loadMatrix(fn2,Isoc,head)
			else if (strcmp(model,strMod3)) then
				call num2str(Z_input,5,Z_inputS)
				fn3=trim(percorso)//trim(Z_inputS)//".dat"
				call loadMatrix(fn3,Isoc,head)
			end if
			!!Check whether logt has equally spaced steps
			allocate(logtiso(size(Isoc,1)))
			if (strcmp(model,strMod3)) then 
				logtiso=log10(Isoc(:,clogt))
			else
				logtiso=Isoc(:,clogt)
			end if
			if (.not.isEq(logt8c,0.D0,2)) then !Activity check done
				ilogt8=1
				do while (abs(logt8-Isoc(ilogt8,clogt))>Inct/2. .and. ilogt8<size(Isoc,1)) 
					ilogt8=ilogt8+1
				end do 
				allocate(IsocTmp(size(Isoc,1)-ilogt8+1,size(Isoc,2)))
				IsocTmp=Isoc(ilogt8:size(Isoc,1),:)
				deallocate(Isoc)
				allocate(Isoc(size(IsocTmp,1),size(IsocTmp,2)))
				Isoc=IsocTmp
				deallocate(IsocTmp)
			end if
			allocate(FeHcol(size(Isoc,1)))
			FeHcol=(log10(Z_input)+costZ)*ones(size(Isoc,1))
			allocate(IsocTmp(size(Isoc,1),size(Isoc,2)+1))
			IsocTmp(:,1:size(IsocTmp,2)-1)=Isoc
			IsocTmp(:,size(IsocTmp,2))=FeHcol
			deallocate(Isoc); deallocate(FeHcol)
			allocate(Isoc(size(IsocTmp,1),size(IsocTmp,2)))
			Isoc=IsocTmp
			deallocate(IsocTmp)
			
			allocate(t_iso(size(Isoc,1))); allocate(logt_iso(size(Isoc,1)))
			if (strcmp(model,strMod3)) then !age stored in linear scale, not as log
				t_iso=Isoc(:,clogt) 
				!the index should be ct. It's still called clogt because of the input of the function
				logt_iso=dnint(log10(t_iso)*100)/100
				t_iso=10.**logt_iso !in this way t_iso comes from logt having 0 or 5 as 2nd decimal digit
			else
				logt_iso=dnint(Isoc(:,clogt)*100)/100
				t_iso=10.**logt_iso
			end if
			allocate(M_iso(size(Isoc,1))); allocate(logL_iso(size(Isoc,1)))
			allocate(L_iso(size(Isoc,1))); allocate(logTeff_iso(size(Isoc,1)))
			allocate(Teff_iso(size(Isoc,1))); allocate(logg_iso(size(Isoc,1)))
			allocate(g_iso(size(Isoc,1))); allocate(R_iso(size(Isoc,1)))
			allocate(rho_iso(size(Isoc,1))); allocate(logrho_iso(size(Isoc,1)))
			 
			M_iso=Isoc(:,cM)
			logL_iso=Isoc(:,clogL)
			L_iso=10.**(logL_iso)
			logTeff_iso=Isoc(:,clogTe)
			Teff_iso=10.**Isoc(:,clogTe)
			logg_iso=Isoc(:,clogg)
			g_iso=10.**Isoc(:,clogg)
			R_iso=sqrt(L_iso/(Teff_iso/TeffSun)**4) !raggi solari
			rho_iso=M_iso/(R_iso**3)*rhoSun !g/cm3
			logrho_iso=log10(rho_iso)
			
			allocate(V_iso(size(Isoc,1))); allocate(BmV_iso(size(Isoc,1)))
			allocate(BC_iso(size(Isoc,1))) 
			V_iso=Isoc(:,cmag0)
			if (isEq(useColor,1.D0,2)) then
				BmV_iso=Isoc(:,cmag2)-Isoc(:,cmag1) !Isoc(:,9)-Isoc(:,11);! !10,11:B-V ->PDisoc 9,11:J-K ->2MASSisoc
			else
				BmV_iso=logTeff_iso
			end if
			!!! 
			BC_iso=Isoc(:,cmbol)-Isoc(:,cmag0)
			!!! 
			

			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			!!!!!!!NEW CODE 12/2/14!!!!!!!!!!!!!!! 
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
			!size equal to any of _iso
			allocate(x_iso(size(logTeff_iso))); allocate(y_iso(size(logTeff_iso)))
			allocate(xx_iso(size(logTeff_iso)))
			if (calibHRD) then 
				x_iso=BmV_iso
				xx_iso=logTeff_iso
				y_iso=V_iso
			end if 
			if (calibNoD) then
				if (caliblogL) then
				 	y_iso=logL_iso
				else if (loggAvail) then 
					y_iso=logg_iso
				else !rhoAvail
					y_iso=logrho_iso
				end if 
				x_iso=BmV_iso
				xx_iso=logTeff_iso
			end if 
			if (calibSPEC) then
				if (caliblogL) then
					y_iso=logL_iso
				else if (loggAvail) then 
					y_iso=logg_iso
				else !rhoAvail
					y_iso=logrho_iso
				end if 
				x_iso=logTeff_iso
				xx_iso=BmV_iso
			end if 
			
			allocate(dist_ndx(size(x_iso))); allocate(ageMinDist(size(logt_iso)))
			call indexx(sqrt((x-x_iso)**2+(y-y_iso)**2),dist_ndx)
			ageMinDist=logt_iso(dist_ndx)
			call uniqueFast(ageMinDist,2,ix,.true.)
			allocate(dist_ndxU(size(ix)))
			dist_ndxU=dist_ndx(ix)
			
			deallocate(dist_ndx); deallocate(ageMinDist)
			deallocate(ix); deallocate(y_iso)
			
			allocate(t_iTmp(size(dist_ndxU)));allocate(distTmp(size(dist_ndxU),2))
			allocate(logTeff_iTmp(size(dist_ndxU)));allocate(logL_iTmp(size(dist_ndxU)))
			allocate(M_iTmp(size(dist_ndxU)));allocate(logg_iTmp(size(dist_ndxU)))
			allocate(logrho_iTmp(size(dist_ndxU)))
			allocate(BmV_iTmp(size(dist_ndxU)));allocate(BC_iTmp(size(dist_ndxU)))
			
			rowAge=1
			do Nage=1,size(dist_ndxU,1)
				ndxDU=dist_ndxU(Nage)
				BmV_i1=BmV_iso(ndxDU)
				V_i1=V_iso(ndxDU)
				
				logTeff_i1=logTeff_iso(ndxDU)
				logL_i1=logL_iso(ndxDU)
				M_i1=M_iso(ndxDU)
				logg_i1=logg_iso(ndxDU)
				logrho_i1=logrho_iso(ndxDU)
				
				call choose_i2Col(ndxDU,x,x_iso,xx_iso,calibSPEC,hstar,hstarlim,V_iso,logL_iso, &
				& M_iso,logg_iso,logrho_iso,logt_iso,BmV_i2,V_i2,logTeff_i2,logL_i2,M_i2,logg_i2,logrho_i2)
				if (calibSPEC) then
					x_i1=logTeff_i1
					x_i2=logTeff_i2
				else
					x_i1=BmV_i1
					x_i2=BmV_i2
				end if
				if (isEq(x_i1,x_i2,4)) then 
					cycle
				end if
				
				t_iTmp(rowAge)=t_iso(ndxDU)
				BC_iTmp(rowAge)=Isoc(ndxDU,cmbol)-Isoc(ndxDU,cmag0) !just rough. 
				!When I recover BC in the logg-BmV plane, I'll make the interpolation using BC_i1 e BC_i2
				!!!
				if (calibHRD) then
					call findX_i(x,y,BmV_i1,BmV_i2,V_i1,V_i2,x_i,dist0)
					BmV_iTmp(rowAge)=x_i
				end if
				if (calibNoD) then
					if (caliblogL) then
						call findX_i(x,y,BmV_i1,BmV_i2,logL_i1,logL_i2,x_i,dist0)
					else if (loggAvail) then 
						call findX_i(x,y,BmV_i1,BmV_i2,logg_i1,logg_i2,x_i,dist0)
					else !rhoAvail
						call findX_i(x,y,BmV_i1,BmV_i2,logrho_i1,logrho_i2,x_i,dist0)
					end if 
					BmV_iTmp(rowAge)=x_i
				end if
				if (calibSPEC) then
					if (caliblogL) then
						call findX_i(x,y,logTeff_i1,logTeff_i2,logL_i1,logL_i2,x_i,dist0)
					else if (loggAvail) then 
						call findX_i(x,y,logTeff_i1,logTeff_i2,logg_i1,logg_i2,x_i,dist0)
					else !rhoAvail
						call findX_i(x,y,logTeff_i1,logTeff_i2,logrho_i1,logrho_i2,x_i,dist0)
					end if 
					logTeff_iTmp(rowAge)=x_i
				end if
				
				distTmp(rowAge,:)=(/ dist0, logt_iso(ndxDU) /)
				
				!!18/3/15 The following if cancel the criterion of perpendicular distance
				!!between a star and the isochrone if the theoretical point to be interpolated
				!!doesn't fall inside the isochrone (i.e. between x_i1 and x_i2), but on its
				!!extension. This possibility, in fact, would build fictitious isochrones that
				!!could be erroneously close to the star
				if (.not.((x_i>=x_i1 .and. x_i<=x_i2) .or. (x_i>=x_i2 .and. x_i<=x_i1))) then 
					BmV_iTmp(rowAge)=BmV_i1
					logg_iTmp(rowAge)=logg_i1
					logTeff_iTmp(rowAge)=logTeff_i1
					logL_iTmp(rowAge)=logL_i1
					M_iTmp(rowAge)=M_i1
					logrho_iTmp(rowAge)=logrho_i1
				else
					if (calibSPEC) then 
						mBV=(BmV_i2-BmV_i1)/(x_i2-x_i1)
						qBV=-mBV*x_i1+BmV_i1
						BmV_iTmp(rowAge)=mBV*x_i+qBV
					else
						mT=(logTeff_i2-logTeff_i1)/(x_i2-x_i1)
						qT=-mT*x_i1+logTeff_i1
						logTeff_iTmp(rowAge)=mT*x_i+qT
					end if 
					mg=(logg_i2-logg_i1)/(x_i2-x_i1)
					qg=-mg*x_i1+logg_i1
					logg_iTmp(rowAge)=mg*x_i+qg
					mL=(logL_i2-logL_i1)/(x_i2-x_i1)
					qL=-mL*x_i1+logL_i1
					logL_iTmp(rowAge)=mL*x_i+qL
					mM=(M_i2-M_i1)/(x_i2-x_i1)
					qM=-mM*x_i1+M_i1
					M_iTmp(rowAge)=mM*x_i+qM
					mrho=(logrho_i2-logrho_i1)/(x_i2-x_i1)
					qrho=-mrho*x_i1+logrho_i1
					logrho_iTmp(rowAge)=mrho*x_i+qrho
				end if 

				rowAge=rowAge+1
			end do
			deallocate(t_iso); deallocate(logt_iso); deallocate(M_iso)
			deallocate(logL_iso); deallocate(L_iso); deallocate(logTeff_iso)
			deallocate(Teff_iso); deallocate(logg_iso); deallocate(g_iso)
			deallocate(R_iso); deallocate(rho_iso); deallocate(logrho_iso)
			deallocate(V_iso); deallocate(BmV_iso); deallocate(BC_iso)
			deallocate(x_iso); deallocate(xx_iso)
			
			allocate(t_i(rowAge-1));allocate(dist(rowAge-1,2))
			allocate(logTeff_i(rowAge-1));allocate(logL_i(rowAge-1))
			allocate(M_i(rowAge-1));allocate(logg_i(rowAge-1))
			allocate(logrho_i(rowAge-1))
			allocate(BmV_i(rowAge-1));allocate(BC_i(rowAge-1))
			t_i=t_iTmp(1:rowAge-1)
			dist=distTmp(1:rowAge-1,:)
			logTeff_i=logTeff_iTmp(1:rowAge-1)
			logL_i=logL_iTmp(1:rowAge-1)
			M_i=M_iTmp(1:rowAge-1)
			logg_i=logg_iTmp(1:rowAge-1)
			logrho_i=logrho_iTmp(1:rowAge-1)
			BmV_i=BmV_iTmp(1:rowAge-1)
			BC_i=BC_iTmp(1:rowAge-1)
			deallocate(t_iTmp);deallocate(distTmp);deallocate(logTeff_iTmp)
			deallocate(logL_iTmp);deallocate(M_iTmp);deallocate(logg_iTmp)
			deallocate(logrho_iTmp);deallocate(BmV_iTmp);deallocate(BC_iTmp)
			
			allocate(Teff_i(1:rowAge-1));allocate(L_i(1:rowAge-1))
			allocate(g_i(1:rowAge-1));allocate(rho_i(1:rowAge-1))
			Teff_i=10.**logTeff_i
			L_i=10.**logL_i
			g_i=10.**logg_i
			rho_i=10.**logrho_i
			
			deallocate(dist_ndxU)
			
			allocate(Isoc_i(rowAge-1,13)) !11 nel T
			do j=1,rowAge-1
				Isoc_i(j,:)=(/Z_input,t_i(j),logTeff_i(j),Teff_i(j),logL_i(j),L_i(j),logg_i(j),g_i(j), &
							& M_i(j),logrho_i(j),rho_i(j),BmV_i(j),BC_i(j) /)
			end do
			!11 nel T
						
			deallocate(t_i); deallocate(dist); deallocate(BC_i); deallocate(BmV_i)
			deallocate(logTeff_i); deallocate(logL_i); deallocate(M_i)
			deallocate(logg_i); deallocate(logrho_i); deallocate(Teff_i)
			deallocate(L_i); deallocate(g_i); deallocate(rho_i)
			
			!End pre-while 
			n=n+1
			Z_input2n=Z_low*qgp**n
			indxZinput=minloc(abs(Z_input2n-Z_iso),1)
			Z_input2=Z_iso(indxZinput);
			do while (abs(Z_input2-Z_input)<1.e-6 .and. Z_input2n<=Z_isoUP*(1+qgp) .and. n<=n_is) 
			!Z to be loaded is equal to the previous. Can happen if Z is widely spaced like in LG
				n=n+1
				Z_input2n=Z_low*qgp**n
				indxZinput=minloc(abs(Z_input2n-Z_iso),1)
				Z_input2=Z_iso(indxZinput)
			end do
			Z_input=Z_input2
				
			do while (Z_input2n<=Z_isoUP*(1+qgp) .and. n<=n_is) !Z_input2nn<=Z_isoUP && Z_input<=Z_isoUP
				if (strcmp(model,strMod1)) then
					call num2str(Z_input,5,Z_inputS)
					fn1=trim(percorso)//trim(Z_inputS)//"M.txt"
					call loadMatrix(fn1,Isoc1,head)
				else if (strcmp(model,strMod2)) then
					write(Z_inputSLG,'(F5.3)') Z_input 
					fn2=trim(percorso)//Ystr//"Z"//Z_inputSLG//".txt"
					call loadMatrix(fn2,Isoc1,head)
				else if (strcmp(model,strMod3)) then
					call num2str(Z_input,5,Z_inputS)
					fn3=trim(percorso)//trim(Z_inputS)//".dat"
					call loadMatrix(fn3,Isoc1,head)
				end if
				!!Check whether logt has equally spaced steps
				allocate(logtiso(size(Isoc1,1)))
				if (strcmp(model,strMod3)) then 
					logtiso=log10(Isoc1(:,clogt))
				else
					logtiso=Isoc1(:,clogt)
				end if 
				if (.not.isEq(logt8c,0.D0,2)) then !Activity check done
					ilogt8=1
					do while (abs(logt8-Isoc1(ilogt8,clogt))>Inct/2. .and. ilogt8<size(Isoc1,1)) 
						ilogt8=ilogt8+1
					end do 
					allocate(IsocTmp(size(Isoc1,1)-ilogt8+1,size(Isoc1,2)))
					IsocTmp=Isoc1(ilogt8:size(Isoc1,1),:)
					deallocate(Isoc1)
					allocate(Isoc1(size(IsocTmp,1),size(IsocTmp,2)))
					Isoc1=IsocTmp
					deallocate(IsocTmp)
				end if
				allocate(FeHcol(size(Isoc1,1)))
				FeHcol=(log10(Z_input)+costZ)*ones(size(Isoc1,1))
				allocate(IsocTmp(size(Isoc1,1),size(Isoc1,2)+1))
				IsocTmp(:,1:size(IsocTmp,2)-1)=Isoc1
				IsocTmp(:,size(IsocTmp,2))=FeHcol
				deallocate(Isoc1); deallocate(FeHcol)
				allocate(Isoc1(size(IsocTmp,1),size(IsocTmp,2)))
				Isoc1=IsocTmp
				deallocate(IsocTmp)

				allocate(t_iso(size(Isoc1,1))); allocate(logt_iso(size(Isoc1,1)))
				if (strcmp(model,strMod3)) then !age stored in linear scale, not as log
					t_iso=Isoc1(:,clogt) 
					!the index should be ct. It's still caled clogt because of the input of the function
					logt_iso=dnint(log10(t_iso)*100)/100
					t_iso=10.**logt_iso !in this way t_iso comes from logt having 0 or 5 as 2nd decimal digit
				else
					logt_iso=dnint(Isoc1(:,clogt)*100)/100
					t_iso=10.**logt_iso
				end if
				allocate(M_iso(size(Isoc1,1))); allocate(logL_iso(size(Isoc1,1)))
				allocate(L_iso(size(Isoc1,1))); allocate(logTeff_iso(size(Isoc1,1)))
				allocate(Teff_iso(size(Isoc1,1))); allocate(logg_iso(size(Isoc1,1)))
				allocate(g_iso(size(Isoc1,1))); allocate(R_iso(size(Isoc1,1)))
				allocate(rho_iso(size(Isoc1,1))); allocate(logrho_iso(size(Isoc1,1)))
				 
				M_iso=Isoc1(:,cM)
				logL_iso=Isoc1(:,clogL)
				L_iso=10.**(logL_iso)
				logTeff_iso=Isoc1(:,clogTe)
				Teff_iso=10.**Isoc1(:,clogTe)
				logg_iso=Isoc1(:,clogg)
				g_iso=10.**Isoc1(:,clogg)
				R_iso=sqrt(L_iso/(Teff_iso/TeffSun)**4) !raggi solari
				rho_iso=M_iso/(R_iso**3)*rhoSun !g/cm3
				logrho_iso=log10(rho_iso)
				
				allocate(V_iso(size(Isoc1,1))); allocate(BmV_iso(size(Isoc1,1)))
				allocate(BC_iso(size(Isoc1,1))) 
				V_iso=Isoc1(:,cmag0)
				if (isEq(useColor,1.D0,2)) then
					BmV_iso=Isoc1(:,cmag2)-Isoc1(:,cmag1) !Isoc(:,9)-Isoc(:,11);! !10,11:B-V ->PDisoc 9,11:J-K ->2MASSisoc
				else
					BmV_iso=logTeff_iso
				end if
				!!! 
				BC_iso=Isoc1(:,cmbol)-Isoc1(:,cmag0)
				!!! 
				
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
				!!!!!!!NEW CODE 12/2/14!!!!!!!!!!!!!!! 
				!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
				!size equal to any of _iso
				allocate(x_iso(size(logTeff_iso))); allocate(y_iso(size(logTeff_iso)))
				allocate(xx_iso(size(logTeff_iso)))
				if (calibHRD) then 
					x_iso=BmV_iso
					xx_iso=logTeff_iso
					y_iso=V_iso
				end if 
				if (calibNoD) then
					if (caliblogL) then
						y_iso=logL_iso 	
					else if (loggAvail) then 
						y_iso=logg_iso
					else !rhoAvail
						y_iso=logrho_iso
					end if 
					x_iso=BmV_iso
					xx_iso=logTeff_iso
				end if 
				if (calibSPEC) then
					if (caliblogL) then
						y_iso=logL_iso 
					else if (loggAvail) then 
						y_iso=logg_iso
					else !rhoAvail
						y_iso=logrho_iso
					end if 
					x_iso=logTeff_iso
					xx_iso=BmV_iso
				end if 
			
				allocate(dist_ndx(size(x_iso))); allocate(ageMinDist(size(logt_iso)))
				call indexx(sqrt((x-x_iso)**2+(y-y_iso)**2),dist_ndx)
				ageMinDist=logt_iso(dist_ndx)
				call uniqueFast(ageMinDist,2,ix,.true.)
				allocate(dist_ndxU(size(ix)))
				dist_ndxU=dist_ndx(ix)
			
				deallocate(dist_ndx); deallocate(ageMinDist)
				deallocate(ix); deallocate(y_iso)
			
				allocate(t_iTmp(size(dist_ndxU)));allocate(distTmp(size(dist_ndxU),2))
				allocate(logTeff_iTmp(size(dist_ndxU)));allocate(logL_iTmp(size(dist_ndxU)))
				allocate(M_iTmp(size(dist_ndxU)));allocate(logg_iTmp(size(dist_ndxU)))
				allocate(logrho_iTmp(size(dist_ndxU)))
				allocate(BmV_iTmp(size(dist_ndxU)));allocate(BC_iTmp(size(dist_ndxU)))
			
				rowAge=1
				do Nage=1,size(dist_ndxU,1) 
					ndxDU=dist_ndxU(Nage)
					BmV_i1=BmV_iso(ndxDU)
					V_i1=V_iso(ndxDU)
				
					logTeff_i1=logTeff_iso(ndxDU)
					logL_i1=logL_iso(ndxDU)
					M_i1=M_iso(ndxDU)
					logg_i1=logg_iso(ndxDU)
					logrho_i1=logrho_iso(ndxDU)
				
					call choose_i2Col(ndxDU,x,x_iso,xx_iso,calibSPEC,hstar,hstarlim,V_iso,logL_iso, &
					& M_iso,logg_iso,logrho_iso,logt_iso,BmV_i2,V_i2,logTeff_i2,logL_i2,M_i2,logg_i2,logrho_i2)
					if (calibSPEC) then
						x_i1=logTeff_i1
						x_i2=logTeff_i2
					else
						x_i1=BmV_i1
						x_i2=BmV_i2
					end if
					
					if (isEq(x_i1,x_i2,4)) then 
						cycle
					end if
				
					t_iTmp(rowAge)=t_iso(ndxDU)
					BC_iTmp(rowAge)=Isoc1(ndxDU,cmbol)-Isoc1(ndxDU,cmag0) !just rough. 
					!When I recover BC in the logg-BmV plane, I'll make the interpolation using BC_i1 e BC_i2
					!!!
					if (calibHRD) then
						call findX_i(x,y,BmV_i1,BmV_i2,V_i1,V_i2,x_i,dist0)
						BmV_iTmp(rowAge)=x_i
					end if
					if (calibNoD) then
						if (caliblogL) then
							call findX_i(x,y,BmV_i1,BmV_i2,logL_i1,logL_i2,x_i,dist0)
						else if (loggAvail) then 
							call findX_i(x,y,BmV_i1,BmV_i2,logg_i1,logg_i2,x_i,dist0)
						else !rhoAvail
							call findX_i(x,y,BmV_i1,BmV_i2,logrho_i1,logrho_i2,x_i,dist0)
						end if 
						BmV_iTmp(rowAge)=x_i
					end if
					if (calibSPEC) then
						if (caliblogL) then
							call findX_i(x,y,logTeff_i1,logTeff_i2,logL_i1,logL_i2,x_i,dist0)
						else if (loggAvail) then 
							call findX_i(x,y,logTeff_i1,logTeff_i2,logg_i1,logg_i2,x_i,dist0)
						else !rhoAvail
							call findX_i(x,y,logTeff_i1,logTeff_i2,logrho_i1,logrho_i2,x_i,dist0)
						end if 
						logTeff_iTmp(rowAge)=x_i
					end if
				
					distTmp(rowAge,:)=(/ dist0, logt_iso(ndxDU) /)
				
					!!18/3/15 The following if cancel the criterion of perpendicular distance
					!!between a star and the isochrone if the theoretical point to be interpolated
					!!doesn't fall inside the isochrone (i.e. between x_i1 and x_i2), but on its
					!!extension. This possibility, in fact, would build fictitious isochrones that
					!!could be erroneously close to the star
					if (.not.((x_i>=x_i1 .and. x_i<=x_i2) .or. (x_i>=x_i2 .and. x_i<=x_i1))) then 
						BmV_iTmp(rowAge)=BmV_i1
						logg_iTmp(rowAge)=logg_i1
						logTeff_iTmp(rowAge)=logTeff_i1
						logL_iTmp(rowAge)=logL_i1
						M_iTmp(rowAge)=M_i1
						logrho_iTmp(rowAge)=logrho_i1
					else
						if (calibSPEC) then 
							mBV=(BmV_i2-BmV_i1)/(x_i2-x_i1)
							qBV=-mBV*x_i1+BmV_i1
							BmV_iTmp(rowAge)=mBV*x_i+qBV
						else
							mT=(logTeff_i2-logTeff_i1)/(x_i2-x_i1)
							qT=-mT*x_i1+logTeff_i1
							logTeff_iTmp(rowAge)=mT*x_i+qT
						end if 
						mg=(logg_i2-logg_i1)/(x_i2-x_i1)
						qg=-mg*x_i1+logg_i1
						logg_iTmp(rowAge)=mg*x_i+qg
						mL=(logL_i2-logL_i1)/(x_i2-x_i1)
						qL=-mL*x_i1+logL_i1
						logL_iTmp(rowAge)=mL*x_i+qL
						mM=(M_i2-M_i1)/(x_i2-x_i1)
						qM=-mM*x_i1+M_i1
						M_iTmp(rowAge)=mM*x_i+qM
						mrho=(logrho_i2-logrho_i1)/(x_i2-x_i1)
						qrho=-mrho*x_i1+logrho_i1
						logrho_iTmp(rowAge)=mrho*x_i+qrho
					end if 

					rowAge=rowAge+1
				end do
				deallocate(Isoc1)
				deallocate(t_iso); deallocate(logt_iso); deallocate(M_iso)
				deallocate(logL_iso); deallocate(L_iso); deallocate(logTeff_iso)
				deallocate(Teff_iso); deallocate(logg_iso); deallocate(g_iso)
				deallocate(R_iso); deallocate(rho_iso); deallocate(logrho_iso)
				deallocate(V_iso); deallocate(BmV_iso); deallocate(BC_iso)
				deallocate(x_iso); deallocate(xx_iso)
			
				allocate(t_i(rowAge-1));allocate(dist(rowAge-1,2))
				allocate(logTeff_i(rowAge-1));allocate(logL_i(rowAge-1))
				allocate(M_i(rowAge-1));allocate(logg_i(rowAge-1))
				allocate(logrho_i(rowAge-1))
				allocate(BmV_i(rowAge-1));allocate(BC_i(rowAge-1))
				t_i=t_iTmp(1:rowAge-1)
				dist=distTmp(1:rowAge-1,:)
				logTeff_i=logTeff_iTmp(1:rowAge-1)
				logL_i=logL_iTmp(1:rowAge-1)
				M_i=M_iTmp(1:rowAge-1)
				logg_i=logg_iTmp(1:rowAge-1)
				logrho_i=logrho_iTmp(1:rowAge-1)
				BmV_i=BmV_iTmp(1:rowAge-1)
				BC_i=BC_iTmp(1:rowAge-1)
				deallocate(t_iTmp);deallocate(distTmp);deallocate(logTeff_iTmp)
				deallocate(logL_iTmp);deallocate(M_iTmp);deallocate(logg_iTmp)
				deallocate(logrho_iTmp);deallocate(BmV_iTmp);deallocate(BC_iTmp)
			
				allocate(Teff_i(1:rowAge-1));allocate(L_i(1:rowAge-1))
				allocate(g_i(1:rowAge-1));allocate(rho_i(1:rowAge-1))
				Teff_i=10.**logTeff_i
				L_i=10.**logL_i
				g_i=10.**logg_i
				rho_i=10.**logrho_i
			
				deallocate(dist_ndxU)
			
				allocate(Isoc_i1(rowAge-1,13)) !11 nel T
				do j=1,rowAge-1
					Isoc_i1(j,:)=(/Z_input,t_i(j),logTeff_i(j),Teff_i(j),logL_i(j),L_i(j),logg_i(j),g_i(j), &
								& M_i(j),logrho_i(j),rho_i(j),BmV_i(j),BC_i(j) /)
				end do
				!11 nel T
						
				deallocate(t_i); deallocate(dist); deallocate(BC_i); deallocate(BmV_i)
				deallocate(logTeff_i); deallocate(logL_i); deallocate(M_i)
				deallocate(logg_i); deallocate(logrho_i); deallocate(Teff_i)
				deallocate(L_i); deallocate(g_i); deallocate(rho_i)
				
				!End 
				
				!!!
				call append(Isoc_i,Isoc_i1)
				!!!
				deallocate(Isoc_i1)
				
			!!!!!!!!!!!!!!!!!!!!!!!!!!
			!!!!!!!!!!!!!!!!!!!
				n=n+1
				Z_input2n=Z_low*qgp**n
				indxZinput=minloc(abs(Z_input2n-Z_iso),1)
				Z_input2=Z_iso(indxZinput);
				do while (abs(Z_input2-Z_input)<1.e-6 .and. Z_input2n<=Z_isoUP*(1+qgp) .and. n<=n_is) 
				!Z to be loaded is equal to the previous. Can happen if Z is widely spaced like in LG
					n=n+1
					Z_input2n=Z_low*qgp**n
					indxZinput=minloc(abs(Z_input2n-Z_iso),1)
					Z_input2=Z_iso(indxZinput)
				end do
				Z_input=Z_input2
				
			!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			!!!!!!!!!!!!!!!!!!! 
			end do !end while
		end if 
	end if !If not entered this if, Isoc_i will remain deallocated
	
	if (strcmp(model,strMod3)) then 
		logt8=log10(logt8) !Re-set the logt8 as logarithm for the output of Stelle
	end if

end subroutine multipleZisocPhSCP

subroutine append(A,B)
	IMPLICIT NONE
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: A
	DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: B
	
	DOUBLE PRECISION, DIMENSION(size(A,1)+size(B,1),size(A,2)) :: Atmp
	
	Atmp(1:size(A,1),:)=A
	Atmp(size(A,1)+1:size(Atmp,1),:)=B
	deallocate(A)
	allocate(A(size(Atmp,1),size(Atmp,2)))
	A=Atmp

end subroutine append

subroutine append1D(A,B)
	IMPLICIT NONE
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE, INTENT(inout) :: A
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: B
	
	DOUBLE PRECISION, DIMENSION(size(A)+size(B)) :: Atmp
	
	Atmp(1:size(A))=A
	Atmp(size(A)+1:size(Atmp))=B
	deallocate(A)
	allocate(A(size(Atmp)))
	A=Atmp

end subroutine append1D

subroutine M2R(logt_iso,mid_n,Isoc,Mi,caliblogL,loggAvail,loggAvailCal,logY,Rv) 
	!First the isochrone whose logt is represented by the mid_n index is chosen
	!Then Radius Rv is the value reported on that line of the isochrone where mass is ~Mi
	!Thus this routine gives the correspondence between stellar mass and stellar radius and
	!it is useful for low MS stars (very low Mass stars) for which R=R(t) is almost constant
	!Actually this subroutine has been tested in this case and assume low logY
	USE SunPD, ONLY: TeffSun,rhoSun
	USE IsoPD
	
	IMPLICIT NONE
	DOUBLE PRECISION Mi,logY,Rv
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: logt_iso
	DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: Isoc
	INTEGER mid_n
	LOGICAL caliblogL,loggAvail,loggAvailCal
	
	DOUBLE PRECISION logYv,xlow,xup
	DOUBLE PRECISION, DIMENSION(2) :: ylow,yup
	DOUBLE PRECISION, DIMENSION(2) :: MR05
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: R05,logrho05
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: mid_Isoc,mid_Isoc05,mid_Isoc05rho
	INTEGER xM,cY
	INTEGER, DIMENSION(:), ALLOCATABLE :: ndxit,ndxft
		
	interface
		subroutine trovaIndici(v,ndxi,ndxf)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: v
			INTEGER, DIMENSION(:), ALLOCATABLE :: ndxi,ndxf
		end subroutine trovaIndici
		subroutine InterpLin_M(fileM,x,colx,OutRange,colY,y,xlow,ylow,xup,yup)
			IMPLICIT NONE
			! Arguments declarations 
			DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: fileM
			DOUBLE PRECISION x,xlow,xup
			INTEGER, DIMENSION(:), INTENT(in) :: colY
			DOUBLE PRECISION, DIMENSION(size(colY)) :: ylow,yup,y
			INTEGER colx,OutRange
		end subroutine InterpLin_M
	end interface
	
	call trovaIndici(logt_iso,ndxit,ndxft) 
	
	allocate(mid_Isoc(ndxft(mid_n)-ndxit(mid_n)+1,size(Isoc,2)))
	mid_Isoc=Isoc(ndxit(mid_n):ndxft(mid_n),:)
	xM=1
	do while (mid_Isoc(xM,cM)<1.11*Mi)
		xM=xM+1
	end do
	allocate(mid_Isoc05(xM-1,size(mid_Isoc,2)))
	allocate(R05(size(mid_Isoc05,1))); allocate(logrho05(size(mid_Isoc05,1))) 
	allocate(mid_Isoc05rho(size(mid_Isoc05,1),size(mid_Isoc05,2)+2))
	
	mid_Isoc05=mid_Isoc(1:xM-1,:)
	R05=sqrt(10.**mid_Isoc05(:,clogL)/(10.**mid_Isoc05(:,clogTe)/TeffSun)**4) !solRadii
	logrho05=log10(mid_Isoc05(:,cM)/R05**3*rhoSun) !g/cm3
	mid_Isoc05rho(:,1:size(mid_Isoc05,2))=mid_Isoc05
	mid_Isoc05rho(:,size(mid_Isoc05,2)+1)=R05
	mid_Isoc05rho(:,size(mid_Isoc05,2)+2)=logrho05
	if (caliblogL) then
		cY=clogL
		logYv=logY
	else if (loggAvail.or.loggAvailCal) then 
		cY=clogg
		logYv=-logY !re-turn to correct value
	else
		cY=size(mid_Isoc05rho,2)
		logYv=-logY !re-turn to correct value
	end if 
	call InterpLin_M(mid_Isoc05rho,logYv,cY,1,(/ cM,size(mid_Isoc05rho,2)-1 /),MR05,xlow,ylow,xup,yup)

	Rv=MR05(2)

end subroutine M2R

subroutine searchTOold(y_iso,logTeff_iso,caliblogL,logt_iso,last_logt,logY_soglia)
	IMPLICIT NONE
	DOUBLE PRECISION last_logt,logY_soglia
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: y_iso,logTeff_iso,logt_iso
	LOGICAL caliblogL
	
	DOUBLE PRECISION mOld
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: logY_isoOld,logTeff_isoOld
	INTEGER kto,kOld
	LOGICAL isEq
	
	kto=1
	do while (isEq(logt_iso(size(logt_iso,1)-kto),last_logt,2).and.kto<size(logt_iso,1)) 
		kto=kto+1
	end do 
	kto=kto-1
	allocate(logY_isoOld(kto+1)); allocate(logTeff_isoOld(kto+1))
	logY_isoOld=y_iso(size(logt_iso,1)-kto:size(logt_iso,1))
	logTeff_isoOld=logTeff_iso(size(logt_iso,1)-kto:size(logt_iso,1))
	kOld=1
	do while (isEq(logTeff_isoOld(kOld+1),logTeff_isoOld(kOld),4).and.kOld+1<size(logTeff_isoOld)) 
		kOld=kOld+1
	end do
	mOld=(logY_isoOld(kOld+1)-logY_isoOld(kOld))/(logTeff_isoOld(kOld+1)-logTeff_isoOld(kOld))
	if (.not.caliblogL) then !19/9/2017. logL not defined. I'm either in the plane logg-logTeff or in logrho-logTeff
		mOld=-mOld !In both the planes MS is decreasing.
		! Changing here the sign, I do not change inequality mOld>0 that looks for increasing isoc
	end if 
	do while (mOld>=0) 
		kOld=kOld+1
		do while (isEq(logTeff_isoOld(kOld+1),logTeff_isoOld(kOld),4) .and. kOld+1<size(logTeff_isoOld)) 
			kOld=kOld+1
		end do
		mOld=(logY_isoOld(kOld+1)-logY_isoOld(kOld))/(logTeff_isoOld(kOld+1)-logTeff_isoOld(kOld))
		if (.not.caliblogL) then !19/9/2017. logL not defined
			mOld=-mOld
		end if 
	end do 
	if (caliblogL) then !19/9/2017
		logY_soglia=logY_isoOld(kOld)-0.2 !decrease arbitrarely of 0.2 to locate a bit under the T.O.
	else !logY=logrho o logY=logg
		logY_soglia=-(logY_isoOld(kOld)+0.2) !With minus to use the inequality logY<logY_soglia,
											 ! which yields for Y=luminosity
	end if 
	
	deallocate(logY_isoOld); deallocate(logTeff_isoOld)
	
end subroutine searchTOold

subroutine calibrateDiag(state,d1lim,x1,x2,x1star,x2star)
	!This subroutine converts one coordinate of the star (x1star e.g. VMag) from a diagram (e.g. CMD)
	!to another diagram (e.g. hybrid logL-BmV diagram) giving the corresponding of x1star that is
	!x2star (in this example logL)  
	!x1, x2: vectors(dim=2), each component represents the point coordinate ON one of the two
	!reference isochrones. In this example:
	!x1(1): Vmag of the 1st isochrone (computed considering the point of the isochrone having the
	!		same BmV of the star
	!x1(2): VMag of the 2nd isochrone (computed...)
	!x2(1): logL of the 1st isochrone (computed...)
	!x2(2): logL of the 2nd isochrone (computed...)
	!
	IMPLICIT NONE
	INTEGER state
	DOUBLE PRECISION, DIMENSION(2) :: x1,x2
	DOUBLE PRECISION d1lim,x1star,x2star
	
	DOUBLE PRECISION, DIMENSION(2) :: d1star_is
	DOUBLE PRECISION r,d1_is,d2_is,dr
	!dr: fraction of d2_is that expresses how much the star is far from an isochrone in the 2nd diagram
	
	d1_is=abs(x1(1)-x1(2)) !distance between the two points on the two reference isochrones
						   !taken in the 1st diagram (e.g. CMD). This distance is always taken
						   !vertically (as in this example, being BmV fixed) or horizontally
	d2_is=abs(x2(1)-x2(2))
	d1star_is=abs(x1star-x1)	!distance between the star and the two reference isochrones in the
								!1st diagram (e.g. CMD). This distance is always taken vertically
								!(as in this example) or horizontally
	select case (state) !
	case (1) !
		if (x1(1)<x1(2)) then 
			r=d1star_is(1)/(d1star_is(1)+d1_is)
			dr=r*d2_is/(1-r)
			x2star=x2(1)+dr
		else
			r=d1star_is(2)/(d1star_is(2)+d1_is)
			dr=r*d2_is/(1-r)
			x2star=x2(2)+dr
		end if 
	case (2) !
		if (x1(1)<x1(2)) then 
			r=d1star_is(2)/(d1star_is(2)+d1_is)
			dr=r*d2_is/(1-r)
			x2star=x2(2)-dr
		else
			r=d1star_is(1)/(d1star_is(1)+d1_is)
			dr=r*d2_is/(1-r)
			x2star=x2(1)-dr
		end if 
	case (3) !
		if (x1(1)<x1(2)) then 
			if (d1star_is(2)>d1lim) then !should be .neq.0, but I put this lim to avoid numerical errors
				r=d1star_is(1)/d1star_is(2)
				dr=r*d2_is/(1+r)
				x2star=x2(1)-dr
			else
				r=d1star_is(2)/d1star_is(1)
				dr=r*d2_is/(1+r)
				x2star=x2(2)+dr
			end if 
		else
			if (d1star_is(1)>d1lim) then  !!=0
				r=d1star_is(2)/d1star_is(1)
				dr=r*d2_is/(1+r)
				x2star=x2(2)-dr
			else
				r=d1star_is(1)/d1star_is(2)
				dr=r*d2_is/(1+r)
				x2star=x2(1)+dr
			end if 
		end if
	case default
		print*,'Check state: could be either 1, 2 or 3' 
	end select
end subroutine calibrateDiag

subroutine selectIsoc(logt_iso,age,logtstep,last_logt,ti0,tf0)
	IMPLICIT NONE
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: logt_iso
	DOUBLE PRECISION age,logtstep,last_logt
	INTEGER ti0,tf0
	
	INTEGER j_iso
	LOGICAL isEq
	
	do j_iso=1,size(logt_iso) 
		if (isEq(logt_iso(j_iso),age,2)) then
			ti0=j_iso
			exit
		end if 
	end do
	if (isEq(age,last_logt,2)) then
		tf0=size(logt_iso)
	else
		do j_iso=ti0,size(logt_iso) 
			if (isEq(logt_iso(j_iso),age+logtstep,2)) then 
				tf0=j_iso-1
				exit
			end if 
		end do
	end if
	
end subroutine

subroutine findX_i(x,y,x_i1,x_i2,y_i1,y_i2,x_i,dist0)
	IMPLICIT NONE
	DOUBLE PRECISION x,y,x_i1,x_i2,y_i1,y_i2,x_i,dist0
	
	DOUBLE PRECISION mY,q,q1

	mY=(y_i2-y_i1)/(x_i2-x_i1)
	q=-mY*x_i1+y_i1
	dist0=abs(y-(mY*x+q))/sqrt(1+mY**2)
	if (mY.ne.0) then 
		q1=x/mY+y
		x_i=(q1-q)*mY/(mY**2+1)
	else
		x_i=x
	end if 

end subroutine findX_i

subroutine choose_i2xy1y2(x,x_isor,y1_isor,y2_isor,i_jd,hstar,hstarlim, & !y1,
	& x_is1,x_is2,y1_is1,y1_is2,y2_is1,y2_is2) !_is1 _is2 two points on the same isoc
	IMPLICIT NONE
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: x_isor,y1_isor,y2_isor
	DOUBLE PRECISION x,hstar,hstarlim,x_is1,x_is2,y1_is1,y1_is2,y2_is1,y2_is2 !y1,
	INTEGER i_jd
		
	INTEGER step,ks
	LOGICAL isEq,incs
	
	interface
		subroutine calibrate12(x,x_isor,y1_isor,y2_isor,id,ks,incs,x_is1,x_is2,y1_is1,y1_is2,y2_is1,y2_is2)
			IMPLICIT NONE
			DOUBLE PRECISION x,x_is1,x_is2,y1_is1,y1_is2,y2_is1,y2_is2
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: x_isor,y1_isor,y2_isor
			INTEGER id,ks
			LOGICAL incs
		end subroutine calibrate12
	end interface

	if (i_jd==1) then 
		step=1
		do while (isEq(x_is1,x_isor(i_jd+step),4)) 
			step=step+1
		end do 
		x_is2=x_isor(i_jd+step)
		y1_is2=y1_isor(i_jd+step)
		y2_is2=y2_isor(i_jd+step)
	else if (i_jd==size(x_isor)) then 
		step=1
		do while (isEq(x_is1,x_isor(i_jd-step),4)) 
			step=step+1
		end do 
		x_is2=x_isor(i_jd-step)
		y1_is2=y1_isor(i_jd-step)
		y2_is2=y2_isor(i_jd-step)
	else
		if (x>x_is1) then
			incs=.true. 
			if (x_isor(i_jd+1)>x_is1) then !increasing isoc. with respect to x
				ks=1
				call calibrate12(x,x_isor,y1_isor,y2_isor,i_jd,ks,incs,x_is1,x_is2,y1_is1,y1_is2,y2_is1,y2_is2)
!				x_is2=x_isor(i_jd+1)
!				y1_is2=y1_isor(i_jd+1)
!				y2_is2=y2_isor(i_jd+1)
			else if (x_isor(i_jd-1)>x_is1) then !decreasing isoc with respect to x
				ks=-1
				call calibrate12(x,x_isor,y1_isor,y2_isor,i_jd,ks,incs,x_is1,x_is2,y1_is1,y1_is2,y2_is1,y2_is2)
!				x_is2=x_isor(i_jd-1)
!				y1_is2=y1_isor(i_jd-1)
!				y2_is2=y2_isor(i_jd-1)
			else
				if (hstar>hstarlim) then !!low-MS star => go forward with the indices
					step=1
					do while (isEq(x_is1,x_isor(i_jd+step),4)) 
						step=step+1
					end do
					if (x_isor(i_jd+step).gt.x_is1) then
						i_jd=i_jd+step-1
						ks=1
					else
						ks=-1
					end if
					call calibrate12(x,x_isor,y1_isor,y2_isor,i_jd,ks,incs,x_is1,x_is2,y1_is1,y1_is2,y2_is1,y2_is2) 
!					x_is2=x_isor(i_jd+step)
!					y1_is2=y1_isor(i_jd+step)
!					y2_is2=y2_isor(i_jd+step)
				else
					step=1
					do while (isEq(x_is1,x_isor(i_jd-step),4)) 
						step=step+1
					end do
					if (x_isor(i_jd-step).gt.x_is1) then
						i_jd=i_jd-step+1
						ks=-1
					else
						ks=1
					end if
					call calibrate12(x,x_isor,y1_isor,y2_isor,i_jd,ks,incs,x_is1,x_is2,y1_is1,y1_is2,y2_is1,y2_is2)
!					x_is2=x_isor(i_jd-step)
!					y1_is2=y1_isor(i_jd-step)
!					y2_is2=y2_isor(i_jd-step)
				end if 
			end if 
		else
			incs=.false.
			if (x_isor(i_jd-1)<x_is1) then  !!increasing isoc. with respect to x
				ks=-1
				call calibrate12(x,x_isor,y1_isor,y2_isor,i_jd,ks,incs,x_is1,x_is2,y1_is1,y1_is2,y2_is1,y2_is2)
!				x_is2=x_isor(i_jd-1)
!				y1_is2=y1_isor(i_jd-1)
!				y2_is2=y2_isor(i_jd-1)
			else if (x_isor(i_jd+1)<x_is1) then !!decreasing isoc with respect to x
				ks=1
				call calibrate12(x,x_isor,y1_isor,y2_isor,i_jd,ks,incs,x_is1,x_is2,y1_is1,y1_is2,y2_is1,y2_is2)
!				x_is2=x_isor(i_jd+1)
!				y1_is2=y1_isor(i_jd+1)
!				y2_is2=y2_isor(i_jd+1)
			else
				if (hstar>hstarlim) then !!low-MS star => go forward with the indices
					step=1
					do while (isEq(x_is1,x_isor(i_jd+step),4)) 
						step=step+1
					end do
					if (x_isor(i_jd+step).lt.x_is1) then
						i_jd=i_jd+step-1
						ks=1
					else
						ks=-1
					end if
					call calibrate12(x,x_isor,y1_isor,y2_isor,i_jd,ks,incs,x_is1,x_is2,y1_is1,y1_is2,y2_is1,y2_is2)
!					x_is2=x_isor(i_jd+step)
!					y1_is2=y1_isor(i_jd+step)
!					y2_is2=y2_isor(i_jd+step)
				else
					step=1
					do while (isEq(x_is1,x_isor(i_jd-step),4)) 
						step=step+1
					end do
					if (x_isor(i_jd-step).lt.x_is1) then
						i_jd=i_jd-step+1
						ks=-1
					else
						ks=1
					end if
					call calibrate12(x,x_isor,y1_isor,y2_isor,i_jd,ks,incs,x_is1,x_is2,y1_is1,y1_is2,y2_is1,y2_is2)
!					x_is2=x_isor(i_jd-step)
!					y1_is2=y1_isor(i_jd-step)
!					y2_is2=y2_isor(i_jd-step)
				end if 
			end if 
		end if 
	end if
end subroutine choose_i2xy1y2

subroutine calibrate12(x,x_isor,y1_isor,y2_isor,id,ks,incs,x_is1,x_is2,y1_is1,y1_is2,y2_is1,y2_is2)
	IMPLICIT NONE
	DOUBLE PRECISION x,x_is1,x_is2,y1_is1,y1_is2,y2_is1,y2_is2
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: x_isor,y1_isor,y2_isor
	INTEGER id,ks
	LOGICAL incs
	
	DOUBLE PRECISION, DIMENSION(20) :: dX
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dXX
	INTEGER step,idX,step1
	
	x_is1=x_isor(id)
	step=1
	if (ks.gt.0) then
		x_is2=x_isor(id+step)
		dX(step)=abs(x-x_is2)
		do while (x_is1.eq.x_is2 .and. id+step.lt.size(x_isor) .and. step.lt.10)
			step=step+1
			x_is2=x_isor(id+step)
			dX(step)=abs(x-x_is2)
		end do
		if (incs) then !I want to increase because x_is1<x. While also x_is2<x, continue increase x_is2 until x_is2>x, so that x_is1<x<x_is2
			do while (x.gt.x_is2 .and. id+step.lt.size(x_isor) .and. step.lt.10)
				step=step+1
				x_is2=x_isor(id+step)
				dX(step)=abs(x-x_is2)
			end do
		else !x_is1>x
			do while (x.lt.x_is2 .and. id+step.lt.size(x_isor) .and. step.lt.10)
				step=step+1
				x_is2=x_isor(id+step)
				dX(step)=abs(x-x_is2)
			end do
		end if
		if (step.gt.1) then
			if ((incs .and. x.le.x_is2) .or. (.not.incs .and. x.ge.x_is2)) then
				id=id+step-1
				ks=1
			else
				!first try to go backwards to the other direction
				step1=1
				x_is2=x_isor(id-step1)
				dX(step+step1)=abs(x-x_is2)
				do while (x_is1.eq.x_is2 .and. id-step1.gt.1 .and. step1.lt.10)
					step1=step1+1
					x_is2=x_isor(id-step1)
					dX(step+step1)=abs(x-x_is2)
				end do
				if (incs) then !x_is1<x
					do while (x.gt.x_is2 .and. id-step1.gt.1 .and. step1.lt.10)
						step1=step1+1
						x_is2=x_isor(id-step1)
						dX(step+step1)=abs(x-x_is2)
					end do
				else !x_is1>x
					do while (x.lt.x_is2 .and. id-step1.gt.1 .and. step1.lt.10)
						step1=step1+1
						x_is2=x_isor(id-step1)
						dX(step+step1)=abs(x-x_is2)
					end do
				end if
				if ((incs.and.x.le.x_is2) .or. (.not.incs.and.x.ge.x_is2)) then
					id=id-step1+1
					ks=-1
				else
					!!end try block
					allocate(dXX(step+step1))
					dXX=dX(1:step+step1)
					idX=minloc(dXX,1)
					if (idX.le.step) then
						id=id+idX-1
						if (idX.eq.1) then !go forward in dXX, so continue going forward in the _isor
							ks=1
						else if (idX.eq.step) then
							ks=-1
						else
							if (dXX(idX+1).lt.dXX(idX-1)) then
								ks=1 !go forward
							else
								ks=-1
							end if
						end if
					else
						id=id-idX+1
						if (idX.eq.step+1) then
							ks=-1
						else if (idX.eq.size(dXX)) then
							ks=1
						else
							if (dXX(idX+1).lt.dXX(idX-1)) then
								ks=-1 !go backward
							else
								ks=1
							end if
						end if
					end if
				end if
			end if
		end if
	else
		x_is2=x_isor(id-step)
		dX(step)=abs(x-x_is2)
		do while (x_is1.eq.x_is2 .and. id-step.gt.1 .and. step.lt.10)
			step=step+1
			x_is2=x_isor(id-step)
			dX(step)=abs(x-x_is2)
		end do
		if (incs) then !x_is1<x
			do while (x.gt.x_is2 .and. id-step.gt.1 .and. step.lt.10)
				step=step+1
				x_is2=x_isor(id-step)
				dX(step)=abs(x-x_is2)
			end do
		else !x_is1>x
			do while (x.lt.x_is2 .and. id-step.gt.1 .and. step.lt.10)
				step=step+1
				x_is2=x_isor(id-step)
				dX(step)=abs(x-x_is2)
			end do
		end if
		if (step.gt.1) then
			if ((incs .and. x.le.x_is2) .or. (.not.incs .and. x.ge.x_is2)) then
				id=id-step+1
				ks=-1
			else
				!first try to go towards the other direction
				step1=1
				x_is2=x_isor(id+step1)
				dX(step+step1)=abs(x-x_is2)
				do while (x_is1.eq.x_is2 .and. id+step1.lt.size(x_isor) .and. step1.lt.10)
					step1=step1+1
					x_is2=x_isor(id+step1)
					dX(step+step1)=abs(x-x_is2)
				end do
				if (incs) then !I want to increase because x_is1<x. While also x_is2<x, continue increase x_is2
						  !until x_is2>x, so that x_is1<x<x_is2
					do while (x.gt.x_is2 .and. id+step1.lt.size(x_isor) .and. step1.lt.10)
						step1=step1+1
						x_is2=x_isor(id+step1)
						dX(step+step1)=abs(x-x_is2)
					end do
				else !x_is1>x
					do while (x.lt.x_is2 .and. id+step1.lt.size(x_isor) .and. step1.lt.10)
						step1=step1+1
						x_is2=x_isor(id+step1)
						dX(step+step1)=abs(x-x_is2)
					end do
				end if
				if ((incs.and.x.le.x_is2) .or. (.not.incs.and.x.ge.x_is2)) then
					id=id+step1-1
					ks=1
				else
				!!end try block
					allocate(dXX(step+step1))
					dXX=dX(1:step+step1)
					idX=minloc(dXX,1)
					if (idX.le.step) then
						id=id-idX+1
						if (idX.eq.1) then !go forward in dXX, so continue going backward in the _isor
							ks=-1
						else if (idX.eq.step) then
							ks=1
						else
							if (dXX(idX+1).lt.dXX(idX-1)) then
								ks=-1 !go backward
							else
								ks=1
							end if
						end if
					else !idX>step !idX is among the indices added thanks to step1
						id=id+idX-1 !follow ks>0 criterium
						if (idX.eq.step+1) then
							ks=1 !among the last dXX values I was going forward in the _isor (step1 criterium)
						else if (idX.eq.size(dXX)) then
							ks=-1
						else
							if (dXX(idX+1).lt.dXX(idX-1)) then
								ks=1 !go forward
							else
								ks=-1
							end if
						end if
					end if
				end if
			end if
		end if
	end if
	x_is1=x_isor(id)
	y1_is1=y1_isor(id)
	y2_is1=y2_isor(id)
	!!
	x_is2=x_isor(id+ks)
	y1_is2=y1_isor(id+ks)
	y2_is2=y2_isor(id+ks)

end subroutine calibrate12

subroutine choose_i2Col(ndxDU,x,x_iso,xx_iso,calibSPEC,hstar,hstarlim,V_iso,logL_iso, & !x_i1,
	& M_iso,logg_iso,logrho_iso,logt_iso,BmV_i2,V_i2,logTeff_i2,logL_i2,M_i2,logg_i2,logrho_i2)
	
	IMPLICIT NONE
	INTEGER ndxDU
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: x_iso
	DOUBLE PRECISION, DIMENSION(size(x_iso)) :: xx_iso
	DOUBLE PRECISION, DIMENSION(size(x_iso)) :: V_iso,logL_iso
	DOUBLE PRECISION, DIMENSION(size(x_iso)) :: M_iso,logg_iso,logrho_iso,logt_iso
	DOUBLE PRECISION x,BmV_i2,V_i2,logTeff_i2,logL_i2,M_i2,logg_i2,logrho_i2
	DOUBLE PRECISION hstar,hstarlim
	LOGICAL calibSPEC
	
	INTEGER step
	DOUBLE PRECISION, DIMENSION(size(x_iso)) :: logTeff_iso,BmV_iso
	DOUBLE PRECISION x_i1
	LOGICAL isEq
	
	x_i1=x_iso(ndxDU)
	
	if (calibSPEC.eqv..true.) then
		logTeff_iso=x_iso
		BmV_iso=xx_iso
	else
		BmV_iso=x_iso
		logTeff_iso=xx_iso
	endif
		
	if (ndxDU==1) then !1st value of the file
		step=1
		do while (isEq(x_i1,x_iso(ndxDU+step),4))
			step=step+1
		end do 
		! photIsocAvail always true for PD isochrones
		BmV_i2=BmV_iso(ndxDU+step)
		V_i2=V_iso(ndxDU+step)
		
		logTeff_i2=logTeff_iso(ndxDU+step)
		logL_i2=logL_iso(ndxDU+step)
		M_i2=M_iso(ndxDU+step)
		logg_i2=logg_iso(ndxDU+step)
		logrho_i2=logrho_iso(ndxDU+step)
	else if (ndxDU==size(x_iso)) then !last value of the file
		step=1
		do while (isEq(x_i1,x_iso(ndxDU-step),4))
			step=step+1
		end do 
		BmV_i2=BmV_iso(ndxDU-step)
		V_i2=V_iso(ndxDU-step)
		logTeff_i2=logTeff_iso(ndxDU-step)
		logL_i2=logL_iso(ndxDU-step)
		M_i2=M_iso(ndxDU-step)
		logg_i2=logg_iso(ndxDU-step)
		logrho_i2=logrho_iso(ndxDU-step)
	else
		if (.not.isEq(logt_iso(ndxDU),logt_iso(ndxDU-1),2)) then !1st value of an isochrone inside the file
			step=1
			do while (isEq(x_i1,x_iso(ndxDU+step),4)) 
				step=step+1
			end do 
			BmV_i2=BmV_iso(ndxDU+step)
			V_i2=V_iso(ndxDU+step)
			logTeff_i2=logTeff_iso(ndxDU+step)
			logL_i2=logL_iso(ndxDU+step)
			M_i2=M_iso(ndxDU+step)
			logg_i2=logg_iso(ndxDU+step)
			logrho_i2=logrho_iso(ndxDU+step)
		else if (.not.isEq(logt_iso(ndxDU),logt_iso(ndxDU+1),2)) then !last value of an isochrone inside the grid
			step=1
			do while (isEq(x_i1,x_iso(ndxDU-step),4)) 
				step=step+1
			end do 
			BmV_i2=BmV_iso(ndxDU-step)
			V_i2=V_iso(ndxDU-step)
			logTeff_i2=logTeff_iso(ndxDU-step)
			logL_i2=logL_iso(ndxDU-step)
			M_i2=M_iso(ndxDU-step)
			logg_i2=logg_iso(ndxDU-step)
			logrho_i2=logrho_iso(ndxDU-step)
		else
			if (x>x_i1) then 
				if (x_iso(ndxDU+1)>x_i1) then !increasing isoc in BmV (or in Teff 29/4/15)
					step=1
					do while(isEq(x_i1,x_iso(ndxDU+step),4))
						if (isEq(logt_iso(ndxDU),logt_iso(ndxDU+step),2)) then 
							step=step+1
						else
							step=step-1 !in this case BmV_i1=BmV_i2, but it should never enter this else
							exit
						end if 
					end do
					BmV_i2=BmV_iso(ndxDU+step)
					V_i2=V_iso(ndxDU+step)
					logTeff_i2=logTeff_iso(ndxDU+step)
					logL_i2=logL_iso(ndxDU+step)
					M_i2=M_iso(ndxDU+step)
					logg_i2=logg_iso(ndxDU+step)
					logrho_i2=logrho_iso(ndxDU+step)
				else if (x_iso(ndxDU-1)>x_i1) then !decreasing isoc in BmV
					step=1
					do while(isEq(x_i1,x_iso(ndxDU-step),4))
						if (isEq(logt_iso(ndxDU),logt_iso(ndxDU-step),2)) then 
							step=step+1
						else
							step=step-1
							exit
						end if 
					end do
					BmV_i2=BmV_iso(ndxDU-step)
					V_i2=V_iso(ndxDU-step)
					logTeff_i2=logTeff_iso(ndxDU-step)
					logL_i2=logL_iso(ndxDU-step)
					M_i2=M_iso(ndxDU-step)
					logg_i2=logg_iso(ndxDU-step)
					logrho_i2=logrho_iso(ndxDU-step)
				else
					if (hstar>hstarlim) then !Low MS --> go forward with the indices
						step=1
						do while (isEq(x_i1,x_iso(ndxDU+step),4)) 
							if (isEq(logt_iso(ndxDU),logt_iso(ndxDU+step),2)) then 
								step=step+1
							else
								step=step-1 !in this case BmV_i1=BmV_i2, but it should never enter this else
								exit
							end if 
						end do 
						BmV_i2=BmV_iso(ndxDU+step)
						V_i2=V_iso(ndxDU+step)
						logTeff_i2=logTeff_iso(ndxDU+step)
						logL_i2=logL_iso(ndxDU+step)
						M_i2=M_iso(ndxDU+step)
						logg_i2=logg_iso(ndxDU+step)
						logrho_i2=logrho_iso(ndxDU+step)
					else
						step=1
						do while (isEq(x_i1,x_iso(ndxDU-step),4)) 
							if (isEq(logt_iso(ndxDU),logt_iso(ndxDU-step),2)) then 
								step=step+1
							else
								step=step-1
								exit
							end if 
						end do 
						BmV_i2=BmV_iso(ndxDU-step)
						V_i2=V_iso(ndxDU-step)
						logTeff_i2=logTeff_iso(ndxDU-step)
						logL_i2=logL_iso(ndxDU-step)
						M_i2=M_iso(ndxDU-step)
						logg_i2=logg_iso(ndxDU-step)
						logrho_i2=logrho_iso(ndxDU-step)
					end if 
				end if 
			else !x<x_i1
				if (x_iso(ndxDU-1)<x_i1) then !increasing isoc in BmV
					step=1
					do while(isEq(x_i1,x_iso(ndxDU-step),4))
						if (isEq(logt_iso(ndxDU),logt_iso(ndxDU-step),2)) then 
							step=step+1
						else
							step=step-1
							exit
						end if
					end do
					BmV_i2=BmV_iso(ndxDU-step)
					V_i2=V_iso(ndxDU-step)
					logTeff_i2=logTeff_iso(ndxDU-step)
					logL_i2=logL_iso(ndxDU-step)
					M_i2=M_iso(ndxDU-step)
					logg_i2=logg_iso(ndxDU-step)
					logrho_i2=logrho_iso(ndxDU-step)
				else if (x_iso(ndxDU+1)<x_i1) then !decreasing isoc in BmV
					step=1
					do while(isEq(x_i1,x_iso(ndxDU+step),4))
						if (isEq(logt_iso(ndxDU),logt_iso(ndxDU+step),2)) then 
							step=step+1
						else
							step=step-1 !in this case BmV_i1=BmV_i2, but it should never enter this else
							exit
						end if
					end do
					BmV_i2=BmV_iso(ndxDU+step)
					V_i2=V_iso(ndxDU+step)
					logTeff_i2=logTeff_iso(ndxDU+step)
					logL_i2=logL_iso(ndxDU+step)
					M_i2=M_iso(ndxDU+step)
					logg_i2=logg_iso(ndxDU+step)
					logrho_i2=logrho_iso(ndxDU+step)
				else
					if (hstar>hstarlim) then !Low MS => go forward with the indices
						step=1
						do while (isEq(x_i1,x_iso(ndxDU+step),4))
							if (isEq(logt_iso(ndxDU),logt_iso(ndxDU+step),2)) then 
								step=step+1
							else
								step=step-1 !in this case BmV_i1=BmV_i2, but it should never enter this else
								exit
							end if 
						end do 
						BmV_i2=BmV_iso(ndxDU+step)
						V_i2=V_iso(ndxDU+step)
						logTeff_i2=logTeff_iso(ndxDU+step)
						logL_i2=logL_iso(ndxDU+step)
						M_i2=M_iso(ndxDU+step)
						logg_i2=logg_iso(ndxDU+step)
						logrho_i2=logrho_iso(ndxDU+step)
					else
						step=1
						do while (isEq(x_i1,x_iso(ndxDU-step),4)) 
							if (isEq(logt_iso(ndxDU),logt_iso(ndxDU-step),2)) then 
								step=step+1
							else
								step=step-1
								exit
							end if 
						end do 
						BmV_i2=BmV_iso(ndxDU-step)
						V_i2=V_iso(ndxDU-step)
						logTeff_i2=logTeff_iso(ndxDU-step)
						logL_i2=logL_iso(ndxDU-step)
						M_i2=M_iso(ndxDU-step)
						logg_i2=logg_iso(ndxDU-step)
						logrho_i2=logrho_iso(ndxDU-step)
					end if 
				end if 
			end if 
		end if 
	end if
	
end subroutine choose_i2Col

subroutine uniqueFast(list,d,indices,first)
	IMPLICIT NONE
	! find "indices" that correspond to the unique numbers in "list"
	! list(indices) displays the unique elements of list in ascending order
	!!
	! first=.true. For each unique element, select the index
	! of the first occurrence in "list"; otherwise the last is picked.
	! first=.true. is a little bit more time consuming than first=.false.
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: list
	INTEGER d
	INTEGER, DIMENSION(:), ALLOCATABLE :: indices
	LOGICAL first 
	
	INTEGER :: j,i,ext,sL,kx,sI
	INTEGER, DIMENSION(size(list)) :: ndxS,ndxS2
	INTEGER, DIMENSION(:), ALLOCATABLE :: ndxi,ndxf,ndxStmp,itmp,indLast,indFirst
	LOGICAL, DIMENSION(size(list)) :: mask
	DOUBLE PRECISION, DIMENSION(size(list)) :: listS
	
	interface
		subroutine indexx(arr,indx)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: arr !(n)
			INTEGER indx(size(arr))
		end subroutine indexx
		subroutine trovaIndici(v,ndxi,ndxf)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: v
			INTEGER, DIMENSION(:), ALLOCATABLE :: ndxi,ndxf
		end subroutine trovaIndici
		subroutine indexI(arr,indx)
			IMPLICIT NONE
			INTEGER, DIMENSION(:), INTENT(in) :: arr
			INTEGER indx(size(arr))
 		end subroutine indexI
		
		function isEq_vv(x1,x2,n)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: x1,x2
			INTEGER, INTENT(in) :: n
			LOGICAL, DIMENSION(size(x1)) :: isEq_vv
		end function isEq_vv
	end interface
	
	!Order list in ascending order
	call indexx(list,ndxS)
	listS=list(ndxS)
!	print*,'uniqueFast: ndxS',ndxS
	!In case multiple elements are actually present in list, order these multiple
	! elements in the same order they appeared in the original vector 
	!(elements - even if equal - are considered distinguishable)
	call trovaIndici(listS,ndxi,ndxf)
	j=1
	do i=1,size(ndxi)
		ext=ndxf(i)-ndxi(i)+1
		if (ext.gt.1) then
			allocate(ndxStmp(ext)); allocate(itmp(ext))
			ndxStmp=ndxS(ndxi(i):ndxf(i))
			call indexI(ndxStmp,itmp)
			ndxStmp=ndxStmp(itmp)
			ndxS2(j:j+ext-1)=ndxStmp
			deallocate(ndxStmp); deallocate(itmp)
			j=j+ext
		else
			ndxS2(j)=ndxS(ndxi(i))
			j=j+1
		end if
	end do
!	print*,'uniqueFast: ndxS2',ndxS2
	if (first) then
		listS=listS(size(listS):1:-1) !flipped (sorted in descending order)
	end if	
	
	sL=size(listS)
	mask=.false.
	mask(1:sL-1)=isEq_vv(listS(1:sL-1),listS(2:sL),d)
	mask=.not.mask
	
	sI=count(mask)
	allocate(indLast(sI)) !indices(sI)
	indLast=pack((/(kx,kx=1,sL)/),mask)
	
	allocate(indices(size(indLast)))
	if (first) then
		allocate(indFirst(size(indLast)))
		indFirst=sL-indLast+1 !indices of first occurrence in the ascending ordered vector
		indices=ndxS2(indFirst) !would pick up the values of list in descending order
		indices=indices(size(indices):1:-1) !now values of list picked up in ascending order
	else
		indices=ndxS2(indLast)
	end if
	!!list(indices) picks up the unique values of list in ascending order
	!!list(sort(indices)) picks up the values of list in the original sequence they were
	
end subroutine uniqueFast

subroutine findYgivenX_v(M,x,xx,ny,y)
	IMPLICIT NONE
	! Arguments declarations
	DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: M
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: y
	DOUBLE PRECISION, DIMENSION(size(M,1)) :: xx
	DOUBLE PRECISION x
	INTEGER ny 
	! Variable declarations 
	DOUBLE PRECISION, DIMENSION(size(M,1)) :: xCol,yCol
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ytmp
	DOUBLE PRECISION, DIMENSION(2,2) :: Mint
	DOUBLE PRECISION, DIMENSION(1) :: ylow,yup,yInt
	DOUBLE PRECISION minxCol,maxxCol
	DOUBLE PRECISION xCol1,xCol2,yCol1,yCol2
	DOUBLE PRECISION xlow,xup
	INTEGER i,kif,ii
	
	interface
		subroutine InterpLin_M(fileM,x,colx,OutRange,colY,y,xlow,ylow,xup,yup)
			IMPLICIT NONE
			! Arguments declarations 
			DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: fileM
			DOUBLE PRECISION x,xlow,xup
			INTEGER, DIMENSION(:), INTENT(in) :: colY
			DOUBLE PRECISION, DIMENSION(size(colY)) :: ylow,yup,y
			INTEGER colx,OutRange
		end subroutine InterpLin_M
	end interface
	! 
	xCol=xx
	yCol=M(:,ny)
	
	minxCol=minval(xCol)
	maxxCol=maxval(xCol)
	
	if (x>=minxCol .and. x<=maxxCol) then
		allocate(ytmp(size(M,1))) !it remains void if x is outside the xCol range
		i=1
		kif=0
		do while (i<=size(xCol,1)-1) 
			xCol1=xCol(i)
			xCol2=xCol(i+1)
			if ((x>=xCol1 .and. x<=xCol2) .or. (x<=xCol1 .and. x>=xCol2)) then
				kif=kif+1
				yCol1=yCol(i)
				yCol2=yCol(i+1)
				Mint=reshape((/ xCol1,xCol2,yCol1,yCol2 /) , (/ 2 , 2/)) !ACCORDING TO COLUMN
				call InterpLin_M(Mint,x,1,1,(/2 /),yInt,xlow,ylow,xup,yup)
				ytmp(kif)=yInt(1)
			end if 
			i=i+1
		end do
		allocate(y(kif))
		do ii=1,kif
			y(ii)=ytmp(ii) !,1
		end do
	end if
	!If the if statement is not satisfied, y remains de-allocated

end subroutine findYgivenX_v



!Consider a matrix M. Search all the x values read at the column xx (if xx is a vector) or at the column of M
! whose index is xx (if xx is a scalar). 
!Give the corrisponding y values that are tabulated in the column whose index is ny 
!If M is the Leo track								|||	If M is PD isochrone 
!nx=	1	Mass			!nx=	12	LC			|||	!nx=	1	Z			!nx=	11	V 
!nx=	2	Age				!nx=	13	Lneutr		|||	!nx=	2	logt		!nx=	12	R 
!nx=	3	logL			!nx=	14	Lgrav		|||	!nx=	3	M_ini		!nx=	13	I 
!nx=	4	logTe			!nx=	15	Rstar		|||	!nx=	4	M_act		!nx=	14	J 
!nx=	5	Xcen			!nx=	16	Xsup		||| !nx=	5	logL		!nx=	15	H 
!nx=	6	Ycen			!nx=	17	Ysup		|||	!nx=	6	logTe		!nx=	16	K 
!nx=	7	XCcen			!nx=	18	XCsup		|||	!nx=	7	logg		!nx=	17	iIMF 
!nx=	8	XOcen			!nx=	19	XOsup		|||	!nx=	8	mbol		!nx=	18	stage 
!nx=	9	XNcen			!nx=	20	XNsup		|||	!nx=	9	U 
!nx=	10	Lx				!nx=	21	XC13sup		|||	!nx=	10	B 
!nx=	11	Ly 

!> 
subroutine findYgivenX_i(M,x,xx,ny,y)
	IMPLICIT NONE
	! Arguments declarations
	DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: M
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: y
	DOUBLE PRECISION x
	INTEGER xx,ny 
	! Variable declarations 
	DOUBLE PRECISION, DIMENSION(size(M,1)) :: xCol,yCol
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ytmp
	DOUBLE PRECISION, DIMENSION(2,2) :: Mint
	DOUBLE PRECISION, DIMENSION(1) :: ylow,yup,yInt
	DOUBLE PRECISION minxCol,maxxCol
	DOUBLE PRECISION xCol1,xCol2,yCol1,yCol2
	DOUBLE PRECISION xlow,xup
	INTEGER i,kif,ii
	
	interface
		subroutine InterpLin_M(fileM,x,colx,OutRange,colY,y,xlow,ylow,xup,yup)
			IMPLICIT NONE
			! Arguments declarations 
			DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: fileM
			DOUBLE PRECISION x,xlow,xup
			INTEGER, DIMENSION(:), INTENT(in) :: colY
			DOUBLE PRECISION, DIMENSION(size(colY)) :: ylow,yup,y
			INTEGER colx,OutRange
		end subroutine InterpLin_M
	end interface
	! 
	xCol=M(:,xx)
	yCol=M(:,ny)
	
	minxCol=minval(xCol)
	maxxCol=maxval(xCol)
	
	if (x>=minxCol .and. x<=maxxCol) then
		allocate(ytmp(size(M,1))) !it remains void if x is outside the xCol range
		i=1
		kif=0
		do while (i<=size(xCol,1)-1) 
			xCol1=xCol(i)
			xCol2=xCol(i+1)
			if ((x>=xCol1 .and. x<=xCol2) .or. (x<=xCol1 .and. x>=xCol2)) then
				kif=kif+1
				yCol1=yCol(i)
				yCol2=yCol(i+1)
				Mint=reshape((/ xCol1,xCol2,yCol1,yCol2 /) , (/ 2 , 2/)) !ACCORDING TO COLUMN
				call InterpLin_M(Mint,x,1,1,(/2 /),yInt,xlow,ylow,xup,yup)
				ytmp(kif)=yInt(1)
			end if 
			i=i+1
		end do
		allocate(y(kif))
		do ii=1,kif
			y(ii)=ytmp(ii) !,1
		end do
	end if
	!If the if statement is not satisfied, y remains de-allocated

end subroutine findYgivenX_i

!Given a table, where the values (xi, yi) are reported; 
!yi may be made of several columns, each representative of a parameter to be interpolated
!colY is the vector containg those column indices corresponding to the parameters to be interpolated
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!IT'S ASSUMED THAT THE xi COLUMN IS ORDERED IN ASCENDING ORDER!!
!! If t's not the case, matrix is ordered according to the xi ascending order
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!Considered a given x input value, the corresponding y is found according to linear interpolation
!between the two points of coordinates (x1,y1) and (x2,y2), having the closest abscissa to x

!If x is outside the range of xi values
! - OutRange=1: give the extreme value (first or last) of xi as the result of interpolation and
!	display a warning
! - OutRange=0: do not interpolate

!> 
subroutine InterpLin_M(fileM,x,colx,OutRange,colY,y,xlow,ylow,xup,yup)
	IMPLICIT NONE
	! Arguments declarations 
	DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: fileM
	DOUBLE PRECISION x,xlow,xup
	INTEGER, DIMENSION(:), INTENT(in) :: colY
	DOUBLE PRECISION, DIMENSION(size(colY)) :: ylow,yup,y
	INTEGER colx,OutRange
	
	! Variable declarations 
	DOUBLE PRECISION, DIMENSION(size(fileM,1)) :: xi
	DOUBLE PRECISION, DIMENSION(size(fileM,1),size(colY)) :: yi
	DOUBLE PRECISION, DIMENSION(size(colY)) :: y1,y2,m
	DOUBLE PRECISION x1,x2
	INTEGER, DIMENSION(size(fileM,1)) :: ndxsorted
	INTEGER, DIMENSION(1) :: ndx
	INTEGER i,ndx2
		
	interface
		subroutine indexx(arr,indx)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: arr !(n)
			INTEGER indx(size(arr))
		end subroutine indexx
	end interface

	xi=fileM(:,colx)
	yi=fileM(:,colY)
	
	call indexx(xi,ndxsorted) 
	xi=xi(ndxsorted) !6/10/17
	yi=yi(ndxsorted,:)
	xlow=xi(1)
	ylow=yi(1,:)
	xup=xi(size(xi,1))
	yup=yi(size(xi,1),:)

	if (x<xi(1)) then 
		if (OutRange==1) then 
			y=ylow
		else
			print*, 'from InterpLin_M: Out of interpolation range'
			write(*,*) x,' < ',xi(1)
			do i=1,size(y)
				call setNaN(y(i))
			end do 
		end if 
	else if (x>xi(size(xi))) then 
		if (OutRange==1) then 
			y=yup
		else
			print*, 'from InterpLin_M: Out of interpolation range'
			write(*,*) x,' > ',xi(size(xi))
			do i=1,size(y)
				call setNaN(y(i))
			end do
		end if 
	else
		!!! Modified on 6/09/2017 
		!!!	[diffsort, ndx]=sort(abs(x-xi)); 

		!!!	ndx1=ndx(1); 
		!!!	ndx2=ndx(2); 
		!!!	x1=xi(ndx1); 
		!!!	y1=yi(ndx1); 
		!!!	x2=xi(ndx2); 
		!!!	y2=yi(ndx2); 
		!!!
		ndx=minloc(abs(x-xi),1) 
		x1=xi(ndx(1))
		y1=yi(ndx(1),:)
		if (x1==xlow) then !ndx=1
			ndx2=ndx(1)+1
		else if (x1==xup) then !ndx=size(xi,1)
			ndx2=ndx(1)-1
		else
			if (x>x1) then 
				ndx2=ndx(1)+1
			else
				ndx2=ndx(1)-1
			end if 
		end if 
		x2=xi(ndx2)
		y2=yi(ndx2,:)

		m=(y2-y1)/(x2-x1) !Surely x1!=x2, otherwise there won't be univoc correspondence x-->y
		y=m*(x-x1)+y1
	end if 
end subroutine InterpLin_M

subroutine selectColExcept1(matrix,colx,colY)
	IMPLICIT NONE
	DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: matrix
	INTEGER, DIMENSION(size(matrix,2)-1) :: colY
	INTEGER colx
	
	INTEGER, DIMENSION(colx-1) :: colY01
	INTEGER, DIMENSION(size(matrix,2)-colx) :: colY02
	INTEGER colM,idum
	
	colM=size(matrix,2)
	if (colx==1) then
		colY=(/ (idum, idum=2,colM) /)
	else if (colx==colM) then
		colY=(/ (idum, idum=1,colM-1) /)
	else
		colY01=(/ (idum, idum=1,colx-1) /)
		colY02=(/ (idum, idum=colx+1,colM) /)
		colY=(/ colY01, colY02 /)
	end if

end subroutine selectColExcept1

!Open the track with given Z and M 
!PD refers to the theoretical models to be used 
!	PD=0 --> Padova models (PARSEC) v1.0 
!	PD=1 --> Padova models (PARSEC) v1.2S

!In output it gives:
!	- the track
!	- the mass of that track (not necessarily coincident with M because the track with M
!	as specified in input might not exist
!   - the Is string of the path, if it hasn't already been specified as an input parameter

!> 
subroutine openTracksPD(Z,M,PD,Tracks,M_in,Is)
	IMPLICIT NONE
	! Arguments declarations 
	DOUBLE PRECISION Z,M,M_in
	INTEGER PD
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Tracks
	CHARACTER(LEN=:), ALLOCATABLE :: Is
	! Variable declarations
	INTEGER xM,x,kdot,idum,i
	INTEGER kDiffZ,f_unit,ex,eof,k
	INTEGER, DIMENSION(:), ALLOCATABLE :: kDiffM
	
	CHARACTER :: TABUL=achar(9)
	CHARACTER(LEN=11) :: Tab
	CHARACTER(LEN=5) :: f7
	CHARACTER(LEN=4) :: outa
	CHARACTER(LEN=7) :: MinS
	CHARACTER(LEN=1000) txt,head
	CHARACTER(LEN=:), ALLOCATABLE :: TrPath,ML,Est,perc,Mstr
	CHARACTER(LEN=:), ALLOCATABLE :: ZinS,YinS
	
	DOUBLE PRECISION Mlow,Mlow1,Mlow2,Msup,Msup1,Msup2,DeltaM,DeltaM1,DeltaM2
	DOUBLE PRECISION Yin,Zin
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Zmod,Ymod
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: vdum,vdum1,vdum2,Mmod,diffM
	
	LOGICAL vuoto,isEq
	 
	interface
		subroutine loadMatrix(fileName,matrix,head)
			IMPLICIT NONE
			CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: fileName
			DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: matrix
			CHARACTER(LEN=1000) :: head
		end subroutine loadMatrix
		subroutine indexx(arr,indx)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: arr
			INTEGER indx(size(arr))
 		end subroutine indexx
 		subroutine num2str(x,d,str)
			IMPLICIT NONE
			DOUBLE PRECISION x
			CHARACTER(LEN=:), ALLOCATABLE :: str !was 40
			INTEGER d
		end subroutine num2str
	end interface

	if (allocated(Is)) then !Is is the input string specifying the path
		call loadMatrix(Is,Tracks,head)
		if (PD==0) then
			xM=scan(trim(Is),"M") !left-most 'M'. Should be the only one except for 'PMS' on its right
		else
			xM=scan(trim(Is),"M",.true.) !infer M value from the path (prima ndxM)
		end if

		x=xM+1
		kdot=0
		do while ((iachar(Is(x:x))>=48 .and. iachar(Is(x:x))<=57) .or. (kdot<1)) 
			if (iachar(Is(x:x))==46) then !found '.'
				kdot=kdot+1
				if (kdot>1) then 
					exit
				end if 
			end if 
			x=x+1
		end do 
		x=x-1

		Mstr=Is(xM+1:x)
		read(Mstr,'(F7.2)') M_in
		if (.not.isEq(M_in,M,2)) then
			print*, 'Whatch out: Input mass of openTracks is inconsistent with the input path'
		end if
		return
	end if 
	
	if (PD==0) then !!Padova v.1.0
		TrPath="/home/bonfanti/Documents/ArticoloTesi/TracceLeo/"
		Tab="tab2_S12D_Z"
		ML="OUTA1.74_F7_M"
		Est=".PMS.hrdat"
		
		allocate(Zmod(16)); allocate(Ymod(16))
		Zmod=(/ 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, &
			& 0.01, 0.014, 0.017, 0.02, 0.03, 0.04, 0.05, 0.06 /)
		Ymod=(/ 0.249, 0.249, 0.249, 0.25, 0.252, 0.256, 0.259, 0.263, 0.267, &
			& 0.273, 0.279, 0.284, 0.302, 0.321, 0.339, 0.356 /)

		kDiffZ=minloc(abs(Z-Zmod),1) 
		Zin=Zmod(kDiffZ)
		call num2str(Zin,4,ZinS)
		Yin=Ymod(kDiffZ)
		call num2str(Yin,3,YinS)
		
		allocate(Mmod(1200))
		Mmod=(/ (dble(idum)/100, idum=1,1200) /) !before 16/10/17 the massStep was 0.05
	else if (PD==1) then !!Padova PARSEC v1.2S
		TrPath="/home/bonfanti/Documents/PostDocLiegi/PARSEC1.2S/Tracce/"
		outa="OUTA" !1.74 o 1.77
		f7="_F7_M"
		Est=".DAT"
		
		allocate(Zmod(15)); allocate(Ymod(15))
		Zmod=(/ 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, &
			& 0.01, 0.014, 0.017, 0.02, 0.03, 0.04, 0.06 /)
		Ymod=(/ 0.249, 0.249, 0.249, 0.25, 0.252, 0.256, 0.259, 0.263, 0.267, &
			& 0.273, 0.279, 0.284, 0.302, 0.321, 0.356 /)

		kDiffZ=minloc(abs(Z-Zmod),1) 
		Zin=Zmod(kDiffZ)
		call num2str(Zin,4,ZinS)
		Yin=Ymod(kDiffZ)
		call num2str(Yin,3,YinS)
		
		Mlow=0.09
		Msup=12
		DeltaM=0.01
		allocate(vdum(nint((Msup-Mlow)/DeltaM)+1))
		vdum=(/ (dble(idum)/100, idum=nint(Mlow*100),nint(Msup*100)) /)
		! 
		DeltaM1=2
		Mlow1=Msup+DeltaM1
		Msup1=30
		allocate(vdum1(nint((Msup1-Mlow1)/DeltaM1)+1))
		vdum1=(/ (dble(idum), idum=nint(Mlow1),nint(Msup1),nint(DeltaM1)) /)
		! 
		DeltaM2=5
		Mlow2=Msup1+DeltaM2
		Msup2=350
		allocate(vdum2(nint((Msup2-Mlow2)/DeltaM2)+1))
		vdum2=(/ (dble(idum), idum=nint(Mlow2),nint(Msup2),nint(DeltaM2)) /)
		allocate(Mmod(size(vdum)+size(vdum1)+size(vdum2)))
		Mmod(1:size(vdum))=vdum
		Mmod(size(vdum)+1:size(vdum)+size(vdum1))=vdum1
		Mmod(size(vdum)+size(vdum1)+1:size(vdum)+size(vdum1)+size(vdum2))=vdum2
	end if 
	
	allocate(diffM(size(Mmod))); allocate(kDiffM(size(Mmod)))
	diffM=abs(M-Mmod)
	call indexx(diffM,kDiffM)
	k=1
	M_in=Mmod(kDiffM(k))
	if (PD==0) then 
		if (M_in<10) then 
			write(MinS,'(F4.2)') M_in !MinS=sprintf("%4.2f",Min);
			if (M_in<1) then 
				MinS=MinS(2:len(MinS))
			end if 
		else
			write(MinS,'(F4.1)') M_in
		end if 
		perc=trim(TrPath)//Tab//trim(ZinS)//"_Y"//trim(YinS)//"/Z"//trim(ZinS)// &
			& "Y"//trim(YinS)//trim(ML)//trim(MinS)//trim(Est)
	else if (PD==1) then 
		write(MinS,'(I3.3,F0.3)') int(M_in),M_in-int(M_in)
		if (M_in.le.0.70 .or. isEq(M_in,0.7D0,2)) then !!anint(M_in*100)/100<=0.70
			ML="1.77"
		else
			ML="1.74"
		end if 
		perc=trim(TrPath)//"Z"//trim(ZinS)//"Y"//trim(YinS)//"/Z"//trim(ZinS)// &
			& "Y"//trim(YinS)//outa//trim(ML)//f7//trim(MinS)//trim(Est)
	end if 

	call get_unit(f_unit)
	OPEN(UNIT=f_unit, FILE=perc, STATUS='OLD', IOSTAT=ex)
	
	vuoto=.true. !initialized even if not mandatory
	if (ex.eq.0) then !file does exist
		!!!!! 
		!!!If file doesn't contain data, it won't be charged!!! 
		READ(f_unit, '(A)', IOSTAT=eof) txt
		i=1
		do while(((txt(i:i).eq.' ').or.(txt(i:i).eq.TABUL)).and.(i.lt.len(trim(txt))))
			i=i+1
		end do
		do while ((txt(i:i)=="#".or.txt(i:i)=="%").and.eof==0) 
			READ(f_unit, '(A)', IOSTAT=eof) txt
		end do
		if (eof==1) then 
			vuoto=.true.
		else
			vuoto=.false.
		end if 
		CLOSE(f_unit)
		!!!End of part checking whether the file contains data!!! 
		!!!!!!! 
	end if
	do while ((ex.ne.0 .or. vuoto) .and. k<size(Mmod)) !file doesn't exist or it's empty
		k=k+1
		M_in=Mmod(kDiffM(k))
		if (PD==0) then 
			if (M_in<10) then 
				write(MinS,'(F4.2)') M_in
				if (M_in<1) then 
					MinS=MinS(2:len(MinS))
				end if 
			else
				write(MinS,'(F4.1)') M_in
			end if
			perc=trim(TrPath)//Tab//trim(ZinS)//"_Y"//trim(YinS)//"/Z"//trim(ZinS)// &
				& "Y"//trim(YinS)//trim(ML)//trim(MinS)//trim(Est) 
		else if (PD==1) then 
			write(MinS,'(I3.3,F0.3)') int(M_in),M_in-int(M_in)
			if (M_in.le.0.70 .or. isEq(M_in,0.7D0,2)) then 
				ML="1.77"
			else
				ML="1.74"
			end if 
			perc=trim(TrPath)//"Z"//trim(ZinS)//"Y"//trim(YinS)//"/Z"//trim(ZinS)// &
				& "Y"//trim(YinS)//outa//trim(ML)//f7//trim(MinS)//trim(Est)
		end if 
	
		call get_unit(f_unit)
		OPEN(UNIT=f_unit, FILE=perc, STATUS='OLD', IOSTAT=ex)
	
		if (ex.eq.0) then !file does exist
			!!!!! 
			!!!If file doesn't contain data, it won't be charged!!!
			READ(f_unit, '(A)', IOSTAT=eof) txt
			i=1
			do while(((txt(i:i).eq.' ').or.(txt(i:i).eq.TABUL)).and.(i.lt.len(trim(txt))))
				i=i+1
			end do
			do while ((txt(i:i)=="#".or.txt(i:i)=="%").and.eof==0) 
				READ(f_unit, '(A)', IOSTAT=eof) txt
			end do
			if (eof==1) then 
				vuoto=.true.
			else
				vuoto=.false.
			end if 
			CLOSE(f_unit)
			!!!End of part checking whether the file contains data!!! 
			!!!!!!! 
		end if 
	end do 
	if (ex.ne.0 .or. vuoto) then 
		print*, "from openTracksPD: File does not exist"
		return
	end if 
	
	call loadMatrix(perc,Tracks,head)
	Is=perc

end subroutine openTracksPD

subroutine checkTracksPD(Z,M,PD,M_in0,M_in,sameFile,Is)
!Given the input M, find out whether the tracks to be opened (M_in)
!is the same as M_in0. If it's the same (sameFile=.true.), this file won't be
!loaded in the main
	IMPLICIT NONE
	! Arguments declarations 
	DOUBLE PRECISION Z,M,M_in0,M_in
	INTEGER PD
	LOGICAL sameFile
	CHARACTER(LEN=:), ALLOCATABLE :: Is
	! Variable declarations
	INTEGER idum,i
	INTEGER kDiffZ,f_unit,ex,eof,k
	INTEGER, DIMENSION(:), ALLOCATABLE :: kDiffM
	
	CHARACTER :: TABUL=achar(9)
	CHARACTER(LEN=11) :: Tab
	CHARACTER(LEN=5) :: f7
	CHARACTER(LEN=4) :: outa
	CHARACTER(LEN=7) :: MinS
	CHARACTER(LEN=1000) txt
	CHARACTER(LEN=:), ALLOCATABLE :: TrPath,ML,Est,perc
	CHARACTER(LEN=:), ALLOCATABLE :: ZinS,YinS
	
	DOUBLE PRECISION Mlow,Mlow1,Mlow2,Msup,Msup1,Msup2,DeltaM,DeltaM1,DeltaM2
	DOUBLE PRECISION Yin,Zin
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Zmod,Ymod
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: vdum,vdum1,vdum2,Mmod,diffM
	
	LOGICAL vuoto,isEq
	 
	interface
		subroutine indexx(arr,indx)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: arr
			INTEGER indx(size(arr))
 		end subroutine indexx
 		subroutine num2str(x,d,str)
			IMPLICIT NONE
			DOUBLE PRECISION x
			CHARACTER(LEN=:), ALLOCATABLE :: str !was 40
			INTEGER d
		end subroutine num2str
	end interface

	if (PD==0) then !!Padova v.1.0
		TrPath="/home/bonfanti/Documents/ArticoloTesi/TracceLeo/"
		Tab="tab2_S12D_Z"
		ML="OUTA1.74_F7_M"
		Est=".PMS.hrdat"
		
		allocate(Zmod(16)); allocate(Ymod(16))
		Zmod=(/ 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, &
			& 0.01, 0.014, 0.017, 0.02, 0.03, 0.04, 0.05, 0.06 /)
		Ymod=(/ 0.249, 0.249, 0.249, 0.25, 0.252, 0.256, 0.259, 0.263, 0.267, &
			& 0.273, 0.279, 0.284, 0.302, 0.321, 0.339, 0.356 /)

		kDiffZ=minloc(abs(Z-Zmod),1) 
		Zin=Zmod(kDiffZ)
		call num2str(Zin,4,ZinS)
		Yin=Ymod(kDiffZ)
		call num2str(Yin,3,YinS)
		
		allocate(Mmod(1200))
		Mmod=(/ (dble(idum)/100, idum=1,1200) /) !before 16/10/17 the massStep was 0.05
	else if (PD==1) then !!Padova PARSEC v1.2S
		TrPath="/home/bonfanti/Documents/PostDocLiegi/PARSEC1.2S/Tracce/"
		outa="OUTA" !1.74 o 1.77
		f7="_F7_M"
		Est=".DAT"
		
		allocate(Zmod(15)); allocate(Ymod(15))
		Zmod=(/ 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, &
			& 0.01, 0.014, 0.017, 0.02, 0.03, 0.04, 0.06 /)
		Ymod=(/ 0.249, 0.249, 0.249, 0.25, 0.252, 0.256, 0.259, 0.263, 0.267, &
			& 0.273, 0.279, 0.284, 0.302, 0.321, 0.356 /)

		kDiffZ=minloc(abs(Z-Zmod),1) 
		Zin=Zmod(kDiffZ)
		call num2str(Zin,4,ZinS)
		Yin=Ymod(kDiffZ)
		call num2str(Yin,3,YinS)
		
		Mlow=0.09
		Msup=12
		DeltaM=0.01
		allocate(vdum(nint((Msup-Mlow)/DeltaM)+1))
		vdum=(/ (dble(idum)/100, idum=nint(Mlow*100),nint(Msup*100)) /)
		! 
		DeltaM1=2
		Mlow1=Msup+DeltaM1
		Msup1=30
		allocate(vdum1(nint((Msup1-Mlow1)/DeltaM1)+1))
		vdum1=(/ (dble(idum), idum=nint(Mlow1),nint(Msup1),nint(DeltaM1)) /)
		! 
		DeltaM2=5
		Mlow2=Msup1+DeltaM2
		Msup2=350
		allocate(vdum2(nint((Msup2-Mlow2)/DeltaM2)+1))
		vdum2=(/ (dble(idum), idum=nint(Mlow2),nint(Msup2),nint(DeltaM2)) /)
		allocate(Mmod(size(vdum)+size(vdum1)+size(vdum2)))
		Mmod(1:size(vdum))=vdum
		Mmod(size(vdum)+1:size(vdum)+size(vdum1))=vdum1
		Mmod(size(vdum)+size(vdum1)+1:size(vdum)+size(vdum1)+size(vdum2))=vdum2
	end if 
	
	allocate(diffM(size(Mmod))); allocate(kDiffM(size(Mmod)))
	diffM=abs(M-Mmod)
	call indexx(diffM,kDiffM)
	k=1
	M_in=Mmod(kDiffM(k))
	if (PD==0) then 
		if (M_in<10) then 
			write(MinS,'(F4.2)') M_in !MinS=sprintf("%4.2f",Min);
			if (M_in<1) then 
				MinS=MinS(2:len(MinS))
			end if 
		else
			write(MinS,'(F4.1)') M_in
		end if 
		perc=trim(TrPath)//Tab//trim(ZinS)//"_Y"//trim(YinS)//"/Z"//trim(ZinS)// &
			& "Y"//trim(YinS)//trim(ML)//trim(MinS)//trim(Est)
	else if (PD==1) then 
		write(MinS,'(I3.3,F0.3)') int(M_in),M_in-int(M_in)
		if (M_in.le.0.70 .or. isEq(M_in,0.7D0,2)) then 
			ML="1.77"
		else
			ML="1.74"
		end if 
		perc=trim(TrPath)//"Z"//trim(ZinS)//"Y"//trim(YinS)//"/Z"//trim(ZinS)// &
			& "Y"//trim(YinS)//outa//trim(ML)//f7//trim(MinS)//trim(Est)
	end if 

	call get_unit(f_unit)
	OPEN(UNIT=f_unit, FILE=perc, STATUS='OLD', IOSTAT=ex)
	
	vuoto=.true. !initialized even if not mandatory
	if (ex.eq.0) then !file does exist
		!!!!! 
		!!!If file doesn't contain data, it won't be charged!!! 
		READ(f_unit, '(A)', IOSTAT=eof) txt
		i=1
		do while(((txt(i:i).eq.' ').or.(txt(i:i).eq.TABUL)).and.(i.lt.len(trim(txt))))
			i=i+1
		end do
		do while ((txt(i:i)=="#".or.txt(i:i)=="%").and.eof==0) 
			READ(f_unit, '(A)', IOSTAT=eof) txt
		end do
		if (eof==1) then 
			vuoto=.true.
		else
			vuoto=.false.
		end if 
		CLOSE(f_unit)
		!!!End of part checking whether the file contains data!!! 
		!!!!!!! 
	end if
	do while ((ex.ne.0 .or. vuoto) .and. k<size(Mmod)) !file doesn't exist or it's empty
		k=k+1
		M_in=Mmod(kDiffM(k))
		if (PD==0) then 
			if (M_in<10) then 
				write(MinS,'(F4.2)') M_in
				if (M_in<1) then 
					MinS=MinS(2:len(MinS))
				end if 
			else
				write(MinS,'(F4.1)') M_in
			end if
			perc=trim(TrPath)//Tab//trim(ZinS)//"_Y"//trim(YinS)//"/Z"//trim(ZinS)// &
				& "Y"//trim(YinS)//trim(ML)//trim(MinS)//trim(Est) 
		else if (PD==1) then 
			write(MinS,'(I3.3,F0.3)') int(M_in),M_in-int(M_in)
			if (M_in.le.0.70 .or. isEq(M_in,0.7D0,2)) then 
				ML="1.77"
			else
				ML="1.74"
			end if 
			perc=trim(TrPath)//"Z"//trim(ZinS)//"Y"//trim(YinS)//"/Z"//trim(ZinS)// &
				& "Y"//trim(YinS)//outa//trim(ML)//f7//trim(MinS)//trim(Est)
		end if 
	
		call get_unit(f_unit)
		OPEN(UNIT=f_unit, FILE=perc, STATUS='OLD', IOSTAT=ex)
	
		if (ex.eq.0) then !file does exist
			!!!!! 
			!!!If file doesn't contain data, it won't be charged!!!
			READ(f_unit, '(A)', IOSTAT=eof) txt
			i=1
			do while(((txt(i:i).eq.' ').or.(txt(i:i).eq.TABUL)).and.(i.lt.len(trim(txt))))
				i=i+1
			end do
			do while ((txt(i:i)=="#".or.txt(i:i)=="%").and.eof==0) 
				READ(f_unit, '(A)', IOSTAT=eof) txt
			end do
			if (eof==1) then 
				vuoto=.true.
			else
				vuoto=.false.
			end if 
			CLOSE(f_unit)
			!!!End of part checking whether the file contains data!!! 
			!!!!!!! 
		end if 
	end do 
	if (ex.ne.0 .or. vuoto) then 
		print*, "from checkTracksPD: File does not exist"
		sameFile=.true. !to prevent loading in the main, where the idea is to load
						!the file if it's different from M_in0
	else
		if (isEq(M_in,M_in0,3)) then
			sameFile=.true.
		else
			sameFile=.false.
			Is=perc !if sameFile=.true. it remains deallocated
		end if
	end if
	
end subroutine checkTracksPD

subroutine findTracksPD(Z,M,PD,M_in,Zin,found)
!Given the mass M, find out which is the nearest mass M_in among the available tracks
	IMPLICIT NONE
	! Arguments declarations 
	DOUBLE PRECISION Z,M,M_in,Zin
	INTEGER PD
	! Variable declarations
	INTEGER idum,i
	INTEGER kDiffZ,f_unit,ex,eof,k
	INTEGER, DIMENSION(:), ALLOCATABLE :: kDiffM
	
	CHARACTER :: TABUL=achar(9)
	CHARACTER(LEN=11) :: Tab
	CHARACTER(LEN=5) :: f7
	CHARACTER(LEN=4) :: outa
	CHARACTER(LEN=7) :: MinS
	CHARACTER(LEN=1000) txt
	CHARACTER(LEN=:), ALLOCATABLE :: TrPath,ML,Est,perc
	CHARACTER(LEN=:), ALLOCATABLE :: ZinS,YinS
	
	DOUBLE PRECISION Mlow,Mlow1,Mlow2,Msup,Msup1,Msup2,DeltaM,DeltaM1,DeltaM2
	DOUBLE PRECISION Yin
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Zmod,Ymod
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: vdum,vdum1,vdum2,Mmod,diffM
	
	LOGICAL vuoto,found,isEq
	 
	interface
		subroutine indexx(arr,indx)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: arr
			INTEGER indx(size(arr))
 		end subroutine indexx
 		subroutine num2str(x,d,str)
			IMPLICIT NONE
			DOUBLE PRECISION x
			CHARACTER(LEN=:), ALLOCATABLE :: str !was 40
			INTEGER d
		end subroutine num2str
	end interface

	if (PD==0) then !!Padova v.1.0
		TrPath="/home/bonfanti/Documents/ArticoloTesi/TracceLeo/"
		Tab="tab2_S12D_Z"
		ML="OUTA1.74_F7_M"
		Est=".PMS.hrdat"
		
		allocate(Zmod(16)); allocate(Ymod(16))
		Zmod=(/ 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, &
			& 0.01, 0.014, 0.017, 0.02, 0.03, 0.04, 0.05, 0.06 /)
		Ymod=(/ 0.249, 0.249, 0.249, 0.25, 0.252, 0.256, 0.259, 0.263, 0.267, &
			& 0.273, 0.279, 0.284, 0.302, 0.321, 0.339, 0.356 /)

		kDiffZ=minloc(abs(Z-Zmod),1) 
		Zin=Zmod(kDiffZ)
		call num2str(Zin,4,ZinS)
		Yin=Ymod(kDiffZ)
		call num2str(Yin,3,YinS)
		
		allocate(Mmod(1200))
		Mmod=(/ (dble(idum)/100, idum=1,1200) /) !before 16/10/17 the massStep was 0.05
	else if (PD==1) then !!Padova PARSEC v1.2S
		TrPath="/home/bonfanti/Documents/PostDocLiegi/PARSEC1.2S/Tracce/"
		outa="OUTA" !1.74 o 1.77
		f7="_F7_M"
		Est=".DAT"
		
		allocate(Zmod(15)); allocate(Ymod(15))
		Zmod=(/ 0.0001, 0.0002, 0.0005, 0.001, 0.002, 0.004, 0.006, 0.008, &
			& 0.01, 0.014, 0.017, 0.02, 0.03, 0.04, 0.06 /)
		Ymod=(/ 0.249, 0.249, 0.249, 0.25, 0.252, 0.256, 0.259, 0.263, 0.267, &
			& 0.273, 0.279, 0.284, 0.302, 0.321, 0.356 /)

		kDiffZ=minloc(abs(Z-Zmod),1) 
		Zin=Zmod(kDiffZ)
		call num2str(Zin,4,ZinS)
		Yin=Ymod(kDiffZ)
		call num2str(Yin,3,YinS)
		
		Mlow=0.09
		Msup=12
		DeltaM=0.01
		allocate(vdum(nint((Msup-Mlow)/DeltaM)+1))
		vdum=(/ (dble(idum)/100, idum=nint(Mlow*100),nint(Msup*100)) /)
		! 
		DeltaM1=2
		Mlow1=Msup+DeltaM1
		Msup1=30
		allocate(vdum1(nint((Msup1-Mlow1)/DeltaM1)+1))
		vdum1=(/ (dble(idum), idum=nint(Mlow1),nint(Msup1),nint(DeltaM1)) /)
		! 
		DeltaM2=5
		Mlow2=Msup1+DeltaM2
		Msup2=350
		allocate(vdum2(nint((Msup2-Mlow2)/DeltaM2)+1))
		vdum2=(/ (dble(idum), idum=nint(Mlow2),nint(Msup2),nint(DeltaM2)) /)
		allocate(Mmod(size(vdum)+size(vdum1)+size(vdum2)))
		Mmod(1:size(vdum))=vdum
		Mmod(size(vdum)+1:size(vdum)+size(vdum1))=vdum1
		Mmod(size(vdum)+size(vdum1)+1:size(vdum)+size(vdum1)+size(vdum2))=vdum2
	end if 
	
	allocate(diffM(size(Mmod))); allocate(kDiffM(size(Mmod)))
	diffM=abs(M-Mmod)
	call indexx(diffM,kDiffM)
	k=1
	M_in=Mmod(kDiffM(k))
	if (PD==0) then 
		if (M_in<10) then 
			write(MinS,'(F4.2)') M_in !MinS=sprintf("%4.2f",Min);
			if (M_in<1) then 
				MinS=MinS(2:len(MinS))
			end if 
		else
			write(MinS,'(F4.1)') M_in
		end if 
		perc=trim(TrPath)//Tab//trim(ZinS)//"_Y"//trim(YinS)//"/Z"//trim(ZinS)// &
			& "Y"//trim(YinS)//trim(ML)//trim(MinS)//trim(Est)
	else if (PD==1) then 
		write(MinS,'(I3.3,F0.3)') int(M_in),M_in-int(M_in)
		if (M_in.le.0.70 .or. isEq(M_in,0.7D0,2)) then 
			ML="1.77"
		else
			ML="1.74"
		end if 
		perc=trim(TrPath)//"Z"//trim(ZinS)//"Y"//trim(YinS)//"/Z"//trim(ZinS)// &
			& "Y"//trim(YinS)//outa//trim(ML)//f7//trim(MinS)//trim(Est)
	end if 

	call get_unit(f_unit)
	OPEN(UNIT=f_unit, FILE=perc, STATUS='OLD', IOSTAT=ex)
	
	vuoto=.true. !initialized even if not mandatory
	if (ex.eq.0) then !file does exist
		!!!!! 
		!!!If file doesn't contain data, it won't be charged!!! 
		READ(f_unit, '(A)', IOSTAT=eof) txt
		i=1
		do while(((txt(i:i).eq.' ').or.(txt(i:i).eq.TABUL)).and.(i.lt.len(trim(txt))))
			i=i+1
		end do
		do while ((txt(i:i)=="#".or.txt(i:i)=="%").and.eof==0) 
			READ(f_unit, '(A)', IOSTAT=eof) txt
		end do
		if (eof==1) then 
			vuoto=.true.
		else
			vuoto=.false.
		end if 
		CLOSE(f_unit)
		!!!End of part checking whether the file contains data!!! 
		!!!!!!! 
	end if
	do while ((ex.ne.0 .or. vuoto) .and. k<size(Mmod)) !file doesn't exist or it's empty
		k=k+1
		M_in=Mmod(kDiffM(k))
		if (PD==0) then 
			if (M_in<10) then 
				write(MinS,'(F4.2)') M_in
				if (M_in<1) then 
					MinS=MinS(2:len(MinS))
				end if 
			else
				write(MinS,'(F4.1)') M_in
			end if
			perc=trim(TrPath)//Tab//trim(ZinS)//"_Y"//trim(YinS)//"/Z"//trim(ZinS)// &
				& "Y"//trim(YinS)//trim(ML)//trim(MinS)//trim(Est) 
		else if (PD==1) then 
			write(MinS,'(I3.3,F0.3)') int(M_in),M_in-int(M_in)
			if (M_in.le.0.70 .or. isEq(M_in,0.7D0,2)) then 
				ML="1.77"
			else
				ML="1.74"
			end if 
			perc=trim(TrPath)//"Z"//trim(ZinS)//"Y"//trim(YinS)//"/Z"//trim(ZinS)// &
				& "Y"//trim(YinS)//outa//trim(ML)//f7//trim(MinS)//trim(Est)
		end if 
	
		call get_unit(f_unit)
		OPEN(UNIT=f_unit, FILE=perc, STATUS='OLD', IOSTAT=ex)
	
		if (ex.eq.0) then !file does exist
			!!!!! 
			!!!If file doesn't contain data, it won't be charged!!!
			READ(f_unit, '(A)', IOSTAT=eof) txt
			i=1
			do while(((txt(i:i).eq.' ').or.(txt(i:i).eq.TABUL)).and.(i.lt.len(trim(txt))))
				i=i+1
			end do
			do while ((txt(i:i)=="#".or.txt(i:i)=="%").and.eof==0) 
				READ(f_unit, '(A)', IOSTAT=eof) txt
			end do
			if (eof==1) then 
				vuoto=.true.
			else
				vuoto=.false.
			end if 
			CLOSE(f_unit)
			!!!End of part checking whether the file contains data!!! 
			!!!!!!! 
		end if 
	end do 
	if (ex.ne.0 .or. vuoto) then
		found=.false.
		print*, "from findTracksPD: File does not exist"
	else
		found=.true.
	end if
	
end subroutine findTracksPD

subroutine selectMass(vrho,vg,I_vrho,I_vg,Teff,I_Teff,FeH,I_FeH,Rinp,I_Rinp,Mabs,idCol,Minf,Msup)
	USE empRel, only : blogg,brho
	USE SunPD
	IMPLICIT NONE
	INTEGER idCol
	DOUBLE PRECISION, DIMENSION(3) :: vrho,I_vrho
	DOUBLE PRECISION, DIMENSION(2) :: vg,I_vg
	DOUBLE PRECISION Teff,I_Teff,FeH,I_FeH,Rinp,I_Rinp,Mabs,Minf,Msup
	
	INTEGER i,j,sz,k
	INTEGER, DIMENSION(:), ALLOCATABLE :: ir,ig
	DOUBLE PRECISION rho,rho1,rho2,I_rho,logrho
	DOUBLE PRECISION g,g1,g2,I_g,logg,M,I_M,XT,Re
	DOUBLE PRECISION, DIMENSION(10) :: vecT
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: MlVec,MuVec
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: MlMat,MuMat
	LOGICAL isEq
	
	interface
		subroutine find(b,ix)
			IMPLICIT NONE
			LOGICAL, DIMENSION(:), INTENT(in) :: b
			INTEGER, DIMENSION(:), ALLOCATABLE :: ix
		end subroutine find
		
		function isEq_v(x1,x2,n)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: x1
			DOUBLE PRECISION, INTENT(in) :: x2
			INTEGER, INTENT(in) :: n
			LOGICAL, DIMENSION(size(x1)) :: isEq_v
		end function isEq_v
	end interface
	
	!vrho=(/ rho_t, deltaNu, rho_i /)
	!vg=(/ logg, numax /)
	call find(.not.isEq_v(vrho,-1.D0,1),ir)
	call find(.not.isEq_v(vg,-1.D0,1),ig)
	if (.not.isEq(Rinp,-1.D0,2)) then
		sz=3
	else
		sz=0
	end if
	k=0
	if (allocated(ir).and.allocated(ig)) then
		allocate(MlMat(size(ir),size(ig)+sz))
		allocate(MuMat(size(ir),size(ig)+sz))
		do i=1,size(ir)
			if (ir(i).eq.2) then
				rho=(vrho(ir(i))/DnuSun)**2 !solar
				rho1=((vrho(ir(i))+I_vrho(ir(i)))/DnuSun)**2
				rho2=((vrho(ir(i))-I_vrho(ir(i)))/DnuSun)**2
				I_rho=(abs(rho1-rho)+abs(rho2-rho))/2.
			else
				rho=vrho(ir(i)) !solar
				I_rho=I_vrho(ir(i))
			end if
			if (sz.eq.3) then
				M=rho*Rinp**3
				I_M=M*(I_rho/rho+3.*I_Rinp/Rinp)
				k=k+1
				MlMat(i,k)=M-I_M
				MuMat(i,k)=M+I_M
			end if
			do j=1,size(ig)
				if (ig(j).eq.2) then
					g=vg(ig(j))/numaxSun*sqrt(Teff/TeffSun) !solar
					g1=(vg(ig(j))+I_vg(ig(j)))/numaxSun*sqrt(Teff/TeffSun)
					g2=(vg(ig(j))-I_vg(ig(j)))/numaxSun*sqrt(Teff/TeffSun)
					I_g=(abs(g1-g)+abs(g2-g))/2.
				else
					g=10**(vg(ig(j))-loggSun) !solar
					I_g=g*I_vg(ig(j))*log(10.)
				end if
				if (sz.eq.3) then
					M=g*Rinp**2
					I_M=M*(I_g/g+2.*I_Rinp/Rinp)
					k=k+1
					MlMat(i,k)=M-I_M
					MuMat(i,k)=M+I_M
				end if
				M=g**3/rho**2
				I_M=M*(3.*I_g/g+2.*I_rho/rho)
				k=k+1
				MlMat(i,k)=M-I_M
				MuMat(i,k)=M+I_M
			end do
		end do
		Minf=minval(MlMat)
		Msup=maxval(MuMat)
	else
		if (allocated(ir).or.allocated(ig)) then
			XT=log10(Teff)-4.1
			if (allocated(ir)) then
				allocate(MlVec(size(ir)))
				allocate(MuVec(size(ir)))
				do i=1,size(ir)
					if (ir(i).eq.2) then
						rho=(vrho(ir(i))/DnuSun)**2 !solar
						rho1=((vrho(ir(i))+I_vrho(ir(i)))/DnuSun)**2
						rho2=((vrho(ir(i))-I_vrho(ir(i)))/DnuSun)**2
						I_rho=(abs(rho1-rho)+abs(rho2-rho))/2.
					else
						rho=vrho(ir(i)) !solar
						I_rho=I_vrho(ir(i))
					end if
					if (sz.eq.3) then
						M=rho*Rinp**3
						I_M=M*(I_rho/rho+3.*I_Rinp/Rinp)
					else
						logrho=log10(rho*rhoSun) !g/cm3
						vecT=(/ 1.D0,XT,XT**2,XT**3,logrho,logrho**2,logrho**3,FeH,FeH**2,FeH**3 /)
						Re=10.**(dot_product(brho,vecT)) !solar. Relative uncertainty Torres09: 3.2%
						M=rho*Re**3
						I_M=M*(I_rho/rho+3.*sqrt(0.032**2+(I_FeH*log(10.))**2+(I_Teff/Teff)**2))
						I_M=min(I_M,max(0.5,0.5*M))
					end if
					MlVec(i)=M-I_M
					MuVec(i)=M+I_M
				end do
			else !allocated(ig)
				allocate(MlVec(size(ig)))
				allocate(MuVec(size(ig)))
				do j=1,size(ig)
					if (ig(j).eq.2) then
						g=vg(ig(j))/numaxSun*sqrt(Teff/TeffSun) !solar
						g1=(vg(ig(j))+I_vg(ig(j)))/numaxSun*sqrt(Teff/TeffSun)
						g2=(vg(ig(j))-I_vg(ig(j)))/numaxSun*sqrt(Teff/TeffSun)
						I_g=(abs(g1-g)+abs(g2-g))/2.
					else
						g=10**(vg(ig(j))-loggSun) !solar
						I_g=g*I_vg(ig(j))*log(10.)
					end if
					if (sz.eq.3) then
						M=g*Rinp**2
						I_M=M*(I_g/g+2.*I_Rinp/Rinp)
					else
						logg=log10(g)+loggSun !cgs
						vecT=(/ 1.D0,XT,XT**2,XT**3,logg,logg**2,logg**3,FeH,FeH**2,FeH**3 /)
						Re=10.**(dot_product(blogg,vecT)) !solar (unc. 3.2%)
						M=g*Re**2
						I_M=M*(I_g/g+2.*sqrt(0.032**2+(I_FeH*log(10.))**2+(I_Teff/Teff)**2))
						I_M=min(I_M,max(0.5,0.5*M))
					end if
					MlVec(j)=M-I_M
					MuVec(j)=M+I_M
				end do
			end if
			Minf=minval(MlVec)
			Msup=maxval(MuVec)
		else
			call MfromEker(Teff,Rinp,Mabs,idCol,M,I_M)
			Minf=M-I_M
			Msup=M+I_M
		end if
	end if
	Minf=max(Minf,0.05)
	
end subroutine selectMass

subroutine MfromEker(Teff,R,Mabs,idCol,M,I_M)
	USE empRel, only : logLsup,aE,bE
	USE SunPD

	IMPLICIT NONE
	INTEGER idCol
	DOUBLE PRECISION Teff,R,Mabs,M,I_M
	
	DOUBLE PRECISION logL,logM
	LOGICAL isEq
	
	if (.not.isEq(R,-1.D0,2)) then
		logL=log10(R**2*(Teff/TeffSun)**4)
	else !necessarily Mabs available (.neq.-100)
		logL=-0.4*(Mabs-MabsSun(idCol))
	end if
	if (logL.lt.logLsup(1)) then
		logM=aE(1)*logL+bE(1)
	else if (logL.lt.logLsup(2)) then
		logM=aE(2)*logL+bE(2)
	else if (logL.lt.logLsup(3)) then
		logM=aE(3)*logL+bE(3)
	else
		logM=aE(4)*logL+bE(4)
	end if
	M=10**logM
	I_M=1.5
	I_M=min(I_M,max(0.5,0.5*M))

end subroutine MfromEker

subroutine tracksNumber(Minf,Msup,Mstep,MlowMS,Zvec,fnd05,Zt,nM)
	IMPLICIT NONE
	
	DOUBLE PRECISION Minf,Msup,Mstep,MlowMS
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: Zvec
	DOUBLE PRECISION, DIMENSION(size(Zvec)) :: Zt
	INTEGER, DIMENSION(size(Zvec)) :: nM
	LOGICAL, DIMENSION(size(Zvec)) :: fnd05
	!
	INTEGER iz
	LOGICAL fnd,found,isEq,sameFile
	DOUBLE PRECISION Ml05,MTr,M_in0,Z_in,M_in
	CHARACTER(LEN=:), ALLOCATABLE :: Isp
	
	interface
		subroutine checkTracksPD(Z,M,PD,M_in0,M_in,sameFile,Is)
		  IMPLICIT NONE
		  DOUBLE PRECISION Z,M,M_in0,M_in
		  INTEGER PD
		  LOGICAL sameFile
		  CHARACTER(LEN=:), ALLOCATABLE :: Is
		end subroutine checkTracksPD
	end interface
	
	fnd=.false.
	!for each metallicity, evaluate number of tracks (nM) to be loaded
	do iz=1,size(Zvec)
		fnd05(iz)=.false.
		call findTracksPD(Zvec(iz),MlowMS,1,Ml05,Zt(iz),found)
		if (abs(MlowMS-Ml05)/Ml05.lt.0.15) then
			fnd=.true.
		end if
		MTr=Minf
		call findTracksPD(Zvec(iz),MTr,1,M_in0,Z_in,found)
		if (found) then
			nM(iz)=1
			if (fnd .and. isEq(M_in0,Ml05,2)) then
				fnd05(iz)=.true.
			end if
		else !no tracks exist in that metallicity folder. That metallicity value should be
			 !removed from the list in check/find/openTracksPD
			nM(iz)=0
			print*,'No tracks in the specified metallicity folder'; cycle
		end if
		do while (MTr.lt.Msup) !+Mstep REMOVED !.le.
			MTr=MTr+Mstep
			call checkTracksPD(Zvec(iz),MTr,1,M_in0,M_in,sameFile,Isp)
			if (.not.sameFile) then
				deallocate(Isp) !if sameFile=.true., Isp is not allocated
				nM(iz)=nM(iz)+1
				M_in0=M_in
				if (fnd .and. isEq(M_in0,Ml05,2)) then
					fnd05(iz)=.true.
				end if
			end if
		end do
		if (.not.fnd05(iz)) then !if not already in, add also MlowMS for gyro
			nM(iz)=nM(iz)+1
		end if
	end do
	
end subroutine tracksNumber

subroutine storeIsoc(Zvec,IsocTab,Zndxi,Zndxf,idCol)
	IMPLICIT NONE
	
	INTEGER idCol
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: Zvec
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: IsocTab
	INTEGER, DIMENSION(:), ALLOCATABLE :: Zndxi,Zndxf
	
	INTEGER iz,row0
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Isoc0
	CHARACTER(LEN=:), ALLOCATABLE :: kz,path
	CHARACTER(LEN=:), ALLOCATABLE :: percorso
	CHARACTER(LEN=1000) head
	
	interface
		subroutine num2str(x,d,str)
			IMPLICIT NONE
			DOUBLE PRECISION x
			CHARACTER(LEN=:), ALLOCATABLE :: str !was 40
			INTEGER d
		end subroutine num2str
		subroutine loadMatrix(fileName,matrix,head)
			IMPLICIT NONE
			CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: fileName
			DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: matrix
			CHARACTER(LEN=1000) :: head
		end subroutine loadMatrix
		subroutine append(A,B)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: A
			DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: B
		end subroutine append
		subroutine avoidLastFlag(Isoc)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: Isoc
		end subroutine avoidLastFlag
	end interface
	
	print*,'begin storeIsoc'
	select case (idCol)
	case (1)
		percorso="/home/bonfanti/Documents/PostDocLiegi/PARSEC1.2S/IsocroneRed/Z" !reduced tables
	case (2)
		percorso="/home/bonfanti/Documents/PostDocLiegi/PARSEC1.2S/IsocroneGaia/Z" !reduced tables
	case default
		percorso="/home/bonfanti/Documents/PostDocLiegi/PARSEC1.2S/IsocroneRed/Z" !reduced tables
	end select
	
	allocate(Zndxi(size(Zvec))); allocate(Zndxf(size(Zvec)))
	row0=0
	Zndxi(1)=1
	do iz=1,size(Zvec)
		call num2str(Zvec(iz),5,kz) !5: max number of decimal digits
						
		path=trim(percorso)//trim(kz)//'.dat'
		deallocate(kz)
		call loadMatrix(path,Isoc0,head)
		call avoidLastFlag(Isoc0)
		if (iz.eq.1) then
			allocate(IsocTab(size(Isoc0,1),size(Isoc0,2)))
			IsocTab=Isoc0
		else
			call append(IsocTab,Isoc0)
		end if
		Zndxf(iz)=row0+size(Isoc0,1)
		row0=Zndxf(iz)
		if (iz.gt.1) then
			Zndxi(iz)=Zndxf(iz-1)+1
		end if
		deallocate(path); deallocate(Isoc0)
	end do
	
	print*,'end storeIsoc'

end subroutine storeIsoc

subroutine storeEvoZ(Minf,Msup,ZtevoTab,Endxi,Endxf,MeAv)

	USE IsoPD, only : MLeo
	
	IMPLICIT NONE
	DOUBLE PRECISION Minf,Msup
	DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: MeAv
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ZtevoTab
	INTEGER, DIMENSION(:), ALLOCATABLE :: Endxi,Endxf
	
	INTEGER i,row0
	INTEGER, DIMENSION(:), ALLOCATABLE :: ix
	DOUBLE PRECISION Me
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Ztevo0
	CHARACTER(LEN=7) Mgstr
	CHARACTER(LEN=62) :: EvoZ="/home/bonfanti/Documents/PostDocLiegi/PARSEC1.2S/Tracce/EvoZ/M"
	CHARACTER(LEN=1000) head
	CHARACTER(LEN=:), ALLOCATABLE :: path
	LOGICAL, DIMENSION(size(MLeo)) :: Mb
	
	interface
		subroutine find(b,ix)
			IMPLICIT NONE
			LOGICAL, DIMENSION(:), INTENT(in) :: b
			INTEGER, DIMENSION(:), ALLOCATABLE :: ix
		end subroutine find
		subroutine loadMatrix(fileName,matrix,head)
			IMPLICIT NONE
			CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: fileName
			DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: matrix
			CHARACTER(LEN=1000) :: head
		end subroutine loadMatrix
		subroutine append(A,B)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: A
			DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: B
		end subroutine append
	end interface
	
	Mb=(MLeo.ge.Minf .and. MLeo.le.Msup)
	call find(Mb,ix)
	
	if (allocated(ix)) then
		allocate(MeAv(size(ix)))
		MeAv=MLeo(ix)
		allocate(Endxi(size(ix))); allocate(Endxf(size(ix)))
		row0=0
		Endxi(1)=1
		do i=1,size(ix)
			Me=MLeo(ix(i))
			write(Mgstr,'(I3.3,F0.3)') int(Me),Me-int(Me)
			path=EvoZ//Mgstr//".txt"
			call loadMatrix(path,Ztevo0,head) !readCSV_h
			if (i.eq.1) then
				allocate(ZtevoTab(size(Ztevo0,1),size(Ztevo0,2)))
				ZtevoTab=Ztevo0
			else
				call append(ZtevoTab,Ztevo0)
			end if
			Endxf(i)=row0+size(Ztevo0,1)
			row0=Endxf(i)
			if (i.gt.1) then
				Endxi(i)=Endxf(i-1)+1
			end if
			deallocate(path); deallocate(Ztevo0)
		end do
	else
		print*,'from storeEvoZ: no mass values for element diffusion'
	end if

end subroutine

subroutine storeTracks_V_ZAMS(Minf,Msup,Mstep,MlowMS,Ztvec,nMt,fnd05t,Mavail,TrackTab,Tndxi,Tndxf, &
	& velTrackTab,Vndxi,Vndxf,ZAMStab,ZAndxi,ZAndxf)
	IMPLICIT NONE
	DOUBLE PRECISION Minf,Msup,Mstep,MlowMS
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: Ztvec
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Mavail,TrackTab,velTrackTab,ZAMStab
	INTEGER, DIMENSION(:), INTENT(in) :: nMt !size of Ztvec
	INTEGER, DIMENSION(:), ALLOCATABLE :: Vndxi,Vndxf,ZAndxi,ZAndxf
	INTEGER, DIMENSION(:,:), ALLOCATABLE :: Tndxi,Tndxf
	LOGICAL, DIMENSION(:), INTENT(in) :: fnd05t !size of Ztvec
	
	INTEGER Mcol,iz
	INTEGER row0T,row0V,row0Z,kM
	DOUBLE PRECISION MTr,M_in0,M_in
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Track0,velTrack0,ZAMS0
	CHARACTER(LEN=1000) head
	CHARACTER(LEN=56) :: tracceLeo="/home/bonfanti/Documents/PostDocLiegi/PARSEC1.2S/Tracce/"
	CHARACTER(LEN=:), ALLOCATABLE :: Isp,kv,velPath,ZAMSpath
	LOGICAL sameFile
	
	interface
		subroutine num2str(x,d,str)
			IMPLICIT NONE
			DOUBLE PRECISION x
			CHARACTER(LEN=:), ALLOCATABLE :: str !was 40
			INTEGER d
		end subroutine num2str
		subroutine openTracksPD(Z,M,PD,Tracks,M_in,Is)
			IMPLICIT NONE
			DOUBLE PRECISION Z,M,M_in
			INTEGER PD
			DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Tracks
			CHARACTER(LEN=:), ALLOCATABLE :: Is
		end subroutine
		subroutine append(A,B)
			IMPLICIT NONE
			DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE, INTENT(inout) :: A
			DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: B
		end subroutine append
		subroutine checkTracksPD(Z,M,PD,M_in0,M_in,sameFile,Is)
		  IMPLICIT NONE
		  DOUBLE PRECISION Z,M,M_in0,M_in
		  INTEGER PD
		  LOGICAL sameFile
		  CHARACTER(LEN=:), ALLOCATABLE :: Is
		end subroutine checkTracksPD
		subroutine loadMatrix(fileName,matrix,head)
			IMPLICIT NONE
			CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: fileName
			DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: matrix
			CHARACTER(LEN=1000) :: head
		end subroutine loadMatrix
	end interface
	
	print*,'begin storeTrack'
	
	Mcol=maxval(nMt) !columns of matrices containing the starting/ending
					!indices of the track tables
	
	allocate(Mavail(size(Ztvec),Mcol)); allocate(Tndxi(size(Ztvec),Mcol))
	allocate(Tndxf(size(Ztvec),Mcol))
	row0T=0
	Tndxi(1,1)=1
	!-
	allocate(Vndxi(size(Ztvec))); allocate(Vndxf(size(Ztvec)))
	row0V=0
	Vndxi(1)=1
	!-
	allocate(ZAndxi(size(Ztvec))); allocate(ZAndxf(size(Ztvec)))
	row0Z=0
	ZAndxi(1)=1
		
	!store relevant tracks in the big matrix TrackTab
	do iz=1,size(Ztvec)					
		MTr=Minf
		kM=0
		do while (MTr.lt.Msup+Mstep) !.le.
			if (MTr.eq.Minf) then
				kM=kM+1
				call openTracksPD(Ztvec(iz),MTr,1,Track0,M_in0,Isp)
				if (iz.eq.1) then
					allocate(TrackTab(size(Track0,1),size(Track0,2)))
					TrackTab=Track0
				else
					call append(TrackTab,Track0)
				end if
				Mavail(iz,kM)=M_in0
				Tndxf(iz,kM)=row0T+size(Track0,1)
				row0T=Tndxf(iz,kM)
				if (.not.(iz.eq.1 .and. kM.eq.1)) then
					if (kM.eq.1) then !beginning of a new line in the matrix
						Tndxi(iz,kM)=Tndxf(iz-1,nMt(iz-1))+1
					else
						Tndxi(iz,kM)=Tndxf(iz,kM-1)+1
					end if
				end if
				deallocate(Track0); deallocate(Isp)
			else
				call checkTracksPD(Ztvec(iz),MTr,1,M_in0,M_in,sameFile,Isp)
				if (.not.sameFile) then
					kM=kM+1
					call loadMatrix(Isp,Track0,head)
					call append(TrackTab,Track0)
					M_in0=M_in
					Mavail(iz,kM)=M_in0
					Tndxf(iz,kM)=row0T+size(Track0,1)
					row0T=Tndxf(iz,kM)
					if (kM.eq.1) then !beginning of a new line in the matrix. Should never happen here
						Tndxi(iz,kM)=Tndxf(iz-1,nMt(iz-1))+1
					else
						Tndxi(iz,kM)=Tndxf(iz,kM-1)+1
					end if
					deallocate(Track0); deallocate(Isp)
				end if
			end if
									
			MTr=MTr+Mstep
		end do
		if (.not.fnd05t(iz)) then
			kM=kM+1
			call openTracksPD(Ztvec(iz),MlowMS,1,Track0,M_in0,Isp)
			call append(TrackTab,Track0)
			Mavail(iz,kM)=M_in0
			Tndxf(iz,kM)=row0T+size(Track0,1)
			row0T=Tndxf(iz,kM)
			if (kM.eq.1) then !beginning of a new line in the matrix. Should never happen here
				Tndxi(iz,kM)=Tndxf(iz-1,nMt(iz-1))+1
			else
				Tndxi(iz,kM)=Tndxf(iz,kM-1)+1
			end if
			deallocate(Track0); deallocate(Isp)
		end if
							
		call num2str(Ztvec(iz),5,kv)
		!Charge also tables of evolutionary speed
		velPath=tracceLeo//"TracksVel/Z"//kv//".txt"
		call loadMatrix(velPath,velTrack0,head)
		if (iz.eq.1) then
			allocate(velTrackTab(size(velTrack0,1),size(velTrack0,2)))
			velTrackTab=velTrack0
		else
			call append(velTrackTab,velTrack0)
		end if
		Vndxf(iz)=row0V+size(velTrack0,1)
		row0V=Vndxf(iz)
		if (iz.gt.1) then
			Vndxi(iz)=Vndxf(iz-1)+1
		end if
		deallocate(velPath); deallocate(velTrack0)
		
		!Charge also ZAMS tables
		ZAMSpath=tracceLeo//"ZAMSdata/Z"//kv//".txt"
		call loadMatrix(ZAMSpath,ZAMS0,head)
		if (iz.eq.1) then
			allocate(ZAMStab(size(ZAMS0,1),size(ZAMS0,2)))
			ZAMStab=ZAMS0
		else
			call append(ZAMStab,ZAMS0)
		end if
		ZAndxf(iz)=row0Z+size(ZAMS0,1)
		row0Z=ZAndxf(iz)
		if (iz.gt.1) then
			ZAndxi(iz)=ZAndxf(iz-1)+1
		end if
		deallocate(ZAMSpath); deallocate(ZAMS0)
	end do
	
	print*,'end storeTrack'

end subroutine

subroutine sort0(arr)
	!deals with INTEGER
	IMPLICIT NONE
	INTEGER M,NSTACK
	INTEGER, DIMENSION(:), INTENT(inout) :: arr
	PARAMETER (M=7,NSTACK=50)
!	Sorts an array arr(1:n) into ascending numerical order using the Quicksort algorithm. n
!	is input; arr is replaced on output by its sorted rearrangement.
!	Parameters: M is the size of subarrays sorted by straight insertion and NSTACK is the required
!	auxiliary storage.
	INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
	INTEGER a,temp
	jstack=0
	l=1
	ir=size(arr) !n
1 	if (ir-l.lt.M) then 	!Insertion sort when subarray small enough.
		do j=l+1,ir
			a=arr(j)
			do i=j-1,l,-1
				if(arr(i).le.a)goto 2
				arr(i+1)=arr(i)
			end do
			i=l-1
2		 	arr(i+1)=a
		end do
		if (jstack.eq.0) return
		ir=istack(jstack) 	!Pop stack and begin a new round of partitioning.
		l=istack(jstack-1)
		jstack=jstack-2
	else
		k=(l+ir)/2 			
		!Choose median of left, center, and right elements as partitioning
		!element a. Also rearrange so that a(l)<=a(l+1)<=a(ir).
		temp=arr(k)
		arr(k)=arr(l+1)
		arr(l+1)=temp
		if(arr(l).gt.arr(ir))then
			temp=arr(l)
			arr(l)=arr(ir)
			arr(ir)=temp
		end if
		if(arr(l+1).gt.arr(ir))then
			temp=arr(l+1)
			arr(l+1)=arr(ir)
			arr(ir)=temp
		end if
		if(arr(l).gt.arr(l+1))then
			temp=arr(l)
			arr(l)=arr(l+1)
			arr(l+1)=temp
		end if
		i=l+1 				!Initialize pointers for partitioning.
		j=ir
		a=arr(l+1) 			!Partitioning element.
3	 	continue 			!Beginning of innermost loop.
		i=i+1 				!Scan up to find element > a.
		if(arr(i).lt.a)goto 3
4	 	continue
		j=j-1 				!Scan down to find element < a.
		if(arr(j).gt.a)goto 4
		if(j.lt.i)goto 5 	!Pointers crossed. Exit with partitioning complete.
		temp=arr(i) 	 	!Exchange elements.
		arr(i)=arr(j)
		arr(j)=temp
		goto 3 				!End of innermost loop.
5	 	arr(l+1)=arr(j) 	!Insert partitioning element.
		arr(j)=a
		jstack=jstack+2
		!Push pointers to larger subarray on stack, process smaller subarray immediately.
		if (jstack.gt.NSTACK) then
			print*, 'NSTACK too small in sort'
		end if
		if (ir-i+1.ge.j-l) then
			istack(jstack)=ir
			istack(jstack-1)=i
			ir=j-1
		else
			istack(jstack)=j-1
			istack(jstack-1)=l
			l=i
		end if
	end if
	goto 1
end subroutine sort0

subroutine sort1(arr)
	!deals with DOUBLE
	IMPLICIT NONE
	INTEGER M,NSTACK
	DOUBLE PRECISION, DIMENSION(:), INTENT(inout) :: arr
	PARAMETER (M=7,NSTACK=50)
!	Sorts an array arr(1:n) into ascending numerical order using the Quicksort algorithm. n
!	is input; arr is replaced on output by its sorted rearrangement.
!	Parameters: M is the size of subarrays sorted by straight insertion and NSTACK is the required
!	auxiliary storage.
	INTEGER i,ir,j,jstack,k,l,istack(NSTACK)
	DOUBLE PRECISION a,temp !REAL
	jstack=0
	l=1
	ir=size(arr) !n
1 	if (ir-l.lt.M) then 	!Insertion sort when subarray small enough.
		do j=l+1,ir
			a=arr(j)
			do i=j-1,l,-1
				if(arr(i).le.a)goto 2
				arr(i+1)=arr(i)
			end do
			i=l-1
2		 	arr(i+1)=a
		end do
		if (jstack.eq.0) return
		ir=istack(jstack) 	!Pop stack and begin a new round of partitioning.
		l=istack(jstack-1)
		jstack=jstack-2
	else
		k=(l+ir)/2 			
		!Choose median of left, center, and right elements as partitioning
		!element a. Also rearrange so that a(l)<=a(l+1)<=a(ir).
		temp=arr(k)
		arr(k)=arr(l+1)
		arr(l+1)=temp
		if(arr(l).gt.arr(ir))then
			temp=arr(l)
			arr(l)=arr(ir)
			arr(ir)=temp
		end if
		if(arr(l+1).gt.arr(ir))then
			temp=arr(l+1)
			arr(l+1)=arr(ir)
			arr(ir)=temp
		end if
		if(arr(l).gt.arr(l+1))then
			temp=arr(l)
			arr(l)=arr(l+1)
			arr(l+1)=temp
		end if
		i=l+1 				!Initialize pointers for partitioning.
		j=ir
		a=arr(l+1) 			!Partitioning element.
3	 	continue 			!Beginning of innermost loop.
		i=i+1 				!Scan up to find element > a.
		if(arr(i).lt.a)goto 3
4	 	continue
		j=j-1 				!Scan down to find element < a.
		if(arr(j).gt.a)goto 4
		if(j.lt.i)goto 5 	!Pointers crossed. Exit with partitioning complete.
		temp=arr(i) 	 	!Exchange elements.
		arr(i)=arr(j)
		arr(j)=temp
		goto 3 				!End of innermost loop.
5	 	arr(l+1)=arr(j) 	!Insert partitioning element.
		arr(j)=a
		jstack=jstack+2
		!Push pointers to larger subarray on stack, process smaller subarray immediately.
		if (jstack.gt.NSTACK) then
			print*, 'NSTACK too small in sort'
		end if
		if (ir-i+1.ge.j-l) then
			istack(jstack)=ir
			istack(jstack-1)=i
			ir=j-1
		else
			istack(jstack)=j-1
			istack(jstack-1)=l
			l=i
		end if
	end if
	goto 1
end subroutine sort1

subroutine indexI(arr,indx) !n,
	IMPLICIT NONE
	INTEGER, DIMENSION(:), INTENT(in) :: arr !(n)
	INTEGER indx(size(arr)),M,NSTACK !n,indx(n)
	
	PARAMETER (M=7,NSTACK=50)
	!Indexes an array arr(1:n), i.e., outputs the array indx(1:n) such that arr(indx(j))
	!is in ascending order for j = 1;2; : : :;N. The input quantities n and arr are not changed
	INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK),n
	DOUBLE PRECISION a
	
	n=size(arr)
	do j=1,n
		indx(j)=j
	end do
	jstack=0
	l=1
	ir=n
1 	if(ir-l.lt.M)then
		do j=l+1,ir
			indxt=indx(j)
			a=arr(indxt)
			do i=j-1,l,-1
				if(arr(indx(i)).le.a)goto 2
				indx(i+1)=indx(i)
			end do
			i=l-1 !!!
2 			indx(i+1)=indxt
		end do
		if(jstack.eq.0)return
		ir=istack(jstack)
		l=istack(jstack-1)
		jstack=jstack-2
	else
		k=(l+ir)/2
		itemp=indx(k)
		indx(k)=indx(l+1)
		indx(l+1)=itemp
		if(arr(indx(l)).gt.arr(indx(ir)))then
			itemp=indx(l)
			indx(l)=indx(ir)
			indx(ir)=itemp
		end if
		if(arr(indx(l+1)).gt.arr(indx(ir)))then
			itemp=indx(l+1)
			indx(l+1)=indx(ir)
			indx(ir)=itemp
		end if
		if(arr(indx(l)).gt.arr(indx(l+1)))then
			itemp=indx(l)
			indx(l)=indx(l+1)
			indx(l+1)=itemp
		end if
		i=l+1
		j=ir
		indxt=indx(l+1)
		a=arr(indxt)
3		continue
			i=i+1
		if(arr(indx(i)).lt.a)goto 3
4 		continue
			j=j-1
		if(arr(indx(j)).gt.a)goto 4
		if(j.lt.i)goto 5
		itemp=indx(i)
		indx(i)=indx(j)
		indx(j)=itemp
		goto 3
5 		indx(l+1)=indx(j)
		indx(j)=indxt
		jstack=jstack+2
		if (jstack.gt.NSTACK) then
			print*, 'NSTACK too small in indexI'
		end if
		if(ir-i+1.ge.j-l)then
			istack(jstack)=ir
			istack(jstack-1)=i
			ir=j-1
		else
			istack(jstack)=j-1
			istack(jstack-1)=l
			l=i
		end if
	end if
	goto 1
end subroutine indexI

subroutine indexx(arr,indx)
	IMPLICIT NONE
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: arr
	INTEGER indx(size(arr)),M,NSTACK
	
	PARAMETER (M=7,NSTACK=50)
	!Indexes an array arr(1:n), i.e., outputs the array indx(1:n) such that arr(indx(j))
	!is in ascending order for j = 1;2; : : :;N. The input quantities n and arr are not changed
	INTEGER i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK),n
	DOUBLE PRECISION a
	
	n=size(arr)
	do j=1,n
		indx(j)=j
	end do
	jstack=0
	l=1
	ir=n
1 	if(ir-l.lt.M)then
		do j=l+1,ir
			indxt=indx(j)
			a=arr(indxt)
			do i=j-1,l,-1
				if(arr(indx(i)).le.a)goto 2
				indx(i+1)=indx(i)
			end do
			i=l-1 !!!
2 			indx(i+1)=indxt
		end do
		if(jstack.eq.0)return
		ir=istack(jstack)
		l=istack(jstack-1)
		jstack=jstack-2
	else
		k=(l+ir)/2
		itemp=indx(k)
		indx(k)=indx(l+1)
		indx(l+1)=itemp
		if(arr(indx(l)).gt.arr(indx(ir)))then
			itemp=indx(l)
			indx(l)=indx(ir)
			indx(ir)=itemp
		end if
		if(arr(indx(l+1)).gt.arr(indx(ir)))then
			itemp=indx(l+1)
			indx(l+1)=indx(ir)
			indx(ir)=itemp
		end if
		if(arr(indx(l)).gt.arr(indx(l+1)))then
			itemp=indx(l)
			indx(l)=indx(l+1)
			indx(l+1)=itemp
		end if
		i=l+1
		j=ir
		indxt=indx(l+1)
		a=arr(indxt)
3		continue
			i=i+1
		if(arr(indx(i)).lt.a)goto 3
4 		continue
			j=j-1
		if(arr(indx(j)).gt.a)goto 4
		if(j.lt.i)goto 5
		itemp=indx(i)
		indx(i)=indx(j)
		indx(j)=itemp
		goto 3
5 		indx(l+1)=indx(j)
		indx(j)=indxt
		jstack=jstack+2
		if (jstack.gt.NSTACK) then
			print*, 'NSTACK too small in indexx'
		end if
		if(ir-i+1.ge.j-l)then
			istack(jstack)=ir
			istack(jstack-1)=i
			ir=j-1
		else
			istack(jstack)=j-1
			istack(jstack-1)=l
			l=i
		end if
	end if
	goto 1
end subroutine indexx

subroutine trovaIndici(v,ndxi,ndxf)

	!Given a vector v, containing ordered (and possible multiple) values such as
	!v=[a a a b b c c c c c c] 
	!it gives two vectors ndxi e ndxf 
	!ndxi=indices corresponding to the first occurrence of x
	!ndxf=indices corresponding to the last occurrence of x

	IMPLICIT NONE
	! Arguments declarations 
	DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: v
	INTEGER, DIMENSION(:), ALLOCATABLE :: ndxi,ndxf
	! Variable declarations 
	INTEGER k,ki,kf
	INTEGER, DIMENSION(:), ALLOCATABLE :: tmpi,tmpf
	DOUBLE PRECISION lastv
	
	interface
		subroutine allocVecI(vTmp,n,v)
			IMPLICIT NONE
			INTEGER, DIMENSION(:), INTENT(in) :: vTmp
			INTEGER, DIMENSION(:), ALLOCATABLE :: v
			INTEGER n
		end subroutine allocVecI
	end interface


	allocate(tmpi(size(v)))
	allocate(tmpf(size(v)))
	!
	ki=1
	kf=1
	if (size(v)>2) then
		lastv=v(size(v))
		tmpi(ki)=1
		k=1
		do while (v(k).ne.lastv) 
			do while (v(k)==v(k+1)) 
				k=k+1
			end do 
			tmpf(kf)=k
			kf=kf+1
			ki=ki+1
			tmpi(ki)=k+1
			k=k+1
		end do 
		tmpf(kf)=size(v)
	else if (size(v)==2) then
		if (v(1)==v(2)) then 
			tmpi(ki)=1
			tmpf(kf)=2
		else
			tmpi(ki)=1
			tmpf(kf)=1
			ki=ki+1
			kf=kf+1
			tmpi(ki)=2
			tmpf(kf)=2
		end if 
	else
		tmpi(ki)=1
		tmpf(kf)=1
	end if

	call allocVecI(tmpi,ki,ndxi)
	call allocVecI(tmpf,ki,ndxf)

end subroutine trovaIndici

subroutine allocVecI(vTmp,n,v)
	!Consider a vTmp vector that has already been allocated with a size that, in general, is >= than the number (n)
	!of elements that are actually hosted (a priori the number of elements entering vTmp was not known)
	!This subroutine writes the elements of vTmp into a new v vector, which has the proper size
	!By definition size(v)<=size(vTmp) 

	IMPLICIT NONE
	INTEGER, DIMENSION(:), INTENT(in) :: vTmp
	INTEGER, DIMENSION(:), ALLOCATABLE :: v
	INTEGER n
	
	INTEGER i
	
	allocate(v(n))
	do i=1,n
		v(i)=vTmp(i)
	end do

end subroutine allocVecI  

subroutine loadMatrix(fileName,matrix,head)
	IMPLICIT NONE
	CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: fileName
	DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: matrix
	CHARACTER(LEN=1000) :: head
	
	CHARACTER(LEN=1000) line
	CHARACTER :: ch1='#',ch2='%'
	INTEGER, DIMENSION(:), ALLOCATABLE :: ndx
	INTEGER f_unit,iol
	INTEGER k,k1
	INTEGER i,kx
	INTEGER nRow,nCol
	LOGICAL skip0,skip
	
	interface
		subroutine countRow(filename,nRow,k1)
			IMPLICIT NONE
			CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: filename
			INTEGER nRow,k1
		end subroutine countRow
	end interface
	
	call countRow(fileName,nRow,k1)
	allocate(ndx(nRow))
	
	call get_unit(f_unit)
	OPEN(UNIT=f_unit, FILE=fileName, STATUS='OLD')
	skip0=.false. !becomes true if enters in a comment block
	k=0		!counts total lines of file
	kx=0
	READ(f_unit, '(A)', IOSTAT=iol) line
	do while(iol.eq.0)
		k=k+1
		call skipLines(line,ch1,ch2,skip) !trim()
		if (.not.skip) then
			kx=kx+1
			ndx(kx)=k
		end if
		if (k1.eq.1) then
			head='No HEADER is present'
		else
			if (k.eq.k1-1) then
				head=line !take header if present
			end if
		end if
		if (k.eq.k1) then !take just once the 1st useful line to detect the number of columns
			call countCol(line,nCol)
		end if
		READ(f_unit, '(A)', IOSTAT=iol) line
	end do
	CLOSE(f_unit)
	
	allocate(matrix(nRow,nCol))
		
	call get_unit(f_unit)
	OPEN(UNIT=f_unit, FILE=fileName, STATUS='OLD')
	kx=1
	do i=1,k
		if (kx.le.nRow) then !.lt.
			if (i.eq.ndx(kx)) then
				READ(f_unit,*) matrix(kx,:)
				kx=kx+1
			else
				READ(f_unit,*) line		
			end if
		end if
	end do
	CLOSE(f_unit)
	
end subroutine loadMatrix

subroutine countRow(filename,nRow,k1)
	IMPLICIT NONE
	CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: filename
	INTEGER nRow,k1
	
	CHARACTER(LEN=1000) :: line
	CHARACTER(1) :: ch1='#',ch2='%'
	INTEGER f_unit,ios,iol
	INTEGER k,kc,kbc
	LOGICAL skip,skip0
	
	
	call get_unit(f_unit)
	OPEN(UNIT=f_unit, FILE=filename, STATUS='OLD', IOSTAT=ios)
	skip0=.false. !becomes true if enters in a comment block
	k=0		!counts total lines of file
	kc=0	!counts number of commented lines
	kbc=0	!counts number of comment blocks
	k1=0	!number of the 1st line without comments
	if (ios.eq.0) then
		READ(f_unit, '(A)', IOSTAT=iol) line
		do while(iol.eq.0)
			k=k+1
			call skipLines(trim(line),ch1,ch2,skip)
			if (skip) then !skip==.true.
				kc=kc+1
			else
				if (k1.eq.0) then !first line without comments
					k1=k !number of the first line without comment
				end if
			end if
			if (skip.neqv.skip0) then !change between comment/not_comment
				if (.not.skip0) then !skip0==false. We're passing from not_comment to comment -> block comment found!
					kbc=kbc+1
				end if
				skip0=.not.skip0
			end if
			READ(f_unit, '(A)', IOSTAT=iol) line
		end do
	else
		print*, "from CountRow: File does not exist"
		CLOSE(f_unit)
		return
	end if
	CLOSE(f_unit)
	
	nRow=k-kc
end subroutine countRow

subroutine countCol(line,nCol)
	IMPLICIT NONE
	CHARACTER(LEN=1000) :: line
	CHARACTER :: TAB=achar(9)
	INTEGER nCol
	INTEGER i,j,endl
	
	i=1
	do while(((line(i:i).eq.' ').or.(line(i:i).eq.TAB)).and.(i.lt.len(line))) !
		i=i+1
	end do
	
	nCol=0
	j=i
	endl=len_trim(line)-(j-1) !line(endl:endl) is for sure a character (no blank)
	do while(j.le.endl) !len_trim(line). Inequality holds if len_trim(line)>=1
		if ((line(j:j).eq.' ').or.(line(j:j).eq.TAB).or.(line(j:j).eq.',')) then
			nCol=nCol+1
			j=j+1
			do while(((line(j:j).eq.' ').or.(line(j:j).eq.TAB)).and.(j.lt.endl)) !
				j=j+1
				if (line(j:j).eq.',') then
					nCol=nCol+1
					j=j+1
				end if
			end do
		end if
		j=j+1
	end do
	nCol=nCol+1
	
end subroutine countCol

subroutine skipLines(line,ch1,ch2,skip)
	IMPLICIT NONE
	CHARACTER ch1,ch2
	CHARACTER(LEN=1000) line
	CHARACTER :: TAB=achar(9)
	INTEGER i
	LOGICAL skip
	
	i=1
	do while(((line(i:i).eq.' ').or.(line(i:i).eq.TAB)).and.(i.lt.len(trim(line))))
		i=i+1
	end do
	
	if ((line(i:i).eq.ch1).or.(line(i:i).eq.ch2)) then
		skip=.true.
	else
		skip=.false.
	end if
	
end subroutine skipLines


subroutine num2str(x,d,str)
	IMPLICIT NONE
	DOUBLE PRECISION x,lx
	CHARACTER(LEN=:), ALLOCATABLE :: str !was 40
	CHARACTER(LEN=10) frmt
	INTEGER i,d,u,t
	
	lx=log10(x)
	if (lx.ge.1) then
		u=floor(lx)+1
	else
		u=1
	end if
	
	t=u+1+d
!	print*, t
	
	str="000000000000000000000000000000000000000000"	
!	str=repeat(achar(48)//achar(48)//achar(48)//achar(48)//achar(48)//achar(48),t)
	!Because of the bug that repeat one single character for 7 times even if t>7,
	!in this way if t>7 str will be made of 42 characters
	!	allocate(CHARACTER(len=lungh) :: str) --> the bug in the compiler doesn't
	!allow this instruction. You could use it only if lungh is a number and not a
	!variable
!	Update 7/3/18: if inside another subroutine, repeat command *sometimes* doesn't
!	work, so forget about it
	
	write(frmt,'(A2,I0,A,I0,A)') "(F",t,".",d,")"
		
	write(str,frmt) x
	
	i=t
	do while(str(i:i).eq.'0')
		i=i-1
	end do
	str=str(1:i)

end subroutine num2str

subroutine setNaN(x)
	IMPLICIT NONE
	DOUBLE PRECISION x
	
	x=0
	x=x/x
	
end subroutine setNaN

subroutine get_unit ( iunit )

!*****************************************************************************80
!
!! GET_UNIT returns a free FORTRAN unit number.
!
!  Discussion:
!
!    A "free" FORTRAN unit number is an integer between 1 and 99 which
!    is not currently associated with an I/O device.  A free FORTRAN unit
!    number is needed in order to open a file with the OPEN command.
!
!    If IUNIT = 0, then no free FORTRAN unit could be found, although
!    all 99 units were checked (except for units 5, 6 and 9, which
!    are commonly reserved for console I/O).
!
!    Otherwise, IUNIT is an integer between 1 and 99, representing a
!    free FORTRAN unit.  Note that GET_UNIT assumes that units 5 and 6
!    are special, and will never return those values.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    26 October 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Output, integer ( kind = 4 ) IUNIT, the free unit number.
!
  implicit none

  integer ( kind = 4 ) i
  integer ( kind = 4 ) ios
  integer ( kind = 4 ) iunit
  logical lopen

  iunit = 0

  do i = 1, 99

    if ( i /= 5 .and. i /= 6 .and. i /= 9 ) then

      inquire ( unit = i, opened = lopen, iostat = ios )

      if ( ios == 0 ) then
        if ( .not. lopen ) then
          iunit = i
          return
        end if
      end if

    end if

  end do

  return
end subroutine get_unit

