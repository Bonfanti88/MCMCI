PROGRAM MCMCI
	
	USE SunPD, sunra=>RSunm, sunmass=>MSunkg
	USE IsoPD, only : BmVinf,BmVsup,Zi_inf,Zi_sup
	USE empRel, only : cJ
  !************************************************************************************
  !1. VARIABLES DECLARATION                                                           *
  !************************************************************************************

  IMPLICIT NONE

  INTEGER :: abs_errorm,abs_errorp,ord_errorm,ord_errorp,ncol_file,abs_val,ord_val
  INTEGER :: nbinh,nu,i,j,ntr,na,nburn,i2,dof,ntiming,k,test,mco,mco2,l,nb,loc,nchain
  INTEGER :: ntimingoc, en_int3sig,nenoch2,limmed(5),nb2,n2,io,nbliss,k3,l3       
  INTEGER :: nrv,nat,lim1,lim2,ce,t,link,isze,statlen,soluce,numbin,l2,nbliss2
  INTEGER :: npla,ordertrendrv,medi,be_int1sig,be_int3sig,en_int1sig,nbox,nbox2
  INTEGER :: medis,limmeds(5)
  INTEGER :: ngroup=0,tmax=0,statlen2,nttvmax,epochtime
  INTEGER :: nsystr,order,gibbscount,nburn2,ppdum,nubin,di,nmpix,nmpiy,k2,dtred
  INTEGER, PARAMETER :: ne3=50                    
  INTEGER :: nfi=0,nphotot=0,nmax=0,nc=5,nd=2,ne=4,njump,njump_rv=0,ncra=0
  INTEGER :: nrvtot=0,nrvmax=0,nphototbin=0,nenoch=265,testrampe=0    
  INTEGER :: nsysglo=0,nsysparmax=81,nrvparglo=3,ne2=33,testpld=0
  INTEGER :: nrvsysparmax=21,tgus=1,nxd=0,nddf=0,testsin,testflare
  INTEGER :: test_kb=0,test_t0=0,test_per=0,test_b=0,test_dF=0,test_tidel=0,test_dur=0
  INTEGER :: row,rowTot,acc,sismo,iz,xZi,xZu,jj,idum,ntra,idCol,d,dN,dar,dimZ_iso,iol

  INTEGER, DIMENSION(:), ALLOCATABLE :: np,nprv,IArr,jumpgibbs,fwhmyorder,nf,testphef      
  INTEGER, DIMENSION(:), ALLOCATABLE :: epoch,epoch_oc,gibbsonoff,skyorder,nep
  INTEGER, DIMENSION(:), ALLOCATABLE :: length,flagtr,ndot,ppforder,grouporder,group
  INTEGER, DIMENSION(:), ALLOCATABLE :: sinusnumber,ramporder,jumporder,fwhmorder,pldcor
  INTEGER, DIMENSION(:), ALLOCATABLE :: colororder,timeorder,pporder,rvloghkorder
  INTEGER, DIMENSION(:), ALLOCATABLE :: rvtimeorder,rvbisorder,rvfwhmorder,fwhmxorder
  INTEGER, DIMENSION(:), ALLOCATABLE :: offsetnumber,rvcontrastorder,tflick,jumpnumber 
  INTEGER, DIMENSION(:), ALLOCATABLE :: flarenumber,nttv
  INTEGER, DIMENSION(:), ALLOCATABLE :: Zndxi,Zndxf,nM,nMt,iZt,indxZt,Vndxi,Vndxf
  INTEGER, DIMENSION(:), ALLOCATABLE :: ZAndxi,ZAndxf,Endxi,Endxf,rowIO
  INTEGER, DIMENSION(:), ALLOCATABLE :: nda
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: nlast,nbin,np_bin,npjh,xbox,ybox,epochtr
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: Tndxi,Tndxf

  DOUBLE PRECISION :: min_x,min_y,max_x,max_y,medi_val,trho,strho,sigddf,redfdur
  DOUBLE PRECISION :: m2,dub,delta_1,feg,fdeg,ratio2,phef,phef1,corbliss,ratiog
  DOUBLE PRECISION :: incli,ave,c1,c2,c3,c4,q1,q2,test2,dum4,lc1,success,rms,dif,mdum
  DOUBLE PRECISION :: starmass,dmass_s,dmass_s_ini,starradius,dstarradius,dstarradius_ini
  DOUBLE PRECISION :: delta_2,delta_3,thet,delta,kk,rs,ud,us,rossi,chisq,eg0,rp_iron
  DOUBLE PRECISION :: fseg,fteg,eg1,ratio,nmeanbin,harvest,r_R,bcor,resave,pond
  DOUBLE PRECISION :: dbeta,x_se,z_se,x_tr,z_tr,nacce,y_tr,y_se,temp_sys,dvsini
  DOUBLE PRECISION :: cons1,phef2,tempu,phef3,vrotsini,nincli,pos_y,pos_z,flaremodel
  DOUBLE PRECISION :: meX,meX2,meY,meY2,tie,thrown,resavebin,rhotep
  DOUBLE PRECISION :: ini_svsinisinbeta,dsvsinisinbeta,ini_svsinicosbeta
  DOUBLE PRECISION :: dsvsinicosbeta,binsize,phase,phaset,phase2,emsrsb
  DOUBLE PRECISION :: dsvsinicosbeta_ini,dsvsinisinbeta_ini,initsvsinisinb
  DOUBLE PRECISION :: Eano_se,Tano_se,resmeanerror,as1,as2,das,dis,dis2
  DOUBLE PRECISION :: ttiming,kepler,initsvsinicosb,tidelmax,redup,reddown,deltared
  DOUBLE PRECISION :: bigbetared,dum_beta_red,Mano,Eano,Tano,pos_x,dtime
  DOUBLE PRECISION :: aic,bic,dic2,bf_1,bf_1b,massmeanerror
  DOUBLE PRECISION :: meancotime,blissmodel,msrsb,dstarlum,bf_val,dage
  DOUBLE PRECISION :: teff,feh,massenoch,maxadapt,prot,burnfrac,test2b
  DOUBLE PRECISION :: rpup,rpdo,time,dtemp_s,dmet_s,dcol_s,cortide,mulimb0g
  DOUBLE PRECISION :: mulimb0t,tdbtra,rvmodplat,sf2,codilu,midx,midy
  DOUBLE PRECISION :: enoch_pmin,enoch_mmin,enoch_mmax,trho2,strho2,starlum
  DOUBLE PRECISION :: profullsecond=0.,prosecond=0.,bestmerit=1.E6,regul=1.         
  DOUBLE PRECISION :: meanposterior=0.,protransit=0.,profulltransit=0.
  DOUBLE PRECISION :: npara,ndata=0.,dic=0.,inclistar=90, phefsave
  DOUBLE PRECISION :: dumnum1,dumnum2,red_dur(ne3),limpd(8)
  DOUBLE PRECISION :: f01,f02,f03,radq,rade,dcon,kePM,Dm2_m1,massTmp,dmassTmp,radiusTmp,dradiusTmp
  DOUBLE PRECISION :: daRdP,daRdW,daRddF,daRdb,daRde,daRdom,Da_R
  DOUBLE PRECISION :: Drho,Drho1,rhoI,DrhoI,loggI,I_loggI,vsiniI,I_vsiniI,protI,I_protI
  DOUBLE PRECISION :: logRHKI,I_logRHKI,YMgI,I_YMgI,numaxI,I_numaxI,deltanuI,I_deltanuI
  DOUBLE PRECISION :: Zinf,Zsup,Mstep,Minf,Msup
  DOUBLE PRECISION :: MlowMS=0.5,rhotepi,Drhotepi,rhoInp,I_rhoInp,Rinp,I_Rinp
  DOUBLE PRECISION :: mag,I_mag,col,I_col,dist,I_dist,Hipf,useColor,Mabs
  DOUBLE PRECISION :: colI,I_colI,TeffI,I_TeffI,colL,TeffL,I_TeffL,Lum,I_Lum
  DOUBLE PRECISION :: mass_Isoch,dmass_Isoch,radius_Isoch,dradius_Isoch
  DOUBLE PRECISION :: start,finish

  DOUBLE PRECISION, PARAMETER :: pi=3.14159265358979,expi=2.71828183
  DOUBLE PRECISION, PARAMETER :: light=299792458.,earthra=6378000.  	
  DOUBLE PRECISION, PARAMETER :: ua=149597.8707E6,jupra=69911000.
  DOUBLE PRECISION, PARAMETER :: consbf=2.5066   !sunra=6.9598E8,
  DOUBLE PRECISION, PARAMETER :: jupmass=1898.6E24,earthmass=5.9742E24
  DOUBLE PRECISION, PARAMETER :: gravi=6.67384E-11 !,sunmass=1.96215E30
  DOUBLE PRECISION, PARAMETER :: FeHinf=log10(Zi_inf)+costZ
  DOUBLE PRECISION, PARAMETER :: FeHsup=log10(Zi_sup)+costZ
  DOUBLE PRECISION, DIMENSION(1) :: theo
  DOUBLE PRECISION, DIMENSION(2) :: rold_ini,qtemp,sqtemp
  DOUBLE PRECISION, DIMENSION(4) :: ntemp,sntemp
  DOUBLE PRECISION, DIMENSION(7) :: spec,sspec,enochcoef
  DOUBLE PRECISION, DIMENSION(13) :: eno
  DOUBLE PRECISION, DIMENSION(29) :: star
  DOUBLE PRECISION, DIMENSION(2) :: vg,I_vg
  DOUBLE PRECISION, DIMENSION(3) :: vrho,I_vrho
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Z_iso

  DOUBLE PRECISION :: gelman_svsinicosbeta,gelman_svsinisinbeta,gelman_rs,gelman_ms
  DOUBLE PRECISION :: gelman_f2,gelmanval,gelman_teff,gelman_met,gelman_col
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: gelman_ttv,ttv_ini
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gelman_dF,gelman_b,gelman_t0
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gelman_dur,gelman_per,gelman_kb
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gelman_secosw,ttvrms
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: gelman_ddf
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gelman_tidel,gelman_dfgroup
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gelman_sesinw
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: gelman_dFsec
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gelman_phampli1,gelman_phampli2
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gelman_phampli3,gelman_phoffset
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: gelman_t1ramp,gelman_t2ramp     
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: gelman_sit0,gelman_sip
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: gelman_flampli,gelman_fltau,gelman_flt0
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: gelman_jumplimb
  CHARACTER :: gelman

  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: permin,permax,kmax
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: temp1D,acf,cortime,testiron
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: edfgroup,dfgroup_ini,edfgroup_ini
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: binoverphot,binoverrv,f2,ktide
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: f2_end,ktide_end,sigusti
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: nsyspar,Mano_tr, Mano_se,ltt
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dilu,edilu,dilution,dummass
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: beta_red_glo,jitterglo,avusti
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: er_om,e_ka,initsecosw,m2_m1,temp3
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: initsesinw,timing,stiming
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: y,z,sig,w,a,timerit
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: timingoc,stimingoc,merit_end
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: omctiming,phampli1_ini,regul2    
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phampli2_ini,phampli3_ini
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: phoffset_ini,met_s,met_s_end
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: sesinw_ini,t0_ini,Tano_tr
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: svsinicosbeta_end,svsinisinbeta_end
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: vsini_end,beta_end,Eano_tr,temp_end 
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: glorvbjd,glorvresi,glorverror,rho_end
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: lum_s,lum_s_end,temp_s,temp_s_end
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: col_s,col_s_end
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Tperi,per_ini,mass_s_end,pepoch,tidel_ini
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ephampli3,ephoffset,ephampli1
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ephampli2,sb,gbeta_red
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: radius_s_end,beta,sdF
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: logg,logg_end,bf_resrms,age,age_end
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: best_beta_red,beta_red,gjitter
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mulimb_qd,z2,rvz,massjitter,massjitter_end
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rossiter,ephampli1_ini
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: fbjd,fphot,ferror,fphotcor,ferrorcor
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: fresi,fmodel,fmodel_tr
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: frvbjd,frv,frverror,frv2
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: frvresi,frvmodel,frvmodel2, trvmodel2
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: tbjd,tphot,terror,tphotcor,terrorcor
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: tresi,tmodel,tmodel_tr,tresibin,trv2
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: muph,t1ramp_ini,et1ramp_ini
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: t2ramp_ini,et2ramp_ini,et1ramp,et2ramp
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: trvbjd,trv,trverror,trvresi,trvmodel  
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: vsini,sper,sdur,st0,dkb,svsinisinbeta,stidel
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: svsinicosbeta,rephotmerit,photmerit
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rervmerit,rvmerit,remerit,merit,kb_ini
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: photmerit_end,rvmerit_end,remerit_end
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ave_presglo,sig_presglo,radius_s
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ave_presbinglo,sig_presbinglo,jitter
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ave_presglorv,sig_presglorv,rho,b_ini 
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: mass_s,srold,dsesinw,dsecosw
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dsesinw_ini,dsecosw_ini,sdF_ini,sb_ini
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: sdur_ini,sper_ini,st0_ini,dkb_ini,stidel_ini
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dF_ini,dur_ini,secosw_ini,nacib,nacinb
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: rvresrms,rvresrms_soluce,best_bred_dur
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ephampli2_ini,ephampli3_ini,bnpave,best_bnpave
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: ephoffset_ini,best_jitter,er_exc  
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: eddf,eddf_ini,tju,dFsec_ini,edFsec,edFsec_ini
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: incli_s,incli_s_end
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: sinincli_s,sinincli_s_end
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: best_bresrmsbin,bresrmsbin,bred_dur
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Zvec,Zt,Ztvec,MeAv

  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: fbjd_bin,fphot_bin,fphotcor_bin,dfgroup
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ferrorcor_bin,ferror_bin,fresi_bin,dfgroup_end
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: sit0_ini,esit0,esit0_ini,logg_p,logg_p_end
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: fltau_ini,efltau,efltau_ini,fltime
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: flampli_ini,eflampli,eflampli_ini
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: flt0_ini,eflt0,eflt0_ini,sttv,sttv_ini
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: sip_ini,esip,esip_ini,teq_p,teq_p_end
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: hillrad_p, hillrad_p_end
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: offsettime,resigroup
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: sysva,sysva_bestfit,ecosw_end,esinw_end
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: rvsysva,phampli2,phampli3
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: rvsysva_bestfit,npave,texp,rvtexp
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: irrad,irrad_end,compttv
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: u,v,cvm,po,octime_end,octime,mulimb1
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: exc_end,sesinw_end,phoffset,ecosw,esinw   
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: dF_end,b_end,t0_end,kb_end,ka_end
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: per_end,dur_end,secosw_end,b3_end,tidel_end
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: omega_end,rr_end,a_R_end,b2_end
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: rhop_end,ql_ini,nl_ini,dF,b,snl
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: prtr,proc,prtr_end,proc_end,syscor 
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: sql,eresibin,phampli1
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: mass_p_end,rvbis,rvloghk
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: bf_resrmsbin,resrmsbin,model_tr2 
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: rvtrendco,safro,safro_end,model_tr2bin
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: ddf,ddf_end,radipla,radipla_end,ddepth,ddepth_end
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: dratio,dratio_end,ttr,ttr_end,dFsec,dFsec_end
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: resrms     
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: rvtrendco_end,djumplimb_ini,rephotbin 
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: photbin,ephotbin,ephotcorbin
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: photcorbin,bjdbin,fobjdbin
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: photcor2bin,ephotcor2bin,mulimb0,mulimb2
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: enochm,enochr,enocht,enochf
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: modelbin,model_trbin,per,tidel,secosw
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: rvmod,rvmod2,radius_p_end
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: mass_p_sini,mass_p_sini_end
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: t1ramp,t2ramp,t1ramp_end,t2ramp_end
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: semi_end,inclian_end,photcor2
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: djumplimb,t0,rvfwhm,phampli1_end   
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ka,kb,a_R,b2,b3,rr,dur
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phampli2_end,phampli3_end
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: phoffset_end,inclian,rvchi2
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: sesinw,radius_p,mass_p,semi
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: roche,a_roche,roche_end,a_roche_end
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: rervchi2,photchi2,rephotchi2,rhop
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: rvcontrast,fobjd,bjd2
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: exc,omega,bjd,phot,error,resi
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: photcor,rold,dX,dY,airmass,fwhmy
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: rvbjd,rv,rv2,rverror,rvresi
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: mulimb_nl,fwhm,sky,fwhmx
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: forvbjd,forvbjd2,model,model_tr
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: fwhmbin,fwhmybin,skybin,fwhmxbin,airmassbin
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: xjh, yjh, fluxjh,dXbin,dYbin
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: cotime, maxt, mint
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: IsocTab,Mavail,TrackTab,ZtevoTab !Isoc0,Track0
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: velTrackTab,GyroTab,ZAMStab !velTrack0,ZAMS0
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: jumplimb_end,ql_end,nl,ql,ttv,ttv_end
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE:: rvmodpla,rv2pla,pixcon
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: jumplimb,resibin  
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: sit0,sip,sip_end,sit0_end
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: flampli,fltau,flampli_end,fltau_end
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: flt0,flt0_end

  CHARACTER :: errorx,errory,ifacf,jump,ironlow,isttv
  CHARACTER :: dum,let1,let2,let3,let4,let5,let6,ok,kesol,istiming,fitmsrs,enoch
  CHARACTER :: priorrad,priorvsini,gibbs_sampler,istimingoc,iftrendrv,fixstellar
  CHARACTER :: maxred,cone,newf,priorteff,priormet
  CHARACTER :: priorlogg,rossitif,isddf,stelincli,dotide,priorf2,testf2,m2op,isoch
  CHARACTER :: massfromr='n',dfcond,msrs,enoch_torres,priorrho,priorlum,priormass
  CHARACTER :: priorprot,priordeltanu,priornumax,priorlogRHK,priorYMg
  CHARACTER :: priormag,priorcol,priordist,MjumpIso,RjumpIso
  CHARACTER(LEN=2) :: limb,tefil,tefil2,temp_a1                                                
  CHARACTER(LEN=3) :: rampmod     
  CHARACTER(LEN=7) :: texunit
  CHARACTER(LEN=8) :: fin                                                                               
  CHARACTER(LEN=12) :: name,name2,name3,name4,name5                                                    
  CHARACTER(LEN=15) :: fil2,inputfile,graphfile2,acffile,graphacffile
  CHARACTER(LEN=16) :: temp_a2,graphfile          
  CHARACTER(LEN=20) :: texb,ylabel,xlabel,xlabel2
  CHARACTER(LEN=37) :: temp2                                                          
  CHARACTER(LEN=47) :: dum2           
  CHARACTER(LEN=58) :: temp
  CHARACTER(LEN=171) :: ltitle
  CHARACTER(LEN=1000) :: head
  CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE :: filter,ddf_filter,tilo,filo
  CHARACTER(LEN=2), DIMENSION(:), ALLOCATABLE :: wfilter
  CHARACTER(LEN=:), ALLOCATABLE :: Barnes
  CHARACTER, DIMENSION(:), ALLOCATABLE :: fitlimb,iffitoc,bliss,iffitph,utc_tdb
  CHARACTER, DIMENSION(:), ALLOCATABLE :: accepted,burned,overphot,overrv,oldform,donegroup        
  CHARACTER, DIMENSION(:,:), ALLOCATABLE :: isjump
  
  LOGICAL isEq,ybesR
  LOGICAL, DIMENSION(:), ALLOCATABLE :: fnd05,fnd05t

  EXTERNAL baseline,rvbaseline,rvglo,enochlaw
  
  interface
  	
  	function polyval(v,x)
		IMPLICIT NONE
		DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: v
		DOUBLE PRECISION x,polyval
	end function polyval
    
	subroutine loadMatrix(fileName,matrix,head)
	  IMPLICIT NONE
	  CHARACTER(LEN=:), ALLOCATABLE, INTENT(in) :: fileName
	  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: matrix
	  CHARACTER(LEN=1000) :: head
	end subroutine loadMatrix
	
	subroutine uniqueFast(list,d,indices,first)
		IMPLICIT NONE
		DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: list
		INTEGER d
		INTEGER, DIMENSION(:), ALLOCATABLE :: indices
		LOGICAL first
	end subroutine uniqueFast
	
	subroutine tracksNumber(Minf,Msup,Mstep,MlowMS,Zvec,fnd05,Zt,nM)
		IMPLICIT NONE
		DOUBLE PRECISION Minf,Msup,Mstep,MlowMS
		DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: Zvec
		DOUBLE PRECISION, DIMENSION(size(Zvec)) :: Zt
		INTEGER, DIMENSION(size(Zvec)) :: nM
		LOGICAL, DIMENSION(size(Zvec)) :: fnd05
	end subroutine tracksNumber
	
	subroutine storeIsoc(Zvec,IsocTab,Zndxi,Zndxf,idCol)
		IMPLICIT NONE
		INTEGER idCol
		DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: Zvec
		DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: IsocTab
		INTEGER, DIMENSION(:), ALLOCATABLE :: Zndxi,Zndxf
	end subroutine storeIsoc
	
	subroutine storeTracks_V_ZAMS(Minf,Msup,Mstep,MlowMS,Ztvec,nMt,fnd05t,Mavail,TrackTab, &
		& Tndxi,Tndxf,velTrackTab,Vndxi,Vndxf,ZAMStab,ZAndxi,ZAndxf)
		IMPLICIT NONE
		DOUBLE PRECISION Minf,Msup,Mstep,MlowMS
		DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: Ztvec
		DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: Mavail,TrackTab,velTrackTab,ZAMStab
		INTEGER, DIMENSION(:), INTENT(in) :: nMt !size of Ztvec
		INTEGER, DIMENSION(:), ALLOCATABLE :: Vndxi,Vndxf,ZAndxi,ZAndxf
		INTEGER, DIMENSION(:,:), ALLOCATABLE :: Tndxi,Tndxf
		LOGICAL, DIMENSION(:), INTENT(in) :: fnd05t !size of Ztvec
	end subroutine storeTracks_V_ZAMS
	
	subroutine storeEvoZ(Minf,Msup,ZtevoTab,Endxi,Endxf,MeAv)
		IMPLICIT NONE
		DOUBLE PRECISION Minf,Msup
		DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: MeAv
		DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: ZtevoTab
		INTEGER, DIMENSION(:), ALLOCATABLE :: Endxi,Endxf
	end subroutine storeEvoZ

	subroutine SCPmcmcPD12S(SCP,Intestaz,IsocTab,Zvec,Zndxi,Zndxf, &
	& TrackTab,nM,Mavail,indxZt,Tndxi,Tndxf,velTrackTab,Vndxi,Vndxf, &
	& GyroTab,ZAMStab,ZAndxi,ZAndxf,ZtevoTab,Endxi,Endxf,MeAv,M_star,I_M_star, &
	& R_star,I_R_star,Teff,I_Teff,t_star,I_t_star,L,I_L,link,row,acc)
	  IMPLICIT NONE
	  DOUBLE PRECISION :: M_star,I_M_star,R_star,I_R_star,Teff,I_Teff,t_star,I_t_star
	  DOUBLE PRECISION :: L,I_L
	  DOUBLE PRECISION, DIMENSION(29) :: SCP
	  DOUBLE PRECISION, DIMENSION(:), INTENT(in) :: Zvec,MeAv
	  DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: IsocTab,TrackTab,Mavail
	  DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: velTrackTab,GyroTab,ZAMStab
	  DOUBLE PRECISION, DIMENSION(:,:), INTENT(in) :: ZtevoTab
	  INTEGER link,row,acc
	  INTEGER, DIMENSION(:), INTENT(in) :: Zndxi,Zndxf,nM,indxZt,Vndxi,Vndxf
	  INTEGER, DIMENSION(:), INTENT(in) :: ZAndxi,ZAndxf,Endxi,Endxf
	  INTEGER, DIMENSION(:,:), INTENT(in) :: Tndxi,Tndxf
	  CHARACTER(LEN=171) :: Intestaz
	end subroutine SCPmcmcPD12S
	
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

  !******************************************************************************
  !2. FORMAT DEFINITION                                                         !
  !******************************************************************************

100 FORMAT(A6,F7.4,A4,F7.4) 
108 FORMAT(A21,1x,F10.4)
109 FORMAT(A21,1x,F10.4,A13)                   
110 FORMAT(I5,7(1x,F11.6))    
111 FORMAT(F15.6,8(1x,F15.6))                    
116 FORMAT(A21,ES14.5)
117 FORMAT(2x,I3,4x,F5.1,3x,F5.1,2x,F6.1,2x,F6.1)               
119 FORMAT(1x,I3,3x,F5.1,3x,F5.1,6(ES9.2))                 
121 FORMAT(I5,1x,F15.6,1x,F13.6,4(1x,F13.6)) 
122 FORMAT(F15.6,13(1x,F15.6))           
129 FORMAT(A21,F13.8,A49)   
132 FORMAT(A8,I2,A1,F12.4,2x,F12.4)

  !******************************************************************************
  !3. INPUT FILE READING + ARRAY ALLOCATION                                     !
  !******************************************************************************

  CALL RANDOM_SEED(SIZE=isze)
  ALLOCATE(IArr(isze))
  CALL RANDOM_SEED(PUT=IArr(1:isze))

  OPEN(UNIT=11,FILE='mcmc.dat')
  DO i=1,3
     READ(11,*) dum 
  ENDDO
  READ(11,*) dum2,nchain
  READ(11,*) dum2,nat  
  na = nat*nchain	  
  ALLOCATE(accepted(na));ALLOCATE(burned(na))
  rowTot=0 !counter of I/O consistency in models
  allocate(rowIO(na)) !vector of I/O consistency
  READ(11,*) dum2,burnfrac
  nburn = NINT((burnfrac/100.)*nat)
  READ(11,*) dum2,statlen
  READ(11,*) dum2,binsize
  READ(11,*) dum2,success
  success=success/100.
  IF(success.GE.1)THEN
     PRINT*, 'Error in mcmc.dat/success rate too large'; STOP
  ENDIF
  READ(11,*) dum2,kesol
  READ(11,*) dum2,kepler
  READ(11,*) dum2,gibbs_sampler,statlen2
  READ(11,*) dum2,maxadapt
  READ(11,*) dum2,gelman
  READ(11,*) dum2,ifacf
  READ(11,*) dum2,rossitif
  READ(11,*) dum2,dotide 
  njump = 3; npara = 3 !Stellar mass, teff & metallicity
  READ(11,*) dum2,fixstellar
  IF(fixstellar.EQ.'y')THEN
    njump=0; npara=0
  ENDIF
  READ(11,*) dum2,fitmsrs
  IF(fitmsrs.EQ.'y')THEN
    njump = njump+1; npara=npara+1
  ENDIF
  IF(fitmsrs.EQ.'y'.AND.fixstellar.EQ.'y')THEN
     PRINT*, 'Error in mcmc.dat/fixstellar and fitmsrs are turned on'; STOP
  ENDIF
  READ(11,*) dum2,msrs,msrsb,emsrsb
  IF(msrs.EQ.'y')THEN
    njump=njump-1; npara=npara-1
  ENDIF
  IF(fitmsrs.EQ.'y'.AND.msrs.EQ.'y')THEN
     PRINT*, 'Error in mcmc.dat/fitmsrs and msrs are turned on'; STOP
  ENDIF
  IF(fixstellar.EQ.'y'.AND.msrs.EQ.'y')THEN
     PRINT*, 'Error in mcmc.dat/fixstellar and msrs are turned on'; STOP
  ENDIF
  READ(11,*) dum2,massfromr
  IF(fitmsrs.EQ.'y'.AND.massfromr.EQ.'y')THEN
     PRINT*, 'Error in mcmc.dat/fitmsrs and massfromr are turned on'; STOP
  ENDIF
  IF(fixstellar.EQ.'y'.AND.massfromr.EQ.'y')THEN
     PRINT*, 'Error in mcmc.dat/fixstellar and massfromr are turned on'; STOP
  ENDIF
  IF(msrs.EQ.'y'.AND.massfromr.EQ.'y')THEN
     PRINT*, 'Error in mcmc.dat/msrs and massfromr are turned on'; STOP
  ENDIF
  READ(11,*) dum2,enoch
  IF(enoch.EQ.'y')THEN
    njump=njump-1; npara=npara-1
  ENDIF
  IF(fitmsrs.EQ.'y'.AND.enoch.EQ.'y')THEN
     PRINT*, 'Error in mcmc.dat/fitmsrs and enoch are turned on'; STOP
  ENDIF
  IF(fixstellar.EQ.'y'.AND.enoch.EQ.'y')THEN
     PRINT*, 'Error in mcmc.dat/fixstellar and enoch are turned on'; STOP
  ENDIF
  IF(massfromr.EQ.'y'.AND.enoch.EQ.'y')THEN
     PRINT*, 'Error in mcmc.dat/massfromr and enoch are turned on'; STOP
  ENDIF
  IF(msrs.EQ.'y'.AND.enoch.EQ.'y')THEN
     PRINT*, 'Error in mcmc.dat/msrs and enoch are turned on'; STOP
  ENDIF
  READ(11,*) dum2,enoch_pmin
  READ(11,*) dum2,enoch_mmin,enoch_mmax
  READ(11,*) dum2,enoch_torres
  IF(enoch.EQ.'y')THEN
    OPEN(UNIT=3,FILE='/home/bonfanti/Documents/PostDocLiegi/MCMC_Andrea/Lib/eb.dat')
    nenoch2=0
    DO i=1,nenoch
      test=1
      READ(3,*) eno(1:12),l2,eno(13)
      IF(eno(13).LT.enoch_pmin) test=0
      IF(eno(1).LT.enoch_mmin) test=0
      IF(eno(1).GT.enoch_mmax) test=0
      IF(enoch_torres.EQ.'y'.AND.l2.NE.1) test=0
      IF(test.EQ.1) nenoch2=nenoch2+1
    ENDDO
    PRINT*, 'Number of EB stars used for calibration law:',nenoch2
    IF(nenoch2.LT.1)THEN
      PRINT*, 'No EB stars meet the requirements'; STOP
    ENDIF
    ALLOCATE(enochm(nenoch2,2));ALLOCATE(enochr(nenoch2,2))
    ALLOCATE(enocht(nenoch2,2));ALLOCATE(enochf(nenoch2,2))
    ALLOCATE(massjitter(na))
    CLOSE(3)
    OPEN(UNIT=3,FILE='/home/bonfanti/Documents/PostDocLiegi/MCMC_Andrea/Lib/eb.dat')
    j=0
    DO i=1,nenoch 
      test=1
      READ(3,*) eno(1:12),l2,eno(13)
      IF(eno(13).LT.enoch_pmin) test=0
      IF(eno(1).LT.enoch_mmin) test=0
      IF(eno(1).GT.enoch_mmax) test=0
      IF(enoch_torres.EQ.'y'.AND.l2.NE.1) test=0
      IF(test.EQ.1)THEN
        j=j+1
        enochm(j,1)=eno(1)
        enochm(j,2)=eno(2)
        enochr(j,1)=eno(3)
        enochr(j,2)=eno(4)
        enocht(j,1)=eno(9)
        enocht(j,2)=eno(10)
        enochf(j,1)=eno(7)
        enochf(j,2)=eno(8)
      ENDIF
    ENDDO
  ENDIF
  MjumpIso='n'
  RjumpIso='n'
  read(11,*) dum2,isoch
  if (isoch.eq.'y') then
  	MjumpIso='y'
  end if
  if(fitmsrs.EQ.'y'.AND.isoch.EQ.'y')then
     print*, 'Error in mcmc.dat/fitmsrs and isoch are turned on'; STOP
  end if
  if(fixstellar.EQ.'y'.AND.isoch.EQ.'y')then
     print*, 'Error in mcmc.dat/fixstellar and isoch are turned on'; STOP
  end if
  if(massfromr.EQ.'y'.AND.isoch.EQ.'y')then
     print*, 'Error in mcmc.dat/massfromr and isoch are turned on'; STOP
  end if
  if(msrs.EQ.'y'.AND.isoch.EQ.'y')then
     print*, 'Error in mcmc.dat/msrs and isoch are turned on'; STOP
  end if
  if (enoch.eq.'y'.and.isoch.eq.'y') then
	print*, 'Error in mcmc.dat/enoch and isoch are turned on'; STOP
  end if
  
  READ(11,*) dum2,m2op
  READ(11,*) dum2,stelincli
  READ(11,*) dum2,dfcond
  READ(11,*) dum2,limb
  IF(limb.NE.'nl'.AND.limb.NE.'qd'.AND.limb.NE.'no')THEN
     PRINT*, 'Error in reading mcmc.dat/LD law'; STOP
  ENDIF
  READ(11,*) dum2,rampmod
  READ(11,*) dum2,isddf,sigddf
  READ(11,*) dum2,isttv
  READ(11,*) dum2,reddown,redup
  deltared=(redup-reddown)/DBLE(ne3-1)
  READ(11,*) dum2,maxred,redfdur
  test2=1E12
  DO i=1,ne3-1
     IF(i.EQ.1)THEN
       red_dur(i)=binsize
     ELSE IF(i.EQ.2)THEN
        red_dur(i)=reddown
     ELSE
       red_dur(i)=red_dur(i-1)+deltared
     ENDIF
  ENDDO
  IF(maxred.EQ.'y')THEN
    red_dur(ne3)=redfdur
    dtred=ne3
  ELSE
    red_dur(ne3)=red_dur(ne3-1)+deltared
  ENDIF  
  red_dur=red_dur/(24*60.)
  READ(11,*) dum2,nbliss
  READ(11,*) dum2,npla
  ALLOCATE(testiron(npla)); ALLOCATE(m2_m1(npla));ALLOCATE(temp3(npla))
  ALLOCATE(permin(npla)); ALLOCATE(permax(npla)); ALLOCATE(kmax(npla))
  READ(11,*) dum2,permin(1:npla)
  READ(11,*) dum2,permax(1:npla)
  READ(11,*) dum2,kmax(1:npla)
  READ(11,*) dum2,rpdo
  READ(11,*) dum2,rpup
  rpdo=rpdo*earthra/jupra;rpup=rpup*earthra/jupra
  READ(11,*) dum2,ironlow
  READ(11,*) dum2,tidelmax
  IF(isttv.EQ.'y')THEN
    ALLOCATE(nttv(npla))
    DO i=1,npla
      let1 = CHAR(48+int(i/10))
      let2 = CHAR(48+mod(i,10))
      name = 'ttv' // let1 // let2 // '.dat'
      nttv(i)=0 !. is an integer
      OPEN(UNIT=3,FILE=name)
      DO  
        READ(3,*,IOSTAT=io) dumnum1,dumnum2
        IF(io < 0)THEN
          WRITE(*,'(A13,I2,A2,I3)') 'N_TTV planet ',i,': ',nttv(i)
          EXIT
        ELSE
          nttv(i)=nttv(i)+1
          njump=njump+1
          npara=npara+1
        ENDIF 
      ENDDO
      CLOSE(3)
    ENDDO
    nttvmax=MAXVAL(nttv)
    ALLOCATE(sttv(npla,nttvmax));ALLOCATE(sttv_ini(npla,nttvmax))
    ALLOCATE(ttv(npla,nttvmax,na));ALLOCATE(ttv_ini(npla,nttvmax))
    ALLOCATE(epochtr(npla,nttvmax));ALLOCATE(ttr(npla,nttvmax,na))
    ALLOCATE(ttvrms(npla));ALLOCATE(compttv(npla,nttvmax))
    IF(gelman.EQ.'y'.AND.nchain.GT.1) ALLOCATE(gelman_ttv(npla,nttvmax))
    DO i=1,npla
      let1 = CHAR(48+int(i/10))
      let2 = CHAR(48+mod(i,10))
      name = 'ttv' // let1 // let2 // '.dat'
      OPEN(UNIT=3,FILE=name)
      DO j=1,nttv(i)
        READ(3,*) epochtr(i,j),ttv_ini(i,j),sttv_ini(i,j)
        sttv(i,j)=sttv_ini(i,j)
      ENDDO
      CLOSE(3)
    ENDDO
  ENDIF
 
  ! Data
  DO i=1,3
     READ(11,*) dum
  ENDDO
  READ(11,*) dum2,ntr 
  READ(11,*) dum
  IF(ntr.GT.0)THEN
     ! Allocation of arrays related to photometry baselines
     ALLOCATE(bf_resrms(ntr));ALLOCATE(oldform(ntr));ALLOCATE(utc_tdb(ntr))
     ALLOCATE(cortime(ntr));ALLOCATE(grouporder(ntr));ALLOCATE(group(ntr))
     ALLOCATE(bf_resrmsbin(ntr,ne3));ALLOCATE(resrmsbin(ntr,ne3))
     ALLOCATE(sysva(ntr,nsysparmax));ALLOCATE(sysva_bestfit(ntr,nsysparmax))
     ALLOCATE(sit0_ini(ntr,4));ALLOCATE(esit0(ntr,4));ALLOCATE(esit0_ini(ntr,4))
     ALLOCATE(sip_ini(ntr,4));ALLOCATE(esip(ntr,4));;ALLOCATE(esip_ini(ntr,4))
     ALLOCATE(fltau_ini(ntr,4));ALLOCATE(efltau(ntr,4));ALLOCATE(efltau_ini(ntr,4))
     ALLOCATE(flt0_ini(ntr,4));ALLOCATE(eflt0(ntr,4));;ALLOCATE(eflt0_ini(ntr,4))
     ALLOCATE(flampli_ini(ntr,4));ALLOCATE(eflampli(ntr,4));;ALLOCATE(eflampli_ini(ntr,4))
     ALLOCATE(nsyspar(ntr));ALLOCATE(sinusnumber(ntr));ALLOCATE(ppforder(ntr))
     ALLOCATE(ramporder(ntr));ALLOCATE(colororder(ntr));ALLOCATE(flarenumber(ntr))
     ALLOCATE(fwhmorder(ntr));ALLOCATE(pporder(ntr));ALLOCATE(fwhmyorder(ntr))
     ALLOCATE(timeorder(ntr));ALLOCATE(offsetnumber(ntr))
     ALLOCATE(tju(ntr,4))
     ALLOCATE(offsettime(ntr,4));ALLOCATE(binoverphot(ntr));ALLOCATE(overphot(ntr))    
     ALLOCATE(pldcor(ntr));ALLOCATE(jumporder(ntr));ALLOCATE(jumpnumber(ntr))
     ALLOCATE(pepoch(ntr))
     ALLOCATE(np(ntr));ALLOCATE(filter(ntr));ALLOCATE(t1ramp_ini(ntr))
     ALLOCATE(et1ramp_ini(ntr));ALLOCATE(t2ramp_ini(ntr));ALLOCATE(et2ramp_ini(ntr))
     ALLOCATE(et1ramp(ntr));ALLOCATE(et2ramp(ntr));ALLOCATE(skyorder(ntr))
     ALLOCATE(wfilter(ntr));ALLOCATE(nbin(ntr,ne3));ALLOCATE(fwhmxorder(ntr))
     ALLOCATE(np_bin(ntr,ne3));ALLOCATE(gbeta_red(ntr));ALLOCATE(bred_dur(ntr))
     ALLOCATE(bnpave(ntr));ALLOCATE(npave(ntr,ne3));ALLOCATE(nlast(ntr,ne3))
     ALLOCATE(photchi2(ntr,na)); ALLOCATE(rephotchi2(ntr,na));ALLOCATE(best_bnpave(ntr))
     ALLOCATE(resrms(ntr,na));ALLOCATE(bresrmsbin(ntr));ALLOCATE(best_bred_dur(ntr))
     ALLOCATE(beta_red(ntr)); ALLOCATE(best_beta_red(ntr));ALLOCATE(flagtr(ntr))
     ALLOCATE(best_bresrmsbin(ntr));ALLOCATE(bliss(ntr));ALLOCATE(ndot(ntr))
     testsin=0
     testflare=0
     testrampe=0
     DO i=1,ntr
        let1 = CHAR(48+int(i/1000))
        i2 = i - INT(i/1000)*1000
        let2 = CHAR(48+int(i2/100))
        i2 = i - INT(i/100)*100
        let3 = CHAR(48+int(i2/10))
        let4 = CHAR(48+mod(i2,10))
        fin =  let1 // let2 // let3 // let4 // '.txt'
        name = 'phot' // fin
        np(i)=0  !. is an integer
        OPEN(UNIT=3,FILE=name)
        DO  
          READ(3,*,IOSTAT=io) dumnum1,dumnum2
          IF(io < 0)THEN
            WRITE(*,'(A3,I2,A2,I6,A7)') 'LC-',i,': ',np(i),' points'
            EXIT
          ELSE
            np(i)=np(i)+1
           ENDIF 
        ENDDO
        CLOSE(3)
        nsyspar(i)=0
        flagtr(i)=0
        READ(11,*) dum,filter(i),timeorder(i),colororder(i),fwhmorder(i),  &
             & fwhmxorder(i),fwhmyorder(i),skyorder(i),pporder(i),ppforder(i),   & 
             & sinusnumber(i),ramporder(i), jumpnumber(i),jumporder(i), offsetnumber(i), &
             & pldcor(i),flarenumber(i),gbeta_red(i),overphot(i),   &
             & binoverphot(i),bliss(i),oldform(i),utc_tdb(i),cortime(i),group(i), &
             & grouporder(i)
        IF(overphot(i).EQ.'y') binoverphot(i)=binoverphot(i)/(24*3600.)
        nsystr=0
        IF(timeorder(i).GE.0)THEN
           IF(timeorder(i).GT.4)THEN
             PRINT*, 'mcmc.dat/time model has a too large order. Set to 4'; timeorder(i)=4
           ENDIF
           nsyspar(i) = nsyspar(i) + 1
           nsystr=nsystr+1
           npara = npara+timeorder(i)+1
        ENDIF 
        IF(grouporder(i).GE.0.AND.group(i).GT.0)THEN
           IF(grouporder(i).GT.4)THEN
             PRINT*, 'mcmc.dat/group model has a too large order. Set to 4'; grouporder(i)=4
           ENDIF
           ngroup = group(i)
           nsyspar(i) = nsyspar(i) + 1
           nsystr=nsystr+1
           npara = npara+grouporder(i)+1
        ENDIF
        IF(colororder(i).GT.0)THEN
           IF(colororder(i).GT.4)THEN
             PRINT*, 'mcmc.dat/color model has a too large order. Set to 4'; colororder(i)=4
           ENDIF
           nsyspar(i) = nsyspar(i) + 1
           npara=npara+colororder(i)+1
           IF(nsystr.GT.0)THEN
              npara=npara-1
           ENDIF
           nsystr=nsystr+1
        ENDIF
        IF(fwhmorder(i).GT.0)THEN
           IF(fwhmorder(i).GT.4)THEN
             PRINT*, 'mcmc.dat/FWHM model has a too large order. Set to 4'; fwhmorder(i)=4
           ENDIF
           nsyspar(i) = nsyspar(i) + 1
           npara=npara+fwhmorder(i)+1
           IF(nsystr.GT.0)THEN
              npara=npara-1
           ENDIF
           nsystr=nsystr+1
        ENDIF
        IF(fwhmxorder(i).GT.0)THEN
           IF(fwhmxorder(i).GT.4)THEN
             PRINT*, 'mcmc.dat/FWHM_x model has a too large order. Set to 4'; fwhmxorder(i)=4
           ENDIF
           nsyspar(i) = nsyspar(i) + 1
           npara=npara+fwhmxorder(i)+1
           IF(nsystr.GT.0)THEN
              npara=npara-1
           ENDIF
           nsystr=nsystr+1
        ENDIF
        IF(fwhmyorder(i).GT.0)THEN
           IF(fwhmyorder(i).GT.4)THEN
             PRINT*, 'mcmc.dat/FWHM_y model has a too large order. Set to 4'; fwhmyorder(i)=4
           ENDIF
           nsyspar(i) = nsyspar(i) + 1
           npara=npara+fwhmyorder(i)+1
           IF(nsystr.GT.0)THEN
              npara=npara-1
           ENDIF
           nsystr=nsystr+1
        ENDIF
        IF(skyorder(i).GT.0)THEN
           IF(skyorder(i).GT.4)THEN
             PRINT*, 'mcmc.dat/background model has a too large order. Set to 4'; skyorder(i)=4
           ENDIF
           nsyspar(i) = nsyspar(i) + 1
           npara=npara+skyorder(i)+1
           IF(nsystr.GT.0)THEN
              npara=npara-1
           ENDIF
           nsystr=nsystr+1
        ENDIF
        IF(pporder(i).GT.0)THEN
           IF(pporder(i).GT.4)THEN
             PRINT*, 'mcmc.dat/x-y model has a too large order. Set to 4'; pporder(i)=4
           ENDIF
           IF(pporder(i).EQ.1)THEN
              ppdum = 3
           ELSE IF(pporder(i).EQ.2)THEN
              ppdum = 6 
           ELSE IF(pporder(i).EQ.3)THEN 
              ppdum = 10
           ELSE IF(pporder(i).EQ.4)THEN
              ppdum = 15
           ELSE
              PRINT*, 'Error in mcmc.dat/order of pp-function'; STOP
           ENDIF
           nsyspar(i) = nsyspar(i) + 1
           npara=npara+ppdum
           IF(nsystr.GT.0)THEN
              npara=npara-1
           ENDIF
           nsystr=nsystr+1
        ENDIF
        IF(ppforder(i).GT.0)THEN
           IF(ppforder(i).GT.2)THEN
             PRINT*, 'mcmc.dat/x-y-FWHM model has a too large order. Set to 1'; ppforder(i)=1
           ENDIF
           IF(ppforder(i).EQ.1)THEN
              ppdum = 6
           ELSE
              PRINT*, 'Error in mcmc.dat/order of x-y-FWHM function'; STOP
           ENDIF
           nsyspar(i) = nsyspar(i) + 1
           npara=npara+ppdum
           IF(nsystr.GT.0)THEN
              npara=npara-1
           ENDIF
           nsystr=nsystr+1
        ENDIF
        IF(sinusnumber(i).GT.0)THEN
           IF(testsin.LT.1)THEN
              ALLOCATE(sit0(ntr,4,na))
              ALLOCATE(sip(ntr,4,na)) 
              testsin=2
           ENDIF
           njump = njump + 2*sinusnumber(i)
           nsyspar(i) = nsyspar(i) + sinusnumber(i)
           npara = npara+3*sinusnumber(i)+1
           IF(nsystr.GT.0)THEN
              npara=npara-1
           ENDIF
           nsystr = nsystr + 1
        ENDIF
        IF(sinusnumber(i).GT.0)THEN
          k=sinusnumber(i)
          READ(11,*) dum,sip_ini(i,1:k),esip_ini(i,1:k)
        ENDIF
        IF(ramporder(i).GT.0)THEN
           IF(ramporder(i).GT.2)THEN
             PRINT*, 'mcmc.dat/ramp model has a too large order. Set to 2'; ramporder(i)=2
           ENDIF
           nsyspar(i) = nsyspar(i) + 1
           IF(ramporder(i).EQ.1.AND.rampmod.EQ.'log') npara=npara+2
           IF(ramporder(i).EQ.1.AND.rampmod.EQ.'exp')THEN
             npara=npara+3 
             READ(11,*) dum,t1ramp_ini(i),et1ramp_ini(i)
             IF(testrampe.LE.1) testrampe=1
           ENDIF
           IF(ramporder(i).EQ.2.AND.rampmod.EQ.'log') npara=npara+3
           IF(ramporder(i).EQ.2.AND.rampmod.EQ.'exp')THEN
             npara=npara+5
             READ(11,*) dum,t1ramp_ini(i),et1ramp_ini(i),t2ramp_ini(i),et2ramp_ini(i)
             testrampe=2
           ENDIF
           IF(nsystr.GT.0)THEN
              npara=npara-1
           ENDIF
           nsystr=nsystr+1
        ENDIF
        IF(i.EQ.ntr.AND.testrampe.GT.0)THEN
          ALLOCATE(t1ramp(ntr,na))
        ENDIF
        IF(i.EQ.ntr.AND.testrampe.GT.1)THEN
          ALLOCATE(t2ramp(ntr,na))
        ENDIF
        IF(jumpnumber(i).GT.0.AND.jumporder(i).GT.0)THEN
          k = jumpnumber(i)
          READ(11,*) dum,tju(i,1:k)
          IF(utc_tdb(i).EQ.'y')THEN
            DO j=1,k
              dum4=tju(i,j)
              tju(i,j)=tdbtra(dum4)
            ENDDO 
          ENDIF
          IF(jumporder(i).GT.3)THEN
           PRINT*, 'mcmc.dat/jump model has a too large order. Set to 3'; jumporder(i)=2
          ENDIF
          npara = npara + jumporder(i) + 1
          nsyspar(i) = nsyspar(i) + 1
          IF(nsystr.GT.0)THEN
            npara=npara-1
          ENDIF
          nsystr=nsystr+1
        ENDIF
        IF(offsetnumber(i).GT.0)THEN
          k=offsetnumber(i)
          READ(11,*) dum,offsettime(i,1:k)
          IF(utc_tdb(i).EQ.'y')THEN
            DO j=1,k
              dum4=offsettime(i,j)
              offsettime(i,j)=tdbtra(dum4)
            ENDDO 
          ENDIF
          nsyspar(i) = nsyspar(i) + 1
          npara = npara + offsetnumber(i)
          nsystr=nsystr+1
        ENDIF
        IF(pldcor(i).GT.0)THEN
          IF(pldcor(i).GT.2)THEN
            WRITE(*,*) 'PLD order greater than 2, set to 2' 
            pldcor(i)=2
          ENDIF
          testpld = 1
          nsyspar(i) = nsyspar(i) + 1
          npara = npara + pldcor(i)
          IF(pldcor(i).GT.1) npara = npara+17
          nsystr=nsystr+1
          IF(nsystr.GT.0)THEN
            npara=npara-1
          ENDIF
        ENDIF   
        IF(flarenumber(i).GT.0)THEN
           IF(testflare.LT.1)THEN
              ALLOCATE(fltau(ntr,4,na))
              ALLOCATE(flampli(ntr,4,na))  
              ALLOCATE(flt0(ntr,4,na)) 
              testflare=2
           ENDIF
           njump = njump + 3*flarenumber(i)
           npara = npara+3*flarenumber(i)
        ENDIF
        IF(flarenumber(i).GT.0)THEN
          k=flarenumber(i)
          READ(11,*) dum,flt0_ini(i,1:k),eflt0_ini(i,1:k),flampli_ini(i,1:k),eflampli_ini(i,1:k), &
            & fltau_ini(i,1:k),efltau_ini(i,1:k)
          IF(utc_tdb(i).EQ.'y')THEN
            DO j=1,k
              dum4=flt0_ini(i,j)
              flt0_ini(i,j)=tdbtra(dum4)
            ENDDO 
          ENDIF
        ENDIF
        IF(nsyspar(i).GT.0) nsysglo=nsysglo+1
        IF(isddf.NE.'n')THEN
          nxd=nxd+1
          DO j=1,i-1
            IF(filter(j).EQ.filter(i))THEN
               nxd=nxd-1
               EXIT
             ENDIF
          ENDDO
        ENDIF
        IF(nmax.LT.np(i)) nmax=np(i)
        nphotot=nphotot + np(i)
        nfi=nfi+1   
        wfilter(nfi)=filter(i)
        IF(i.GT.1)THEN
           DO j=1,i-1
              IF(filter(i).EQ.filter(j))THEN
                 nfi=nfi-1
                 EXIT
              ENDIF
           ENDDO
        ENDIF
     ENDDO
     IF(ngroup.GT.0)THEN
       ALLOCATE(dfgroup(ngroup,na));ALLOCATE(dfgroup_ini(ngroup))
       ALLOCATE(edfgroup(ngroup));ALLOCATE(edfgroup_ini(ngroup))
       IF(gelman.EQ.'y'.AND.nchain.GT.1) ALLOCATE(gelman_dfgroup(ngroup))
       DO i=1,ngroup
         njump=njump+1; npara=npara+1
         l=0
         DO j=1,ntr
           IF(group(j).EQ.i)THEN
             l=l+1
             IF(l.EQ.1)THEN
               t=np(j)
               IF(t.GT.tmax) tmax=t
             ELSE
               IF(np(j).NE.t)THEN
                 WRITE(*,'(A6,I2,A29)') 'Group ',i,': number of measurements vary'
                 STOP
               ENDIF
             ENDIF 
           ENDIF
         ENDDO
       ENDDO
       ALLOCATE(resigroup(ngroup,tmax));ALLOCATE(donegroup(ngroup))
       resigroup=0
       donegroup='n'
     ENDIF
     IF(ntr.GT.0.AND.nchain.GT.1.AND.testsin.GT.0)THEN
        ALLOCATE(gelman_sip(ntr,4));ALLOCATE(gelman_sit0(ntr,4))
     ENDIF
     IF(ntr.GT.0.AND.nchain.GT.1.AND.testflare.GT.0)THEN
        ALLOCATE(gelman_flampli(ntr,4));ALLOCATE(gelman_fltau(ntr,4))
        ALLOCATE(gelman_flt0(ntr,4))
     ENDIF
     IF(ntr.GT.0.AND.nchain.GT.1.AND.testrampe.GT.0)THEN
        ALLOCATE(gelman_t1ramp(ntr))
     ENDIF
     IF(ntr.GT.0.AND.nchain.GT.1.AND.testrampe.GT.1)THEN
        ALLOCATE(gelman_t2ramp(ntr))
     ENDIF
     IF(nxd.GT.1.AND.ntr.GT.0.AND.isddf.NE.'n')THEN
        nddf=nxd-1
        njump=njump+nddf*npla
        npara=npara+nddf*npla
        IF(isddf.EQ.'p') ndata=ndata+nddf*npla
        ALLOCATE(ddf(npla,nddf,na));ALLOCATE(eddf(npla,nddf))
        ALLOCATE(radipla(npla,nddf,na));ALLOCATE(ddepth(npla,nddf,na))
        ALLOCATE(dratio(npla,nddf,na))
        ALLOCATE(eddf_ini(npla,nddf));ALLOCATE(ddf_filter(nddf))
        IF(gelman.EQ.'y'.AND.nchain.GT.1) ALLOCATE(gelman_ddf(npla,nddf))
        k=0
        l=0
        DO i=1,ntr
           IF(nddf.GT.0.AND.isddf.NE.'n')THEN
              k=k+1
              test=0
              IF(k.EQ.1)THEN
                 test=1
              ELSE
                 DO j=1,i-1
                    IF(filter(i).EQ.filter(j))THEN
                       test=1
                       flagtr(i)=flagtr(j)
                    ENDIF
                 ENDDO
              ENDIF
              IF(test.EQ.0)THEN
                 l=l+1
                 flagtr(i)=l
                 ddf_filter(l)=filter(i)
              ENDIF
           ENDIF
        ENDDO
     ENDIF
     IF(testpld.GT.0) ALLOCATE(pixcon(ntr,nmax,9))
     ALLOCATE(phot(ntr,nmax)); ALLOCATE(bjd(ntr,nmax))
     ALLOCATE(xbox(ntr,nmax)); ALLOCATE(ybox(ntr,nmax))
     ALLOCATE(texp(ntr,nmax)); ALLOCATE(photcor2(ntr,nmax))
     ALLOCATE(bjd2(ntr,nmax)); ALLOCATE(photcor(ntr,nmax))
     ALLOCATE(fobjd(ntr,nmax)); ALLOCATE(error(ntr,nmax))
     ALLOCATE(resi(ntr,nmax));ALLOCATE(model(ntr,nmax))
     ALLOCATE(model_tr(ntr,nmax));ALLOCATE(fwhmy(ntr,nmax))
     ALLOCATE(model_tr2(ntr,nmax));ALLOCATE(testphef(nfi))
     ALLOCATE(dX(ntr,nmax)); ALLOCATE(dY(ntr,nmax))
     ALLOCATE(airmass(ntr,nmax));ALLOCATE(fwhm(ntr,nmax))
     ALLOCATE(dFsec(nfi,na,npla)); ALLOCATE(edFsec(nfi,npla))
     ALLOCATE(edFsec_ini(nfi,npla));ALLOCATE(sky(ntr,nmax))
     ALLOCATE(phampli1(nfi,na)); ALLOCATE(ephampli1(nfi))
     ALLOCATE(phampli2(nfi,na)); ALLOCATE(ephampli2(nfi))
     ALLOCATE(phampli3(nfi,na)); ALLOCATE(ephampli3(nfi))
     ALLOCATE(phoffset(nfi,na)); ALLOCATE(ephoffset(nfi))
     ALLOCATE(ephampli1_ini(nfi));ALLOCATE(ephampli2_ini(nfi))
     ALLOCATE(dilu(nfi));ALLOCATE(fwhmx(ntr,nmax))
     ALLOCATE(ephampli3_ini(nfi));ALLOCATE(ephoffset_ini(nfi))
     ALLOCATE(edilu(nfi));ALLOCATE(dilution(nfi))
     ALLOCATE(iffitoc(nfi));ALLOCATE(iffitph(nfi))
     ALLOCATE(dFsec_ini(nfi,npla))
     ALLOCATE(phampli1_ini(nfi));ALLOCATE(phampli2_ini(nfi))
     ALLOCATE(phampli3_ini(nfi));ALLOCATE(phoffset_ini(nfi))
     IF(nchain.GT.1.AND.gelman.EQ.'y')THEN
       ALLOCATE(gelman_dFsec(nfi,npla))
       ALLOCATE(gelman_phampli1(nfi));;ALLOCATE(gelman_phampli2(nfi))
       ALLOCATE(gelman_phampli3(nfi));ALLOCATE(gelman_phoffset(nfi))
     ENDIF
  ENDIF
  READ(11,*) dum2,nrv,ordertrendrv
  IF(dotide.EQ.'y'.AND.nrv.GT.0)THEN 
     testf2 = 'y'
  ELSE IF(dotide.EQ.'y')THEN
     testf2='n'
     WRITE(*,*) 'mcmc.dat/tides model not included as no RV'
  ENDIF
  IF(ordertrendrv.GT.0.AND.nrv.GT.0)THEN
    IF(ordertrendrv.GT.2)THEN
      PRINT*, 'mcmc.dat/RV global trend order has a too large order. Set to 2'; ordertrendrv=2
    ENDIF
    npara=npara+ordertrendrv
    iftrendrv='y'
    ALLOCATE(rvtrendco(na,nrvparglo))
  ENDIF
  READ(11,*) dum
  IF(nrv.GT.0)THEN
     ALLOCATE(rvtimeorder(nrv));ALLOCATE(rvfwhmorder(nrv));ALLOCATE(rvbisorder(nrv))
     ALLOCATE(rvcontrastorder(nrv));ALLOCATE(binoverrv(nrv));ALLOCATE(overrv(nrv))
     ALLOCATE(rvchi2(nrv,na)); ALLOCATE(rervchi2(nrv,na));ALLOCATE(rvresrms(nrv))
     ALLOCATE(jitter(nrv));ALLOCATE(best_jitter(nrv));ALLOCATE(rvresrms_soluce(nrv))
     ALLOCATE(nprv(nrv));ALLOCATE(gjitter(nrv));ALLOCATE(rvloghkorder(nrv))
     OPEN(UNIT=447,FILE='rvsyste.res')
     DO i=1,nrv
        let1 = CHAR(48+int(i/1000))
        i2 = i - INT(i/1000)*1000
        let2 = CHAR(48+int(i2/100))
        i2 = i - INT(i/100)*100
        let3 = CHAR(48+int(i2/10))
        let4 = CHAR(48+mod(i2,10))
        fin =  let1 // let2 // let3 // let4 // '.txt'
        name = 'rv' // fin
        nprv(i)=0 !. is an integer
        OPEN(UNIT=3,FILE=name)
        DO  
          READ(3,*,IOSTAT=io) dumnum1,dumnum2
          IF(io < 0)THEN
            WRITE(*,'(A3,I2,A2,I6,A7)') 'RV-',i,': ',nprv(i),' points'
            EXIT
          ELSE
            nprv(i)=nprv(i)+1
           ENDIF 
        ENDDO
        CLOSE(3)
        READ(11,*) dum,gjitter(i),rvtimeorder(i),rvfwhmorder(i),rvbisorder(i), &
        & rvcontrastorder(i),rvloghkorder(i),overrv(i),binoverrv(i)  

        IF(overrv(i).EQ.'y') binoverrv(i)=binoverrv(i)/(24*3600.)
        npara = npara+1
        IF(rvtimeorder(i).GT.0)THEN
          IF(rvtimeorder(i).GT.4)THEN
            PRINT*, 'mcmc.dat/rv time correction has a too large order. Set to 4'; rvtimeorder(i)=4
          ENDIF
          npara = npara+rvtimeorder(i)
        ENDIF
        IF(rvfwhmorder(i).GT.0)THEN
          IF(rvfwhmorder(i).GT.4)THEN
             PRINT*, 'mcmc.dat/rv fwhm correction has a too large order. Set to 4'; rvfwhmorder(i)=4
          ENDIF
          npara = npara+rvfwhmorder(i)
        ENDIF
        IF(rvbisorder(i).GT.0)THEN
          IF(rvbisorder(i).GT.4)THEN
            PRINT*, 'mcmc.dat/rv BIS correction has a too large order. Set to 4'; rvbisorder(i)=4
          ENDIF
          npara = npara+rvbisorder(i)
        ENDIF
        IF(rvloghkorder(i).GT.0)THEN
          IF(rvloghkorder(i).GT.4)THEN
            PRINT*, 'mcmc.dat/rv RHK correction has a too large order. Set to 4'; rvloghkorder(i)=4
          ENDIF
          npara = npara+rvloghkorder(i)
        ENDIF
        IF(rvcontrastorder(i).GT.0)THEN
          IF(rvcontrastorder(i).GT.4)THEN
             PRINT*, 'mcmc.dat/rv contrast correction has a too large order. Set to 4'; rvcontrastorder(i)=4
          ENDIF
          npara = npara+rvcontrastorder(i)
        ENDIF
        IF(nrvmax.LT.nprv(i)) nrvmax=nprv(i)
        nrvtot = nrvtot + nprv(i)
     ENDDO
  ENDIF
  ndata=ndata+DBLE(nrvtot+nphotot)
  READ(11,*) dum2, istiming
  IF(istiming.EQ.'y')THEN
     ntiming=0
     OPEN(UNIT=3,FILE='timing.dat')
     DO
       READ(3,*,IOSTAT=io) dumnum1,dumnum2
       IF(io < 0)THEN
         WRITE(*,'(A13,I3)') 'Ntrtimings = ',ntiming
         EXIT
       ELSE
         ntiming=ntiming+1
       ENDIF
     ENDDO
     CLOSE(3)
     ndata=ndata+DBLE(ntiming)
     ALLOCATE(timerit(na))
     ALLOCATE(timing(ntiming))
     ALLOCATE(stiming(ntiming))
     ALLOCATE(epoch(ntiming))
     ALLOCATE(omctiming(ntiming))
     OPEN(UNIT=3,FILE='timing.dat')
     DO i=1,ntiming
        READ(3,*) timing(i),stiming(i)
     ENDDO
     CLOSE(3)
  ENDIF
  READ(11,*) dum2,istimingoc
  IF(istimingoc.EQ.'y')THEN
     ntimingoc=0
     OPEN(UNIT=3,FILE='timingoc.dat')
     DO
       READ(3,*,IOSTAT=io) dumnum1,dumnum2
       IF(io < 0)THEN
         WRITE(*,'(A13,I3)') 'Noctimings = ',ntimingoc
         EXIT
       ELSE
         ntimingoc=ntimingoc+1
       ENDIF
     ENDDO
     CLOSE(3)
     ALLOCATE(timingoc(ntimingoc))
     ALLOCATE(stimingoc(ntimingoc))
     ALLOCATE(epoch_oc(ntimingoc))
     OPEN(UNIT=3,FILE='timingoc.dat')
     DO i=1,ntimingoc
        READ(3,*) timingoc(i),stimingoc(i)
     ENDDO
     CLOSE(3)
  ENDIF

  ALLOCATE(isjump(npla,9));ALLOCATE(kb_ini(npla)) 
  ALLOCATE(logg_p(npla,na));ALLOCATE(teq_p(npla,na))
  ALLOCATE(hillrad_p(npla,na))
  ALLOCATE(dur_ini(npla));ALLOCATE(b_ini(npla));ALLOCATE(dF_ini(npla))
  ALLOCATE(rold(na,nd)); ALLOCATE(srold(nd))
  ALLOCATE(vsini(na));ALLOCATE(secosw_ini(npla));ALLOCATE(sesinw_ini(npla))
  ALLOCATE(per_ini(npla));ALLOCATE(t0_ini(npla));ALLOCATE(tidel_ini(npla))
  ALLOCATE(nacib(nchain)); ALLOCATE(nacinb(nchain))
  ALLOCATE(beta(na));ALLOCATE(octime(npla,na))
  ALLOCATE(svsinisinbeta(na)); ALLOCATE(svsinicosbeta(na))
  ALLOCATE(beta_red_glo(na));ALLOCATE(jitterglo(na))
  ALLOCATE(ql(nfi,nd,na));ALLOCATE(ql_ini(nfi,nd)); ALLOCATE(sql(nfi,nd))
  ALLOCATE(jumplimb(nfi,nd,na)); ALLOCATE(djumplimb(nfi,nd))
  ALLOCATE(fitlimb(nfi)); ALLOCATE(djumplimb_ini(nfi,nd))
  ALLOCATE(nl(nfi,ne,na));ALLOCATE(nl_ini(nfi,ne));ALLOCATE(snl(nfi,ne))
  ALLOCATE(dF(npla,na));ALLOCATE(sdF(npla));ALLOCATE(sdF_ini(npla))
  ALLOCATE(prtr(npla,na));ALLOCATE(proc(npla,na))
  ALLOCATE(mass_s(na)); ALLOCATE(lum_s(na));ALLOCATE(mass_p(npla,na))
  ALLOCATE(mass_p_sini(npla,na));ALLOCATE(temp_s(na));ALLOCATE(met_s(na))
  ALLOCATE(radius_s(na)); ALLOCATE(radius_p(npla,na));allocate(col_s(na))
  ALLOCATE(logg(na)); ALLOCATE(irrad(npla,na));allocate(age(na))
  ALLOCATE(inclian(npla,na)); ALLOCATE(semi(npla,na))
  ALLOCATE(roche(npla,na));ALLOCATE(a_roche(npla,na))
  ALLOCATE(sb_ini(npla));ALLOCATE(st0_ini(npla))
  ALLOCATE(sdur_ini(npla));ALLOCATE(sper_ini(npla));ALLOCATE(stidel_ini(npla))
  ALLOCATE(b(npla,na)); ALLOCATE(sb(npla));ALLOCATE(rr(npla,na))
  ALLOCATE(b2(npla,na));ALLOCATE(b3(npla,na));ALLOCATE(a_R(npla,na))
  ALLOCATE(ka(npla,na));ALLOCATE(e_ka(npla));ALLOCATE(kb(npla,na))
  ALLOCATE(dkb(npla));ALLOCATE(dkb_ini(npla))
  ALLOCATE(dur(npla,na)); ALLOCATE(sdur(npla))
  ALLOCATE(safro(npla,na))
  ALLOCATE(t0(npla,na)); ALLOCATE(st0(npla))
  ALLOCATE(per(npla,na)); ALLOCATE(sper(npla))
  ALLOCATE(tidel(npla,na)); ALLOCATE(stidel(npla))
  ALLOCATE(secosw(npla,na)); ALLOCATE(dsecosw(npla))
  ALLOCATE(sesinw(npla,na)); ALLOCATE(dsesinw(npla))
  ALLOCATE(esinw(npla,na)); ALLOCATE(ecosw(npla,na))
  ALLOCATE(dsesinw_ini(npla)); ALLOCATE(dsecosw_ini(npla))
  ALLOCATE(rephotmerit(na)); ALLOCATE(photmerit(na))
  ALLOCATE(rervmerit(na)); ALLOCATE(rvmerit(na))
  ALLOCATE(remerit(na)); ALLOCATE(merit(na))
  ALLOCATE(exc(npla,na)); ALLOCATE(omega(npla,na))
  ALLOCATE(rho(na)); ALLOCATE(rhop(npla,na))
  ALLOCATE(ave_presglo(na)); ALLOCATE(sig_presglo(na))
  ALLOCATE(ave_presbinglo(na)); ALLOCATE(sig_presbinglo(na))
  ALLOCATE(ave_presglorv(na));ALLOCATE(sig_presglorv(na))
  ALLOCATE(mulimb_nl(nmax,nc))
  ALLOCATE(mulimb_qd(nmax))
  ALLOCATE(muph(nmax))
  ALLOCATE(z2(1));ALLOCATE(rvz(nrvmax)); ALLOCATE(rossiter(nrvmax))
  ALLOCATE(er_exc(npla));ALLOCATE(er_om(npla))
  ALLOCATE(initsecosw(npla));ALLOCATE(initsesinw(npla))
  ALLOCATE(Tano_tr(npla));ALLOCATE(Eano_tr(npla))
  ALLOCATE(Tperi(npla))
  ALLOCATE(Mano_tr(npla));ALLOCATE(Mano_se(npla));ALLOCATE(ltt(npla))
  IF(ntr.GT.0)THEN
     ALLOCATE(syscor(ntr,nmax));ALLOCATE(mulimb0(ntr,nmax))
     ALLOCATE(mulimb1(ntr,nmax))
     ALLOCATE(tbjd(nphotot));ALLOCATE(tphot(nphotot))
     ALLOCATE(terror(nphotot));ALLOCATE(tresi(nphotot))
     ALLOCATE(tphotcor(nphotot));ALLOCATE(terrorcor(nphotot))
     ALLOCATE(tmodel(nphotot));ALLOCATE(tmodel_tr(nphotot))
     ALLOCATE(fbjd(nphotot)); ALLOCATE(fphot(nphotot))
     ALLOCATE(ferror(nphotot));ALLOCATE(fmodel_tr(nphotot))
     ALLOCATE(fphotcor(nphotot)); ALLOCATE(ferrorcor(nphotot))
     ALLOCATE(fresi(nphotot));ALLOCATE(fmodel(nphotot))
     ALLOCATE(tilo(nphotot));ALLOCATE(filo(nphotot))
     ALLOCATE(tflick(nphotot))
  ENDIF
  IF(nrv.GT.0)THEN
     ALLOCATE(trvbjd(nrvtot)); ALLOCATE(trv(nrvtot))
     ALLOCATE(trverror(nrvtot)); ALLOCATE(trvresi(nrvtot))
     ALLOCATE(trvmodel(nrvtot)); ALLOCATE(trvmodel2(nrvtot))
     ALLOCATE(glorvresi(nrvtot))
     ALLOCATE(glorvbjd(nrvtot));ALLOCATE(glorverror(nrvtot))
     ALLOCATE(frvbjd(nrvtot)); ALLOCATE(frv(nrvtot))
     ALLOCATE(frv2(nrvtot)); ALLOCATE(trv2(nrvtot))
     ALLOCATE(frverror(nrvtot)); ALLOCATE(frvresi(nrvtot))
     ALLOCATE(frvmodel(nrvtot)); ALLOCATE(frvmodel2(nrvtot))
     ALLOCATE(rvmodpla(npla,nrv,nrvmax));ALLOCATE(rv2pla(npla,nrv,nrvmax))
     ALLOCATE(rvbjd(nrv,nrvmax)); ALLOCATE(rv(nrv,nrvmax))
     ALLOCATE(rvfwhm(nrv,nrvmax));ALLOCATE(rvbis(nrv,nrvmax));ALLOCATE(rvcontrast(nrv,nrvmax))
     ALLOCATE(rv2(nrv,nrvmax));ALLOCATE(rverror(nrv,nrvmax)); ALLOCATE(rvresi(nrv,nrvmax))
     ALLOCATE(rvtexp(nrv,nrvmax));ALLOCATE(rvloghk(nrv,nrvmax))
     ALLOCATE(rvmod(nrv,nrvmax)); ALLOCATE(forvbjd(nrv,nrvmax)); ALLOCATE(forvbjd2(nrv,nrvmax))
     ALLOCATE(rvmod2(nrv,nrvmax));ALLOCATE(rvsysva(nrv,nrvsysparmax))
     ALLOCATE(rvsysva_bestfit(nrv,nrvsysparmax))
  ENDIF
  IF(ngroup.GT.0)THEN
    ALLOCATE(mulimb2(ntr,nmax))
  ENDIF
  IF(stelincli.EQ.'y')THEN
     ALLOCATE(incli_s(na));ALLOCATE(sinincli_s(na))
  ENDIF
  ALLOCATE(f2(na));ALLOCATE(ktide(na))

  DO i=1,3
     READ(11,*) dum
  ENDDO
  READ(11,*) dum2,starmass,dmass_s_ini,priormass
  dmass_s = dmass_s_ini
  IF(priormass.NE.'p'.AND.msrs.NE.'y'.AND.enoch.NE.'y'.AND.fitmsrs.NE.'y' &
    & .AND.fixstellar.NE.'y'.AND.massfromr.NE.'y'.and.isoch.ne.'y')THEN
    PRINT*, 'Stellar mass is an unconstrained free parameter. Do you still want to proceed?'
    READ(*,*) dum      
    IF(dum.NE.'y') STOP
  ENDIF
  IF(priormass.EQ.'p'.AND.ABS(dmass_s_ini).LT.1.E-10)THEN
    PRINT*, 'The normal prior PDF for M* has a zero width!'; STOP
  ENDIF
  if (isoch.eq.'y'.and.priormass.eq.'n') then
  	print*, 'Stellar mass must be set as jump parameter if isochrones are used'; stop
  endif
  READ(11,*) dum2,starradius,dstarradius_ini,priorrad
  dstarradius = dstarradius_ini
  IF(priorrad.NE.'p'.AND.massfromr.EQ.'y')THEN
    PRINT*, 'Stellar mass from stellar radius without prior. Do you still want to proceed?'
    READ(*,*) dum      
    IF(dum.NE.'y') STOP
  ENDIF
  IF(priorrad.EQ.'p'.AND.ABS(dstarradius_ini).LT.1.E-10)THEN
    PRINT*, 'The normal prior PDF for R* has a zero width!'; STOP
  ENDIF
  if (priorrad.eq.'p') then
	Rinp=starradius
	I_Rinp=dstarradius_ini
  else
	Rinp=-1.
	I_Rinp=-1.
  end if
  if(isoch.eq.'y'.and.ntr.eq.0.and.nrv.eq.0.and.priorrad.eq.'n')then !only perform the IsochPlacement
    print*, 'R* must be set as a jump parameter, either p or y'; stop
  end if
  READ(11,*) dum2,trho2,strho2,priorrho
  IF(priorrho.EQ.'p'.AND.ABS(strho2).LT.1.E-10)THEN
    PRINT*, 'The normal prior PDF for rho* has a zero width!'; STOP
  ENDIF
  trho = starmass/(starradius**3)
  strho = SQRT((dmass_s/(starradius**3))**2 + &
       & (3*starmass*dstarradius/(starradius**4))**2)
  IF(priormass.EQ.'p') ndata=ndata+1
  IF(priorrad.EQ.'p') ndata=ndata+1
  IF(priorrho.EQ.'p') ndata=ndata+1
  if (priorrho.ne.'n') then
  	rhoInp=trho2
  	I_rhoInp=strho2
  else
  	rhoInp=-1.
  	I_rhoInp=-1.
  end if
  READ(11,*) dum2,starlum,dstarlum,priorlum
  IF(priorlum.EQ.'p') ndata=ndata+1
  IF(priorlum.EQ.'p'.AND.ABS(dstarlum).LT.1.E-10)THEN
    PRINT*, 'The normal prior PDF for L* has a zero width!'; STOP
  ENDIF
  READ(11,*) dum2,spec(2),sspec(2),priorlogg
  IF(priorlogg.EQ.'p') ndata=ndata+1
  IF(priorlogg.EQ.'p'.AND.ABS(sspec(2)).LT.1.E-10)THEN
    PRINT*, 'The normal prior PDF for logg has a zero width!'; STOP
  ENDIF
  if (priorlogg.ne.'n') then
  	loggI=spec(2)
  	I_loggI=sspec(2)
  	if (I_loggI.gt.0.15) then !more than 35%
  		I_loggI=0.15
  		print*,'logg uncertainty too high. Set to 0.15 dex'
  	end if
  else
  	loggI=-1.
  	I_loggI=-1.
  end if
  READ(11,*) dum2,spec(1),sspec(1),priorteff
  IF(priorteff.EQ.'p'.AND.ABS(sspec(1)).LT.1.E-10)THEN
    PRINT*, 'The normal prior PDF for Teff has a zero width!'; STOP
  ENDIF
  READ(11,*) dum2,spec(3),sspec(3),priormet
  IF(priormet.EQ.'p'.AND.ABS(sspec(3)).LT.1.E-10)THEN
    PRINT*, 'The normal prior PDF for [Fe/H] a zero width!'; STOP
  ENDIF
  if (isoch.eq.'y'.and.priormet.eq.'n') then
  	print*, 'At least Teff (or color) and [Fe/H] must be set as p or y to compute mass from isochrones'
  	stop
  end if
!  FeHinf=log10(Zi_inf)+costZ
!  FeHsup=log10(Zi_sup)+costZ
  READ(11,*) dum2,spec(4),sspec(4)
  READ(11,*) dum2,spec(5),sspec(5),priorvsini
  IF(priorvsini.EQ.'p'.AND.ABS(sspec(5)).LT.1.E-10)THEN
    PRINT*, 'The normal prior PDF for VsinI* a zero width!'; STOP
  ENDIF
  IF(priorvsini.EQ.'p') ndata=ndata+1
  READ(11,*) dum2,spec(6),sspec(6),priorprot
  if (priorprot.eq.'y') then
  	protI=spec(6)
  	I_protI=sspec(6)
  else
  	protI=-1.
  	I_protI=-1.
  end if
  READ(11,*) dum2,spec(7),sspec(7),priorf2
  IF(priorf2.EQ.'p'.AND.ABS(sspec(7)).LT.1.E-10)THEN
    PRINT*, 'The normal prior PDF for F2 a zero width!'; STOP
  ENDIF
  IF(testf2.EQ.'y')THEN
     IF(priorf2.NE.'n')THEN
        njump=njump+1
        npara=npara+1
        njump_rv=njump_rv+1
     ENDIF
     IF(priorf2.EQ.'p') ndata=ndata+1
  ENDIF
  vsini(1)=spec(5);dvsini=sspec(5)
  read(11,*) dum2,logRHKI,I_logRHKI,priorlogRHK
  if (priorlogRHK.eq.'n') then
  	logRHKI=0.
  	I_logRHKI=0.
  end if
  read(11,*) dum2,YMgI,I_YMgI,priorYMg
  if (priorYMg.eq.'n') then
  	YMgI=-100.
  	I_YmgI=-100.
  end if
  read(11,*) dum2,numaxI,I_numaxI,priornumax
  if (priornumax.eq.'n') then
  	numaxI=-1.
  	I_numaxI=-1.
  end if
  read(11,*) dum2,deltanuI,I_deltanuI,priordeltanu
  if (priordeltanu.eq.'n') then
  	deltanuI=-1.
  	I_deltanuI=-1.
  end if
  if (priornumax.eq.'n'.and.priordeltanu.eq.'n') then
  	sismo=0
  else
  	sismo=1
  end if
  read(11,*) dum2,mag,I_mag,priormag
  if (priormag.eq.'n') then
    mag=-100.
    I_mag=-100.
  end if
  read(11,*) dum2,col,I_col,idCol,priorcol
  if (idCol.lt.1.or.idCol.gt.2) then
    print*,'Photometry different from (V,B-V) and (G,G_BP-G_RP) not yet implemented'; stop
  end if
  if (priorcol.eq.'n') then
    col=-100.
    I_col=-100.
    if (priorteff.ne.'p') then
      print*, 'Specify either Teff or color'; stop
    end if
  else
    if (isoch.ne.'y') then
      print*,'Isochrones are needed to calibrate Teff from color'; stop
    end if
    if (I_col.lt.1.e-10) then
      print*,'Specify uncertainty for color'; stop
    end if
  end if
  if (priorteff.eq.'p'.and.priorcol.eq.'p') then
    print*, 'Specify just one prior among Teff and color'; stop
  end if
  read(11,*) dum2,dist,I_dist,priordist
  if (.not.((priormag.eq.'n'.and.priordist.eq.'n').or. &
  	& (priormag.eq.'y'.and.priordist.eq.'y').or. &
  	& (priormag.eq.'p'.and.priordist.eq.'p'))) then
  	print*,'Magnitude and distance must be set both to n, y or p'; stop
  end if
  if (priordist.eq.'n') then
    dist=-1.
    I_dist=-1.
    if (priorcol.eq.'p') then 
      mag=0.
      I_mag=0.
    end if
    Hipf=-1.
  else
    Hipf=0.
  end if
  if (priormag.ne.'n'.and.priordist.ne.'n'.and.priorcol.eq.'n') then
    useColor=0.
  else
    useColor=1.
  end if
  if (priormag.eq.'p'.and.priordist.eq.'p'.and.priorlum.eq.'p') then
  	print*,'Error: both (mag,dist) and luminosity are set as priors'; stop
  end if
  if (priormag.eq.'p'.and.priordist.eq.'p') then
  	ndata=ndata+1
  end if
  if (isoch.eq.'y') then
  	if (priorlogg.ne.'n'.or.(priormag.ne.'n'.and.priordist.ne.'n') &
		& .or.priorlum.ne.'n'.or.sismo.eq.1) then
	  ybesR=.true.
  	else
	  ybesR=.false.
	endif
  endif
  DO i=1,3
     READ(11,*) dum
  ENDDO

  d=0
  dN=0
  DO i=1,npla
     READ(11,*) dum2,dF(i,1),sdF_ini(i),isjump(i,1)
     dF_ini(i)=dF(i,1)
     IF(nddf.GT.0.AND.isddf.NE.'n')THEN
        DO j=1,nddf
           ddf(i,j,1)=0.
           eddf_ini(i,j)=sigddf
           eddf(i,j)=sigddf
        ENDDO
     ENDIF
     IF(isjump(i,1).NE.'n')THEN
        njump=njump+1
        npara=npara+1
        IF(gelman.EQ.'y'.AND.nchain.GT.1.AND.test_dF.EQ.0)THEN 
          ALLOCATE(gelman_dF(npla));test_dF=1
        ENDIF
     ENDIF
     IF(isjump(i,1).EQ.'p') ndata=ndata+1
     READ(11,*) dum2,tidel(i,1),stidel_ini(i),isjump(i,9)
     IF(tidel(i,1).LT.1) tidel(i,1)=1.
     IF(tidel(i,1).GT.tidelmax) tidel(i,1)=tidelmax
     tidel_ini(i)=tidel(i,1)
     IF(isjump(i,9).NE.'n')THEN
        njump=njump+1
        npara=npara+1
        IF(gelman.EQ.'y'.AND.nchain.GT.1.AND.test_tidel.EQ.0)THEN 
          ALLOCATE(gelman_tidel(npla));test_tidel=1
        ENDIF
     ENDIF
     IF(isjump(i,9).EQ.'p') ndata=ndata+1
     READ(11,*) dum2,b(i,1),sb_ini(i), isjump(i,2)
     b_ini(i)=b(i,1)
     IF(isjump(i,2).NE.'n')THEN
        njump=njump+1
        npara=npara+1
        IF(gelman.EQ.'y'.AND.nchain.GT.1.AND.test_b.EQ.0)THEN
          ALLOCATE(gelman_b(npla));test_b=1
        ENDIF
     ENDIF
     IF(isjump(i,2).EQ.'p') ndata=ndata+1
     READ(11,*) dum2,dur(i,1),sdur_ini(i), isjump(i,3)
     dur_ini(i)=dur(i,1)
     IF(fitmsrs.EQ.'y'.OR.fixstellar.EQ.'y')THEN
       dur_ini(i)=0.
       sdur_ini(i)=0.
       isjump(i,3)='n'
     ENDIF
     if(isoch.eq.'y'.and.ntr.eq.0.and.nrv.eq.0)then !only perform the IsochPlacement
       dur(i,1)=0.
       sdur_ini(i)=0.
	   isjump(i,3)='n'
     end if
     if (fitmsrs.eq.'n'.and.fixstellar.eq.'n') then
	  	if (isEq(dur(i,1),0.D0,5)) then
	  	  if (isoch.eq.'y') then
	  	      !just to be sure everything=0
	  	      dur_ini(i)=0.
		      sdur_ini(i)=0.
		      isjump(i,3)='n'
	  	  end if
	  	else
	  		d=1
	  		dN=dN+1
	  	end if
	 end if
     IF(isjump(i,3).NE.'n')THEN
        njump=njump+1
        npara=npara+1      
        IF(gelman.EQ.'y'.AND.nchain.GT.1.AND.test_dur.EQ.0)THEN
          ALLOCATE(gelman_dur(npla));test_dur=1
        ENDIF
     ENDIF
     IF(isjump(i,3).EQ.'p') ndata=ndata+1
     READ(11,*) dum2,t0(i,1),st0_ini(i),isjump(i,4)
     IF(isjump(i,4).NE.'n')THEN
        njump=njump+1; njump_rv=njump_rv+1; npara=npara+1
        IF(gelman.EQ.'y'.AND.nchain.GT.1.AND.test_t0.EQ.0)THEN
          ALLOCATE(gelman_t0(npla));test_t0=1
        ENDIF
     ENDIF
     IF(isjump(i,4).EQ.'p') ndata=ndata+1
     t0(i,1) = t0(i,1) - 0.
     t0_ini(i) = t0(i,1)
     READ(11,*) dum2,per(i,1),sper_ini(i),isjump(i,5)
     if(isoch.eq.'y'.and.ntr.eq.0.and.nrv.eq.0)then !only perform the IsochPlacement
     	per(i,1)=10. !just a fake
     end if
     IF(isjump(i,5).NE.'n')THEN
        njump=njump+1; njump_rv=njump_rv+1; npara=npara+1
        IF(gelman.EQ.'y'.AND.nchain.GT.1.AND.test_per.EQ.0)THEN
          ALLOCATE(gelman_per(npla));test_per=1
        ENDIF
     ENDIF
     IF(isjump(i,5).EQ.'p') ndata=ndata+1
     per_ini(i) = per(i,1)
     READ(11,*) dum2,exc(i,1),er_exc(i),isjump(i,6)  
     IF(isjump(i,6).NE.'n')THEN
        njump=njump+2; njump_rv=njump_rv+2; npara=npara+2
        IF(gelman.EQ.'y'.AND.nchain.GT.1.AND.tgus.EQ.1)THEN 
           ALLOCATE(gelman_secosw(npla));ALLOCATE(gelman_sesinw(npla))
           tgus=2
        ENDIF
     ENDIF
     IF(isjump(i,6).EQ.'p') ndata=ndata+2
     READ(11,*) dum2,omega(i,1),er_om(i)
     secosw(i,1) = SQRT(exc(i,1))*DCOS(omega(i,1)*pi/180)
     sesinw(i,1) = SQRT(exc(i,1))*DSIN(omega(i,1)*pi/180)
     secosw_ini(i)=secosw(i,1);sesinw_ini(i)=sesinw(i,1)
     READ(11,*) dum2, initsecosw(i), initsesinw(i)
     IF(initsecosw(i).LT.1.E-6.AND.initsesinw(i).LT.1.E-6)THEN
        dsecosw_ini(i) = SQRT(er_exc(i))
        dsesinw_ini(i) = SQRT(er_exc(i))
     ELSE
        dsecosw_ini(i) = initsecosw(i)
        dsesinw_ini(i) = initsesinw(i)
     ENDIF
     READ(11,*) dum2,ka(i,1),e_ka(i),isjump(i,7)
     IF(isjump(i,7).NE.'n')THEN
        njump=njump+1; njump_rv=njump_rv+1+nrv; npara=npara+1
        IF(gelman.EQ.'y'.AND.nchain.GT.1.AND.test_kb.EQ.0)THEN
          ALLOCATE(gelman_kb(npla));test_kb=1
        ENDIF
     ENDIF
     IF(isjump(i,7).EQ.'p') ndata=ndata+1
     kb(i,1) = ka(i,1)*per(i,1)**(1./3.)
     kb(i,1) = kb(i,1)*SQRT(1.-exc(i,1)**2) 
     kb_ini(i)=kb(i,1)
     dkb_ini(i) = e_ka(i)
          ! SQRT((e_ka(i)*SQRT(1-exc(i,1)**2)*per(i,1)**(1./3.))**2 + &
          ! & ((ka(i,1)/(3.*per(i,1)**(2./3.)))*SQRT(1.-exc(i,1)**2)*sper_ini(i))**2 + &
          ! & ((ka(i,1)*per(i,1)**(1./3.))*exc(i,1)*er_exc(i)/SQRT(1.-exc(i,1)**2))**2) 
     IF(i.EQ.1)THEN
        READ(11,*) dum2,beta(1),dbeta,isjump(i,8)
        IF(isjump(i,8).NE.'n')THEN
           njump=njump+2; njump_rv=njump_rv+2; npara=npara+2
        ENDIF
        IF(isjump(i,8).EQ.'p') ndata=ndata+2
        READ(11,*) dum2,initsvsinicosb, initsvsinisinb
        svsinisinbeta(1) = SQRT(vsini(1))*DSIN(beta(1)*pi/180.)
        svsinicosbeta(1) = SQRT(vsini(1))*DCOS(beta(1)*pi/180.)
        ini_svsinisinbeta =   svsinisinbeta(1)
        ini_svsinicosbeta =   svsinicosbeta(1)
        IF(initsvsinicosb.LT.1.E-10.AND.initsvsinisinb.LT.1.E-10)THEN
           dsvsinisinbeta_ini=(0.5*dvsini*DSIN(beta(1)*pi/180.)/SQRT(vsini(1)))**2+ & 
           & (SQRT(vsini(1))*DCOS(beta(1)*pi/180.)*dbeta*pi/180)**2
           dsvsinisinbeta_ini=SQRT(dsvsinisinbeta_ini)
           dsvsinicosbeta_ini=(0.5*dvsini*DCOS(beta(1)*pi/180.)/SQRT(vsini(1)))**2+ &
           & (SQRT(vsini(1))*DSIN(beta(1)*pi/180.)*dbeta*pi/180)**2
           dsvsinicosbeta_ini=SQRT(dsvsinicosbeta_ini)
        ELSE
          dsvsinisinbeta_ini=initsvsinisinb
          dsvsinicosbeta_ini=initsvsinicosb
        ENDIF
     ENDIF
  ENDDO
  if (priorrad.eq.'y'.and..not.ybesR) then
  	Rinp=starradius
  	I_Rinp=dstarradius_ini
  endif
  if (isoch.eq.'y'.and.(d.eq.0)) then !.or.(d.eq.1.and.priorrad.eq.'p'.and.nrv.eq.0)
  	njump=njump+1
  	npara=npara+1
  	RjumpIso='y'
  	if (priorrad.eq.'n') then
  		print*, 'Set radius in the input form'; stop
  	endif
  endif
  DO i=1,4
     READ(11,*) dum
  ENDDO

  !******************************************************************************
  !4. LIMB-DARKENING                                                            !
  !******************************************************************************

  IF(nrv.GT.0)THEN
     tefil = 'JV'
     CALL qdinter(spec,sspec,tefil,qtemp,sqtemp)
     rold(1,1) = qtemp(1)
     rold(1,2) = qtemp(2)
     srold(1) = sqtemp(1)
     srold(2) = sqtemp(2)
     PRINT*, '  Interpolated quadratic LD coefficients for Rossiter:'
     WRITE(*,100) '   u1= ',rold(1,1),' +- ',srold(1)
     WRITE(*,100) '   u2= ',rold(1,2),' +- ',srold(2)
     PRINT*, ' '
     rold_ini(1)=rold(1,1)
     rold_ini(2)=rold(1,2)
  ENDIF
  IF(ntr.GT.0)THEN
     ce=0
     SELECT CASE(limb)
     CASE('qd')
        PRINT*, 'Limb-darkening model: quadratic'
        DO i=1,nfi
           loc=0
           temp = wfilter(i) // '-filter'
           tefil = wfilter(i)
           WRITE(*,'(I3,1x,A10)') i,temp
           IF(wfilter(i).NE."Su".AND.wfilter(i).NE."Sv".AND. &
                &   wfilter(i).NE."Sb".AND.wfilter(i).NE."Sy".AND. &
                &   wfilter(i).NE."JU".AND.wfilter(i).NE."JB".AND.wfilter(i).NE."JV".AND. &
                &   wfilter(i).NE."JR".AND.wfilter(i).NE."JI".AND.wfilter(i).NE."JJ".AND. &
                &   wfilter(i).NE."JH".AND.wfilter(i).NE."JK".AND.wfilter(i).NE."u'".AND. &
                &   wfilter(i).NE."g'".AND.wfilter(i).NE."r'".AND.wfilter(i).NE."i'".AND. &
                &   wfilter(i).NE."z'".AND.wfilter(i).NE."Ke".AND.wfilter(i).NE."Co".AND. &
                &   wfilter(i).NE."S1".AND.wfilter(i).NE."S2".AND.wfilter(i).NE."S3".AND. &
                &   wfilter(i).NE."S4".AND.wfilter(i).NE."Ch".AND.wfilter(i).NE."TE")THEN
              temp2 = '  ' // wfilter(i) // '-filter is not in tables.'
              PRINT*, temp2
              READ(11,*) dum,tefil2,dFsec(i,1,1),edFsec_ini(i,1),dilu(i),edilu(i),phampli1_ini(i), &
               & ephampli1_ini(i),phampli2_ini(i),ephampli2_ini(i),phampli3_ini(i),ephampli3_ini(i), &
               & phoffset_ini(i),ephoffset_ini(i),iffitoc(i),& 
               & fitlimb(i),iffitph(i),ql(i,1,1),sql(i,1), ql(i,2,1),sql(i,2)
              IF(iffitoc(i).NE.'n')THEN
                njump=njump+npla
                npara=npara+npla
                dFsec(i,1,:)=dFsec(i,1,1)
                edFsec_ini(i,:)=edFsec_ini(i,1)   
              ENDIF
              IF(iffitph(i).NE.'n')THEN
                IF(ephampli1_ini(i).GT.1E-10)THEN
                  njump=njump+1
                  npara=npara+1
                  IF(iffitph(i).EQ.'p') ndata=ndata+1
                ENDIF 
                IF(ephampli2_ini(i).GT.1E-10)THEN
                  njump=njump+1
                  npara=npara+1
                  IF(iffitph(i).EQ.'p') ndata=ndata+1
                ENDIF  
                IF(ephampli3_ini(i).GT.1E-10)THEN
                  njump=njump+1
                  npara=npara+1
                  IF(iffitph(i).EQ.'p') ndata=ndata+1
                ENDIF 
                IF(ephoffset_ini(i).GT.1E-10)THEN
                  njump=njump+1
                  npara=npara+1
                  IF(iffitph(i).EQ.'p') ndata=ndata+1
                ENDIF   
              ENDIF
              IF(iffitoc(i).EQ.'p') ndata=ndata+npla
              dFsec_ini(i,:)=dFsec(i,1,:)
              phampli1(i,1)=phampli1_ini(i)
              phampli2(i,1)=phampli2_ini(i)
              phampli3(i,1)=phampli3_ini(i)
              phoffset(i,1)=phoffset_ini(i)
              IF(tefil.NE.tefil2)THEN
                PRINT*, 'Error in mcmc.dat: ',tefil,' /= ',tefil2
                STOP
              ENDIF
              PRINT*, '  Quadratic LD coefficients read in mcmc.dat:'
              WRITE(*,100) '   u1= ',ql(i,1,1),' +- ',sql(i,1)
              WRITE(*,100) '   u2= ',ql(i,2,1),' +- ',sql(i,2)
              loc=1
           ENDIF
           IF(loc.NE.1)THEN
              CALL qdinter(spec,sspec,tefil,qtemp,sqtemp)
              ql(i,1,1) = qtemp(1)
              ql(i,2,1) = qtemp(2)
              sql(i,1) = sqtemp(1)
              sql(i,2) = sqtemp(2)
              PRINT*, '  Interpolated quadratic LD coefficients:'
              WRITE(*,100) '   u1= ',ql(i,1,1),' +- ',sql(i,1)
              WRITE(*,100) '   u2= ',ql(i,2,1),' +- ',sql(i,2)
              READ(11,*) dum,tefil2,dFsec(i,1,1),edFsec_ini(i,1),dilu(i),edilu(i),phampli1_ini(i), &
               & ephampli1_ini(i),phampli2_ini(i),ephampli2_ini(i),phampli3_ini(i),ephampli3_ini(i), &
               & phoffset_ini(i),ephoffset_ini(i),iffitoc(i),&
               & fitlimb(i),iffitph(i)
!              print*,'val letti',dum,tefil2,dFsec(i,1,1)
              IF(iffitoc(i).NE.'n')THEN
                njump=njump+1
                npara=npara+1      
              ENDIF
              IF(iffitph(i).NE.'n')THEN
                IF(ephampli1_ini(i).GT.1E-10)THEN
                  njump=njump+1
                  npara=npara+1
                  IF(iffitph(i).EQ.'p') ndata=ndata+1
                ENDIF 
                IF(ephampli2_ini(i).GT.1E-10)THEN
                  njump=njump+1
                  npara=npara+1
                  IF(iffitph(i).EQ.'p') ndata=ndata+1
                ENDIF  
                IF(ephampli3_ini(i).GT.1E-10)THEN
                  njump=njump+1
                  npara=npara+1
                  IF(iffitph(i).EQ.'p') ndata=ndata+1
                ENDIF 
                IF(ephoffset_ini(i).GT.1E-10)THEN
                  njump=njump+1
                  npara=npara+1
                  IF(iffitph(i).EQ.'p') ndata=ndata+1
                ENDIF   
              ENDIF
              IF(iffitoc(i).EQ.'p') ndata=ndata+npla
              dFsec_ini(i,:)=dFsec(i,1,:)
              phampli1(i,1)=phampli1_ini(i)
              phampli2(i,1)=phampli2_ini(i)
              phampli3(i,1)=phampli3_ini(i)
              phoffset(i,1)=phoffset_ini(i)
              IF(tefil.NE.tefil2)THEN
                PRINT*, 'Error in mcmc.dat: ',tefil,' /= ',tefil2
                STOP
              ENDIF
           ENDIF
           IF(fitlimb(i).NE.'n')THEN
              njump = njump+2; npara=npara+2
              IF(gelman.EQ.'y'.AND.nchain.GT.1.AND.ce.LT.1)THEN
                 ALLOCATE(gelman_jumplimb(nfi,2))
                 ce=2
              ENDIF
              jumplimb(i,1,1) = 2*ql(i,1,1) + ql(i,2,1)
              jumplimb(i,2,1) = ql(i,1,1) - 2*ql(i,2,1)
              djumplimb_ini(i,1) = SQRT(4*sql(i,1)**2 + sql(i,2)**2)
              djumplimb_ini(i,2) = SQRT(sql(i,1)**2 + 4*sql(i,2)**2)
              djumplimb(i,1)=djumplimb_ini(i,1)
              djumplimb(i,2)=djumplimb_ini(i,2)
           ENDIF
           IF(fitlimb(i).EQ.'p') ndata=ndata+2
           ql_ini(i,1)=ql(i,1,1)
           ql_ini(i,2)=ql(i,2,1)
           PRINT*, 
        ENDDO
     CASE('no')
        PRINT*, 'Limb-darkening model: none'
        DO i=1,nfi
           loc=0
           temp = wfilter(i) // '-filter'
           tefil = wfilter(i)
           READ(11,*) dum,tefil2,dFsec(i,1,1),edFsec_ini(i,1),dilu(i),edilu(i),phampli1_ini(i), &
                & ephampli1_ini(i),phampli2_ini(i),ephampli2_ini(i),phampli3_ini(i),ephampli3_ini(i), &
                & phoffset_ini(i),ephoffset_ini(i),iffitoc(i),&
                & fitlimb(i),iffitph(i)
            IF(iffitoc(i).NE.'n')THEN
              njump=njump+1
              npara=npara+1      
           ENDIF
           IF(iffitph(i).NE.'n')THEN
             IF(ephampli1_ini(i).GT.1E-10)THEN
               njump=njump+1
               npara=npara+1
               IF(iffitph(i).EQ.'p') ndata=ndata+1
             ENDIF 
             IF(ephampli2_ini(i).GT.1E-10)THEN
               njump=njump+1
               npara=npara+1
               IF(iffitph(i).EQ.'p') ndata=ndata+1
             ENDIF  
             IF(ephampli3_ini(i).GT.1E-10)THEN
               njump=njump+1
               npara=npara+1
               IF(iffitph(i).EQ.'p') ndata=ndata+1
             ENDIF 
             IF(ephoffset_ini(i).GT.1E-10)THEN
               njump=njump+1
               npara=npara+1
               IF(iffitph(i).EQ.'p') ndata=ndata+1
             ENDIF   
           ENDIF
           IF(iffitoc(i).EQ.'p') ndata=ndata+npla
           dFsec_ini(i,:)=dFsec(i,1,:)
           phampli1(i,1)=phampli1_ini(i)
           phampli2(i,1)=phampli2_ini(i)
           phampli3(i,1)=phampli3_ini(i)
           phoffset(i,1)=phoffset_ini(i)
           IF(tefil.NE.tefil2)THEN
             PRINT*, 'Error in mcmc.dat: ',tefil,' /= ',tefil2
             STOP
           ENDIF
           ql(i,1,1) = 0.
           ql(i,2,1) = 0.
           sql(i,1) = 0.
           sql(i,2) = 0.
        ENDDO
     CASE('nl')
        PRINT*, 'Limb-darkening model: non-linear'
        DO i=1,nfi
           loc=0
           temp = wfilter(i) // '-filter'
           tefil = wfilter(i)
           PRINT*, temp
           IF(wfilter(i).NE."Su".AND.wfilter(i).NE."Sv".AND. &
                &   wfilter(i).NE."Sb".AND.wfilter(i).NE."Sy".AND. &
                &   wfilter(i).NE."JU".AND.wfilter(i).NE."JB".AND.wfilter(i).NE."JV".AND. &
                &   wfilter(i).NE."JR".AND.wfilter(i).NE."JI".AND.wfilter(i).NE."JJ".AND. &
                &   wfilter(i).NE."JH".AND.wfilter(i).NE."JK".AND.wfilter(i).NE."u'".AND. &
                &   wfilter(i).NE."g'".AND.wfilter(i).NE."r'".AND.wfilter(i).NE."i'".AND. &
                &   wfilter(i).NE."z'".AND.wfilter(i).NE."Ke".AND.wfilter(i).NE."Co".AND. &
                &   wfilter(i).NE."S1".AND.wfilter(i).NE."S2".AND.wfilter(i).NE."S3".AND. &
                &   wfilter(i).NE."S4".AND.wfilter(i).NE."Ch".AND.wfilter(i).NE."TE")THEN
              temp2 = wfilter(i) // '-filter is not in tables.'
              PRINT*, temp2
              READ(11,*) dum,tefil2,dFsec(i,1,1),edFsec_ini(i,1),dilu(i),edilu(i),phampli1_ini(i), &
                   & ephampli1_ini(i),phampli2_ini(i),ephampli2_ini(i),phampli3_ini(i),ephampli3_ini(i), &
                   & phoffset_ini(i),ephoffset_ini(i),iffitoc(i),&
                   & fitlimb(i),iffitph(i),nl(i,1,1),snl(i,1), nl(i,2,1),snl(i,2),nl(i,3,1),snl(i,3), nl(i,4,1),snl(i,4)
              IF(iffitoc(i).NE.'n')THEN
                njump=njump+1
                npara=npara+1      
              ENDIF
              IF(iffitph(i).NE.'n')THEN
                IF(ephampli1_ini(i).GT.1E-10)THEN
                  njump=njump+1
                  npara=npara+1
                  IF(iffitph(i).EQ.'p') ndata=ndata+1
                ENDIF 
                IF(ephampli2_ini(i).GT.1E-10)THEN
                  njump=njump+1
                  npara=npara+1
                  IF(iffitph(i).EQ.'p') ndata=ndata+1
                ENDIF  
                IF(ephampli3_ini(i).GT.1E-10)THEN
                  njump=njump+1
                  npara=npara+1
                  IF(iffitph(i).EQ.'p') ndata=ndata+1
                ENDIF 
                IF(ephoffset_ini(i).GT.1E-10)THEN
                  njump=njump+1
                  npara=npara+1
                  IF(iffitph(i).EQ.'p') ndata=ndata+1
                ENDIF   
              ENDIF
              IF(iffitoc(i).EQ.'p') ndata=ndata+npla
              dFsec_ini(i,:)=dFsec(i,1,:)
              phampli1(i,1)=phampli1_ini(i)
              phampli2(i,1)=phampli2_ini(i)
              phampli3(i,1)=phampli3_ini(i)
              phoffset(i,1)=phoffset_ini(i)
              IF(tefil.NE.tefil2)THEN
                PRINT*, 'Error in mcmc.dat: ',tefil,' /= ',tefil2
                STOP
              ENDIF
              PRINT*, '  Non-linear LD coefficients read in mcmc.dat:'
              WRITE(*,100) '   a1= ',nl(i,1,1),' +- ',snl(i,1)
              WRITE(*,100) '   a2= ',nl(i,2,1),' +- ',snl(i,2)
              WRITE(*,100) '   a3= ',nl(i,3,1),' +- ',snl(i,3)
              WRITE(*,100) '   a3= ',nl(i,4,1),' +- ',snl(i,4)
              loc=1
           ENDIF
           IF(loc.NE.1)THEN
              CALL nlinter(spec,sspec,tefil,ntemp,sntemp)
              nl(i,1,1) = ntemp(1)
              nl(i,2,1) = ntemp(2)
              nl(i,3,1) = ntemp(3)
              nl(i,4,1) = ntemp(4)
              snl(i,1) = sntemp(1)
              snl(i,2) = sntemp(2)
              snl(i,3) = sntemp(3)
              snl(i,4) = sntemp(4)
              PRINT*, '  Interpolated non-linear LD coefficients:'
              WRITE(*,100) '  a1= ',nl(i,1,1),' +- ',snl(i,1)
              WRITE(*,100) '  a2= ',nl(i,2,1),' +- ',snl(i,2)
              WRITE(*,100) '  a3= ',nl(i,3,1),' +- ',snl(i,3)
              WRITE(*,100) '  a4= ',nl(i,4,1),' +- ',snl(i,4)
              READ(11,*) dum,tefil2,dFsec(i,1,1),edFsec_ini(i,1),dilu(i),edilu(i),phampli1_ini(i), &
                   & ephampli1_ini(i),phampli2_ini(i),ephampli2_ini(i),phampli3_ini(i),ephampli3_ini(i), &
                   & phoffset_ini(i),ephoffset_ini(i),iffitoc(i),fitlimb(i),iffitph(i)
              IF(iffitoc(i).NE.'n')THEN
                njump=njump+1
                npara=npara+1      
              ENDIF
              IF(iffitph(i).NE.'n')THEN
                IF(ephampli1_ini(i).GT.1E-10)THEN
                  njump=njump+1
                  npara=npara+1
                  IF(iffitph(i).EQ.'p') ndata=ndata+1
                ENDIF 
                IF(ephampli2_ini(i).GT.1E-10)THEN
                  njump=njump+1
                  npara=npara+1
                  IF(iffitph(i).EQ.'p') ndata=ndata+1
                ENDIF  
                IF(ephampli3_ini(i).GT.1E-10)THEN
                  njump=njump+1
                  npara=npara+1
                  IF(iffitph(i).EQ.'p') ndata=ndata+1
                ENDIF 
                IF(ephoffset_ini(i).GT.1E-10)THEN
                  njump=njump+1
                  npara=npara+1
                  IF(iffitph(i).EQ.'p') ndata=ndata+1
                ENDIF   
              ENDIF
              IF(iffitoc(i).EQ.'p') ndata=ndata+npla
              dFsec_ini(i,:)=dFsec(i,1,:)
              phampli1(i,1)=phampli1_ini(i)
              phampli2(i,1)=phampli2_ini(i)
              phampli3(i,1)=phampli3_ini(i)
              phoffset(i,1)=phoffset_ini(i)
              IF(tefil.NE.tefil2)THEN
                PRINT*, 'Error in mcmc.dat: ',tefil,' /= ',tefil2
                STOP 
              ENDIF
           ENDIF
           nl_ini(i,1)=nl(i,1,1)
           nl_ini(i,2)=nl(i,2,1)
           nl_ini(i,3)=nl(i,3,1)
           nl_ini(i,4)=nl(i,4,1)
        ENDDO
     END SELECT
  ENDIF
  IF(gibbs_sampler.EQ.'y')THEN
    WRITE(*,'(A21,I6)') 'N jump parameters =  ',njump
    WRITE(*,'(A21,I6)') 'N GS cycles/param =  ',nburn/(statlen2*njump)
    PRINT*, ' '
  ENDIF

  IF(ngroup.GT.0)THEN
    DO i=1,3
      READ(11,*) dum
    ENDDO
    DO i=1,ngroup
      READ(11,*) dum,dfgroup_ini(i),edfgroup_ini(i)
    ENDDO
  ENDIF
  CLOSE(11)

  !******************************************************************************
  !5. OPENING OF OUTPUT FILES                                                   !
  !******************************************************************************

  OPEN(UNIT=51,FILE='mcmc_chi2.res')
  OPEN(UNIT=444,FILE='mcmc_bf.res')
  OPEN(UNIT=445,FILE='mcmc_med.res')
  IF(ntr.GT.0)THEN
     OPEN(UNIT=446,FILE='mcm_sysphot.res')
  ENDIF

  !******************************************************************************
  !6. INPUT FILE PARAMETERS WRITING                                             !
  !******************************************************************************
 
  IF(gibbs_sampler.EQ.'y')THEN
     nburn2 = nburn*0.75
     nburn2 = nburn2 - MOD(nburn2-1,njump)
  ENDIF
  IF(nburn.GE.nat)THEN
     PRINT*, 'Error in mcmc.dat/nburn too large'; STOP
  ENDIF
  k=0
  DO i=1,nchain
     DO j=1,nat
        k=k+1
        IF(j.LE.nburn)THEN
           burned(k)='y'
        ELSE
           burned(k)='n'
        ENDIF
     ENDDO
  ENDDO
  
  !******************************************************************************
  !7. READING OF THE PHOTOMETRIC TIME-SERIES                                    !
  !   FORMAT: BJD, PHOT, ERROR, X, Y, AIRMASS                                   !
  !******************************************************************************

  IF(ntr.GT.0)THEN
     DO i=1,ntr         
        let1 = CHAR(48+int(i/1000))
        i2 = i - INT(i/1000)*1000
        let2 = CHAR(48+int(i2/100))
        i2 = i - INT(i/100)*100
        let3 = CHAR(48+int(i2/10))
        let4 = CHAR(48+mod(i2,10))
        fin =  let1 // let2 // let3 // let4 // '.txt'
        name = 'phot' // fin
        OPEN(UNIT=3,FILE=name)
        nb=np(i)
        meX = 0.
        meY = 0.
        DO j=1,nb
          IF(pldcor(i).GT.0)THEN
            READ(3,*) bjd(i,j),phot(i,j),error(i,j),dX(i,j),dY(i,j),fwhm(i,j),fwhmx(i,j),fwhmy(i,j), &
               & sky(i,j),airmass(i,j),texp(i,j),pixcon(i,j,1:9)
          ! fwhmx(i,j)=DCOS(fwhmx(i,j)*pi/180.)
          ELSE
            IF(oldform(i).EQ.'y')THEN
              READ(3,*) bjd(i,j),phot(i,j),error(i,j),dX(i,j),dY(i,j),airmass(i,j),fwhm(i,j), &
               & sky(i,j),texp(i,j)
              fwhmy(i,j)=fwhm(i,j)
              fwhmx(i,j)=fwhm(i,j)
            ELSE
              READ(3,*) bjd(i,j),phot(i,j),error(i,j),dX(i,j),dY(i,j),fwhm(i,j),fwhmx(i,j),fwhmy(i,j), &
               & sky(i,j),airmass(i,j),texp(i,j)
            ENDIF
          ENDIF
           texp(i,j)=texp(i,j)/(3600*24.)              ! assumes Texp in s
           bjd(i,j) = bjd(i,j) - 0.                    ! assumes date = BJD-2450000.
           IF(utc_tdb(i).EQ.'y')THEN
              dum4 = bjd(i,j)
              bjd(i,j) = tdbtra(dum4)
              ! print*, bjd(i,j)
           ENDIF
           bjd(i,j) = bjd(i,j) + cortime(i)
           meX = meX + dX(i,j)/DBLE(nb)
           meY = meY + dY(i,j)/DBLE(nb)
        ENDDO
        IF(sinusnumber(i).GT.0)THEN
           k = sinusnumber(i)
           DO j=1,k
             sit0_ini(i,j)=bjd(i,1)
             esit0_ini(i,j)=sip_ini(i,j)/10.
           ENDDO
        ENDIF
        tie = bjd(i,1) - red_dur(1)
        DO j=1,nb
           bjd2(i,j) = bjd(i,j) - tie
        ENDDO
        meX2 = DBLE(NINT(meX))
        meY2 = DBLE(NINT(meY))
        IF(meX.GE.meX2)THEN
           meX2 = meX2 + 0.5
        ELSE
           meX2 = meX2 - 0.5
        ENDIF
        IF(meY.GE.meY2)THEN
           meY2 = meY2 + 0.5
        ELSE
           meY2 = meY2 - 0.5
        ENDIF
        CLOSE(3)
        DO j=1,nb
           dX(i,j) = dX(i,j) - meX2
           dY(i,j) = dY(i,j) - meY2
        ENDDO
        name = 'phoc' // fin
        OPEN(UNIT=3,FILE=name)
        DO j=1,nb
           WRITE(3,*) bjd(i,j),phot(i,j),error(i,j),dX(i,j),dY(i,j), &
            & fwhm(i,j),fwhmx(i,j),fwhmy(i,j), &
            & sky(i,j),airmass(i,j),texp(i,j)*3600.*24.
        ENDDO
        CLOSE(3)
     ENDDO
     DO i=1,ntr
        nb = np(i)
        DO j=1,ne3
           nbin(i,j) = NINT((bjd(i,nb)-bjd(i,1))/red_dur(j))
           IF(nbin(i,j).LT.3) nbin(i,j)=3
           IF(nbin(i,j).GE.nb)THEN
              nbin(i,j)=nb
              np_bin(i,j)=1
              nlast(i,j)=1
           ELSE
              np_bin(i,j) = NINT(DBLE(nb)/DBLE(nbin(i,j)))
              test=0
              DO WHILE(test.LT.1)
                 nlast(i,j) = nb-(nbin(i,j)-1)*np_bin(i,j)
                 IF(nlast(i,j).LT.np_bin(i,j))THEN
                    nbin(i,j)=nbin(i,j)-1
                 ELSE IF(nlast(i,j).GT.2*np_bin(i,j))THEN
                    nbin(i,j)=nbin(i,j)+1
                 ELSE
                   test=2
                 ENDIF
              ENDDO
           ENDIF
           IF(j.EQ.1.AND.nbin(i,j).GT.ncra) ncra = nbin(i,j)
        ENDDO
      ENDDO
      DO i=1,ntr
         DO j=1,ne3
           IF(j.EQ.1)THEN
              nphototbin = nphototbin + nbin(i,j)
              nmeanbin = nmeanbin + np_bin(i,j)/DBLE(ntr)
           ENDIF
           npave(i,j) = np_bin(i,j)        
           IF(npave(i,j).LT.1) npave(i,j) = nlast(i,j)
        ENDDO
     ENDDO
     ALLOCATE(tresibin(nphototbin)); ALLOCATE(bjdbin(ntr,ncra))
     ALLOCATE(resibin(ntr,ncra,ne3)); ALLOCATE(eresibin(ntr,ncra))
     ALLOCATE(photbin(ntr,ncra)); ALLOCATE(ephotbin(ntr,ncra))
     ALLOCATE(modelbin(ntr,ncra)); ALLOCATE(fobjdbin(ntr,ncra))
     ALLOCATE(photcorbin(ntr,ncra)); ALLOCATE(ephotcorbin(ntr,ncra))
     ALLOCATE(photcor2bin(ntr,ncra)); ALLOCATE(ephotcor2bin(ntr,ncra))
     ALLOCATE(model_trbin(ntr,ncra));ALLOCATE(dXbin(ntr,ncra))
     ALLOCATE(fwhmbin(ntr,ncra));ALLOCATE(fwhmybin(ntr,ncra)) 
     ALLOCATE(fwhmxbin(ntr,ncra));ALLOCATE(airmassbin(ntr,ncra))
     ALLOCATE(dYbin(ntr,ncra));ALLOCATE(skybin(ntr,ncra))
     ALLOCATE(rephotbin(ntr,ncra));ALLOCATE(model_tr2bin(ntr,ncra))
  ENDIF

  !******************************************************************************
  !8. READING OF THE RV TIME-SERIES                                             !
  !   FORMAT: BJD, RV, ERROR, FWHM, C, BIS, BISERROR                            !
  !******************************************************************************

  IF(nrv.GT.0)THEN
     DO i=1,nrv       
        let1 = CHAR(48+int(i/1000))
        i2 = i - INT(i/1000)*1000
        let2 = CHAR(48+int(i2/100))
        i2 = i - INT(i/100)*100
        let3 = CHAR(48+int(i2/10))
        let4 = CHAR(48+mod(i2,10))
        fin =  let1 // let2 // let3 // let4 // '.txt'
        name2 = 'rv' // fin
        OPEN(UNIT=3,FILE=name2)
        nb=nprv(i)
        DO j=1,nb
           READ(3,*) rvbjd(i,j),rv(i,j),rverror(i,j),rvfwhm(i,j), &
                & rvcontrast(i,j),rvbis(i,j),rvloghk(i,j),rvtexp(i,j)
           rvtexp(i,j)=rvtexp(i,j)/(24*3600.)          ! assumes Texp in s 
           rv(i,j)=rv(i,j)*1000.                       ! assumes RV in km/s
           rverror(i,j)=rverror(i,j)*1000.             ! assumes dRV in km/s 
           rvbjd(i,j) = rvbjd(i,j) - 0.                ! assumes date = BJD-2450000.
        ENDDO
        CLOSE(3) 
        WRITE(*,'(A3,I2,A11,F8.3)') 'RV-',i, ' mean JD = ', SUM(rvbjd(i,:))/DBLE(nb)
     ENDDO
     PRINT*, ' '
  ENDIF

  IF(gibbs_sampler.EQ.'y')THEN
     ALLOCATE(gibbsonoff(njump));ALLOCATE(regul2(njump));ALLOCATE(jumpgibbs(na))
  ENDIF

  !******************************************************************************  
  !10. MCMC                                                                     !
  !******************************************************************************

  DO link = 1,na 
     IF(nfi.GT.0) testphef=0
     accepted(link)='y'
     IF(nrv.GT.0.AND.link.EQ.1)THEN
        DO i=1,nrv 
           nb = nprv(i)
           DO j=1,nb
              rverror(i,j) = SQRT(rverror(i,j)**2 + gjitter(i)**2)
           ENDDO
        ENDDO
     ENDIF
     IF(ntr.GT.0.AND.link.EQ.1)THEN
        DO i=1,ntr  
           nb = np(i)
           DO j=1,nb
              error(i,j) = error(i,j) * gbeta_red(i)
           ENDDO
        ENDDO
     ENDIF
     ! IF BEGINNING OF A NEW CHAIN, RANDOMIZE THE STARTING POINT
     IF(MOD(link-1,nat).EQ.0.AND.link.GE.1)THEN
        regul=1.
        IF(gibbs_sampler.EQ.'y')THEN
           DO i=1,njump-1
              gibbsonoff(i)=0
              regul2(i)=1.
           ENDDO
           gibbsonoff(njump)=1
           regul2(njump)=1.
        ENDIF
        CALL gasdev_s(harvest)
        mass_s(link) = starmass + dmass_s_ini*harvest
        CALL gasdev_s(harvest)
        if (priorteff.eq.'p') then
			temp_s(link) = spec(1) + sspec(1)*harvest
			dtemp_s = sspec(1)
        else
		    col_s(link)=col+I_col*harvest
		    dcol_s=I_col
        end if
        CALL gasdev_s(harvest)
        met_s(link) = spec(3) + sspec(3)*harvest
        dmet_s = sspec(3)
        IF(fitmsrs.EQ.'y'.OR.massfromr.EQ.'y'.OR.fixstellar.EQ.'y'.or.RjumpIso.eq.'y')THEN
          CALL gasdev_s(harvest)
          radius_s(link) = starradius + dstarradius_ini*harvest
          if (isoch.eq.'y'.and.((priorrad.eq.'y'.and..not.ybesR).or.priorrad.eq.'p')) then
		      Rinp=radius_s(link)
		  endif
        ENDIF
        IF(ngroup.GT.0)THEN
          DO i=1,ngroup
            CALL gasdev_s(harvest)
            edfgroup(i)=edfgroup_ini(i)
            dfgroup(i,link) = dfgroup_ini(i)+ edfgroup(i)*harvest
          ENDDO
        ENDIF
        DO i=1,npla
           test=0
           DO WHILE(test.LT.1)
             CALL gasdev_s(harvest)
             dsesinw(i)=dsesinw_ini(i)
             sesinw(i,link) = sesinw_ini(i)+ dsesinw(i)*harvest
             IF(ABS(sesinw(i,link)).LT.1) test=2
           ENDDO
           test=0
           DO WHILE(test.LT.1)
             CALL gasdev_s(harvest)
             dsecosw(i)=dsecosw_ini(i)
             secosw(i,link) = secosw_ini(i)+ dsecosw(i)*harvest
             IF(ABS(secosw(i,link)).LT.1) test=2
           ENDDO
           test=0
           exc(i,link)=secosw(i,link)**2 +sesinw(i,link)**2
           CALL gasdev_s(harvest)
           sper(i)=sper_ini(i)
           per(i,link) = per_ini(i) + sper(i)*harvest
           IF(per(i,link).GT.permax(i)) per(i,link)=permax(i)
           IF(per(i,link).LT.permin(i)) per(i,link)=permin(i)
           CALL gasdev_s(harvest)
           stidel(i)=stidel_ini(i)
           tidel(i,link) = tidel_ini(i) + stidel(i)*harvest
           IF(tidel(i,link).GT.tidelmax) tidel(i,link)=tidelmax
           IF(tidel(i,link).LT.1.) tidel(i,link)=1.
           sdF(i)=sdF_ini(i)
           CALL gasdev_s(harvest)
           sdF(i)=sdF_ini(i)
           dF(i,link) = dF_ini(i) + sdF(i)*harvest
           IF(dF(i,link).LT.0.) dF(i,link)=1.E-6 
           CALL gasdev_s(harvest)
           sdur(i)=sdur_ini(i)
           dur(i,link) = dur_ini(i) + sdur(i)*harvest
           IF(dur(i,link).LE.0.) dur(i,link)=1.E-8
           CALL gasdev_s(harvest)
           dkb(i)=dkb_ini(i)
           kb(i,link) = kb_ini(i) + dkb(i)*harvest
           IF(kb(i,link).LE.0.) kb(i,link)=1.E-6
           ka(i,link) = kb(i,link)/(SQRT(1.-exc(i,link)**2)*per(i,link)**(1./3.))
           IF(ka(i,link).GT.kmax(i))THEN 
             ka(i,link)=kmax(i)
             kb(i,link)=kmax(i)*(SQRT(1.-exc(i,link)**2)*per(i,link)**(1./3.))
           ENDIF 
           sb(i)=sb_ini(i)
           test=0
           DO WHILE(test.LT.1) 
              CALL gasdev_s(harvest)
              b(i,link) = b_ini(i) + sb(i)*harvest
              IF(b(i,link).GE.0.) test=2
           ENDDO
           IF(ntr.GT.0.AND.nddf.GT.0.AND.isddf.NE.'n')THEN
             DO j=1,nddf         
               CALL gasdev_s(harvest)
               eddf(i,j)=eddf_ini(i,j)
               ddf(i,j,link) = ddf(i,j,link) + eddf(i,j)*harvest
               IF(ddf(i,j,link).LT.-dF(i,link)) ddf(i,j,link)=0.
             ENDDO
           ENDIF
           CALL gasdev_s(harvest)
           st0(i)=st0_ini(i)
           t0(i,link) = t0_ini(i) + st0(i)*harvest
           IF(i.EQ.1.AND.ntiming.GT.0.AND.istiming.EQ.'y')THEN
              timerit(link) = 0.
              DO j=1,ntiming
                 test2 = (timing(j)-t0(i,link))/per(i,link)
                 epoch(i) = NINT(test2)
                 ttiming = t0(i,link) + epoch(i)*per(i,link)
                 omctiming(j) = timing(j) - ttiming
                 timerit(link) = timerit(link) + (omctiming(j)/stiming(j))**2
              ENDDO
           ENDIF
        ENDDO
        IF(limb.EQ.'qd'.AND.ntr.GT.0)THEN
           DO i=1,nfi
              k=0
              DO WHILE(k.LT.1)
                 ! CALL gasdev_s(harvest)
                 ql(i,1,link) = ql_ini(i,1) ! + sql(i,1)*harvest
                 ! CALL gasdev_s(harvest)
                 ql(i,2,link) = ql_ini(i,2) ! + sql(i,2)*harvest
                 IF(fitlimb(i).NE.'n')THEN
                    jumplimb(i,1,link) = 2*ql(i,1,link) + ql(i,2,link)
                    jumplimb(i,2,link) = ql(i,1,link) - 2*ql(i,2,link)
                    djumplimb(i,1)=djumplimb_ini(i,1)
                    djumplimb(i,2)=djumplimb_ini(i,2)
                 ENDIF
                 IF((3*jumplimb(i,1,link)-jumplimb(i,2,link)).LT.5)THEN
                    k=2
                 ENDIF
              ENDDO
           ENDDO
        ENDIF
        IF(limb.EQ.'nl'.AND.ntr.GT.0)THEN
           DO i=1,nfi           
              ! CALL gasdev_s(harvest)
              nl(i,1,link) = nl_ini(i,1) ! + harvest * snl(i,1)
              ! CALL gasdev_s(harvest)
              nl(i,2,link) = nl_ini(i,2) ! + harvest * snl(i,2)
              ! CALL gasdev_s(harvest)
              nl(i,3,link) = nl_ini(i,3) ! + harvest * snl(i,3)
              ! CALL gasdev_s(harvest)
              nl(i,4,link) = nl_ini(i,4) ! + harvest * snl(i,4)
           ENDDO
        ENDIF
        IF(ntr.GT.0.)THEN
         DO j=1,npla
           DO i=1,nfi         
              edFsec(i,j)=edFsec_ini(i,j)
              test=2
              DO WHILE(test.GT.1)
                CALL gasdev_s(harvest)
                dFsec(i,link,j) = dFsec_ini(i,j) + edFsec(i,j)*harvest
                IF(dFsec(i,link,j).GE.0.) test=0
              ENDDO
           ENDDO
          ENDDO
        ENDIF
        IF(ntr.GT.0.)THEN
           DO i=1,nfi     
              CALL gasdev_s(harvest)
              ephampli1(i)=ephampli1_ini(i)
              ephampli2(i)=ephampli2_ini(i)
              ephampli3(i)=ephampli3_ini(i)
              ephoffset(i)=ephoffset_ini(i)
              phoffset(i,link) = phoffset_ini(i) + ephoffset(i)*harvest
              IF(phoffset(i,link).LT.-180.) phoffset(i,link)=phoffset(i,link)+360.
              IF(phoffset(i,link).GT.180.)  phoffset(i,link)=phoffset(i,link)-360.
              test=2
              DO WHILE(test.GT.1)
                CALL gasdev_s(harvest)
                phampli1(i,link) = phampli1_ini(i) + ephampli1(i)*harvest
                IF(dfcond.EQ.'y')THEN
                  IF(phampli1(i,link).GE.0.AND.phampli1(i,link).LE.dFsec(i,link,1)) test=0
                ELSE
                  IF(phampli1(i,link).GE.0.) test=0
                ENDIF
              ENDDO
              test=2
              DO WHILE(test.GT.1)
                CALL gasdev_s(harvest)
                phampli2(i,link) = phampli2_ini(i) + ephampli2(i)*harvest
                IF(phampli2(i,link).GE.0.) test=0
              ENDDO
              test=2
              DO WHILE(test.GT.1)
                CALL gasdev_s(harvest)
                phampli3(i,link) = phampli3_ini(i) + ephampli3(i)*harvest
                IF(phampli3(i,link).GE.0.) test=0
              ENDDO            
           ENDDO
        ENDIF
        IF(ntr.GT.0)THEN
           DO i=1,ntr
              IF(sinusnumber(i).GT.0.AND.testsin.GT.0)THEN
                 k=sinusnumber(i)
                 DO j=1,k
                   CALL gasdev_s(harvest)
                   esit0(i,j)=esit0_ini(i,j) 
                   sit0(i,j,link) = sit0_ini(i,j) + esit0(i,j)*harvest
                   CALL gasdev_s(harvest)
                   esip(i,j)=esip_ini(i,j)
                   sip(i,j,link) = sip_ini(i,j) + esip(i,j)*harvest
                   IF(sip(i,j,link).LT.0.) sip(i,j,link)=sip_ini(i,j)
                 ENDDO
              ENDIF
              IF(flarenumber(i).GT.0.AND.testflare.GT.0)THEN
                 k=flarenumber(i)
                 DO j=1,k
                   CALL gasdev_s(harvest)
                   eflampli(i,j)=eflampli_ini(i,j) 
                   flampli(i,j,link) = flampli_ini(i,j) + eflampli(i,j)*harvest
                   CALL gasdev_s(harvest)
                   efltau(i,j)=efltau_ini(i,j)
                   fltau(i,j,link) = fltau_ini(i,j) + efltau(i,j)*harvest
                   CALL gasdev_s(harvest)
                   eflt0(i,j)=eflt0_ini(i,j)
                   flt0(i,j,link) = flt0_ini(i,j) + eflt0(i,j)*harvest
                 ENDDO
              ENDIF
              IF(ramporder(i).GT.0.AND.rampmod.EQ.'exp')THEN
                 k=ramporder(i)
                 test=0
                 DO WHILE(test.LT.1)
                   CALL gasdev_s(harvest)
                   et1ramp(i)=et1ramp_ini(i) 
                   t1ramp(i,link) = t1ramp_ini(i) + et1ramp(i)*harvest
                   IF(t1ramp(i,link).GT.0) test=2
                 ENDDO
                 IF(k.GT.1)THEN
                   test=0
                   DO WHILE(test.LT.1)
                     CALL gasdev_s(harvest)
                     et2ramp(i)=et2ramp_ini(i) 
                     t2ramp(i,link) = t2ramp_ini(i) + et2ramp(i)*harvest
                     IF(t2ramp(i,link).GT.0) test=2
                 ENDDO
                 ENDIF
              ENDIF
           ENDDO
        ENDIF

        IF(nrv.GT.0)THEN     
           CALL gasdev_s(harvest)
           rold(link,1) = rold_ini(1) + srold(1)*harvest
           CALL gasdev_s(harvest)
           rold(link,2) = rold_ini(2) + srold(2)*harvest
           test=0
           DO WHILE(test.EQ.0)
              CALL gasdev_s(harvest)
              dsvsinisinbeta=dsvsinisinbeta_ini
              svsinisinbeta(link) = ini_svsinisinbeta +dsvsinisinbeta*harvest
              CALL gasdev_s(harvest)
              dsvsinicosbeta=dsvsinicosbeta_ini
              svsinicosbeta(link) = ini_svsinicosbeta +dsvsinicosbeta*harvest
              vsini(link)= svsinicosbeta(link)**2+ svsinisinbeta(link)**2
              IF(vsini(link).GT.0.)test=1 
           ENDDO
           IF(ABS(svsinicosbeta(link)).GT.1.E-13)THEN
              beta(link) = DATAN(svsinisinbeta(link)/svsinicosbeta(link))
           ELSE
              beta(link) = DASIN(svsinisinbeta(link)/SQRT(vsini(link)))
           ENDIF
           IF(svsinicosbeta(link).LT.0)THEN
              IF(svsinisinbeta(link).LT.0)THEN
                 beta(link)=beta(link)+pi
              ELSE
                 beta(link)=beta(link)-pi
              ENDIF
           ENDIF
           beta(link)=beta(link)*180/pi
        ELSE
          CALL gasdev_s(harvest)
          vsini(link)=spec(5) + harvest*sspec(5)    
        ENDIF
        IF(testf2.EQ.'y')THEN
           test=0
           sf2 = sspec(7) 
           DO WHILE(test.LT.1)
              CALL gasdev_s(harvest)
              f2(link) = spec(7)+harvest*sspec(7) 
              IF(f2(link).GT.0) test=2
           ENDDO
        ENDIF

       ! TTV
        IF(nttvmax.GT.0.AND.isttv.EQ.'y')THEN
          DO i=1,npla
            DO j=1,nttv(i)
              CALL gasdev_s(harvest)
              ttv(i,j,link)=ttv_ini(i,j)+sttv_ini(i,j)*harvest
              sttv(i,j)=sttv_ini(i,j)
              ttr(i,j,link)=t0(i,link)+per(i,link)*epochtr(i,j)+ttv(i,j,link)
            ENDDO
          ENDDO
        ENDIF

  ! ELSE, MAKE A NEW JUMP OF THE PARAMETER STATE

     ELSE IF(link.GT.1)THEN
        IF(gibbs_sampler.EQ.'y'.AND.MOD(link,nat).LE.nburn2)THEN
          gibbscount=0
          DO i=1,njump-1
            IF(gibbsonoff(i).EQ.1) jumpgibbs(link)=i+1
          ENDDO
          IF(gibbsonoff(njump).EQ.1) jumpgibbs(link)=1 
        ENDIF
        IF(fixstellar.EQ.'y')THEN
          CALL gasdev_s(harvest)
          mass_s(link) = starmass + dmass_s_ini*harvest
          CALL gasdev_s(harvest)
          temp_s(link) = spec(1) + sspec(1)*harvest
          !col_s is only available with isoch that is incompatible with fixstellar
          CALL gasdev_s(harvest)
          met_s(link) = spec(3) + sspec(3)*harvest
          CALL gasdev_s(harvest)
          radius_s(link) = starradius + dstarradius_ini*harvest
        ELSE
          IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
            CALL gasdev_s(harvest)
            if (priorteff.eq.'p') then
		        temp_s(link) = temp_s(link-1) + regul*dtemp_s*harvest
		        IF(temp_s(link).LT.0.)THEN
		          temp_s(link)=temp_s(link-1)
		          accepted(link)='n'
		        ENDIF
            else
            	col_s(link)=col_s(link-1)+regul*dcol_s*harvest
            	if (col_s(link).lt.BmVinf.or.col_s(link).gt.BmVsup) then
            	  col_s(link)=col_s(link-1)
            	  accepted(link)='n'
            	end if
            end if
          ELSE
            gibbscount=gibbscount+1
            IF(jumpgibbs(link).EQ.gibbscount)THEN
              CALL gasdev_s(harvest)
              if (priorteff.eq.'p') then
			    temp_s(link) = temp_s(link-1)+ regul2(gibbscount)*sspec(1)*harvest 
			    IF(temp_s(link).LT.0.)THEN
			      temp_s(link)=temp_s(link-1)
			      accepted(link)='n'
			    ENDIF
			    dtemp_s=regul2(gibbscount)*sspec(1)
		      else
		      	col_s(link)=col_s(link-1)+regul2(gibbscount)*I_col*harvest
		      	if (col_s(link).lt.BmVinf.or.col_s(link).gt.BmVsup) then
            	  col_s(link)=col_s(link-1)
            	  accepted(link)='n'
            	end if
            	dcol_s=regul2(gibbscount)*I_col
		      end if
              gibbsonoff(gibbscount) = 1
            ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
              if (priorteff.eq.'p') then
              	temp_s(link) = temp_s(link-1)
              else
              	col_s(link)=col_s(link-1)
              end if
              gibbsonoff(gibbscount)=0
            ELSE
              if (priorteff.eq.'p') then
              	temp_s(link) = temp_s(link-1)
              else
                col_s(link)=col_s(link-1)
              end if
            ENDIF
          ENDIF
          IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
            CALL gasdev_s(harvest)
            met_s(link) = met_s(link-1) + regul*dmet_s*harvest
            if (priormet.eq.'y'.and.(met_s(link).lt.FeHinf.or.met_s(link).gt.FeHsup)) then
              met_s(link)=met_s(link-1)
              accepted(link)='n'
            end if
          ELSE
            gibbscount=gibbscount+1
            IF(jumpgibbs(link).EQ.gibbscount)THEN
              CALL gasdev_s(harvest)
              met_s(link) = met_s(link-1)+ regul2(gibbscount)*sspec(3)*harvest 
              if (priormet.eq.'y'.and.(met_s(link).lt.FeHinf.or.met_s(link).gt.FeHsup)) then
		        met_s(link)=met_s(link-1)
		        accepted(link)='n'
		      end if
              dmet_s=regul2(gibbscount)*sspec(3) 
              gibbsonoff(gibbscount) = 1
            ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
              met_s(link) = met_s(link-1)
              gibbsonoff(gibbscount)=0
            ELSE
              met_s(link) = met_s(link-1)
            ENDIF
          ENDIF
          IF(msrs.EQ.'y'.OR.massfromr.EQ.'y'.OR.enoch.EQ.'y')THEN !!.or.isoch.eq.'y'
            CALL gasdev_s(harvest)
            mass_s(link) = starmass + dmass_s_ini*harvest
            IF(mass_s(link).LT.0.)THEN
              mass_s(link)=mass_s(link-1)
              accepted(link)='n'
            ENDIF
          ELSE
            IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
              CALL gasdev_s(harvest)
              mass_s(link) = mass_s(link-1) + regul*dmass_s*harvest
              IF(mass_s(link).LT.0.)THEN
                 mass_s(link)=mass_s(link-1)
                accepted(link)='n'
              ENDIF
            ELSE
              gibbscount=gibbscount+1
              IF(jumpgibbs(link).EQ.gibbscount)THEN
                CALL gasdev_s(harvest)
                mass_s(link) = mass_s(link-1)+ regul2(gibbscount)*dmass_s_ini*harvest 
                IF(mass_s(link).LT.0.)THEN
                  mass_s(link)=mass_s(link-1)
                  accepted(link)='n'
                ENDIF              
                dmass_s=regul2(gibbscount)*dmass_s_ini
                gibbsonoff(gibbscount) = 1
              ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                mass_s(link) = mass_s(link-1)
                gibbsonoff(gibbscount)=0
              ELSE
                mass_s(link) = mass_s(link-1)
              ENDIF
            ENDIF
          ENDIF
        ENDIF
        IF(fitmsrs.EQ.'y'.OR.massfromr.EQ.'y'.or.RjumpIso.eq.'y')THEN
          IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
            CALL gasdev_s(harvest)
            radius_s(link) = radius_s(link-1) + regul*dstarradius*harvest
            IF(radius_s(link).LT.0.)THEN
               radius_s(link)=radius_s(link-1)
               accepted(link)='n'
            ENDIF
          ELSE
            gibbscount=gibbscount+1
            IF(jumpgibbs(link).EQ.gibbscount)THEN
              CALL gasdev_s(harvest)
              radius_s(link) = radius_s(link-1)+ regul2(gibbscount)*dstarradius_ini*harvest 
              IF(radius_s(link).LT.0.)THEN
                 radius_s(link)=radius_s(link-1)
                 accepted(link)='n'
              ENDIF
              dstarradius=regul2(gibbscount)*dstarradius_ini
              gibbsonoff(gibbscount) = 1
            ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
              radius_s(link) = radius_s(link-1)
              gibbsonoff(gibbscount)=0
            ELSE
              radius_s(link) = radius_s(link-1)
            ENDIF
          ENDIF
          if (isoch.eq.'y'.and.((priorrad.eq.'y'.and..not.ybesR).or.priorrad.eq.'p')) then
		      Rinp=radius_s(link)
		      I_Rinp=dstarradius
		  endif
        ENDIF
        IF(nrv.GT.0)THEN
           rold(link,1) = rold_ini(1) 
           rold(link,2) = rold_ini(2) 
        ENDIF
        DO i=1,npla
           IF(isjump(i,6).EQ.'n')THEN
              sesinw(i,link) = sesinw_ini(i)   
              secosw(i,link) = secosw_ini(i)   
              exc(i,link) = secosw(i,link)**2 + sesinw(i,link)**2
           ELSE
              IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                CALL gasdev_s(harvest)
                sesinw(i,link) = sesinw(i,link-1)+ regul*dsesinw(i)*harvest
              ELSE 
                 gibbscount=gibbscount+1
                 IF(jumpgibbs(link).EQ.gibbscount)THEN
                    CALL gasdev_s(harvest)
                    sesinw(i,link) = sesinw(i,link-1)+ regul2(gibbscount)*dsesinw_ini(i)*harvest 
                    dsesinw(i)=regul2(gibbscount)*dsesinw_ini(i)
                    gibbsonoff(gibbscount) = 1
                 ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                    sesinw(i,link) = sesinw(i,link-1)
                    gibbsonoff(gibbscount)=0
                 ELSE
                    sesinw(i,link) = sesinw(i,link-1)
                 ENDIF
              ENDIF
              IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                 CALL gasdev_s(harvest)
                 secosw(i,link) = secosw(i,link-1)+ regul*dsecosw(i)*harvest
              ELSE 
                 gibbscount=gibbscount+1
                 IF(jumpgibbs(link).EQ.gibbscount)THEN
                    CALL gasdev_s(harvest)
                    secosw(i,link) = secosw(i,link-1)+ regul2(gibbscount)*dsecosw_ini(i)*harvest
                    dsecosw(i)=regul2(gibbscount)* dsecosw_ini(i)
                    gibbsonoff(gibbscount) = 1
                 ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                    secosw(i,link) = secosw(i,link-1)
                    gibbsonoff(gibbscount)=0
                 ELSE
                    secosw(i,link) = secosw(i,link-1)
                 ENDIF
              ENDIF
              exc(i,link) = secosw(i,link)**2 + sesinw(i,link)**2
!!              IF(exc(i,link).GE.1.)exc(i,link)=0.999
           ENDIF
           IF(ABS(sesinw(i,link)).GE.1)THEN
             accepted(link)='n'    
             sesinw(i,link)=sesinw(i,link-1) 
           ENDIF
           IF(ABS(secosw(i,link)).GE.1)THEN
             accepted(link)='n'
             secosw(i,link)=secosw(i,link-1) 
           ENDIF
           if (exc(i,link).lt.0.or.exc(i,link).ge.1) then
	       	 accepted(link)='n'
	       	 exc(i,link)=exc(i,link-1)
	       endif
        ENDDO
        IF(nrv.GT.0)THEN
           IF(isjump(1,8).EQ.'n')THEN
              svsinisinbeta(link) = ini_svsinisinbeta
              svsinicosbeta(link) = ini_svsinicosbeta
              CALL gasdev_s(harvest)
              vsini(link)=spec(5) + harvest*sspec(5)    
           ELSE         
              IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                 CALL gasdev_s(harvest)
                 svsinisinbeta(link) = svsinisinbeta(link-1) + regul*dsvsinisinbeta*harvest
              ELSE
                 gibbscount=gibbscount+1
                 IF(jumpgibbs(link).EQ.gibbscount)THEN
                    CALL gasdev_s(harvest)
                    svsinisinbeta(link) = svsinisinbeta(link-1) +regul2(gibbscount)*dsvsinisinbeta_ini*harvest
                    dsvsinisinbeta=regul2(gibbscount)*dsvsinisinbeta_ini
                    gibbsonoff(gibbscount) = 1
                 ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                    svsinisinbeta(link) = svsinisinbeta(link-1)
                    gibbsonoff(gibbscount) = 0
                 ELSE
                    svsinisinbeta(link) = svsinisinbeta(link-1)
                 ENDIF
              ENDIF
              IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                 CALL gasdev_s(harvest)
                 svsinicosbeta(link) = svsinicosbeta(link-1)+regul*dsvsinicosbeta*harvest
              ELSE
                 gibbscount=gibbscount+1
                 IF(jumpgibbs(link).EQ.gibbscount)THEN
                    CALL gasdev_s(harvest)
                    svsinicosbeta(link) = svsinicosbeta(link-1) + regul2(gibbscount)*dsvsinicosbeta_ini*harvest
                    dsvsinicosbeta=regul2(gibbscount)*dsvsinicosbeta_ini
                    gibbsonoff(gibbscount) = 1
                 ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                    svsinicosbeta(link) = svsinicosbeta(link-1)
                    gibbsonoff(gibbscount) = 0
                 ELSE
                    svsinicosbeta(link) = svsinicosbeta(link-1)
                 ENDIF
              ENDIF
              vsini(link)= svsinicosbeta(link)**2+ svsinisinbeta(link)**2
           ENDIF
           IF(ABS(svsinicosbeta(link)).GT.1.E-13)THEN
              beta(link) = DATAN(svsinisinbeta(link)/svsinicosbeta(link))
           ELSE
              beta(link) = DASIN(svsinisinbeta(link)/SQRT(vsini(link)))
           ENDIF
           IF(svsinicosbeta(link).LT.0)THEN
              IF(svsinisinbeta(link).LT.0)THEN
                 beta(link)=beta(link)+pi
              ELSE
                 beta(link)=beta(link)-pi
              ENDIF
           ENDIF
           beta(link)=beta(link)*180/pi
        ELSE
          CALL gasdev_s(harvest)
          vsini(link)=spec(5) + harvest*sspec(5)    
        ENDIF
        IF(ngroup.GT.0)THEN
          DO i=1,ngroup
            IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
              CALL gasdev_s(harvest)
              dfgroup(i,link)=dfgroup_ini(i)+regul*edfgroup(i)*harvest
              IF(dfgroup(i,link).LT.0)THEN
                dfgroup(i,link)=dfgroup(i,link-1)
                accepted(link)='n'
              ENDIF
            ELSE
              gibbscount=gibbscount+1
              IF(jumpgibbs(link).EQ.gibbscount)THEN
                CALL gasdev_s(harvest)
                dfgroup(i,link) = dfgroup(i,link-1) + regul2(gibbscount)*edfgroup_ini(i)*harvest
                IF(dfgroup(i,link).LT.0.)THEN
                  dfgroup(i,link)=dfgroup(i,link-1)
                  accepted(link)='n'
                ENDIF  
                edfgroup(i)=regul2(gibbscount)*edfgroup_ini(i)
                gibbsonoff(gibbscount) = 1
              ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                dfgroup(i,link) = dfgroup(i,link-1) 
                gibbsonoff(gibbscount) = 0
              ELSE
               dfgroup(i,link) = dfgroup(i,link-1) 
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        DO i=1,npla
           IF(isjump(i,1).EQ.'n')THEN
             CALL gasdev_s(harvest)
             dF(i,link) = dF_ini(i)+harvest*sdF_ini(i)
             IF(dF(i,link).LT.0.)THEN
               dF(i,link)=dF(i,link-1)
               accepted(link)='n'
             ENDIF
           ELSE
              IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                 CALL gasdev_s(harvest)
                 dF(i,link) = dF(i,link-1) + regul*sdF(i)*harvest
                 IF(dF(i,link).LT.0.)THEN
                   dF(i,link)=dF(i,link-1)
                   accepted(link)='n'
                 ENDIF       
              ELSE
                 gibbscount=gibbscount+1
                 IF(jumpgibbs(link).EQ.gibbscount)THEN
                    CALL gasdev_s(harvest)
                    dF(i,link) = dF(i,link-1) + regul2(gibbscount)*sdF_ini(i)*harvest
                    IF(dF(i,link).LT.0.)THEN
                      dF(i,link)=dF(i,link-1)
                      accepted(link)='n'
                    ENDIF  
                    sdF(i)=regul2(gibbscount)*sdF_ini(i)
                    gibbsonoff(gibbscount) = 1
                 ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                    dF(i,link) = dF(i,link-1) 
                    gibbsonoff(gibbscount) = 0
                 ELSE
                    dF(i,link) = dF(i,link-1) 
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
        IF(ntr.GT.0)THEN
          DO i=1,nfi
            DO j=1,npla
              IF(iffitoc(i).EQ.'n')THEN
                 CALL gasdev_s(harvest)
                 dFsec(i,link,j) = dFsec_ini(i,j) + harvest*edFsec_ini(i,j) 
              ELSE  
                 IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                    CALL gasdev_s(harvest)
                    dFsec(i,link,j) = dFsec(i,link-1,j) + regul*edFsec(i,j)*harvest
                    IF(dFsec(i,link,j).LT.0.)THEN
                      dFsec(i,link,j)=dFsec(i,link-1,j)
                      accepted(link)='n'
                    ENDIF  
                 ELSE
                    gibbscount=gibbscount+1
                    IF(jumpgibbs(link).EQ.gibbscount)THEN
                       CALL gasdev_s(harvest)
                       dFsec(i,link,j) = dFsec(i,link-1,j) + regul2(gibbscount)*edFsec_ini(i,j)*harvest
                       IF(dFsec(i,link,j).LT.0.)THEN
                         dFsec(i,link,j)=dFsec(i,link-1,j)
                         accepted(link)='n'
                       ENDIF  
                       edFsec(i,j)=regul2(gibbscount)*edFsec_ini(i,j)
                       gibbsonoff(gibbscount) = 1
                    ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                       gibbsonoff(gibbscount) = 0
                       dFsec(i,link,j) = dFsec(i,link-1,j)
                    ELSE
                       dFsec(i,link,j) = dFsec(i,link-1,j)
                    ENDIF
                 ENDIF
              ENDIF
            ENDDO
              IF(iffitph(i).EQ.'n'.OR.ephampli1_ini(i).LE.1E-10)THEN
                 CALL gasdev_s(harvest)
                 phampli1(i,link) = phampli1_ini(i) + harvest*ephampli1_ini(i) 
                 IF(phampli1(i,link).LT.0.)THEN
                   phampli1(i,link)=phampli1(i,link-1)
                   accepted(link)='n'
                 ENDIF                                
                 IF(dfcond.EQ.'y'.AND.phampli1(i,link).GT.dFsec(i,link,1))THEN
                    phampli1(i,link)=phampli1(i,link-1)
                    accepted(link)='n'
                 ENDIF
              ELSE   
                 IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN  
                    CALL gasdev_s(harvest)
                    phampli1(i,link) = phampli1(i,link-1) + regul*ephampli1(i)*harvest
                    IF(phampli1(i,link).LT.0.)THEN
                      phampli1(i,link)=phampli1(i,link-1)
                      accepted(link)='n'
                    ENDIF  
                    IF(dfcond.EQ.'y'.AND.phampli1(i,link).GT.dFsec(i,link,1))THEN
                       phampli1(i,link)=phampli1(i,link-1)
                       accepted(link)='n'
                    ENDIF
                 ELSE
                    gibbscount=gibbscount+1
                    IF(jumpgibbs(link).EQ.gibbscount)THEN
                       CALL gasdev_s(harvest)
                       phampli1(i,link) = phampli1(i,link-1) + regul2(gibbscount)*ephampli1_ini(i)*harvest
                       IF(phampli1(i,link).LT.0.)THEN
                         phampli1(i,link)=phampli1(i,link-1)
                         accepted(link)='n'
                       ENDIF  
                       IF(dfcond.EQ.'y'.AND.phampli1(i,link).GT.dFsec(i,link,1))THEN
                         phampli1(i,link)=phampli1(i,link-1)
                         accepted(link)='n'
                       ENDIF
                       ephampli1(i)=regul2(gibbscount)*ephampli1_ini(i)
                       gibbsonoff(gibbscount) = 1
                    ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                       gibbsonoff(gibbscount) = 0
                       phampli1(i,link) = phampli1(i,link-1)
                    ELSE
                       phampli1(i,link) = phampli1(i,link-1)
                    ENDIF
                 ENDIF
              ENDIF
              IF(iffitph(i).EQ.'n'.OR.ephampli2_ini(i).LT.1E-10)THEN
                 CALL gasdev_s(harvest)
                 phampli2(i,link) = phampli2_ini(i) + harvest*ephampli2_ini(i) 
                 IF(phampli2(i,link).LT.0.)THEN
                   phampli2(i,link)=phampli2(i,link-1)
                   accepted(link)='n'
                 ENDIF  
              ELSE  
                 IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN   
                    CALL gasdev_s(harvest)
                    phampli2(i,link) = phampli2(i,link-1) + regul*ephampli2(i)*harvest
                    IF(phampli2(i,link).LT.0.)THEN
                      phampli2(i,link)=phampli2(i,link-1)
                      accepted(link)='n'
                    ENDIF  
                 ELSE
                    gibbscount=gibbscount+1
                    IF(jumpgibbs(link).EQ.gibbscount)THEN
                       CALL gasdev_s(harvest)
                       phampli2(i,link) = phampli2(i,link-1) + regul2(gibbscount)*ephampli2_ini(i)*harvest
                       IF(phampli2(i,link).LT.0.)THEN
                         phampli2(i,link)=phampli2(i,link-1)
                         accepted(link)='n'
                       ENDIF  
                       ephampli2(i)=regul2(gibbscount)*ephampli2_ini(i)
                       gibbsonoff(gibbscount) = 1
                    ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                       gibbsonoff(gibbscount) = 0
                       phampli2(i,link) = phampli2(i,link-1)
                    ELSE
                       phampli2(i,link) = phampli2(i,link-1)
                    ENDIF
                 ENDIF
              ENDIF
              IF(iffitph(i).EQ.'n'.OR.ephampli3_ini(i).LT.1E-10)THEN
                 CALL gasdev_s(harvest)
                 phampli3(i,link) = phampli3_ini(i) + harvest*ephampli3_ini(i) 
                 IF(phampli3(i,link).LT.0.)THEN
                   phampli3(i,link)=phampli3(i,link-1)
                   accepted(link)='n'
                 ENDIF  
              ELSE  
                 IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN  
                    CALL gasdev_s(harvest)
                    phampli3(i,link) = phampli3(i,link-1) + regul*ephampli3(i)*harvest
                    IF(phampli3(i,link).LT.0.)THEN
                      phampli3(i,link)=phampli3(i,link-1)
                      accepted(link)='n'
                    ENDIF 
                 ELSE
                    gibbscount=gibbscount+1
                    IF(jumpgibbs(link).EQ.gibbscount)THEN
                       CALL gasdev_s(harvest)
                       phampli3(i,link) = phampli3(i,link-1) + regul2(gibbscount)*ephampli3_ini(i)*harvest
                       IF(phampli3(i,link).LT.0.)THEN
                         phampli3(i,link)=phampli3(i,link-1)
                         accepted(link)='n'
                       ENDIF 
                       ephampli3(i)=regul2(gibbscount)*ephampli3_ini(i)
                       gibbsonoff(gibbscount) = 1
                    ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                       gibbsonoff(gibbscount) = 0
                       phampli3(i,link) = phampli3(i,link-1)
                    ELSE
                       phampli3(i,link) = phampli3(i,link-1)
                    ENDIF
                 ENDIF
              ENDIF
              IF(iffitph(i).EQ.'n'.OR.ephoffset_ini(i).LT.1E-10)THEN
                 CALL gasdev_s(harvest)
                 phoffset(i,link) = phoffset_ini(i) + harvest*ephoffset_ini(i) 
              ELSE  
                 IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                    CALL gasdev_s(harvest)
                    phoffset(i,link) = phoffset(i,link-1) + regul*ephoffset(i)*harvest
                    IF(phoffset(i,link).GT.180.) phoffset(i,link)=phoffset(i,link)-360.
                    IF(phoffset(i,link).LT.-180.) phoffset(i,link)=phoffset(i,link)+360.
                 ELSE
                    gibbscount=gibbscount+1
                    IF(jumpgibbs(link).EQ.gibbscount)THEN   
                      CALL gasdev_s(harvest)
                      phoffset(i,link) = phoffset(i,link-1) + regul2(gibbscount)*ephoffset_ini(i)*harvest
                      IF(phoffset(i,link).GT.180.) phoffset(i,link)=phoffset(i,link)-360.
                      IF(phoffset(i,link).LT.-180.) phoffset(i,link)=phoffset(i,link)+360.  
                      ephoffset(i)=regul2(gibbscount)*ephoffset_ini(i)
                      gibbsonoff(gibbscount) = 1
                    ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                       gibbsonoff(gibbscount) = 0
                       phoffset(i,link) = phoffset(i,link-1)
                    ELSE
                       phoffset(i,link) = phoffset(i,link-1)
                    ENDIF
                 ENDIF
              ENDIF
           ENDDO
        ENDIF
        DO i=1,npla
          IF(nddf.GT.0.AND.ntr.GT.0.AND.isddf.NE.'n')THEN
           DO j=1,nddf
              IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                 CALL gasdev_s(harvest)
                 ddf(i,j,link) = ddf(i,j,link-1) + regul*eddf(i,j)*harvest
                 IF(ddf(i,j,link).LT.-dF(i,link))THEN
                   ddf(i,j,link)=ddf(i,j,link-1)
                   accepted(link)='n'
                 ENDIF              
              ELSE
                 gibbscount=gibbscount+1
                 IF(jumpgibbs(link).EQ.gibbscount)THEN
                    CALL gasdev_s(harvest)
                    ddf(i,j,link) = ddf(i,j,link-1) + regul2(gibbscount)*eddf_ini(i,j)*harvest
                    IF(ddf(i,j,link).LT.-dF(i,link))THEN
                      ddf(i,j,link)=ddf(i,j,link-1)
                      accepted(link)='n'
                    ENDIF           
                    eddf(i,j)=regul2(gibbscount)*eddf_ini(i,j)
                    gibbsonoff(gibbscount) = 1
                 ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                    gibbsonoff(gibbscount) = 0
                    ddf(i,j,link) = ddf(i,j,link-1)
                 ELSE
                    ddf(i,j,link) = ddf(i,j,link-1)
                 ENDIF
              ENDIF
           ENDDO
          ENDIF
           IF(isjump(i,5).EQ.'n')THEN
              per(i,link) = per_ini(i) 
              IF(per(i,link).GT.permax(i)) per(i,link)=permax(i)
              IF(per(i,link).LT.permin(i)) per(i,link)=permin(i)
           ELSE 
              IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                 CALL gasdev_s(harvest)
                 per(i,link) = per(i,link-1) + regul*sper(i)*harvest 
                 IF(per(i,link).GT.permax(i).OR.per(i,link).LT.permin(i))THEN
                   per(i,link)=per(i,link-1)
                   accepted(link)='n'
                 ENDIF               
              ELSE
                 gibbscount=gibbscount+1  
                 IF(jumpgibbs(link).EQ.gibbscount)THEN
                    CALL gasdev_s(harvest)
                    per(i,link) = per(i,link-1) + regul2(gibbscount)*sper_ini(i)*harvest 
                    IF(per(i,link).GT.permax(i).OR.per(i,link).LT.permin(i))THEN
                       per(i,link)=per(i,link-1)
                       accepted(link)='n'
                    ENDIF     
                    sper(i)=regul2(gibbscount)*sper_ini(i)
                    gibbsonoff(gibbscount) = 1
                 ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                    gibbsonoff(gibbscount) = 0
                    per(i,link) = per(i,link-1)
                 ELSE
                    per(i,link) = per(i,link-1)
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
        DO i=1,npla
           IF(isjump(i,7).EQ.'n')THEN
              kb(i,link) = kb_ini(i)
              IF(kb(i,link).LE.0.) kb(i,link)=1.E-6
              ka(i,link) = kb(i,link)/(SQRT(1.-exc(i,link)**2)*per(i,link)**(1./3.))
              IF(ka(i,link).GT.kmax(i))THEN
                ka(i,link)=kmax(i)
                kb(i,link)=kmax(i)*(SQRT(1.-exc(i,link)**2)*per(i,link)**(1./3.))
             ENDIF 
           ELSE
              IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                 CALL gasdev_s(harvest)
                 kb(i,link) = kb(i,link-1) + regul*dkb(i)*harvest
                 IF(kb(i,link).LT.0.)THEN
                   kb(i,link)=kb(i,link-1)
                   accepted(link)='n'
                 ENDIF   
                 ka(i,link) = kb(i,link)/(SQRT(1.-exc(i,link)**2)*per(i,link)**(1./3.))
                 IF(ka(i,link).GT.kmax(i))THEN
                   ka(i,link)=ka(i,link-1)
                   kb(i,link)=kb(i,link-1)
                   accepted(link)='n'
                 ENDIF 
              ELSE
                 gibbscount=gibbscount+1
                 IF(jumpgibbs(link).EQ.gibbscount)THEN
                   CALL gasdev_s(harvest)
                   kb(i,link) = kb(i,link-1)+regul2(gibbscount)*dkb_ini(i)*harvest
                   IF(kb(i,link).LT.0.)THEN
                     kb(i,link)=kb(i,link-1)
                     accepted(link)='n'
                   ENDIF   
                   ka(i,link) = kb(i,link)/(SQRT(1.-exc(i,link)**2)*per(i,link)**(1./3.))
                   IF(ka(i,link).GT.kmax(i))THEN
                     ka(i,link)=ka(i,link-1)
                     kb(i,link)=kb(i,link-1)
                     accepted(link)='n'
                   ENDIF
                   dkb(i)=regul2(gibbscount)*dkb_ini(i)
                   gibbsonoff(gibbscount) = 1
                 ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                   gibbsonoff(gibbscount) = 0
                   kb(i,link) = kb(i,link-1)
                   ka(i,link) = kb(i,link)/(SQRT(1.-exc(i,link)**2)*per(i,link)**(1./3.))
                 ELSE
                   kb(i,link) = kb(i,link-1)
                   ka(i,link) = kb(i,link)/(SQRT(1.-exc(i,link)**2)*per(i,link)**(1./3.))
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
        DO i=1,npla
           IF(isjump(i,9).EQ.'n')THEN
              tidel(i,link) = tidel_ini(i) 
              IF(tidel(i,link).GT.tidelmax) tidel(i,link)=tidelmax
              IF(tidel(i,link).LT.1.) tidel(i,link)=1.
           ELSE 
              IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                 CALL gasdev_s(harvest)
                 tidel(i,link) = tidel(i,link-1) + regul*stidel(i)*harvest 
                 IF(tidel(i,link).GT.tidelmax.OR.tidel(i,link).LT.1.)THEN
                   tidel(i,link)=tidel(i,link-1)
                   accepted(link)='n'
                 ENDIF               
              ELSE
                 gibbscount=gibbscount+1  
                 IF(jumpgibbs(link).EQ.gibbscount)THEN
                    CALL gasdev_s(harvest)
                    tidel(i,link) = tidel(i,link-1) + regul2(gibbscount)*stidel_ini(i)*harvest 
                    IF(tidel(i,link).GT.tidelmax.OR.tidel(i,link).LT.1.)THEN
                       tidel(i,link)=tidel(i,link-1)
                       accepted(link)='n'
                    ENDIF     
                    stidel(i)=regul2(gibbscount)*stidel_ini(i)
                    gibbsonoff(gibbscount) = 1
                 ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                    gibbsonoff(gibbscount) = 0
                    tidel(i,link) = tidel(i,link-1)
                 ELSE
                    tidel(i,link) = tidel(i,link-1)
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
        DO i=1,npla
           IF(isjump(i,3).EQ.'n')THEN
              dur(i,link) = dur_ini(i) 
           ELSE          
              IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                 CALL gasdev_s(harvest)
                 dur(i,link) = dur(i,link-1) + regul*sdur(i)*harvest
                 IF(dur(i,link).LT.0)THEN
                   dur(i,link)=dur(i,link-1)
                   accepted(link)='n'
                 ENDIF     
              ELSE
                 gibbscount=gibbscount+1
                 IF(jumpgibbs(link).EQ.gibbscount)THEN
                    CALL gasdev_s(harvest)
                    dur(i,link) = dur(i,link-1) + regul2(gibbscount)*sdur_ini(i)*harvest
                    IF(dur(i,link).LT.0)THEN
                      dur(i,link)=dur(i,link-1)
                      accepted(link)='n'
                    ENDIF    
                    sdur(i)=regul2(gibbscount)*sdur_ini(i)
                    gibbsonoff(gibbscount) = 1
                 ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                    gibbsonoff(gibbscount) = 0
                    dur(i,link) = dur(i,link-1)
                 ELSE
                    dur(i,link) = dur(i,link-1)
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
        DO i=1,npla
           IF(isjump(i,2).EQ.'n')THEN
             CALL gasdev_s(harvest)
             b(i,link) = b_ini(i) + sb_ini(i)*harvest
             IF(b(i,link).LT.0)THEN
                b(i,link)=b(i,link-1)
                accepted(link)='n'
             ENDIF    
           ELSE
              IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                 CALL gasdev_s(harvest)
                 b(i,link) = b(i,link-1) + regul*sb(i)*harvest
                 IF(b(i,link).LT.0)THEN
                   b(i,link)=b(i,link-1)
                   accepted(link)='n'
                 ENDIF  
              ELSE
                 gibbscount=gibbscount+1
                 IF(jumpgibbs(link).EQ.gibbscount)THEN
                    CALL gasdev_s(harvest)
                    b(i,link) = b(i,link-1) + regul2(gibbscount)*sb_ini(i)*harvest
                    IF(b(i,link).LT.0)THEN
                      b(i,link)=b(i,link-1)
                      accepted(link)='n'
                    ENDIF  
                    sb(i)=regul2(gibbscount)*sb_ini(i)
                    gibbsonoff(gibbscount) = 1
                 ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                    gibbsonoff(gibbscount) = 0 
                    b(i,link) = b(i,link-1)
                 ELSE
                    b(i,link) = b(i,link-1)
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
        DO i=1,npla
           IF(isjump(i,4).EQ.'n')THEN
              t0(i,link) = t0_ini(i) 
           ELSE
              IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                 CALL gasdev_s(harvest)
                 t0(i,link) = t0(i,link-1) + regul*st0(i)*harvest
                 dif = (t0(i,link)-t0_ini(i))
                 IF((dif/per(i,link)).GT.0.5) t0(i,link)=t0(i,link)-per(i,link)
                 IF((dif/per(i,link)).LT.-0.5) t0(i,link)=t0(i,link)+per(i,link)
              ELSE
                 gibbscount=gibbscount+1
                 IF(jumpgibbs(link).EQ.gibbscount)THEN
                    CALL gasdev_s(harvest)
                    t0(i,link) = t0(i,link-1) + regul2(gibbscount)*st0_ini(i)*harvest
                    st0(i)=regul2(gibbscount)*st0_ini(i)
                    gibbsonoff(gibbscount) = 1
                    dif = (t0(i,link)-t0_ini(i))
                    IF((dif/per(i,link)).GT.0.5) t0(i,link)=t0(i,link)-per(i,link)
                    IF((dif/per(i,link)).LT.-0.5) t0(i,link)=t0(i,link)+per(i,link)
                 ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                    gibbsonoff(gibbscount) = 0
                    t0(i,link) = t0(i,link-1)
                 ELSE
                    t0(i,link) = t0(i,link-1)
                 ENDIF
              ENDIF
           ENDIF
        ENDDO
        IF(ntiming.GT.0.AND.istiming.EQ.'y')THEN
           timerit(link) = 0.
           DO i=1,ntiming
              test2 = (timing(i)-t0(1,link))/per(1,link)
              epoch(i) = NINT(test2)
              ttiming = t0(1,link) + epoch(i)*per(1,link)
              omctiming(i) = timing(i) - ttiming
              timerit(link) = timerit(link) + (omctiming(i)/stiming(i))**2
           ENDDO
        ENDIF
        IF(limb.EQ.'qd'.AND.ntr.GT.0)THEN
           DO i=1,nfi
              IF(fitlimb(i).EQ.'n')THEN
                 ql(i,1,link) = ql_ini(i,1)
                 ql(i,2,link) = ql_ini(i,2)
              ELSE
                 IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                    CALL gasdev_s(harvest)
                    jumplimb(i,1,link) = jumplimb(i,1,link-1)+ &
                            & regul* djumplimb(i,1)*harvest
                 ELSE
                    gibbscount=gibbscount+1
                    IF(jumpgibbs(link).EQ.gibbscount)THEN
                       CALL gasdev_s(harvest)
                       jumplimb(i,1,link) = jumplimb(i,1,link-1)+ &
                            & regul2(gibbscount)*djumplimb_ini(i,1)*harvest
                       djumplimb(i,1)=regul2(gibbscount)*djumplimb_ini(i,1)
                       gibbsonoff(gibbscount) = 1
                    ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                       gibbsonoff(gibbscount) = 0
                       jumplimb(i,1,link) = jumplimb(i,1,link-1)
                    ELSE
                       jumplimb(i,1,link) = jumplimb(i,1,link-1)
                    ENDIF
                 ENDIF
                 IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                    CALL gasdev_s(harvest)
                    jumplimb(i,2,link) = jumplimb(i,2,link-1)+ &
                         & regul*djumplimb(i,2)*harvest
                 ELSE
                    gibbscount=gibbscount+1
                    IF(jumpgibbs(link).EQ.gibbscount)THEN
                       CALL gasdev_s(harvest)
                       jumplimb(i,2,link) = jumplimb(i,2,link-1)+ &
                            &regul2(gibbscount)* djumplimb_ini(i,2)*harvest
                       djumplimb(i,2)=regul2(gibbscount)*djumplimb_ini(i,2)
                       gibbsonoff(gibbscount) = 1
                    ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                       gibbsonoff(gibbscount) = 0
                       jumplimb(i,2,link) = jumplimb(i,2,link-1)
                    ELSE
                       jumplimb(i,2,link) = jumplimb(i,2,link-1)
                    ENDIF
                 ENDIF
                 IF((3*jumplimb(i,1,link)-jumplimb(i,2,link)).GE.5)THEN
                    jumplimb(i,1,link)=jumplimb(i,1,link-1)
                    jumplimb(i,2,link)=jumplimb(i,2,link-1)
                    accepted(link)='n'
                 ENDIF
                 ql(i,1,link) = (2*jumplimb(i,1,link)+jumplimb(i,2,link))/5.
                 ql(i,2,link) = (jumplimb(i,1,link)-2*jumplimb(i,2,link))/5.
              ENDIF
           ENDDO
        ENDIF
        IF(limb.EQ.'nl'.AND.ntr.GT.0)THEN
           DO i=1,nfi           
              nl(i,1,link) = nl_ini(i,1)
              nl(i,2,link) = nl_ini(i,2)
              nl(i,3,link) = nl_ini(i,3)
              nl(i,4,link) = nl_ini(i,4)
           ENDDO
        ENDIF
        IF(ntr.GT.0)THEN
           DO i=1,ntr
              IF(ramporder(i).GT.0.AND.rampmod.EQ.'exp')THEN
                 k=ramporder(i)
                 IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                    DO j=1,k
                      IF(j.EQ.1)THEN
                         CALL gasdev_s(harvest)
                         t1ramp(i,link)=t1ramp(i,link-1)+regul*et1ramp(i)*harvest
                         IF(t1ramp(i,link).LT.0)THEN
                           t1ramp(i,link)=t1ramp(i,link-1) 
                           accepted(link)='n'
                         ENDIF
                       ELSE IF(j.EQ.2)THEN
                         CALL gasdev_s(harvest)
                         t2ramp(i,link)=t2ramp(i,link-1)+regul*et2ramp(i)*harvest
                         IF(t2ramp(i,link).LT.0)THEN
                           t2ramp(i,link)=t2ramp(i,link-1) 
                           accepted(link)='n'
                         ENDIF
                      ENDIF
                   ENDDO
                 ELSE  
                    DO j=1,k
                      IF(j.EQ.1)THEN
                        gibbscount=gibbscount+1
                        IF(jumpgibbs(link).EQ.gibbscount)THEN
                           CALL gasdev_s(harvest)
                           t1ramp(i,link) = t1ramp(i,link-1)+regul2(gibbscount)*et1ramp_ini(i)*harvest
                           IF(t1ramp(i,link).LT.0)THEN
                              t1ramp(i,link)=t1ramp(i,link-1) 
                              accepted(link)='n'
                           ENDIF
                           et1ramp(i)=regul2(gibbscount)*et1ramp_ini(i)
                           gibbsonoff(gibbscount) = 1
                        ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                           gibbsonoff(gibbscount) = 0
                           t1ramp(i,link) = t1ramp(i,link-1)
                        ELSE
                           t1ramp(i,link) = t1ramp(i,link-1)
                        ENDIF
                      ELSE IF(j.EQ.2)THEN
                        gibbscount=gibbscount+1
                        IF(jumpgibbs(link).EQ.gibbscount)THEN
                           CALL gasdev_s(harvest)
                           t2ramp(i,link) = t2ramp(i,link-1)+regul2(gibbscount)*et2ramp_ini(i)*harvest
                           IF(t2ramp(i,link).LT.0)THEN
                              t2ramp(i,link)=t2ramp(i,link-1) 
                              accepted(link)='n'
                           ENDIF
                           et2ramp(i)=regul2(gibbscount)*et2ramp_ini(i)
                           gibbsonoff(gibbscount) = 1
                        ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                           gibbsonoff(gibbscount) = 0
                           t2ramp(i,link) = t2ramp(i,link-1)
                        ELSE
                           t2ramp(i,link) = t2ramp(i,link-1)
                        ENDIF
                      ENDIF
                    ENDDO
                 ENDIF
              ENDIF
              IF(sinusnumber(i).GT.0.AND.testsin.GT.0)THEN
                 k=sinusnumber(i)
                 IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                    DO j=1,k
                      CALL gasdev_s(harvest)
                      sip(i,j,link)=sip(i,j,link-1)+regul*esip(i,j)*harvest
                      IF(sip(i,j,link).LT.0)THEN
                        sip(i,j,link)=sip(i,j,link-1) 
                        accepted(link)='n'
                      ENDIF
                   ENDDO
                 ELSE  
                    DO j=1,k
                      gibbscount=gibbscount+1
                      IF(jumpgibbs(link).EQ.gibbscount)THEN
                         CALL gasdev_s(harvest)
                         sip(i,j,link) = sip(i,j,link-1)+regul2(gibbscount)*esip_ini(i,j)*harvest
                         IF(sip(i,j,link).LT.0)THEN
                           sip(i,j,link)=sip(i,j,link-1) 
                           accepted(link)='n'
                         ENDIF
                         esip(i,j)=regul2(gibbscount)*esip_ini(i,j)
                         gibbsonoff(gibbscount) = 1
                      ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                         gibbsonoff(gibbscount) = 0
                         sip(i,j,link) = sip(i,j,link-1)
                      ELSE
                         sip(i,j,link) = sip(i,j,link-1)
                      ENDIF
                    ENDDO
                 ENDIF
              ENDIF
              IF(sinusnumber(i).GT.0.AND.testsin.GT.0)THEN
                 k=sinusnumber(i)
                 IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                    DO j=1,k
                      test=0
                      CALL gasdev_s(harvest)
                      sit0(i,j,link)=sit0(i,j,link-1)+regul*esit0(i,j)*harvest
                    ENDDO
                 ELSE  
                    DO j=1,k
                      gibbscount=gibbscount+1
                      IF(jumpgibbs(link).EQ.gibbscount)THEN
                         CALL gasdev_s(harvest)
                         sit0(i,j,link) = sit0(i,j,link-1)+regul2(gibbscount)*esit0_ini(i,j)*harvest
                         esit0(i,j)=regul2(gibbscount)*esit0_ini(i,j)
                         gibbsonoff(gibbscount) = 1
                      ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                         gibbsonoff(gibbscount) = 0
                         sit0(i,j,link) = sit0(i,j,link-1)
                      ELSE
                         sit0(i,j,link) = sit0(i,j,link-1)
                      ENDIF
                    ENDDO
                 ENDIF
              ENDIF
              IF(flarenumber(i).GT.0.AND.testflare.GT.0)THEN
                 k=flarenumber(i)
                 IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                    DO j=1,k
                      CALL gasdev_s(harvest)
                      flt0(i,j,link)=flt0(i,j,link-1)+regul*eflt0(i,j)*harvest
                      IF(flt0(i,j,link).LT.0)THEN
                        flt0(i,j,link)=flt0(i,j,link-1) 
                        accepted(link)='n'
                      ENDIF
                   ENDDO
                 ELSE  
                    DO j=1,k
                      gibbscount=gibbscount+1
                      IF(jumpgibbs(link).EQ.gibbscount)THEN
                         CALL gasdev_s(harvest)
                         flt0(i,j,link) = flt0(i,j,link-1)+regul2(gibbscount)*eflt0_ini(i,j)*harvest
                         IF(flt0(i,j,link).LT.0)THEN
                           flt0(i,j,link)=flt0(i,j,link-1) 
                           accepted(link)='n'
                         ENDIF
                         eflt0(i,j)=regul2(gibbscount)*eflt0_ini(i,j)
                         gibbsonoff(gibbscount) = 1
                      ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                         gibbsonoff(gibbscount) = 0
                         flt0(i,j,link) = flt0(i,j,link-1)
                      ELSE
                         flt0(i,j,link) = flt0(i,j,link-1)
                      ENDIF
                    ENDDO
                 ENDIF
              ENDIF
              IF(flarenumber(i).GT.0.AND.testflare.GT.0)THEN
                 k=flarenumber(i)
                 IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                    DO j=1,k
                      CALL gasdev_s(harvest)
                      flampli(i,j,link)=flampli(i,j,link-1)+regul*eflampli(i,j)*harvest
                      IF(flampli(i,j,link).LT.0)THEN
                        flampli(i,j,link)=flampli(i,j,link-1) 
                        accepted(link)='n'
                      ENDIF
                   ENDDO
                 ELSE  
                    DO j=1,k
                      gibbscount=gibbscount+1
                      IF(jumpgibbs(link).EQ.gibbscount)THEN
                         CALL gasdev_s(harvest)
                         flampli(i,j,link) = flampli(i,j,link-1)+regul2(gibbscount)*eflampli_ini(i,j)*harvest
                         IF(flampli(i,j,link).LT.0)THEN
                           flampli(i,j,link)=flampli(i,j,link-1) 
                           accepted(link)='n'
                         ENDIF
                         eflampli(i,j)=regul2(gibbscount)*eflampli_ini(i,j)
                         gibbsonoff(gibbscount) = 1
                      ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                         gibbsonoff(gibbscount) = 0
                         flampli(i,j,link) = flampli(i,j,link-1)
                      ELSE
                         flampli(i,j,link) = flampli(i,j,link-1)
                      ENDIF
                    ENDDO
                 ENDIF
              ENDIF
              IF(flarenumber(i).GT.0.AND.testflare.GT.0)THEN
                 k=flarenumber(i)
                 IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
                    DO j=1,k
                      test=0
                      CALL gasdev_s(harvest)
                      fltau(i,j,link)=fltau(i,j,link-1)+regul*efltau(i,j)*harvest
                    ENDDO
                 ELSE  
                    DO j=1,k
                      gibbscount=gibbscount+1
                      IF(jumpgibbs(link).EQ.gibbscount)THEN
                         CALL gasdev_s(harvest)
                         fltau(i,j,link) = fltau(i,j,link-1)+regul2(gibbscount)*efltau_ini(i,j)*harvest
                         efltau(i,j)=regul2(gibbscount)*efltau_ini(i,j)
                         gibbsonoff(gibbscount) = 1
                      ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                         gibbsonoff(gibbscount) = 0
                         fltau(i,j,link) = fltau(i,j,link-1)
                      ELSE
                         fltau(i,j,link) = fltau(i,j,link-1)
                      ENDIF
                    ENDDO
                 ENDIF
              ENDIF
         ENDDO
      ENDIF
      IF(nttvmax.GT.0.AND.isttv.EQ.'y')THEN
        DO i=1,npla
          DO j=1,nttv(i)
            IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
              CALL gasdev_s(harvest)
              ttv(i,j,link)=ttv(i,j,link-1)+regul*sttv(i,j)*harvest
            ELSE  
              gibbscount=gibbscount+1
              IF(jumpgibbs(link).EQ.gibbscount)THEN
                CALL gasdev_s(harvest)
                ttv(i,j,link) = ttv(i,j,link-1)+regul2(gibbscount)*sttv_ini(i,j)*harvest
                sttv(i,j)=regul2(gibbscount)*sttv_ini(i,j)
                gibbsonoff(gibbscount) = 1
              ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                gibbsonoff(gibbscount) = 0
                ttv(i,j,link) = ttv(i,j,link-1)
              ELSE
                ttv(i,j,link) = ttv(i,j,link-1)
              ENDIF
            ENDIF
            ttr(i,j,link)=t0(i,link)+per(i,link)*epochtr(i,j)+ttv(i,j,link) 
          ENDDO
        ENDDO
      ENDIF
    ENDIF
    IF(ntr.GT.0)THEN
      DO i=1,nfi
        CALL gasdev_s(harvest)
        dilution(i) = 1+0.01*(dilu(i)+harvest*edilu(i))
      ENDDO
    ENDIF
    IF(testf2.EQ.'y')THEN
       IF(priorf2.EQ.'n')THEN
          test=0
          DO WHILE(test.LT.1)
             CALL gasdev_s(harvest)
             f2(link) = spec(7)+harvest*sspec(7) 
             IF(f2(link).GT.0) test = 2
          ENDDO
       ELSE
          IF(gibbs_sampler.EQ.'n'.OR.MOD(link,nat).GT.nburn2)THEN
             CALL gasdev_s(harvest)
             f2(link) = f2(link-1) + regul*sf2*harvest
             IF(f2(link).LT.0)THEN
               f2(link)=f2(link-1)
               accepted(link)='n'
             ENDIF
          ELSE
             gibbscount=gibbscount+1
             IF(jumpgibbs(link).EQ.gibbscount)THEN
                CALL gasdev_s(harvest)
                f2(link) = f2(link-1) + regul2(gibbscount)*sspec(7)*harvest
                IF(f2(link).LT.0)THEN
                  f2(link)=f2(link-1)
                  accepted(link)='n'
                ENDIF   
                sf2=regul2(gibbscount)*sspec(7)
                gibbsonoff(gibbscount) = 1     
             ELSE IF(gibbsonoff(gibbscount).EQ.1)THEN
                gibbsonoff(gibbscount) = 0
                f2(link) = f2(link-1)
             ELSE
                f2(link) = f2(link-1)
             ENDIF
          ENDIF
       ENDIF
    ENDIF
    
    ntra=0
    dar=0

! MODEL COMPUTATION 
    DO i=1,npla                
      IF(exc(i,link).EQ.0.) THEN
        CALL gasdev_s(harvest)
        if (link.eq.1) then
        	omega(i,link) = omega(i,1) + harvest*er_om(i) !link-
        else
        	omega(i,link) = omega(i,link-1) + harvest*er_om(i)
        end if
      ELSE
        IF(ABS(secosw(i,link)).GT.1.E-13)THEN
          omega(i,link) = DATAN(sesinw(i,link)/secosw(i,link))
        ELSE
          omega(i,link) = DASIN(sesinw(i,link)/SQRT(exc(i,link)))
        ENDIF 
        IF(omega(i,link).LT.0) omega(i,link)=omega(i,link)+2*pi
        IF(secosw(i,link).LT.0)THEN
          IF(sesinw(i,link).LT.0)THEN
            omega(i,link)=omega(i,link)+pi
          ELSE
            omega(i,link)=omega(i,link)-pi
          ENDIF
        ENDIF
        omega(i,link)=omega(i,link)*180/pi
      ENDIF
      ecosw(i,link)=exc(i,link)*DCOS(omega(i,link)*pi/180)
      esinw(i,link)=exc(i,link)*DSIN(omega(i,link)*pi/180)   
      Tano_tr(i) = (pi/2.) - omega(i,link)*pi/180
      Eano_tr(i) = 2.* DATAN(SQRT((1-exc(i,link))/(1+exc(i,link)))*DTAN(Tano_tr(i)/2.))
      rr(i,link) = SQRT(dF(i,link)) 
      temp3(i) = b(i,link)*(1-exc(i,link)**2)/(1.+exc(i,link)*DSIN(omega(i,link)*pi/180))
      
      rade=SQRT(1-exc(i,link)**2)
      IF(m2op.EQ.'y')THEN  
        m2 = (mass_s(link)*sunmass)**(2./3.)*kb(i,link)   
        m2 = m2/(2.*pi*gravi)**(1./3.)                    
        m2 = m2/(jupmass*(24*3600.)**(-1./3.))           
        m2_m1(i) = (m2*jupmass)/(mass_s(link)*sunmass)  
        m2 = m2*((1+m2_m1(i))**(2./3.))
        m2_m1(i) = (m2*jupmass)/(mass_s(link)*sunmass)  
        m2_m1(i) = (1+m2_m1(i))**(1./3.)
        kePM=ka(i,link)*rade*((per(i,link)*24.*3600.)/(2.*pi*gravi*mass_s(link)*sunmass))**(1./3.)
        f02=1./3.*m2_m1(i)**(-2./3.)*(kePM+5./3.*kePM**2)/(1.+kePM)**(1./3.)
        Dm2_m1=f02*sqrt((dkb(i)/ka(i,link))**2+(sper(i)/per(i,link))**2+(dmass_s/mass_s(link))**2+ &
        	& (exc(i,link)*er_exc(i)/rade**2)**2)
      ELSE                                             
        m2_m1(i) = 1.
        f02=0.
        Dm2_m1=0.
      ENDIF
 
      IF(i.EQ.1.AND.(fitmsrs.EQ.'y'.OR.fixstellar.EQ.'y'))THEN
        rhotep = mass_s(link)/(radius_s(link)**3)
        rho(link)=rhotep
      ENDIF
      IF(fitmsrs.EQ.'y'.OR.fixstellar.EQ.'y')THEN   
        a_R(i,link) = (((rhotep)/0.0134235)*per(i,link)**2)**(1./3.) 
        a_R(i,link) = a_R(i,link)*m2_m1(i)
        IF(b(i,link).GE.a_R(i,link))THEN
          b(i,link)=b(i,link-1)
          IF(b(i,link).GE.a_R(i,link)) b(i,link)=b(i,link)-2*(b(i,link)-a_R(i,link))
          accepted(link)='n'
        ENDIF
      ENDIF
      IF(fitmsrs.EQ.'n'.AND.fixstellar.EQ.'n')THEN
        IF(i.EQ.1) then
         rho(link)=0.
         Drho=0.
        end if
        IF(((1+rr(i,link))**2 - temp3(i)**2).GT.0.AND.RjumpIso.eq.'n'.and.(isoch.eq.'y'.or.d.eq.1))THEN
          if (.not.isEq(dur(i,link),0.D0,5)) then !!dur(i,link).GT.0.
          ! Approximation in Andrew's powerpoint 
              dar=dar+1
              a_R(i,link) = (per(i,link)/dur(i,link))*(SQRT((1+rr(i,link))**2 - temp3(i)**2)/pi)* &
              & (SQRT(1-exc(i,link)**2))/(1.+exc(i,link)*DSIN(omega(i,link)*pi/180))
              f01=per(i,link)*SQRT(1-exc(i,link)**2)/(2*pi*dur(i,link)*(1.+exc(i,link)*DSIN(omega(i,link)*pi/180)))
              radq=SQRT((1+rr(i,link))**2 - temp3(i)**2)
              dcon=1.+exc(i,link)*DSIN(omega(i,link)*pi/180)
              daRdP=a_R(i,link)/per(i,link)
              daRdW=a_R(i,link)/dur(i,link)
              daRddF=f01*(1+rr(i,link))/(rr(i,link)*radq)
              daRdb=f01*2.*b(i,link)*(rade**2/dcon)**2/radq !ok
              !temp3(i)/b(i,link) could lead to 0/0
              daRde=per(i,link)/(pi*dur(i,link))*abs(-b(i,link)**2*(2*exc(i,link)+DSIN(omega(i,link)*pi/180)* &
              & (exc(i,link)**2+1))/radq*rade**3/dcon**4 - radq*(exc(i,link)+DSIN(omega(i,link)*pi/180))/ &
              & (dcon**2*rade)) !ok
              daRdom=per(i,link)/(pi*dur(i,link))*rade*exc(i,link)*abs(DCOS(omega(i,link)*pi/180)/dcon**2* &
              & (b(i,link)**2*rade**4/(dcon**2*radq)-radq)) !ok
              
              Da_R=SQRT((daRdP*sper(i))**2+(daRdW*sdur(i))**2+(daRddF*sdF(i))**2+(daRdb*sb(i))**2+ &
              & (daRde*er_exc(i))**2+(daRdom*er_om(i))**2) !ok
          ! Actual formula (Seager & Mallen-Ornelas 2002)
          ! a_R(i,link) = SQRT(((1+rr(i,link))**2 - (1-(DSIN(dur(i,link)*pi/per(i,link)))**2) & 
          ! & *temp3(i)**2)/DSIN(dur(i,link)*pi/per(i,link))**2) * &
          ! & (SQRT(1-exc(i,link)**2))/(1.+exc(i,link)*DSIN(omega(i,link)*pi/180))
          else
          	  a_R(i,link)=0.
              Da_R=0.
          endif
        else if (RjumpIso.eq.'y') then !(shallow transit|posterior_e)!isoch.eq.'y'.and.isEq(dur(i,link),0.D0,5)
          !elseifStatem implies logg/sismo/R/L
	      if (priorteff.eq.'p') then
	        TeffI=temp_s(link)
	        I_TeffI=dtemp_s
	        colI=col
	        I_colI=I_col
	      else !priorcol. TeffI used just in selectMass
	        TeffI=10.**polyval(cJ(:,idCol),col_s(link))
	        !print*,'TeffI-polyval',TeffI
	        !(3.908-0.234*col_s(link)) !logTe-BmV relation according to Johnson (1966)
	        I_TeffI=0.02*log(10.)*TeffI
	        colI=col_s(link)
	        I_colI=dcol_s
	      end if
	      if (priorvsini.ne.'n') then
		    vsiniI=vsini(link)
		    if (vsiniI.lt.0) then
		    	vsiniI=spec(5)
		    end if
			I_vsiniI=dvsini
		  else
			vsiniI=-1.
			I_vsiniI=-1.
		  end if
	      if (link.eq.1.and.ntra.lt.1) then !!load isoch
	        if (priorrho.ne.'n') then
			  rhoI=rhoInp*rhoSun !conversion into g/cm3
			  DrhoI=I_rhoInp*rhoSun
			else
			  rhoI=-1.
			  DrhoI=-1.
			end if
			
	        ltitle="#row FeH I_FeH Teff I_Teff rho I_rho logg I_logg vsini I_vsini Prot I_Prot logRHK YMg"// &
					& " numax I_numax Dnu I_Dnu sismo Rinput I_Rinput mag color dist I_dist Hipf useCol idCol"
	        call system('mkdir -p Isoch')
			OPEN(UNIT=3,FILE='./Isoch/Star.txt')
			!write parameters to be used by SCP just once to check
			write(3,*) ltitle
			write(3,*) "1 ",met_s(link),dmet_s,TeffI,I_TeffI,rhoI,DrhoI,loggI,I_loggI,vsiniI,I_vsiniI, &
						& protI,I_protI,logRHKI,YMgI,numaxI,I_numaxI,deltanuI,I_deltanuI,sismo,Rinp,I_Rinp, &
						& mag,colI,dist,I_dist,Hipf,useColor,idCol
			close(3)
									    	
			OPEN(UNIT=4,FILE='/home/bonfanti/Documents/TopCat/Z_iso.txt',STATUS='OLD')
			iz=0
			READ(4, '(A)', IOSTAT=iol) head
			do while (iol.eq.0)
				iz=iz+1
				READ(4, '(A)', IOSTAT=iol) head
			end do
			print*, 'dimZ_iso',iz
			dimZ_iso=iz
			close(4)
			allocate(Z_iso(dimZ_iso))
			OPEN(UNIT=4,FILE='/home/bonfanti/Documents/TopCat/Z_iso.txt',STATUS='OLD')
			do iz=1,dimZ_iso
				read(4,*) Z_iso(iz)
			end do
			close(4)
			if (.not.isEq(Z_iso(1),Zi_inf,5).or..not.isEq(Z_iso(size(Z_iso)),Zi_sup,3)) then
				print*,'Specified Zi_inf or Zi_sup inconsistent with available Z values'
				print*,'Check either Ziso.txt or IsoPD.f90'
				stop
			end if
			
			Mstep=0.01
			vrho=(/ -1.D0,deltanuI,rhoInp /)
			I_vrho=(/ -1.D0,I_deltanuI,I_rhoInp /)
			vg=(/ loggI,numaxI /)
			I_vg=(/ I_loggI,I_numaxI /)
			if (priormag.ne.'n'.and.priordist.ne.'n') then
				Mabs=mag-5.*log10(dist)+5.
			else
				Mabs=-100.
			end if
			call selectMass(vrho,vg,I_vrho,I_vg,TeffI,I_TeffI,met_s(link),dmet_s,Rinp,I_Rinp,Mabs,idCol,Minf,Msup)
			
			if (priormet.eq.'p') then    		
				Zinf=10.**(spec(3)-3*sspec(3)-costZ)
				Zsup=10.**(spec(3)+3*sspec(3)-costZ)
				xZi=minloc(abs(Zinf-Z_iso),1)
				xZu=minloc(abs(Zsup-Z_iso),1)
			else !priormet.eq.'y' -> load every isoch
				xZi=1
				xZu=size(Z_iso)
				Zinf=Z_iso(xZi)
				Zsup=Z_iso(xZu)
			end if
			allocate(Zvec(xZu-xZi+1))
			Zvec=Z_iso(xZi:xZu)
						
			allocate(nM(size(Zvec)))
			allocate(fnd05(size(Zvec)))
			allocate(Zt(size(Zvec)))
			!for each metallicity, evaluate number of tracks (nM) to be loaded
			call tracksNumber(Minf,Msup,Mstep,MlowMS,Zvec,fnd05,Zt,nM)
						
			call uniqueFast(Zt,4,iZt,.true.)
			allocate(Ztvec(size(iZt)))
			Ztvec=Zt(iZt)
			allocate(indxZt(size(Zvec)))
			!nth position of indxZt reports the index of the tracks whose metallicity is the closest to the
			!one of the isochrone that has been opened (the opened isochrone is the nth in terms of Z)
			do jj=1,size(iZt)-1
				indxZt(iZt(jj):iZt(jj+1)-1)=(/ (jj, idum=1,iZt(jj+1)-iZt(jj)) /)
			end do
			indxZt(iZt(size(iZt)):size(indxZt))=(/ (size(iZt), idum=1,size(indxZt)-iZt(size(iZt))+1) /)
			
			allocate(nMt(size(Ztvec)))
			nMt=nM(iZt)
			deallocate(nM)
			
			allocate(fnd05t(size(Ztvec)))
			fnd05t=fnd05(iZt)
			deallocate(fnd05)
			
			print*,'Storing isochrones and tracks...'
			!store relevant isochrones in the big matrix IsocTab
			call cpu_time(start)
			call storeIsoc(Zvec,IsocTab,Zndxi,Zndxf,idCol)
			call cpu_time(finish)
			write(*,'(A20,F5.2,A2)') 'Isochrones storage ',finish-start,' s'
			
			call cpu_time(start)
			call storeTracks_V_ZAMS(Minf,Msup,Mstep,MlowMS,Ztvec,nMt,fnd05t,Mavail,TrackTab,Tndxi,Tndxf, &
					& velTrackTab,Vndxi,Vndxf,ZAMStab,ZAndxi,ZAndxf)
			call cpu_time(finish)
			write(*,'(A16,F5.2,A2)') 'Tracks storage ',finish-start,' s'				
						
			call storeEvoZ(Minf,Msup,ZtevoTab,Endxi,Endxf,MeAv)
			
			Barnes="/home/bonfanti/Documents/Dottorato/Articoli/BrownSWPAgeGiroIsoc/tauBarnesKim2010.txt"
			call loadMatrix(Barnes,GyroTab,head)
			print*,'Storing ended'
!			!Printing test
!			print*,'First 3 Z',Zvec(1:3)
!			print*,'First 3 indices init',Zndxi(1:3)
!			print*,'First 3 indices fin',Zndxf(1:3)
!			print*,'1st row, 3rd isoc',IsocTab(Zndxi(3),:)
!			print*,'last row, 3rd isoc',IsocTab(Zndxf(3),:)
!			print*,'------'

	      end if
	      if (ntra.lt.1) then
	      	  ntra=ntra+1
	      	  if (link.eq.1.and.priormag.eq.'p'.and.priordist.eq.'p') then !just to obtain L=L(d,v)
	      	  	if (priorteff.eq.'p') then
	      	  		TeffL=spec(1)
	      	  		I_TeffL=sspec(1)
	      	  		colL=-100.
	      	  	else
	      	  		TeffL=0.
	      	  		I_TeffL=0.
	      	  		colL=col
	      	  	end if
	      	  	star= (/ dble(1),spec(3),sspec(3),TeffL,I_TeffL,rhoI,DrhoI,loggI,I_loggI,vsiniI,I_vsiniI, &
					& protI,I_protI,logRHKI,YMgI,numaxI,I_numaxI,deltanuI,I_deltanuI,dble(sismo),-1.D0,-1.D0, &
					& mag,colL,dist,I_dist,Hipf,useColor,dble(idCol) /)
				call SCPmcmcPD12S(star,ltitle,IsocTab,Zvec,Zndxi,Zndxf,TrackTab,nMt,Mavail,indxZt,Tndxi,Tndxf,velTrackTab, &
						& Vndxi,Vndxf,GyroTab,ZAMStab,ZAndxi,ZAndxf,ZtevoTab,Endxi,Endxf,MeAv,massTmp,dmassTmp, &
						& radiusTmp,dradiusTmp,temp_s(link),dtemp_s,age(link),dage,starlum,dstarlum,0,row,acc)
				print*,'calibL',starlum,dstarlum
	      	  end if
	      	  
	      	  star= (/ dble(1),met_s(link),dmet_s,TeffI,I_TeffI,rhoI,DrhoI,loggI,I_loggI,vsiniI,I_vsiniI, &
					& protI,I_protI,logRHKI,YMgI,numaxI,I_numaxI,deltanuI,I_deltanuI,dble(sismo),Rinp,I_Rinp, &
					& mag,colI,dist,I_dist,Hipf,useColor,dble(idCol) /)
								
			  call SCPmcmcPD12S(star,ltitle,IsocTab,Zvec,Zndxi,Zndxf,TrackTab,nMt,Mavail,indxZt,Tndxi,Tndxf,velTrackTab, &
							& Vndxi,Vndxf,GyroTab,ZAMStab,ZAndxi,ZAndxf,ZtevoTab,Endxi,Endxf,MeAv,mass_Isoch,dmass_Isoch, &
							& radius_Isoch,dradius_Isoch,temp_s(link),dtemp_s,age(link),dage,Lum,I_Lum,link,row,acc)
	          age(link)=age(link)/1.E9 !from yrs to Gyr
	          dage=dage/1.E9
	          if (acc.eq.0) then !Algorithm hasn't converged
				print*,'No convergence reached in theoretical models.'
				print*,'chain step: ',link,' chain ',(link-1)/nat+1
				accepted(link)='n'
			  else
			  	rowTot=rowTot+row
			  	rowIO(link)=row
			  end if
			  rhotepi=mass_s(link)/radius_s(link)**3
			  Drhotepi=rhotepi*sqrt((dmass_s/mass_s(link))**2+(3*dstarradius/radius_s(link))**2) !sigma !ok
		  end if
		  a_R(i,link)=(((rhotepi)/0.0134235)*per(i,link)**2)**(1./3.)
		  f03=a_R(i,link)
	      a_R(i,link) = a_R(i,link)*m2_m1(i)
	      Drho1=Drhotepi
	      !rho is not inferred from a_R. It's the opposite
	    ELSE
          dar=dar+1
          t=0
          DO WHILE(t.LT.1)
            CALL gasdev_s(harvest)
            rhotep= trho+harvest*strho
            IF(rhotep.GT.0)THEN
              t=2
              a_R(i,link) = (((rhotep)/0.0134235)*per(i,link)**2)**(1./3.)
              f03=a_R(i,link)
              a_R(i,link) = a_R(i,link)*m2_m1(i)
              Drho1=abs(harvest*strho)
              !rho not inferred from a_R
            ENDIF
          ENDDO
        ENDIF
        rho(link)=rho(link)+(0.0134235)*((a_R(i,link)/m2_m1(i))**3.)*(1./per(i,link)**2.)
		if (((1+rr(i,link))**2 - temp3(i)**2).GT.0.AND.RjumpIso.eq.'n') then !.not.isEq(dur(i,link),0.D0,5)
		    if (.not.isEq(dur(i,link),0.D0,5)) then
		    	Drho=Drho+rho(link)*sqrt((3.*Da_R/a_R(i,link))**2+(3.*Dm2_m1/m2_m1(i))**2+(2.*sper(i)/per(i,link))**2) !ok
		    endif
		else	
			Drho=Drho+Drho1
		end if
        IF(i.EQ.npla) then
         if (RjumpIso.eq.'n') then
           rho(link)=rho(link)/DBLE(dar)
           Drho=Drho/DBLE(dar)
         else
           rho(link)=rho(link)/DBLE(npla)
           Drho=Drho/DBLE(npla)
         endif
        end if
      ENDIF
    ENDDO
    DO i=1,npla
      a_R(i,link) = (((rho(link))/0.0134235)*per(i,link)**2)**(1./3.) 
      a_R(i,link) = a_R(i,link)*m2_m1(i)
      if(.not.(isoch.eq.'y'.and.ntr.eq.0.and.nrv.eq.0))then !not in case of the IsochPlacement only
		  IF((a_R(i,link)*(1-exc(i,link))).LE.(1+rr(i,link)))THEN
	  !       PRINT*, 'Boum!'
		     accepted(link)='n'
		  ENDIF
	  end if
!      if (mod(link,100).eq.0) print*,'Dur condition',((1+rr(i,link))**2 - temp3(i)**2)
      IF(((1+rr(i,link))**2 - temp3(i)**2).GT.0.)THEN
        ! Approximation in Andrew's powerpoint
          dur(i,link) = (per(i,link)/a_R(i,link))*(SQRT((1+rr(i,link))**2 - temp3(i)**2)/pi)* &
          & (SQRT(1-exc(i,link)**2))/(1.+exc(i,link)*DSIN(omega(i,link)*pi/180))
        ! Actual formula (Seager & Mallen-Ornelas 2002)
        !  dur(i,link)=(per(i,link)/pi)*ASIN((1./a_R(i,link))*SQRT(((1+rr(i,link))**2 - & 
        !  & temp3(i)**2)/(1. - (b(i,link)/a_R(i,link))**2))) *    & 
        !  & (SQRT(1-exc(i,link)**2))/(1.+exc(i,link)*DSIN(omega(i,link)*pi/180))
      ELSE
        dur(i,link)=0.
      ENDIF
      IF(i.EQ.1)THEN
        if (priorteff.eq.'p') teff = temp_s(link)
        feh = met_s(link)
        IF(enoch.EQ.'y')THEN
          order=7
          ALLOCATE(po(nenoch2,3));ALLOCATE(y(nenoch2));ALLOCATE(sig(nenoch2))
          ALLOCATE(w(order));ALLOCATE(v(order,order));ALLOCATE(u(nenoch2,order))
          ALLOCATE(dummass(nenoch2))
          DO k=1,nenoch2
            CALL gasdev_s(harvest)
            po(k,1) = enocht(k,1) + harvest*enocht(k,2)
            po(k,1) = LOG10(po(k,1))-4.1
            CALL gasdev_s(harvest)
            dumnum1 = enochm(k,1) + harvest*enochm(k,2)
            dummass(k) = dumnum1
            CALL gasdev_s(harvest)
            dumnum2 = enochr(k,1) + harvest*enochr(k,2)
            po(k,2) = LOG10(dumnum1/dumnum2**3)
            CALL gasdev_s(harvest)
            po(k,3) = enochf(k,1) + harvest*enochf(k,2)
            y(k)=LOG10(dumnum1)
            sig(k)=0.5*(LOG10(dumnum1+enochm(k,2))-LOG10(dumnum1-enochm(k,2)))
          ENDDO
          CALL ebfit(po,y,sig,nenoch2,enochcoef,order,u,v,w,nenoch2,order,chisq,enochlaw)
          rms = 0.
          massmeanerror = 0.
          DO k=1,nenoch2
            dub = enochcoef(1)+enochcoef(2)*po(k,1)+enochcoef(3)*po(k,1)**2 + &
             & enochcoef(4)*po(k,2) + enochcoef(5)*po(k,2)**2 + enochcoef(6)*po(k,2)**3 + &
             & enochcoef(7)*po(k,3)
            rms = rms + (10**(dub)-dummass(k))**2
            massmeanerror=massmeanerror+enochm(k,2)/DBLE(nenoch2)
          ENDDO
          rms = SQRT(rms/nenoch2)
          IF(rms.GE.massmeanerror)THEN
            massjitter(link) = SQRT(rms**2 - massmeanerror**2)
          ELSE
            massjitter(link) = 0.
          ENDIF
          CALL enochfor(teff,feh,rho(link),enochcoef,massenoch)
          CALL gasdev_s(harvest)
          mass_s(link)=massenoch + harvest*massjitter(link)
          radius_s(link)=(mass_s(link)/rho(link))**(1./3.)
          DEALLOCATE(po);DEALLOCATE(y);DEALLOCATE(sig);DEALLOCATE(w);DEALLOCATE(v);DEALLOCATE(u)
          DEALLOCATE(dummass)
        ENDIF
        if (isoch.eq.'y'.and.ntra.lt.1) then !.not.isEq(dur_ini(i),0.D0,5)
        	ntra=ntra+1 !but shouldn't be necessary
        	if (priorteff.eq.'p') then
		        TeffI=temp_s(link)
		        I_TeffI=dtemp_s
		        colI=col
		        I_colI=I_col
		    else !priorcol. TeffI used just in selectMass
		        TeffI=10.**polyval(cJ(:,idCol),col_s(link))
		        !print*,'TeffI-polyval',TeffI
		        !(3.908-0.234*col_s(link)) !logTe-BmV relation according to Johnson (1966)
		        I_TeffI=0.02*log(10.)*TeffI
		        colI=col_s(link)
		        I_colI=dcol_s
		    end if
        	if (priorvsini.ne.'n') then
				vsiniI=vsini(link)
				if (vsiniI.lt.0) then
			    	vsiniI=spec(5)
			    end if
				I_vsiniI=dvsini
	    	else
	    		vsiniI=-1.
	    		I_vsiniI=-1.
	    	end if
	    	if (priorrho.ne.'n') then
	    		rhoI=(rho(link)/Drho**2+rhoInp/I_rhoInp**2)/(1./Drho**2+1./I_rhoInp**2)*rhoSun
	    		DrhoI=1./sqrt(1./Drho**2+1./I_rhoInp**2)*rhoSun
	    	else
				rhoI=rho(link)*rhoSun !conversion into g/cm3
				DrhoI=Drho*rhoSun
			end if
	    	if (link.eq.1) then
	    		ltitle="#row FeH I_FeH Teff I_Teff rho I_rho logg I_logg vsini I_vsini Prot I_Prot logRHK YMg"// &
	    			& " numax I_numax Dnu I_Dnu sismo Rinput I_Rinput mag color dist I_dist Hipf useCol idCol"
	    		call system('mkdir -p Isoch')
	    		OPEN(UNIT=3,FILE='./Isoch/Star.txt')
	    		!write parameters to be used by SCP just once to check
				write(3,*) ltitle
				write(3,*) "1 ",met_s(link),dmet_s,TeffI,I_TeffI,rhoI,DrhoI,loggI,I_loggI,vsiniI,I_vsiniI, &
							& protI,I_protI,logRHKI,YMgI,numaxI,I_numaxI,deltanuI,I_deltanuI,sismo,Rinp,I_Rinp, &
							& mag,colI,dist,I_dist,Hipf,useColor,idCol
				close(3)
						    			    	
	    		OPEN(UNIT=4,FILE='/home/bonfanti/Documents/TopCat/Z_iso.txt',STATUS='OLD')
	    		iz=0
				READ(4, '(A)', IOSTAT=iol) head
				do while (iol.eq.0)
					iz=iz+1
					READ(4, '(A)', IOSTAT=iol) head
				end do
				print*, 'dimZ_iso',iz
				dimZ_iso=iz
				close(4)
				allocate(Z_iso(dimZ_iso))
				OPEN(UNIT=4,FILE='/home/bonfanti/Documents/TopCat/Z_iso.txt',STATUS='OLD')
	    		do iz=1,dimZ_iso
	    			read(4,*) Z_iso(iz)
	    		end do
	    		close(4)
	    		if (.not.isEq(Z_iso(1),Zi_inf,5).or..not.isEq(Z_iso(size(Z_iso)),Zi_sup,3)) then
					print*,'Specified Zi_inf or Zi_sup inconsistent with available Z values'
					stop
				end if		
	    		
	    		Mstep=0.01
				vrho=(/ rho(link),deltanuI,rhoInp /)
				I_vrho=(/ 3.*Drho,I_deltanuI,I_rhoInp /)
				vg=(/ loggI,numaxI /)
				I_vg=(/ I_loggI,I_numaxI /)
				if (priormag.ne.'n'.and.priordist.ne.'n') then
					Mabs=mag-5.*log10(dist)+5.
				else
					Mabs=-100.
				end if
				call selectMass(vrho,vg,I_vrho,I_vg,TeffI,I_TeffI,met_s(link),dmet_s,Rinp,I_Rinp,Mabs,idCol,Minf,Msup)
	    		
				if (priormet.eq.'p') then
					Zinf=10.**(spec(3)-3.*sspec(3)-costZ)
					Zsup=10.**(spec(3)+3.*sspec(3)-costZ)
					xZi=minloc(abs(Zinf-Z_iso),1)
					xZu=minloc(abs(Zsup-Z_iso),1)
				else !priormet.eq.'y' -> load every isoch
					xZi=1
					xZu=size(Z_iso)
					Zinf=Z_iso(xZi)
					Zsup=Z_iso(xZu)
				end if
				allocate(Zvec(xZu-xZi+1))
				Zvec=Z_iso(xZi:xZu)
				
				allocate(nM(size(Zvec)))
				allocate(fnd05(size(Zvec)))
				allocate(Zt(size(Zvec)))
				!for each metallicity, evaluate number of tracks (nM) to be loaded
				call tracksNumber(Minf,Msup,Mstep,MlowMS,Zvec,fnd05,Zt,nM)
								
				call uniqueFast(Zt,4,iZt,.true.)
				allocate(Ztvec(size(iZt)))
				Ztvec=Zt(iZt)
				allocate(indxZt(size(Zvec)))
				!nth position of indxZt reports the index of the tracks whose metallicity is the closest to the
				!one of the isochrone that has been opened (the opened isochrone is the nth in terms of Z)
				do jj=1,size(iZt)-1
					indxZt(iZt(jj):iZt(jj+1)-1)=(/ (jj, idum=1,iZt(jj+1)-iZt(jj)) /)
				end do
				indxZt(iZt(size(iZt)):size(indxZt))=(/ (size(iZt), idum=1,size(indxZt)-iZt(size(iZt))+1) /)
				
				allocate(nMt(size(Ztvec)))
				nMt=nM(iZt)
				deallocate(nM)
				
				allocate(fnd05t(size(Ztvec)))
				fnd05t=fnd05(iZt)
				deallocate(fnd05)
				
				print*,'Storing isochrones and tracks...'
				!store relevant isochrones in the big matrix IsocTab
				call cpu_time(start)
				call storeIsoc(Zvec,IsocTab,Zndxi,Zndxf,idCol)
				call cpu_time(finish)
				write(*,'(A20,F5.2,A2)') 'Isochrones storage ',finish-start,' s'
				
				call cpu_time(start)
!				print*,'Minf',Minf
!				print*,'Msup',Msup
!				print*,'Mstep',Mstep
!				print*,'MlowMS',MlowMS
!				print*,'Ztvec',Ztvec
!				print*,'nMt',nMt
!				print*,'fnd05t',fnd05t
				call storeTracks_V_ZAMS(Minf,Msup,Mstep,MlowMS,Ztvec,nMt,fnd05t,Mavail,TrackTab,Tndxi,Tndxf, &
					& velTrackTab,Vndxi,Vndxf,ZAMStab,ZAndxi,ZAndxf)
				call cpu_time(finish)
				write(*,'(A16,F5.2,A2)') 'Tracks storage ',finish-start,' s'
				
				call storeEvoZ(Minf,Msup,ZtevoTab,Endxi,Endxf,MeAv)
				
				Barnes="/home/bonfanti/Documents/Dottorato/Articoli/BrownSWPAgeGiroIsoc/tauBarnesKim2010.txt"
				call loadMatrix(Barnes,GyroTab,head)
				print*,'Storing ended'
!			!Printing test
!				print*,'First 3 Z',Zvec(1:3)
!				print*,'First 3 indices init',Zndxi(1:3)
!				print*,'First 3 indices fin',Zndxf(1:3)
!				print*,'1st row, 3rd isoc',IsocTab(Zndxi(3),:)
!				print*,'last row, 3rd isoc',IsocTab(Zndxf(3),:)
!				print*,'------'
!				print*,'Mass init',Minf
!				print*,'Mass fin',MTr-Mstep
!				print*,'size(nMt)',size(nMt)
!				print*,'nMt',nMt
!				print*,'All Zvec'
!				print*,Zvec
!				print*,'Z of tracks Ztvec',Ztvec
!				print*,'indxZt',indxZt
!				print*,'2nd Z, first 3 masses',Mavail(2,1:3)
!				print*,'1st Z T init',Tndxi(1,:)
!				print*,'1st Z T fin',Tndxf(1,:)
!				print*,'2nd Z, first 3 indices T init',Tndxi(2,1:3)
!				print*,'2nd Z, first 3 indices T fin',Tndxf(2,1:3)
!				print*,'first row, 2nd Z, 3rd track',TrackTab(Tndxi(2,3),:)
!				print*,'last row, 2nd Z, 3rd track',TrackTab(Tndxf(2,3),:)
!				print*,'last Z, indici05',Ztvec(size(Ztvec)),Tndxi(size(Ztvec),nMt(size(nMt))),Tndxf(size(Ztvec),nMt(size(nMt)))
!				print*,'first row, last Z, traccia05',TrackTab(Tndxi(size(Ztvec),nMt(size(nMt))),:)
!				print*,'last row, last Z, traccia05',TrackTab(Tndxf(size(Ztvec),nMt(size(nMt))),:)
!				print*,'first 2 rows GyroTab',GyroTab(1:2,:)
!				print*,'3rd and 4th row, 2nd Tab ZAMS',ZAMStab(ZAndxi(2)+2:ZAndxi(2)+3,:)
!				print*,'M evolution',MeAv
!				print*,'1st and 2nd row, 2nd Tab EvoZ',ZtevoTab(Endxi(2):Endxi(2)+1,:)
!				!!
			end if
	    	!!!!
	    	!!!!
			if (link.eq.1.and.priormag.eq.'p'.and.priordist.eq.'p') then !just to obtain L=L(d,v)
		  	  	if (priorteff.eq.'p') then
		  	  		TeffL=spec(1)
		  	  		I_TeffL=sspec(1)
		  	  		colL=-100.
		  	  	else
		  	  		TeffL=0.
		  	  		I_TeffL=0.
		  	  		colL=col
		  	  	end if
		  	  	star= (/ dble(1),spec(3),sspec(3),TeffL,I_TeffL,rhoI,DrhoI,loggI,I_loggI,vsiniI,I_vsiniI, &
					& protI,I_protI,logRHKI,YMgI,numaxI,I_numaxI,deltanuI,I_deltanuI,dble(sismo),-1.D0,-1.D0, &
					& mag,colL,dist,I_dist,Hipf,useColor,dble(idCol) /)
				call SCPmcmcPD12S(star,ltitle,IsocTab,Zvec,Zndxi,Zndxf,TrackTab,nMt,Mavail,indxZt,Tndxi,Tndxf,velTrackTab, &
						& Vndxi,Vndxf,GyroTab,ZAMStab,ZAndxi,ZAndxf,ZtevoTab,Endxi,Endxf,MeAv,massTmp,dmassTmp, &
						& radiusTmp,dradiusTmp,temp_s(link),dtemp_s,age(link),dage,starlum,dstarlum,0,row,acc)
				print*,'calibL',starlum,dstarlum
		  	end if
	    	
	    	star= (/ dble(1),met_s(link),dmet_s,TeffI,I_TeffI,rhoI,DrhoI,loggI,I_loggI,vsiniI,I_vsiniI, &
	    			& protI,I_protI,logRHKI,YMgI,numaxI,I_numaxI,deltanuI,I_deltanuI,dble(sismo),Rinp,I_Rinp, &
	    			& mag,colI,dist,I_dist,Hipf,useColor,dble(idCol) /)
	    			    	
	    	call SCPmcmcPD12S(star,ltitle,IsocTab,Zvec,Zndxi,Zndxf,TrackTab,nMt,Mavail,indxZt,Tndxi,Tndxf,velTrackTab, &
	    					& Vndxi,Vndxf,GyroTab,ZAMStab,ZAndxi,ZAndxf,ZtevoTab,Endxi,Endxf,MeAv,mass_Isoch,dmass_Isoch, &
	    					& radius_Isoch,dradius_Isoch,temp_s(link),dtemp_s,age(link),dage,Lum,I_Lum,link,row,acc)
	    	age(link)=age(link)/1.E9 !from yrs to Gyr
            dage=dage/1.E9
            radius_s(link)=(mass_s(link)/(rhoI/rhoSun))**(1./3.)
	    	
	    	if (acc.eq.0) then !Algorithm hasn't converged
	    		print*,'No convergence reached in theoretical models.'
	    		print*,'chain step: ',link,' chain ',(link-1)/nat+1
	    		accepted(link)='n'
	    	else
			  	rowTot=rowTot+row
			  	rowIO(link)=row
	    	end if
	    end if
	    teff=temp_s(link)
        IF(msrs.EQ.'y')THEN
          t=0
          DO WHILE(t.LT.1)
            CALL gasdev_s(harvest) 
            dum4 = msrsb + harvest*emsrsb
            IF(dum4.GT.0) t=2
          ENDDO
          mass_s(link)=rho(link)**(1./(1-3*dum4))
          radius_s(link) = (mass_s(link)/rho(link))**(1./3.)
        ENDIF
        IF(massfromr.EQ.'y')THEN
          mass_s(link)=rho(link)*radius_s(link)**3
        ENDIF
        IF(massfromr.NE.'y'.AND.enoch.NE.'y'.AND.msrs.NE.'y'.AND. &
           & fitmsrs.NE.'y'.AND.fixstellar.NE.'y'.and.isoch.ne.'y')THEN
          radius_s(link)=(mass_s(link)/rho(link))**(1./3.)
        ENDIF
        IF(stelincli.EQ.'y')THEN
          test=0
          nincli=0.
          DO WHILE(test.LT.1)
            CALL gasdev_s(harvest)
            vrotsini = spec(5) + sspec(5)*harvest  
            CALL gasdev_s(harvest)
            prot = spec(6) + sspec(6)*harvest  
            sinincli_s(link)=1000*vrotsini*prot*24*3600./(2*pi*radius_s(link)*sunra)
            IF(sinincli_s(link).GE.0.)THEN
              test=2
              IF(sinincli_s(link).LE.1)THEN
                incli_s(link)=DASIN(sinincli_s(link))*180/pi 
              ELSE
                incli_s(link)=90.
              ENDIF
            ENDIF
          ENDDO
        ENDIF
        lum_s(link)=(radius_s(link)**2)*(teff/TeffSun)**4
        logg(link) = LOG10(gravi*1.E6*(mass_s(link)*sunmass)/(radius_s(link)*sunra*1.E2)**2)
      ENDIF
      semi(i,link)=a_R(i,link)*radius_s(link)*sunra/ua
      irrad(i,link)=lum_s(link)/semi(i,link)**2
      radius_p(i,link)=rr(i,link)*radius_s(link)*sunra/jupra
      safro(i,link)=(semi(i,link)/radius_p(i,link))*(ua/jupra)
      r_R = a_R(i,link)*(1-exc(i,link)*DCOS(Eano_tr(i)))
      bcor = b(i,link)*(1-exc(i,link)*DCOS(Eano_tr(i)))
      x_tr = r_R*DSIN(Tano_tr(i) + omega(i,link)*pi/180 - pi/2.)
      incli = DACOS(bcor/r_R)
      IF(incli.LT.1e-6) incli=1e-6
      inclian(i,link) = incli*180/pi
      mass_p(i,link) = (mass_s(link)*sunmass)**(2/3.)*kb(i,link)/DSIN(incli)
      mass_p(i,link) = mass_p(i,link)/(2.*pi*gravi)**(1./3.)
      mass_p(i,link) = mass_p(i,link)/(jupmass*(24*3600.)**(-1/3.)) 
      m2_m1(i) = (mass_p(i,link)*jupmass)/(mass_s(link)*sunmass) 
      mass_p(i,link) = mass_p(i,link)*((1+m2_m1(i))**(2/3.))
      mdum = mass_p(i,link)*jupmass/earthmass
      IF(ironlow.EQ.'y')THEN
        CALL MR(mdum,rp_iron)
        testiron(i)=radius_p(i,link)*(jupra/earthra)*(1./rp_iron)
      ENDIF
      safro(i,link)=safro(i,link)*(mass_p(i,link)/mass_s(link))*(jupmass/sunmass)
      roche(i,link)=2.46*radius_p(i,link)*jupra*((mass_s(link)*sunmass)/(mass_p(i,link)*jupmass))**(1/3.)
      roche(i,link)=roche(i,link)/ua
      a_roche(i,link)=semi(i,link)/roche(i,link)
      mass_p_sini(i,link) = mass_p(i,link)*DSIN(incli)
      logg_p(i,link) = LOG10(gravi*1.E6*(mass_p(i,link)*jupmass)/(radius_p(i,link)*jupra*1.E2)**2)
      teq_p(i,link) = teff*SQRT((radius_s(link)*sunra)/(2*semi(i,link)*ua))
      hillrad_p(i,link)=semi(i,link)*ua/(radius_p(i,link)*jupra)
      hillrad_p(i,link)=hillrad_p(i,link)*(mass_p(i,link)*jupmass/(3*mass_s(link)*sunmass))**(1/3.)
      rhop(i,link) = mass_p(i,link)/radius_p(i,link)**3
      y_tr = r_R*DCOS(Tano_tr(i) + omega(i,link)*pi/180 - pi/2.)*DSIN(incli)
      z_tr =  -bcor*DCOS(Tano_tr(i) + omega(i,link)*pi/180 - pi/2.)
      b2(i,link) = SQRT(x_tr**2 + z_tr**2)
      Tano_se = (3*pi/2.) - omega(i,link)*pi/180
      Eano_se = 2.* DATAN(SQRT((1-exc(i,link))/(1+exc(i,link)))*DTAN(Tano_se/2.))
      r_R = a_R(i,link)*(1-exc(i,link)*DCOS(Eano_se))
      bcor = b(i,link)*(1-exc(i,link)*DCOS(Eano_se))
      x_se = r_R*DSIN(Tano_se + omega(i,link)*pi/180 - pi/2.)
      y_se = r_R*DCOS(Tano_se + omega(i,link)*pi/180 - pi/2.)*DSIN(incli)
      z_se =  -bcor*DCOS(Tano_se + omega(i,link)*pi/180 - pi/2.)
      b3(i,link) = SQRT(x_se**2 + z_se**2)
      ltt(i) = ABS(y_tr - y_se)*radius_s(link)*sunra/(light*24*3600)
      Mano_tr(i) = Eano_tr(i) - exc(i,link)*DSIN(Eano_tr(i))
      Mano_tr(i)=MOD(Mano_tr(i),2*pi)
      IF(Mano_tr(i).LT.0.)      Mano_tr(i) = Mano_tr(i)+2*pi
      Mano_se(i) = Eano_se - exc(i,link)*DSIN(Eano_se)
      Mano_se(i)=MOD(Mano_se(i),2*pi)
      IF(Mano_se(i).LT.0.)      Mano_se(i) = Mano_se(i)+2*pi
      Tperi(i) = t0(i,link) - (Mano_tr(i)*per(i,link)/(2.*pi))
      octime(i,link) = Tperi(i) + (Mano_se(i)*per(i,link)/(2.*pi))
      IF(octime(i,link).GT.t0(i,link)) octime(i,link) = octime(i,link) - per(i,link)
      prtr(i,link)=(radius_s(link)*sunra)
      prtr(i,link)=prtr(i,link)/(semi(i,link)*ua)
      proc(i,link)=prtr(i,link)
      prtr(i,link)=prtr(i,link)*(1+exc(i,link)*DCOS((pi/2)-(omega(i,link)*pi/180)))
      prtr(i,link)=prtr(i,link)/(1-exc(i,link)**2)
      proc(i,link)=proc(i,link)*(1+exc(i,link)*DCOS((3*pi/2)-(omega(i,link)*pi/180)))
      proc(i,link)=proc(i,link)/(1-exc(i,link)**2)
    ENDDO
! PHOTOMETRY
    IF(ntr.GT.0)THEN 
      photmerit(link)=0.
      DO i=1,ntr   
         photchi2(i,link)=0.
         nb = np(i)
         DO j=1,nb
            mulimb0(i,j)=0.
            IF(group(i).GT.0) mulimb2(i,j)=0.
            IF(overphot(i).EQ.'n')THEN
               nubin=1
            ELSE
               nubin=CEILING(texp(i,j)/binoverphot(i))
               nubin=nubin-MOD(nubin,2)+1
            ENDIF
            DO l2=1,nubin
               IF(l2.EQ.1)THEN
                  time = bjd(i,j)
               ELSE IF(l2.GT.1.AND.l2.LE.((nubin-1)/2+1))THEN
                  time = bjd(i,j)+DBLE(l2-1)*binoverphot(i)
               ELSE
                  time = bjd(i,j)-DBLE(l2-((nubin-1)/2+1))*binoverphot(i)
               ENDIF
               DO k=1,npla
                  epochtime = NINT((bjd(i,j)-t0(k,link))/per(k,link))
                  Mano = ((2*pi)/per(k,link))*(time-Tperi(k))
                  IF(isttv.EQ.'y'.AND.nttvmax.GT.0)THEN
                    DO l=1,nttv(k)
                      IF(epochtime.EQ.epochtr(k,l))THEN
                         Mano = ((2*pi)/per(k,link))*(time-Tperi(k)-ttv(k,l,link))
                      ENDIF 
                    ENDDO
                  ENDIF
                  Mano = MOD(Mano,2*pi)
                  IF(Mano.LT.0.) Mano = Mano+2*pi
                  IF((ABS(Mano-Mano_se(k)).LT.pi/2.).OR.(ABS(Mano-Mano_se(k)-2*pi).LT.pi/2.).OR. &
                  & (ABS(Mano-Mano_se(k)+2*pi).LT.pi/2.))THEN 
                    Mano = ((2*pi)/per(k,link))*(time-ltt(k)-Tperi(k))
                  ENDIF
                  Mano = MOD(Mano,2*pi)
                  IF(Mano.LT.0.) Mano = Mano+2*pi
                  IF(kesol.EQ.'s')THEN         ! Series solution of Kepler equation
                                               ! Not valid if e > 0.6627434
                     Eano = Mano + exc(k,link)*DSIN(Mano) + &
                     & (exc(k,link)**2)*(0.5*DSIN(2.*Mano)) + &
                     & (exc(k,link)**3)*((3/8.)*DSIN(3.*Mano)-(1/8.)*DSIN(Mano)) + &
                     & (exc(k,link)**4)*((1/3.)*DSIN(4.*Mano)-(1/6.)*DSIN(2.*Mano))
                  ELSE                         ! Numerical solution of Kepler equation
                     IF(Mano.GT.2*pi)THEN
                        test = 0 !. is an integer.
                        DO WHILE(test.LT.1)
                           Mano=Mano-2.*pi
                           IF(Mano.LE.2*pi) test=2 !. is an integer.
                        ENDDO 
                     ENDIF
                     IF(Mano.LT.0*pi)THEN
                        test = 0 !. is an integer.
                        DO WHILE(test.LT.1)
                           Mano=Mano+2.*pi
                           IF(Mano.GE.0*pi) test=2 !. is an integer.
                        ENDDO 
                     ENDIF
                     test2 = 1.
                     dum4 = 1.
                     eg0 = Mano ! + SIGN(dum4,test2)*0.85*exc(k,link)
                     DO WHILE(test2.GT.kepler)
                        feg = eg0 - exc(k,link)*DSIN(eg0) - Mano
                        fdeg = 1. - exc(k,link)*DCOS(eg0) 
                        fseg = exc(k,link)*DSIN(eg0)
                        fteg = exc(k,link)*DCOS(eg0)
                        delta_1 = - feg/fdeg
                        delta_2 = - feg/(fdeg+(0.5*delta_1*fseg))
                        delta_3 = - feg/(fdeg+(0.5*delta_2*fseg)+((delta_2**2)*fteg/6.))
                        eg1 = eg0 + delta_3
                        test2 = ABS(eg1 - eg0)
                        eg0 = eg1
                     ENDDO 
                     Eano = eg0
                  ENDIF
                  Tano = 2. * DATAN(SQRT((1+exc(k,link))/(1-exc(k,link)))*DTAN(Eano/2.))
                  ! test = (COS(Eano)-exc(k,link))/(1-exc(k,link)*DCOS(Eano))                
                  r_R = a_R(k,link)*(1-exc(k,link)*DCOS(Eano))
                  bcor = b(k,link)*(1-exc(k,link)*DCOS(Eano))
                  pos_x = r_R*DSIN(Tano + omega(k,link)*pi/180 - pi/2.)
                  pos_y = SQRT(r_R**2 - bcor**2)*DCOS(Tano + omega(k,link)*pi/180 - pi/2.)
                  pos_z = -bcor*DCOS(Tano + omega(k,link)*pi/180 - pi/2.)
                  z2(1)=SQRT(pos_x**2+pos_z**2)
!                  WRITE(733,*) bjd(i,j),pos_x,pos_y,pos_z         
                  IF(limb.EQ.'nl'.AND.pos_y.GT.0.)THEN
                     phef=1.
                     DO l=1,nfi
                        IF(filter(i).EQ.wfilter(l))THEN
                           c1 = nl(l,1,link)
                           c2 = nl(l,2,link)
                           c3 = nl(l,3,link)
                           c4 = nl(l,4,link)
                           codilu = dilution(l)
                           phase  = 2*pi*(time-t0(1,link))/per(1,link)
                           phase2 = phase - pi*phoffset(l,link)/180.
                           phaset = -DSIN(inclian(1,link)*pi/180.)*DCOS(phase2)
                           phef1  = -phampli1(l,link)*(DSIN(DACOS(phaset)) + (pi-DACOS(phaset))*phaset)/pi  ! REFL+OFFSET
                           cortide = SQRT((DCOS(phase))**2 + (tidel(k,link)*DSIN(phase))**2)
                           phef2 = 0.5*(1-tidel(1,link))*0.5*(dFsec(l,link,1)+dFsec(l,link,1)-phampli2(l,link)) & 
                             & *cortide                                          !ELLI
                           phef3  = phampli3(l,link)*DSIN(phase)                 !BEAM
                           phef   = -(phef1+phef2-phef3)
                        ENDIF
                     ENDDO 
                     IF(flagtr(i).GT.0)THEN
                        l=flagtr(i)
                        phase  = 2*pi*(time-t0(k,link))/per(k,link)
                        cortide = SQRT((DCOS(phase))**2 + (tidel(k,link)*DSIN(phase))**2)
                        ratio=SQRT((dF(k,link)+ddf(k,l,link))*cortide)
                        radipla(k,l,link)=SQRT(dF(k,link)+ddf(k,l,link))*radius_s(link)*sunra/jupra
                        dratio(k,l,link)=SQRT(dF(k,link)+ddf(k,l,link)) 
                        ddepth(k,l,link)=dratio(k,l,link)**2
                     ELSE
                       phase  = 2*pi*(time-t0(k,link))/per(k,link)
                       cortide = SQRT((DCOS(phase))**2 + (tidel(k,link)*DSIN(phase))**2)
                       ratio =  rr(k,link)*SQRT(cortide)
                     ENDIF
                     CALL occultnl(ratio,c1,c2,c3,c4,z2,theo,mulimb_nl,2) 
                     IF(k.EQ.1)THEN
                        mulimb0t=1+(theo(1)*(1+phef)-1.)/codilu
                     ELSE
                        mulimb0t=1+(theo(1)*(((mulimb0t-1)*codilu)+1)-1.)/codilu
                     ENDIF
                     IF(group(i).GT.0)THEN
                       IF(k.EQ.1)THEN
                         phase  = 2*pi*(time-t0(k,link))/per(k,link)
                         cortide = SQRT((DCOS(phase))**2 + (tidel(k,link)*DSIN(phase))**2)
                         ratiog=SQRT(dfgroup(group(i),link)*cortide)
                         CALL occultnl(ratiog,c1,c2,c3,c4,z2,theo,mulimb_nl,2) 
                         mulimb0g=1+(theo(1)*(1+phef)-1.)/codilu
                       ELSE
                         mulimb0g=mulimb0t
                       ENDIF
                     ENDIF
                  ELSE IF(limb.EQ.'qd'.AND.pos_y.GT.0.)THEN
                     phef=1.
                     DO l=1,nfi
                        IF(filter(i).EQ.wfilter(l))THEN
                           q1 = ql(l,1,link)
                           q2 = ql(l,2,link)
                           codilu = dilution(l)
                           phase  = 2*pi*(time-t0(1,link))/per(1,link)
                           phase2 = phase - pi*phoffset(l,link)/180.
                           phaset = -DSIN(inclian(1,link)*pi/180.)*DCOS(phase2)
                           phef1  = -phampli1(l,link)*(DSIN(DACOS(phaset)) + (pi-DACOS(phaset))*phaset)/pi  ! REFL+OFFSET
                            
                           cortide = SQRT((DCOS(phase))**2 + (tidel(k,link)*DSIN(phase))**2)
                           phef2 = 0.5*(1-tidel(1,link))*0.5*(dFsec(l,link,1)+dFsec(l,link,1)-phampli2(l,link)) &
                             & *cortide                                         !ELLI
                           phef3  = phampli3(l,link)*DSIN(phase)                !BEAM
                           phef   = -(phef1+phef2-phef3)
                        ENDIF
                     ENDDO
                     IF(flagtr(i).GT.0)THEN  
                        l=flagtr(i)
                        phase  = 2*pi*(time-t0(k,link))/per(k,link)
                        cortide = SQRT((DCOS(phase))**2 + (tidel(k,link)*DSIN(phase))**2)
                        ratio=SQRT((dF(k,link)+ddf(k,l,link))*cortide)
                        radipla(k,l,link)=SQRT(dF(k,link)+ddf(k,l,link))*radius_s(link)*sunra/jupra
                        dratio(k,l,link)=SQRT(dF(k,link)+ddf(k,l,link))
                        ddepth(k,l,link)=dratio(k,l,link)**2
                     ELSE
                       phase  = 2*pi*(time-t0(k,link))/per(k,link)
                       cortide = SQRT((DCOS(phase))**2 + (tidel(k,link)*DSIN(phase))**2)
                       ratio =  rr(k,link)*SQRT(cortide)
                     ENDIF
                     CALL occultquad(z2,q1,q2,ratio,theo,mulimb_qd,1) 
                     IF(k.EQ.1)THEN
                        mulimb0t=1+(theo(1)*(1+phef)-1.)/codilu
                     ELSE
                        mulimb0t=1+(theo(1)*(((mulimb0t-1)*codilu)+1)-1.)/codilu
                     ENDIF
                     IF(group(i).GT.0)THEN
                       IF(k.EQ.1)THEN
                         phase  = 2*pi*(time-t0(k,link))/per(k,link)
                         cortide = SQRT((DCOS(phase))**2 + (tidel(k,link)*DSIN(phase))**2)
                         ratiog=SQRT(dfgroup(group(i),link)*cortide)
                         CALL occultquad(z2,q1,q2,ratiog,theo,mulimb_qd,1) 
                         mulimb0g=1+(theo(1)*(1+phef)-1.)/codilu
                       ELSE
                         mulimb0g=mulimb0t
                       ENDIF
                     ENDIF
                  ELSE IF(limb.EQ.'no'.AND.pos_y.GT.0.)THEN
                    q1 = 0.
                    q2 = 0.
                    phef=1.
                    DO l=1,nfi
                       IF(filter(i).EQ.wfilter(l))THEN
                          codilu = dilution(l)
                          phase  = 2*pi*(time-t0(1,link))/per(1,link)
                          phase2 = phase - pi*phoffset(l,link)/180.
                          phaset = -DSIN(inclian(1,link)*pi/180.)*DCOS(phase2)
                          phef1  = -phampli1(l,link)*(DSIN(DACOS(phaset)) + (pi-DACOS(phaset))*phaset)/pi  ! REFL+OFFSET
                          cortide = SQRT((DCOS(phase))**2 + (tidel(k,link)*DSIN(phase))**2)
                          phef2 = 0.5*(1-tidel(1,link))*0.5*(dFsec(l,link,1)+dFsec(l,link,1)-phampli2(l,link)) &
                            & *cortide                                         !ELLI
                          phef3  = phampli3(l,link)*DSIN(phase)                !BEAM
                          phef   = -(phef1+phef2+phef3)
                       ENDIF
                    ENDDO
                    IF(flagtr(i).GT.0)THEN  
                       l=flagtr(i)
                       phase  = 2*pi*(time-t0(k,link))/per(k,link)
                       cortide = SQRT((DCOS(phase))**2 + (tidel(k,link)*DSIN(phase))**2)
                       ratio=SQRT((dF(k,link)+ddf(k,l,link))*cortide)
                       radipla(k,l,link)=SQRT(dF(k,link)+ddf(k,l,link))*radius_s(link)*sunra/jupra
                       dratio(k,l,link)=SQRT(dF(k,link)+ddf(k,l,link))
                       ddepth(k,l,link)=dratio(k,l,link)**2
                     ELSE
                       phase  = 2*pi*(time-t0(k,link))/per(k,link)
                       cortide = SQRT((DCOS(phase))**2 + (tidel(k,link)*DSIN(phase))**2)
                       ratio =  rr(k,link)*SQRT(cortide)
                    ENDIF
                    CALL occultquad(z2,q1,q2,ratio,theo,mulimb_qd,1) 
                    IF(k.EQ.1)THEN
                        mulimb0t=1+(theo(1)*(1+phef)-1.)/codilu
                     ELSE
                        mulimb0t=1+(theo(1)*(((mulimb0t-1)*codilu)+1)-1.)/codilu
                    ENDIF
                    IF(group(i).GT.0)THEN
                      IF(k.EQ.1)THEN
                        phase  = 2*pi*(time-t0(k,link))/per(k,link)
                        cortide = SQRT((DCOS(phase))**2 + (tidel(k,link)*DSIN(phase))**2)
                        ratiog=SQRT(dfgroup(group(i),link)*cortide)
                        CALL occultquad(z2,q1,q2,ratiog,theo,mulimb_qd,1) 
                        mulimb0g=1+(theo(1)*(1+phef)-1.)/codilu
                      ELSE
                        mulimb0g=mulimb0t
                      ENDIF
                    ENDIF
                 ELSE IF(pos_y.LE.0.)THEN !OCCULTATION
                    q1 = 0.
                    q2 = 0.
                    phef=1.
                    DO l=1,nfi
                       IF(filter(i).EQ.wfilter(l))THEN
                          codilu = dilution(l)
                          phase  = 2*pi*(time-t0(1,link))/per(1,link)
                          phase2 = phase - pi*phoffset(l,link)/180.
                          phaset = -DSIN(inclian(1,link)*pi/180.)*DCOS(phase2)
                          phef1  = -phampli1(l,link)*(DSIN(DACOS(phaset)) + (pi-DACOS(phaset))*phaset)/pi  ! REFL+OFFSET
                          IF(z2(1).GT.1.AND.pos_x.GE.0.AND.testphef(l).EQ.0)THEN
                            phefsave = phef1
                          ELSE IF(z2(1).LE.1)THEN
                            phef1 = phefsave
                            testphef(l) = 1
                          ENDIF
                          cortide = SQRT((DCOS(phase))**2 + (tidel(k,link)*DSIN(phase))**2)
                          phef2 = 0.5*(1-tidel(1,link))*0.5*(dFsec(l,link,1)+dFsec(l,link,1)-phampli2(l,link)) &
                           & *cortide                                          !ELLI
                          phef3  = phampli3(l,link)*DSIN(phase)                !BEAM
                          phef   = -(phef1+phef2-phef3)
                       ENDIF
                    ENDDO
                    DO l=1,nfi
                      IF(filter(i).EQ.wfilter(l))THEN
                        ratio = rr(k,link)
                        ratio2 = SQRT(dFsec(l,link,k))
                      ENDIF
                    ENDDO
                    CALL occultquad(z2,q1,q2,ratio,theo,mulimb_qd,1)
                    IF(k.EQ.1)THEN
                       IF(ratio.GT.0.)THEN
                          theo(1)=((theo(1)-1)*(ratio2/ratio)**2)+1.
                          mulimb0t=1+(theo(1)*(1+phef)-1.)/codilu
                          IF(group(i).GT.0) mulimb0g=mulimb0t
                       ELSE
                          theo(1)=1.
                          mulimb0t=1+(theo(1)*(1+phef)-1.)/codilu
                          IF(group(i).GT.0) mulimb0g=mulimb0t
                       ENDIF
                    ELSE
                       IF(ratio.GT.0.)THEN
                          theo(1)=((theo(1)-1)*(ratio2/ratio)**2)+1
                          mulimb0t=1+(theo(1)*(((mulimb0t-1)*codilu)+1)-1.)/codilu  
                          IF(group(i).GT.0) mulimb0g=mulimb0t
                       ELSE
                          theo(1)=1.
                          mulimb0t=1+(theo(1)*(((mulimb0t-1)*codilu)+1)-1.)/codilu  
                          IF(group(i).GT.0) mulimb0g=mulimb0t
                       ENDIF
                    ENDIF
                 ENDIF
              ENDDO
              mulimb0(i,j)=mulimb0(i,j)+mulimb0t/DBLE(nubin)
              IF(group(i).GT.0) mulimb2(i,j)=mulimb2(i,j)+mulimb0g/DBLE(nubin)
           ENDDO
        ENDDO
        dof = nb - njump
        DO j=1,nb
            flaremodel=1.
            IF(flarenumber(i).GT.0)THEN
             k = flarenumber(i)
             DO l = 1,k
              IF(bjd(i,j).GT.flt0(i,l,link))THEN
               dtime = bjd(i,j)-flt0(i,l,link)
               flaremodel = flaremodel + flampli(i,l,link)*expi**(-dtime/fltau(i,l,link))
              ENDIF   
             ENDDO
            ENDIF 
            mulimb0(i,j)=mulimb0(i,j)*flaremodel
            mulimb1(i,j)=mulimb0(i,j)
            syscor(i,j)=1.
        ENDDO
     ENDDO
     DO i=1,ntr   
        photchi2(i,link)=0.
        nb = np(i)
        IF(group(i).GT.0)THEN
          t=group(i)
          IF(donegroup(t).EQ.'n')THEN
            donegroup(t)='y'
            DO j=1,nb
              pond = 0.
              DO l=1,ntr
                IF(group(l).EQ.t)THEN
                  pond = pond + 1./(error(l,j)**2)
                ENDIF
              ENDDO
              DO l=1,ntr
                IF(group(l).EQ.t)THEN
                   resigroup(t,j)=resigroup(t,j)+phot(l,j)/(pond*mulimb2(l,j)*error(l,j)**2)
!                    resigroup(t,j)=resigroup(t,j)+phot(l,j)/mulimb2(l,j)
                ENDIF
              ENDDO
            ENDDO 
          ENDIF
        ENDIF
        IF(nsyspar(i).GT.0)THEN
          ALLOCATE(po(nb,ne2));ALLOCATE(y(nb));ALLOCATE(sig(nb));ALLOCATE(z(nb))
          ALLOCATE(u(nb,nsysparmax));ALLOCATE(v(nsysparmax,nsysparmax))
          ALLOCATE(w(nsysparmax));ALLOCATE(cvm(nsysparmax,nsysparmax))
          ALLOCATE(a(nsysparmax));ALLOCATE(nep(15))
          nep=0
          nep(1)=timeorder(i)
          nep(2)=colororder(i)
          nep(3)=fwhmorder(i)
          nep(4)=fwhmxorder(i)
          nep(5)=fwhmyorder(i)
          nep(6)=skyorder(i)
          nep(7)=pporder(i)
          nep(8)=ppforder(i)
          nep(9)=sinusnumber(i)
          nep(10)=ramporder(i) 
          nep(11)=jumpnumber(i)
          nep(12)=jumporder(i)
          nep(13)=offsetnumber(i)
          IF(group(i).GT.0)THEN
            nep(14)=grouporder(i)
          ENDIF
          IF(pldcor(i).GT.0)THEN
            nep(15)=pldcor(i)
          ENDIF
          temp_sys=bjd(i,1)-10./(60.*24.)
          po=0.
          DO j=1,nb
            po(j,1)=bjd(i,j)-temp_sys
            po(j,2)=10**(airmass(i,j))
            po(j,3)=fwhm(i,j)
            po(j,4)=fwhmx(i,j)
            po(j,5)=fwhmy(i,j)
            po(j,6)=sky(i,j)
            po(j,7)=dX(i,j)
            po(j,8)=dY(i,j)
            IF(sinusnumber(i).GT.0)THEN
              po(j,9)=(bjd(i,j)-sit0(i,1,link))/sip(i,1,link)
              po(j,9)=DSIN(po(j,9)*2*pi)
            ENDIF
            IF(sinusnumber(i).GT.1)THEN
              po(j,10)=(bjd(i,j)-sit0(i,2,link))/sip(i,2,link)
              po(j,10)=DSIN(po(j,10)*2*pi)
            ENDIF
            IF(sinusnumber(i).GT.2)THEN
              po(j,11)=(bjd(i,j)-sit0(i,3,link))/sip(i,3,link)
              po(j,11)=DSIN(po(j,11)*2*pi)
            ENDIF
            IF(sinusnumber(i).GT.3)THEN
              po(j,12)=(bjd(i,j)-sit0(i,4,link))/sip(i,4,link)
              po(j,12)=DSIN(po(j,12)*2*pi)
            ENDIF
            IF(offsetnumber(i).GT.0)THEN
              po(j,13)=bjd(i,j)-offsettime(i,1)
            ENDIF
            IF(offsetnumber(i).GT.1)THEN
              po(j,14)=bjd(i,j)-offsettime(i,2)
            ENDIF
            IF(offsetnumber(i).GT.2)THEN
              po(j,15)=bjd(i,j)-offsettime(i,3)
            ENDIF
            IF(offsetnumber(i).GT.3)THEN
              po(j,16)=bjd(i,j)-offsettime(i,4)
            ENDIF
            IF(ramporder(i).GT.0)THEN
              IF(rampmod.EQ.'log')THEN
                po(j,17)=LOG(po(j,1))
              ELSE
                po(j,17)=expi**(-po(j,1)/t1ramp(i,link))
              ENDIF
            ENDIF
            IF(ramporder(i).GT.1)THEN
              IF(rampmod.EQ.'log')THEN
                po(j,18)=LOG(po(j,1))**2
              ELSE
                po(j,18)=expi**(-po(j,1)/t2ramp(i,link))
              ENDIF
            ENDIF
            IF(group(i).GT.0)THEN
              t=group(i)
              IF(grouporder(i).GT.0)THEN
                po(j,19)=resigroup(t,j)
              ENDIF
              IF(grouporder(i).GT.1)THEN
                po(j,20)=resigroup(t,j)**2
              ENDIF  
              IF(grouporder(i).GT.2)THEN
                po(j,21)=resigroup(t,j)**3
              ENDIF
              IF(grouporder(i).GT.3)THEN
                po(j,22)=resigroup(t,j)**4
              ENDIF
            ENDIF
            IF(pldcor(i).GT.0)THEN
              DO l=1,9
                po(j,22+l)= pixcon(i,j,l) /SUM(pixcon(i,j,1:9))
              ENDDO
            ENDIF
            if(jumpnumber(i).GT.0.)THEN
             k = jumpnumber(i)
             DO l=1,k
               IF(l.EQ.1)THEN
                 IF(bjd(i,j).GE.tju(i,l).AND.bjd(i,j).LT.tju(i,l+1))THEN
                   po(j,32) = LOG(bjd(i,j)-(tju(i,l)-5./(60.*24.)))
                 ENDIF
               ELSE IF(l.LT.k)THEN
                 IF(bjd(i,j).GE.tju(i,l).AND.bjd(i,j).LT.tju(i,l+1))THEN
                   po(j,33) = LOG(bjd(i,j)-(tju(i,l)-5./(60.*24.)))
                 ENDIF
               ELSE
                 IF(bjd(i,j).GE.tju(i,l))THEN
                   po(j,33) = LOG(bjd(i,j)-(tju(i,l)-5./(60.*24.)))
                 ENDIF
               ENDIF
             ENDDO
            ENDIF
            y(j)=phot(i,j)/mulimb1(i,j)
            sig(j)=error(i,j)
          ENDDO
          CALL svdfit(po,y,sig,nb,a,nsysparmax,u,v,w,chisq,baseline,nep,ne2)
          CALL svdvar(v,nsysparmax,nsysparmax,w,cvm,nsysparmax)
          DO j=1,nsysparmax 
             CALL gasdev_s(harvest)
             sysva(i,j)=a(j)+harvest*SQRT(cvm(j,j))
          ENDDO
          DO j=1,nb
             z(j)=a(1)                                                      ! Normalization scalar
             IF(nep(1).GE.1) z(j)=z(j)+a(2)*po(j,1)                         ! time polynomial order 1     
             IF(nep(1).GE.2) z(j)=z(j)+a(3)*po(j,1)**2                      !                 order 2 
             IF(nep(1).GE.3) z(j)=z(j)+a(4)*po(j,1)**3                      !                 order 3
             IF(nep(1).GE.4) z(j)=z(j)+a(5)*po(j,1)**4                      !                 order 4
             muph(j)=z(j)
             IF(nep(2).GE.1) z(j)=z(j)+a(6)*po(j,2)                         ! 10**(airmass) pol order 1
             IF(nep(2).GE.2) z(j)=z(j)+a(7)*po(j,2)**2                      !                   order 2 
             IF(nep(2).GE.3) z(j)=z(j)+a(8)*po(j,2)**3                      !                   order 3 
             IF(nep(2).GE.4) z(j)=z(j)+a(9)*po(j,2)**4                      !                   order 4 

             IF(nep(3).GE.1) z(j)=z(j)+a(10)*po(j,3)                        ! FWHM polynomial order 1
             IF(nep(3).GE.2) z(j)=z(j)+a(11)*po(j,3)**2                     !                 order 2 
             IF(nep(3).GE.3) z(j)=z(j)+a(12)*po(j,3)**3                     !                 order 3 
             IF(nep(3).GE.4) z(j)=z(j)+a(13)*po(j,3)**4                     !                 order 4 

             IF(nep(4).GE.1) z(j)=z(j)+a(14)*po(j,4)                        ! FWHM_x polynomial order 1 
             IF(nep(4).GE.2) z(j)=z(j)+a(15)*po(j,4)**2                     !                   order 2 
             IF(nep(4).GE.3) z(j)=z(j)+a(16)*po(j,4)**3                     !                   order 3 
             IF(nep(4).GE.4) z(j)=z(j)+a(17)*po(j,4)**4                     !                   order 4 

             IF(nep(5).GE.1) z(j)=z(j)+a(18)*po(j,5)                        ! FWHM_y polynomial order 1 
             IF(nep(5).GE.2) z(j)=z(j)+a(19)*po(j,5)**2                     !                   order 2 
             IF(nep(5).GE.3) z(j)=z(j)+a(20)*po(j,5)**3                     !                   order 3 
             IF(nep(5).GE.4) z(j)=z(j)+a(21)*po(j,5)**4                     !                   order 4 

             IF(nep(6).GE.1) z(j)=z(j)+a(22)*po(j,6)                        ! background polynomial order 1 
             IF(nep(6).GE.2) z(j)=z(j)+a(23)*po(j,6)**2                     !                       order 2 
             IF(nep(6).GE.3) z(j)=z(j)+a(24)*po(j,6)**3                     !                       order 3 
             IF(nep(6).GE.4) z(j)=z(j)+a(25)*po(j,6)**4                     !                       order 4 

             IF(nep(7).GE.1)THEN                                            ! x- and y-position polynomial order 1
               z(j)=z(j)+a(26)*po(j,7)
               z(j)=z(j)+a(27)*po(j,8)
             ENDIF
             IF(nep(7).GE.2)THEN                                            ! x- and y-position polynomial order 2
               z(j)=z(j)+a(28)*po(j,7)**2
               z(j)=z(j)+a(29)*po(j,8)**2
               z(j)=z(j)+a(30)*po(j,7)*po(j,8)
             ENDIF
             IF(nep(7).GE.3)THEN                                            ! x- and y-position polynomial order 3
               z(j)=z(j)+a(31)*po(j,7)**3
               z(j)=z(j)+a(32)*po(j,8)**3
               z(j)=z(j)+a(33)*po(j,7)*po(j,8)**2
               z(j)=z(j)+a(34)*po(j,8)*po(j,7)**2
             ENDIF
             IF(nep(7).GE.4)THEN                                            ! x- and y-position polynomial order 4
               z(j)=z(j)+a(35)*po(j,7)**4
               z(j)=z(j)+a(36)*po(j,8)**4
               z(j)=z(j)+a(37)*po(j,8)*po(j,7)**3
               z(j)=z(j)+a(38)*po(j,7)*po(j,8)**3
               z(j)=z(j)+a(39)*(po(j,7)**2)*(po(j,8)**2)
             ENDIF 

             IF(nep(8).GE.1)THEN                                            ! cross-terms position and FWHM order 2
               z(j)=z(j)+a(40)*po(j,7)*po(j,3)**2
               z(j)=z(j)+a(41)*po(j,8)*po(j,3)**2
               z(j)=z(j)+a(42)*po(j,7)*po(j,4)**2
               z(j)=z(j)+a(43)*po(j,8)*po(j,4)**2
               z(j)=z(j)+a(44)*po(j,7)*po(j,5)**2
               z(j)=z(j)+a(45)*po(j,8)*po(j,5)**2
             ENDIF

             IF(nep(9).GE.1) z(j)=z(j)+a(46)*po(j,9)                        ! sinus function 1
             IF(nep(9).GE.2) z(j)=z(j)+a(47)*po(j,10)                       ! sinus function 2
             IF(nep(9).GE.3) z(j)=z(j)+a(48)*po(j,11)                       ! sinus function 3
             IF(nep(9).GE.4) z(j)=z(j)+a(49)*po(j,12)                       ! sinus function 4

             IF(nep(10).GE.1) z(j)=z(j)+a(50)*po(j,17)                       ! ramp order 1
             IF(nep(10).GE.2) z(j)=z(j)+a(51)*po(j,18)                       ! ramp order 2

             IF(nep(13).GE.1.AND.po(j,13).GT.0.) z(j)=z(j)+a(52)             ! offset 1
             IF(nep(13).GE.2.AND.po(j,14).GT.0.) z(j)=z(j)+a(53)             ! offset 2
             IF(nep(13).GE.3.AND.po(j,15).GT.0.) z(j)=z(j)+a(54)             ! offset 3
             IF(nep(13).GE.4.AND.po(j,16).GT.0.) z(j)=z(j)+a(55)             ! offset 4

             IF(nep(14).GE.1) z(j)=z(j)+a(56)*po(j,19)                       ! group model
             IF(nep(14).GE.2) z(j)=z(j)+a(57)*po(j,20)                      
             IF(nep(14).GE.3) z(j)=z(j)+a(58)*po(j,21)
             IF(nep(14).GE.4) z(j)=z(j)+a(59)*po(j,22)

             IF(nep(15).GE.0)THEN                                            ! PLD model
               z(j) = z(j)+a(60)*po(j,23)
               z(j) = z(j)+a(61)*po(j,24)
               z(j) = z(j)+a(62)*po(j,25)
               z(j) = z(j)+a(63)*po(j,26)
               z(j) = z(j)+a(64)*po(j,27)
               z(j) = z(j)+a(65)*po(j,28)
               z(j) = z(j)+a(66)*po(j,29)
               z(j) = z(j)+a(67)*po(j,30)
               z(j) = z(j)+a(68)*po(j,31)
             ENDIF
             IF(nep(15).GE.1)THEN                                            
               z(j) = z(j)+a(69)*po(j,23)**2
               z(j) = z(j)+a(70)*po(j,24)**2
               z(j) = z(j)+a(71)*po(j,25)**2
               z(j) = z(j)+a(72)*po(j,26)**2
               z(j) = z(j)+a(73)*po(j,27)**2
               z(j) = z(j)+a(74)*po(j,28)**2
               z(j) = z(j)+a(75)*po(j,29)**2
               z(j) = z(j)+a(76)*po(j,30)**2
               z(j) = z(j)+a(77)*po(j,31)**2
             ENDIF
             IF(nep(11).GE.1.AND.nep(12).GT.0.)THEN   ! Jump 1
              z(j)=z(j)+a(78)*po(j,32)
             ENDIF
             IF(nep(11).GE.1.AND.nep(12).GT.1)THEN 
              z(j)=z(j)+a(79)*po(j,32)**2  
             ENDIF   
             IF(nep(11).GE.2.AND.nep(12).GT.0.)THEN   ! Jump 2
              z(j)=z(j)+a(80)*po(j,33)
             ENDIF
             IF(nep(11).GE.2.AND.nep(12).GT.1)THEN 
              z(j)=z(j)+a(81)*po(j,33)**2  
             ENDIF  
             syscor(i,j)=syscor(i,j)*z(j)
             mulimb1(i,j)=mulimb1(i,j)*syscor(i,j)
          ENDDO
          DEALLOCATE(po);DEALLOCATE(y);DEALLOCATE(sig);DEALLOCATE(z)
          DEALLOCATE(u);DEALLOCATE(v);DEALLOCATE(w);DEALLOCATE(cvm)
          DEALLOCATE(a);DEALLOCATE(nep)
        ELSE
          DO j=1,nb
            muph(j)=1.
          ENDDO
        ENDIF
        IF(bliss(i).NE.'n')THEN
          IF(link.EQ.1)THEN
            test=0
            corbliss=11.
            DO WHILE(test.LT.1)                
! STEP 1: tuning of the box size to get in the box with the largest number of points at least nbliss points
              corbliss=corbliss-1.
              IF(nbliss.GT.NINT(0.1*nb))THEN
                nbliss2 = NINT(0.1*nb)
              ELSE
                nbliss2 = nbliss
              ENDIF
              ndot(i)=NINT(corbliss*SQRT(DBLE(nb)/DBLE(nbliss2)))
              nmpix = ndot(i)
              nmpiy = ndot(i)
              ALLOCATE(xjh(nmpix,nmpiy));ALLOCATE(yjh(nmpix,nmpiy))
              ALLOCATE(npjh(nmpix,nmpiy))
              midx = (MAXVAL(dX(i,1:nb))-MINVAL(dX(i,1:nb)))/DBLE(nmpix)
              midy = (MAXVAL(dY(i,1:nb))-MINVAL(dY(i,1:nb)))/DBLE(nmpiy)
              DO l=1,nmpiy
                DO k=1,nmpix
                  xjh(k,l)=MINVAL(dX(i,1:nb))+DBLE(k-1)*midx+midx/2.
                  yjh(k,l)=MINVAL(dY(i,1:nb))+DBLE(l-1)*midy+midy/2.
                  npjh(k,l)=0
                ENDDO
              ENDDO
              DO j=1,nb
                DO l=1,nmpiy
                  DO k=1,nmpix
                    IF((ABS(dX(i,j)-xjh(k,l)).LE.midx/2.+EPSILON(midx)).AND. &
                    & (ABS(dY(i,j)-yjh(k,l)).LE.midy/2.+EPSILON(midy)))THEN
                      npjh(k,l)=npjh(k,l)+1
                    ENDIF
                  ENDDO
                ENDDO
              ENDDO 
              IF(corbliss.LT.2.) test=2
              IF(MAXVAL(npjh).GT.(nbliss2*3.)) test=2
              DEALLOCATE(xjh);DEALLOCATE(yjh);DEALLOCATE(npjh)
            ENDDO
          ENDIF
          nmpix = ndot(i)
          nmpiy = ndot(i)
          ALLOCATE(maxt(nmpix,nmpiy));ALLOCATE(mint(nmpix,nmpiy))
          ALLOCATE(xjh(nmpix,nmpiy));ALLOCATE(yjh(nmpix,nmpiy))
          ALLOCATE(fluxjh(nmpix,nmpiy));ALLOCATE(npjh(nmpix,nmpiy))
          ALLOCATE(cotime(nmpix,nmpiy))
          midx = (MAXVAL(dX(i,1:nb))-MINVAL(dX(i,1:nb)))/DBLE(nmpix)
          midy = (MAXVAL(dY(i,1:nb))-MINVAL(dY(i,1:nb)))/DBLE(nmpiy)
          DO l=1,nmpiy
            DO k=1,nmpix
              xjh(k,l)=MINVAL(dX(i,1:nb))+DBLE(k-1)*midx+midx/2.
              yjh(k,l)=MINVAL(dY(i,1:nb))+DBLE(l-1)*midy+midy/2.
              fluxjh(k,l)=0.
              npjh(k,l)=0
            ENDDO
          ENDDO 
          IF(link.EQ.1)THEN
            DO j=1,nb       
              DO l=1,nmpiy
                DO k=1,nmpix
                  IF((ABS(dX(i,j)-xjh(k,l)).LE.midx/2.+EPSILON(midx)).AND. &
                  & (ABS(dY(i,j)-yjh(k,l)).LE.midy/2.+EPSILON(midy)))THEN
                  fluxjh(k,l)=fluxjh(k,l)+phot(i,j)/mulimb1(i,j)
                   npjh(k,l)=npjh(k,l)+1
                   xbox(i,j)=k;ybox(i,j)=l
                  ENDIF
                ENDDO
              ENDDO
            ENDDO
! STEP 2: shuffling of the measurements and redefinition of the boxes
            DO j=1,nb 
              k=xbox(i,j);l=ybox(i,j)    
              k3=k;l3=l
              dis2=1E12
              IF(npjh(k,l).LT.nbliss2)THEN 
                DO k2=1,nmpix
                  DO l2=1,nmpiy
                    dis=SQRT((xjh(k2,l2)-xjh(k,l))**2 + (yjh(k2,l2)-yjh(k,l))**2)
                    IF((dis.LT.dis2).AND. &
                      & (npjh(k2,l2).GT.nbliss2))THEN 
                      dis2=dis
                      k3=k2;l3=l2
                    ENDIF
                  ENDDO
                ENDDO              
                npjh(k,l)=npjh(k,l)-1
                npjh(k3,l3)=npjh(k3,l3)+1
                xbox(i,j)=k3;ybox(i,j)=l3  
              ENDIF
            ENDDO
          ENDIF
          xjh=0.
          yjh=0.
          npjh=0
          maxt=0. 
          mint=1E12
          fluxjh=0.
          DO j=1,nb
            k=xbox(i,j);l=ybox(i,j) 
            fluxjh(k,l)=fluxjh(k,l)+phot(i,j)/mulimb1(i,j)
            xjh(k,l)=xjh(k,l)+dX(i,j)
            yjh(k,l)=yjh(k,l)+dY(i,j)
            npjh(k,l)=npjh(k,l)+1
            IF(bjd(i,j).GT.maxt(k,l)) maxt(k,l)=bjd(i,j)
            IF(bjd(i,j).LT.mint(k,l)) mint(k,l)=bjd(i,j)
          ENDDO    
          DO k=1,nmpix
            DO l=1,nmpiy   
              xjh(k,l)=xjh(k,l)/DBLE(npjh(k,l))
              yjh(k,l)=yjh(k,l)/DBLE(npjh(k,l))
            ENDDO  
          ENDDO
          meancotime=0.
          nbox=0 
          nbox2=MAXVAL(npjh)
          DO l=1,nmpiy
            DO k=1,nmpix
              IF(npjh(k,l).GT.0)THEN 
                fluxjh(k,l)=fluxjh(k,l)/DBLE(npjh(k,l))  
                IF(npjh(k,l).LT.nbox2) nbox2 = npjh(k,l)
              ELSE
                fluxjh(k,l)=0.
              ENDIF
              cotime(k,l)=maxt(k,l)-mint(k,l)
              IF(fluxjh(k,l).GT.0.)THEN
                nbox=nbox+1
                meancotime=meancotime+cotime(k,l)
              ENDIF 
            ENDDO
          ENDDO
          IF(nbox.GT.0) meancotime=meancotime/nbox
          IF(link.EQ.1.AND.i.EQ.1) PRINT*, '--------------------------------------------'
!          IF(link.EQ.1) PRINT*, 'LC ',i,': <BLISS nbox>: ',nbox,nbox2
          IF(link.EQ.1) WRITE(*,'(A2,I3,A24,3(1x,I4))') 'LC ',i,': Nbox, Np/box(MIN,MAX): ', &    
          & nbox,nbox2,MAXVAL(npjh)
          IF(link.EQ.1.AND.i.EQ.ntr) PRINT*, '--------------------------------------------'
          DO j=1,nb
            blissmodel=1.  
            nbox2=0      
            k=xbox(i,j);l=ybox(i,j) 
            test=0 !. is an integer
            IF((dX(i,j)-xjh(k,l)).GE.0.)THEN
              k2=k+1
            ELSE
              k2=k-1
            ENDIF
            IF((dY(i,j)-yjh(k,l)).GE.0.)THEN
              l2=l+1
            ELSE
              l2=l-1
            ENDIF   
            IF((k2.GE.1).AND.(k2.LE.nmpix).AND. & 
             &(l2.GE.1).AND.(l2.LE.nmpiy))THEN
               IF((fluxjh(k,l).GT.0).AND. &
                &(fluxjh(k2,l).GT.0).AND. & 
                &(fluxjh(k,l2).GT.0).AND. &
               &(fluxjh(k2,l2).GT.0))THEN
                test=2
              ENDIF
            ENDIF
            IF(test.LT.1)THEN ! Nearest-Neighbor Interpolation (NNI) 
              blissmodel=fluxjh(k,l)
              nbox2=nbox2+1
            ELSE ! Bilinear Interpolation (BLI)
              blissmodel=fluxjh(k,l)*((xjh(k2,l2)-dX(i,j))/(xjh(k2,l2)-xjh(k,l)))* &
              & ((yjh(k2,l2)-dY(i,j))/(yjh(k2,l2)-yjh(k,l))) 
              blissmodel=blissmodel+fluxjh(k2,l)*((dX(i,j)-xjh(k,l))/&
              & (xjh(k2,l2)-xjh(k,l)))*((yjh(k2,l2)-dY(i,j))/(yjh(k2,l2)-yjh(k,l)))
              blissmodel=blissmodel+fluxjh(k,l2)*((xjh(k2,l2)-dX(i,j))/&
              & (xjh(k2,l2)-xjh(k,l)))*((dY(i,j)-yjh(k,l))/(yjh(k2,l2)-yjh(k,l)))
              blissmodel=blissmodel+fluxjh(k2,l2)*((dX(i,j)-xjh(k,l))/&
              & (xjh(k2,l2)-xjh(k,l)))*((dY(i,j)-yjh(k,l))/(yjh(k2,l2)-yjh(k,l)))
             ENDIF
            syscor(i,j)=syscor(i,j)*blissmodel
            mulimb1(i,j)=mulimb1(i,j)*blissmodel
          ENDDO
          DEALLOCATE(maxt);DEALLOCATE(mint);DEALLOCATE(xjh)
          DEALLOCATE(yjh);DEALLOCATE(fluxjh);DEALLOCATE(npjh)
          DEALLOCATE(cotime)
        ENDIF
        DO j=1,nb
          photcor(i,j)=phot(i,j)/syscor(i,j)
          photcor2(i,j)=photcor(i,j)*muph(j)
          photchi2(i,link) = photchi2(i,link)+((phot(i,j)-mulimb1(i,j))/error(i,j))**2
          resi(i,j) = phot(i,j)-mulimb1(i,j)
          model(i,j) = mulimb1(i,j)
          model_tr(i,j) = mulimb0(i,j) 
          model_tr2(i,j) = mulimb0(i,j)*muph(j)
        ENDDO 
        rephotchi2(i,link) = photchi2(i,link)/DBLE(dof)
        photmerit(link) = photmerit(link) + photchi2(i,link)  
      ENDDO
      rephotmerit(link) = photmerit(link)/DBLE(nphotot-njump)
      IF(link.EQ.1)THEN
        PRINT*, 'Statistics for the initial phot. residuals'
        PRINT*, '        P_R  Np/bin   Sig    <error>   Sig_N    Beta_W   Beta_R    CorF'
        WRITE(444,*) 'Statistics for the initial phot. residuals'
        WRITE(444,*) '        P_R  Np/bin   Sig    <error>   Sig_N    Beta_W   Beta_R    CorF'
        WRITE(445,*) 'Statistics for the initial phot. residuals'
        WRITE(445,*) '        P_R  Np/bin   Sig    <error>   Sig_N    Beta_W   Beta_R    CorF'
      ENDIF
      IF(link.EQ.na)THEN
        PRINT*, '--------------------------------------------' 
        WRITE(444,*) '--------------------------------------------' 
        WRITE(445,*) '--------------------------------------------' 
        PRINT*, 'Statistics for the best-fit phot. residuals'
        PRINT*, '        P_R  Np/bin   Sig    <error>   Sig_N    Beta_W   Beta_R    CorF'
        WRITE(444,*) 'Statistics for the first chain phot. residuals'
        WRITE(444,*) '        P_R  Np/bin   Sig    <error>   Sig_N    Beta_W   Beta_R    CorF'
        WRITE(445,*) 'Statistics for the first chain phot. residuals'
        WRITE(445,*) '        P_R  Np/bin   Sig    <error>   Sig_N    Beta_W   Beta_R    CorF'
      ENDIF
      DO i=1,ntr
      	nb = np(i)
        resave = 0.
        resmeanerror = 0.
        resrms(i,link) = 0.
        DO j=1,nb
          resave = resave + resi(i,j)/DBLE(nb)
          resmeanerror = resmeanerror + error(i,j)/DBLE(nb)
        ENDDO
        DO j=1,nb
          resrms(i,link) = resrms(i,link) + (resi(i,j)-resave)**2
        ENDDO
        resrms(i,link) = SQRT(resrms(i,link)/DBLE(nb))
        bigbetared=0.
        DO l=1,ne3
           IF(l.EQ.1)THEN
              k=0
              ce = nbin(i,l)
              DO j=1,ce-1
                 resibin(i,j,l)=0.
                 photbin(i,j)=0.
                 photcorbin(i,j)=0.
                 photcor2bin(i,j)=0.
                 eresibin(i,j)=0.
                 ephotbin(i,j)=0.
                 rephotbin(i,j)=0.
                 ephotcorbin(i,j)=0.
                 ephotcor2bin(i,j)=0.
                 bjdbin(i,j)=0.
                 dXbin(i,j)=0.
                 dYbin(i,j)=0.
                 fwhmbin(i,j)=0.
                 fwhmxbin(i,j)=0.
                 fwhmybin(i,j)=0.
                 airmassbin(i,j)=0.
                 skybin(i,j)=0.
                 modelbin(i,j)=0.
                 model_trbin(i,j)=0.
                 model_tr2bin(i,j)=0.
                 lim1=k+1
                 lim2=k+np_bin(i,l)
                 k=k+np_bin(i,l)
                 DO t=lim1,lim2
                   resibin(i,j,l)=resibin(i,j,l)+resi(i,t)/DBLE(np_bin(i,l))
                   photbin(i,j)=photbin(i,j)+phot(i,t)/DBLE(np_bin(i,l))
                   rephotbin(i,j)=rephotbin(i,j)+error(i,t)/DBLE(np_bin(i,l))
                   photcorbin(i,j)=photcorbin(i,j)+photcor(i,t)/DBLE(np_bin(i,l))
                   photcor2bin(i,j)=photcor2bin(i,j)+photcor2(i,t)/DBLE(np_bin(i,l))
                   bjdbin(i,j)=bjdbin(i,j)+bjd(i,t)/DBLE(np_bin(i,l))
                   dXbin(i,j)=dXbin(i,j)+dX(i,t)/DBLE(np_bin(i,l))
                   dYbin(i,j)=dYbin(i,j)+dY(i,t)/DBLE(np_bin(i,l))
                   fwhmbin(i,j)=fwhmbin(i,j)+fwhm(i,t)/DBLE(np_bin(i,l))
                   fwhmxbin(i,j)=fwhmxbin(i,j)+fwhmx(i,t)/DBLE(np_bin(i,l))
                   fwhmybin(i,j)=fwhmybin(i,j)+fwhmy(i,t)/DBLE(np_bin(i,l))
                   airmassbin(i,j)=airmassbin(i,j)+airmass(i,t)/DBLE(np_bin(i,l))
                   skybin(i,j)=skybin(i,j)+sky(i,t)/DBLE(np_bin(i,l))
                   modelbin(i,j)=modelbin(i,j)+model(i,t)/DBLE(np_bin(i,l))
                   model_trbin(i,j)=model_trbin(i,j)+model_tr(i,t)/DBLE(np_bin(i,l)) 
                   model_tr2bin(i,j)=model_tr2bin(i,j)+model_tr2(i,t)/DBLE(np_bin(i,l)) 
                 ENDDO
                 rephotbin(i,j)=rephotbin(i,j)/SQRT(DBLE(np_bin(i,l)))
                 DO t=lim1,lim2
                   eresibin(i,j)=eresibin(i,j)+(resi(i,t)-resibin(i,j,l))**2
                   ephotbin(i,j)=ephotbin(i,j)+(phot(i,t)-photbin(i,j))**2
                   ephotcorbin(i,j)=ephotcorbin(i,j)+(photcor(i,t)-photcorbin(i,j))**2
                   ephotcor2bin(i,j)=ephotcor2bin(i,j)+(photcor2(i,t)-photcor2bin(i,j))**2
                 ENDDO
                 eresibin(i,j)=SQRT(eresibin(i,j))/np_bin(i,l)
                 ephotbin(i,j)=SQRT(ephotbin(i,j))/np_bin(i,l)
                 ephotcorbin(i,j)=SQRT(ephotcorbin(i,j))/np_bin(i,l)
                 ephotcor2bin(i,j)=SQRT(ephotcor2bin(i,j))/np_bin(i,l)
              ENDDO
              resibin(i,ce,l)=0.
              photbin(i,ce)=0.
              photcorbin(i,ce)=0.
              photcor2bin(i,ce)=0.
              eresibin(i,ce)=0.
              ephotbin(i,ce)=0.
              rephotbin(i,ce)=0.
              ephotcorbin(i,ce)=0.
              ephotcor2bin(i,ce)=0.
              bjdbin(i,ce)=0.
              dXbin(i,ce)=0.
              dYbin(i,ce)=0.
              fwhmbin(i,ce)=0.
              fwhmxbin(i,ce)=0.
              fwhmybin(i,ce)=0.
              airmassbin(i,ce)=0.
              skybin(i,ce)=0.
              modelbin(i,ce)=0.
              model_trbin(i,ce)=0.
              model_tr2bin(i,ce)=0.
              lim1=k+1
              lim2=np(i)
              DO t=lim1,lim2
                resibin(i,ce,l)=resibin(i,ce,l)+resi(i,t)/DBLE(nlast(i,l))
                photbin(i,ce)=photbin(i,ce)+phot(i,t)/DBLE(nlast(i,l))
                rephotbin(i,ce)=rephotbin(i,ce)+error(i,t)/SQRT(DBLE(nlast(i,l)))
                photcorbin(i,ce)=photcorbin(i,ce)+photcor(i,t)/DBLE(nlast(i,l))
                photcor2bin(i,ce)=photcor2bin(i,ce)+photcor2(i,t)/DBLE(nlast(i,l))
                bjdbin(i,ce)=bjdbin(i,ce)+bjd(i,t)/DBLE(nlast(i,l))
                dXbin(i,j)=dXbin(i,ce)+dX(i,t)/DBLE(nlast(i,l))
                dYbin(i,j)=dYbin(i,ce)+dY(i,t)/DBLE(nlast(i,l))
                fwhmbin(i,j)=fwhmbin(i,ce)+fwhm(i,t)/DBLE(nlast(i,l))
                fwhmxbin(i,j)=fwhmxbin(i,ce)+fwhmx(i,t)/DBLE(nlast(i,l))
                fwhmybin(i,j)=fwhmybin(i,ce)+fwhmy(i,t)/DBLE(nlast(i,l))
                airmassbin(i,j)=airmassbin(i,ce)+airmass(i,t)/DBLE(nlast(i,l))
                skybin(i,j)=skybin(i,ce)+sky(i,t)/DBLE(nlast(i,l))
                modelbin(i,ce)=modelbin(i,ce)+model(i,t)/DBLE(nlast(i,l))
                model_trbin(i,ce)=model_trbin(i,ce)+model_tr(i,t)/DBLE(nlast(i,l))
                model_tr2bin(i,ce)=model_tr2bin(i,ce)+model_tr2(i,t)/DBLE(nlast(i,l))
              ENDDO
              rephotbin(i,ce)=rephotbin(i,ce)/SQRT(DBLE(nlast(i,l)))
              DO t=lim1,lim2
                eresibin(i,ce)=eresibin(i,ce)+(resi(i,t)-resibin(i,j,l))**2
                ephotbin(i,ce)=ephotbin(i,ce)+(phot(i,t)-photbin(i,j))**2
                ephotcorbin(i,ce)=ephotcorbin(i,ce)+(photcor(i,t)-photcorbin(i,j))**2
                ephotcor2bin(i,ce)=ephotcor2bin(i,ce)+(photcor2(i,t)-photcor2bin(i,j))**2
              ENDDO
              eresibin(i,ce)=SQRT(eresibin(i,ce))/DBLE(nlast(i,l))
              ephotbin(i,ce)=SQRT(ephotbin(i,ce))/DBLE(nlast(i,l))
              ephotcorbin(i,ce)=SQRT(ephotcorbin(i,ce))/DBLE(nlast(i,l))
              ephotcor2bin(i,ce)=SQRT(ephotcor2bin(i,ce))/DBLE(nlast(i,l))
          ELSE
              k=0
              ce = nbin(i,l)
              DO j=1,ce-1
                 resibin(i,j,l)=0.
                 lim1=k+1
                 lim2=k+np_bin(i,l)
                 k=k+np_bin(i,l)
                 DO t=lim1,lim2
                   resibin(i,j,l)=resibin(i,j,l)+resi(i,t)/DBLE(np_bin(i,l))
                 ENDDO
              ENDDO
              resibin(i,ce,l)=0.
              lim1=k+1
              lim2=np(i)
              DO t=lim1,lim2
                resibin(i,ce,l)=resibin(i,ce,l)+resi(i,t)/DBLE(nlast(i,l))
              ENDDO
          ENDIF
          resavebin = 0.
          resrmsbin(i,l) = 0.
          DO j=1,ce
             resavebin = resavebin + resibin(i,j,l)/DBLE(ce)
          ENDDO
          DO j=1,ce
            resrmsbin(i,l) = resrmsbin(i,l) + (resibin(i,j,l)-resavebin)**2
          ENDDO
          resrmsbin(i,l) = SQRT(resrmsbin(i,l)/DBLE(ce-1))
          dum_beta_red=SQRT(DBLE(npave(i,l)))*(resrmsbin(i,l)/resrms(i,link))
          IF(dum_beta_red.LT.1) dum_beta_red=1.
          IF(maxred.EQ.'y')THEN
            IF(l.EQ.dtred)THEN
              bigbetared = dum_beta_red
              beta_red(i)=dum_beta_red
              bnpave(i)=npave(i,l)
              bresrmsbin(i)=resrmsbin(i,l)
              bred_dur(i)= red_dur(l)*24*60
            ENDIF
          ELSE
            IF(l.GT.1.AND.dum_beta_red.GT.bigbetared.AND.np_bin(i,l).GE.3)THEN
              bigbetared = dum_beta_red
              beta_red(i)=dum_beta_red
              bnpave(i)=npave(i,l)
              bresrmsbin(i)=resrmsbin(i,l)
              bred_dur(i)= red_dur(l)*24*60
            ENDIF
          ENDIF
        ENDDO
        IF(link.EQ.1)THEN 
            WRITE(*,119) i,bred_dur(i),bnpave(i),resrms(i,link),resmeanerror,bresrmsbin(i),resrms(i,link)/resmeanerror, &
            & beta_red(i),beta_red(i)*resrms(i,link)/resmeanerror
            WRITE(444,119) i,bred_dur(i),bnpave(i),resrms(i,link),resmeanerror,bresrmsbin(i),resrms(i,link)/resmeanerror, &
            & beta_red(i),beta_red(i)*resrms(i,link)/resmeanerror
            WRITE(445,119) i,bred_dur(i),bnpave(i),resrms(i,link),resmeanerror,bresrmsbin(i),resrms(i,link)/resmeanerror, &
            & beta_red(i),beta_red(i)*resrms(i,link)/resmeanerror
        ELSE IF(link.EQ.na)THEN
           WRITE(*,119) i,best_bred_dur(i),best_bnpave(i),resrms(i,soluce),resmeanerror,best_bresrmsbin(i), & 
             & resrms(i,soluce)/resmeanerror,best_beta_red(i),best_beta_red(i)*resrms(i,soluce)/resmeanerror  
           WRITE(444,119) i,best_bred_dur(i),best_bnpave(i),resrms(i,soluce),resmeanerror,best_bresrmsbin(i), &
            & resrms(i,soluce)/resmeanerror,best_beta_red(i),best_beta_red(i)*resrms(i,soluce)/resmeanerror 
           WRITE(445,119) i,best_bred_dur(i),best_bnpave(i),resrms(i,soluce),resmeanerror,best_bresrmsbin(i), &
            & resrms(i,soluce)/resmeanerror,best_beta_red(i),best_beta_red(i)*resrms(i,soluce)/resmeanerror 
        ENDIF
      ENDDO
      IF(link.EQ.1.OR.link.EQ.na)THEN
        PRINT*, '--------------------------------------------' 
        WRITE(444,*) '--------------------------------------------' 
        WRITE(445,*) '--------------------------------------------' 
      ENDIF
    ENDIF
! RV
    IF(nrv.GT.0)THEN
      i2=0
      IF(testf2.EQ.'y')THEN     
         ktide(link) = f2(link)*(1.5/(a_R(1,link)**3))*(2*pi/(per(1,link)*3600.*24.))
         ktide(link) = ktide(link)*(DSIN(inclian(1,link)*pi/180.))**2
         ktide(link) = ktide(link)*(mass_p(1,link)*jupmass)/(mass_s(link)*sunmass)
         ktide(link) = ktide(link)*radius_s(link)*sunra
      ENDIF
      DO i=1,nrv
        rvchi2(i,link)=0. 
        nb = nprv(i)
        DO j=1,nb
           rvmod(i,j)=0.
           IF(overrv(i).EQ.'n')THEN
              nubin=1
           ELSE
              nubin=CEILING(rvtexp(i,j)/binoverrv(i))
              nubin=nubin-MOD(nubin,2)+1
           ENDIF
           DO k=1,npla  
              rvmodpla(k,i,j)=0.
              DO l2=1,nubin
                 IF(l2.EQ.1)THEN
                    time = rvbjd(i,j)
                 ELSE IF(l2.GT.1.AND.l2.LE.((nubin-1)/2+1))THEN   
                    time = rvbjd(i,j)+DBLE(l2-1)*binoverrv(i)
                 ELSE
                    time = rvbjd(i,j)-DBLE(l2-((nubin-1)/2+1))*binoverrv(i)
                 ENDIF
                 Mano = ((2*pi)/per(k,link))*(time-Tperi(k))
                 Mano = MOD(Mano,2*pi)
                 IF(Mano.LT.0.) Mano = Mano+2*pi     
                 IF(kesol.EQ.'s')THEN         ! Series solution of Kepler equation
                                              ! Not valid if e > 0.6627434
                    Eano = Mano + exc(k,link)*DSIN(Mano) + &
                    & (exc(k,link)**2)*(0.5*DSIN(2.*Mano)) + &
                    & (exc(k,link)**3)*((3/8.)*DSIN(3.*Mano)-(1/8.)*DSIN(Mano)) + &
                    & (exc(k,link)**4)*((1/3.)*DSIN(4.*Mano)-(1/6.)*DSIN(2.*Mano))
                 ELSE                         ! Numerical solution of Kepler equation
                    IF(Mano.GT.2*pi)THEN
                       test = 0 !. is an integer
                       DO WHILE(test.LT.1)
                          Mano=Mano-2.*pi
                          IF(Mano.LE.2*pi) test=2 !. is an integer  
                       ENDDO 
                    ENDIF
                    IF(Mano.LT.0*pi)THEN
                       test = 0 !. is an integer
                       DO WHILE(test.LT.1)
                          Mano=Mano+2.*pi
                          IF(Mano.GE.0*pi) test=2 !. is an integer
                       ENDDO 
                    ENDIF
                    test2 = 1.
                    dum4 = 1.
                    eg0 = Mano  ! + SIGN(dum4,test2)*0.85*exc(k,link)
                    DO WHILE(test2.GT.kepler)
                       feg = eg0 - exc(k,link)*DSIN(eg0) - Mano
                       fdeg = 1. - exc(k,link)*DCOS(eg0) 
                       fseg = exc(k,link)*DSIN(eg0)
                       fteg = exc(k,link)*DCOS(eg0)
                       delta_1 = - feg/fdeg
                       delta_2 = - feG/(fdeg+(0.5*delta_1*fseg))
                       delta_3 = - feG/(fdeg+(0.5*delta_2*fseg)+((delta_2**2)*fteg/6.))
                       eg1 = eg0 + delta_3
                       test2 = ABS(eg1 - eg0)
                       eg0 = eg1
                    ENDDO 
                    Eano = eg0   
                 ENDIF     
                 Tano = 2. * DATAN(SQRT((1+exc(k,link))/(1-exc(k,link)))*DTAN(Eano/2.))
                 ! test = (COS(Eano)-exc(k,link))/(1-exc(k,link)*DCOS(Eano)) 
                 r_R = a_R(k,link)*(1-exc(k,link)*DCOS(Eano))
                 bcor = b(k,link)*(1-exc(k,link)*DCOS(Eano))
                 pos_x = r_R*DSIN(Tano + omega(k,link)*pi/180 - pi/2.)
                 pos_y = SQRT(r_R**2 - bcor**2)*DCOS(Tano + omega(k,link)*pi/180 - pi/2.)
                 pos_z = -bcor*DCOS(Tano + omega(k,link)*pi/180 - pi/2.)
                 rvz(j)=SQRT(pos_x**2+pos_z**2) 
                 rossiter(j) = 0.
                 delta = rvz(j)/a_R(k,link)
                 thet = Tano - Tano_tr(k)
                 kk = rr(k,link)
                 rs = (1.+ kk)/a_R(k,link)
                 us = rold(link,1) + rold(link,2)
                 ud = rold(link,1) - rold(link,2)
                 IF(k.EQ.1.AND.rossitif.EQ.'y')THEN
                    CALL RVC(delta,rs,kk,us,ud,rossi)
                    CALL LUMC(delta,rs,kk,us,ud,lc1)
                    IF(pos_y.GT.0.)THEN
                       rossiter(j) = (rossi/(delta*lc1))*(((vsini(link)*DSIN(beta(link)*pi/180))* &   
                       & (b(k,link)/a_R(k,link))*DCOS(thet))-(vsini(link)*DCOS(beta(link)*pi/180)*DSIN(thet)))
                       rossiter(j)=rossiter(j)*1000
                    ENDIF
                 ENDIF
                 IF(k.EQ.1)THEN
                    IF(testf2.EQ.'y')THEN
                      rvmodplat =rossiter(j) + ka(k,link)*(exc(k,link)*DCOS(omega(k,link)*pi/180) + &
                      &  DCOS(Tano+omega(k,link)*pi/180))-ktide(link)*DCOS(2*(Tano+omega(k,link)*pi/180)) 
                    ELSE
                      rvmodplat =rossiter(j) + ka(k,link)*(exc(k,link)*DCOS(omega(k,link)*pi/180) + &
                      &  DCOS(Tano+omega(k,link)*pi/180)) 
                    ENDIF
                 ELSE
                    rvmodplat = ka(k,link)*(exc(k,link)*DCOS(omega(k,link)*pi/180) + &
                    &  DCOS(Tano+omega(k,link)*pi/180))
                 ENDIF
                 rvmod(i,j)=rvmod(i,j)+rvmodplat/DBLE(nubin)
                 rvmodpla(k,i,j)=rvmodpla(k,i,j)+rvmodplat/DBLE(nubin)
              ENDDO
           ENDDO
           rvresi(i,j) = rv(i,j) - rvmod(i,j)   
           rvmod2(i,j)=rvmod(i,j)
        ENDDO
        ALLOCATE(po(nb,5));ALLOCATE(y(nb));ALLOCATE(sig(nb))
        ALLOCATE(u(nb,nrvsysparmax));ALLOCATE(v(nrvsysparmax,nrvsysparmax))
        ALLOCATE(w(nrvsysparmax));ALLOCATE(cvm(nrvsysparmax,nrvsysparmax))
        ALLOCATE(a(nrvsysparmax));ALLOCATE(nep(5))
        nep(1)=rvtimeorder(i)
        nep(2)=rvfwhmorder(i)
        nep(3)=rvbisorder(i)
        nep(4)=rvcontrastorder(i)
        nep(5)=rvloghkorder(i)
        temp_sys=rvbjd(i,1)-10./(60.*24.)
        DO j=1,nb
          po(j,1)=(rvbjd(i,j)-temp_sys)/365.25
          po(j,2)=rvfwhm(i,j)
          po(j,3)=rvbis(i,j)
          po(j,4)=rvcontrast(i,j)
          po(j,5)=rvloghk(i,j)
          y(j)=rvresi(i,j)
          sig(j)=rverror(i,j)
        ENDDO 
        CALL svdfit(po,y,sig,nb,a,nrvsysparmax,u,v,w,chisq,rvbaseline,nep,5)
        CALL svdvar(v,nrvsysparmax,nrvsysparmax,w,cvm,nrvsysparmax)
        DO j=1,nrvsysparmax 
          CALL gasdev_s(harvest)
          rvsysva(i,j)=a(j)+harvest*SQRT(cvm(j,j))
        ENDDO
        WRITE(447,*) a(1)
        DO j=1,nb
          i2=i2+1
          tempu = a(1)
          IF(nep(1).GE.1) tempu=tempu+a(2)*po(j,1)
          IF(nep(1).GE.2) tempu=tempu+a(3)*po(j,1)**2
          IF(nep(1).GE.3) tempu=tempu+a(4)*po(j,1)**3
          IF(nep(1).GE.4) tempu=tempu+a(5)*po(j,1)**4
          IF(nep(2).GE.1) tempu=tempu+a(6)*po(j,2)
          IF(nep(2).GE.2) tempu=tempu+a(7)*po(j,2)**2
          IF(nep(2).GE.3) tempu=tempu+a(8)*po(j,2)**3
          IF(nep(2).GE.4) tempu=tempu+a(9)*po(j,2)**4
          IF(nep(3).GE.1) tempu=tempu+a(10)*po(j,3)
          IF(nep(3).GE.2) tempu=tempu+a(11)*po(j,3)**2
          IF(nep(3).GE.3) tempu=tempu+a(12)*po(j,3)**3
          IF(nep(3).GE.4) tempu=tempu+a(13)*po(j,3)**4
          IF(nep(4).GE.1) tempu=tempu+a(14)*po(j,4)
          IF(nep(4).GE.2) tempu=tempu+a(15)*po(j,4)**2
          IF(nep(4).GE.3) tempu=tempu+a(16)*po(j,4)**3
          IF(nep(4).GE.4) tempu=tempu+a(17)*po(j,4)**4
          IF(nep(5).GE.1) tempu=tempu+a(18)*po(j,5)
          IF(nep(5).GE.2) tempu=tempu+a(19)*po(j,5)**2
          IF(nep(5).GE.3) tempu=tempu+a(20)*po(j,5)**3
          IF(nep(5).GE.4) tempu=tempu+a(21)*po(j,5)**4
          rvmod(i,j)=rvmod(i,j)+tempu
          rv2(i,j)=rv(i,j)-tempu
          rvresi(i,j)=rv(i,j)-rvmod(i,j)
          IF(iftrendrv.EQ.'y')THEN
            glorvbjd(i2) = rvbjd(i,j)
            glorvresi(i2) = rvresi(i,j)
            glorverror(i2) = rverror(i,j) 
          ENDIF
        ENDDO
        DEALLOCATE(po);DEALLOCATE(y);DEALLOCATE(sig)
        DEALLOCATE(u);DEALLOCATE(v);DEALLOCATE(w);DEALLOCATE(cvm)
        DEALLOCATE(a);DEALLOCATE(nep)
      ENDDO
      IF(iftrendrv.EQ.'y')THEN
        ALLOCATE(po(nrvtot,1));ALLOCATE(y(nrvtot));ALLOCATE(sig(nrvtot))
        ALLOCATE(u(nrvtot,nrvparglo));ALLOCATE(v(nrvparglo,nrvparglo));ALLOCATE(w(nrvparglo))
        ALLOCATE(cvm(nrvparglo,nrvparglo));ALLOCATE(a(nrvparglo));ALLOCATE(nep(1))
        nep(1)=ordertrendrv
        temp_sys=glorvbjd(1)-10./(60.*24.)
        DO j=1,nrvtot
          po(j,1)=(glorvbjd(j)-temp_sys)/365.25
          y(j)=glorvresi(j)
          sig(j)=glorverror(j)
        ENDDO
        CALL svdfit(po,y,sig,nrvtot,a,nrvparglo,u,v,w,chisq,rvglo,nep,1)
        CALL svdvar(v,nrvparglo,nrvparglo,w,cvm,nrvparglo)
        DO j=1,nrvparglo
          CALL gasdev_s(harvest)
          rvtrendco(link,j)=a(j)+harvest*SQRT(cvm(j,j))
        ENDDO
        i2=0
        DO i=1,nrv
          nb = nprv(i)
          DO j=1,nb
            i2=i2+1
            rvmod(i,j)=rvmod(i,j)+a(1)
            IF(ordertrendrv.GE.1) rvmod(i,j)=rvmod(i,j)+a(2)*po(i2,1)
            IF(ordertrendrv.GE.2) rvmod(i,j)=rvmod(i,j)+a(3)*po(i2,1)*po(i2,1)
            rvresi(i,j)=rv(i,j)-rvmod(i,j)
          ENDDO  
        ENDDO
        DEALLOCATE(po);DEALLOCATE(y);DEALLOCATE(sig)
        DEALLOCATE(u);DEALLOCATE(v);DEALLOCATE(w);DEALLOCATE(cvm)
        DEALLOCATE(a);DEALLOCATE(nep)
      ENDIF
      DO i=1,nrv
        nb=nprv(i)
        DO j=1,nb
          IF(npla.GT.1)THEN
             DO l=1,npla
               rv2pla(l,i,j)=rv2(i,j)
               DO i2=1,npla
                 IF(l.NE.i2)THEN
                    rv2pla(l,i,j)=rv2pla(l,i,j)-rvmodpla(i2,i,j)
                 ENDIF
               ENDDO
             ENDDO        
           ELSE
             rv2pla(1,i,j)=rv2(i,j)
           ENDIF
        ENDDO
        DO j=1,nb
          rvchi2(i,link) = rvchi2(i,link) + ((rvmod(i,j)-rv(i,j))/rverror(i,j))**2
        ENDDO
        dof = nb - njump_rv
        IF(njump_rv.GE.nb)THEN
          dof = 1 !. is an integer
        ENDIF
        rervchi2(i,link) = rvchi2(i,link)/DBLE(dof)
        rvmerit(link) = rvmerit(link) + rvchi2(i,link)
      ENDDO
      IF(njump_rv.GE.nrvtot)THEN
        rervmerit(link) = rvmerit(link)
      ELSE
        rervmerit(link) = rvmerit(link)/DBLE(nrvtot-njump_rv)
      ENDIF
      IF(link.EQ.1)THEN
        PRINT*, 'Statistics for the initial RV residuals'
        PRINT*, '         Sig    <error>  jitter  RChi2'
      ENDIF
      IF(link.EQ.na)THEN
        PRINT*, 'Statistics for the best-fit RV residuals'
        PRINT*, '         Sig    <error>  jitter  RChi2'
        WRITE(444,*) 'Statistics for the first chain RV. residuals'
        WRITE(444,*) '         Sig    <error>  jitter  RChi2'
        WRITE(445,*) 'Statistics for the first chain RV. residuals'
        WRITE(445,*) '         Sig    <error>  jitter  RChi2'
      ENDIF
      DO i=1,nrv 
        nb=nprv(i)
        resave = 0.
        resmeanerror = 0.
        rvresrms(i) = 0.
        DO j=1,nb
          resave = resave + rvresi(i,j)/DBLE(nb)
          resmeanerror = resmeanerror + rverror(i,j)/DBLE(nb)
        ENDDO
        DO j=1,nb
          rvresrms(i) = rvresrms(i) + (rvresi(i,j)-resave)**2
        ENDDO
        rvresrms(i) = SQRT(rvresrms(i)/DBLE(nb))
        IF((rvresrms(i)**2 - resmeanerror**2).GE.0.)THEN
          jitter(i) = SQRT(rvresrms(i)**2 - resmeanerror**2)
        ELSE
          jitter(i) = 0.
        ENDIF
        IF(link.EQ.1)THEN
          WRITE(*,117) i,rvresrms(i),resmeanerror,jitter(i),rervchi2(i,link)
        ENDIF
        IF(link.EQ.na)THEN
          WRITE(*,117) i,rvresrms_soluce(i),resmeanerror,best_jitter(i),rervchi2(i,soluce)
          WRITE(444,117) i,rvresrms_soluce(i),resmeanerror,best_jitter(i),rervchi2(i,soluce)
          WRITE(445,117) i,rvresrms_soluce(i),resmeanerror,best_jitter(i),rervchi2(i,soluce)
        ENDIF
      ENDDO 
      IF(link.EQ.1.OR.link.EQ.na)THEN
        PRINT*, '---------------------------------------------'
        WRITE(444,*) '---------------------------------------------'
        WRITE(445,*) '---------------------------------------------'
      ENDIF
    ENDIF
  
    dof = nphotot + nrvtot - njump
    merit(link) = photmerit(link) + rvmerit(link) 

! TRANSIT TIMING PENALTY
   IF(ntiming.GT.0.AND.istiming.EQ.'y')THEN                
     merit(link) = merit(link) + timerit(link)
   ENDIF

! DDF PENALTY
   IF(ntr.GT.0.AND.nddf.GT.0.AND.isddf.EQ.'p')THEN
     DO k=1,npla
       DO i=1,nddf         
         merit(link)=merit(link)+(ddf(k,i,link)/sigddf)**2
       ENDDO
     ENDDO
   ENDIF

! GROUP WHITE LC DEPTH PENALTY
   IF(ngroup.GT.0)THEN
     DO i=1,ngroup
        merit(link)=merit(link)+((dfgroup(i,link)-dfgroup_ini(i))/edfgroup_ini(i))**2
     ENDDO
   ENDIF

! SINUS PERIOD PENALTY
   IF(ntr.GT.0)THEN
     DO i=1,ntr
        IF(sinusnumber(i).GT.0)THEN
           k=sinusnumber(i)
           DO l=1,k
             merit(link)=merit(link)+((sip_ini(i,l)-sip(i,l,link))/esip_ini(i,l))**2
           ENDDO
        ENDIF
     ENDDO
   ENDIF

! OCCULTATION TIMING PENALTY 
    IF(ntimingoc.GT.0.AND.istimingoc.EQ.'y')THEN              
      IF(link.EQ.1) ndata =ndata+DBLE(ntimingoc)
      DO i=1,ntimingoc     
        test2 = (timingoc(i)-octime(1,link))/per(1,link)
        epoch_oc(i) = NINT(test2)
        ttiming = octime(1,link) + epoch_oc(i)*per(1,link)
        merit(link) = merit(link) + ((timing(i) - ttiming)/stimingoc(i))**2
      ENDDO
    ENDIF

! STELLAR MASS NORMAL PRIOR
    IF(priormass.eq.'p')THEN     
      merit(link) = merit(link) + ((mass_s(link)-starmass)/dmass_s_ini)**2
    ENDIF
    if (MjumpIso.eq.'y') then
      merit(link) = merit(link) + ((mass_s(link)-mass_Isoch)/dmass_Isoch)**2
    endif

! STELLAR TEFF NORMAL PRIOR
    IF(priorteff.eq.'p')THEN 
      merit(link) = merit(link) + ((temp_s(link)-spec(1))/sspec(1))**2
    ENDIF
    
! Stellar color normal prior
	IF(priorcol.eq.'p')THEN 
      merit(link) = merit(link) + ((col_s(link)-col)/I_col)**2
    ENDIF
  
! STELLAR METALLICITY PRIOR
    IF(priormet.eq.'p')THEN   
      merit(link) = merit(link) + ((met_s(link)-spec(3))/sspec(3))**2
    ENDIF

! STELLAR SIZE NORMAL PRIOR
    IF(priorrad.eq.'p')THEN     
      merit(link) = merit(link) + ((radius_s(link)-starradius)/dstarradius_ini)**2
    ENDIF
    if (RjumpIso.eq.'y') then
      merit(link) = merit(link) + ((radius_s(link)-radius_Isoch)/dradius_Isoch)**2
    endif

! STELLAR DENSITY NORMAL PRIOR
    IF(priorrho.eq.'p')THEN     
      merit(link) = merit(link) + ((rho(link)-trho2)/strho2)**2
    ENDIF

! STELLAR LUMINOSITY NORMAL PRIOR
    IF(priorlum.eq.'p'.or.(priormag.eq.'p'.and.priordist.eq.'p'))THEN     
      merit(link) = merit(link) + ((lum_s(link)-starlum)/dstarlum)**2
    ENDIF

! VSINI NORMAL PRIOR
    IF(isjump(1,8).EQ.'y'.AND.nrv.GT.0.AND.priorvsini.eq.'p')THEN  
      merit(link) = merit(link) + ((vsini(link)-spec(5))/sspec(5))**2
    ENDIF 

! LOGG NORMAL PRIOR
    IF(ntr.GT.0.AND.priorlogg.eq.'p')THEN  
      merit(link) = merit(link) + ((logg(link)-spec(2))/sspec(2))**2
    ENDIF 
 
! F2 NORMAL PRIOR
    IF(testf2.EQ.'y'.AND.priorf2.eq.'p')THEN  
      merit(link) = merit(link) + ((f2(link)-spec(7))/sspec(7))**2
    ENDIF 

! JUMP PARAMETERS NORMAL PRIORS
    DO k=1,npla
      IF(isjump(k,1).EQ.'p'.AND.sdF(k).GT.1E-10)THEN
        merit(link) = merit(link) + ((dF(k,link)-dF_ini(k))/sdF_ini(k))**2
      ENDIF
      IF(isjump(k,2).EQ.'p'.AND.sb(k).GT.1E-10)THEN
        merit(link) = merit(link) + ((b(k,link)-b_ini(k))/sb_ini(k))**2
      ENDIF
      IF(isjump(k,3).EQ.'p'.AND.sdur(k).GT.1E-10)THEN
        merit(link) = merit(link) + ((dur(k,link)-dur_ini(k))/sdur_ini(k))**2
      ENDIF
      IF(isjump(k,4).EQ.'p'.AND.st0(k).GT.1E-10)THEN
        IF(link.EQ.1) ndata=ndata+1
        merit(link) = merit(link) + ((t0(k,link)-t0_ini(k))/st0_ini(k))**2 
      ENDIF
      IF(isjump(k,5).EQ.'p'.AND.sper(k).GT.1E-10)THEN
        merit(link) = merit(link) + ((per(k,link)-per_ini(k))/sper_ini(k))**2
      ENDIF
      IF(isjump(k,9).EQ.'p'.AND.stidel(k).GT.1E-10)THEN
        merit(link) = merit(link) + ((tidel(k,link)-tidel_ini(k))/stidel_ini(k))**2
      ENDIF
      IF(isjump(k,6).EQ.'p'.AND.ABS(dsecosw(k)).GT.1E-10.AND.ABS(dsesinw(k)).GT.1E-10)THEN
        merit(link) = merit(link) + ((secosw(k,link)-secosw_ini(k))/dsecosw_ini(k))**2
        merit(link) = merit(link) + ((sesinw(k,link)-sesinw_ini(k))/dsesinw_ini(k))**2
      ENDIF    
      IF(isjump(k,7).EQ.'p'.AND.dkb(k).GT.1E-10)THEN
        merit(link) = merit(link) + ((kb(k,link)-kb_ini(k))/dkb_ini(k))**2 
      ENDIF
    ENDDO
    IF(isjump(1,8).EQ.'p'.AND.ABS(dsvsinisinbeta).GT.1E-10.AND.ABS(dsvsinisinbeta).GT.1E-10)THEN
        merit(link)=merit(link)+((svsinicosbeta(link)-ini_svsinicosbeta)/dsvsinicosbeta_ini)**2
        merit(link)=merit(link)+((svsinisinbeta(link)-ini_svsinisinbeta)/dsvsinisinbeta_ini)**2
    ENDIF
    IF(limb.EQ.'qd'.AND.ntr.GT.0)THEN
      DO i=1,nfi
         IF(fitlimb(i).EQ.'p')THEN
           merit(link) = merit(link) + ((ql(i,1,link)-ql_ini(i,1))/sql(i,1))**2
           merit(link) = merit(link) + ((ql(i,2,link)-ql_ini(i,2))/sql(i,2))**2
        ENDIF
      ENDDO
    ENDIF
    IF(ntr.GT.0)THEN
      DO i=1,nfi
        DO j=1,npla
          IF(iffitoc(i).EQ.'p')THEN
            merit(link) = merit(link) + ((dFsec(i,link,j)-dFsec_ini(i,j))/edFsec_ini(i,j))**2
          ENDIF
        ENDDO
      ENDDO
      DO i=1,nfi
        IF(iffitph(i).EQ.'p')THEN
          merit(link) = merit(link) + ((phampli1(i,link)-phampli1_ini(i))/ephampli1_ini(i))**2
          merit(link) = merit(link) + ((phampli2(i,link)-phampli2_ini(i))/ephampli2_ini(i))**2
          merit(link) = merit(link) + ((phampli3(i,link)-phampli3_ini(i))/ephampli3_ini(i))**2
          merit(link) = merit(link) + ((phoffset(i,link)-phoffset_ini(i))/ephoffset_ini(i))**2
        ENDIF
      ENDDO
    ENDIF
    remerit(link) = merit(link)/DBLE(dof)

    ok='n'
    IF(ntr.GT.0)THEN
      k=0  
      ave_presglo(link) = 0.
      sig_presglo(link) = 0.
      ave_presbinglo(link) = 0.
      sig_presbinglo(link) = 0.
      DO i=1,ntr   
        nb = np(i)
        DO j=1,nb
          k=k+1
          tbjd(k) = bjd(i,j)
          tilo(k) = filter(i)
          tmodel(k) = model(i,j)
          tmodel_tr(k) = model_tr(i,j)
          tphot(k) = phot(i,j) 
          tphotcor(k) = photcor(i,j)
          tresi(k) = resi(i,j)
          ave_presglo(link) = ave_presglo(link)+tresi(k)/DBLE(nphotot)
          terror(k) = error(i,j)
          terrorcor(k) = error(i,j)
        ENDDO
      ENDDO
      k=0
      DO i=1,ntr
        nb = nbin(i,1)
        DO j=1,nb
          k = k + 1 
          tresibin(k) = resibin(i,j,1)
          ave_presbinglo(link) = ave_presbinglo(link)+tresibin(k)/DBLE(nphototbin)
        ENDDO
      ENDDO
      k=0
      DO i=1,ntr
        nb = np(i)
        DO j=1,nb
          k=k+1
          sig_presglo(link) = sig_presglo(link)+(tresi(k)-ave_presglo(link))**2
        ENDDO
      ENDDO
      sig_presglo(link) = SQRT(sig_presglo(link)/DBLE(nphotot))
      k=0
      DO i=1,ntr
        nb = nbin(i,1)
        DO j=1,nb
          k = k + 1 
          sig_presbinglo(link) = sig_presbinglo(link)+(tresibin(k))**2
        ENDDO
      ENDDO
      sig_presbinglo(link) = SQRT(sig_presbinglo(link)/DBLE(nphototbin))
      beta_red_glo(link)=SQRT(DBLE(nphototbin-1)/DBLE(nphototbin))*SQRT(DBLE(nmeanbin))
      beta_red_glo(link)= beta_red_glo(link)*sig_presbinglo(link)/sig_presglo(link)
      DO i=1,nphotot
        test2 = (tbjd(i)-t0(1,link))/per(1,link)
        test = NINT(test2)
        fbjd(i) = (test2-test)*per(1,link)
        IF(fbjd(i).LT.-0.25*per(1,link)) fbjd(i)=fbjd(i)+per(1,link)
        IF(isttv.EQ.'y'.AND.nttvmax.GT.0)THEN       
          epochtime = NINT((tbjd(i)-t0(1,link))/per(1,link))
          DO l=1,nttv(1)
            IF(epochtime.EQ.epochtr(1,l))THEN
              fbjd(i)=fbjd(i)-ttv(1,l,link)
            ENDIF 
          ENDDO
        ENDIF
        tbjd(i) = fbjd(i)
      ENDDO
      CALL SORT(nphotot,fbjd)
      tflick=0
      DO i=1,nphotot
        test = 0
        j=0
        DO WHILE(test.LT.1)
          j=j+1
          IF((ABS(fbjd(i)-tbjd(j)).lt.1.E-8).AND.tflick(j).EQ.0)THEN
            fphot(i) = tphot(j)
            filo(i) = tilo(j)
            fphotcor(i)=tphotcor(j)
            ferror(i) = terror(j)
            ferrorcor(i) = terrorcor(j)
            fmodel(i) = tmodel(j)
            fmodel_tr(i) = tmodel_tr(j)
            fresi(i) = tresi(j)
            tflick(j) = tflick(j)+1
            test = 1
          ENDIF
          IF(j.GT.nphotot) test = 1
        ENDDO
      ENDDO
      numbin =  FLOOR((MAXVAL(fbjd)-MINVAL(fbjd))*24*60./(binsize))
      ! Binning
      ALLOCATE(nf(nfi+1))
      ALLOCATE(fbjd_bin(numbin,nfi+1)); ALLOCATE(fphot_bin(numbin,nfi+1))
      ALLOCATE(ferror_bin(numbin,nfi+1)); ALLOCATE(fresi_bin(numbin,nfi+1))
      ALLOCATE(fphotcor_bin(numbin,nfi+1));ALLOCATE(ferrorcor_bin(numbin,nfi+1))
      DO l=1,nfi+1
        as1 = fbjd(1) - 1.5*binsize/(24*60) 
        nf(l)=0
        DO i=1,numbin
          as1 = as1 + binsize/(24.*60.) 
          as2 = as1 + binsize/(24.*60.)
          k=0
          DO j=1,nphotot
            IF(fbjd(j).GT.as1.AND.fbjd(j).LE.as2)THEN
              IF(l.LE.nfi)THEN
                IF(filo(j).EQ.wfilter(l)) k=k+1
              ELSE
                k=k+1
              ENDIF
            ENDIF
          ENDDO
          IF(k.GT.0)THEN
            nf(l)=nf(l)+1
          ENDIF
        ENDDO
        as1 = fbjd(1) - 1.5*binsize/(24*60)
        test=0 
        DO i=1,numbin
          as1 = as1 + binsize/(24.*60.) 
          as2 = as1 + binsize/(24.*60.)
          k=0
          DO j=1,nphotot
            IF(fbjd(j).GT.as1.AND.fbjd(j).LE.as2)THEN
              IF(l.LE.nfi)THEN
                IF(filo(j).EQ.wfilter(l)) k=k+1
              ELSE             
                k=k+1
               ENDIF
            ENDIF
          ENDDO
          IF(k.GT.0)THEN
            test=test+1
            fbjd_bin(test,l)=0.
            fphot_bin(test,l)=0.
            fphotcor_bin(test,l)=0.
            ferror_bin(test,l)=0.
            ferrorcor_bin(test,l)=0.
            fresi_bin(test,l)=0.
            IF(k.EQ.1)THEN
              DO j=1,nphotot
                IF(fbjd(j).GT.as1.AND.fbjd(j).LE.as2)THEN
                  IF(l.GT.nfi.OR.filo(j).EQ.wfilter(l))THEN
                    fbjd_bin(test,l)=fbjd(j)
                    fphot_bin(test,l)=fphot(j)
                    fphotcor_bin(test,l)=fphotcor(j)
                    fresi_bin(test,l)=fresi(j)
                    ferror_bin(test,l)=ferror(j)
                    ferrorcor_bin(test,l)=ferror(j)
                  ENDIF
                ENDIF
              ENDDO
            ELSE
              DO j=1,nphotot
                IF(fbjd(j).GT.as1.AND.fbjd(j).LE.as2)THEN
                  IF(l.GT.nfi.OR.filo(j).EQ.wfilter(l))THEN
                    fbjd_bin(test,l)=fbjd_bin(test,l)+fbjd(j)/DBLE(k)
                    fphot_bin(test,l)=fphot_bin(test,l)+fphot(j)/DBLE(k)
                    fphotcor_bin(test,l)=fphotcor_bin(test,l)+fphotcor(j)/DBLE(k)
                    fresi_bin(test,l)=fresi_bin(test,l)+fresi(j)/DBLE(k)
                  ENDIF
                ENDIF
              ENDDO
              DO j=1,nphotot
                IF(fbjd(j).GT.as1.AND.fbjd(j).LE.as2)THEN
                  IF(l.GT.nfi.OR.filo(j).EQ.wfilter(l))THEN
                    ferror_bin(test,l)=ferror_bin(test,l)+ &
                      & (fphot(j)-fphot_bin(test,l))**2
                    ferrorcor_bin(test,l)=ferrorcor_bin(test,l)+ & 
                      & (fresi(j)-fresi_bin(test,l))**2
                   ENDIF
                ENDIF
              ENDDO
              ferror_bin(test,l)=SQRT(ferror_bin(test,l))/DBLE(k)
              ferrorcor_bin(test,l)=SQRT(ferrorcor_bin(test,l))/DBLE(k)
            ENDIF
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    IF(nrv.GT.0)THEN
      k=0
      ave_presglorv(link)=0.
      sig_presglorv(link)=0.
      resmeanerror=0.
      DO i=1,nrv
        nb = nprv(i)
        DO j=1,nb
          k=k+1
          trvbjd(k) = rvbjd(i,j)
          trvmodel(k) = rvmod(i,j)
          trvmodel2(k) = rvmod2(i,j)
          trv(k) = rv(i,j)
          trv2(k) = rv2(i,j)
          trvresi(k) = rvresi(i,j)
          ave_presglorv(link) = ave_presglorv(link)+trvresi(k)/DBLE(nrvtot)
          trverror(k) = rverror(i,j)
          resmeanerror = resmeanerror + trverror(k)/DBLE(nrvtot)
        ENDDO
      ENDDO
      k=0
      DO i=1,nrv
        nb = nprv(i)
        DO j=1,nb
          k=k+1
          sig_presglorv(link) = sig_presglorv(link)+(trvresi(k)-ave_presglorv(link))**2
        ENDDO
      ENDDO
      sig_presglorv(link) = SQRT(sig_presglorv(link)/DBLE(nrvtot))
      IF((sig_presglorv(link)**2 - resmeanerror**2).GT.0.)THEN
        jitterglo(link) = SQRT(sig_presglorv(link)**2 - resmeanerror**2)
      ELSE
        jitterglo(link) = 0.
      ENDIF   
      DO i=1,nrvtot
        test2 = (trvbjd(i)-t0(1,link))/per(1,link)
        test = NINT(test2)
        frvbjd(i) = (test2-test)*per(1,link)
        trvbjd(i) = frvbjd(i)
      ENDDO
      CALL SORT(nrvtot,frvbjd)
      DO i=1,nrvtot
        DO j=1,nrvtot
          IF(ABS(frvbjd(i)-trvbjd(j)).lt.0.000001)THEN
            frv(i) = trv(j)
            frv2(i) = trv2(j)
            frverror(i) = trverror(j)
            frvmodel(i) = trvmodel(j)
            frvmodel2(i) = trvmodel2(j)
            frvresi(i) = trvresi(j)
          ENDIF
        ENDDO
      ENDDO 
    ENDIF

    IF(link.EQ.1)THEN
      PRINT*,'Initial deduced parameters'
      IF(limb.eq.'qd'.AND.ntr.GT.0)THEN
        DO i=1,nfi
          IF(fitlimb(i).EQ.'y')THEN
            temp = '  U1 ' // wfilter(i) // '-filter    '
            WRITE(*,108) temp, ql(i,1,link)
            temp = '  U2 ' // wfilter(i) // '-filter    '
            WRITE(*,108) temp, ql(i,2,link)
          ENDIF
        ENDDO
      ENDIF
      DO k=1,npla
        WRITE(*,108) ' e =                ',exc(k,link)
        WRITE(*,109) ' omega =            ',omega(k,link), '    deg'  
        WRITE(*,108) ' Rp/R* =            ',rr(k,link)
        WRITE(*,108) ' a/R* =             ',a_R(k,link)
        WRITE(*,108) ' b_tr =             ',b2(k,link)
        WRITE(*,108) ' b_oc =             ',b3(k,link)
        WRITE(*,109) ' K =                ',ka(k,link),'   ms-1'  
      ENDDO
      WRITE(*,109) ' rho* =             ',rho(link),'rho_sun'
      WRITE(*,'(A14,I3)') ' N_parameters= ',INT(npara)
      PRINT*, '--------------------------------------------'
    ENDIF
    

! METROPOLIS TEST  
    IF(MOD(link-1,nat).EQ.0.OR.link.EQ.1.AND.accepted(link).EQ.'y')THEN
      ok='y'                                           !If beginning of a new chain, accepted
    ELSE IF(accepted(link).EQ.'y')THEN
      thrown = expi**((merit(link-1)-merit(link))/2.)
      DO k=1,npla  
        IF(isjump(k,5).EQ.'j')THEN !Jeffrey's prior for P
          thrown = thrown * (per(k,link-1)/per(k,link))
        ENDIF
        IF(isjump(k,7).EQ.'j')THEN !Modified Jeffrey's prior for K
          thrown = thrown * ((ka(k,link-1)+0.1)/(ka(k,link)+0.1))
        ENDIF
!        if (isjump(k,3).ne.'n') then !Carter+2008eq36. Do not disfavor b=0, when stepping in W
!          thrown=thrown*sqrt(1-b(k,link)**2)/sqrt(1-b(k,link-1)**2)
!        endif
      ENDDO
      IF(thrown.GT.1.AND.accepted(link).EQ.'y')THEN
        ok='y'
      ELSE       
        CALL RANDOM_NUMBER(harvest)
        IF(harvest.LE.thrown) ok='y'
      ENDIF
      if(.not.(isoch.eq.'y'.and.ntr.eq.0.and.nrv.eq.0))then !not in case of the IsochPlacement only
		  DO k=1,npla
		    IF(radius_p(k,link).LT.rpdo) ok='n'
		    IF(radius_p(k,link).GT.rpup) ok='n'
		    IF(ironlow.EQ.'y'.AND.testiron(k).LT.1)THEN
		      ok='n'
		    ENDIF
		  ENDDO
	  end if
    ELSE
      ok='n'
    ENDIF
    IF(ok.EQ.'y')THEN
      accepted(link)='y'
    ELSE
      accepted(link)='n'
    ENDIF
    
!    if (mod(link,100).eq.0) print*,'radius',radius_s(link)
!    if (mod(link,100).eq.0) print*,'radius_Isoch',radius_Isoch
!    if (mod(link,100).eq.0) print*,'Inc radius_Isoch',dradius_Isoch
!    if (mod(link,100).eq.0) print*,'mass',mass_s(link)
!    if (mod(link,100).eq.0) print*,'mass_Isoch',mass_Isoch
!    if (mod(link,100).eq.0) print*,'Inc mass_Isoch',dmass_Isoch
!    if (mod(link,100).eq.0) print*,'teff',temp_s(link)
!    if (mod(link,100).eq.0) print*,'teff-1',temp_s(link-1)
!    if (mod(link,100).eq.0) print*,'FeH',met_s(link)
!    if (mod(link,100).eq.0) print*,'FeH-1',met_s(link-1)
!    if (mod(link,100).eq.0) print*,'MeritMass',((mass_s(link)-mass_Isoch)/dmass_Isoch)**2
!    if (mod(link,100).eq.0) print*,'MeritRadius',((radius_s(link)-radius_Isoch)/dradius_Isoch)**2
!    if (mod(link,100).eq.0) print*,'MeritTeff',((temp_s(link)-spec(1))/sspec(1))**2
!    if (mod(link,100).eq.0) print*,'MeritFeH',((met_s(link)-spec(3))/sspec(3))**2
!    if (mod(link,100).eq.0) print*,'Merit',merit(link)
!    if (mod(link,100).eq.0) print*,'Merit-1',merit(link-1)
!    if (mod(link,100).eq.0) print*,'accepted AFTER',accepted(link)
    

    IF(ok.EQ.'n')THEN    !If jump fails, go back to the previous values for all parameters
      DO k=1,npla
        prtr(k,link)=prtr(k,link-1)
        proc(k,link)=proc(k,link-1)
        rhop(k,link)=rhop(k,link-1)
        dF(k,link)=dF(k,link-1)
        b(k,link)=b(k,link-1)
        b2(k,link)=b2(k,link-1)
        b3(k,link)=b3(k,link-1)
        roche(k,link)=roche(k,link-1)
        a_roche(k,link)=a_roche(k,link-1)
        dur(k,link)=dur(k,link-1)
        t0(k,link)=t0(k,link-1)
        octime(k,link)=octime(k,link-1)
        semi(k,link)=semi(k,link-1)
        irrad(k,link)=irrad(k,link-1)
        inclian(k,link)=inclian(k,link-1)
        logg_p(k,link)=logg_p(k,link-1)
        teq_p(k,link)=teq_p(k,link-1)
        hillrad_p(k,link)=hillrad_p(k,link-1)
        per(k,link)=per(k,link-1)
        tidel(k,link)=tidel(k,link-1)
        secosw(k,link)=secosw(k,link-1)
        sesinw(k,link)=sesinw(k,link-1)
        ecosw(k,link)=ecosw(k,link-1)
        esinw(k,link)=esinw(k,link-1)
        exc(k,link)=exc(k,link-1)
        rr(k,link)=rr(k,link-1)
        a_R(k,link)=a_R(k,link-1)
        omega(k,link)=omega(k,link-1)
        mass_p(k,link)=mass_p(k,link-1)
        mass_p_sini(k,link)=mass_p_sini(k,link-1)
        radius_p(k,link)=radius_p(k,link-1)
        kb(k,link)=kb(k,link-1)
        ka(k,link)=ka(k,link-1)
        safro(k,link)=safro(k,link-1)
      ENDDO
      f2(link)=f2(link-1) 
      vsini(link)=vsini(link-1)
      beta(link)=beta(link-1)
      svsinisinbeta(link)=svsinisinbeta(link-1)
      svsinicosbeta(link)=svsinicosbeta(link-1)
      met_s(link)=met_s(link-1)
      temp_s(link)=temp_s(link-1)
      if (priorteff.ne.'p') col_s(link)=col_s(link-1)
      rho(link)=rho(link-1)
      mass_s(link)=mass_s(link-1) 
      radius_s(link)=radius_s(link-1)
      lum_s(link)=lum_s(link-1)
      !!
      logg(link)=logg(link-1)
      if (isoch.eq.'y') age(link)=age(link-1)
      !!
      ktide(link)=ktide(link-1)
      IF(stelincli.EQ.'y')THEN
        incli_s(link)=incli_s(link-1)
        sinincli_s(link)=sinincli_s(link-1)
      ENDIF
      merit(link)=merit(link-1)
      photmerit(link)=photmerit(link-1)
      rvmerit(link)=rvmerit(link-1)
      remerit(link)=remerit(link-1)
      IF(ntr.GT.0)THEN
        DO i=1,nfi
          DO j=1,npla
            dFsec(i,link,j)=dFsec(i,link-1,j)
          ENDDO
          phampli1(i,link)=phampli1(i,link-1)
          phampli2(i,link)=phampli2(i,link-1)
          phampli3(i,link)=phampli3(i,link-1)
          phoffset(i,link)=phoffset(i,link-1)
        ENDDO
      ENDIF 
      IF(ntr.GT.0.AND.nddf.GT.0.AND.isddf.NE.'n')THEN
        DO k=1,npla
          DO i=1,nddf
            ddf(k,i,link)=ddf(k,i,link-1)
            radipla(k,i,link)=radipla(k,i,link-1)
            dratio(k,i,link)=dratio(k,i,link-1)
            ddepth(k,i,link)=ddepth(k,i,link-1)
          ENDDO
        ENDDO
      ENDIF  
      IF(ntr.GT.0.AND.testsin.GT.0)THEN
        DO i=1,ntr
          IF(sinusnumber(i).GT.0)THEN
            k=sinusnumber(i)
            DO j=1,k
              sip(i,j,link)=sip(i,j,link-1)
              sit0(i,j,link)=sit0(i,j,link-1)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      IF(ntr.GT.0.AND.rampmod.EQ.'exp')THEN
        DO i=1,ntr
          IF(ramporder(i).GT.0)THEN
            t1ramp(i,link)=t1ramp(i,link-1)
            IF(ramporder(i).GT.1) t2ramp(i,link)=t2ramp(i,link-1)
          ENDIF
        ENDDO
      ENDIF
      IF(ntr.GT.0.AND.testflare.GT.0)THEN
        DO i=1,ntr
          IF(flarenumber(i).GT.0)THEN
            k=flarenumber(i)
            DO j=1,k
              flampli(i,j,link)=flampli(i,j,link-1)
              fltau(i,j,link)=fltau(i,j,link-1)
              flt0(i,j,link)=flt0(i,j,link-1)
            ENDDO
          ENDIF
        ENDDO
      ENDIF
      IF(nttvmax.GT.0.AND.isttv.EQ.'y')THEN    
        DO i=1,npla
          DO j=1,nttv(i)
            ttv(i,j,link)=ttv(i,j,link-1)
            ttr(i,j,link)=ttr(i,j,link-1)    
          ENDDO
        ENDDO
      ENDIF
      IF(ntr.GT.0)THEN
        IF(limb.EQ.'qd'.OR.limb.EQ.'no')THEN
          DO i=1,nfi
            IF(fitlimb(i).NE.'n')THEN
              DO j=1,2
                jumplimb(i,j,link)= jumplimb(i,j,link-1)
              ENDDO
            ENDIF
          ENDDO
        ENDIF
      ENDIF
    ENDIF

    IF(link.EQ.1)THEN
      PRINT*,      ' STEP          CHI2            RVCHI2         PHOTCHI2          BIC'
      WRITE(444,*) ' STEP          CHI2            RVCHI2         PHOTCHI2          BIC'
      WRITE(445,*) ' STEP          CHI2            RVCHI2         PHOTCHI2          BIC'
    ENDIF

    bf_1=-0.5*merit(link)
    IF(burned(link).EQ.'n') dic = dic-2*bf_1/DBLE((nat-nburn)*nchain)
    IF(merit(link).LT.bestmerit)THEN   
      bic = -2*bf_1+npara*LOG(ndata)
      WRITE(*,'(I7,4(1x,F15.4))') link,merit(link),rvmerit(link),photmerit(link),bic
      bf_1b = bf_1
      soluce=link
      bestmerit=merit(link)
      dic2 = - 2*bf_1
      aic = -2*bf_1+2*npara
      IF(ntr.GT.0)THEN 
        IF(nsysglo.GT.0) sysva_bestfit=sysva
        best_beta_red=beta_red
        best_bred_dur = bred_dur
        best_bnpave = bnpave
        best_bresrmsbin = bresrmsbin
        bf_resrmsbin = resrmsbin
        bf_resrms = resrms(:,link)
      ENDIF
      IF(nrv.GT.0)THEN
        rvsysva_bestfit=rvsysva
        best_jitter=jitter
        rvresrms_soluce=rvresrms
      ENDIF
      IF(ntiming.GT.0.AND.istiming.EQ.'y')THEN
        OPEN(UNIT=3,FILE='timing.res')
        ave = 0.
        DO i=1,ntiming
          test2 = (timing(i)-t0(1,link))/per(1,link)
          epoch(i) = NINT(test2)
          ttiming = t0(1,link) + epoch(i)*per(1,link)
          omctiming(i) = (timing(i) - ttiming)*24.*60.
          ave = ave + omctiming(i)/DBLE(ntiming)   
          stiming(i)=stiming(i)*24.*60.
          WRITE(3,*) epoch(i), omctiming(i), stiming(i)
        ENDDO
        CLOSE(3)
        fil2 = 'timing.res'; graphfile = 'timing.sm'
        min_x = MINVAL(epoch)-(MAXVAL(epoch)-MINVAL(epoch))*0.1
        max_x = MAXVAL(epoch)+(MAXVAL(epoch)-MINVAL(epoch))*0.1
        min_y = MINVAL(omctiming)-(MAXVAL(omctiming)-MINVAL(omctiming))*0.1-MAXVAL(stiming)
        max_y = MAXVAL(omctiming)+(MAXVAL(omctiming)-MINVAL(omctiming))*0.1+MAXVAL(stiming)
        ncol_file=3;abs_val=1;ord_val=2
        xlabel2='        EPOCH       ';ylabel='     O-C (min)      '
        abs_errorm=0;abs_errorp=0;ord_errorm=3;ord_errorp=3
        errory='y';errorx='n';cone='n';newf='y'
        CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
         & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
         & ord_errorp,1,cone,newf,2)
        DO i=1,ntiming
          stiming(i)=stiming(i)/(24*60.)
        ENDDO
      ENDIF 

      IF(ntr.GT.0)THEN
        DO i=1,ntr  
          let1 = CHAR(48+INT(i/1000))
          i2 = i - int(i/1000)*1000
          let2 = CHAR(48+INT(i2/100))
          i2 = i - int(i/100)*100
          let3 = CHAR(48+INT(i2/10))
          let4 = CHAR(48+MOD(i2,10))
          fin =  let1 // let2 // let3 // let4 // '.res'
          name = 'phot' // fin
          OPEN(UNIT=3,FILE=name)
          nb=np(i)
          DO j=1,nb
            test2 = (bjd(i,j)-t0(1,link))/per(1,link)
            test2b = (bjd(i,j)-t0(1,link))/per(1,link)
            test = NINT(test2)
            fobjd(i,j) = (test2b-test)*per(1,link)
            IF(fobjd(i,j).LT.-0.25*per(1,link))THEN
              fobjd(i,j)=fobjd(i,j)+per(1,link)
            ENDIF
            WRITE(3,111) fobjd(i,j),bjd(i,j),phot(i,j),photcor(i,j), &
             & error(i,j),resi(i,j),model(i,j),model_tr(i,j), &
             & model(i,j)-model(i,j)
          ENDDO
          CLOSE(3)
          name4 = 'phpt' // fin
          OPEN(UNIT=3,FILE=name4)
          nb=np(i)
          DO j=1,nb
            test2 = (bjd(i,j)-t0(1,link))/per(1,link)
            test2b = (bjd(i,j)-t0(1,link))/per(1,link)
            test = NINT(test2)
            fobjd(i,j) = (test2b-test)*per(1,link)
            IF(fobjd(i,j).LT.-0.25*per(1,link))THEN
              fobjd(i,j)=fobjd(i,j)+per(1,link)
            ENDIF
            WRITE(3,111) fobjd(i,j),bjd(i,j),phot(i,j),photcor2(i,j), &
             & error(i,j),resi(i,j),model(i,j),model_tr2(i,j), &
             & model(i,j)-model(i,j)
          ENDDO
          CLOSE(3)
          name2 = 'phob' // fin
          OPEN(UNIT=3,FILE=name2)
          IF(link.EQ.1) OPEN(UNIT=333,FILE='phot'//let1//let2//let3//let4//'.bin')
          nb2=nbin(i,1)
          DO j=1,nb2
            test2 = (bjdbin(i,j)-t0(1,link))/per(1,link)
            test2b = (bjdbin(i,j)-t0(1,link))/per(1,link)
            test = NINT(test2)
            fobjdbin(i,j) = (test2b-test)*per(1,link)
            IF(fobjdbin(i,j).LT.-0.25*per(1,link))THEN
              fobjdbin(i,j)=fobjdbin(i,j)+per(1,link)
            ENDIF
            WRITE(3,122) fobjdbin(i,j),bjdbin(i,j),photbin(i,j),ephotbin(i,j), &
              & photcorbin(i,j),ephotcorbin(i,j),resibin(i,j,1),eresibin(i,j), &
              & modelbin(i,j),model_trbin(i,j),dXbin(i,j),dYbin(i,j),fwhmbin(i,j),skybin(i,j)
            IF(link.EQ.1)THEN
              WRITE(333,*) bjdbin(i,j),photbin(i,j),rephotbin(i,j),dXbin(i,j),dYbin(i,j),fwhmbin(i,j), &
              & fwhmxbin(i,j),fwhmybin(i,j),skybin(i,j),airmassbin(i,j),binsize*60
            ENDIF
          ENDDO
          CLOSE(3)
          IF(link.EQ.1) CLOSE(333)
          name5 = 'phpb' // fin
          OPEN(UNIT=3,FILE=name5)
          DO j=1,nb2
            WRITE(3,122) fobjdbin(i,j),bjdbin(i,j),photbin(i,j),ephotbin(i,j), &
              & photcor2bin(i,j),ephotcor2bin(i,j),resibin(i,j,1),eresibin(i,j), &
              & modelbin(i,j),model_trbin(i,j),dXbin(i,j),dYbin(i,j),fwhmbin(i,j),skybin(i,j)
          ENDDO
          CLOSE(3)
          ALLOCATE(temp1D(ne3))
          DO l=1,ne3
            temp1D(l)=bf_resrmsbin(i,l)*SQRT(npave(i,l))/bf_resrms(i)
          ENDDO
          name3 = 'rms'//fin
          OPEN(UNIT=3,FILE=name3)
          DO l=2,ne3
            IF(l.EQ.2.AND.red_dur(l).GE.red_dur(1))THEN
              WRITE(3,*) red_dur(1)*24*60., 1E6*bf_resrms(i)/SQRT(npave(i,1)),1E6*bf_resrmsbin(i,1), &
               & temp1D(1), ' 1. '
              WRITE(3,*) red_dur(l)*24*60., 1E6*bf_resrms(i)/SQRT(npave(i,l)),1E6*bf_resrmsbin(i,l), &
               & temp1D(l), ' 1. '
            ELSE IF(l.LT.ne3.AND.red_dur(l).LE.red_dur(1).AND.red_dur(l+1).GT.red_dur(1))THEN 
              WRITE(3,*) red_dur(l)*24*60., 1E6*bf_resrms(i)/SQRT(npave(i,l)),1E6*bf_resrmsbin(i,l), &
               & temp1D(l), ' 1. '
              WRITE(3,*) red_dur(1)*24*60., 1E6*bf_resrms(i)/SQRT(npave(i,1)),1E6*bf_resrmsbin(i,1), &
               & temp1D(1), ' 1. '
            ELSE
              WRITE(3,*) red_dur(l)*24*60.,1E6*bf_resrms(i)/SQRT(npave(i,l)),1E6*bf_resrmsbin(i,l), & 
                & temp1D(l), ' 1. '
            ENDIF
          ENDDO
          CLOSE(3)

          fil2 = name3; graphfile = 'rms' //  let1 // let2 // let3 // let4 // '.sm'
          min_x = 0.
          max_x = MAXVAL(red_dur)*1.05*24*60.
          min_y = MINVAL(bf_resrmsbin(i,1:ne3))-(MAXVAL(bf_resrmsbin(i,1:ne3)) - &
           & MINVAL(bf_resrmsbin(i,1:ne3)))*0.1
          max_y = MAXVAL(bf_resrmsbin(i,1:ne3))+(MAXVAL(bf_resrmsbin(i,1:ne3)) - &
            & MINVAL(bf_resrmsbin(i,1:ne3)))*0.1
          min_y = min_y*1E6;max_y=max_y*1E6
          ncol_file=5;abs_val=1;ord_val=3
          xlabel2='';ylabel=''
          abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
          errory='n';errorx='n';cone='n';newf='y'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,0)
          abs_val=1;ord_val=2
          xlabel2='    Binning [min]   ';ylabel='      rms [ppm]     '
          abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
          errory='n';errorx='n';cone='y';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,2)
          min_x = 0.
          max_x = MAXVAL(red_dur)*1.05*24*60.
          min_y = MINVAL(temp1D)-(MAXVAL(temp1D)-MINVAL(temp1D))*0.1
          max_y = MAXVAL(temp1D)+(MAXVAL(temp1D)-MINVAL(temp1D))*0.1
          abs_val=1;ord_val=4
          xlabel2='';ylabel=''
          abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
          errory='n';errorx='n';cone='n';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,0)
          abs_val=1;ord_val=5
          xlabel2='    Binning [min]   ';ylabel='        Beta        '
          abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
          errory='n';errorx='n';cone='y';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,0)
          DEALLOCATE(temp1D)

          fil2 = name; graphfile = 'phot' //  let1 // let2 // let3 // let4 // '.sm'
          min_x = MINVAL(bjd(i,1:nb))-(MAXVAL(bjd(i,1:nb))-MINVAL(bjd(i,1:nb)))*0.05
          max_x = MAXVAL(bjd(i,1:nb))+(MAXVAL(bjd(i,1:nb))-MINVAL(bjd(i,1:nb)))*0.05
          min_y = MINVAL(phot(i,1:nb))-(MAXVAL(phot(i,1:nb))-MINVAL(phot(i,1:nb)))*0.1
          max_y = MAXVAL(phot(i,1:nb))+(MAXVAL(phot(i,1:nb))-MINVAL(phot(i,1:nb)))*0.1
          ncol_file=9;abs_val=2;ord_val=3
          xlabel2='';ylabel=''
          abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
          errory='n';errorx='n';cone='n';newf='y'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,1)
          abs_val=2;ord_val=7
          xlabel2='    JD - 2450000    ';ylabel='        FLUX       '
          errory='n';errorx='n';cone='y';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,0)
          fil2 = name2
          ncol_file=14;abs_val=2;ord_val=3
          xlabel2='';ylabel=''
          abs_errorm=0;abs_errorp=0;ord_errorm=4;ord_errorp=4
          errory='y';errorx='n';cone='n';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,2)
          fil2 = name
          min_y = MINVAL(resi(i,1:nb))-(MAXVAL(resi(i,1:nb))-MINVAL(resi(i,1:nb)))*0.1
          max_y = MAXVAL(resi(i,1:nb))+(MAXVAL(resi(i,1:nb))-MINVAL(resi(i,1:nb)))*0.1
          ncol_file=9;abs_val=2;ord_val=6
          xlabel2='    JD - 2450000    ';ylabel='        O-C         '
          errory='n';errorx='n';cone='n';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,1)
          abs_val=2;ord_val=9
          xlabel2='';ylabel=''
          errory='n';errorx='n';cone='y';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,0)
          fil2 = name2
          ncol_file=14;abs_val=2;ord_val=7
          xlabel2='';ylabel=''
          abs_errorm=0;abs_errorp=0;ord_errorm=8;ord_errorp=8
          errory='y';errorx='n';cone='n';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,2)

          fil2 = name; graphfile = 'phoc' //  let1 // let2 // let3 // let4 // '.sm'
          min_x = MINVAL(bjd(i,1:nb))-(MAXVAL(bjd(i,1:nb))-MINVAL(bjd(i,1:nb)))*0.05
          max_x = MAXVAL(bjd(i,1:nb))+(MAXVAL(bjd(i,1:nb))-MINVAL(bjd(i,1:nb)))*0.05
          min_y = MINVAL(photcor(i,1:nb))-(MAXVAL(photcor(i,1:nb))-MINVAL(photcor(i,1:nb)))*0.1
          max_y = MAXVAL(photcor(i,1:nb))+(MAXVAL(photcor(i,1:nb))-MINVAL(photcor(i,1:nb)))*0.1
          ncol_file=9;abs_val=2;ord_val=4
          xlabel2='';ylabel=''
          abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
          errory='n';errorx='n';cone='n';newf='y'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,1)
          abs_val=2;ord_val=8
          xlabel2='    JD - 2450000    ';ylabel='        FLUX       '
           errory='n';errorx='n';cone='y';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,0)
          fil2 = name2
          ncol_file=14;abs_val=2;ord_val=5
          xlabel2='';ylabel=''
          abs_errorm=0;abs_errorp=0;ord_errorm=6;ord_errorp=6
          errory='y';errorx='n';cone='n';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,2)
          fil2 = name
          min_y = MINVAL(resi(i,1:nb))-(MAXVAL(resi(i,1:nb))-MINVAL(resi(i,1:nb)))*0.1
          max_y = MAXVAL(resi(i,1:nb))+(MAXVAL(resi(i,1:nb))-MINVAL(resi(i,1:nb)))*0.1
          ncol_file=9;abs_val=2;ord_val=6
          xlabel2='    JD - 2450000    ';ylabel='        O-C         '
          errory='n';errorx='n';cone='n';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,1)
          abs_val=2;ord_val=9
          xlabel2='';ylabel=''
          errory='n';errorx='n';cone='y';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,0)
          fil2 = name2
          ncol_file=14;abs_val=2;ord_val=7
          xlabel2='';ylabel=''
          abs_errorm=0;abs_errorp=0;ord_errorm=8;ord_errorp=8
          errory='y';errorx='n';cone='n';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,2)

          fil2 = name4; graphfile = 'phop' //  let1 // let2 // let3 // let4 // '.sm'
          min_x = MINVAL(bjd(i,1:nb))-(MAXVAL(bjd(i,1:nb))-MINVAL(bjd(i,1:nb)))*0.05
          max_x = MAXVAL(bjd(i,1:nb))+(MAXVAL(bjd(i,1:nb))-MINVAL(bjd(i,1:nb)))*0.05
          min_y = MINVAL(photcor2(i,1:nb))-(MAXVAL(photcor2(i,1:nb))-MINVAL(photcor2(i,1:nb)))*0.1
          max_y = MAXVAL(photcor2(i,1:nb))+(MAXVAL(photcor2(i,1:nb))-MINVAL(photcor2(i,1:nb)))*0.1
          ncol_file=9;abs_val=2;ord_val=4
          xlabel2='';ylabel=''
          abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
          errory='n';errorx='n';cone='n';newf='y'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,1)
          abs_val=2;ord_val=8
          xlabel2='    JD - 2450000    ';ylabel='        FLUX       '
           errory='n';errorx='n';cone='y';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,0)
          fil2 = name5
          ncol_file=14;abs_val=2;ord_val=5
          xlabel2='';ylabel=''
          abs_errorm=0;abs_errorp=0;ord_errorm=6;ord_errorp=6
          errory='y';errorx='n';cone='n';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,2)
          fil2 = name
          min_y = MINVAL(resi(i,1:nb))-(MAXVAL(resi(i,1:nb))-MINVAL(resi(i,1:nb)))*0.1
          max_y = MAXVAL(resi(i,1:nb))+(MAXVAL(resi(i,1:nb))-MINVAL(resi(i,1:nb)))*0.1
          ncol_file=9;abs_val=2;ord_val=6
          xlabel2='    JD - 2450000    ';ylabel='        O-C         '
          errory='n';errorx='n';cone='n';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,1)
          abs_val=2;ord_val=9
          xlabel2='';ylabel=''
          errory='n';errorx='n';cone='y';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,0)
          fil2 = name2
          ncol_file=14;abs_val=2;ord_val=7
          xlabel2='';ylabel=''
          abs_errorm=0;abs_errorp=0;ord_errorm=8;ord_errorp=8
          errory='y';errorx='n';cone='n';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,2)

          fil2 = name; graphfile = 'fophot' //  let1 // let2 // let3 // let4 // '.sm'
          min_x = MINVAL(fobjd(i,1:nb))-(MAXVAL(fobjd(i,1:nb))-MINVAL(fobjd(i,1:nb)))*0.05
          max_x = MAXVAL(fobjd(i,1:nb))+(MAXVAL(fobjd(i,1:nb))-MINVAL(fobjd(i,1:nb)))*0.05
          min_y = MINVAL(photcor(i,1:nb))-(MAXVAL(photcor(i,1:nb))-MINVAL(photcor(i,1:nb)))*0.1
          max_y = MAXVAL(photcor(i,1:nb))+(MAXVAL(photcor(i,1:nb))-MINVAL(photcor(i,1:nb)))*0.1
          ncol_file=9;abs_val=1;ord_val=4
          xlabel2='';ylabel=''
          abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
          errory='n';errorx='n';cone='n';newf='y'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,1)
          abs_val=1;ord_val=8
          xlabel2='       dT (d)       ';ylabel='       FLUX        '
           errory='n';errorx='n';cone='y';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,0)
          fil2 = name2
          ncol_file=14;abs_val=1;ord_val=5
          xlabel2='';ylabel=''
          abs_errorm=0;abs_errorp=0;ord_errorm=6;ord_errorp=6
          errory='y';errorx='n';cone='n';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,2)
          fil2 = name
          min_y = MINVAL(resi(i,1:nb))-(MAXVAL(resi(i,1:nb))-MINVAL(resi(i,1:nb)))*0.1
          max_y = MAXVAL(resi(i,1:nb))+(MAXVAL(resi(i,1:nb))-MINVAL(resi(i,1:nb)))*0.1
          ncol_file=9;abs_val=1;ord_val=6
          xlabel2='       dT (d)       ';ylabel='        O-C         '
          errory='n';errorx='n';cone='n';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,1)
          abs_val=1;ord_val=9
          xlabel2='';ylabel=''
          errory='n';errorx='n';cone='y';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,0)
          fil2 = name2
          ncol_file=14;abs_val=1;ord_val=7
          xlabel2='';ylabel=''
          abs_errorm=0;abs_errorp=0;ord_errorm=8;ord_errorp=8
          errory='y';errorx='n';cone='n';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,2)

        ENDDO

        DO k=1,nfi+1
          IF(k.LE.nfi)THEN
            OPEN(UNIT=3,FILE='mcmc_phot'//wfilter(k)//'.res')
          ELSE 
            OPEN(UNIT=3,FILE='mcmc_phot.res')
          ENDIF
          DO i=1,nphotot
            IF(filo(i).EQ.wfilter(k).OR.k.EQ.nfi+1)THEN
              WRITE(3,110) i,fbjd(i),fphot(i),fphotcor(i),ferror(i), &
               & fresi(i),fmodel_tr(i),fmodel(i)-fmodel(i)
            ENDIF
          ENDDO
          CLOSE(3)
          IF(k.LE.nfi)THEN
            OPEN(UNIT=3,FILE='mcmc_phob'//wfilter(k)//'.res') 
          ELSE
            OPEN(UNIT=3,FILE='mcmc_phob.res') 
          ENDIF
          n2 = nf(k)
          DO i=1,n2
            WRITE(3,*) i,fbjd_bin(i,k),fphot_bin(i,k),fphotcor_bin(i,k),ferrorcor_bin(i,k), &
             & fresi_bin(i,k)
          ENDDO
          CLOSE(3)
          IF(k.LE.nfi)THEN
            fil2 = 'mcmc_phot'//wfilter(k)//'.res';graphfile = 'mcmc_phot'//wfilter(k)//'.sm'
          ELSE
            fil2 = 'mcmc_phot.res';graphfile = 'mcmc_phot.sm'
          ENDIF
          min_x = MINVAL(fbjd)-(MAXVAL(fbjd)-MINVAL(fbjd))*0.05
          max_x = MAXVAL(fbjd)+(MAXVAL(fbjd)-MINVAL(fbjd))*0.05
          min_y = MINVAL(fphotcor)-(MAXVAL(fphotcor)-MINVAL(fphotcor))*0.1
          max_y = MAXVAL(fphotcor)+(MAXVAL(fphotcor)-MINVAL(fphotcor))*0.1
          ncol_file=8;abs_val=2;ord_val=4
          xlabel2='';ylabel=''
          abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
          errory='n';errorx='n';cone='n';newf='y'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,1)
          abs_val=2;ord_val=7
          xlabel2='       dT (d)       ';ylabel='        FLUX        '
          errory='n';errorx='n';cone='y';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,0)
          fil2 = 'mcmc_phob.res' 
          IF(k.LE.nfi) fil2 = 'mcmc_phob'//wfilter(k)//'.res'  
          ncol_file=6;abs_val=2;ord_val=4
          xlabel2='';ylabel=''
          abs_errorm=0;abs_errorp=0;ord_errorm=5;ord_errorp=5
          errory='y';errorx='n';cone='n';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,2)
          fil2 = 'mcmc_phot.res'
          IF(k.LE.nfi) fil2 = 'mcmc_phot'//wfilter(k)//'.res'
          min_y = MINVAL(fresi)-(MAXVAL(fresi)-MINVAL(fresi))*0.1
          max_y = MAXVAL(fresi)+(MAXVAL(fresi)-MINVAL(fresi))*0.1
          ncol_file=8;abs_val=2;ord_val=6
          xlabel2='       dT (d)       ';ylabel='        O-C         '
          errory='n';errorx='n';cone='n';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,1)
          abs_val=2;ord_val=8
          xlabel2='';ylabel=''
          errory='n';errorx='n';cone='y';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,0)
          fil2 = 'mcmc_phob.res'
          IF(k.LE.nfi) fil2 = 'mcmc_phob'//wfilter(k)//'.res'
          ncol_file=6;abs_val=2;ord_val=6
          xlabel2='';ylabel=''
          abs_errorm=0;abs_errorp=0;ord_errorm=5;ord_errorp=5
          errory='y';errorx='n';cone='n';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,2)
         ENDDO
      ENDIF

      IF(nrv.GT.0)THEN
        OPEN(UNIT=21,FILE='mcmc_rv.res')
        DO i=1,nrvtot
          WRITE(21,121) i,frvbjd(i),frv2(i),frverror(i),frvresi(i),frvmodel2(i), &
           & frvresi(i)-frvresi(i)
        ENDDO
        CLOSE(21)
        fil2 = 'mcmc_rv.res';graphfile = 'mcmc_rv.sm'
        min_x = MINVAL(frvbjd)-(MAXVAL(frvbjd)-MINVAL(frvbjd))*0.05
        max_x = MAXVAL(frvbjd)+(MAXVAL(frvbjd)-MINVAL(frvbjd))*0.05
        min_y = MINVAL(frv2)-(MAXVAL(frv2)-MINVAL(frv2))*0.1
        max_y = MAXVAL(frv2)+(MAXVAL(frv2)-MINVAL(frv2))*0.1
        ncol_file=7;abs_val=2;ord_val=3
        xlabel2='';ylabel=''
        abs_errorm=0;abs_errorp=0;ord_errorm=4;ord_errorp=4
        errory='y';errorx='n';cone='n';newf='y'
        CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
         & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
         & ord_errorp,2,cone,newf,0)
        abs_val=2;ord_val=6
        xlabel2='       dT (d)       ';ylabel='      RV (m/s)      '
        errory='n';errorx='n';cone='y';newf='n'
        CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
         & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
         & ord_errorp,2,cone,newf,0)
        min_y = MINVAL(frvresi)-(MAXVAL(frvresi)-MINVAL(frvresi))
        max_y = MAXVAL(frvresi)+(MAXVAL(frvresi)-MINVAL(frvresi))
        abs_val=2;ord_val=5
        xlabel2='       dT (d)       ';ylabel='      O-C (m/s)     '
        errory='y';errorx='n';cone='n';newf='n'
        CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
         & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
         & ord_errorp,3,cone,newf,0)
        abs_val=2;ord_val=7
        xlabel2='';ylabel=''
        errory='n';errorx='n';cone='y';newf='n'
        CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
         & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
         & ord_errorp,3,cone,newf,0)
        DO i=1,nrv   
          let1 = CHAR(48+INT(i/1000))
          i2 = i - INT(i/1000)*1000
          let2 = CHAR(48+INT(i2/100))
          i2 = i - INT(i/100)*100
          let3 = CHAR(48+INT(i2/10))
          let4 = CHAR(48+MOD(i2,10))
          fil2 =  'rv' // let1 // let2 // let3 // let4 // '.res'
          OPEN(UNIT=3,FILE=fil2)
          nb=nprv(i)
          DO j=1,nb
            WRITE(3,121) j,rvbjd(i,j),rv2(i,j),rverror(i,j),rvresi(i,j), &
             & rvmod2(i,j),rvresi(i,j)-rvresi(i,j)
          ENDDO
          CLOSE(3)
          graphfile = 'rv' // let1 // let2 // let3 // let4 // '.sm'     
          min_x = MINVAL(rvbjd(i,1:nb))-(MAXVAL(rvbjd(i,1:nb))-MINVAL(rvbjd(i,1:nb)))*0.05
          max_x = MAXVAL(rvbjd(i,1:nb))+(MAXVAL(rvbjd(i,1:nb))-MINVAL(rvbjd(i,1:nb)))*0.05
          min_y = MINVAL(rv2(i,1:nb))-(MAXVAL(rv2(i,1:nb))-MINVAL(rv2(i,1:nb)))*0.1
          max_y = MAXVAL(rv2(i,1:nb))+(MAXVAL(rv2(i,1:nb))-MINVAL(rv2(i,1:nb)))*0.1
          ncol_file=7;abs_val=2;ord_val=3
          xlabel2='';ylabel=''
          abs_errorm=0;abs_errorp=0;ord_errorm=4;ord_errorp=4
          errory='y';errorx='n';cone='n';newf='y'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,0)
          abs_val=2;ord_val=6
          xlabel2='     JD-2450000     ';ylabel='      RV (m/s)      '
          errory='n';errorx='n';cone='y';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,2,cone,newf,0)
          min_y = MINVAL(rvresi(i,1:nb))-(MAXVAL(rvresi(i,1:nb))-MINVAL(rvresi(i,1:nb)))
          max_y = MAXVAL(rvresi(i,1:nb))+(MAXVAL(rvresi(i,1:nb))-MINVAL(rvresi(i,1:nb)))
          abs_val=2;ord_val=5
          xlabel2='     JD-245000     ';ylabel='      O-C (m/s)     '
          errory='y';errorx='n';cone='n';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,0)
          abs_val=2;ord_val=7
          xlabel2='';ylabel=''
          errory='n';errorx='n';cone='y';newf='n'
          CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,3,cone,newf,0)
          DO k=1,npla
            let5 = CHAR(48+INT(k/10))
            let6 = CHAR(48+MOD(k,10))
            fil2 =   'forv'// let1 // let2 // let3 // let4 // 'p' // let5 // let6 // '.res'
            OPEN(UNIT=3,FILE=fil2)
            nb=nprv(i)
            DO j=1,nb
              test2 = (rvbjd(i,j)-t0(k,link))/per(k,link)
              test = NINT(test2)
              forvbjd2(i,j) = (test2-test)*per(k,link)
              forvbjd(i,j) = forvbjd2(i,j)
            ENDDO
            CALL SORT(nb,forvbjd2(i,:))
            DO j=1,nb
              DO l=1,nb
                IF(ABS(forvbjd2(i,j)-forvbjd(i,l)).LT.1E-8)THEN
                  WRITE(3,121) l,forvbjd(i,l),rv2pla(k,i,l),rverror(i,l),rvresi(i,l), &
                   & rvmodpla(k,i,l),rvresi(i,l)-rvresi(i,l)
                ENDIF
              ENDDO
            ENDDO
            CLOSE(3)
            graphfile = 'forv'// let1 // let2 // let3 // let4 // 'p' // let5 // let6 // '.sm'
            min_x = MINVAL(forvbjd(i,1:nb))-(MAXVAL(forvbjd(i,1:nb))-MINVAL(forvbjd(i,1:nb)))*0.05
            max_x = MAXVAL(forvbjd(i,1:nb))+(MAXVAL(forvbjd(i,1:nb))-MINVAL(forvbjd(i,1:nb)))*0.05
            min_y = MINVAL(rv2pla(k,i,1:nb))-(MAXVAL(rv2pla(k,i,1:nb))-MINVAL(rv2pla(k,i,1:nb)))*0.1
            max_y = MAXVAL(rv2pla(k,i,1:nb))+(MAXVAL(rv2pla(k,i,1:nb))-MINVAL(rv2pla(k,i,1:nb)))*0.1
            ncol_file=7;abs_val=2;ord_val=3
            xlabel2='';ylabel=''
            abs_errorm=0;abs_errorp=0;ord_errorm=4;ord_errorp=4
            errory='y';errorx='n';cone='n';newf='y'
            CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
             & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
             & ord_errorp,2,cone,newf,0)
            abs_val=2;ord_val=6
            xlabel2='       dT (d)       ';ylabel='      RV (m/s)      '
            errory='n';errorx='n';cone='y';newf='n'
            CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
             & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
             & ord_errorp,2,cone,newf,0)
            min_y = MINVAL(rvresi(i,1:nb))-(MAXVAL(rvresi(i,1:nb))-MINVAL(rvresi(i,1:nb)))
            max_y = MAXVAL(rvresi(i,1:nb))+(MAXVAL(rvresi(i,1:nb))-MINVAL(rvresi(i,1:nb)))
            abs_val=2;ord_val=5
            xlabel2='       dT (d)       ';ylabel='      O-C (m/s)     '
            errory='y';errorx='n';cone='n';newf='n'
            CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
             & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
             & ord_errorp,3,cone,newf,0)
            abs_val=2;ord_val=7
            xlabel2='';ylabel=''
            errory='n';errorx='n';cone='y';newf='n'
            CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
             & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
             & ord_errorp,3,cone,newf,0)
          ENDDO
        ENDDO
      ENDIF

      IF(ntr.GT.0.AND.nrv.GT.0)THEN
        fil2 = 'mcmc_phot.res';graphfile = 'mcmc.sm'
        min_x = MINVAL(fbjd)-(MAXVAL(fbjd)-MINVAL(fbjd))*0.05
        max_x = MAXVAL(fbjd)+(MAXVAL(fbjd)-MINVAL(fbjd))*0.05
        min_y = MINVAL(fphotcor)-(MAXVAL(fphotcor)-MINVAL(fphotcor))*0.1
        max_y = MAXVAL(fphotcor)+(MAXVAL(fphotcor)-MINVAL(fphotcor))*0.1
        ncol_file=8;abs_val=2;ord_val=4
        xlabel2='';ylabel=''
        abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
        errory='n';errorx='n';cone='n';newf='y'
        CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
         & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
         & ord_errorp,2,cone,newf,1)
        abs_val=2;ord_val=7
        xlabel2='       dT (d)       ';ylabel='        FLUX        '
        errory='n';errorx='n';cone='y';newf='n'
        CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
         & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
         & ord_errorp,2,cone,newf,0)
        fil2 = 'mcmc_phob.res'
        ncol_file=6;abs_val=2;ord_val=4
        xlabel2='';ylabel=''
        abs_errorm=0;abs_errorp=0;ord_errorm=5;ord_errorp=5
        errory='y';errorx='n';cone='n';newf='n'
        CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
         & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
         & ord_errorp,2,cone,newf,2)
        fil2 = 'mcmc_rv.res'
        min_x = MINVAL(frvbjd)-(MAXVAL(frvbjd)-MINVAL(frvbjd))*0.05
        max_x = MAXVAL(frvbjd)+(MAXVAL(frvbjd)-MINVAL(frvbjd))*0.05
        min_y = MINVAL(frv2)-(MAXVAL(frv2)-MINVAL(frv2))*0.1
        max_y = MAXVAL(frv2)+(MAXVAL(frv2)-MINVAL(frv2))*0.1
        ncol_file=7;abs_val=2;ord_val=3
        xlabel2='';ylabel=''
        abs_errorm=0;abs_errorp=0;ord_errorm=4;ord_errorp=4
        errory='y';errorx='n';cone='n';newf='n'
        CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
         & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
         & ord_errorp,3,cone,newf,0)
        abs_val=2;ord_val=6
        xlabel2='       dT (d)       ';ylabel='      RV (m/s)      '
        errory='n';errorx='n';cone='y';newf='n'
        CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
         & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
         & ord_errorp,3,cone,newf,0)
      ENDIF
    ENDIF
    IF(ntr.GT.0)THEN
      DEALLOCATE(fbjd_bin); DEALLOCATE(fphot_bin)
      DEALLOCATE(ferror_bin); DEALLOCATE(fresi_bin)
      DEALLOCATE(ferrorcor_bin); DEALLOCATE(fphotcor_bin)
      DEALLOCATE(nf)
    ENDIF

    WRITE(51,*) link,merit(link),photmerit(link),rvmerit(link)

    IF(link.EQ.1)THEN
      PRINT*, 'STOP HERE?  [y/n]'
      READ(*,*) dum
      PRINT*, '--------------------------------------------' 
      IF(dum.EQ.'y') STOP
    ENDIF

! Computation of the standard deviation on the N last steps for all the parameters
    IF(gibbs_sampler.EQ.'n')THEN
      IF(MOD(link,statlen).EQ.0)THEN
        nacce = 0.
        DO i=link-statlen+1,link
          IF(accepted(i).EQ.'y')THEN
            nacce = nacce+1.
          ENDIF
        ENDDO      
        nacce = nacce/DBLE(statlen)
        IF(nacce.GT.0.)THEN
          regul = regul*nacce/success
        ELSE
          regul=regul*0.5
        ENDIF
        IF(regul.LT.(1./maxadapt)) regul=1./maxadapt
        IF(regul.GT.maxadapt) regul=maxadapt
      ENDIF
    ENDIF
    IF(gibbs_sampler.EQ.'y'.AND.MOD(link,nat).GT.nburn2)THEN
      IF(MOD(link-nburn2,statlen).EQ.0)THEN
        nacce = 0.
        DO i=link-statlen+1,link
          IF(accepted(i).EQ.'y')THEN
            nacce = nacce+1.
          ENDIF
        ENDDO      
        nacce = nacce/DBLE(statlen)
        IF(nacce.GT.0.)THEN
          regul = regul*nacce/success
        ELSE
          regul=regul*0.5
        ENDIF
      ENDIF
    ENDIF 
    IF(gibbs_sampler.EQ.'y'.AND.MOD(link,nat).LE.nburn2.AND.MOD(link,nat).GT.1 &
     &.AND.MOD(MOD(link,nat)-1,statlen2*njump).EQ.0)THEN
      DO i=1,njump
        nacce=0.
        DO j=link-(statlen2*njump)+1,link
          IF(accepted(j).EQ.'y'.AND.jumpgibbs(j).EQ.i)THEN
            nacce = nacce+1.
          ENDIF
        ENDDO 
        nacce = nacce/DBLE(statlen2)        
        IF(nacce.GT.0.)THEN
          regul2(i) = regul2(i)*nacce/0.4
        ELSE
          regul2(i)=regul2(i)*0.5
        ENDIF
        IF(regul2(i).LT.(1./maxadapt)) regul2(i)=1./maxadapt
        IF(regul2(i).GT.maxadapt) regul2(i)=maxadapt
      ENDDO
    ENDIF
  ENDDO
  CLOSE(51)

  ! Figure of the evolution of the merit function
  fil2 = 'mcmc_chi2.res';graphfile = 'mcmc_chi2.sm'
  min_x = 0; max_x = na + 1
  min_y = MINVAL(merit)-(MAXVAL(merit)-MINVAL(merit))*0.1
  max_y = MAXVAL(merit)+(MAXVAL(merit)-MINVAL(merit))*0.1
  ncol_file=4;abs_val=1;ord_val=2
  xlabel2='   Merit function   ';ylabel='        Step        '
  abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
  errory='n';errorx='n';cone='n';newf='y'
  CALL GRAPH2D(fil2,graphfile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
   & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
   & ord_errorp,1,cone,newf,2)
   
   !Writing I/O consistency in models
   if (isoch.eq.'y') then
   	write(*,'(A30,F4.2)') 'IO consistency rate in models ',dble(rowTot)/dble(na)
   	open(unit=7,file='rowIO.res')
   	do j=1,na
   		write(7,*) rowIO(j)
   	end do
   	close(7)
   end if

!**********************
!11. FINAL STATISTICS *
!**********************

  IF(gelman.EQ.'y'.AND.nchain.GT.1)THEN
    ALLOCATE(avusti(nchain));ALLOCATE(sigusti(nchain));ALLOCATE(length(nchain))
    IF(fixstellar.NE.'y')THEN
      CALL gelmancomp(na,nchain,nburn,nat,temp_s,gelman_teff)
      if (priorcol.eq.'p') then
        call gelmancomp(na,nchain,nburn,nat,col_s,gelman_col)
      end if
      CALL gelmancomp(na,nchain,nburn,nat,met_s,gelman_met)
      IF(massfromr.NE.'y'.AND.enoch.NE.'y'.AND.msrs.NE.'y'.and.isoch.eq.'y')THEN !ne
        CALL gelmancomp(na,nchain,nburn,nat,mass_s,gelman_ms)
      ENDIF
      IF(fitmsrs.EQ.'y'.OR.massfromr.EQ.'y'.or.RjumpIso.eq.'y')THEN 
        CALL gelmancomp(na,nchain,nburn,nat,radius_s,gelman_rs)
      ENDIF
    ENDIF
    IF(limb.EQ.'qd'.AND.ntr.GT.0)THEN  
      DO j=1,nfi
        IF(fitlimb(j).NE.'n')THEN
          DO l=1,nd
            ALLOCATE(temp1D(na));temp1D=jumplimb(j,l,:)
            CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)
            gelman_jumplimb(j,l)=gelmanval;DEALLOCATE(temp1D)
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    IF(ngroup.GT.0)THEN
      DO j=1,ngroup
        ALLOCATE(temp1D(na));temp1D=dfgroup(j,:)
        CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)
        gelman_dfgroup(j)=gelmanval;DEALLOCATE(temp1D)
      ENDDO
    ENDIF
    DO j=1,npla
      IF(isjump(j,1).NE.'n')THEN
        ALLOCATE(temp1D(na));temp1D=dF(j,:)
        CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)
        gelman_dF(j)=gelmanval;DEALLOCATE(temp1D)
      ENDIF
      IF(isjump(j,2).NE.'n')THEN
        ALLOCATE(temp1D(na));temp1D=b(j,:)
        CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
        gelman_b(j)=gelmanval;DEALLOCATE(temp1D)
      ENDIF
      IF(isjump(j,3).NE.'n')THEN
        ALLOCATE(temp1D(na));temp1D=dur(j,:)
        CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
        gelman_dur(j)=gelmanval;DEALLOCATE(temp1D)        
      ENDIF
      IF(isjump(j,4).NE.'n')THEN
        ALLOCATE(temp1D(na));temp1D=t0(j,:)
        CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
        gelman_t0(j)=gelmanval;DEALLOCATE(temp1D)     
      ENDIF
      IF(isjump(j,5).NE.'n')THEN
        ALLOCATE(temp1D(na));temp1D=per(j,:)
        CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
        gelman_per(j)=gelmanval;DEALLOCATE(temp1D)   
      ENDIF 
      IF(isjump(j,9).NE.'n')THEN
        ALLOCATE(temp1D(na));temp1D=tidel(j,:)
        CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
        gelman_tidel(j)=gelmanval;DEALLOCATE(temp1D)               
      ENDIF
      IF(isjump(j,6).NE.'n')THEN
        ALLOCATE(temp1D(na));temp1D=secosw(j,:)
        CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
        gelman_secosw(j)=gelmanval
        temp1D=sesinw(j,:)
        CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
        gelman_sesinw(j)=gelmanval;DEALLOCATE(temp1D)
      ENDIF
      IF(isjump(j,7).NE.'n')THEN
        ALLOCATE(temp1D(na));temp1D=kb(j,:)
        CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
        gelman_kb(j)=gelmanval;DEALLOCATE(temp1D)                      
      ENDIF
      IF(nddf.GT.0.AND.ntr.GT.0.AND.isddf.NE.'n')THEN
        DO l=1,nddf
          ALLOCATE(temp1D(na));temp1D=ddf(j,l,:)
          CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
          gelman_ddf(j,l)=gelmanval;DEALLOCATE(temp1D)              
        ENDDO
      ENDIF
      IF(nttvmax.GT.0.AND.isttv.EQ.'y')THEN
        DO k=1,npla
          DO l=1,nttv(k)
            ALLOCATE(temp1D(na));temp1D=ttv(k,l,:)
            CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
            gelman_ttv(k,l)=gelmanval;DEALLOCATE(temp1D)                      
          ENDDO
        ENDDO
      ENDIF
      IF(j.EQ.1)THEN
        IF(isjump(j,8).NE.'n')THEN
          ALLOCATE(temp1D(na));temp1D=svsinicosbeta
          CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
          gelman_svsinicosbeta=gelmanval  
          temp1D=svsinisinbeta
          CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
          gelman_svsinisinbeta=gelmanval;DEALLOCATE(temp1D)            
        ENDIF
      ENDIF  
      IF(ntr.GT.0)THEN
        DO l=1,nfi
          ALLOCATE(temp1D(na));temp1D=dFsec(l,:,j)
          CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
          gelman_dFsec(l,j)=gelmanval;DEALLOCATE(temp1D)
        ENDDO
      ENDIF
      IF(j.EQ.1.AND.ntr.GT.0)THEN
        DO l=1,nfi
          ALLOCATE(temp1D(na));temp1D=phampli1(l,:)
          CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
          gelman_phampli1(l)=gelmanval
          temp1D=phampli2(l,:)
          CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
          gelman_phampli2(l)=gelmanval
          temp1D=phampli3(l,:)
          CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
          gelman_phampli3(l)=gelmanval
          temp1D=phoffset(l,:)
          CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
          gelman_phoffset(l)=gelmanval;DEALLOCATE(temp1D)           
        ENDDO
      ENDIF
    ENDDO
    IF(testf2.EQ.'y'.AND.priorf2.NE.'n')THEN
      ALLOCATE(temp1D(na));temp1D=f2
      CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
      gelman_f2=gelmanval;DEALLOCATE(temp1D)              
    ENDIF
    IF(ntr.GT.0.AND.testsin.GT.0)THEN
      DO l=1,ntr
        IF(sinusnumber(l).GT.0)THEN
          di = sinusnumber(l)
          DO j=1,di
            ALLOCATE(temp1D(na));temp1D=sit0(l,j,:)
            CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
            gelman_sit0(l,j)=gelmanval;DEALLOCATE(temp1D) 
            ALLOCATE(temp1D(na));temp1D=sip(l,j,:)
            CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
            gelman_sip(l,j)=gelmanval;DEALLOCATE(temp1D) 
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    IF(ntr.GT.0.AND.testflare.GT.0)THEN
      DO l=1,ntr
        IF(flarenumber(l).GT.0)THEN
          di = flarenumber(l)
          DO j=1,di
            ALLOCATE(temp1D(na));temp1D=flampli(l,j,:)
            CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
            gelman_flampli(l,j)=gelmanval;DEALLOCATE(temp1D) 
            ALLOCATE(temp1D(na));temp1D=fltau(l,j,:)
            CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
            gelman_fltau(l,j)=gelmanval;DEALLOCATE(temp1D) 
            ALLOCATE(temp1D(na));temp1D=flt0(l,j,:)
            CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
            gelman_flt0(l,j)=gelmanval;DEALLOCATE(temp1D) 
          ENDDO
        ENDIF
      ENDDO
    ENDIF
    IF(ntr.GT.0.AND.testrampe.GT.0)THEN
      DO l=1,ntr
        IF(ramporder(l).GT.0.AND.rampmod.EQ.'exp')THEN
          di=ramporder(l)
          DO j=1,di
            SELECT CASE(j)
              CASE(1)
                ALLOCATE(temp1D(na));temp1D=t1ramp(l,:)
                CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
                gelman_t1ramp(l)=gelmanval;DEALLOCATE(temp1D)       
              CASE(2)
                ALLOCATE(temp1D(na));temp1D=t2ramp(l,:)
                CALL gelmancomp(na,nchain,nburn,nat,temp1D,gelmanval)  
                gelman_t2ramp(l)=gelmanval;DEALLOCATE(temp1D)      
            END SELECT
          ENDDO
        ENDIF
      ENDDO
    ENDIF
  ENDIF

  k=0
  nacib=0.
  nacinb=0.
  OPEN(UNIT=945,FILE='transit.res')
  DO i=1,na
    IF(burned(i).EQ.'y'.AND.accepted(i).EQ.'y')THEN
      test = CEILING(DBLE(i)/DBLE(nat))
      nacib(test)=nacib(test)+1.
    ENDIF
    IF(burned(i).EQ.'n')THEN
      meanposterior = meanposterior+merit(i)
      IF(b3(1,i).LT.(1.+SQRT(dF(1,i)))) prosecond=prosecond+1.
      IF(b2(1,i).LT.(1.+SQRT(dF(1,i)))) protransit=protransit+1.
      IF(b3(1,i).LT.(1.-SQRT(dF(1,i)))) profullsecond=profullsecond+1.
      IF(b2(1,i).LT.(1.-SQRT(dF(1,i)))) profulltransit=profulltransit+1. 
      IF(accepted(i).EQ.'y')THEN
        test = CEILING(DBLE(i)/DBLE(nat))
        nacinb(test)=nacinb(test)+1.
      ENDIF
      k=k+1
    ENDIF
  ENDDO
  CLOSE(945)
  mco = k
  meanposterior=meanposterior/DBLE(mco)
  prosecond=prosecond/DBLE(mco)
  profullsecond=profullsecond/DBLE(mco)
  protransit=protransit/DBLE(mco)
  profulltransit=profulltransit/DBLE(mco)
  k=0
  PRINT*, 'STATISTICS BASED ON ',mco,' POINTS'
  WRITE(444,*) 'STATISTICS BASED ON ',mco,' POINTS'
  WRITE(445,*) 'STATISTICS BASED ON ',mco,' POINTS'

  ALLOCATE(merit_end(mco)); ALLOCATE(remerit_end(mco))
  ALLOCATE(rvmerit_end(mco));ALLOCATE(photmerit_end(mco))
  ALLOCATE(dur_end(npla,mco));ALLOCATE(octime_end(npla,mco))
  ALLOCATE(t0_end(npla,mco));ALLOCATE(per_end(npla,mco))
  ALLOCATE(tidel_end(npla,mco))
  ALLOCATE(dF_end(npla,mco)); ALLOCATE(b_end(npla,mco))
  ALLOCATE(secosw_end(npla,mco)); ALLOCATE(sesinw_end(npla,mco))
  ALLOCATE(exc_end(npla,mco)); ALLOCATE(omega_end(npla,mco))
  ALLOCATE(ecosw_end(npla,mco)); ALLOCATE(esinw_end(npla,mco))
  ALLOCATE(b2_end(npla,mco)); ALLOCATE(b3_end(npla,mco))
  ALLOCATE(rho_end(mco)); ALLOCATE(rhop_end(npla,mco))
  ALLOCATE(rr_end(npla,mco)); ALLOCATE(a_R_end(npla,mco)); ALLOCATE(temp_end(mco))
  ALLOCATE(radius_s_end(mco));ALLOCATE(radius_p_end(npla,mco))
  ALLOCATE(mass_s_end(mco));ALLOCATE(temp_s_end(mco))
  ALLOCATE(met_s_end(mco));ALLOCATE(lum_s_end(mco))
  ALLOCATE(logg_p_end(npla,mco));ALLOCATE(teq_p_end(npla,mco))
  ALLOCATE(hillrad_p_end(npla,mco))
  ALLOCATE(mass_p_end(npla,mco));ALLOCATE(mass_p_sini_end(npla,mco))
  ALLOCATE(prtr_end(npla,mco));ALLOCATE(proc_end(npla,mco))
  ALLOCATE(semi_end(npla,mco));ALLOCATE(inclian_end(npla,mco))
  ALLOCATE(roche_end(npla,mco));ALLOCATE(a_roche_end(npla,mco))
  ALLOCATE(kb_end(npla,mco)); ALLOCATE(ka_end(npla,mco))
  ALLOCATE(dFsec_end(nfi,mco,npla));ALLOCATE(phampli1_end(nfi,mco));
  ALLOCATE(phampli2_end(nfi,mco));ALLOCATE(phampli3_end(nfi,mco));
  ALLOCATE(phoffset_end(nfi,mco));ALLOCATE(logg_end(mco))
  ALLOCATE(safro_end(npla,mco));ALLOCATE(irrad_end(npla,mco))
  allocate(age_end(mco));allocate(col_s_end(mco))
  IF(ngroup.GT.0)THEN
    ALLOCATE(dfgroup_end(ngroup,mco))
  ENDIF
  IF(enoch.EQ.'y')THEN
    ALLOCATE(massjitter_end(mco))
  ENDIF
  IF(testf2.EQ.'y')THEN
    ALLOCATE(f2_end(mco));ALLOCATE(ktide_end(mco))
  ENDIF
  IF(nddf.GT.0.AND.ntr.GT.0.AND.isddf.NE.'n')THEN
    ALLOCATE(ddf_end(npla,nddf,mco))
    ALLOCATE(ddepth_end(npla,nddf,mco))
    ALLOCATE(dratio_end(npla,nddf,mco))
    ALLOCATE(radipla_end(npla,nddf,mco))
  ENDIF
  IF(nttvmax.GT.0.AND.isttv.EQ.'y')THEN
    ALLOCATE(ttv_end(npla,nttvmax,mco));ALLOCATE(ttr_end(npla,nttvmax,mco))
  ENDIF
  IF(ntr.GT.0.AND.testsin.GT.0)THEN
    ALLOCATE(sip_end(ntr,4,mco));ALLOCATE(sit0_end(ntr,4,mco))
  ENDIF
  IF(ntr.GT.0.AND.testflare.GT.0)THEN
    ALLOCATE(fltau_end(ntr,4,mco));ALLOCATE(flampli_end(ntr,4,mco))
    ALLOCATE(flt0_end(ntr,4,mco))
  ENDIF
  IF(ntr.GT.0.AND.testrampe.GT.0)THEN
    ALLOCATE(t1ramp_end(ntr,mco))
  ENDIF
  IF(ntr.GT.0.AND.testrampe.GT.1)THEN
    ALLOCATE(t2ramp_end(ntr,mco))
  ENDIF
  IF(limb.EQ.'qd')THEN
    ALLOCATE(jumplimb_end(nfi,nd,mco))
    ALLOCATE(ql_end(nfi,nd,mco))
  ENDIF
  ALLOCATE(svsinicosbeta_end(mco)); ALLOCATE(svsinisinbeta_end(mco))
  ALLOCATE(vsini_end(mco)); ALLOCATE(beta_end(mco))
  IF(nrv.GT.0.AND.iftrendrv.EQ.'y') ALLOCATE(rvtrendco_end(mco,nrvparglo))
  IF(stelincli.EQ.'y')THEN
    ALLOCATE(incli_s_end(mco));ALLOCATE(sinincli_s_end(mco))
  ENDIF
  k=0
  DO i=1,na
    IF(burned(i).EQ.'n')THEN
      k=k+1
      remerit_end(k)=remerit(i)
      merit_end(k)=merit(i)
      rvmerit_end(k)=rvmerit(i)
      photmerit_end(k)=photmerit(i)
      radius_s_end(k)=radius_s(i)
      logg_end(k)=logg(i)
      if (isoch.eq.'y') age_end(k)=age(i)
      mass_s_end(k)=mass_s(i)
      met_s_end(k)=met_s(i)
      temp_s_end(k)=temp_s(i)
      if (priorcol.eq.'p') col_s_end(k)=col_s(i)
      lum_s_end(k)=lum_s(i)
      rho_end(k)=rho(i)
      IF(stelincli.EQ.'y')THEN
        incli_s_end(k)=incli_s(i)
        sinincli_s_end(k)=sinincli_s(i)
      ENDIF
      IF(ngroup.GT.0)THEN
        DO j=1,ngroup
          dfgroup_end(j,k)=dfgroup(j,i)
        ENDDO
      ENDIF
      DO j=1,npla
        IF(nddf.GT.0.AND.ntr.GT.0.AND.isddf.NE.'n')THEN
          DO l=1,nddf
            ddf_end(j,l,k)=ddf(j,l,i)
            dratio_end(j,l,k)=dratio(j,l,i)
            ddepth_end(j,l,k)=ddepth(j,l,i)
            radipla_end(j,l,k)=radipla(j,l,i)
          ENDDO
        ENDIF
        radius_p_end(j,k)=radius_p(j,i)
        mass_p_end(j,k)=mass_p(j,i)
        mass_p_sini_end(j,k)=mass_p_sini(j,i)
        semi_end(j,k)=semi(j,i)
        irrad_end(j,k)=irrad(j,i)
        roche_end(j,k)=roche(j,i)
        a_roche_end(j,k)=a_roche(j,i)
        prtr_end(j,k)=prtr(j,i)
        proc_end(j,k)=proc(j,i)
        inclian_end(j,k)=inclian(j,i)
        b_end(j,k)=b(j,i)
        octime_end(j,k)=octime(j,i)
        dF_end(j,k)=dF(j,i)
        per_end(j,k)=per(j,i)
        tidel_end(j,k)=tidel(j,i)
        dur_end(j,k)=dur(j,i)
        t0_end(j,k)=t0(j,i)
        secosw_end(j,k)=secosw(j,i)
        sesinw_end(j,k)=sesinw(j,i)
        ecosw_end(j,k)=ecosw(j,i)
        esinw_end(j,k)=esinw(j,i)
        exc_end(j,k)=exc(j,i)
        omega_end(j,k)=omega(j,i)
        IF(omega_end(j,k).LT.omega(j,soluce)-180.) omega_end(j,k)=omega_end(j,k)+360.
        IF(omega_end(j,k).GT.omega(j,soluce)+180.) omega_end(j,k)=omega_end(j,k)-360.
        b2_end(j,k)=b2(j,i)
        b3_end(j,k)=b3(j,i)
        rhop_end(j,k)=rhop(j,i)
        a_R_end(j,k)=a_R(j,i)
        rr_end(j,k)=rr(j,i)
        logg_p_end(j,k)=logg_p(j,i)
        teq_p_end(j,k)=teq_p(j,i)
        hillrad_p_end(j,k)=hillrad_p(j,i)
        safro_end(j,k)=safro(j,i)
      ENDDO
      IF(enoch.EQ.'y')THEN
        massjitter_end(k)=massjitter(i)
      ENDIF
      IF(testf2.EQ.'y')THEN
        f2_end(k)=f2(i)
        ktide_end(k)=ktide(i)
      ENDIF
      IF(ntr.GT.0)THEN
        DO j=1,nfi 
          DO l=1,npla
            dFsec_end(j,k,l)=dFsec(j,i,l)
          ENDDO
          phampli1_end(j,k)=phampli1(j,i)
          phampli2_end(j,k)=phampli2(j,i)
          phampli3_end(j,k)=phampli3(j,i)
          phoffset_end(j,k)=phoffset(j,i)
        ENDDO
      ENDIF
      IF(nttvmax.GT.0.AND.isttv.EQ.'y')THEN
        DO j=1,npla
          DO l=1,nttv(j)
            ttv_end(j,l,k)=ttv(j,l,i)
            ttr_end(j,l,k)=ttr(j,l,i)
          ENDDO
        ENDDO
      ENDIF
      IF(testsin.GT.0.AND.ntr.GT.0)THEN
        DO j=1,ntr
          di = sinusnumber(j)
          DO l=1,di
            sit0_end(j,l,k)=sit0(j,l,i)
            sip_end(j,l,k)=sip(j,l,i)
          ENDDO
        ENDDO
      ENDIF
      IF(testflare.GT.0.AND.ntr.GT.0)THEN
        DO j=1,ntr
          di = flarenumber(j)
          DO l=1,di
            flt0_end(j,l,k)=flt0(j,l,i)
            flampli_end(j,l,k)=flampli(j,l,i)
            fltau_end(j,l,k)=fltau(j,l,i)
          ENDDO
        ENDDO
      ENDIF
      IF(testrampe.GT.0.AND.ntr.GT.0)THEN
        DO j=1,ntr
          di = ramporder(j)
          DO l=1,di
            IF(l.EQ.1) t1ramp_end(j,k)=t1ramp(j,i)
            IF(l.EQ.2) t2ramp_end(j,k)=t2ramp(j,i)
          ENDDO
        ENDDO
      ENDIF
      IF(limb.EQ.'qd')THEN
        DO j=1,nfi
          jumplimb_end(j,1,k)=jumplimb(j,1,i)
          jumplimb_end(j,2,k)=jumplimb(j,2,i)
          ql_end(j,1,k)=ql(j,1,i)
          ql_end(j,2,k)=ql(j,2,i)
        ENDDO
      ENDIF
      DO j=1,npla
        kb_end(j,k)=kb(j,i)
        ka_end(j,k)=ka(j,i)
      ENDDO
      svsinicosbeta_end(k)=svsinicosbeta(i)
      svsinisinbeta_end(k)=svsinisinbeta(i)
      vsini_end(k)=vsini(i)
      beta_end(k)=beta(i)
      IF(beta_end(k).LT.beta(soluce)-180.) beta_end(k)=beta_end(k)+360.
      IF(beta_end(k).GT.beta(soluce)+180.) beta_end(k)=beta_end(k)-360. 
      IF(nrv.GT.0.AND.iftrendrv.EQ.'y')THEN
        DO j=1,nrvparglo
          rvtrendco_end(k,j)=rvtrendco(i,j)
        ENDDO
      ENDIF
    ENDIF
  ENDDO

  medi = INT(mco/2)
  limmed(5) = medi
  limmed(1) = CEILING((mco*(1.-0.683))/2.)
  limmed(2) = limmed(1) + FLOOR(mco*0.683)
  limmed(3) = CEILING((mco*(1.-0.997))/2.)
  limmed(4) = limmed(3) + FLOOR(mco*0.997)
  be_int1sig = CEILING((mco*(1.-0.683))/2.)
  en_int1sig = be_int1sig + FLOOR(mco*0.683)
  be_int3sig = CEILING((mco*(1.-0.997))/2.)
  en_int3sig = be_int3sig + FLOOR(mco*0.997)
  PRINT*, 'Median=',medi
  PRINT*, '1sig_b=',be_int1sig
  PRINT*, '1sig_t=',en_int1sig
  PRINT*, '3sig_b=',be_int3sig
  PRINT*, '3sig_t=',en_int3sig
  WRITE(444,*) 'Median=',medi
  WRITE(444,*) '1sig_b=',be_int1sig
  WRITE(444,*) '1sig_t=',en_int1sig
  WRITE(444,*) '3sig_b=',be_int3sig
  WRITE(444,*) '3sig_t=',en_int3sig
  WRITE(445,*) 'Median=',medi
  WRITE(445,*) '1sig_b=',be_int1sig
  WRITE(445,*) '1sig_t=',en_int1sig
  WRITE(445,*) '3sig_b=',be_int3sig
  WRITE(445,*) '3sig_t=',en_int3sig

  DO i=1,nchain
    nacib(i)=nacib(i)/nburn
    nacinb(i)=nacinb(i)/(nat-nburn)
  ENDDO
  
  PRINT*,      '------------------------------------------------------------------------------------------'
  WRITE(444,*) '------------------------------------------------------------------------------------------'
  WRITE(445,*) '------------------------------------------------------------------------------------------'
  PRINT*,      'Acceptance rate for the burn-in phase and the rest'
  WRITE(444,*) 'Acceptance rate for the burn-in phase and the rest'
  WRITE(445,*) 'Acceptance rate for the burn-in phase and the rest'
  DO i=1,nchain
    WRITE(*,132)   ' Chain ',i,':',nacib(i),nacinb(i)
    WRITE(444,132) ' Chain ',i,':',nacib(i),nacinb(i)
    WRITE(445,132) ' Chain ',i,':',nacib(i),nacinb(i)
  ENDDO
  PRINT*,      '------------------------------------------------------------------------------------------'
  WRITE(444,*) '------------------------------------------------------------------------------------------'
  WRITE(445,*) '------------------------------------------------------------------------------------------'

  PRINT*, 'RV global SVD parameters: median of the MPDF + 68.3 and 99.7% limits + Unit' 
  WRITE(444,*) 'RV global SVD parameters: median of the MPDF + 68.3 and 99.7% limits + Unit' 
  WRITE(445,*) 'RV global SVD parameters: median of the MPDF + 68.3 and 99.7% limits + Unit' 

  ALLOCATE(temp1D(mco))
  ylabel = ' Number of MCMC steps'; nbinh = 20  
  jump='n'

  IF(nrv.GT.0)THEN
    IF(iftrendrv.EQ.'y')THEN
      DO i=1,nrvparglo
         test = 0
         SELECT CASE(i)
           CASE(1)            
              texb =' RV g-zero =       ';texunit='   ms-1';bf_val=rvtrendco(soluce,1);& 
               & temp1D=rvtrendco_end(:,1);test=1
              fil2 = 'rv_gzero.res'; graphfile = 'rv_gzero_h.sm'
              xlabel = '     RV-0 [m/s]     '
           CASE(2)         
              texb =' RV g-slope =      ';texunit=' ms-1/d';bf_val=rvtrendco(soluce,2);& 
              & temp1D=rvtrendco_end(:,2);test=1
              fil2 = 'rv_gslope.res'; graphfile = 'rv_gslope_h.sm'
              xlabel = '    RV-1 [m/s/d]    '
           CASE(3)
             IF(ordertrendrv.EQ.2)THEN        
               texb =' RV g-curvature = ';texunit='ms-1/d2';bf_val=rvtrendco(soluce,3);& 
               & temp1D=rvtrendco_end(:,3);test=1
              fil2 = 'rv_gcurve.res'; graphfile = 'rv_gcurve_h.sm'
              xlabel = '   RV-2 [m/s/d^2]   '
             ENDIF
         END SELECT
         IF(test.EQ.1)THEN
           OPEN(UNIT=124,FILE=fil2)
           DO l=1,mco
             WRITE(124,*) l,temp1D(l)
           ENDDO
           CLOSE(124)
           CALL SORT(mco,temp1D)
           CALL PDLIM(mco,temp1D,bf_val,limmed,limpd)
           medi_val = temp1d(limmed(5))
           CALL HISTO(fil2,graphfile,bf_val,medi_val,xlabel,nbinh)
           CALL PDFSCREEN(texb,texunit,bf_val,medi_val,limpd,gelman,gelmanval,nchain,jump)
         ENDIF
      ENDDO
    ENDIF
  ENDIF
          
  PRINT*,      '------------------------------------------------------------------------------------------'
  WRITE(444,*) '------------------------------------------------------------------------------------------'
  WRITE(445,*) '------------------------------------------------------------------------------------------'

  PRINT*, 'MCMC jump parameters: median of the MPDF + 68.3 and 99.7% limits + Unit (+ GR)' 
  WRITE(444,*) 'MCMC jump parameters: median of the MPDF + 68.3 and 99.7% limits + Unit (+ GR)' 
  WRITE(445,*) 'MCMC jump parameters: median of the MPDF + 68.3 and 99.7% limits + Unit (+ GR)' 

  jump='y'

  WRITE(*,'(A4)') 'LCs'
  WRITE(444,'(A4)') 'LCs'  
  WRITE(445,'(A4)') 'LCs'
  IF(ntr.GT.0)THEN
    DO i=1,ntr
      IF(i.LT.10)THEN
        let1=CHAR(i+48)
      ELSE IF(i.LT.100)THEN
        let1=CHAR(48+INT(i/10)); let2=CHAR(48+MOD(i,10))
      ELSE
        let1 = CHAR(48+INT(i/100))
        i2 = i - INT(i/100)*100
        let2 = CHAR(48+INT(i2/10))
        let3 = CHAR(48+MOD(i2,10))
      ENDIF
      IF(testsin.GT.0.AND.sinusnumber(i).GT.0)THEN              ! Sinusoids t0 and P 
        di=sinusnumber(i)
        DO j=1,di
          DO k=1,2
            SELECT CASE (k)
              CASE (1)
                IF(i.LT.10)THEN 
                  texb =' SinT0-'//let1//'_'//CHAR(j+48)//' =       '
                  fil2 = 'SinT0_'//let1//'.res'
                  graphfile = 'SinT0_'//let1//'_h.sm'
                  graphfile2 = 'SinT0_'//let1//'.sm'
                ELSE IF(i.LT.100)THEN
                  texb =' SinT0-'//let1//let2//'-'//CHAR(j+48)//' =      '
                  fil2 = 'SinT0_'//let1//let2//'.res' 
                  graphfile = 'SinT0_'//let1//let2//'_h.sm'
                  graphfile2 = 'SinT0_'//let1//let2//'.sm'
                ELSE
                  texb =' SinT0-'//let1//let2//let3//'-'//CHAR(j+48)//' =     '
                  fil2 = 'SinT0_'//let1//let2//let3//'.res' 
                  graphfile = 'SinT0_'//let1//let2//let3//'_h.sm'
                  graphfile2 = 'SinT0_'//let1//let2//let3//'.sm'
                ENDIF
                xlabel = '      SinT0 [JD]    '
                texunit = '     JD';bf_val=sit0(i,j,soluce);temp1D=sit0_end(i,j,:)
                IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_sit0(i,j)
                IF(ifacf.EQ.'y'.AND.i.LT.10) THEN
                  acffile= 'sinT0_'//let1 //'.acf'
                  graphacffile='sinT0_'//let1//'_acf.sm'
                ELSE IF(ifacf.EQ.'y'.AND.i.LT.100) THEN
                  acffile= 'sinT0_'//let1//let2//'.acf'
                  graphacffile='sinT0_'//let1//let2//'_acf.sm'
                ELSE IF(ifacf.EQ.'y')THEN
                  acffile= 'sinT0_'//let1//let2//let3//'.acf'
                  graphacffile='sinT0_'//let1//let2//let3//'_acf.sm'
                ENDIF
              CASE (2)
                IF(i.LT.10)THEN 
                  texb =' SinP-'//let1//'-'//CHAR(j+48)//' =        '
                  fil2 = 'SinP_'//let1//'.res'
                  graphfile = 'SinP_'//let1//'_h.sm'
                  graphfile2 = 'SinP_'//let1//'.sm'
                ELSE IF(i.LT.100)THEN
                  texb =' SinP-'//let1//let2//'-'//CHAR(j+48)//' =       '
                  fil2 = 'SinP_'//let1//let2//'.res'
                  graphfile = 'SinP_'//let1//let2//'_h.sm'
                  graphfile2 = 'SinP_'//let1//let2//'.sm' 
                ELSE
                  texb =' SinP-'//let1//let2//let3//'-'//CHAR(j+48)//' =      '
                  fil2 = 'SinP_'//let1//let2//let3//'.res'
                  graphfile = 'SinP_'//let1//let2//let3//'_h.sm'
                  graphfile2 = 'SinP_'//let1//let2//let3//'.sm' 
                ENDIF
                xlabel = '     SinP [days]    '
                texunit = '   days';bf_val=sip(i,j,soluce);temp1D=sip_end(i,j,:)
                IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_sip(i,j)
                IF(ifacf.EQ.'y'.AND.i.LT.10) THEN
                  acffile= 'sinP_'//let1//'.acf'
                  graphacffile='sinP_'//let1//'_acf.sm'
                ELSE IF(ifacf.EQ.'y'.AND.i.LT.100) THEN
                  acffile= 'sinP_'//let1//let2//'.acf'
                  graphacffile='sinP_'//let1//let2//'_acf.sm'
                ELSE IF(ifacf.EQ.'y') THEN
                  acffile= 'sinP_'//let1//let2//let3//'.acf'
                  graphacffile='sinP_'//let1//let2//let3//'_acf.sm'
                ENDIF
            END SELECT
            OPEN(UNIT=124,FILE=fil2)
            DO l=1,mco
              WRITE(124,*) l,temp1D(l)
            ENDDO
            CLOSE(124)
           IF(ifacf.EQ.'y')THEN
             nu = (mco/(nchain*2))
             ALLOCATE(acf(nu))
             OPEN(UNIT=135,FILE=acffile)
             DO l=1,nu
               CALL AUTOCOR(mco,temp1D,nchain,l-1,acf(l))
               WRITE(135,*) l-1,acf(l)
             ENDDO
             CLOSE(135)
             min_x = 0-nu*0.05;max_x=nu*1.05
             min_y = MINVAL(acf)-(MAXVAL(acf)-MINVAL(acf))*0.05
             max_y = MAXVAL(acf)+(MAXVAL(acf)-MINVAL(acf))*0.05
             ncol_file=2;abs_val=1;ord_val=2
             xlabel2='    Lag (Nsteps)    ';ylabel='   Autocorrelation  '
             abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
             errory='n';errorx='n';cone='n';newf='y'
             CALL GRAPH2D(acffile,graphacffile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
               & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
               & ord_errorp,1,cone,newf,0)
             DEALLOCATE(acf)
           ENDIF
           min_x = 0-mco*0.05;max_x=mco*1.05
           min_y = MINVAL(temp1D)-ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
           max_y = MAXVAL(temp1D)+ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
           ncol_file=2;abs_val=1;ord_val=2
           xlabel2='       Steps        ';ylabel=xlabel
           abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
           errory='n';errorx='n';cone='n';newf='y'
           CALL GRAPH2D(fil2,graphfile2,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
            & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
            & ord_errorp,1,cone,newf,0)
           CALL SORT(mco,temp1D)
           CALL PDLIM(mco,temp1D,bf_val,limmed,limpd)
           medi_val = temp1d(limmed(5))
           CALL HISTO(fil2,graphfile,bf_val,medi_val,xlabel,nbinh)
           CALL PDFSCREEN(texb,texunit,bf_val,medi_val,limpd,gelman,gelmanval, & 
            & nchain,jump)
         ENDDO
       ENDDO
     ENDIF
      IF(testflare.GT.0.AND.flarenumber(i).GT.0)THEN              ! Flares Ampli and Tau
        di=flarenumber(i)
        DO j=1,di
          DO k=1,3
            SELECT CASE (k)
              CASE (1)
                IF(i.LT.10)THEN 
                  texb =' Flat0-'//let1//'_'//CHAR(j+48)//' =       '
                  fil2 = 'Flat0_'//let1//'.res'
                  graphfile = 'Flat0_'//let1//'_h.sm'
                  graphfile2 = 'Flat0_'//let1//'.sm'
                ELSE IF(i.LT.100)THEN
                  texb =' Flat0-'//let1//let2//'-'//CHAR(j+48)//' =      '
                  fil2 = 'Flat0_'//let1//let2//'.res' 
                  graphfile = 'Flat0_'//let1//let2//'_h.sm'
                  graphfile2 = 'Flat0_'//let1//let2//'.sm'
                ELSE
                  texb =' Flat0-'//let1//let2//let3//'-'//CHAR(j+48)//' =     '
                  fil2 = 'Flat0_'//let1//let2//let3//'.res' 
                  graphfile = 'Flat0_'//let1//let2//let3//'_h.sm'
                  graphfile2 = 'Flat0_'//let1//let2//let3//'.sm'
                ENDIF
                xlabel = '     Flare t0       '
                texunit = '  BJD  ';bf_val=flt0(i,j,soluce);temp1D=flt0_end(i,j,:)
                IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_flt0(i,j)
                IF(ifacf.EQ.'y'.AND.i.LT.10) THEN
                  acffile= 'Flat0_'//let1 //'.acf'
                  graphacffile='Flat0_'//let1//'_acf.sm'
                ELSE IF(ifacf.EQ.'y'.AND.i.LT.100) THEN
                  acffile= 'Flat0_'//let1//let2//'.acf'
                  graphacffile='Flat0_'//let1//let2//'_acf.sm'
                ELSE IF(ifacf.EQ.'y')THEN
                  acffile= 'Flat0_'//let1//let2//let3//'.acf'
                  graphacffile='Flat0_'//let1//let2//let3//'_acf.sm'
                ENDIF
              CASE (2)
                IF(i.LT.10)THEN 
                  texb =' FlaAm-'//let1//'_'//CHAR(j+48)//' =       '
                  fil2 = 'FlaAm_'//let1//'.res'
                  graphfile = 'FlaAm_'//let1//'_h.sm'
                  graphfile2 = 'FlaAm_'//let1//'.sm'
                ELSE IF(i.LT.100)THEN
                  texb =' FlaAm-'//let1//let2//'-'//CHAR(j+48)//' =      '
                  fil2 = 'FlaAm_'//let1//let2//'.res' 
                  graphfile = 'FlaAm_'//let1//let2//'_h.sm'
                  graphfile2 = 'FlaAm_'//let1//let2//'.sm'
                ELSE
                  texb =' FlaAm-'//let1//let2//let3//'-'//CHAR(j+48)//' =     '
                  fil2 = 'FlaAm_'//let1//let2//let3//'.res' 
                  graphfile = 'FlaAm_'//let1//let2//let3//'_h.sm'
                  graphfile2 = 'FlaAm_'//let1//let2//let3//'.sm'
                ENDIF
                xlabel = '     Flare Amp      '
                texunit = '       ';bf_val=flampli(i,j,soluce);temp1D=flampli_end(i,j,:)
                IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_flampli(i,j)
                IF(ifacf.EQ.'y'.AND.i.LT.10) THEN
                  acffile= 'FlaAm_'//let1 //'.acf'
                  graphacffile='FlaAm_'//let1//'_acf.sm'
                ELSE IF(ifacf.EQ.'y'.AND.i.LT.100) THEN
                  acffile= 'FlaAm_'//let1//let2//'.acf'
                  graphacffile='FlaAm_'//let1//let2//'_acf.sm'
                ELSE IF(ifacf.EQ.'y')THEN
                  acffile= 'FlaAm_'//let1//let2//let3//'.acf'
                  graphacffile='FlaAm_'//let1//let2//let3//'_acf.sm'
                ENDIF
              CASE (3)
                IF(i.LT.10)THEN 
                  texb =' FlaT-'//let1//'-'//CHAR(j+48)//' =        '
                  fil2 = 'FlaT_'//let1//'.res'
                  graphfile = 'FlaT_'//let1//'_h.sm'
                  graphfile2 = 'FlaT_'//let1//'.sm'
                ELSE IF(i.LT.100)THEN
                  texb =' FlaT-'//let1//let2//'-'//CHAR(j+48)//' =       '
                  fil2 = 'Flat_'//let1//let2//'.res'
                  graphfile = 'FlaT_'//let1//let2//'_h.sm'
                  graphfile2 = 'FlaT_'//let1//let2//'.sm' 
                ELSE
                  texb =' FlaT-'//let1//let2//let3//'-'//CHAR(j+48)//' =      '
                  fil2 = 'FlaT_'//let1//let2//let3//'.res'
                  graphfile = 'FlaT_'//let1//let2//let3//'_h.sm'
                  graphfile2 = 'FlaT_'//let1//let2//let3//'.sm' 
                ENDIF
                xlabel = '     FlaT [days]    '
                texunit = '   days';bf_val=fltau(i,j,soluce);temp1D=fltau_end(i,j,:)
                IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_fltau(i,j)
                IF(ifacf.EQ.'y'.AND.i.LT.10) THEN
                  acffile= 'FlaT_'//let1//'.acf'
                  graphacffile='FlaT_'//let1//'_acf.sm'
                ELSE IF(ifacf.EQ.'y'.AND.i.LT.100) THEN
                  acffile= 'FlaT_'//let1//let2//'.acf'
                  graphacffile='FlaT_'//let1//let2//'_acf.sm'
                ELSE IF(ifacf.EQ.'y') THEN
                  acffile= 'FlaT_'//let1//let2//let3//'.acf'
                  graphacffile='FlaT_'//let1//let2//let3//'_acf.sm'
                ENDIF
            END SELECT
            OPEN(UNIT=124,FILE=fil2)
            DO l=1,mco
              WRITE(124,*) l,temp1D(l)
            ENDDO
            CLOSE(124)
           IF(ifacf.EQ.'y')THEN
             nu = (mco/(nchain*2))
             ALLOCATE(acf(nu))
             OPEN(UNIT=135,FILE=acffile)
             DO l=1,nu
               CALL AUTOCOR(mco,temp1D,nchain,l-1,acf(l))
               WRITE(135,*) l-1,acf(l)
             ENDDO
             CLOSE(135)
             min_x = 0-nu*0.05;max_x=nu*1.05
             min_y = MINVAL(acf)-(MAXVAL(acf)-MINVAL(acf))*0.05
             max_y = MAXVAL(acf)+(MAXVAL(acf)-MINVAL(acf))*0.05
             ncol_file=2;abs_val=1;ord_val=2
             xlabel2='    Lag (Nsteps)    ';ylabel='   Autocorrelation  '
             abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
             errory='n';errorx='n';cone='n';newf='y'
             CALL GRAPH2D(acffile,graphacffile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
               & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
               & ord_errorp,1,cone,newf,0)
             DEALLOCATE(acf)
           ENDIF
           min_x = 0-mco*0.05;max_x=mco*1.05
           min_y = MINVAL(temp1D)-ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
           max_y = MAXVAL(temp1D)+ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
           ncol_file=2;abs_val=1;ord_val=2
           xlabel2='       Steps        ';ylabel=xlabel
           abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
           errory='n';errorx='n';cone='n';newf='y'
           CALL GRAPH2D(fil2,graphfile2,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
            & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
            & ord_errorp,1,cone,newf,0)
           CALL SORT(mco,temp1D)
           CALL PDLIM(mco,temp1D,bf_val,limmed,limpd)
           medi_val = temp1d(limmed(5))
           CALL HISTO(fil2,graphfile,bf_val,medi_val,xlabel,nbinh)
           CALL PDFSCREEN(texb,texunit,bf_val,medi_val,limpd,gelman,gelmanval, & 
            & nchain,jump)
         ENDDO
       ENDDO
     ENDIF
     IF(ramporder(i).GT.0.AND.rampmod.EQ.'exp')THEN  ! Ramp double exponential T1 and T2
        DO k=1,2
          SELECT CASE (k)
            CASE (1)
              IF(i.LT.10)THEN 
                texb =' T1-'//let1//' =             '
                fil2 = 'Er_T1_'//let1//'.res'
                graphfile = 'Er_T1_'//let1//'_h.sm'
                graphfile2 = 'Er_T1_'//let1//'.sm'
              ELSE IF(i.LT.100)THEN
                texb =' T1-'//let1//let2//' =            '
                fil2 = 'Er_T1_'//let1//let2//'.res'
                graphfile = 'Er_T1_'//let1//let2//'_h.sm'
                graphfile2 = 'Er_T1_'//let1//let2//'.sm'
              ELSE
                texb =' T1-'//let1//let2//let3//' =           '
                fil2 = 'Er_T1_'//let1//let2//let3//'.res'
                graphfile = 'Er_T1_'//let1//let2//let3//'_h.sm'
                graphfile2 = 'Er_T1_'//let1//let2//let3//'.sm'
              ENDIF  
              xlabel = '   ErampT1 [days]   '
              texunit = '   days';bf_val=t1ramp(i,soluce);temp1D=t1ramp_end(i,:)
              IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_t1ramp(i)
              IF(ifacf.EQ.'y'.AND.i.LT.10) THEN
                acffile= 'ErampT1_'//let1 //'.acf'
                graphacffile='ErampT1_'//let1//'_acf.sm'
              ELSE IF(ifacf.EQ.'y'.AND.i.LT.100) THEN
                acffile= 'ErampT1_'//let1 //let2//'.acf'
                graphacffile='ErampT1_'//let1//let2//'_acf.sm'
              ELSE IF(ifacf.EQ.'y') THEN
                acffile= 'ErampT1_'//let1//let2//let3//'.acf'
                graphacffile='ErampT1_'//let1//let2//let3//'_acf.sm'
              ENDIF
            CASE (2)
              IF(i.LT.10)THEN 
                texb =' T2-'//let1//' =             '
                fil2 = 'Er_T2_'//let1//'.res'
                graphfile = 'Er_T2_'//let1//'_h.sm'
                graphfile2 = 'Er_T2_'//let1//'.sm'
              ELSE IF(i.LT.100)THEN
                texb =' T2-'//let1//let2//' =            '
                fil2 = 'Er_T2_'//let1//let2//'.res'
                graphfile = 'Er_T2_'//let1//let2//'_h.sm'
                graphfile2 = 'Er_T2_'//let1//let2//'.sm'
              ELSE
                texb =' T2-'//let1//let2//let3//' =           '
                fil2 = 'Er_T2_'//let1//let2//let3//'.res'
                graphfile = 'Er_T2_'//let1//let2//let3//'_h.sm'
                graphfile2 = 'Er_T2_'//let1//let2//let3//'.sm'
              ENDIF
              xlabel = '   ErampT2 [days]   '
              texunit = '   days';bf_val=t2ramp(i,soluce);temp1D=t2ramp_end(i,:)
              IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_t2ramp(i)
              IF(ifacf.EQ.'y'.AND.i.LT.10) THEN
                acffile= 'Er_T2_'//let1//'.acf'
                graphacffile='Er_T2_'//let1 //'_acf.sm'
              ELSE IF(ifacf.EQ.'y'.AND.i.LT.100) THEN
                acffile= 'Er_T2_'//let1 //let2//'.acf'
                graphacffile='Er_T2_'//let1//let2//'_acf.sm'
              ELSE IF(ifacf.EQ.'y') THEN
                acffile= 'Er_T2_'//let1 //let2//let3//'.acf'
                graphacffile='Er_T2_'//let1//let2//let3//'_acf.sm'
              ENDIF
          END SELECT
          OPEN(UNIT=124,FILE=fil2)
          DO l=1,mco
            WRITE(124,*) l,temp1D(l)
          ENDDO
          CLOSE(124)
          IF(ifacf.EQ.'y')THEN
            nu = (mco/(nchain*2))
            ALLOCATE(acf(nu))
            OPEN(UNIT=135,FILE=acffile)
            DO l=1,nu
              CALL AUTOCOR(mco,temp1D,nchain,l-1,acf(l))
              WRITE(135,*) l-1,acf(l)
            ENDDO
            CLOSE(135)
            min_x = 0-nu*0.05;max_x=nu*1.05
            min_y = MINVAL(acf)-(MAXVAL(acf)-MINVAL(acf))*0.05
            max_y = MAXVAL(acf)+(MAXVAL(acf)-MINVAL(acf))*0.05
            ncol_file=2;abs_val=1;ord_val=2
            xlabel2='    Lag (Nsteps)    ';ylabel='   Autocorrelation  '
            abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
            errory='n';errorx='n';cone='n';newf='y'
            CALL GRAPH2D(acffile,graphacffile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
              & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
              & ord_errorp,1,cone,newf,0)
            DEALLOCATE(acf)
          ENDIF
          min_x = 0-mco*0.05;max_x=mco*1.05
          min_y = MINVAL(temp1D)-ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
          max_y = MAXVAL(temp1D)+ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
          ncol_file=2;abs_val=1;ord_val=2
          xlabel2='       Steps        ';ylabel=xlabel
          abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
          errory='n';errorx='n';cone='n';newf='y'
          CALL GRAPH2D(fil2,graphfile2,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,1,cone,newf,0)
          CALL SORT(mco,temp1D)
          CALL PDLIM(mco,temp1D,bf_val,limmed,limpd)
          medi_val = temp1d(limmed(5))
          CALL HISTO(fil2,graphfile,bf_val,medi_val,xlabel,nbinh)
          CALL PDFSCREEN(texb,texunit,bf_val,medi_val,limpd,gelman,gelmanval, & 
           nchain,jump)
        ENDDO
      ENDIF
    ENDDO
  ENDIF

  WRITE(*,'(A5)') 'STAR'
  WRITE(444,'(A5)') 'STAR'  
  WRITE(445,'(A5)') 'STAR'

  DO i=1,6
    test=0
    SELECT CASE(i)
      CASE(1) ! stellar mass
        IF(fixstellar.NE.'y'.AND.massfromr.NE.'y'.AND.enoch.NE.'y'.AND.msrs.NE.'y'.and.isoch.eq.'y')THEN !ne
          fil2 = 'ms.res'; graphfile = 'ms_h.sm'
          graphfile2 = 'ms.sm'
          xlabel = '     M* [Msun]      '
          texb =' M* =               ';texunit='  M_sun'  
          bf_val=mass_s(soluce);temp1D=mass_s_end;test=1
          IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_ms
          IF(ifacf.EQ.'y') THEN
            acffile= 'ms.acf'
            graphacffile='ms_acf.sm'
          ENDIF
        ENDIF
      CASE(2) ! stellar radius
        IF(fixstellar.NE.'y'.AND.(fitmsrs.EQ.'y'.OR.massfromr.EQ.'y'.or.RjumpIso.eq.'y'))THEN !()
          fil2 = 'rs.res'; graphfile = 'rs_h.sm'
          graphfile2 = 'rs.sm'
          xlabel = '     R* [Rsun]      '
          texb =' R* =               ';texunit='  R_sun'  
          bf_val=radius_s(soluce);temp1D=radius_s_end;test=1
          IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_rs
          IF(ifacf.EQ.'y') THEN
            acffile= 'rs.acf'
            graphacffile='rs_acf.sm'
          ENDIF
        ENDIF
      CASE(3) ! Fe/H
        IF(fixstellar.NE.'y')THEN
          fil2 = 'feh.res'; graphfile = 'feh_h.sm'
          graphfile2 = 'feh.sm'
          xlabel = '       [FeH]        '
          texb =' [Fe/H] =           ';texunit='       '
          bf_val=met_s(soluce); temp1D=met_s_end(:); test=1
          IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_met
          IF(ifacf.EQ.'y') THEN
            acffile= 'feh.acf'
            graphacffile='feh_acf.sm'
          ENDIF
        ENDIF
      CASE(4) ! Teff 
        IF(fixstellar.NE.'y')THEN
          fil2 = 'teff.res'; graphfile = 'teff_h.sm'
          graphfile2 = 'teff.sm'
          xlabel = '        Teff [K]        '
          texb =' Teff =             ';texunit='      K'
          bf_val=temp_s(soluce); temp1D=temp_s_end(:); test=1
          IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_teff
          IF(ifacf.EQ.'y') THEN
            acffile= 'teff.acf'
            graphacffile='teff_acf.sm'
          ENDIF
        ENDIF
      case(5) !color
        if (priorcol.eq.'p') then
          fil2 = 'col.res'; graphfile = 'col_h.sm'
          graphfile2 = 'col.sm'
          xlabel = '      color [mag]       '
          texb =' color =            ';texunit='    mag'
          bf_val=col_s(soluce); temp1D=col_s_end(:); test=1
          IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_col
          IF(ifacf.EQ.'y') THEN
            acffile= 'col.acf'
            graphacffile='col_acf.sm'
          ENDIF
        end if
      CASE(6) ! F2 tidal factor
        IF(testf2.EQ.'y')THEN
          fil2 = 'F2.res'; graphfile = 'F2_h.sm'
          graphfile2 = 'F2.sm'
          xlabel = '         F2        '
          texb =' F2 =               ';texunit='       '
          bf_val=f2(soluce);temp1D=f2_end;test=1
          IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_f2
          IF(ifacf.EQ.'y') THEN
            acffile= 'F2.acf'
            graphacffile='F2_acf.sm'
          ENDIF
        ENDIF
    END SELECT
    IF(test.EQ.1)THEN
      OPEN(UNIT=124,FILE=fil2)
      DO l=1,mco  
        WRITE(124,*) l,temp1D(l)
      ENDDO
      CLOSE(124)
      IF(ifacf.EQ.'y')THEN
        nu = (mco/(nchain*2))
        ALLOCATE(acf(nu))
        OPEN(UNIT=135,FILE=acffile)
        DO l=1,nu
          CALL AUTOCOR(mco,temp1D,nchain,l-1,acf(l))
          WRITE(135,*) l-1,acf(l)
        ENDDO
        CLOSE(135)
        min_x = 0-nu*0.05;max_x=nu*1.05
        min_y = MINVAL(acf)-(MAXVAL(acf)-MINVAL(acf))*0.05
        max_y = MAXVAL(acf)+(MAXVAL(acf)-MINVAL(acf))*0.05
        ncol_file=2;abs_val=1;ord_val=2
        xlabel2='    Lag (Nsteps)    ';ylabel='   Autocorrelation  '
        abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
        errory='n';errorx='n';cone='n';newf='y'
        CALL GRAPH2D(acffile,graphacffile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
          & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
          & ord_errorp,1,cone,newf,0)
        DEALLOCATE(acf)
      ENDIF
      min_x = 0-mco*0.05;max_x=mco*1.05
      min_y = MINVAL(temp1D)-ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
      max_y = MAXVAL(temp1D)+ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
      ncol_file=2;abs_val=1;ord_val=2
      xlabel2='       Steps        ';ylabel=xlabel
      abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
      errory='n';errorx='n';cone='n';newf='y'
      CALL GRAPH2D(fil2,graphfile2,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
       & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
       & ord_errorp,1,cone,newf,0)
      CALL SORT(mco,temp1D)
      CALL PDLIM(mco,temp1D,bf_val,limmed,limpd)
      medi_val = temp1d(limmed(5))
      CALL HISTO(fil2,graphfile,bf_val,medi_val,xlabel,nbinh)
      CALL PDFSCREEN(texb,texunit,bf_val,medi_val,limpd,gelman,gelmanval,&
       & nchain,jump)
    ENDIF
  ENDDO

  IF(limb.EQ.'qd'.AND.ntr.GT.0)THEN   
    DO i=1,nfi
      IF(fitlimb(i).NE.'n')THEN                     ! quadratic-LD jump parameters 
        DO l=1,nd 
          temp_a1 = '-C'
          temp_a2 = ' =           '
          let1 = CHAR(l+48)
          texb = ' '//wfilter(i)//temp_a1//let1//temp_a2
          fil2 = 'C'//let1//'_'//wfilter(i)//'.res'
          graphfile =  'C'//let1//'_'//wfilter(i)//'_h.sm'
          graphfile2 =  'C'//let1//'_'//wfilter(i)//'.sm'
          texunit = '       ';bf_val=jumplimb(i,l,soluce);temp1D=jumplimb_end(i,l,:)
          xlabel = '          C'//let1//'        '
          IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_jumplimb(i,l)
          IF(ifacf.EQ.'y') THEN
            acffile= 'C' // let1 // '_'//wfilter(i)//'.acf'
            graphacffile='C' // let1 // '_'//wfilter(i)//'_acf.sm'
          ENDIF
          OPEN(UNIT=124,FILE=fil2)
          DO k=1,mco
            WRITE(124,*) k,temp1D(k)
          ENDDO
          CLOSE(124)
          IF(ifacf.EQ.'y')THEN
            nu = (mco/(nchain*2))
            ALLOCATE(acf(nu))
            OPEN(UNIT=135,FILE=acffile)
            DO l2=1,nu
              CALL AUTOCOR(mco,temp1D,nchain,l2-1,acf(l2))
              WRITE(135,*) l2-1,acf(l2)
            ENDDO
            CLOSE(135)
            min_x = 0-nu*0.05;max_x=nu*1.05
            min_y = MINVAL(acf)-(MAXVAL(acf)-MINVAL(acf))*0.05
            max_y = MAXVAL(acf)+(MAXVAL(acf)-MINVAL(acf))*0.05
            ncol_file=2;abs_val=1;ord_val=2
            xlabel2='    Lag (Nsteps)    ';ylabel='   Autocorrelation  '
            abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
            errory='n';errorx='n';cone='n';newf='y'
            CALL GRAPH2D(acffile,graphacffile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
              & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
              & ord_errorp,1,cone,newf,0)
            DEALLOCATE(acf)
          ENDIF
          min_x = 0-mco*0.05;max_x=mco*1.05
          min_y = MINVAL(temp1D)-ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
          max_y = MAXVAL(temp1D)+ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
          ncol_file=2;abs_val=1;ord_val=2
          xlabel2='       Steps        ';ylabel=xlabel
          abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
          errory='n';errorx='n';cone='n';newf='y'
          CALL GRAPH2D(fil2,graphfile2,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
           & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
           & ord_errorp,1,cone,newf,0)
          CALL SORT(mco,temp1D)
          CALL PDLIM(mco,temp1D,bf_val,limmed,limpd)
          medi_val = temp1d(limmed(5))
          CALL HISTO(fil2,graphfile,bf_val,medi_val,xlabel,nbinh)
          CALL PDFSCREEN(texb,texunit,bf_val,medi_val,limpd,gelman,gelmanval, & 
           & nchain,jump)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
 
  DO j=1,npla
    if(.not.(isoch.eq.'y'.and.ntr.eq.0.and.nrv.eq.0))then !not in case of the IsochPlacement only
      WRITE(*,'(A8,1x,I2)') 'PLANET ',j
      WRITE(444,'(A8,1x,I2)') 'PLANET ',j
      WRITE(445,'(A8,1x,I2)') 'PLANET ',j
    end if
    let1=CHAR(48+j)
    DO k=1,11
      test=0
      SELECT CASE (k)
        CASE(1)  ! Transit depth dF
          IF(isjump(j,1).NE.'n')THEN
            fil2 = 'df' // let1 // '.res'; graphfile = 'df'//let1//'_h.sm'
            graphfile2 =  'df'//let1//'.sm'
            xlabel = '          dF        '
            texb =' dF =               ';texunit='       ';bf_val=dF(j,soluce);&
              & temp1D=dF_end(j,:); test=1
            IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_dF(j)
            IF(ifacf.EQ.'y') THEN
              acffile= 'df' // let1 // '.acf'
              graphacffile='df' // let1 // '_acf.sm'
            ENDIF
          ENDIF   
        CASE(2)  ! Tidel elongation TE
          IF(isjump(j,9).NE.'n')THEN 
            fil2 = 'tel' // let1 // '.res'; graphfile = 'tel'//let1//'_h.sm'
            graphfile2 = 'tel'//let1//'.sm'
            xlabel = '         TE         '
            texb =' TE =               ';texunit='       ';bf_val=tidel(j,soluce);& 
              & temp1D=tidel_end(j,:);test=1
            IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_tidel(j)
            IF(ifacf.EQ.'y') THEN
              acffile= 'tel' // let1 // '.acf'
              graphacffile='tel' // let1 // '_acf.sm'
            ENDIF
          ENDIF   
        CASE(3)  ! Transit circular impact parameter
          IF(isjump(j,2).NE.'n')THEN
            fil2 = 'b' // let1 // '.res'; graphfile = 'b'//let1//'_h.sm'
            graphfile2 = 'b'//let1//'.sm'
            xlabel = '        b [R*]      '
            texb =' b =                ';texunit='     R*';bf_val=b(j,soluce);&
              & temp1D=b_end(j,:);test=1 
            IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_b(j)   
            IF(ifacf.EQ.'y') THEN
              acffile= 'b' // let1 // '.acf'
              graphacffile='b' // let1 // '_acf.sm'
            ENDIF  
          ENDIF
        CASE(4)  ! Transit duration W
          IF(isjump(j,3).NE.'n')THEN
            fil2 = 'dur' // let1 // '.res'; graphfile = 'dur'//let1//'_h.sm'
            graphfile2 = 'dur'//let1//'.sm'
            xlabel = '        W [d]       '
            texb =' W =                ';texunit='   days';bf_val=dur(j,soluce);&
              & temp1D=dur_end(j,:);test=1
            IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_dur(j) 
            IF(ifacf.EQ.'y') THEN
              acffile= 'dur' // let1 // '.acf'
              graphacffile='dur' // let1 // '_acf.sm'
            ENDIF
          ENDIF 
        CASE(5)  ! Inferior conjunction time T0
          IF(isjump(j,4).NE.'n')THEN
            fil2 = 't0_' // let1 // '.res'; graphfile = 't0_'//let1//'_h.sm'
            graphfile2 = 't0_'//let1//'.sm'
            xlabel = '       T0 [JD]      '
            texb =' T0 =               ';texunit='     JD';bf_val=t0(j,soluce);&
              & temp1D=t0_end(j,:);test=1
            IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_t0(j) 
            IF(ifacf.EQ.'y') THEN
              acffile= 't0_' // let1 // '.acf'
              graphacffile='t0_' // let1 // '_acf.sm'
            ENDIF
          ENDIF
        CASE(6)  ! Orbital period P
          IF(isjump(j,5).NE.'n')THEN 
            fil2 = 'per' // let1 // '.res'; graphfile = 'per'//let1//'_h.sm'
            graphfile2 = 'per'//let1//'.sm'
            xlabel = '        P [d]       '
            texb =' P =                ';texunit='   days';bf_val=per(j,soluce);& 
              & temp1D=per_end(j,:);test=1
            IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_per(j)
            IF(ifacf.EQ.'y') THEN
              acffile= 'per' // let1 // '.acf'
              graphacffile='per' // let1 // '_acf.sm'
            ENDIF
          ENDIF 
        CASE(7)  ! sqrt(e)cosw 
          IF(isjump(j,6).NE.'n')THEN
            fil2 = 'secosw' // let1 // '.res'; graphfile = 'secosw'//let1//'_h.sm'
            graphfile2 = 'secosw'//let1//'.sm'
            xlabel = '     sqrt(e)cosw    '
            texb =' sqrt(e)cosw =      ';texunit='       ';bf_val=secosw(j,soluce);&
              & temp1D=secosw_end(j,:);test=1
            IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_secosw(j)
            IF(ifacf.EQ.'y') THEN
              acffile= 'secosw' // let1 // '.acf'
              graphacffile='secosw' // let1 // '_acf.sm'
            ENDIF
          ENDIF   
        CASE(8)  ! sqrt(e)sinw
          IF(isjump(j,6).NE.'n')THEN 
            fil2 = 'sesinw' // let1 // '.res'; graphfile = 'sesinw'//let1//'_h.sm'
            graphfile2 = 'sesinw'//let1//'.sm'
            xlabel = '     sqrt(e)sinw    '
            texb =' sqrt(e)sinw =      ';texunit = '       '
            bf_val=sesinw(j,soluce);temp1D=sesinw_end(j,:);test=1
            IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_sesinw(j)
            IF(ifacf.EQ.'y') THEN
              acffile= 'sesinw' // let1 // '.acf'
              graphacffile='sesinw' // let1 // '_acf.sm'
            ENDIF
          ENDIF    
        CASE(9)  ! K2 
          IF(isjump(j,7).NE.'n')THEN
            fil2 = 'k2_' // let1 // '.res'; graphfile = 'k2_'//let1//'_h.sm'
            graphfile2 = 'k2_'//let1//'.sm'
            xlabel = '         K2         '
            texb =' K2 =               ';texunit = '       ';bf_val=kb(j,soluce);&
              & temp1D=kb_end(j,:);test=1
            IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_kb(j) 
            IF(ifacf.EQ.'y') THEN
              acffile= 'k2_' // let1 // '.acf'
              graphacffile='k2_' // let1 // '_acf.sm'
            ENDIF
          ENDIF
        CASE(10)  ! sqrt(VsinI)cosB 
          IF(isjump(j,8).NE.'n'.AND.j.EQ.1.AND.nrv.GT.0)THEN
            fil2 = 'svsinicosb' // let1 // '.res'
            graphfile = 'svsinicosb'//let1//'_h.sm'
            graphfile2 = 'svsinicosb'//let1//'.sm'
            xlabel = '   sqrt(VsinI)cosB  '
            texb =' sqrt(VsinI)cosB =  ';texunit = '       '
            bf_val=svsinicosbeta(soluce);temp1D=svsinicosbeta_end;test=1
            IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_svsinicosbeta   
            IF(ifacf.EQ.'y') THEN
              acffile= 'svsinicosb' // let1 // '.acf'
              graphacffile='svsinicosb' // let1 // '_acf.sm'
            ENDIF
          ENDIF
        CASE(11)  ! sqrt(VsinI)sinB 
          IF(isjump(j,8).NE.'n'.AND.j.EQ.1.AND.nrv.GT.0)THEN
            fil2 = 'svsinisinb' // let1 // '.res'
            graphfile = 'svsinisinb'//let1//'_h.sm'
            graphfile2 = 'svsinisinb'//let1//'.sm'
            xlabel = '   sqrt(VsinI)sinB  '
            texb =' sqrt(VsinI)sinB =  ';texunit = '       '
            bf_val=svsinisinbeta(soluce);temp1D=svsinisinbeta_end;test=1
            IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_svsinisinbeta    
            IF(ifacf.EQ.'y') THEN
              acffile= 'svsinisinb' // let1 // '.acf'
              graphacffile='svsinisinb' // let1 // '_acf.sm'
            ENDIF
          ENDIF
      END SELECT
      IF(test.EQ.1)THEN
        OPEN(UNIT=124,FILE=fil2)
        DO l=1,mco
          WRITE(124,*) l,temp1D(l)
        ENDDO
        CLOSE(124)
        IF(ifacf.EQ.'y')THEN
          nu = (mco/(nchain*2))
          ALLOCATE(acf(nu))
          OPEN(UNIT=135,FILE=acffile)
          DO l=1,nu
            CALL AUTOCOR(mco,temp1D,nchain,l-1,acf(l))
            WRITE(135,*) l-1,acf(l)
          ENDDO
          CLOSE(135)
          min_x = 0-nu*0.05;max_x=nu*1.05
          min_y = MINVAL(acf)-(MAXVAL(acf)-MINVAL(acf))*0.05
          max_y = MAXVAL(acf)+(MAXVAL(acf)-MINVAL(acf))*0.05
          ncol_file=2;abs_val=1;ord_val=2
          xlabel2='    Lag (Nsteps)    ';ylabel='   Autocorrelation  '
          abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
          errory='n';errorx='n';cone='n';newf='y'
          CALL GRAPH2D(acffile,graphacffile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
            & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
            & ord_errorp,1,cone,newf,0)
          DEALLOCATE(acf)
        ENDIF
        min_x = 0-mco*0.05;max_x=mco*1.05
        min_y = MINVAL(temp1D)-ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
        max_y = MAXVAL(temp1D)+ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
        ncol_file=2;abs_val=1;ord_val=2
        xlabel2='       Steps        ';ylabel=xlabel
        abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
        errory='n';errorx='n';cone='n';newf='y'
        CALL GRAPH2D(fil2,graphfile2,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
         & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
         & ord_errorp,1,cone,newf,0)
        CALL SORT(mco,temp1D)
        CALL PDLIM(mco,temp1D,bf_val,limmed,limpd)
        medi_val = temp1d(limmed(5))
        CALL HISTO(fil2,graphfile,bf_val,medi_val,xlabel,nbinh)
        CALL PDFSCREEN(texb,texunit,bf_val,medi_val,limpd,gelman,gelmanval, & 
         & nchain,jump)
      ENDIF
    ENDDO
    IF(ntr.GT.0)THEN
      let1 = CHAR(48+INT(j/10))
      DO i=1,nfi
        DO k=1,5 
           test=0
           SELECT CASE(k)
             CASE(1) ! Occultation depth
               IF(iffitoc(i).NE.'n')THEN
                 fil2 = 'dFoc_' // let1 // wfilter(i) // '.res'
                 graphfile = 'dFoc_' // let1 // wfilter(i)//'_h.sm'
                 graphfile2 = 'dFoc_' // let1 //wfilter(i)//'.sm'
                 xlabel = '       '//wfilter(i)//'-dFoc      '
                 texb =' '//wfilter(i)//'-dFoc =          ';texunit='       ';&
                   & bf_val=dFsec(i,soluce,j);temp1D=dFsec_end(i,:,j);test=1
                 IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_dFsec(i,j)
                 IF(ifacf.EQ.'y') THEN
                   acffile = 'dFoc_'//let1//wfilter(i)//'.acf'
                   graphacffile = 'dFoc_'//let1//wfilter(i)//'_acf.sm'
                 ENDIF
               ENDIF  
             CASE(2) ! Phase curve amplitude 1
               IF(iffitph(i).NE.'n'.AND.j.EQ.1)THEN
                 fil2 = 'PhA1_' // wfilter(i) // '.res'
                 graphfile = 'PhA1_'//wfilter(i)//'_h.sm'
                 graphfile2 = 'PhA1_'//wfilter(i)//'.sm'
                 xlabel = '       '//wfilter(i)//'-PhA1      '
                 texb =' '//wfilter(i)//'-Ph_A1 =         ';texunit='       ';&
                   & bf_val=phampli1(i,soluce);temp1D=phampli1_end(i,:);test=1
                 IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_phampli1(i)
                 IF(ifacf.EQ.'y') THEN
                   acffile = 'PhA1_'//wfilter(i)//'.acf'
                   graphacffile = 'PhA1_'//wfilter(i)//'_acf.sm'
                 ENDIF
               ENDIF  
             CASE(3) ! Phase curve amplitude 2
               IF(iffitph(i).NE.'n'.AND.j.EQ.1)THEN
                 fil2 = 'PhA2_' // wfilter(i) // '.res'
                 graphfile = 'PhA2_'//wfilter(i)//'_h.sm'
                 graphfile2 = 'PhA2_'//wfilter(i)//'.sm'
                 xlabel = '       '//wfilter(i)//'-PhA2      '
                 texb =' '//wfilter(i)//'-Ph_A2 =         ';texunit='       ';&
                   & bf_val=phampli2(i,soluce);temp1D=phampli2_end(i,:);test=1
                 IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_phampli2(i)
                 IF(ifacf.EQ.'y') THEN
                   acffile = 'PhA2_'//wfilter(i)//'.acf'
                   graphacffile = 'PhA2_'//wfilter(i)//'_acf.sm'
                 ENDIF
               ENDIF  
             CASE(4) ! Phase curve amplitude 3
               IF(iffitph(i).NE.'n'.AND.j.EQ.1)THEN
                 fil2 = 'PhA3_' // wfilter(i) // '.res'
                 graphfile = 'PhA3_'//wfilter(i)//'_h.sm'
                 graphfile2 = 'PhA3_'//wfilter(i)//'.sm'
                 xlabel = '       '//wfilter(i)//'-PhA3      '
                 texb =' '//wfilter(i)//'-Ph_A3 =         ';texunit='       ';&
                   & bf_val=phampli3(i,soluce);temp1D=phampli3_end(i,:);test=1
                 IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_phampli3(i)  
                 IF(ifacf.EQ.'y') THEN
                   acffile = 'PhA3_'//wfilter(i)//'.acf'
                   graphacffile = 'PhA3_'//wfilter(i)//'_acf.sm'
                 ENDIF
               ENDIF
             CASE(5) ! Phase curve offset
               IF(iffitph(i).NE.'n'.AND.j.EQ.1)THEN
                 fil2 = 'PhOf_' // wfilter(i) // '.res'
                 graphfile = 'PhOf_'//wfilter(i)//'_h.sm'
                 graphfile2 = 'PhOf_'//wfilter(i)//'.sm'
                 xlabel = '    '//wfilter(i)//'-PhOf [deg]   '
                 texb =' '//wfilter(i)//'-PH-Of =         ';texunit='    deg';&
                   & bf_val=phoffset(i,soluce);temp1D=phoffset_end(i,:);test=1
                 IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_phoffset(i)
                 IF(ifacf.EQ.'y') THEN
                   acffile = 'PhOf_'//wfilter(i)//'.acf'
                   graphacffile = 'PhOf_'//wfilter(i)//'_acf.sm'
                 ENDIF
               ENDIF  
           END SELECT
           IF(test.EQ.1)THEN
             OPEN(UNIT=124,FILE=fil2)
             DO l=1,mco
               WRITE(124,*) l,temp1D(l)
             ENDDO
             CLOSE(124)
             IF(ifacf.EQ.'y')THEN
               nu = (mco/(nchain*2))
               ALLOCATE(acf(nu))
               OPEN(UNIT=135,FILE=acffile)
               DO l=1,nu
                 CALL AUTOCOR(mco,temp1D,nchain,l-1,acf(l))
                 WRITE(135,*) l-1,acf(l)
               ENDDO
               CLOSE(135)
               min_x = 0-nu*0.05;max_x=nu*1.05
               min_y = MINVAL(acf)-(MAXVAL(acf)-MINVAL(acf))*0.05
               max_y = MAXVAL(acf)+(MAXVAL(acf)-MINVAL(acf))*0.05
               ncol_file=2;abs_val=1;ord_val=2
               xlabel2='    Lag (Nsteps)    ';ylabel='   Autocorrelation  '
               abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
               errory='n';errorx='n';cone='n';newf='y'
               CALL GRAPH2D(acffile,graphacffile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
                 & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
                 & ord_errorp,1,cone,newf,0)
               DEALLOCATE(acf)
             ENDIF
             min_x = 0-mco*0.05;max_x=mco*1.05
             min_y = MINVAL(temp1D)-ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
             max_y = MAXVAL(temp1D)+ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
             ncol_file=2;abs_val=1;ord_val=2
             xlabel2='       Steps        ';ylabel=xlabel
             abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
             errory='n';errorx='n';cone='n';newf='y'
             CALL GRAPH2D(fil2,graphfile2,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
              & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
              & ord_errorp,1,cone,newf,0)
             CALL SORT(mco,temp1D)
             CALL PDLIM(mco,temp1D,bf_val,limmed,limpd)
             medi_val = temp1d(limmed(5))
             CALL HISTO(fil2,graphfile,bf_val,medi_val,xlabel,nbinh)
             CALL PDFSCREEN(texb,texunit,bf_val,medi_val,limpd,gelman,gelmanval, & 
              & nchain,jump)
          ENDIF
        ENDDO
      ENDDO
    ENDIF

    IF(nddf.GT.0.AND.ntr.GT.0.AND.isddf.NE.'n')THEN  ! Transit depth differences
      DO i=1,nddf
        texb =' '//ddf_filter(i)//'-DdF =           ';texunit='       ';&
          & bf_val=ddf(j,i,soluce);temp1D=ddf_end(j,i,:)
        fil2 = 'Ddf_' // let1// '_' //ddf_filter(i) // '.res'
        graphfile = 'Ddf_'// let1 // '_' // ddf_filter(i)//'_h.sm'
        graphfile2 = 'Ddf_'// let1 // '_' // ddf_filter(i)//'.sm'
        xlabel = '         DdF        '
        IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_ddf(j,i)  
        IF(ifacf.EQ.'y') THEN
          acffile = 'Ddf_'// let1 // '_' // ddf_filter(i)//'.acf'
          graphacffile = 'Ddf_'// let1 // '_' // ddf_filter(i)//'_acf.sm'
        ENDIF
        OPEN(UNIT=124,FILE=fil2)
        DO l=1,mco
          WRITE(124,*) l,temp1D(l)
        ENDDO
        CLOSE(124)
        IF(ifacf.EQ.'y')THEN
          nu = (mco/(nchain*2))
          ALLOCATE(acf(nu))
          OPEN(UNIT=135,FILE=acffile)
          DO l=1,nu
            CALL AUTOCOR(mco,temp1D,nchain,l-1,acf(l))
            WRITE(135,*) l-1,acf(l)
          ENDDO
          CLOSE(135)
          min_x = 0-nu*0.05;max_x=nu*1.05
          min_y = MINVAL(acf)-(MAXVAL(acf)-MINVAL(acf))*0.05
          max_y = MAXVAL(acf)+(MAXVAL(acf)-MINVAL(acf))*0.05
          ncol_file=2;abs_val=1;ord_val=2
          xlabel2='    Lag (Nsteps)    ';ylabel='   Autocorrelation  '
          abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
          errory='n';errorx='n';cone='n';newf='y'
          CALL GRAPH2D(acffile,graphacffile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
            & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
            & ord_errorp,1,cone,newf,0)
          DEALLOCATE(acf)
        ENDIF
        min_x = 0-mco*0.05;max_x=mco*1.05
        min_y = MINVAL(temp1D)-ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
        max_y = MAXVAL(temp1D)+ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
        ncol_file=2;abs_val=1;ord_val=2
        xlabel2='       Steps        ';ylabel=xlabel
        abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
        errory='n';errorx='n';cone='n';newf='y'
        CALL GRAPH2D(fil2,graphfile2,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
         & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
         & ord_errorp,1,cone,newf,0)
        CALL SORT(mco,temp1D)
        CALL PDLIM(mco,temp1D,bf_val,limmed,limpd)
        medi_val = temp1d(limmed(5))
        CALL HISTO(fil2,graphfile,bf_val,medi_val,xlabel,nbinh)
        CALL PDFSCREEN(texb,texunit,bf_val,medi_val,limpd,gelman,gelmanval, & 
         & nchain,jump)
      ENDDO
    ENDIF

    IF(ngroup.GT.0)THEN  ! Transit depth white LCs for each group
      DO i=1,ngroup
        let1 = CHAR(48+INT(i/1000))
        i2 = i - INT(i/1000)*1000
        let2 = CHAR(48+INT(i2/100))
        i2 = i - INT(i/100)*100
        let3 = CHAR(48+INT(i2/10))
        let4 = CHAR(48+MOD(i2,10))
        texb ='G'//let4//'-GdF =           ';texunit='       ';&
          & bf_val=dfgroup(i,soluce);temp1D=dfgroup_end(i,:)
        fil2 = 'Gdf_' // let4 // '.res'
        graphfile = 'Gdf_'//let4//'_h.sm'
        graphfile2 = 'Gdf_'//let4//'.sm'
        xlabel = '          dF        '
        IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_dfgroup(i)  
        IF(ifacf.EQ.'y') THEN
          acffile = 'Gdf_'//let4//'.acf'
          graphacffile = 'Gdf_'//let4//'_acf.sm'
        ENDIF
        OPEN(UNIT=124,FILE=fil2)
        DO l=1,mco
          WRITE(124,*) l,temp1D(l)
        ENDDO
        CLOSE(124)
        IF(ifacf.EQ.'y')THEN
          nu = (mco/(nchain*2))
          ALLOCATE(acf(nu))
          OPEN(UNIT=135,FILE=acffile)
          DO l=1,nu
            CALL AUTOCOR(mco,temp1D,nchain,l-1,acf(l))
            WRITE(135,*) l-1,acf(l)
          ENDDO
          CLOSE(135)
          min_x = 0-nu*0.05;max_x=nu*1.05
          min_y = MINVAL(acf)-(MAXVAL(acf)-MINVAL(acf))*0.05
          max_y = MAXVAL(acf)+(MAXVAL(acf)-MINVAL(acf))*0.05
          ncol_file=2;abs_val=1;ord_val=2
          xlabel2='    Lag (Nsteps)    ';ylabel='   Autocorrelation  '
          abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
          errory='n';errorx='n';cone='n';newf='y'
          CALL GRAPH2D(acffile,graphacffile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
            & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
            & ord_errorp,1,cone,newf,0)
          DEALLOCATE(acf)
        ENDIF
        min_x = 0-mco*0.05;max_x=mco*1.05
        min_y = MINVAL(temp1D)-ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
        max_y = MAXVAL(temp1D)+ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
        ncol_file=2;abs_val=1;ord_val=2
        xlabel2='       Steps        ';ylabel=xlabel
        abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
        errory='n';errorx='n';cone='n';newf='y'
        CALL GRAPH2D(fil2,graphfile2,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
         & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
         & ord_errorp,1,cone,newf,0)
        CALL SORT(mco,temp1D)
        CALL PDLIM(mco,temp1D,bf_val,limmed,limpd)
        medi_val = temp1d(limmed(5))
        CALL HISTO(fil2,graphfile,bf_val,medi_val,xlabel,nbinh)
        CALL PDFSCREEN(texb,texunit,bf_val,medi_val,limpd,gelman,gelmanval, & 
         & nchain,jump)
      ENDDO
    ENDIF

    IF(nttvmax.GT.0.AND.isttv.EQ.'y')THEN  ! Transit timing variations 
      let1 = CHAR(48+int(j/10))
      let2 = CHAR(48+mod(j,10))
      name = 'ttv' // let1 // let2 // '.res'
      OPEN(UNIT=125,FILE=name)
      ttvrms(j)=0.; ave=0.; ce=0
      DO l=1,nttv(j)
        IF(l.LT.10)THEN  
           let3=CHAR(48+l)
           texb =' TTV ' // let3 // ' =   '
           fil2 = 'ttv' // let1 // let2 // '_' // let3 // '.res'
           graphfile = 'ttv' // let1 // let2 // '_' // let3 // '_h.sm'
           graphfile2 = 'ttv' // let1 // let2 // '_' // let3 // '.sm'
        ELSE
           let3 = CHAR(48+INT(l/10)); let4 = CHAR(48+MOD(l,10))
           texb =' TTV ' // let3 // let4 // ' =   '
           fil2 = 'ttv' // let1 // let2 // '_' // let3 // let4 // '.res'
           graphfile = 'ttv' // let1 // let2 // '_' // let3 // let4 // '_h.sm'
           graphfile2 = 'ttv' // let1 // let2 // '_' // let3 // let4 // '.sm'
        ENDIF
        xlabel = '      TTV [min]     '
        texunit='    min';bf_val=ttv(j,l,soluce)*24*60.;temp1D=ttv_end(j,l,:)*24*60.
        IF(gelman.EQ.'y'.AND.nchain.GT.1) gelmanval=gelman_ttv(j,l)
        IF(ifacf.EQ.'y'.AND.l.LT.10) THEN
          acffile = 'ttv'// let1 // let2 // '_' // let3 //'.acf'
          graphacffile = 'ttv'// let1 // let2 // '_' // let3 //'_acf.sm'
        ELSE IF(ifacf.EQ.'y'.AND.l.GE.10) THEN
          acffile = 'ttv'// let1 // let2 // '_' // let3 // let4 //'.acf'
          graphacffile = 'ttv'// let1 // let2 // '_' // let3 // let4 // '_acf.sm'
        ENDIF
        OPEN(UNIT=124,FILE=fil2)
        DO k=1,mco
          WRITE(124,*) k,temp1D(k)
        ENDDO
        CLOSE(124)
        IF(ifacf.EQ.'y')THEN
          nu = (mco/(nchain*2))
          ALLOCATE(acf(nu))
          OPEN(UNIT=135,FILE=acffile)
          DO k=1,nu
            CALL AUTOCOR(mco,temp1D,nchain,k-1,acf(k))
            WRITE(135,*) k-1,acf(k)
          ENDDO
          CLOSE(135)
          min_x = 0-nu*0.05;max_x=nu*1.05
          min_y = MINVAL(acf)-(MAXVAL(acf)-MINVAL(acf))*0.05
          max_y = MAXVAL(acf)+(MAXVAL(acf)-MINVAL(acf))*0.05
          ncol_file=2;abs_val=1;ord_val=2
          xlabel2='    Lag (Nsteps)    ';ylabel='   Autocorrelation  '
          abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
          errory='n';errorx='n';cone='n';newf='y'
          CALL GRAPH2D(acffile,graphacffile,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
            & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
            & ord_errorp,1,cone,newf,0)
          DEALLOCATE(acf)
        ENDIF
        min_x = 0-mco*0.05;max_x=mco*1.05
        min_y = MINVAL(temp1D)-ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
        max_y = MAXVAL(temp1D)+ABS(MAXVAL(temp1D)-MINVAL(temp1D))*0.05
        ncol_file=2;abs_val=1;ord_val=2
        xlabel2='       Steps        ';ylabel=xlabel
        abs_errorm=0;abs_errorp=0;ord_errorm=0;ord_errorp=0
        errory='n';errorx='n';cone='n';newf='y'
        CALL GRAPH2D(fil2,graphfile2,ncol_file,abs_val,ord_val,xlabel2,ylabel, & 
         & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
         & ord_errorp,1,cone,newf,0)
        CALL SORT(mco,temp1D)
        CALL PDLIM(mco,temp1D,bf_val,limmed,limpd)
        compttv(j,l)=temp1d(limmed(5))
        ave=ave+compttv(j,l); ce=ce+1
        WRITE(125,*) epochtr(j,l),compttv(j,l),limpd(5),limpd(6)
        medi_val = temp1d(limmed(5))
        CALL HISTO(fil2,graphfile,bf_val,medi_val,xlabel,nbinh)
        CALL PDFSCREEN(texb,texunit,bf_val,medi_val,limpd,gelman,gelmanval, & 
           & nchain,jump)
      ENDDO
      ave=ave/DBLE(ce)
      DO l=1,nttv(j)
        ttvrms(j)=ttvrms(j)+(compttv(j,l)-ave)**2
      ENDDO
      ttvrms(j)=SQRT(ttvrms(j)/ce)
      CLOSE(125)
      das = (MAXVAL(epochtr(j,:))-MINVAL(epochtr(j,:)))/10.
      min_x = MINVAL(epochtr(j,:))-das; max_x = MAXVAL(epochtr(j,:))+das
      das = 5*ttvrms(j)
      min_y = MINVAL(compttv(j,:))-das; max_y = MAXVAL(compttv(j,:))+das
      inputfile='ttv'//let1//let2//'.res';graphfile='ttv'//let1//let2//'.sm';ncol_file=4 
      abs_val=1;ord_val=2;xlabel='Epoch';ylabel='O-C (min)'
      abs_errorm=0;abs_errorp=0;ord_errorm=3;ord_errorp=4
      errory='y';errorx='n';cone='n';newf='y'
      CALL GRAPH2D(inputfile,graphfile,ncol_file,abs_val,ord_val,xlabel,ylabel, & 
        & errorx,errory,min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm, & 
        & ord_errorp,1,cone,newf,0)
    ENDIF
  
  ENDDO

  PRINT*,      '------------------------------------------------------------------------------------------'
  WRITE(444,*) '------------------------------------------------------------------------------------------'
  WRITE(445,*) '------------------------------------------------------------------------------------------'

  PRINT*, 'Derived parameters: median of the MPDF + 68.3 and 99.7% limits' 
  WRITE(444,*) 'Derived parameters: median of the MPDF + 68.3 and 99.7% limits' 
  WRITE(445,*) 'Derived parameters: median of the MPDF + 68.3 and 99.7% limits'
  jump='n'
 
  WRITE(*,'(A5)') 'STAR'
  WRITE(444,'(A5)') 'STAR'  
  WRITE(445,'(A5)') 'STAR'

  DO k=1,13
    test=0
      SELECT CASE (k)
        CASE(1)  ! logg
          fil2 = 'logg.res'; graphfile = 'logg_h.sm'
          xlabel = '      log(g)        '
          texb =' logg =             ';texunit='       '
          bf_val=logg(soluce); temp1D=logg_end(:); test=1
        CASE(2)  ! Ktide
          IF(testf2.EQ.'y')THEN
            fil2 = 'ktide.res'; graphfile = 'ktide_h.sm'
            xlabel = '     Ktide [m/s]    '
            texb =' Ktide =            ';texunit='   ms-1'
            bf_val=ktide(soluce); temp1D=ktide_end(:); test=1
          ENDIF
        CASE(3)  ! rho*
          fil2 = 'rhos.res'; graphfile = 'rhos_h.sm'
          xlabel = '  rho* [rho(sun)]   '
          texb =' rho* =             ';texunit='rho_sun'
          bf_val=rho(soluce); temp1D=rho_end(:); test=1
        CASE(4)  ! M*
          IF(fixstellar.EQ.'y'.OR.enoch.EQ.'y'.OR.msrs.EQ.'y'.OR.massfromr.EQ.'y'.or.isoch.ne.'y')THEN !eq
            fil2 = 'ms.res'; graphfile = 'ms_h.sm'
            xlabel = '     M* [Msun]      '
            texb =' M* =               ';texunit='  M_sun'
            bf_val=mass_s(soluce); temp1D=mass_s_end(:); test=1
          ENDIF
        CASE(5)  ! R*
          IF(fixstellar.EQ.'y'.OR.(massfromr.NE.'y'.AND.fitmsrs.NE.'y'.and.RjumpIso.eq.'n'))THEN
            fil2 = 'rs.res'; graphfile = 'rs_h.sm'
            xlabel = '     R* [Rsun]      '
            texb =' R* =               ';texunit='  R_sun'
            bf_val=radius_s(soluce); temp1D=radius_s_end(:); test=1
          ENDIF
        CASE(6)  ! Teff
          IF(fixstellar.EQ.'y')THEN
            fil2 = 'teff.res'; graphfile = 'teff_h.sm'
            xlabel = '      Teff [K]      '
            texb =' Teff =             ';texunit='      K'
            bf_val=temp_s(soluce); temp1D=temp_s_end(:); test=1
          ENDIF
        CASE(7)  ! Fe/H
          IF(fixstellar.EQ.'y')THEN
            fil2 = 'feh.res'; graphfile = 'feh_h.sm'
            xlabel = '       [Fe/H]       '
            texb =' [Fe/H] =           ';texunit='       '
            bf_val=met_s(soluce); temp1D=met_s_end(:); test=1
          ENDIF
        CASE(8)  ! L*
          fil2 = 'lums.res'; graphfile = 'lums_h.sm'
          xlabel = '     L* [Lsun]      '
          texb =' L* =                 ';texunit='  L_sun'
          bf_val=lum_s(soluce); temp1D=lum_s_end(:); test=1
        CASE(9)  ! Age
          if (isoch.eq.'y') then
			fil2 = 'age.res'; graphfile = 'age_h.sm'
			xlabel = '      age           '
			texb =' Age =              ';texunit='    Gyr'
			bf_val=age(soluce)
			call find(.not.isEq_v(age_end(:),0.D0,1),nda) !skip flagged values
			if (allocated(nda)) then
			  print*,'size(nda)',size(nda)
			  if (size(nda).gt.size(age_end)/2.) then !flagged ages are a minority
			    !then show the age. Otherwise the star is to faint and no reasonable
			    !isochronal age may be provided
			    temp1D=age_end(nda)
			    test=1
			  end if
			end if
			print*,'testAge',test
		  end if
        CASE(10)  ! VsinI*
          fil2 = 'vsini.res'; graphfile = 'vsini_h.sm'
          xlabel = '   VsinI* [km/s]    '
          texb =' VsinI* =           ';texunit='   km/s'
          bf_val=vsini(soluce); temp1D=vsini_end(:); test=1
        CASE(11) ! SinI*
          IF(stelincli.EQ.'y')THEN
            fil2 = 'sininclis.res'; graphfile = 'sininclis_h.sm'
            xlabel = '        SinI*       '
            texb =' SinI* =            ';texunit='       '
            bf_val=sinincli_s(soluce); temp1D=sinincli_s_end(:); test=1
          ENDIF
        CASE(12) ! I*
          IF(stelincli.EQ.'y')THEN
            fil2 = 'inclis.res'; graphfile = 'inclis_h.sm'
            xlabel = '      I* [deg]      '
            texb =' I* =               ';texunit='    deg'
            bf_val=incli_s(soluce); temp1D=incli_s_end(:); test=1
          ENDIF
        CASE(13) ! M*_jitter
          IF(enoch.EQ.'y')THEN
            fil2 = 'mjitter.res'; graphfile = 'mjitters_h.sm'
            xlabel = '  M* jitter [Msun]  '
            texb =' M* jitter =        ';texunit='  M_sun'
            bf_val=massjitter(soluce); temp1D=massjitter_end(:); test=1
          ENDIF
     END SELECT
     IF(test.EQ.1)THEN
       if (k.eq.9) then !age
          mco2=size(nda)
       else
          mco2=mco
       end if
       medis = INT(mco2/2)
	   limmeds(5) = medis
	   limmeds(1) = CEILING((mco2*(1.-0.683))/2.)
	   limmeds(2) = limmeds(1) + FLOOR(mco2*0.683)
	   limmeds(3) = CEILING((mco2*(1.-0.997))/2.)
	   limmeds(4) = limmeds(3) + FLOOR(mco2*0.997)
       OPEN(UNIT=124,FILE=fil2)
       DO l=1,mco2
          WRITE(124,*) l,temp1D(l)
       ENDDO
       CLOSE(124)    
       CALL SORT(mco2,temp1D)
       CALL PDLIM(mco2,temp1D,bf_val,limmeds,limpd)
       medi_val = temp1d(limmeds(5))
       CALL HISTO(fil2,graphfile,bf_val,medi_val,xlabel,nbinh)
       CALL PDFSCREEN(texb,texunit,bf_val,medi_val,limpd,gelman,gelmanval, & 
        & nchain,jump)
    ENDIF
  ENDDO

  IF(limb.EQ.'qd'.AND.ntr.GT.0)THEN   
    DO i=1,nfi
      IF(fitlimb(i).NE.'n')THEN                     ! quadratic-LD coefficients'
        DO l=1,nd 
          temp_a1 = '-U'
          temp_a2 = ' =           '
          let1 = CHAR(l+48)
          texb = ' '//wfilter(i)//temp_a1//let1//temp_a2
          fil2 = 'U'//let1//'_'//wfilter(i)//'.res'
          graphfile =  'U'//let1//'_'//wfilter(i)//'_h.sm'
          texunit = '       ';bf_val=ql(i,l,soluce);temp1D=ql_end(i,l,:)
          xlabel = '         U'//let1//'         '
          OPEN(UNIT=124,FILE=fil2)
          DO k=1,mco
            WRITE(124,*) k,temp1D(k)
          ENDDO
          CLOSE(124)
          CALL SORT(mco,temp1D)
          CALL PDLIM(mco,temp1D,bf_val,limmed,limpd)
          medi_val = temp1d(limmed(5))
          CALL HISTO(fil2,graphfile,bf_val,medi_val,xlabel,nbinh)
          CALL PDFSCREEN(texb,texunit,bf_val,medi_val,limpd,gelman,gelmanval, & 
           & nchain,jump)
        ENDDO
      ENDIF
    ENDDO
  ENDIF
  
  if (ntr.gt.0.or.nrv.gt.0) then
	  DO j=1,npla
		WRITE(*,'(A8,1x,I2)') 'PLANET ',j
		WRITE(444,'(A8,1x,I2)') 'PLANET ',j
		WRITE(445,'(A8,1x,I2)') 'PLANET ',j
		let1=CHAR(48+j)
		DO k=1,42
		  test=0
		  SELECT CASE (k)
		    CASE(1) ! Transit depth
		      IF(isjump(j,1).EQ.'n')THEN
		        fil2 = 'df' // let1 // '.res'; graphfile = 'df'//let1//'_h.sm'
		        xlabel = '         dF         '
		        texb =' dF =               ';texunit='       ';bf_val=dF(j,soluce);&
		          & temp1D=dF_end(j,:); test=1
		      ENDIF   
		    CASE(2) ! Tidal elongation
		      IF(isjump(j,1).EQ.'n')THEN
		        fil2 = 'tel' // let1 // '.res'; graphfile = 'tel'//let1//'_h.sm'
		        xlabel = '         TE         '
		        texb =' TE =               ';texunit='       ';bf_val=tidel(j,soluce);&
		          & temp1D=tidel_end(j,:); test=1
		      ENDIF        
		    CASE(3)  ! Transit circular impact parameter
		      IF(isjump(j,2).EQ.'n')THEN
		        fil2 = 'b' // let1 // '.res'; graphfile = 'b'//let1//'_h.sm'
		        xlabel = '        b [R*]      '
		        texb =' b =                ';texunit='     R*';bf_val=b(j,soluce);&
		          & temp1D=b_end(j,:);test=1 
		      ENDIF
		    CASE(4)  ! Transit duration W
		      IF(isjump(j,3).EQ.'n')THEN
		        fil2 = 'dur' // let1 // '.res'; graphfile = 'dur'//let1//'_h.sm'
		        xlabel = '        W [d]       '
		        texb =' W =                ';texunit='   days';bf_val=dur(j,soluce);&
		          & temp1D=dur_end(j,:);test=1
		      ENDIF 
		    CASE(5)  ! Inferior conjunction time
		      IF(isjump(j,4).EQ.'n')THEN
		        fil2 = 't0_' // let1 // '.res'; graphfile = 't0_'//let1//'_h.sm'
		        xlabel = '       T0 [JD]      '
		        texb =' T0 =               ';texunit='     JD';bf_val=t0(j,soluce);&
		          & temp1D=t0_end(j,:);test=1
		      ENDIF
		    CASE(6)  ! Orbital period P
		      IF(isjump(j,5).EQ.'n')THEN 
		        fil2 = 'per' // let1 // '.res'; graphfile = 'per'//let1//'_h.sm'
		        xlabel = '        P [d]       '
		        texb =' P =                ';texunit='   days';bf_val=per(j,soluce);& 
		          & temp1D=per_end(j,:);test=1
		      ENDIF 
		    CASE(7)  ! sqrt(e)cosw 
		      IF(isjump(j,6).EQ.'n')THEN
		        fil2 = 'secosw' // let1 // '.res'; graphfile = 'secosw'//let1//'_h.sm'
		        xlabel = '     sqrt(e)cosw    '
		        texb =' sqrt(e)cosw =      ';texunit='       ';bf_val=secosw(j,soluce);&
		          & temp1D=secosw_end(j,:);test=1
		      ENDIF   
		    CASE(8)  ! sqrt(e)sinw
		      IF(isjump(j,6).EQ.'n')THEN 
		        fil2 = 'sesinw' // let1 // '.res'; graphfile = 'sesinw'//let1//'_h.sm'
		        xlabel = '     sqrt(e)sinw    '
		        texb =' sqrt(e)sinw =      ';texunit = '       '
		        bf_val=sesinw(j,soluce);temp1D=sesinw_end(j,:);test=1
		      ENDIF    
		    CASE(9)  ! K2 
		      IF(isjump(j,7).EQ.'n')THEN
		        fil2 = 'k2_' // let1 // '.res'; graphfile = 'k2_'//let1//'_h.sm'
		        xlabel = '         K2         '
		        texb =' K2 =               ';texunit = '       ';bf_val=kb(j,soluce);&
		          & temp1D=kb_end(j,:);test=1
		      ENDIF
		    CASE(10)  ! sqrt(VsinI)cosB 
		      IF(isjump(j,8).EQ.'n'.AND.j.EQ.1.AND.nrv.GT.0)THEN
		        fil2 = 'svsinicosb' // let1 // '.res'
		        graphfile = 'svsinicosb'//let1//'_h.sm'
		        xlabel = '   sqrt(VsinI)cosB  '
		        texb =' sqrt(VsinI)cosB =  ';texunit = '       '
		        bf_val=svsinicosbeta(soluce);temp1D=svsinicosbeta_end;test=1
		      ENDIF
		    CASE(11)  ! sqrt(VsinI)sinB 
		      IF(isjump(j,8).EQ.'n'.AND.j.EQ.1.AND.nrv.GT.0)THEN
		        fil2 = 'svsinisinb' // let1 // '.res'
		        graphfile = 'svsinisinb'//let1//'_h.sm'
		        xlabel = '   sqrt(VsinI)sinB  '
		        texb =' sqrt(VsinI)sinB =  ';texunit = '       '
		        bf_val=svsinisinbeta(soluce);temp1D=svsinisinbeta_end;test=1
		      ENDIF
		    CASE(12)  ! Superior conjunction time
		      fil2 = 'oct' // let1 // '.res'; graphfile = 'oct'//let1//'_h.sm'
		      xlabel = '     Toc [JD]       '
		      texb =' Tocc =             ';texunit='     JD';bf_val=octime(j,soluce);&
		        & temp1D=octime_end(j,:); test=1  
		    CASE(13)  ! Rp/R*
		      fil2 = 'rr' // let1 // '.res'; graphfile = 'rr'//let1//'_h.sm'
		      xlabel = '       Rp/Rs        '
		      texb =' Rp/Rs =            ';texunit='       ';bf_val=rr(j,soluce);&
		        & temp1D=rr_end(j,:); test=1
		    CASE(14)  ! a/R*
		      fil2 = 'ar' // let1 // '.res'; graphfile = 'ar'//let1//'_h.sm'
		      xlabel = '        a/Rs        '
		      texb =' a/Rs =             ';texunit='       ';bf_val=a_R(j,soluce);&
		        & temp1D=a_R_end(j,:); test=1
		    CASE(15)  ! RV semi-amplitude
		      fil2 = 'k' // let1 // '.res'; graphfile = 'k'//let1//'_h.sm'
		      xlabel = '      K [m/s]       '
		      texb =' K =                ';texunit='    m/s';bf_val=ka(j,soluce);&
		        & temp1D=ka_end(j,:); test=1
		    CASE(16)  ! Excentricity
		      fil2 = 'e' // let1 // '.res'; graphfile = 'e'//let1//'_h.sm'
		      xlabel = '         e          '
		      texb =' e =                ';texunit='       ';bf_val=exc(j,soluce);&
		        & temp1D=exc_end(j,:); test=1
		    CASE(17)  ! Argument of pericenter
		      fil2 = 'omega' // let1 // '.res'; graphfile = 'omega'//let1//'_h.sm'
		      xlabel = '      w [deg]       '
		      texb =' omega =            ';texunit='    deg';bf_val=omega(j,soluce);&
		        & temp1D=omega_end(j,:); test=1
		    CASE(18)  ! ecosw
		      fil2 = 'ecosw' // let1 // '.res'; graphfile = 'ecosw'//let1//'_h.sm'
		      xlabel = '       ecosw        '
		      texb =' ecosw =            ';texunit='       ';bf_val=ecosw(j,soluce);&
		        & temp1D=ecosw_end(j,:); test=1
		    CASE(19)  ! esinw
		      fil2 = 'esinw' // let1 // '.res'; graphfile = 'esinw'//let1//'_h.sm'
		      xlabel = '       esinw        '
		      texb =' esinw =            ';texunit='       ';bf_val=esinw(j,soluce);&
		        & temp1D=esinw_end(j,:); test=1
		    CASE(20)  ! Semi-major axis
		      fil2 = 'a' // let1 // '.res'; graphfile = 'a'//let1//'_h.sm'
		      xlabel = '       a [AU]       '
		      texb =' a =                ';texunit='     AU';bf_val=semi(j,soluce);&
		        & temp1D=semi_end(j,:); test=1
		    CASE(21)  ! Roche limit
		      fil2 = 'roche' // let1 // '.res'; graphfile = 'roche'//let1//'_h.sm'
		      xlabel = '    a_roche [AU]    '
		      texb =' a_roche =          ';texunit='     AU';bf_val=roche(j,soluce);&
		        & temp1D=roche_end(j,:); test=1
		    CASE(22)  ! Ratio semi-major axis on Roche limit
		      fil2 = 'a_on_roche' // let1 // '.res'; graphfile = 'a_on_roche'//let1//'_h.sm'
		      xlabel = '     a/a_roche      '
		      texb =' a/a_roche =        ';texunit='       ';bf_val=a_roche(j,soluce);&
		        & temp1D=a_roche_end(j,:); test=1
		    CASE(23)  ! Orbital inclination
		      fil2 = 'i' // let1 // '.res'; graphfile = 'i'//let1//'_h.sm'
		      xlabel = '      i [deg]       '
		      texb =' i =                ';texunit='    deg';bf_val=inclian(j,soluce);&
		        & temp1D=inclian_end(j,:); test=1
		    CASE(24)  ! Projected obliquity
		      IF(j.EQ.1)THEN
		        fil2 = 'beta' // let1 // '.res'; graphfile = 'beta'//let1//'_h.sm'
		        xlabel = '     Beta [deg]     '
		        texb =' beta =             ';texunit='    deg';bf_val=beta(soluce);&
		          & temp1D=beta_end(:); test=1
		      ENDIF
		    CASE(25)  ! Transit impact parameter
		      fil2 = 'btr' // let1 // '.res'; graphfile = 'btr'//let1//'_h.sm'
		      xlabel = '      btr [R*]      '
		      texb =' b_tr =             ';texunit='     R*';bf_val=b2(j,soluce);&
		        & temp1D=b2_end(j,:); test=1
		    CASE(26)  ! Occultation impact parameter
		      fil2 = 'boc' // let1 // '.res'; graphfile = 'boc'//let1//'_h.sm'
		      xlabel = '      boc [R*]      '
		      texb =' b_oc =             ';texunit='     R*';bf_val=b3(j,soluce);&
		        & temp1D=b3_end(j,:); test=1
		    CASE(27)  ! Prior transit probability
		      fil2 = 'prtr' // let1 // '.res'; graphfile = 'prtr'//let1//'_h.sm'
		      xlabel = '       Ptr [%]      '
		      texb =' P_tr =             ';texunit='      %';bf_val=prtr(j,soluce)*100.;&
		        & temp1D=prtr_end(j,:)*100.; test=1
		    CASE(28)  ! Prior occultation probability
		      fil2 = 'proc' // let1 // '.res'; graphfile = 'proc'//let1//'_h.sm'
		      xlabel = '       Poc [%]      '
		      texb =' P_oc =             ';texunit='      %';bf_val=proc(j,soluce)*100.;&
		        & temp1D=proc_end(j,:)*100.; test=1
		    CASE(29)  ! Planetary density (Jupiter unit)
		      fil2 = 'rhop_j' // let1 // '.res'; graphfile = 'rhop_j'//let1//'_h.sm'
		      xlabel = '     rho [rho_J]    '
		      texb =' rho_p =            ';texunit='rho_jup';bf_val=rhop(j,soluce);&
		        & temp1D=rhop_end(j,:); test=1
		    CASE(30)  ! Planetary density (Earth unit)
		      cons1 = (earthmass/jupmass)/(earthra/jupra)**3
		      fil2 = 'rhop_e' // let1 // '.res'; graphfile = 'rhop_e'//let1//'_h.sm'
		      xlabel = '     rho [rho_E]    '
		      texb ='                    ';texunit='rho_ear';bf_val=rhop(j,soluce)/cons1;&
		        & temp1D=rhop_end(j,:)/cons1; test=1  
		    CASE(31)  ! Planetary density (cgs)
		      fil2 = 'rhop' // let1 // '.res'; graphfile = 'rhop'//let1//'_h.sm'
		      xlabel = '    rho [g.cm-3]    '
		      texb ='                    ';texunit=' g.cm-3';bf_val=rhop(j,soluce)*1.33;&
		        & temp1D=rhop_end(j,:)*1.33; test=1
		    CASE(32)  ! Planetary gravity (dex)
		      fil2 = 'loggp' // let1 // '.res'; graphfile = 'loggp'//let1//'_h.sm'
		      xlabel = '      log(g)_p      '
		      texb =' logg_p =           ';texunit='       ';bf_val=logg_p(j,soluce);&
		        & temp1D=logg_p_end(j,:); test=1 
		    CASE(33)  ! Equilibrium temperature 
		      fil2 = 'teqp' // let1 // '.res'; graphfile = 'teqp'//let1//'_h.sm'
		      xlabel = '       Teq [K]      '
		      texb =' Teq_p =            ';texunit='      K';bf_val=teq_p(j,soluce);&
		        & temp1D=teq_p_end(j,:); test=1
		    CASE(34)  ! Hill radius
		      fil2 = 'hill' // let1 // '.res'; graphfile = 'hill'//let1//'_h.sm'
		      xlabel = '      HilR [Rp]     '
		      texb =' Hill radius =         ';texunit='     Rp';bf_val=hillrad_p(j,soluce);&
		        & temp1D=hillrad_p_end(j,:); test=1
		    CASE(35)  ! Irradiation 
		      fil2 = 'irra' // let1 // '.res'; graphfile = 'irra'//let1//'_h.sm'
		      xlabel = ' Irradation [Earth] '
		      texb =' Irradation =       ';texunit='  Earth';bf_val=irrad(j,soluce);&
		        & temp1D=irrad_end(j,:); test=1
		    CASE(36)  ! Planet mass (Jupiter unit) 
		      fil2 = 'mp' // let1 // '.res'; graphfile = 'mp'//let1//'_h.sm'
		      xlabel = '      M_p [M_J]     '
		      texb =' Mp =               ';texunit='  M_jup';bf_val=mass_p(j,soluce);&
		        & temp1D=mass_p_end(j,:); test=1
		    CASE(37)  ! Planet mass (Earth unit)
		      cons1 = (earthmass/jupmass)
		      fil2 = 'mpe' // let1 // '.res'; graphfile = 'mp_e'//let1//'_h.sm'
		      xlabel = '      M_p [M_E]     '
		      texb ='                    ';texunit='  M_ear';bf_val=mass_p(j,soluce)/cons1;&
		        & temp1D=mass_p_end(j,:)/cons1; test=1
		    CASE(38)  ! Planet msini (Jupiter unit) 
		      fil2 = 'mpsini' // let1 // '.res'; graphfile = 'mpsini'//let1//'_h.sm'
		      xlabel = '    Msini [M_J]     '
		      texb =' Mp sini =          ';texunit='  M_jup';bf_val=mass_p_sini(j,soluce);&
		        & temp1D=mass_p_sini_end(j,:); test=1
		    CASE(39)  ! Planet msini (Earth unit)
		      cons1 = (earthmass/jupmass)
		      fil2 = 'mpsinie' // let1 // '.res'; graphfile = 'mpsini_e'//let1//'_h.sm'
		      xlabel = '    Msini [M_E]     '
		      texb ='                    ';texunit='  M_ear';bf_val=mass_p_sini(j,soluce)/cons1;&
		        & temp1D=mass_p_sini_end(j,:)/cons1; test=1
		    CASE(40)  ! Planet radius (Jupiter unit) 
		      fil2 = 'rp' // let1 // '.res'; graphfile = 'rp'//let1//'_h.sm'
		      xlabel = '      R_p [R_J]     '
		      texb =' Rp =               ';texunit='  R_jup';bf_val=radius_p(j,soluce);&
		        & temp1D=radius_p_end(j,:); test=1
		    CASE(41)  ! Planet radius (Earth unit)
		      cons1 = (earthra/jupra)
		      fil2 = 'rpe' // let1 // '.res'; graphfile = 'rp_e'//let1//'_h.sm'
		      xlabel = '      R_p [R_E]     '
		      texb ='                    ';texunit='  R_ear';bf_val=radius_p(j,soluce)/cons1;&
		        & temp1D=radius_p_end(j,:)/cons1; test=1
		    CASE(42)  ! Planet Safronov number
		      fil2 = 'safro' // let1 // '.res'; graphfile = 'safro'//let1//'_h.sm'
		      xlabel = '  Safronov number   '
		      texb =' Safronov number =  ';texunit='       ';bf_val=safro(j,soluce);&
		        & temp1D=safro_end(j,:); test=1
		  END SELECT
		  IF(test.EQ.1)THEN
		    OPEN(UNIT=124,FILE=fil2)
		    DO l=1,mco
		      WRITE(124,*) l,temp1D(l)
		    ENDDO
		    CLOSE(124)
		    CALL SORT(mco,temp1D)
		    CALL PDLIM(mco,temp1D,bf_val,limmed,limpd)
		    medi_val = temp1d(limmed(5))
		    CALL HISTO(fil2,graphfile,bf_val,medi_val,xlabel,nbinh)
		    CALL PDFSCREEN(texb,texunit,bf_val,medi_val,limpd,gelman,gelmanval, & 
		     & nchain,jump)
		  ENDIF
		ENDDO    
		IF(nddf.GT.0.AND.ntr.GT.0.AND.isddf.NE.'n')THEN  ! Planet size spectrum
		  DO i=1,nddf
		    DO k=1,3
		      SELECT CASE(k)
		        CASE(1)
		          texb =' '//ddf_filter(i)//'-R_p =            ';texunit='  R_jup';&
		           & bf_val=radipla(j,i,soluce);temp1D=radipla_end(j,i,:)
		          fil2 = 'rp_' // let1 // '_' // ddf_filter(i) // '.res'
		          graphfile = 'rp_'// let1 // '_' // ddf_filter(i)//'_h.sm'
		          xlabel = '      R_p [R_J]     '
		        CASE(2)
		          texb =' '//ddf_filter(i)//'-Rp/Rs =          ';texunit='       ';&
		           & bf_val=dratio(j,i,soluce);temp1D=dratio_end(j,i,:)
		          fil2 = 'rr_' // let1 // '_' // ddf_filter(i) // '.res'
		          graphfile = 'rr_'// let1 // '_' // ddf_filter(i)//'_h.sm'
		          xlabel = '       R_p/R_s      '
		        CASE(3)
		          texb =' '//ddf_filter(i)//'-dF =             ';texunit='       ';&
		           & bf_val=ddepth(j,i,soluce);temp1D=ddepth_end(j,i,:)
		          fil2 = 'dF_' // let1 // '_' // ddf_filter(i) // '.res'
		          graphfile = 'dF_'// let1 // '_' // ddf_filter(i)//'_h.sm'
		          xlabel = '         dF         '
		      END SELECT
		      OPEN(UNIT=124,FILE=fil2)
		      DO l=1,mco
		        WRITE(124,*) l,temp1D(l)
		      ENDDO
		      CLOSE(124)
		      CALL SORT(mco,temp1D)
		      CALL PDLIM(mco,temp1D,bf_val,limmed,limpd)
		      medi_val = temp1d(limmed(5))
		      CALL HISTO(fil2,graphfile,bf_val,medi_val,xlabel,nbinh)
		      CALL PDFSCREEN(texb,texunit,bf_val,medi_val,limpd,gelman,gelmanval, & 
		       & nchain,jump)
		    ENDDO
		  ENDDO
		ENDIF
		IF(nttvmax.GT.0.AND.isttv.EQ.'y')THEN   !Transit timings
		  let1 = CHAR(48+INT(j/10)); let2 = CHAR(48+MOD(j,10))
		  DO i=1,nttv(j)
		    IF(i.LT.10)THEN
		      let3=CHAR(48+i); texb =' Ttr-'//let3//' = '
		      fil2 = 'ttr' // let1 // let2//'_'//let3// '.res'
		      graphfile = 'ttr'//let1//let2//'_'//let3//'_h.sm'
		    ELSE
		      let3 = CHAR(48+INT(i/10)); let4 = CHAR(48+MOD(i,10))
		      texb =' Ttr-'//let3//let4// ' = '
		      fil2 = 'ttr' // let1 // let2 // '_' //let3 // let4 // '.res'
		      graphfile = 'ttr'//let1//let2//'_'//let3//let4//'_h.sm'
		   ENDIF
		   texunit='     JD'; bf_val=ttr(j,i,soluce); temp1D=ttr_end(j,i,:)
		    xlabel = '      T_tr [JD]     '
		    OPEN(UNIT=124,FILE=fil2)
		    DO l=1,mco
		      WRITE(124,*) l,temp1D(l)
		    ENDDO
		    CLOSE(124)
		    CALL SORT(mco,temp1D)
		    CALL PDLIM(mco,temp1D,bf_val,limmed,limpd)
		    medi_val = temp1d(limmed(5))
		    CALL HISTO(fil2,graphfile,bf_val,medi_val,xlabel,nbinh)
		    CALL PDFSCREEN(texb,texunit,bf_val,medi_val,limpd,gelman,gelmanval, & 
		     & nchain,jump)
		  ENDDO
		  WRITE(*,129)   '  rms TTV =           ',ttvrms(j),'   min' 
		  WRITE(444,129) '  rms TTV  =          ',ttvrms(j),'   min' 
		  WRITE(445,129) '  rms TTV  =          ',ttvrms(j),'   min' 
		ENDIF

	  ENDDO
  end if

  PRINT*,      '------------------------------------------------------------------------------------------'
  WRITE(444,*) '------------------------------------------------------------------------------------------'
  WRITE(445,*) '------------------------------------------------------------------------------------------'

  DO k=1,4
     test=0
     SELECT CASE (k)
        CASE(1) ! Reduced Chi2
           fil2 = 'red_merit.res'; graphfile = 'red_merit_h.sm'
           xlabel = '        Rchi2       '
           texb =' RChi2 =            ';texunit='       ';bf_val=remerit(soluce);&
            & temp1D=remerit_end; test=1
        CASE(2) ! Chi2
           fil2 = 'merit.res'; graphfile = 'merit_h.sm'
           xlabel = '        Chi2        '
           texb =' Chi2 =             ';texunit='       ';bf_val=merit(soluce);&
            & temp1D=merit_end; test=1
        CASE(3) ! Phot Chi2
           IF(ntr.GT.0)THEN
              fil2 = 'photmerit.res'; graphfile = 'photmerit_h.sm'
              xlabel = '     Phot_Chi2      '
              texb =' Phot_Chi2 =        ';texunit='       ';bf_val=photmerit(soluce);&
               & temp1D=photmerit_end; test=1
           ENDIF
        CASE(4) ! RV Chi2
           IF(nrv.GT.0)THEN
              fil2 = 'rvmerit.res'; graphfile = 'rvmerit_h.sm'
              xlabel = '      RV_Chi2       '
              texb =' RV_Chi2 =          ';texunit='       ';bf_val=rvmerit(soluce);&
               & temp1D=rvmerit_end; test=1
           ENDIF
     END SELECT
     IF(test.EQ.1)THEN
       OPEN(UNIT=124,FILE=fil2)
       DO l=1,mco
          WRITE(124,*) l,temp1D(l)
       ENDDO
       CLOSE(124)    
       CALL SORT(mco,temp1D)
       CALL PDLIM(mco,temp1D,bf_val,limmed,limpd)
       medi_val = temp1d(limmed(5))
       CALL HISTO(fil2,graphfile,bf_val,medi_val,xlabel,nbinh)
       CALL PDFSCREEN(texb,texunit,bf_val,medi_val,limpd,gelman,gelmanval, & 
        & nchain,jump)
    ENDIF
  ENDDO

  IF(ntr.GT.0)THEN
    WRITE(*,116)   ' BF Beta_red =      ',beta_red_glo(soluce)
    WRITE(444,116) ' BF Beta_red =      ',beta_red_glo(soluce)
    WRITE(445,116) ' BF Beta_red =      ',beta_red_glo(soluce)
  ENDIF
  IF(nrv.GT.0)THEN
    WRITE(*,129)   ' BF jitter =        ',jitterglo(soluce),'    m/s' 
    WRITE(444,129) ' BF jitter =        ',jitterglo(soluce),'    m/s' 
    WRITE(445,129) ' BF jitter =        ',jitterglo(soluce),'    m/s' 
  ENDIF
  WRITE(*,116)     ' Prob. TR =         ', protransit
  WRITE(*,116)     ' Prob. full TR =    ', profulltransit
  WRITE(*,116)     ' Prob. OC =         ', prosecond
  WRITE(*,116)     ' Prob. full OC =    ', profullsecond
  WRITE(444,116)   ' Prob. TR =         ', protransit
  WRITE(444,116)   ' Prob. full TR =    ', profulltransit
  WRITE(444,116)   ' Prob. OC =         ', prosecond
  WRITE(444,116)   ' Prob. full OC =    ', profullsecond
  WRITE(445,116)   ' Prob. TR =         ', protransit
  WRITE(445,116)   ' Prob. full TR =    ', profulltransit
  WRITE(445,116)   ' Prob. OC =         ', prosecond
  WRITE(445,116)   ' Prob. full OC =    ', profullsecond
  WRITE(*,116)     ' AIC =              ', aic
  WRITE(*,116)     ' BIC =              ', bic
  WRITE(*,116)     ' DIC =              ', dic+dic-dic2
  WRITE(444,116)   ' AIC =              ', aic
  WRITE(444,116)   ' BIC =              ', bic
  WRITE(444,116)   ' DIC =              ', dic+dic-dic2
  WRITE(445,116)   ' AIC =              ', aic
  WRITE(445,116)   ' BIC =              ', bic
  WRITE(445,116)   ' DIC =              ', dic+dic-dic2
  
  PRINT*,      '------------------------------------------------------------------------------------------'
  WRITE(444,*) '------------------------------------------------------------------------------------------'
  WRITE(445,*) '------------------------------------------------------------------------------------------'

  CLOSE(444)
  CLOSE(445)

END PROGRAM MCMCI

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE occultnl(rl,c1,c2,c3,c4,b0,mulimb0,mulimbf,nb)
! Please cite Mandel & Agol (2002) if making use of this routine.
  IMPLICIT NONE

  INTEGER, PARAMETER :: nmax=523

  INTEGER i,j,nb,nr,i1,i2
  DOUBLE PRECISION :: mulimbfk(nb,5),pi,c1,c2,c3,c4,rl,bt0(nb),b0k(nb),b0(1), &
     &       mulimb0k(nb),mulimb(nb),mulimbp(nb),dt,t(nmax),th(nmax),r(nmax),&
     &       mulimbf(5),mulimb0(1),&
     &       sig,mulimb1(nb),mulimbhalf(nb),mulimb3half(nb),mulimb2(nb),&
     &       sig1,sig2,omega,dmumax,fac,mu(nb),f1,f2
     PI=DACOS(-1.d0)
!  This routine uses the results for a uniform source to
!  compute the lightcurve for a limb-darkened source
!  (5-1-02 notes)
! Input:
!   rl        radius of the lens   in units of the source radius
!   c1-c4     limb-darkening coefficients
!   b0        impact parameter normalized to source radius
! Output:
!  mulimb0 limb-darkened magnification
!  mulimbf lightcurves for each component
!  
!  First, make grid in radius:
!  Call magnification of uniform source:

!  Mettre b0 de cote, le garder pour la fin 
     b0k(1)=0.1d0
     b0k(2)=b0(1)

   CALL occultuniform(b0k,rl,mulimb0k,nb)
   i1=nb
   i2=1
   fac=0.d0
   DO i=1,nb
     bt0(i)=b0k(i)
     mulimbfk(i,1)=1.d0
     mulimbfk(i,2)=0.8d0
     mulimbfk(i,3)=2.d0/3.d0
     mulimbfk(i,4)=4.d0/7.d0
     mulimbfk(i,5)=0.5d0
     mulimb(i)=mulimb0k(i)
     IF(mulimb0k(i).ne.1.d0)THEN
       i1=MIN(i1,i)
       i2=MAX(i2,i)
     ENDIF
     fac=MAX(fac,ABS(mulimb0k(i)-1.d0))
   ENDDO
!   PRINT*, rl
   omega=4.*((1.d0-c1-c2-c3-c4)/4.+c1/5.+c2/6.+c3/7.+c4/8.)
   nr=2
   dmumax=1.d0
!   WRITE(*,*) 'i1,i2 ',i1,i2
   DO WHILE(dmumax.GT.fac*1.d-3)
     DO i=i1,i2
       mulimbp(i)=mulimb(i)
     ENDDO
     nr=nr*2
!     WRITE(*,*) 'nr ',nr
     dt=0.5d0*pi/DBLE(nr)
     DO j=1,nr+1
       t(j) =dt*DBLE(j-1)
       th(j)=t(j)+0.5d0*dt
       r(j)=DSIN(t(j))
     ENDDO
     sig=SQRT(DCOS(th(nr)))
     DO i=i1,i2
       mulimbhalf(i) =sig**3*mulimb0k(i)/(1.d0-r(nr))
       mulimb1(i)    =sig**4*mulimb0k(i)/(1.d0-r(nr))
       mulimb3half(i)=sig**5*mulimb0k(i)/(1.d0-r(nr))
       mulimb2(i)    =sig**6*mulimb0k(i)/(1.d0-r(nr))
     ENDDO
     DO j=2,nr
       DO i=1,nb
         b0k(i)=bt0(i)/r(j)
       ENDDO
!  Calculate uniform magnification at intermediate radii:
       CALL occultuniform(b0k,rl/r(j),mu,nb)
!  Equation (29):
       sig1=SQRT(DCOS(th(j-1)))
       sig2=SQRT(DCOS(th(j)))
       dmumax=0.d0
       DO i=i1,i2
         f1=r(j)*r(j)*mu(i)/(r(j)-r(j-1))
         f2=r(j)*r(j)*mu(i)/(r(j+1)-r(j))
         mulimbhalf(i) =mulimbhalf(i) +f1*sig1**3-f2*sig2**3
         mulimb1(i)    =mulimb1(i)    +f1*sig1**4-f2*sig2**4
         mulimb3half(i)=mulimb3half(i)+f1*sig1**5-f2*sig2**5
         mulimb2(i)    =mulimb2(i)    +f1*sig1**6-f2*sig2**6
         mulimb(i)=((1.d0-c1-c2-c3-c4)*mulimb0k(i)+c1*mulimbhalf(i)*dt &
     &        +c2*mulimb1(i)*dt+c3*mulimb3half(i)*dt+c4*mulimb2(i)*dt) &
     &        /omega
         IF(ABS(mulimb(i)+mulimbp(i)).GT.1.E-13)THEN 
           dmumax=max(dmumax,abs(mulimb(i)-mulimbp(i))/(mulimb(i)+ &
     &               mulimbp(i)))
         ENDIF
       ENDDO
     ENDDO
   ENDDO
   DO i=i1,i2
     mulimbfk(i,1)=mulimb0k(i)
     mulimbfk(i,2)=mulimbhalf(i)*dt
     mulimbfk(i,3)=mulimb1(i)*dt
     mulimbfk(i,4)=mulimb3half(i)*dt
     mulimbfk(i,5)=mulimb2(i)*dt
     mulimb0k(i)=mulimb(i)
   ENDDO
   DO i=1,nb
     b0k(i)=bt0(i)
   ENDDO
   mulimb0(1)=mulimb0k(2)
   b0(1)     =b0k(2)
   DO i=1,5
      mulimbf(i)=mulimbfk(2,i)
   ENDDO
   RETURN
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE occultnlold(rl,c1,c2,c3,c4,b0,mulimb0,mulimbf,nb)
! Please cite Mandel & Agol (2002) if making use of this routine.
  IMPLICIT NONE

  INTEGER, PARAMETER :: nmax=523

  INTEGER i,j,nb,nr,i1,i2
  DOUBLE PRECISION :: mulimbf(nb,5),pi,c1,c2,c3,c4,rl,bt0(nb),b0(nb), &
     &       mulimb0(nb),mulimb(nb),mulimbp(nb),dt,t(nmax),th(nmax),r(nmax),&
     &       sig,mulimb1(nb),mulimbhalf(nb),mulimb3half(nb),mulimb2(nb),&
     &       sig1,sig2,omega,dmumax,fac,mu(nb),f1,f2
     PI=DACOS(-1.d0)
!  This routine uses the results for a uniform source to
!  compute the lightcurve for a limb-darkened source
!  (5-1-02 notes)
! Input:
!   rl        radius of the lens   in units of the source radius
!   c1-c4     limb-darkening coefficients
!   b0        impact parameter normalized to source radius
! Output:
!  mulimb0 limb-darkened magnification
!  mulimbf lightcurves for each component
!  
!  First, make grid in radius:
!  Call magnification of uniform source:
   CALL occultuniform(b0,rl,mulimb0,nb)
   i1=nb
   i2=1
   fac=0.d0
   DO i=1,nb
     bt0(i)=b0(i)
     mulimbf(i,1)=1.d0
     mulimbf(i,2)=0.8d0
     mulimbf(i,3)=2.d0/3.d0
     mulimbf(i,4)=4.d0/7.d0
     mulimbf(i,5)=0.5d0
     mulimb(i)=mulimb0(i)
     IF(mulimb0(i).ne.1.d0)THEN
       i1=MIN(i1,i)
       i2=MAX(i2,i)
     ENDIF
     fac=MAX(fac,ABS(mulimb0(i)-1.d0))
   ENDDO
!   PRINT*, rl
   omega=4.*((1.d0-c1-c2-c3-c4)/4.+c1/5.+c2/6.+c3/7.+c4/8.)
   nr=2
   dmumax=1.d0
!   WRITE(*,*) 'i1,i2 ',i1,i2
   DO WHILE(dmumax.GT.fac*1.d-3)
     print*, 't65'
     DO i=i1,i2
       mulimbp(i)=mulimb(i)
     ENDDO
     nr=nr*2
!     WRITE(*,*) 'nr ',nr
     dt=0.5d0*pi/DBLE(nr)
     DO j=1,nr+1
       t(j) =dt*DBLE(j-1)
       th(j)=t(j)+0.5d0*dt
       r(j)=DSIN(t(j))
     ENDDO
     sig=SQRT(DCOS(th(nr)))
     DO i=i1,i2
       mulimbhalf(i) =sig**3*mulimb0(i)/(1.d0-r(nr))
       mulimb1(i)    =sig**4*mulimb0(i)/(1.d0-r(nr))
       mulimb3half(i)=sig**5*mulimb0(i)/(1.d0-r(nr))
       mulimb2(i)    =sig**6*mulimb0(i)/(1.d0-r(nr))
     ENDDO
     DO j=2,nr
       DO i=1,nb
         b0(i)=bt0(i)/r(j)
       ENDDO
!  Calculate uniform magnification at intermediate radii:
       CALL occultuniform(b0,rl/r(j),mu,nb)
!  Equation (29):
       sig1=SQRT(DCOS(th(j-1)))
       sig2=SQRT(DCOS(th(j)))
       dmumax=0.d0
       DO i=i1,i2
         f1=r(j)*r(j)*mu(i)/(r(j)-r(j-1))
         f2=r(j)*r(j)*mu(i)/(r(j+1)-r(j))
         mulimbhalf(i) =mulimbhalf(i) +f1*sig1**3-f2*sig2**3
         mulimb1(i)    =mulimb1(i)    +f1*sig1**4-f2*sig2**4
         mulimb3half(i)=mulimb3half(i)+f1*sig1**5-f2*sig2**5
         mulimb2(i)    =mulimb2(i)    +f1*sig1**6-f2*sig2**6
         mulimb(i)=((1.d0-c1-c2-c3-c4)*mulimb0(i)+c1*mulimbhalf(i)*dt &
     &        +c2*mulimb1(i)*dt+c3*mulimb3half(i)*dt+c4*mulimb2(i)*dt) &
     &        /omega
         IF(ABS(mulimb(i)+mulimbp(i)).GT.1.E-13)THEN 
           dmumax=max(dmumax,abs(mulimb(i)-mulimbp(i))/(mulimb(i)+ &
     &               mulimbp(i)))
         ENDIF
       ENDDO
     ENDDO
   ENDDO
   DO i=i1,i2
     mulimbf(i,1)=mulimb0(i)
     mulimbf(i,2)=mulimbhalf(i)*dt
     mulimbf(i,3)=mulimb1(i)*dt
     mulimbf(i,4)=mulimb3half(i)*dt
     mulimbf(i,5)=mulimb2(i)*dt
     mulimb0(i)=mulimb(i)
   ENDDO
   DO i=1,nb
     b0(i)=bt0(i)
   ENDDO
   RETURN
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      
SUBROUTINE occultuniform(b0,w,muo1,nb)
  
  IMPLICIT NONE

  INTEGER :: i,nb
  DOUBLE PRECISION :: muo1(nb),w,b0(nb),z,pi,lambdae,kap0,kap1

  IF(abs(w-0.5d0).lt.1.d-3) w=0.5d0
  pi=DACOS(-1.d0)
!  This routine computes the lightcurve for occultation
!  of a uniform source without microlensing  (Mandel & Agol 2002).
! Input:
! 
!  rs   radius of the source (set to unity)
!  b0   impact parameter in units of rs
!  w    occulting star size in units of rs
! 
! Output:
!  muo1 fraction of flux at each b0 for a uniform source
! 
!  Now, compute pure occultation curve:
   DO i=1,nb
!  substitute z=b0(i) to shorten expressions
     z=b0(i)
!  the source is unocculted:
!  Table 3, I.
     IF(z.ge.1.d0+w)THEN
       muo1(i)=1.d0
       GOTO 1
     ENDIF
!  the  source is completely occulted:
!  Table 3, II.
     IF(w.ge.1.d0.and.z.le.w-1.d0)THEN
       muo1(i)=0.d0
       GOTO 1
    ENDIF
!  the source is partly occulted and the occulting object crosses the limb:
!  Equation (26):
     IF(z.GE.abs(1.d0-w).AND.z.LE.1.d0+w)THEN
       kap1=DACOS(MIN((1.d0-w*w+z*z)/2.d0/z,1.d0))
       kap0=DACOS(MIN((w*w+z*z-1.d0)/2.d0/w/z,1.d0))
       lambdae=w*w*kap0+kap1
       lambdae=(lambdae-0.5d0*SQRT(MAX(4.d0*z*z-(1.d0+z*z-w*w)**2, &
     &            0.d0)))/pi
       muo1(i)=1.d0-lambdae
     ENDIF
!  the occulting object transits the source star (but doesn't
!  completely cover it):
     IF(z.LE.1.d0-w) muo1(i)=1.d0-w*w
 1   CONTINUE
   ENDDO
! muo1=1.d0-lambdae
   RETURN
END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE occultquad(z0,u1,u2,p,muo1,mu0,nz)
!  This routine computes the lightcurve for occultation
!  of a quadratically limb-darkened source without microlensing.
!  Please cite Mandel & Agol (2002) if you make use of this routine
!  in your research.  Please report errors or bugs to agol@tapir.caltech.edu
  
  IMPLICIT NONE

  INTEGER :: i,nz

  DOUBLE PRECISION ::  z0(nz),u1,u2,p,muo1(nz),mu0(nz)
  DOUBLE PRECISION ::  lambdad(nz),etad(nz),lambdae(nz),lam
  DOUBLE PRECISION ::  pi,x1,x2,x3,z,omega,kap0,kap1,q,Kk,Ek,Pk,n
  DOUBLE PRECISION :: ellec,ellk,rj

  IF(ABS(p-0.5d0).LT.1.d-3) p=0.5d0

! Input:
!
! rs   radius of the source (set to unity)
! z0   impact parameter in units of rs
! p    occulting star size in units of rs
! u1   linear    limb-darkening coefficient (gamma_1 in paper)
! u2   quadratic limb-darkening coefficient (gamma_2 in paper)
!
! Output:
!
! muo1 fraction of flux at each z0 for a limb-darkened source
! mu0  fraction of flux at each z0 for a uniform source
!
! Limb darkening has the form:
!  I(r)=[1-u1*(1-sqrt(1-(r/rs)^2))-u2*(1-sqrt(1-(r/rs)^2))^2]/(1-u1/3-u2/6)/pi
! 
! To use this routine
!
! Now, compute pure occultation curve:
  omega=1.d0-u1/3.d0-u2/6.d0
  pi=DACOS(-1.d0)
! Loop over each impact parameter:
  DO i=1,nz
! substitute z=z0(i) to shorten expressions
    z=z0(i)
    x1=(p-z)**2
    x2=(p+z)**2
    x3=p**2-z**2
! the source is unocculted:
! Table 3, I.
    IF(z.GE.1.d0+p)THEN
      lambdad(i)=0.d0
      etad(i)=0.d0
      lambdae(i)=0.d0
      GOTO 10
    ENDIF
! the  source is completely occulted:
! Table 3, II.
    IF(p.GE.1.d0.AND.z.LE.p-1.d0)THEN
      lambdad(i)=1.d0
      etad(i)=1.d0
      lambdae(i)=1.d0
      GOTO 10
    ENDIF
! the source is partly occulted and the occulting object crosses the limb:
! Equation (26):
    IF(z.GE.ABS(1.d0-p).AND.z.LE.1.d0+p)THEN
      kap1=DACOS(MIN((1.d0-p*p+z*z)/2.d0/z,1.d0))
      kap0=DACOS(MIN((p*p+z*z-1.d0)/2.d0/p/z,1.d0))
      lambdae(i)=p*p*kap0+kap1
      lambdae(i)=(lambdae(i)-0.5d0*SQRT(MAX(4.d0*z*z- &
     &               (1.d0+z*z-p*p)**2,0.d0)))/pi
    ENDIF
! the occulting object transits the source star (but doesn't
! completely cover it):
    IF(z.LE.1.d0-p) lambdae(i)=p*p
! the edge of the occulting star lies at the origin- special 
! expressions in this case:
    IF(ABS(z-p).LT.1.d-4*(z+p))THEN
! Table 3, Case V.:
      IF(z.GE.0.5d0)THEN
        lam=0.5d0*pi
        q=0.5d0/p
        Kk=ellk(q)
        Ek=ellec(q)
! Equation 34: lambda_3
        lambdad(i)=1.d0/3.d0+16.d0*p/9.d0/pi*(2.d0*p*p-1.d0)*Ek- &
     &                 (32.d0*p**4-20.d0*p*p+3.d0)/9.d0/pi/p*Kk
! Equation 34: eta_1
        etad(i)=1.d0/2.d0/pi*(kap1+p*p*(p*p+2.d0*z*z)*kap0- &
     &              (1.d0+5.d0*p*p+z*z)/4.d0*sqrt((1.d0-x1)*(x2-1.d0)))
        IF(p.EQ.0.5d0)THEN
! Case VIII: p=1/2, z=1/2
          lambdad(i)=1.d0/3.d0-4.d0/pi/9.d0
          etad(i)=3.d0/32.d0
        ENDIF
        GOTO 10
      ELSE
! Table 3, Case VI.:
        lam=0.5d0*pi
        q=2.d0*p
        Kk=ellk(q)
        Ek=ellec(q)
! Equation 34: lambda_4
        lambdad(i)=1.d0/3.d0+2.d0/9.d0/pi*(4.d0*(2.d0*p*p-1.d0)*Ek+ &
     &                 (1.d0-4.d0*p*p)*Kk)
! Equation 34: eta_2
        etad(i)=p*p/2.d0*(p*p+2.d0*z*z)
        GOTO 10
      ENDIF
    ENDIF
! the occulting star partly occults the source and crosses the limb:
! Table 3, Case III:
    IF((z.GT.0.5d0+ABS(p-0.5d0).AND.z.LT.1.d0+p).OR.(p.GT.0.5d0 &
     &      .AND.z.GT.ABS(1.d0-p)*1.0001d0.AND.z.LT.p))THEN
      lam=0.5d0*pi
      q=SQRT((1.d0-(p-z)**2)/4.d0/z/p)
      Kk=ellk(q)
      Ek=ellec(q)
      n=1.d0/x1-1.d0
      Pk=Kk-n/3.d0*rj(0.d0,1.d0-q*q,1.d0,1.d0+n)
! Equation 34, lambda_1:
      lambdad(i)=1.d0/9.d0/pi/sqrt(p*z)*(((1.d0-x2)*(2.d0*x2+ &
     &        x1-3.d0)-3.d0*x3*(x2-2.d0))*Kk+4.d0*p*z*(z*z+ &
     &        7.d0*p*p-4.d0)*Ek-3.d0*x3/x1*Pk)
      IF(z.LT.p) lambdad(i)=lambdad(i)+2.d0/3.d0
! Equation 34, eta_1:
      etad(i)=1.d0/2.d0/pi*(kap1+p*p*(p*p+2.d0*z*z)*kap0- &
     &          (1.d0+5.d0*p*p+z*z)/4.d0*SQRT((1.d0-x1)*(x2-1.d0)))
      GOTO 10
    ENDIF
! the occulting star transits the source:
! Table 3, Case IV.:
    IF(p.LE.1.d0.AND.z.LE.(1.d0-p)*1.0001d0)THEN
      lam=0.5d0*pi
      q=SQRT((x2-x1)/(1.d0-x1))
      Kk=ellk(q)
      Ek=ellec(q)
      n=x2/x1-1.d0
      Pk=Kk-n/3.d0*rj(0.d0,1.d0-q*q,1.d0,1.d0+n)
! Equation 34, lambda_2:
      lambdad(i)=2.d0/9.d0/pi/SQRT(1.d0-x1)*((1.d0-5.d0*z*z+p*p+ &
     &         x3*x3)*Kk+(1.d0-x1)*(z*z+7.d0*p*p-4.d0)*Ek-3.d0*x3/x1*Pk)
      IF(z.LT.p) lambdad(i)=lambdad(i)+2.d0/3.d0
      IF(ABS(p+z-1.d0).LE.1.d-4)THEN
        lambdad(i)=2/3.d0/pi*DACOS(1.d0-2.d0*p)-4.d0/9.d0/pi* &
     &            SQRT(p*(1.d0-p))*(3.d0+2.d0*p-8.d0*p*p)
      ENDIF
! Equation 34, eta_2:
          etad(i)=p*p/2.d0*(p*p+2.d0*z*z)
    ENDIF
 10 CONTINUE
! Now, using equation (33):
    muo1(i)=1.d0-((1.d0-u1-2.d0*u2)*lambdae(i)+(u1+2.d0*u2)* &
     &      lambdad(i)+u2*etad(i))/omega
! Equation 25:
    mu0(i)=1.d0-lambdae(i)
  ENDDO
  RETURN
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE mr(ms,rs)
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: k1=-0.20945,k2=0.0804,k3=0.304
  DOUBLE PRECISION, PARAMETER :: rearth=6378.,rsun=6.96342e5
  DOUBLE PRECISION, PARAMETER :: m1_fe=4.34,r1_fe=2.23
  DOUBLE PRECISION, PARAMETER :: m1_si=7.38,r1_si=3.58
  DOUBLE PRECISION, PARAMETER :: m1_vo=8.16,r1_vo=4.73
  DOUBLE PRECISION :: ms,rs,m1,r1
  DOUBLE PRECISION :: fra_fe,fra_si,fra_vo
  fra_fe=1.;fra_vo=0.
  fra_si=1.-fra_vo-fra_fe
  m1=fra_fe*m1_fe + fra_vo*m1_vo + fra_si*m1_si
  r1=fra_fe*r1_fe + fra_vo*r1_vo + fra_si*r1_si
  rs = k1+(1./3.)*LOG10(ms/m1)-k2*(ms/m1)**k3
  rs = 10**rs
  rs = rs*r1
  RETURN
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE qdinter(spectro,espectro,filter,a,ea)

  IMPLICIT NONE
  
  INTEGER :: i,k,test2
  INTEGER, PARAMETER :: np=147300 !135516

  DOUBLE PRECISION, DIMENSION(19), PARAMETER :: mec = (/ -5.0,-4.5,-4.0, &
  & -3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,-0.3,-0.2,-0.1,0.,0.1,0.2,  &
  &   0.3,0.5,1.0 /)
  DOUBLE PRECISION, DIMENSION(5), PARAMETER :: vit = (/0., 1., 2., 4., 8./)
  DOUBLE PRECISION, DIMENSION(27) :: hot = (/ 3500,3750,4000,4250,4500,4750, &
  & 5000,5250,5500,5750,6000,6250,6500,6750,7000,7250,7500,7750,8000, &
  & 8250,8500,8750,9000,9250,9500,9750,10000 /)
  DOUBLE PRECISION, DIMENSION(11) :: grav =  (/ 0.0,0.5,1.0,1.5,2.0, &
  & 2.5,3.0,3.5,4.0,4.5,5.0 /)
  DOUBLE PRECISION, DIMENSION(4) ::  spectro,espectro
  DOUBLE PRECISION, DIMENSION(2,4) :: der
  DOUBLE PRECISION, DIMENSION(2) ::  a,ea
  DOUBLE PRECISION ::  dif
  DOUBLE PRECISION ::  teff,steff,teff_2
  DOUBLE PRECISION ::  logg,slogg,logg_2
  DOUBLE PRECISION ::  metal,smetal,metal_2
  DOUBLE PRECISION ::  vturb,svturb,vturb_2
  DOUBLE PRECISION ::  errtol = 0.01
  DOUBLE PRECISION :: rturb,rteff,rlogg,rmet,ru1,ru2
  DOUBLE PRECISION, DIMENSION(2,6) :: coeff

  CHARACTER(LEN=2) :: filter,rfilt

  dif=1000.
  k=0
  DO i=1,27
    IF(ABS(hot(i)-spectro(1)).LT.dif)THEN
      teff = hot(i)
      k=i
      steff = spectro(1)-teff
      dif = ABS(steff)
    ENDIF
  ENDDO
  IF(k.EQ.1)THEN
    teff_2 = hot(2)
  ELSE IF(k.EQ.27)THEN
    teff_2 = hot(26)
  ELSE IF(steff.GT.0.)THEN
    teff_2 = hot(k+1)
  ELSE
    teff_2 = hot(k-1)
  ENDIF
  
  dif = 1000.
  DO i=1,11
    IF(ABS(grav(i)-spectro(2)).LT.dif)THEN
      logg = grav(i)
      k=i
      slogg = spectro(2)-logg
      dif = ABS(slogg)
    ENDIF
  ENDDO
  IF(k.EQ.1)THEN
    logg_2 = grav(2)
  ELSE IF(k.EQ.11)THEN
    logg_2 = grav(10)
  ELSE IF(steff.GT.0.)THEN
    logg_2 = grav(k+1)
  ELSE
    logg_2 = grav(k-1)
  ENDIF
  dif = 1000.
  DO i=1,19
    IF(ABS(mec(i)-spectro(3)).LT.dif)THEN
      metal = mec(i)
      k=i
      smetal = spectro(3)-metal
      dif = ABS(smetal)
    ENDIF
  ENDDO 
  IF(k.EQ.1)THEN
    metal_2 = mec(2)
  ELSE IF(k.EQ.19)THEN
    metal_2 = mec(18)
  ELSE IF(smetal.GT.0.)THEN
    metal_2 = mec(k+1)
  ELSE
    metal_2 = mec(k-1)
  ENDIF
  vturb = 2.
  IF(spectro(4).LT.vturb)THEN
     vturb_2 = 1
     svturb = spectro(4) - vturb
  ENDIF
  IF(spectro(4).GE.vturb)THEN
     vturb_2 = 1
     svturb = spectro(4) - vturb
  ENDIF

  ! print*, 'Teff: ',spectro(1),teff,teff_2,steff
  ! print*, 'logg: ',spectro(2),logg,logg_2,slogg
  ! print*, 'fe/h: ',spectro(3),metal,metal_2,smetal
  ! print*, 'vturb: ',spectro(4),vturb,vturb_2,svturb

  OPEN(UNIT=12,FILE='/home/bonfanti/Documents/PostDocLiegi/MCMC_Andrea/Lib/quadratic.dat')
 
  test2=0
  k=0
  DO WHILE(test2.LT.6)
    k=k+1
    IF(k.GE.np)THEN
       test2=6
       PRINT*, 'Quadratic coefficients not found in table';! STOP
    ENDIF
    READ(12,*) rlogg,rteff,rmet,rturb,ru1,ru2,rfilt
    IF(filter.EQ.rfilt)THEN
      IF(ABS(rturb-vturb).LT.errtol &
      & .AND.ABS(rlogg-logg).LT.errtol &
      & .AND.ABS(rteff-teff).LT.errtol &
      & .AND.ABS(rmet-metal).LT.errtol)THEN
        coeff(1,1)=ru1
        coeff(2,1)=ru2
        test2=test2+1
      ENDIF
      IF(ABS(rturb-vturb).LT.errtol &
      & .AND.ABS(rlogg-logg).LT.errtol &
      & .AND.ABS(rteff-teff_2).LT.errtol &
      & .AND.ABS(rmet-metal).LT.errtol)THEN
        coeff(1,2)=ru1
        coeff(2,2)=ru2
        test2=test2+1
      ENDIF
      IF(ABS(rturb-vturb).LT.errtol &
      & .AND.ABS(rlogg-logg_2).LT.errtol &
      & .AND.ABS(rteff-teff).LT.errtol &
      & .AND.ABS(rmet-metal).LT.errtol)THEN
        coeff(1,3)=ru1
        coeff(2,3)=ru2
        test2=test2+1
      ENDIF
      IF(ABS(rturb-vturb).LT.errtol &
      & .AND.ABS(rlogg-logg).LT.errtol &
      & .AND.ABS(rteff-teff).LT.errtol &
      & .AND.ABS(rmet-metal_2).LT.errtol)THEN
        coeff(1,4)=ru1
        coeff(2,4)=ru2
        test2=test2+1
      ENDIF
      IF(ABS(rturb-vturb_2).LT.errtol &
      & .AND.ABS(rlogg-logg).LT.errtol &
      & .AND.ABS(rteff-teff).LT.errtol &
      & .AND.ABS(rmet).LT.errtol)THEN
        coeff(1,5)=ru1
        coeff(2,5)=ru2
        test2=test2+1
      ENDIF
      IF(ABS(rturb-vturb).LT.errtol &
      & .AND.ABS(rlogg-logg).LT.errtol &
      & .AND.ABS(rteff-teff).LT.errtol &
      & .AND.ABS(rmet).LT.errtol)THEN
        coeff(1,6)=ru1
        coeff(2,6)=ru2
        test2=test2+1
      ENDIF
    ENDIF   
  ENDDO  
  CLOSE(12)

  DO i=1,2
      der(i,1) = (coeff(i,1)-coeff(i,2))/(teff-teff_2) ! Derivee Teff
      der(i,2) = (coeff(i,1)-coeff(i,3))/(logg-logg_2) ! Derivee logg
      der(i,3) = (coeff(i,1)-coeff(i,4))/(metal-metal_2) ! Derivee Fe/H
      der(i,4) = (coeff(i,6)-coeff(i,5))/(vturb-vturb_2) ! Derivee vturb
      a(i) = coeff(i,1) + der(i,1)*steff + der(i,2)*slogg + &
    &    der(i,3)*smetal + der(i,4)*svturb 
      ea(i) = SQRT((der(i,1)*espectro(1))**2 + (der(i,2)*espectro(2))**2 + &
    &    (der(i,3)*espectro(3))**2+ (der(i,4)*espectro(4))**2)
  ENDDO

  RETURN
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE nlinter(spectro,espectro,filter,a,ea)

  IMPLICIT NONE
  
  INTEGER :: i,k,test2
  INTEGER, PARAMETER :: np=147300 !135516

  DOUBLE PRECISION, DIMENSION(19), PARAMETER :: mec = (/ -5.0,-4.5,-4.0, &
  & -3.5,-3.0,-2.5,-2.0,-1.5,-1.0,-0.5,-0.3,-0.2,-0.1,0.,0.1,0.2,  &
  &   0.3,0.5,1.0 /)
  DOUBLE PRECISION, DIMENSION(5), PARAMETER :: vit = (/0., 1., 2., 4., 8./)
  DOUBLE PRECISION, DIMENSION(27) :: hot =  (/ 3500,3750,4000,4250,4500,4750, &
  &  5000,5250,5500,5750,6000,6250,6500,6750,7000,7250,7500,7750,8000, &
  & 8250,8500,8750,9000,9250,9500,9750,10000 /)
  DOUBLE PRECISION, DIMENSION(10) :: grav=(/ 0.5,1.0,1.5,2.0,2.5,3.0,&
    & 3.5,4.0,4.5,5.0 /)
  DOUBLE PRECISION, DIMENSION(4) ::  spectro,espectro
  DOUBLE PRECISION, DIMENSION(4,4) :: der
  DOUBLE PRECISION, DIMENSION(4) ::  a,ea
  DOUBLE PRECISION ::  dif
  DOUBLE PRECISION ::  teff,steff,teff_2
  DOUBLE PRECISION ::  logg,slogg,logg_2
  DOUBLE PRECISION ::  metal,smetal,metal_2
  DOUBLE PRECISION ::  vturb,svturb,vturb_2
  DOUBLE PRECISION ::  errtol = 0.01
  DOUBLE PRECISION ::  rteff,rlogg,rturb,rmet
  DOUBLE PRECISION ::  ru1,ru2,ru3,ru4
  DOUBLE PRECISION, DIMENSION(4,6) :: coeff

  CHARACTER(LEN=2) :: filter,rfilt

  dif=1000.
  DO i=1,27
    IF(ABS(hot(i)-spectro(1)).LT.dif)THEN
      teff = hot(i)
      k=i
      steff = spectro(1)-teff
      dif = ABS(steff)
    ENDIF
  ENDDO 
  IF(k.EQ.1)THEN
    teff_2 = hot(2)
  ELSE IF(k.EQ.27)THEN
    teff_2 = hot(26)
  ELSE IF(steff.GT.0.)THEN
    teff_2 = hot(k+1)
  ELSE
    teff_2 = hot(k-1)
  ENDIF
  dif = 1000.
  DO i=1,10
    IF(ABS(grav(i)-spectro(2)).LT.dif)THEN
      logg = grav(i)
      k=i
      slogg = spectro(2)-logg
      dif = ABS(slogg)
    ENDIF
  ENDDO 
  IF(k.EQ.1)THEN
    logg_2 = grav(2)
  ELSE IF(k.EQ.10)THEN
    logg_2 = grav(9)
  ELSE IF(steff.GT.0.)THEN
    logg_2 = grav(k+1)
  ELSE
    logg_2 = grav(k-1)
  ENDIF
  dif = 1000.
  DO i=1,19
    IF(ABS(mec(i)-spectro(3)).LT.dif)THEN
      metal = mec(i)
      k=i
      smetal = spectro(3)-metal
      dif = ABS(smetal)
    ENDIF
  ENDDO 
  IF(k.EQ.1)THEN
    metal_2 = mec(2)
  ELSE IF(k.EQ.19)THEN
    metal_2 = mec(18)
  ELSE IF(smetal.GT.0.)THEN
    metal_2 = mec(k+1)
  ELSE
    metal_2 = mec(k-1)
  ENDIF
  vturb = 2.
  IF(spectro(4).le.vturb)THEN
     vturb_2 = 1
     svturb = spectro(4) - vturb
  ENDIF
  IF(spectro(4).gt.vturb)THEN
     vturb_2 = 4
     svturb = spectro(4) - vturb
  ENDIF
!  print*, 'Teff: ',spectro(1),teff,teff_2,steff
!  print*, 'logg: ',spectro(2),logg,logg_2,slogg
!  print*, 'fe/h: ',spectro(3),metal,metal_2,smetal
!  print*, 'vturb: ',spectro(4),vturb,vturb_2,svturb

  OPEN(UNIT=12,FILE='/home/bonfanti/Documents/PostDocLiegi/MCMC_Andrea/Lib/nonlinear.dat')

  test2=0
  k=0
  DO WHILE(test2.LT.6)
    k=k+1
    IF(k.GE.np)THEN
       test2=6
       PRINT*, 'Non-linear coefficients not found in table';STOP
    ENDIF
    READ(12,*) rlogg,rteff,rmet,rturb,ru1,ru2,ru3,ru4,rfilt
    IF(filter.EQ.rfilt)THEN
      IF(ABS(rturb-vturb).LT.errtol &
      & .AND.ABS(rlogg-logg).LT.errtol &
      & .AND.ABS(rteff-teff).LT.errtol &
      & .AND.ABS(rmet-metal).LT.errtol)THEN
        coeff(1,1)=ru1
        coeff(2,1)=ru2
        coeff(3,1)=ru3
        coeff(4,1)=ru4
        test2=test2+1
      ENDIF
      IF(ABS(rturb-vturb).LT.errtol &
      & .AND.ABS(rlogg-logg).LT.errtol &
      & .AND.ABS(rteff-teff_2).LT.errtol &
      & .AND.ABS(rmet-metal).LT.errtol)THEN
        coeff(1,2)=ru1
        coeff(2,2)=ru2
        coeff(3,2)=ru3
        coeff(4,2)=ru4
        test2=test2+1
      ENDIF
      IF(ABS(rturb-vturb).LT.errtol &
      & .AND.ABS(rlogg-logg_2).LT.errtol &
      & .AND.ABS(rteff-teff).LT.errtol &
      & .AND.ABS(rmet-metal).LT.errtol)THEN
        coeff(1,3)=ru1
        coeff(2,3)=ru2
        coeff(3,3)=ru3
        coeff(4,3)=ru4
        test2=test2+1
      ENDIF
      IF(ABS(rturb-vturb).LT.errtol &
      & .AND.ABS(rlogg-logg).LT.errtol &
      & .AND.ABS(rteff-teff).LT.errtol &
      & .AND.ABS(rmet-metal_2).LT.errtol)THEN
        coeff(1,4)=ru1
        coeff(2,4)=ru2
        coeff(3,4)=ru3
        coeff(4,4)=ru4
        test2=test2+1
      ENDIF
      IF(ABS(rturb-vturb_2).LT.errtol &
      & .AND.ABS(rlogg-logg).LT.errtol &
      & .AND.ABS(rteff-teff).LT.errtol &
      & .AND.ABS(rmet).LT.errtol)THEN
        coeff(1,5)=ru1
        coeff(2,5)=ru2
        coeff(3,5)=ru3
        coeff(4,5)=ru4
        test2=test2+1
      ENDIF
      IF(ABS(rturb-vturb).LT.errtol &
      & .AND.ABS(rlogg-logg).LT.errtol &
      & .AND.ABS(rteff-teff).LT.errtol &
      & .AND.ABS(rmet).LT.errtol)THEN
        coeff(1,6)=ru1
        coeff(2,6)=ru2
        coeff(3,6)=ru3
        coeff(4,6)=ru4
        test2=test2+1
      ENDIF
    ENDIF   
  ENDDO  
  CLOSE(12)

  DO i=1,4
      der(i,1) = (coeff(i,1)-coeff(i,2))/(teff-teff_2) ! Derivee Teff
      der(i,2) = (coeff(i,1)-coeff(i,3))/(logg-logg_2) ! Derivee logg
      der(i,3) = (coeff(i,1)-coeff(i,4))/(metal-metal_2) ! Derivee Fe/H
      der(i,4) = (coeff(i,6)-coeff(i,5))/(vturb-vturb_2) ! Derivee vturb
      a(i) = coeff(i,1) + der(i,1)*steff + der(i,2)*slogg + &
    &    der(i,3)*smetal + der(i,4)*svturb 
      ea(i) = SQRT((der(i,1)*espectro(1))**2 + (der(i,2)*espectro(2))**2 + &
    &    (der(i,3)*espectro(3))**2+ (der(i,4)*espectro(4))**2)
  ENDDO

  RETURN
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

SUBROUTINE sort(n,arr)

  IMPLICIT NONE

  INTEGER :: n
  INTEGER, PARAMETER :: m=7,nstack=50
  DOUBLE PRECISION :: arr(n)
  INTEGER :: i,ir,j,jstack,k,l,istack(nstack)
  DOUBLE PRECISION :: a,temp
  
  jstack=0
  l=1
  ir=n
1 IF((ir-l).LT.M)THEN
    DO j=l+1,ir
      a=arr(j)
      DO i=j-1,l,-1
        IF(arr(i).LE.a)GOTO 2
        arr(i+1)=arr(i)
      ENDDO
      i=l-1
2     arr(i+1)=a
    ENDDO
    IF(jstack.EQ.0) RETURN
    ir=istack(jstack)
    l=istack(jstack-1)
    jstack=jstack-2
  ELSE
    k=(l+ir)/2
    temp=arr(k)
    arr(k)=arr(l+1)
    arr(l+1)=temp
    IF(arr(l).GT.arr(ir))THEN
      temp=arr(l)
      arr(l)=arr(ir)
      arr(ir)=temp
    ENDIF
    IF(arr(l+1).GT.arr(ir))THEN
      temp=arr(l+1)
      arr(l+1)=arr(ir)
      arr(ir)=temp
    ENDIF
    IF(arr(l).GT.arr(l+1))THEN
      temp=arr(l)
      arr(l)=arr(l+1)
      arr(l+1)=temp
    ENDIF
    i=l+1
    j=ir
    a=arr(l+1)
3   CONTINUE
    i=i+1
    IF(arr(i).LT.a)GOTO 3
4   CONTINUE
    j=j-1
    IF(arr(j).GT.a)GOTO 4
    IF(j.LT.i)GOTO 5
    temp=arr(i)
    arr(i)=arr(j)
    arr(j)=temp
    GOTO 3
5   arr(l+1)=arr(j)
    arr(j)=a
    jstack=jstack+2
    IF(jstack.GT.NSTACK) PRINT*, 'NSTACK too small in sort'
    IF((ir-i+1).GE.j-l)THEN
      istack(jstack)=ir
      istack(jstack-1)=i
      ir=j-1
    ELSE
      istack(jstack)=j-1
      istack(jstack-1)=l
      l=i
    ENDIF
  ENDIF
  GOTO 1
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE LUMC(delta,rr,kk,us,ud,lc)
  
  IMPLICIT NONE

  DOUBLE PRECISION :: a,c,alf(0:3),alfat,p,q,x,valor,alfa,beta,z,us,ud,nu,lc
  DOUBLE PRECISION :: suma,rr,u1,u2,delta,kk,limb
  INTEGER :: k,n,j
  DOUBLE PRECISION :: d(0:2000),e(0:2000)
  DOUBLE PRECISION :: gamma_log,numero,gamma

  n=1000
  u1=(us+ud)/2.d00
  u2=(us-ud)/2.d00
  u1=u1+2.d00*u2
  u2=-u2
  a=1.d00/(1.d00+kk)
  c=delta/rr
  IF(delta.GT.rr) GOTO 100
  DO k=0,2
    nu=(DBLE(K)+2.D00)/2.D00
    alf(k)=(1.D00-a)*(1.D00-a)*((1.D00-c*c)**(nu+1.D00))/(nu*GAMMA(nu+1.D00))
    p=2.D00+nu
    q=1.D00
    beta = p - q
    alfa = q - 1.D00
    z=c*c
    x=1.D00 - 2.D00 * z 
    CALL jacobi_POLY(n,alfa,beta,x,d)
    p=2.D00+nu
    q=1.D00+nu
    z=a
    beta = p - q
    alfa = q - 1.D00
    x=1.D00 - 2.D00 * z
    CALL jacobi_POLY(N,alfa,beta,X,e)
    suma=0.0
    DO j=0,n
      numero=gamma_log(DBLE(j)+nu+1.d00)-gamma_log(DBLE(j)+2.d00)
      valor=((-1)**DBLE(j))*(2.d00+2.d00*DBLE(j)+nu)*DEXP(numero)
      numero = gamma_log(DBLE(j+1)) + gamma_log(alfa+1.D00) - gamma_log(DBLE(j+1)+alfa)
      e(j) = e(j) * DEXP(numero)
      valor = valor * d(j)*e(j)*e(j)
      suma=suma+valor
    ENDDO
    alf(k)=alf(k)*suma
  ENDDO 
  limb=1.0-(u1/3.0+u2/2.0)
  alfat=(alf(0)*(1.0-u1-u2)/limb)+(alf(1)*u1/limb)+(alf(2)*u2/limb)
  lc=1.0-alfat
  GOTO 200
100 LC = 1.d00
200 RETURN
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE RVC(delta,rs,kk,us,ud,lc)

  IMPLICIT NONE

  DOUBLE PRECISION :: a,c,alf(0:3),alfat,p,q,x,valor,alfa,beta,z,us,ud,nu,lc
  DOUBLE PRECISION :: suma,l2,u1,u2,delta,rs,kk,limb
  INTEGER :: k,n,j
  DOUBLE PRECISION :: d(0:2000),e(0:2000)
  DOUBLE PRECISION :: gamma_log,numero,gamma

  l2=0.d00
  n=1000
  u1=(us+ud)/2.d00
  u2=(us-ud)/2.d00
  u1=u1+2.d00*u2
  u2=-u2
  a=1.0d00/(1.0d00+kk)
  c=delta/rs
  IF(delta.LE.rs) GOTO 400
  GOTO 200
400 DO k=0,2
      nu=(DBLE(K)+2.D00)/2.D00
      alf(k)=a*c*(1.d00-a)*(1.D00-a)*GAMMA(nu)* & 
      & ((1.D00-c*c)**(nu+1.D00))/(GAMMA(nu+2.D00)*GAMMA(nu+2.d00))
      p=3.D00+nu
      q=2.D00
      beta = p - q
      alfa = q - 1.D00
      z=c*c
      x=1.D00 - 2.D00 * z 
      CALL jacobi_POLY(n,alfa,beta,x,d)
      DO j=0,n
	numero = gamma_log(DBLE(j+1)) + gamma_log(alfa+1.D00) - gamma_log(DBLE(j+1)+alfa)
	d(j) = d(j) * dexp(numero)
      ENDDO
      p=3.D00+nu
      q=2.D00+nu
      z=a
      beta = p - q
      alfa = q - 1.D00
      x=1.D00 - 2.D00 * z
      CALL jacobi_POLY(n,alfa,beta,x,e)
      suma=0.0
      DO j=0,n
	numero=gamma_log(DBLE(j)+nu+3.d00)-gamma_log(DBLE(j)+1.d00)
	valor=((-1)**DBLE(j))*(3.d00+2.d00*DBLE(j)+nu)*dexp(numero)
	numero = gamma_log(DBLE(j+1)) + gamma_log(alfa+1.D00) &
        & - gamma_log(DBLE(j+1)+alfa)	
	e(j) = e(j) * dexp(numero)	
	valor = valor * d(j)*e(j)*e(j)
	suma=suma+valor
      ENDDO
      alf(k)=alf(k)*suma
    ENDDO
    limb=1.0-(u1/3.0+u2/2.0)
    alfat=(alf(0)*(1.0-u1-u2)/limb)+(alf(1)*u1/limb)+(alf(2)*u2/limb)
    lc=alfat
    GOTO 300
200 CONTINUE
    lc=0.0E00
300 CONTINUE
  RETURN
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE jacobi_poly(n,alpha,beta,x,cx)

  IMPLICIT NONE

  INTEGER :: n
  DOUBLE PRECISION :: alpha
  DOUBLE PRECISION :: beta
  DOUBLE PRECISION :: cx(0:n)
  DOUBLE PRECISION :: c1
  DOUBLE PRECISION :: c2
  DOUBLE PRECISION :: c3
  DOUBLE PRECISION :: c4
  INTEGER :: i
  DOUBLE PRECISION :: r_i
  DOUBLE PRECISION :: x
  IF(alpha <= -1.0D+00)THEN
    WRITE( *, '(a)' ) ' '
    WRITE( *, '(a)' ) 'JACOBI_POLY - Fatal error!'
    WRITE( *, '(a,g14.6)' ) '  Illegal input value of ALPHA = ', alpha
    WRITE( *, '(a)' ) '  But ALPHA must be greater than -1.'
    STOP
  ENDIF
  IF(beta <= -1.0D+00 )THEN
    WRITE( *, '(a)' ) ' '
    WRITE( *, '(a)' ) 'JACOBI_POLY - Fatal error!'
    WRITE( *, '(a,g14.6)' ) '  Illegal input value of BETA = ', beta
    WRITE( *, '(a)' ) '  But BETA must be greater than -1.'
    STOP
  ENDIF
  IF( n < 0 )THEN
    RETURN
  ENDIF
  cx(0) = 1.0D+00
  IF(n == 0)THEN
    RETURN
  ENDIF
  cx(1) = ( 1.0D+00 + 0.5D+00 * ( alpha + beta ) ) * x &
  &  + 0.5D+00 * ( alpha - beta )
  DO i=2,n
    r_i = REAL(i,KIND=8) 
    c1 = 2.0D+00 * r_i * ( r_i + alpha + beta ) &
    &  * ( 2.0D+00 * r_i - 2.0D+00 + alpha + beta )
    c2 = ( 2.0D+00 * r_i - 1.0D+00 + alpha + beta ) &
    &  * ( 2.0D+00 * r_i  + alpha + beta ) &
    &  * ( 2.0D+00 * r_i - 2.0D+00 + alpha + beta )
    c3 = ( 2.0D+00 * r_i - 1.0D+00 + alpha + beta ) &
    &  * ( alpha + beta ) * ( alpha - beta )
    c4 = - 2.0D+00 * ( r_i - 1.0D+00 + alpha ) &
    &  * ( r_i - 1.0D+00 + beta )  * ( 2.0D+00 * r_i + alpha + beta )
    cx(i) = ( ( c3 + c2 * x ) * cx(i-1) + c4 * cx(i-2) ) / c1
  ENDDO
  RETURN
END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE gasdev_s(harvest)
 
  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(out) :: harvest
  DOUBLE PRECISION :: rsq,v1,v2
  DOUBLE PRECISION, SAVE :: g
  LOGICAL, SAVE :: gaus_stored=.false.
  IF(gaus_stored)THEN
    harvest=g
    gaus_stored=.false.
  ELSE
    DO
      CALL RANDOM_NUMBER(v1)
      CALL RANDOM_NUMBER(v2)
      v1=2.0*v1-1.0
      v2=2.0*v2-1.0
      rsq=v1**2+v2**2
      IF(rsq.GT.0..AND.rsq.LT.1.) EXIT
    ENDDO
    rsq=SQRT(-2.0*LOG(rsq)/rsq)
    harvest=v1*rsq
    g=v2*rsq
    gaus_stored=.true.
  ENDIF  

END SUBROUTINE gasdev_s

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE pdlim(n,x,bf,l,limpd)
  
  IMPLICIT NONE
  INTEGER :: n,k,i,test,l(5)
  DOUBLE PRECISION :: x(n),bf,limpd(8)
  k=0
  test=0
  i=0
  DO WHILE(test.LT.1.AND.i.LT.n)
    i=i+1
    IF(x(i).LT.bf)THEN
      k=k+1
    ELSE
      test=2
    ENDIF 
  ENDDO
  test = k - FLOOR(0.683*k)
  limpd(1) = bf-x(test)
  test = FLOOR(0.683*(n-k))+k
  limpd(2) = x(test)-bf
  test = k - FLOOR(0.997*k)
  limpd(3) = bf-x(test)
  test = FLOOR(0.997*(n-k))+k
  limpd(4) = x(test)-bf
  limpd(5) = x(l(5))-x(l(1))
  limpd(6) = x(l(2))-x(l(5))
  limpd(7) = x(l(5))-x(l(3))
  limpd(8) = x(l(4))-x(l(5))

END SUBROUTINE pdlim

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE svdfit(x,y,sig,ndata,a,ma,u,v,w,chisq,funcs,poc,nug)

  INTEGER ma,ndata,NMAX,MMAX,nug
  DOUBLE PRECISION :: chisq,a(ma),sig(ndata),u(ndata,ma),v(ma,ma),w(ma)
  DOUBLE PRECISION :: x(ndata,nug),y(ndata),TOL,x2(nug)
  EXTERNAL funcs 
  PARAMETER (NMAX=20000,MMAX=110,TOL=1.e-13)
  INTEGER  i,j,poc(nug)
  DOUBLE PRECISION :: sum,thresh,tmp,wmax,afunc(MMAX),b(NMAX)
  a = 0.
  DO i=1,ndata
    x2 = x(i,1:nug)
    CALL funcs(x2,afunc,ma,poc)
    tmp=1./sig(i)
    DO j=1,ma
      u(i,j)=afunc(j)*tmp
    ENDDO
    b(i)=y(i)*tmp
  ENDDO
  CALL svdcmp(u,ndata,ma,ndata,ma,w,v)
  wmax=0.
  DO j=1,ma
    IF(w(j).GT.wmax) wmax=w(j)
  ENDDO
  thresh=TOL*wmax
  DO j=1,ma
    IF(w(j).LT.thresh) w(j)=0.
  ENDDO
  CALL svbksb(u,w,v,ndata,ma,ndata,ma,b,a)
  chisq=0.
  DO i=1,ndata
    x2 = x(i,1:nug)
    CALL funcs(x2,afunc,ma,poc)
    sum=0.
    DO j=1,ma
      sum=sum+a(j)*afunc(j)
    ENDDO
    chisq=chisq+((y(i)-sum)/sig(i))**2
  ENDDO
  RETURN
END 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE ebfit(x,y,sig,ndata,a,ma,u,v,w,mp,np,chisq,funcs)
  INTEGER ma,mp,ndata,np,NMAX,MMAX
  REAL*8 chisq,a(ma),sig(ndata),u(mp,np),v(np,np),w(np)
  REAL*8 x(ndata,3),y(ndata),TOL
  EXTERNAL funcs 
  PARAMETER (NMAX=10000,MMAX=50,TOL=1.e-13)
  INTEGER  i,j
  REAL*8 sum,thresh,tmp,wmax,afunc(MMAX),b(NMAX)
  DO i=1,ndata
    CALL funcs(x(i,1),x(i,2),x(i,3),afunc,ma)
    tmp=1./sig(i)
    DO j=1,ma
      u(i,j)=afunc(j)*tmp
    ENDDO
    b(i)=y(i)*tmp
  ENDDO
  CALL svdcmp(u,ndata,ma,mp,np,w,v)
  wmax=0.
  DO j=1,ma
    IF(w(j).GT.wmax) wmax=w(j)
  ENDDO
  thresh=TOL*wmax
  DO j=1,ma
    IF(w(j).LT.thresh) w(j)=0.
  ENDDO
  CALL svbksb(u,w,v,ndata,ma,mp,np,b,a)
  chisq=0.
  DO i=1,ndata
    CALL funcs(x(i,1),x(i,2),x(i,3),afunc,ma)
    sum=0.
    DO j=1,ma
      sum=sum+a(j)*afunc(j)
    ENDDO
    chisq=chisq+((y(i)-sum)/sig(i))**2
  ENDDO
  RETURN
END 


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE svdvar(v,ma,np,w,cvm,ncvm)
  INTEGER :: ma,ncvm,np,MMAX
  DOUBLE PRECISION :: cvm(ncvm,ncvm),v(np,np),w(np)
  PARAMETER (MMAX=102)
  INTEGER :: i,j,k
  DOUBLE PRECISION :: sum,wti(MMAX)
  DO i=1,ma
    wti(i)=0.
    IF(ABS(w(i)).GT.1.E-13) wti(i)=1./(w(i)*w(i))
  ENDDO 
  DO i=1,ma
    DO j=1,i
      sum=0.
      DO k=1,ma
        sum=sum+v(i,k)*v(j,k)*wti(k)
      ENDDO 
      cvm(i,j)=sum
      cvm(j,i)=sum
    ENDDO
  ENDDO
  RETURN
END 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE svbksb(u,w,v,m,n,mp,np,b,x)
  INTEGER :: m,mp,n,np,NMAX
  DOUBLE PRECISION :: b(mp),u(mp,np),v(np,np),w(np),x(np)
  PARAMETER (NMAX=20000)
  INTEGER :: i,j,jj
  DOUBLE PRECISION :: s,tmp(NMAX)
  DO j=1,n
    s=0.
    IF(ABS(w(j)).GT.1.E-13)THEN
      DO i=1,m
        s=s+u(i,j)*b(i)
      ENDDO
      s=s/w(j)
    ENDIF
    tmp(j)=s
  ENDDO
  DO j=1,n
    s=0.
    DO jj=1,n
      s=s+v(j,jj)*tmp(jj)
    ENDDO
    x(j)=s
  ENDDO
  RETURN
END 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE svdcmp(a,m,n,mp,np,w,v)

  INTEGER::  m,mp,n,np,NMAX
  DOUBLE PRECISION :: a(mp,np),v(np,np),w(np)
  PARAMETER (NMAX=20000)
  INTEGER :: i,its,j,jj,k,l,nm
  DOUBLE PRECISION :: anorm,c,f,g,h,s,scale,x,y,z,rv1(NMAX),pythag
  g=0.0
  scale=0.0
  anorm=0.0
  DO i=1,n
     l=i+1
     rv1(i)=scale*g
     g=0.0
     s=0.0
     scale=0.0
     IF(i <= m)THEN
        scale=SUM(ABS(a(i:m,i)))
        IF(scale /= 0.0)THEN
           a(i:m,i)=a(i:m,i)/scale
           DO k=i,m
              s=s+a(k,i)*a(k,i)
           ENDDO
           f=a(i,i)
           g=-SIGN(SQRT(s),f)
           h=f*g-s
           a(i,i)=f-g
           DO j=l,n
              s=0.0
              DO k=i,m
                 s=s+a(k,i)*a(k,j)
              ENDDO
              f=s/h
              DO k=i,m
                 a(k,j)=a(k,j)+f*a(k,i)
              ENDDO
           ENDDO
           a(i:m,i)=scale*a(i:m,i)
        ENDIF
     ENDIF
     w(i)=scale*g
     g=0.0
     s=0.0
     scale=0.0
     IF((i.LE.m).AND.(i.NE.n))THEN
        DO k=l,n
           scale=scale+abs(a(i,k))
        ENDDO
        IF(ABS(scale).GT.1.E-13)THEN
           DO k=l,n
              a(i,k)=a(i,k)/scale
              s=s+a(i,k)*a(i,k)
           ENDDO
           f=a(i,l)
           g=-SIGN(SQRT(s),f)
           h=f*g-s
           a(i,l)=f-g
           DO k=l,n
              rv1(k)=a(i,k)/h
           ENDDO
           DO j=l,m
              s=0.0
              DO k=l,n
                 s=s+a(j,k)*a(i,k)
              ENDDO
              DO k=l,n
                 a(j,k)=a(j,k)+s*rv1(k)
              ENDDO
           ENDDO
           DO k=l,n
              a(i,k)=scale*a(i,k)
           ENDDO
        ENDIF
     ENDIF
     anorm=MAX(anorm,(ABS(w(i))+ABS(rv1(i))))
  ENDDO
  DO i=n,1,-1
     IF(i.LT.n)THEN
        IF(ABS(g).GT.1.E-13)THEN
           DO j=l,n
              v(j,i)=(a(i,j)/a(i,l))/g
           ENDDO
           DO j=l,n
              s=0.0
              DO k=l,n
                 s=s+a(i,k)*v(k,j)
              ENDDO
              DO k=l,n
                 v(k,j)=v(k,j)+s*v(k,i)
              ENDDO
           ENDDO
        ENDIF
        DO j=l,n
           v(i,j)=0.0
           v(j,i)=0.0
        ENDDO
     ENDIF
     v(i,i)=1.0
     g=rv1(i)
     l=i
  ENDDO
  DO i=MIN(m,n),1,-1
     l=i+1
     g=w(i)
     DO j=l,n
        a(i,j)=0.0
     ENDDO
     IF(ABS(g).GT.1.E-13)THEN
        g=1.0/g
        DO j=l,n
           s=0.0
           DO k=l,m
              s=s+a(k,i)*a(k,j)
           ENDDO
           f=(s/a(i,i))*g
           DO k=i,m
              a(k,j)=a(k,j)+f*a(k,i)
           ENDDO
        ENDDO
        DO j=i,m
           a(j,i)=a(j,i)*g
        ENDDO
     ELSE
        DO j=i,m
           a(j,i)=0.0
        ENDDO
     ENDIF
     a(i,i)=a(i,i)+1.0
  ENDDO
  DO k=n,1,-1
     DO its=1,30
        DO l=k,1,-1
           nm=l-1
           IF((ABS(rv1(l))+anorm).EQ.anorm) GOTO 2
           IF((ABS(w(nm))+anorm).EQ.anorm) GOTO 1
        ENDDO
1       c=0.0
        s=1.0
        DO i=l,k
           f=s*rv1(i)
           rv1(i)=c*rv1(i)
           IF((ABS(f)+anorm).EQ.anorm) GOTO 2
           g=w(i)
           h=pythag(f,g)
           w(i)=h
           h=1.0/h
           c= (g*h)
           s=-(f*h)
           DO j=1,m
              y=a(j,nm)
              z=a(j,i)
              a(j,nm)=(y*c)+(z*s)
              a(j,i)=-(y*s)+(z*c)
           ENDDO
        ENDDO
2       z=w(k)
        IF(l.EQ.k)THEN
           IF(z.LT.0.0)THEN
              w(k)=-z
              DO j=1,n
                 v(j,k)=-v(j,k)
              ENDDO
           ENDIF
           GOTO 3
        ENDIF
!        IF(its.EQ.30)THEN
!           PRINT*, 'no convergence in svdcmp';STOP
!        ENDIF
        x=w(l)
        nm=k-1
        y=w(nm)
        g=rv1(nm)
        h=rv1(k) 
        f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y)
        g=pythag(f,DBLE(1.0))
        f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x
        c=1.0
        s=1.0
        DO j=l,nm
           i=j+1
           g=rv1(i) 
           y=w(i)
           h=s*g
           g=c*g
           z=pythag(f,h)
           rv1(j)=z
           c=f/z
           s=h/z
           f= (x*c)+(g*s)
           g=-(x*s)+(g*c)
           h=y*s
           y=y*c
           DO jj=1,n
              x=v(jj,j)
              z=v(jj,i)
              v(jj,j)= (x*c)+(z*s)
              v(jj,i)=-(x*s)+(z*c)
           ENDDO
           z=pythag(f,h)
           w(j)=z
           IF(ABS(z).GT.1.E-13)THEN
              z=1.0/z
              c=f*z
              s=h*z
           ENDIF
           f= (c*g)+(s*y)
           x=-(s*g)+(c*y)
           DO jj=1,m
              y=a(jj,j)
              z=a(jj,i)
              a(jj,j)= (y*c)+(z*s)
              a(jj,i)=-(y*s)+(z*c)
           ENDDO
        ENDDO
        rv1(l)=0.0
        rv1(k)=f
        w(k)=x  
     ENDDO
3    CONTINUE
  ENDDO
  RETURN
END SUBROUTINE svdcmp

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE rvglo(x,p,np,n)
  INTEGER ::  np,n(4)
  DOUBLE PRECISION :: x(4),p(np)
  p(1)=1.
  IF(n(1).GE.1) p(2)=x(1)
  IF(n(1).GE.2) p(3)=x(1)**2
  IF(n(2).GE.1) p(4)=x(2)
  IF(n(3).GE.1) p(5)=x(3)
  IF(n(4).GE.1) p(6)=x(4)
  RETURN
END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE baseline(x,p,np,n) 

  INTEGER :: np,n(16)
  DOUBLE PRECISION :: x(33),p(np)

  p = 0.
  p(1) = 1.
  IF(n(1).GE.1) p(2)=x(1)                             ! TIME
  IF(n(1).GE.2) p(3)=x(1)**2
  IF(n(1).GE.3) p(4)=x(1)**3
  IF(n(1).GE.4) p(5)=x(1)**4

  IF(n(2).GE.1) p(6)=x(2)                             ! LOG10(AIRMASS)
  IF(n(2).GE.2) p(7)=x(2)**2
  IF(n(2).GE.3) p(8)=x(2)**3
  IF(n(2).GE.4) p(9)=x(2)**4

  IF(n(3).GE.1) p(10)=x(3)                             ! FWHM
  IF(n(3).GE.2) p(11)=x(3)**2
  IF(n(3).GE.3) p(12)=x(3)**3
  IF(n(3).GE.4) p(13)=x(3)**4

  IF(n(4).GE.1) p(14)=x(4)                             ! FWHM_x
  IF(n(4).GE.2) p(15)=x(4)**2
  IF(n(4).GE.3) p(16)=x(4)**3
  IF(n(4).GE.4) p(17)=x(4)**4 

  IF(n(5).GE.1) p(18)=x(5)                             ! FWHM_y
  IF(n(5).GE.2) p(19)=x(5)**2
  IF(n(5).GE.3) p(20)=x(5)**3
  IF(n(5).GE.4) p(21)=x(5)**4 

  IF(n(6).GE.1) p(22)=x(6)                             ! BACKGROUND
  IF(n(6).GE.2) p(23)=x(6)**2
  IF(n(6).GE.3) p(24)=x(6)**3
  IF(n(6).GE.4) p(25)=x(6)**4 

  IF(n(7).GE.1)THEN                                     ! POSITION
    p(26)=x(7)
    p(27)=x(8)
  ENDIF
  IF(n(7).GE.2)THEN
    p(28)=x(7)**2
    p(29)=x(8)**2
    p(30)=x(7)*x(8)
  ENDIF
  IF(n(7).GE.3)THEN
    p(31)=x(7)**3
    p(32)=x(8)**3
    p(33)=x(7)*x(8)**2
    p(34)=x(8)*x(7)**2
  ENDIF
  IF(n(7).GE.4)THEN
    p(35)=x(7)**4
    p(36)=x(8)**4
    p(37)=x(8)*x(7)**3
    p(38)=x(7)*x(8)**3
    p(39)=(x(7)**2)*(x(8)**2)
  ENDIF 

  IF(n(8).GE.1)THEN                                       ! CROSS-TERMS FWHM+POSITION
    p(40)=x(7)*x(3)**2
    p(41)=x(8)*x(3)**2
    p(42)=x(7)*x(4)**2
    p(43)=x(8)*x(4)**2
    p(44)=x(7)*x(5)**2
    p(45)=x(8)*x(5)**2
  ENDIF

  IF(n(9).GE.1) p(46)=x(9)                                ! SINUS
  IF(n(9).GE.2) p(47)=x(10) 
  IF(n(9).GE.3) p(48)=x(11)
  IF(n(9).GE.4) p(49)=x(12)

  IF(n(10).GE.1) p(50)=x(17)                                ! RAMP
  IF(n(10).GE.2) p(51)=x(18)

  IF(n(13).GE.1.AND.x(13).GT.0.) p(52)=1.                  ! OFFSETS
  IF(n(13).GE.2.AND.x(14).GT.0.) p(53)=1.
  IF(n(13).GE.3.AND.x(15).GT.0.) p(54)=1.
  IF(n(13).GE.4.AND.x(16).GT.0.) p(55)=1.

  IF(n(14).GE.1) p(56)=x(19)                                ! GROUP MODEL
  IF(n(14).GE.2) p(57)=x(20)
  IF(n(14).GE.3) p(58)=x(21)
  IF(n(14).GE.4) p(59)=x(22)

  IF(n(15).GT.0)THEN                                        ! PLD MODEL
    p(60:68)=x(23:31)
  ENDIF
  IF(n(15).GT.1)THEN 
    p(69:77)=x(23:31)**2
  ENDIF

  IF(n(11).GE.1)THEN                                       ! JUMP
    IF(n(12).GE.1.) p(78)=x(32)
    IF(n(12).GE.2.) p(79)=x(32)**2 
  ENDIF
  IF(n(11).GE.2)THEN                                       ! JUMP
    IF(n(12).GE.1.) p(80)=x(33)
    IF(n(12).GE.2.) p(81)=x(33)**2 
  ENDIF

  RETURN

END SUBROUTINE baseline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE rvbaseline(x,p,np,n) 

  INTEGER :: np,n(5)
  DOUBLE PRECISION :: x(5),p(np)
  p = 0.
  p(1) = 1.
  IF(n(1).GE.1) p(2)=x(1)
  IF(n(1).GE.2) p(3)=x(1)*x(1)
  IF(n(1).GE.3) p(4)=x(1)*x(1)*x(1)
  IF(n(1).GE.4) p(5)=x(1)*x(1)*x(1)*x(1)
  IF(n(2).GE.1) p(6)=x(2)
  IF(n(2).GE.2) p(7)=x(2)*x(2)
  IF(n(2).GE.3) p(8)=x(2)*x(2)*x(2)
  IF(n(2).GE.4) p(9)=x(2)*x(2)*x(2)*x(2)
  IF(n(3).GE.1) p(10)=x(3)
  IF(n(3).GE.2) p(11)=x(3)*x(3)
  IF(n(3).GE.3) p(12)=x(3)*x(3)*x(3)
  IF(n(3).GE.4) p(13)=x(3)*x(3)*x(3)*x(3)
  IF(n(4).GE.1) p(14)=x(4)
  IF(n(4).GE.2) p(15)=x(4)*x(4)
  IF(n(4).GE.3) p(16)=x(4)*x(4)*x(4)
  IF(n(4).GE.4) p(17)=x(4)*x(4)*x(4)*x(4)
  IF(n(5).GE.1) p(18)=x(5)
  IF(n(5).GE.2) p(19)=x(5)*x(5)
  IF(n(5).GE.3) p(20)=x(5)*x(5)*x(5)
  IF(n(5).GE.4) p(21)=x(5)*x(5)*x(5)*x(5)
  RETURN

END SUBROUTINE rvbaseline

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE enochfor(temp,met,dens2,a,mass)

  IMPLICIT NONE
  DOUBLE PRECISION :: temp,met,dens,dens2,mass
  DOUBLE PRECISION, DIMENSION(7) :: a
  DOUBLE PRECISION :: x,logm

  dens = LOG10(dens2)
  x = LOG10(temp)-4.1
  logm=a(1)+a(2)*x+a(3)*x*x+a(4)*dens+a(5)*dens*dens+ &
       & a(6)*dens*dens*dens+a(7)*met
  mass=10**(logm)
  RETURN

END SUBROUTINE enochfor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE gelmancomp(na,nchain,nburn,nat,x,gelmanval)
  IMPLICIT NONE
  INTEGER :: i2,nchain,na,nburn,nat,k,t,lim1,lim2
  DOUBLE PRECISION :: avegtot,meanlength,wz,bz,vario
  DOUBLE PRECISION :: x(na),gelmanval
  DOUBLE PRECISION :: avusti(nchain),sigusti(nchain),length(nchain)
  avegtot=0.;meanlength=0.;wz=0.;bz=0.
  DO i2=1,nchain 
     k=0
     lim1 = nat*(i2-1)+nburn+1
     lim2 = i2*nat
     avusti(i2) = 0.
     sigusti(i2)= 0.
     DO t=lim1,lim2
       k=k+1
       avusti(i2)=avusti(i2)+x(t)
     ENDDO
     avusti(i2)=avusti(i2)/DBLE(k)
     avegtot=avegtot+avusti(i2)/DBLE(nchain)
     DO t=lim1,lim2
       sigusti(i2)=sigusti(i2)+(x(t)-avusti(i2))**2
     ENDDO
     sigusti(i2)=sigusti(i2)/DBLE(k-1)
     wz=wz+sigusti(i2)/DBLE(nchain)
     length(i2)=k 
     meanlength=meanlength+length(i2)/DBLE(nchain)
  ENDDO
  DO i2=1,nchain
     bz=bz+(DBLE(length(i2))/DBLE(nchain-1))*(avegtot-avusti(i2))**2
  ENDDO
  vario=((meanlength-1.)/meanlength)*wz + bz/meanlength
  gelmanval=SQRT(vario/wz)
  RETURN

END SUBROUTINE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE AUTOCOR(mco,param,nchain,k,r_auto)
  IMPLICIT NONE
  INTEGER :: mco,k,i,nu,nchain,l,k2,k3,npas
  DOUBLE PRECISION :: param(mco),t1,t2,r_auto
  k2=mco/nchain
  r_auto=0.
  DO l=1,nchain
    k3=(l-1)*k2
    nu=k3+k2-k
    t1=0.; t2=0.
    npas = nu-k3
    DO i=k3+1,nu
      t2 = t2 + param(i)/DBLE(npas)
    ENDDO
    DO i=k3+1,nu
      t1 = t1 + (param(i)-t2)*(param(i+k)-t2)
    ENDDO
    r_auto = r_auto + t1
  ENDDO 
  r_auto = r_auto/DBLE(nchain)
  RETURN
END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE PDFSCREEN(texb,texunit,bf_val,medi_val,limpd,gelman,gelmanval, & 
  & nchain,jump)
  IMPLICIT NONE
  INTEGER :: nchain
  DOUBLE PRECISION :: bf_val,medi_val,limpd(8),gelmanval
  CHARACTER :: gelman, jump
  CHARACTER(LEN=7) :: texunit
  CHARACTER(LEN=20) :: texb
328 FORMAT(A21,F13.8,A3,E8.2,A2,E8.2,A2,E8.2,A2,E8.2,A8)
528 FORMAT(A21,F13.8,A3,E8.2,A2,E8.2,A2,E8.2,A2,E8.2,A8,F8.4)
  IF(gelman.EQ.'y'.AND.nchain.GT.1.AND.jump.EQ.'y')THEN
     WRITE(444,528) texb,bf_val,'-',limpd(1),'+',limpd(2),'-',limpd(3),'+', &    
      &  limpd(4),texunit,gelmanval
     WRITE(445,528) texb,medi_val,'-',limpd(5),'+',limpd(6),'-', &
      &  limpd(7),'+',limpd(8),texunit,gelmanval
     WRITE(*,528) texb,medi_val,'-',limpd(5),'+',limpd(6),'-', &
      &  limpd(7),'+',limpd(8),texunit,gelmanval
   ELSE
     WRITE(444,328) texb,bf_val,'-',limpd(1),'+',limpd(2),'-',limpd(3),'+', &
      &  limpd(4),texunit
     WRITE(445,328) texb,medi_val,'-',limpd(5),'+',limpd(6),'-',&
      &  limpd(7),'+',limpd(8),texunit
     WRITE(*,328) texb,medi_val,'-',limpd(5),'+',limpd(6),'-',&
     &  limpd(7),'+',limpd(8),texunit
  ENDIF
  RETURN
END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE HISTO(inputfile,graphfile,bf_val,medi_val,xlabel,nbin)
  IMPLICIT NONE
  INTEGER :: nbin
  DOUBLE PRECISION :: bf_val, medi_val
  CHARACTER(LEN=15) :: inputfile, graphfile
  CHARACTER(LEN=20) :: xlabel
  OPEN(UNIT=925,FILE=graphfile)
  WRITE(925,*) 'plot      erase'
  WRITE(925,*) '          data '//inputfile
  WRITE(925,*) '          read c1 2'
  WRITE(925,*) '          vecminmax c1 min max'
  WRITE(925,*) '          set bins = $min, $max, ($max - $min)/',nbin
  WRITE(925,*) '          set c1m = c1 - ($max-$min)/',nbin
  WRITE(925,*) '          set c1p = c1 + ($max-$min)/',nbin
  WRITE(925,*) '          set dum= {',medi_val,' ',medi_val,'}'
  WRITE(925,*) '          set dum2 =  {0 0}'
  WRITE(925,*) '          set dum3= {',bf_val,' ',bf_val,'}'
  WRITE(925,*) '          set dum4=   {0 0}'
  WRITE(925,*) '          vecminmax c1m min2 max2'
  WRITE(925,*) '          vecminmax c1p min3 max3'
  WRITE(925,*) '          set bins2 = $min2,$max3, ($max3 - $min2)/',nbin
  WRITE(925,*) '          expand 1.1'
  WRITE(925,*) '          window 1 1 1 1'
  WRITE(925,*) '          set ac1 = histogram(c1:bins)'
  WRITE(925,*) '          vecminmax ac1 min4 max4'
  WRITE(925,*) '          set dum2[1] = $max4*2'
  WRITE(925,*) '          set dum4[1] = $max4*2'
  WRITE(925,*) '          limits bins2 ac1'
  WRITE(925,*) '          box'
  WRITE(925,*) '          ctype 5'
  WRITE(925,*) '          shade histogram 0 bins ac1 '
  WRITE(925,*) '          ctype 3'
  WRITE(925,*) '          connect dum dum2'
  WRITE(925,*) '          ctype 4'
  WRITE(925,*) '          ltype 2'
  WRITE(925,*) '          connect dum3 dum4'
  WRITE(925,*) '          ctype 2'
  WRITE(925,*) '          ltype 0'
  WRITE(925,*) '          expand 1.5'
  WRITE(925,*) '          xlabel '//xlabel
  WRITE(925,*) '          ylabel N'
  CLOSE(925)
  RETURN
END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE GRAPH2D(inputfile,graphfile,ncol_file,abs_val,ord_val,xlabel,ylabel,errorx,errory, &
   & min_x,min_y,max_x,max_y,abs_errorm,abs_errorp,ord_errorm,ord_errorp,window,cone,newf,small)
  IMPLICIT NONE
  INTEGER :: ncol_file,abs_val,ord_val,i,window,small
  INTEGER :: abs_errorm,abs_errorp,ord_errorm,ord_errorp
  DOUBLE PRECISION :: min_x,min_y,max_x,max_y
  CHARACTER(LEN=15) :: inputfile,graphfile
  CHARACTER(LEN=20) :: xlabel,ylabel
  CHARACTER(LEN=3) :: xabs,yord,e_xp,e_xm,e_yp,e_ym
  CHARACTER :: errorx,errory,cone,newf

  IF(abs_val.LT.10)THEN
    xabs='c'//CHAR(48+abs_val)
  ELSE
    xabs='c'//CHAR(48+INT(abs_val/10))//CHAR(48+MOD(abs_val,10))
  ENDIF
  IF(ord_val.LT.10)THEN
    yord='c'//CHAR(48+ord_val)
  ELSE
    yord='c'//CHAR(48+INT(ord_val/10))//CHAR(48+MOD(ord_val,10))
  ENDIF
  IF(errorx.EQ.'y')THEN
    IF(abs_errorm.LT.10)THEN
      e_xm='c'//CHAR(48+abs_errorm)
    ELSE
      e_xm='c'//CHAR(48+INT(abs_errorm/10))//CHAR(48+MOD(abs_errorm,10))
    ENDIF
    IF(abs_errorp.LT.10)THEN
      e_xp='c'//CHAR(48+abs_errorp)
    ELSE
      e_xp='c'//CHAR(48+INT(abs_errorp/10))//CHAR(48+MOD(abs_errorp,10))
    ENDIF
  ENDIF
  IF(errory.EQ.'y')THEN
    IF(ord_errorm.LT.10)THEN
      e_ym='c'//CHAR(48+ord_errorm)
    ELSE
      e_ym='c'//CHAR(48+INT(ord_errorm/10))//CHAR(48+MOD(ord_errorm,10))
    ENDIF
    IF(ord_errorp.LT.10)THEN
      e_yp='c'//CHAR(48+ord_errorp)
    ELSE
      e_yp='c'//CHAR(48+INT(ord_errorp/10))//CHAR(48+MOD(ord_errorp,10))
    ENDIF
  ENDIF
  IF(newf.EQ.'y')THEN 
    OPEN(UNIT=925,FILE=graphfile)
    WRITE(925,*) 'plot      erase'
  ELSE
    OPEN(UNIT=925,FILE=graphfile,ACCESS='append')
  ENDIF
  WRITE(925,*) '          data '//inputfile
  DO i=1,ncol_file
    IF(i.LT.10)THEN
      WRITE(925,*) '          read c'//CHAR(48+i)//' '//CHAR(48+i)
    ELSE
      WRITE(925,*) '          read c'//CHAR(48+INT(i/10))//CHAR(48+MOD(i,10))// &
       & ' '//CHAR(48+INT(i/10))//CHAR(48+MOD(i,10))
    ENDIF
  ENDDO
  WRITE(925,*) '          expand 1.1'
  IF(window.EQ.1) WRITE(925,*) '          window 1 1 1 1'
  IF(window.EQ.2) WRITE(925,*) '          window 1 2 1 2'
  IF(window.EQ.3) WRITE(925,*) '          window 1 2 1 1'
  WRITE(925,*) '          limits ',min_x,max_x,min_y,max_y
  WRITE(925,*) '          box'
  IF(small.EQ.1)THEN
    WRITE(925,*) '          expand 0.5'
    WRITE(925,*) '          ctype 6'
  ENDIF
  IF(small.EQ.2)THEN
    WRITE(925,*) '          expand 1.1'
    WRITE(925,*) '          ctype 2'
  ENDIF
  WRITE(925,*) '          ptype 10 3'
  WRITE(925,*) '          ctype 2'
  WRITE(925,*) '          expand 1.1'
  IF(cone.NE.'y')THEN
    IF(small.EQ.1)THEN
      WRITE(925,*) '          expand 0.5'
      WRITE(925,*) '          ctype 6'
    ENDIF
    IF(small.EQ.2)THEN
      WRITE(925,*) '          expand 1.1'
    ENDIF
    WRITE(925,*) '          points '//xabs//' '//yord 
  ELSE
    WRITE(925,*) '          expand 1.1'
    WRITE(925,*) '          ctype 3'
    WRITE(925,*) '          connect '//xabs//' '//yord 
    WRITE(925,*) '          ctype 2'
  ENDIF
  IF(errorx.EQ.'y')THEN
    WRITE(925,*) '          errorbar '//xabs//' '//yord//' '//e_xm//' 3'
    WRITE(925,*) '          errorbar '//xabs//' '//yord//' '//e_xp//' 1'
  ENDIF
  IF(errory.EQ.'y')THEN
    WRITE(925,*) '          errorbar '//xabs//' '//yord//' '//e_ym//' 4'
    WRITE(925,*) '          errorbar '//xabs//' '//yord//' '//e_yp//' 2'
  ENDIF
  WRITE(925,*) '          expand 1.1'
  WRITE(925,*) '          ctype 2'
  WRITE(925,*) '          ylabel '//ylabel
  WRITE(925,*) '          xlabel '//xlabel
  CLOSE(925)
  RETURN
END SUBROUTINE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION tdbtra(a)
  IMPLICIT NONE 
  INTEGER, PARAMETER :: n=40
  INTEGER :: i
  DOUBLE PRECISION :: a, tdbtra
  DOUBLE PRECISION :: jd(n),leap(n),inter(n),slope(n)

  a = a + 2450000.
  OPEN(UNIT=12,FILE='/home/bonfanti/Documents/PostDocLiegi/MCMC_Andrea/Lib/tai-utc.dat')
  DO i=1,n
     READ(12,*) jd(i),leap(i),inter(i),slope(i)
  ENDDO
  CLOSE(12)
  IF(a.LT.jd(1))THEN
     tdbtra = a + (32.184 + leap(1) + (a-2400000.5 - inter(1))*slope(1))/(24.*3600.)
  ELSE IF(a.GE.(jd(n)))THEN
     tdbtra = a + (32.184 + leap(n))/(24.*3600.)
  ELSE
     DO i=2,n
        IF(a.GE.jd(i-1).AND.a.LT.jd(i))THEN
           tdbtra = a + (32.184 + leap(i-1) + (a-2400000.5 - inter(i-1))*slope(i-1))/(24.*3600.)
        ENDIF
     ENDDO
  ENDIF
  tdbtra=tdbtra-2450000.
  RETURN
END FUNCTION tdbtra

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

SUBROUTINE enochlaw(x,y,t,p,np)
  INTEGER  np
  REAL*8 x,y,t,p(np)
  p(1)=1.
  p(2)=x
  p(3)=x*x
  p(4)=y
  p(5)=y*y
  p(6)=y*y*y
  p(7)=t
  RETURN
END

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION pythag(a,b)
  IMPLICIT NONE
  DOUBLE PRECISION :: a, b, pythag
  DOUBLE PRECISION :: absa, absb
  absa = ABS(a)
  absb = ABS(b)
  IF(absa.GT.absb)THEN
     pythag = absa*SQRT(1.0+(absb/absa)**2)
  ELSE
     IF(absb.EQ.0)THEN
        pythag = 0.0
     ELSE
        pythag = absb*SQRT(1.0+(absa/absb)**2)
     ENDIF
  ENDIF
  RETURN
END FUNCTION pythag

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION rc(x,y)

  DOUBLE PRECISION :: rc,x,y,errtol,tiny,sqrtny,tnbg,comp1
  DOUBLE PRECISION :: comp2,third,c1,c2,c3,c4

  PARAMETER (errtol=.04d0,tiny=1.69d-38,sqrtny=1.3d-19,big=3.d37, &
       & tnbg=tiny*big,comp1=2.236d0/sqrtny,comp2=tnbg*tnbg/25.d0, &
       & third=1.d0/3.d0,c1=.3d0,c2=1.d0/7.d0,c3=.375d0,c4=9.d0/22.d0)

  DOUBLE PRECISION :: alamb,ave,s,w,xt,yt
  IF(x.LT.0..OR.y.EQ.0..OR.(x+ABS(y)).lt.tiny.or.(x+ &
       & ABS(y)).GT.big.OR.(y.LT.-comp1.AND.x.GT.0. &
       & .AND.x.LT.comp2))THEN
     PRINT*, 'invalid arguments in rc'
  ENDIF
  IF(y.GT.0.d0)THEN
     xt=x
     yt=y
     w=1.
  ELSE
     xt=x-y
     yt=-y
     w=SQRT(x)/SQRT(xt)
  ENDIF
1 CONTINUE
  alamb=2.d0*sqrt(xt)*sqrt(yt)+yt
  xt=.25d0*(xt+alamb)
  yt=.25d0*(yt+alamb)
  ave=third*(xt+yt+yt)
  s=(yt-ave)/ave
  IF(abs(s).gt.errtol) GOTO 1
  rc=w*(1.d0+s*s*(c1+s*(c2+s*(c3+s*c4))))/SQRT(ave)
  RETURN
END FUNCTION rc
!  (C) Copr. 1986-92 Numerical Recipes Software

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION rj(x,y,z,p)

  DOUBLE PRECISION :: rj,p,x,y,z,errtol,tiny,big,c1,c2,c3,c4,c5,c6,c7,c8

  PARAMETER (errtol=.05d0,tiny=2.5d-13,big=9.d11,c1=3.d0/14.d0, &
       & c2=1.d0/3.d0,c3=3.d0/22.d0,c4=3.d0/26.d0,c5=.75d0*c3, &
       & c6=1.5d0*c4,c7=.5d0*c2,c8=c3+c3)
  !    USES rc,rf

  DOUBLE PRECISION :: a,alamb,alpha,ave,b,beta,delp,delx,dely
  DOUBLE PRECISION :: delz,ea,eb,ec,ed,ee,fac,pt,rcx,rho,sqrtx
  DOUBLE PRECISION :: sqrty,sqrtz,sum,tau,xt,yt,zt,rc,rf

  IF(MIN(x,y,z).LT.0..OR.MIN(x+y,x+z,y+z,ABS(p)).LT.tiny.OR.MAX(x,y, &
       & z,ABS(p)).GT.big) PRINT*, 'invalid arguments in rj'
  sum=0.d0
  fac=1.d0
  IF(p.gt.0.d0)THEN
     xt=x
     yt=y
     zt=z
     pt=p
  ELSE
     xt=min(x,y,z)
     zt=max(x,y,z)
     yt=x+y+z-xt-zt
     a=1.d0/(yt-p)
     b=a*(zt-yt)*(yt-xt)
     pt=yt+b
     rho=xt*zt/yt
     tau=p*pt/yt
     rcx=rc(rho,tau)
  ENDIF
1 CONTINUE
  sqrtx=SQRT(xt)
  sqrty=SQRT(yt)
  sqrtz=SQRT(zt)
  alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
  alpha=(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz)**2
  beta=pt*(pt+alamb)**2
  sum=sum+fac*rc(alpha,beta)
  fac=.25d0*fac
  xt=.25d0*(xt+alamb)
  yt=.25d0*(yt+alamb)
  zt=.25d0*(zt+alamb)
  pt=.25d0*(pt+alamb)
  ave=.2d0*(xt+yt+zt+pt+pt)
  delx=(ave-xt)/ave
  dely=(ave-yt)/ave
  delz=(ave-zt)/ave
  delp=(ave-pt)/ave
  IF(MAX(ABS(delx),ABS(dely),ABS(delz),ABS(delp)).GT.errtol) GOTO 1
  ea=delx*(dely+delz)+dely*delz
  eb=delx*dely*delz
  ec=delp**2
  ed=ea-3.d0*ec
  ee=eb+2.d0*delp*(ea-ec)
  rj=3.d0*sum+fac*(1.d0+ed*(-c1+c5*ed-c6*ee)+eb*(c7+delp* &
       & (-c8+delp*c4))+delp*ea*(c2-delp*c3)-c2*delp*ec)/(ave*SQRT(ave))
  IF(p.LE.0.d0) rj=a*(b*rj+3.d0*(rcx-rf(xt,yt,zt)))
  RETURN
END FUNCTION rj
!  (C) Copr. 1986-92 Numerical Recipes Software

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION ellec(k)

  IMPLICIT NONE

  DOUBLE PRECISION :: k,m1,a1,a2,a3,a4,b1,b2,b3,b4,ee1,ee2,ellec

  ! Computes polynomial approximation for the complete elliptic
  ! integral of the second kind (Hasting's approximation):
  m1=1.d0-k*k
  a1=0.44325141463d0
  a2=0.06260601220d0
  a3=0.04757383546d0
  a4=0.01736506451d0
  b1=0.24998368310d0
  b2=0.09200180037d0
  b3=0.04069697526d0
  b4=0.00526449639d0
  ee1=1.d0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
  ee2=m1*(b1+m1*(b2+m1*(b3+m1*b4)))*log(1.d0/m1)
  ellec=ee1+ee2
  RETURN
END FUNCTION ellec

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION ellk(k)

  IMPLICIT NONE

  DOUBLE PRECISION :: a0,a1,a2,a3,a4,b0,b1,b2,b3,b4,ellk, &
       &       ek1,ek2,k,m1

  ! Computes polynomial approximation for the complete elliptic
  ! integral of the first kind (Hasting's approximation):
  m1=1.d0-k*k
  a0=1.38629436112d0
  a1=0.09666344259d0
  a2=0.03590092383d0
  a3=0.03742563713d0
  a4=0.01451196212d0
  b0=0.5d0
  b1=0.12498593597d0
  b2=0.06880248576d0
  b3=0.03328355346d0
  b4=0.00441787012d0
  ek1=a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)))
  ek2=(b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*log(m1)
  ellk=ek1-ek2
  RETURN
END FUNCTION ellk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION rf(x,y,z)

  DOUBLE PRECISION :: rf,x,y,z,errtol,tiny,big,third,c1,c2,c3,c4

  PARAMETER (errtol=.08d0,tiny=1.5d-38,big=3.d37,third=1.d0/3.d0, &
       & c1=1.d0/24.d0,c2=.1d0,c3=3.d0/44.d0,c4=1.d0/14.d0)

  DOUBLE PRECISION :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty, &
       & sqrtz,xt,yt,zt

  IF(MIN(x,y,z).LT.0.d0.OR.MIN(x+y,x+z,y+z).LT.tiny.OR. &
       & MAX(x,y,z).GT.big) PRINT*, 'invalid arguments in rf'
  xt=x
  yt=y
  zt=z
1 CONTINUE
  sqrtx=SQRT(xt)
  sqrty=SQRT(yt)
  sqrtz=SQRT(zt)
  alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
  xt=.25d0*(xt+alamb)
  yt=.25d0*(yt+alamb)
  zt=.25d0*(zt+alamb)
  ave=third*(xt+yt+zt)
  delx=(ave-xt)/ave
  dely=(ave-yt)/ave
  delz=(ave-zt)/ave
  IF(MAX(ABS(delx),ABS(dely),ABS(delz)).GT.errtol) GOTO 1
  e2=delx*dely-delz**2
  e3=delx*dely*delz
  rf=(1.d0+(c1*e2-c2-c3*e3)*e2+c4*e3)/SQRT(ave)
  RETURN 
END FUNCTION rf
!  (C) Copr. 1986-92 Numerical Recipes Software 0NL&WR2.

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION gamma(x)

  IMPLICIT NONE

  DOUBLE PRECISION :: gamma
  DOUBLE PRECISION :: gamma_log
  DOUBLE PRECISION :: x

  gamma = dexp (gamma_log(x )) 
  RETURN
END FUNCTION gamma

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

FUNCTION gamma_log ( x )

  IMPLICIT NONE
  DOUBLE PRECISION, DIMENSION(7), PARAMETER :: c = (/ &
       &  -1.910444077728D-03, &
       &   8.4171387781295D-04, &
       &  -5.952379913043012D-04, &
       &   7.93650793500350248D-04, &
       &  -2.777777777777681622553D-03, &
       &   8.333333333333333331554247D-02, &
       &   5.7083835261D-03 /)
  DOUBLE PRECISION :: corr
  DOUBLE PRECISION, PARAMETER :: d1 = - 5.772156649015328605195174D-01
  DOUBLE PRECISION, PARAMETER :: d2 =   4.227843350984671393993777D-01
  DOUBLE PRECISION, PARAMETER :: d4 =   1.791759469228055000094023D+00
  DOUBLE PRECISION ::  eps
  DOUBLE PRECISION, PARAMETER :: frtbig = 1.42D+09
  INTEGER :: i
  DOUBLE PRECISION :: gamma_log
  DOUBLE PRECISION, DIMENSION(8), PARAMETER :: p1 = (/ &
       &  4.945235359296727046734888D+00, &
       &  2.018112620856775083915565D+02, &
       &  2.290838373831346393026739D+03, &
       &  1.131967205903380828685045D+04, &
       &  2.855724635671635335736389D+04, &
       &  3.848496228443793359990269D+04, &
       &  2.637748787624195437963534D+04, &
       &  7.225813979700288197698961D+03 /)
  DOUBLE PRECISION, DIMENSION(8), PARAMETER :: p2 = (/ &
       &  4.974607845568932035012064D+00, &
       &  5.424138599891070494101986D+02, &
       &  1.550693864978364947665077D+04, &
       &  1.847932904445632425417223D+05, &
       &  1.088204769468828767498470D+06, &
       &  3.338152967987029735917223D+06, &
       &  5.106661678927352456275255D+06, &
       &  3.074109054850539556250927D+06 /)
  DOUBLE PRECISION, DIMENSION(8), PARAMETER :: p4 = (/ &
       &  1.474502166059939948905062D+04, &
       &  2.426813369486704502836312D+06, &
       &  1.214755574045093227939592D+08, &
       &  2.663432449630976949898078D+09, &
       &  2.940378956634553899906876D+10, &
       &  1.702665737765398868392998D+11, &
       &  4.926125793377430887588120D+11, &
       &  5.606251856223951465078242D+11 /)
  DOUBLE PRECISION, PARAMETER :: pnt68 = 0.6796875D+00
  DOUBLE PRECISION, DIMENSION(8), PARAMETER :: q1 = (/ &
       &  6.748212550303777196073036D+01, &
       &  1.113332393857199323513008D+03, &
       &  7.738757056935398733233834D+03, &
       &  2.763987074403340708898585D+04, &
       &  5.499310206226157329794414D+04, &
       &  6.161122180066002127833352D+04, &
       &  3.635127591501940507276287D+04, &
       &  8.785536302431013170870835D+03 /)
  DOUBLE PRECISION, DIMENSION(8), PARAMETER :: q2 = (/ &
       &  1.830328399370592604055942D+02, &
       &  7.765049321445005871323047D+03, &
       &  1.331903827966074194402448D+05, &
       &  1.136705821321969608938755D+06, &
       &  5.267964117437946917577538D+06, &
       &  1.346701454311101692290052D+07, &
       &  1.782736530353274213975932D+07, &
       &  9.533095591844353613395747D+06 /)
  DOUBLE PRECISION, DIMENSION(8), PARAMETER :: q4 = (/ &
       &  2.690530175870899333379843D+03, &
       &  6.393885654300092398984238D+05, &
       &  4.135599930241388052042842D+07, &
       &  1.120872109616147941376570D+09, &
       &  1.488613728678813811542398D+10, &
       &  1.016803586272438228077304D+11, &
       &  3.417476345507377132798597D+11, &
       &  4.463158187419713286462081D+11 /)
  DOUBLE PRECISION :: res
  DOUBLE PRECISION, PARAMETER :: sqrtpi = 0.9189385332046727417803297D+00
  DOUBLE PRECISION ::  x
  DOUBLE PRECISION, PARAMETER :: xbig = 4.08D+36
  DOUBLE PRECISION :: xden
  DOUBLE PRECISION :: xm1
  DOUBLE PRECISION :: xm2
  DOUBLE PRECISION :: xm4
  DOUBLE PRECISION :: xnum
  DOUBLE PRECISION :: xsq
  !
  !  Return immediately if the argument is out of range.
  !
  IF(x <= 0.0D+00 .OR. xbig < x)THEN
     gamma_log = huge(gamma_log)
     RETURN
  ENDIF
  eps = epsilon (eps)
  IF(x <= eps)THEN
     res = -log(x)
  ELSE IF(x <= 1.5D+00)THEN
     IF(x < pnt68)THEN
        corr = -log( x )
        xm1 = x
     ELSE
        corr = 0.0D+00
        xm1 = ( x - 0.5D+00 ) - 0.5D+00
     ENDIF
     IF(x <= 0.5D+00 .OR. pnt68 <= x)THEN
        xden = 1.0D+00
        xnum = 0.0D+00
        DO i=1,8
           xnum = xnum * xm1 + p1(i)
           xden = xden * xm1 + q1(i)
        ENDDO
        res = corr + ( xm1 * ( d1 + xm1 * ( xnum / xden ) ) )
     ELSE
        xm2 = ( x - 0.5D+00 ) - 0.5D+00
        xden = 1.0D+00
        xnum = 0.0D+00
        DO i=1,8
           xnum = xnum * xm2 + p2(i)
           xden = xden * xm2 + q2(i)
        ENDDO
        res = corr + xm2 * ( d2 + xm2 * ( xnum / xden ) )
     ENDIF
  ELSE IF(x <= 4.0D+00)THEN
     xm2 = x - 2.0D+00
     xden = 1.0D+00
     xnum = 0.0D+00
     DO i = 1,8
        xnum = xnum * xm2 + p2(i)
        xden = xden * xm2 + q2(i)
     ENDDO
     res = xm2 * ( d2 + xm2 * ( xnum / xden ) )
  ELSE IF(x <= 12.0D+00)THEN
     xm4 = x - 4.0D+00
     xden = - 1.0D+00
     xnum = 0.0D+00 
     DO i=1,8
        xnum = xnum * xm4 + p4(i)
        xden = xden * xm4 + q4(i)
     ENDDO
     res = d4 + xm4 * ( xnum / xden )
  ELSE 
     res = 0.0D+00
     IF(x <= frtbig)THEN
        res = c(7)
        xsq = x * x
        DO i=1,6
           res = res / xsq + c(i)
        ENDDO
     ENDIF
     res = res/x
     corr = log(x)
     res = res + sqrtpi - 0.5D+00 * corr
     res = res + x * ( corr - 1.0D+00 )
  ENDIF
  gamma_log = res
  RETURN
END FUNCTION gamma_log
