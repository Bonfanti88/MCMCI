|-----------------------------------------------------------------------------|
|                                 README FILE                   v1.1_20200206 |
|                MCMCI code. Bonfanti & Gillon (2020) [BG20]                  |
|                  andrea.bonfanti@oeaw.ac.at                                 |
|                  michael.gillon@uliege.be                                   |
| This is just a succint guide to deal with our code. Please refer to BG20    |
| for an extensive discussion.                                                |
| Please cite BG20 if you use our code.                                       |
|-----------------------------------------------------------------------------|
| Table of contents:                                                          |
| 1) Required files                                                           |
| 2) Compilation - gfortran compilers troubleshooting                         |
| 3) How to fill in the mcmc.dat input form                                   |
| 4) Output files                                                             |
| 5) Examples                                                                 |
| 6) Last updates                                                             |
|-----------------------------------------------------------------------------|

|-----------------------------------------------------------------------------|
| 1) Required files                                                           |
|-----------------------------------------------------------------------------|
Source code:
- SunPD.f90
- IsoPD.f90
- empRel.f90
- SCPmcmcPD12S.f90
- mcmcI.f90
Input files and folder:
- mcmc.dat input form
- phot____.txt and/or RV____.txt  ____ to be filled with progressive numbers
   starting with 0001
  phot____.txt contains the photometric time-series and it is usually composed
   by 11 columns, which are:
   JD-2450000, flux, error, x, y, fwhm, fwhm_x, fwhm_y, sky, airmass, texp
   1. JD-2450000: shifted julian date
   2. flux: normalized stellar flux
   3. error: error on flux
   4. x: x location (px)
   5. y: y location (px)
   6. fwhm: PSF's FWHM (px)
   7. fwhm_x: PSF's FWHM along the x-direction (px)
   8. fwhm_y: PSF's FWHM along the y-direction (px)
   9. sky: sky signal (electrons)
   10. airmass
   11. texp: exposure time (s)
   Essential columns are the first three. If some of the other parameters are
   not available, set their corresponding column to 0.
  RV____.txt contains the radial velocity time-series and it is composed by 8
   columns, which are:
   JD-2450000, rv, error, fwhm, contrast, bis, logR'HK, texp
   1. JD-2450000: shifted julian date
   2. rv: radial velocity (km/s)
   3. error: error on rv (km/s)
   4. fwhm: FWHM of the Cross Correlation Function CCF (km/s)
   5. contrast: contrast of the CCF (%)
   6. bis: bisector velocity span (km/s)
   7. logR'HK
   8. texp: exposure time (s)
- Lib folder
Additional files and folders:
- Ziso.txt 
- Tables of isochrones: http://stev.oapd.inaf.it/cgi-bin/cmd
   to be downloaded according to logtstep that is specified in IsoPD.f90 and
   with the Z values (metallicities) that are specified in Ziso.txt
   Each grid of isochrones has to be named as "Z"+Z numerical value+".dat"
   (e.g. Z0.01.dat or Z0.017.dat)
   So far mcmcI can deal only Johnson-Cousin and Gaia photometric systems
- Tables of tracks: https://people.sissa.it/~sbressan/CAF09_V1.2S_M36_LT/
- EvoZ folder (to be put in the same folder containing Tracks)
- TracksVel folder (to be put in the same folder containing Tracks)
- ZAMSdata folder (to be put in the same folder containing Tracks)
- tauBarnesKim2010.txt

Please update all the directories that are specified in the SCPmcmcPD12S.f90
 and mcmcI.f90 source files with yours. Search for 'bonfanti' to easily find
 them. In case the CHARACTER length of some string to be modified is specified,
 update it
 '/home/bonfanti/Documents/ArticoloTesi/TracceLeo/' refers to very old Padova
  models and can be ignored


|-----------------------------------------------------------------------------|
| 2) Compilation - gfortran compilers troubleshooting                         |
|-----------------------------------------------------------------------------|
gfortran SunPD.f90 IsoPD.f90 empRel.f90 SCPmcmcPD12S.f90 -c
gfortran mcmcI.f90 SunPD.o IsoPD.o empRel.o SCPmcmcPD12S.o -o mcmcI.exe

Program has been successfully tested on machines with gfortran version >=5.3.1

Here follows a list of problems according to our experience and testing:
- Exponent at (1) must be integer for an initialization expression
 It happens with gfortran 4.2 and it has been solved since gfortran 4.5, as it
  is reported here: https://github.com/oorb/oorb/issues/1
 Possible workaround without updating the compiler is avoiding computations
  that involve non-integer exponents in the modules. In  our case, directly
  write the numerical value of Msun (in SunPD.f90) and tZAMS (in empRel.f90)
- Syntax error in CHARACTER declaration
 Errors of this kind are drawn when one tries to allocate CHARACTER that have
  been declared as ALLOCATABLE. This bug is present in gfortran 4.4.7, it’s not
  present in ifort, and it has been solved after updating to 4.8.1, as it is
  stated in
  https://stackoverflow.com/questions/20908053/allocatable-character-variables-in-fortran
- gfortran 4.6 doesn't correctly compare two equal strings if they are declared
  ALLOCATABLE. Thus the functions strcmp and strncmp do not work properly with
  that version of the compiler. Workaround for this is declaring those 
  CHARACTER as POINTER, but this could imply several parts of the code to be
  fixed(!), thus it’s strongly suggested to update the compiler. We guarantee
  that since gfortran 5.3.1 these functions behave correctly; we haven’t tested
  the behaviour in case 4.6 < gfortran < 5.3
- ONLY for MacOS user: suffix or operand invalid for ‘movq’
 This is an assembly error, which is not related to the specific gfortran
  compiler. It is likely due to the homebrew/Macport installation. Following
  this answer
  https://stackoverflow.com/questions/35298202/reason-for-suffix-or-operands-invalid-for-movq
  likely solves the problem.

Summing up, as far as we know, compilation is successful with gfortran >= 5.3.1
 (possibly considering the hacking, that is suggested at the previous point for
 MacOS users).


|-----------------------------------------------------------------------------|
| 3) How to fill in the mcmc.dat input form                                   |
|-----------------------------------------------------------------------------|
mcmc.dat input form
It's an extension of the previous input form, that was compatible with the MCMC
 code by M. Gillon.
Efforts were made to preserve backward compatibility

|--------
|OPTIONS:
|--------
NUMBER_OF_CHAINS:
	Number of chains to be considered during the MCMC
LENGTH_OF_EACH_CHAIN_(steps): 	
	Number of steps of each chain. In the following, 'link' will indicate the
	 step's index 
BURN-IN_PHASE_FRACTION_(%):		
	Fraction of each chain whose steps will be discarded because of the burn-in
STAT_LENGTH_(steps):
	Each time, after this number of steps, the acceptance rate (see e.g. Ford,
	2005) is computed
BIN_SIZE_FOR_LCs_(min):	
	Temporal binning width of photometric data. Just related to the superMongo
	plotting routines
DESIRED_SUCCESS_RATE_(%)
	Desired acceptance rate (see e.g. Ford, 2005)
SERIES_OR_NUMERICAL_SOLUTION_TO_KE_[s-n]: 
	's' eccentric anomaly from Kepler equation is retrieved through a series
	 expansion up to 4-degree in eccentricity e (not valid if e>0.6627434)
    'n' eccentric anomaly from Kepler equation is computed numerically
PRECISION_OF_KE_NUMERICAL_SOLUTION:	
	Maximum absolute difference (in radians) that is allowed between two
	consecutive steps of the iterative process which numerically solves Kepler
	equation; it is the convergence criterion of the numerical algorithm
GIBBS_SAMPLER_[y-n]+STAT_LENGTH	
	'y' implementation of the Gibbs sampler (Casella&George, 1992) during the
	 burn-in phase
	It is followed by the number of steps over which the acceptance rate is
	 computed
MAXADAPT(min=1)	
	Related to the jump width. Any factor by which perturbing jump parameter's
	uncertainties is within [1/MAXADAPT, MAXADAPT]
GELMAN-RUBIN_TEST_[y-n]	
	'y' performs the test by Gelman&Rubin (1992): it checks mutual convergence
	if several chains are launched
COMPUTATION_OF_AUTOCORRELATIONS_[y-n]	
	'y' produces files to plot the autocorrelation function (Ford, 2006)
MODEL_ROSSITER_EFFECT_[y-n]	
	'y' considers the Rossiter-McLaughlin effect model by Giménez 2006 (only if
	 RVs are available during transit)
MODEL_RV_TIDE_[y-n]	
	'y' models the tidal amplitude induced by the exoplanet on its host.
	Magnitude of the tidally induced RV signal is estimated as in Arras+ (2012)
FIX_STELLAR_PARAMS[y-n]	
	'y' stellar parameters are perturbed always starting from the value that is
	 specified in the form.
	Stellar parameters are not jump parameters
FIT_Ms+Rs[y-n]
	'y' Ms and Rs are both jump parameters and they constrain the mean stellar
	 density rhos and the scaled semi-major axis a/Rs, without considering the
	 transit
ASSUME_MsRs_RELATION_[y-n+beta+error]	
	'y' Ms is computed empirically Ms=rhos^(1/(1-3*beta)), consistently with
	eq. (10) of Seager & Mallén-Ornelas (2003), which comes from the empirical
	relation Rs=Ms^beta. 
	For example beta=0.8 for F-K MS-stars (Cox, 2000), but here any beta value
	is accepted and then perturbed at each step, thanks to error. Rs is then
	computed from Ms and rhos. Thus, both Ms and Rs aren't jump parameters
DERIVE_Ms_FROM_RHOs+Rs[y-n]	
	'y' Ms is computed from the jump radius and from rhos that is inferred from
	 transit. Rs is a jump parameter, while Ms is not
USE_RHOs+EB_LAW_FOR_Ms+Rs[y-n]	
	'y' computes Ms from the eclipsing binary (EB) law (like in Torres+ 2010,
	 but revised by Gillon+ 2011)
	Thanks to rhos from transit, then also Rs is computed. 
	Both Ms and Rs are not jump parameters
IF_EB_LOWER_LIMIT_ON_P[d]	
	Minimum allowed orbital period (in days) to apply the EB law
IF_EB_MASS_RANGE[Msun]	
	Specify the minimum and maximum mass values, so to pick up just the binary
	systems belonging to this mass range before applying the EB law
IF_EB_USE_ONLY_TORRES_STARS_[y-n]	
	'y' picks up from the list only those stars considered by Torres+ 2010
COMPUTE_MASS_FROM_ISOCHRONES_[y-n]	
	'y' considers stellar evolutionary models to constrain Ms (and possibly
	 also Rs if no info about rhos comes from transit). 
	It is the only option that implies the use of the Isochrone placement
	 algorithm and it's not compatible with the previous options.
TAKE_INTO_ACCOUNT_M2[y-n]	
	'y' considers the contribution of the exoplanetary mass as non negligible
	 when computing rhos from transit
DEDUCE_STELLAR_INCLINATION_[y-n]	
	'y' infers stellar spin axis inclination by considering both vsini and
	 rotational period (that should be both available)
PH_AMPLI_UPPER_LIMIT[y-n]	
	'y' accepts the jump A_refl (see later) only if it's not greater than the
	 occultation depth. In other words, the amplitude of the phase curve can't
	 exceed the occultation depth.
LIMB-DARKENING_LAW_[qd-nl-no]	
	'qd' quadratic model for limb-darkening (two coefficients)
	'nl' not-linear model for limb-darkening (four coefficients)
	'no' do not consider limb-darkening
RAMP_MODEL_[log-exp]	
	'log' ramp effect is modelled as in Knutson+ 2008, i.e. the time-serie is
	 fitted as a function of ln(dt), where dt is the change in time from the 
	 start of observations. This function can be either linear or quadratic,
	 according to the set ramp order (see later)
	'exp' the time-serie is fitted as a function of exp(-dt/t_ramp), where dt
	 is defined as before, while t_ramp is a normalization factor. If the exp
	 function has degree=1 only one normalization factor with its uncertainty
	 will have to be specified, while for a degree=2 exp-function two
	 normalization factors with their uncertainties will have to be specified
	 (see later t1_ramp and t2_ramp)
DDFs_[y-n-p]+SIGMA	
	'y' considers the transit depth as a function of wavelength; perturbation
	 of DDF is based on SIGMA
	'p' DDF also enters the merit function
TTVs_[y-n] 
	'y' considers also TTVs. If it is the case, a file called ttv01.dat (if,
	 for example, planet 1 is showing TTVs) has to be created, in which there's
	 a line for each transit. Each line has to report the epoch N (i.e. the
	 integer N in the formula T = T0 + N*P), an initial value for the TTV (in
	 days, e.g. 0.), and an initial value for the error (in days, e.g. 0.001).
	 Moreover, both T0 and PERIOD P have to be fixed with null errors in the
	 form (see later)
RED_NOISE_TIMESCALE_MIN_MAX(min)	
	Minimum (MIN) and maximum (MAX) values (in minutes) to be specified. The
	rms of the fit is computed at different timescales, which are established
	by binning the [MIN, MAX] range into 50 parts
FIX_RED_NOISE_TIMESCALE[y-n]+VALUE(min)	
	'y' computes the fit rms just at the timescale that is specified by VALUE
	 (to be expressed in minutes)
BLISS_MIN_N_DATA_PER_BOX	
	If the bliss mapping option (BM) is active (see later), every pixel is
	divided in a grid. Here we specify the minimum number of desired points to
	be present in a single grid's box, so to check and possibly remove the
	different intra-pixel sensitivity (Knutson+ 2008)
NUMBER_OF_PLANETS	
	Number of exoplanets (Np) around the star for which we provide at least one
	transit LC or RV time-series
MINIMAL_PERIOD_(DAYS)	
	Vector - whose length is equal to NUMBER_OF_PLANETS - containing the lower
	acceptable limits for the planetary periods. 
MAXIMAL_PERIOD_(DAYS)	
	Vector - whose length is equal to NUMBER_OF_PLANETS - containing the upper
	 acceptable limits for the planetary periods.
	If the jump P is outside the [MINIMAL_PERIOD, MAXIMAL_PERIOD] range, the
	 jump is not accepted
MAXIMAL_K_(ms)	
	Vector - whose length is equal to NUMBER_OF_PLANETS - containing the upper
	 acceptable limits for the RV semi-amplitude (in m/s)
	If the jump K is greater than MAXIMAL_K, the jump is not accepted
MINIMAL_PLANET_RADIUS_(R_earth)	
	If the planetary radius Rp is lower than this value, the jump is not
	accepted
MAXIMAL_PLANET_RADIUS_(R_earth)	
	If the planetary radius Rp is greater than this value, the jump is not
	accepted
IRON_LOWER_LIMIT_FOR_PLANETS_[y-n]	
	'y' considers the existence of a minimum allowable planetary radius Rp_iron
	 assuming that the planet is all made of iron, and rejects the jump if
	 Rp < Rp_iron
MAXIMAL_TIDAL_ELONGATION	
	Maximum ratio between equatorial and polar radius of the star. If the jump
	TIDAL_ELONGATION (see later) is greater then MAXIMAL_TIDAL_ELONGATION, the
	jump is not accepted

|--------
|DATA:
|--------
NUMBER_OF_LCs:	
	Number of light-curves that are considered
9, sinus numbers	
	If set to k.GT.0 (and up to 4), add a line that reports k sinus periods 
	followed by k errors, which enter the baseline trend model; the code will
	explore the phase and amplitude by itself.
10, ramp order		
	If set to 1 and RAMP_MODEL='exp', add a line that reports t1_ramp and its
	 error
	If set to 2 and RAMP_MODEL='exp', add a line that reports t1_ramp + error,
	 and t2_ramp + error
11, jump numbers	
	If set to k.GT.0 (up to 4), add a line that reports k values of jump time
13, offset numbers	
	If set to k.GT.0 (up to 4), add a line that reports k values of time offset
14, PLD	
	Pixel Level Decorrelation (ONLY for Spitzer Space Telescope). The integer
	 (up to 2) to be specified is the degree of the polynomial that enters the
	 baseline trend model aiming at removing the systematic noise caused by the
	 spacecraft jitter.
	If PLD.GT.0, 9 further columns are expected in the corresponding phot file
15, flare numbers	
	If set to k.GT.0 (up to 4), add a line that reports k flt0 + k errors, k
	flampli + k errors, and k fltau + k errors, so to implement the following
	flaremodel for a given LC and for a given flare:
	flaremodel = flampli*exp[-(bjd-flt0)/fltau]
CF
	To be set to 1.00 during the 1st MCMC run. Then it has to be set equal to
	the CorF value reported in the output mcmc_med.res file, before launching
	the 2nd MCMC run
O	
	'y' if temporal cadence is too long, it considers further time values to
	 build the modelled LC. These further time values are added at steps of BO
	 (in seconds)
BO
	Temporal step's width (in seconds) in case O='y'
BM	
	'y' performs the bliss-mapping
OF	
	'y' accepts phot files according to old format; they were made of 9 columns
	 which were:
	 JD-2450000, flux, error, x, y, airmass, fwhm, sky, texp

NUMBER_OF_RV_TIME-SERIES_+_TREND_ORDER:	
	Number of RV time-series that are considered, followed by the global trend
	order. This global trend aims to have all RV time-series sharing the same
	trend, while the trend options that follow "Jit   1 2 3 4 5 O BO UT Rem"
	are individual, i.e. there is one line for each time-series.

NUMBER_OF_TRANSIT_TIMINGS:	
	'y' then a file timing.dat listing the transit timings and errors must be
	 provided. If it is the case, the Bayesian Penalty involving the provided
	 timings and T0 will enter the merit function
NUMBER_OF_OCCULTATION_TIMINGS:	
	'y' then a file timingoc.dat listing the occultation timings and errors
	 must be provided. If it is the case, the Bayesian Penalty involving the
	 provided timings and the occultation timings inferred from T0 will enter
	 the merit function

Please refer to the readme section in the mcmc.dat for further explanations.

|--------
|STELLAR PARAMETERS VALUES + PRIOR OPTIONS
|--------
If any of the following parameter is set to 'n', the parameter's values at that
 line are not considered
MASS_+_ERROR_(Msun)_+[n-y-p]	
	'y' is the obliged choice if COMPUTE_MASS_FROM_ISOCHRONES is set to 'y': at
	 each step, the jump mass value is driven by the isochronal mass thanks to
	 the Bayesian penalty that enters the merit function
	 'y' option is only compatible with COMPUTE_MASS_FROM_ISOCHRONES='y'
	'p' to be set like this if there is no way to constrain the stellar mass;
	 at each step, the jump mass value is driven by the mass prior by entering
	 the merit function
RADIUS_+_ERROR_(Rsun)_+[n-y-p]	
	'n' is the default choice if the transit is deep, so that rhos can be
	 inferred observationally. But if DERIVE_Ms_FROM_RHOs+Rs='y', Rs should be
	 set as 'p', while Ms should be set as 'n'
	'y' is the default choice if COMPUTE_MASS_FROM_ISOCHRONES='y' and the
	 transit is shallow (rhos cannot be retrieved from transit). Thus, also Rs
	 is a jump parameter and at each step it is driven by the isochronal radius
	 thanks to the Bayesian penalty that enters the merit function. If no
	 gravity proxy nor luminosity are provided, Rs is also part of the input of
	 the Isochrone placement
	 'y' option is only compatible with COMPUTE_MASS_FROM_ISOCHRONES='y'
	'p' is the default choice if COMPUTE_MASS_FROM_ISOCHRONES='n' and the transit
	 is shallow (rhos cannot be retrieved from transit); at each step, the jump
	 radius value is driven by the radius prior. 'p' option may be set also in
	 case COMPUTE_MASS_FROM_ISOCHRONES='y': in this case, Rs is also part of
	 the input of the Isochrone placement.
RHO_+_ERROR_(rho_sun)_+[n-p-y]	
	'y' uses the specified RHO value as input of the Isochrone placement (this
	 is in addition to rhos coming from transit, if it is possible to retrieve
	 it)
	'p' the jump rho value is driven by its prior value by entering the merit
	 function. 
	 Moreover, if COMPUTE_MASS_FROM_ISOCHRONES='y', rho is also part of the
	 input of the Isochrone placement (and this is in addition to rhos coming
	 from transit, if it is possible to retrieve it).
	 In case of a transit-only fit of a Jupiter-sized planet whose transit is
	 very well-defined, it's worth to try to constrain the eccentricity through
	 the photo-eccentric effect (Dawson & Johnson, 2012), by setting RHO as 'p'
	 and EXCENTRICITY as 'y' (see BG20 for an extensive discussion).
LUMINOSITY_+_ERROR_(Lsun)_+[n-p] 
	'n' is the obliged choice if COMPUTE_MASS_FROM_ISOCHRONES='y' (isoch='y').
	 To set Ls as prior if isoch='y', it is possible to provide (mag, distance)
	 as priors, then Ls will be computed within the Isochrone placement and it
	 will drive the jump luminosity value by entering the merit function.
	 Or, instead, if Ls is directly available, another possibility is computing
	 Rs from (Ls, Teff) and set Rs as prior in the input.
	 Do not set mutual dependent parameters as priors
	'p' (optionally and only if COMPUTE_MASS_FROM_ISOCHRONES='n'): the jump
	 luminosity is driven by its prior value by entering the merit function
LOGG_+_ERROR_+[n-p-y]	
	'y' logg is part of the input of the Isochrone placement
	'p' the jump logg value is driven by its prior value by entering the merit
	 function. Moreover, if COMPUTE_MASS_FROM_ISOCHRONES='y', logg is also part
	 of the input of the Isochrone placement
TEFF_+_ERROR_(K)_+[n-p]	
	'p' is always the obliged choice, unless COMPUTE_MASS_FROM_ISOCHRONES='y'
	 and COLOR is available and set to 'p'
[FEH]_+_ERROR_+[p-y] 
	'p' is the default choice
	'y' may be set in case COMPUTE_MASS_FROM_ISOCHRONES='y' when one wants to
	 load all the possible metallic grids of theoretical models. In this case
	 strong priors on other stellar parameters are expected and final
	 convergence is not guaranteed.
VTURB_+_ERROR_(kms-1)	
	Values that are considered for the limb-darkening computation
VSINI_+_ERROR_(kms-1)_+[n-p-y]	
	'y' vsini is part of the input of the Isochrone placement
	'p' the jump vsini value is driven by its prior value by entering the merit
	 function. Moreover, if COMPUTE_MASS_FROM_ISOCHRONES='y', the vsini is also
	 part of the input of the Isochrone placement
PROT_+_ERROR_(days)_+[n-y]	
	Even if set to 'n', it is considered if DEDUCE_STELLAR_INCLINATION='y' and
	 it gives stellar spin axis angle (combined with vsini).
	Moreover, if set to 'y' and COMPUTE_MASS_FROM_ISOCHRONES='y', it is part of
	 the input of the Isochrone placement
F2_+_ERROR_+[y-n-p]	
	It plays a role only if MODEL_RV_TIDE='y' and radial velocity time-series
	 are available
	'y' F2 is simply a jump parameter
	'p' F2 jump value is also driven by its prior value by entering the merit
	 function
LOGRHK_+_ERROR_+_[n-y] 
	'y' logR'HK is part of the input of the Isochrone placement (if isoch='y')
[YMg]_+_ERROR_+_[n-y]	
	'y' [Y/Mg] is part of the input of the Isochrone placement (if isoch='y')
NuMax_+_ERROR_(uHz)_+_[n-y]		
	'y' the asteroseismic nuMax is part of the input of the Isochrone placement
	 (if isoch='y') giving logg through the scaling relation by Chaplin+ 2014
DeltaNu_+_ERROR_(uHz)_+_[n-y]	
	'y' the asteroseismic deltaNu is part of the input of the Isochrone 
	 placement (if isoch='y') giving rho through the scaling relation by
	 Chaplin+ 2014
APPARENT_MAGNITUDE_+_[n-p-y], DISTANCE_(pc)_+_[n-p-y] 
	'y', 'y' They are part of the input of the Isochrone placement (if
	 isoch='y') giving stellar luminosity Ls
	'p', 'p' Like in the previous scenario, but in addition Ls is considered as
	 prior and it enters the merit function
	Only Johnson Vmag and Gaia Gmag have been tested so far
COLOR_INDEX_+_ERR_+_ID_+_[n-p]	
	'n' if Teff is set as 'p'
	'p' if Teff is set as 'n' and isoch='y', so that Teff is inferred from the
	 colour index
	ID is an integer identifier: 1 for Johnson B-V colour index; 2 for Gaia
	 G_BP-G_RP colour index

|--------
|JUMP PARAMETERS: INITIAL VALUES + ERRORS + FIT OPTIONS[y,n,p,j for P]
|--------
If any of the following parameter is set as 'n', it is not considered a jump
 parameter
If any of the following parameter is set as 'p', in addition to what is
 described, at each step the jump value is driven by its corresponding prior
 value by entering the merit function
DF_+_ERROR	
	'y' is the default choice if a lightcurve is available
TIDAL_ELONGATION_+_ERROR	
	Differs from 'n' only if interested in the treatment of tidal elongation
IMPACT_PARAMETER_+_ERROR	
	'y' is the default choice if a lightcurve is available
ECLIPSE_DURATION_+_ERROR_(day)	
	'n' if the transit is shallow and/or noisy: rhos will be inferred from
	 stellar properties and the eclipse duration W will be deduced a posteriori
	'y' if the transit is well defined: rhos will be inferred from transit
T0_+_ERROR_(BJD) 
	'y' is the default choice if a lightcurve is available
PERIOD_+_ERROR_(day)	
	'y' is the default choice if more than one transit is available, otherwise
	 set it to 'n' with the nominal value found in the literature, as there's
	 no way to constrain it
	'j' the stepping probability is weighted by the factor P(link-1)/P(link)
EXCENTRICITY_+_ERROR			
	0. 0. 'n' is the deafault choice if only transit LCs are available. If
	 interested in the photo-eccentric effect (Dawson & Johnson, 2012), it is
	 possible to try e.g. 0.1 0.1 'y' (refer also to RHO_+_ERROR_(rho_sun)
	 description and to BG20)
	If also RV are available, it is worth to try also e.g. 0.1 0.1 'y' and
	 check if the BIC is sensibly lower than the one coming from the circular
	 orbit case
OMEGA_+_ERROR_(degree)			
	0. 0. is the default choice in case circular orbit is assumed
ERROR_SQRT(E)COSW+SQRT(E)SINW	
	If EXCENTRICITY e is set as 'y' or 'p', actual jump parameters are
	 sqrt(e)*cos(w) and sqrt(e)*sin(w) (where w indicates the argument of
	 periastron OMEGA) and their errors have to be specified at this line. If
	 set to 0. 0., then sqrt(EXCENTRICITY_ERROR) are assumed as errors
K_+_ERROR_(ms-1)	
	'y' is the default choice if an RV time-series is available
	'n' if no RV time-series is used in the current analysis. But set anyway
	 the K value and its error, if K is available from the literature. In this
	 case it is still possible to compute the planetary mass
	'j' the stepping probability is weighted by the factor
	 (K(link-1)+0.1)/(K(link)+0.1)
BETA_+_ERROR_(degree)	
	Only for the first exoplanet in the list. 
	It plays a role in case the Rossiter-McLaughlin effect is observed.
	'y' or 'p': actual jump parameters are sqrt(vsini)*cos(BETA) and
	 sqrt(vsini)*sin(BETA), otherwise the only vsini jumps
ERROR_SQRT(VSINI)COSB+SINB		
	In case sqrt(vsini)*cos(BETA), sqrt(vsini)*sin(BETA) are the actual jump
	parameters, their errors have to be specified at this line. If set to 0. 0.
	then error propagation involving vsini and BETA is performed to retrieve
	the error on sqrt(vsini)*cos(BETA) and sqrt(vsini)*sin(BETA)
If Np greater than 1, further blocks of parameters from DF_+_ERROR to
 K_+_ERROR_(ms-1) must be appended (the parameters' lines BETA_+_ERROR_(degree)
 and ERROR_SQRT(VSINI)COSB+SINB must appear only for the first planet)

|--------
|FILTERS: OCC_DEPTH+ERROR, DIL+ERROR, PH_A_(LEB)_O+ERRORS
|OCC+LD+PH_FIT_OPTION[y-n-p] (+ LD VALUES/ERRORS)
|--------
For each kind of LC filter, a line is present. It is composed of 16 elements (+
 further optional elements specifiying the LD coefficients with their errors,
 if the LC filter is not among the available ones). It contains the following
 elements:
1) Filter
2) occultation depth DFsec
3) error on DFsec
4) photometric dilution (Wang+ 2014)
5) error on photometric dilution
6) amplitude of the reflection effect A_refl (For+ 2010)
7) error on A_refl
8) amplitude of the ellipsoidal effect A_elli (Mazeh 2008)
9) error on A_elli
10) amplitude of the relativistic beaming (aka Doppler boosting) effect A_beam
    (Maxted+ 2000)
11) error on A_beam
12) phase offset (in degrees)
13) error on phase offset (in degrees)
14) fit occultation		
	'p' DFsec jumps and also enters the merit function
	'n' DFsec is not a jump parameter
15) fit limb darkening	
	'p' recommended choice so to consider LD in the merit function and avoid
	 fixing LD coefficients
	'n' LD coeffcients are not jump parameters
16) fit phase	
	'p' A_refl, A_elli, A_beam jump and enter the merit function
	'n' A_refl, A_elli, A_beam are not jump parameters
(optionally): pairs of (LD coefficient, error_LD)	
	2 pairs if LIMB-DARKENING_LAW='qd'
	4 pairs if LIMB-DARKENING_LAW='nl'
	not present if LIMB-DARKENING_LAW='no'
												

|-----------------------------------------------------------------------------|
| 4) Output files                                                             |
|-----------------------------------------------------------------------------|
- mcmc_med.res It synthesizes the median values of the parameters of interest
   as inferred from their respective PDFs. It's the reference output file
- mcmc_bf.res  It synthesizes the best-fit values of the parameters, which
   have been used to build the modelled LCs
- phot____.res  It refers to the LC, whose number is specified in ____
   Relevant columns are:
   1. fobjd: folded BJD (modulo 1-orbital-period, zero-point=T0); transformed
    in TDB if UT='y' in the mcmc.dat
   2. bjd: BJD as reported in the phot____.txt, but transformed in TDB if
    UT='y' in the mcmc.dat
   3. phot: original photometric observations as in the phot____.txt
   4. photcor: corrected phot after baseline (and, if set, BM) application
   5. error: photometric errors multiplied by CF
   6. resi: residuals (observations - model, i.e. O-C)
   7. model: model-LC without baseline application
   8. model_tr: model-LC with baseline application
   9. column of zeroes
- phob____.res  Same as phot____.res, but for binned data; bin size is
   according to the BIN_SIZE_FOR_LCs_(min) parameter.
   Here follows the column correspondence phob <--> phot:
   1. fobjdbin <--> 1.
   2. bjdbin <--> 2.
   3. photbin <--> 3.
   4. error on photbin
   5. photcorbin <--> 4.
   6. error on photcorbin
   7. resibin <--> 6.
   8. error on resibin
   9. modelbin <--> 7.
   10. model_trbin <--> 8.
- mcmc_phot.res  Analogue to the phot____.res file, but it contains data
   related to all the phot____.txt files. Here follows the column
   correspondence mcmc_phot <--> phot:
   1. line counter
   2. fbjd <--> 1. and sorted in ascending order
   3. fphot <--> 3.
   4. fphotcor <--> 4.
   5. ferror <--> 5.
   6. fresi <--> 6.
   7. fmodel_tr <--> 8.
   8. column of zeroes
- mcmc_phot__.res  Same as mcmc_phot.res, but it contains only data referring
   to the filter that is specified by __
- mcmc_phob__.res  Same as mcmc_phot__.res, but for binned data; bin size is
   according to the BIN_SIZE_FOR_LCs_(min) parameter.
- rv____.res  It refers to the RV curve, whose number is specified in ____
   RVs are in m/s. Relevant columns are:
   1. line counter
   2. rvbjd: BJD as from rv____.txt
   3. rv2: RV signal after RV-baseline application
   4. rverror: RV errors as from rv____.txt, but summed in quadrature with the
    Jit jitter term as specified in the mcmc.dat
   5. rvresi: residuals (observations - model, i.e. O-C)
   6. rvmod2: global RV model
   7. column of zeroes
- forv____p__.res  Radial velocity signal caused by planet 'p__', as extracted
   from rv____.txt. RVs are in m/s. Relevant columns are:
   1. line counter as in rv____.res
   2. forvbjd: folded BJD (modulo 1-orbital-period, zero-point=T0)
   3. rv2pla: RV signal caused by planet 'p__', after RV-baseline application
   4. rverror: same as in rv____.res
   5. rvresi: same as in rv____.res
   6. rvmodpla: model-RV signal caused by planet 'p__'
   7. column of zeroes
- mcmc_rv.res  Analogue to the rv____.res file, but it contains data
   related to all the rv____.txt files. Columns have the same meaning as in
   each of the rv____.res files
- rms____.res  rms of the O-C residuals (referred to the LC specified by ____)
   that are computed by binning data on different time-scales. For a given
   temporal scale, binning is performed by slicing the LC so to have the same
   number of points in each bin. If gaps are present, it's like ignoring them.
   Relevant columns are:
   1. binning (min): 50 values by default; they are 49 values of an arithmetic
    progression starting from RED_NOISE_TIMESCALE_MIN with common difference
    (MAX-MIN)/49, besides the BIN_SIZE_FOR_LCs value.
   2. unbinned rms [ppm] of all the O-C values divided by the square root of
    the number of points contained in that given bin
   3. binned rms [ppm]: it's computed by averaging the O-C values in each bin
    and considering as mean value the mean of the averaged values
   4. beta: it's the ratio between column 3. and 2.
- mcmc_chi2.res  Evolution of the merit function
   1. link: step of the chain(s), included burn-in phase
   2. merit: global merit function
   3. photmerit: photometric merit function
   4. rvmerit: RV merit function
   
Filenames containing the name of a stellar/planetary parameter may have two
 different extensions:
- *.res: these files list the parameter's values (column 2) that have been 
   assumed step by step (column 1), after the burn-in rejection. If more than
   one chain have been considered, the stepping values coming from each chain
   are appended sequentially
- *.sm: superMongo routines producing plots of the specified parameter
   *_h.sm refers to histogram plots
   Inside the superMongo environment, after opening a graphical device, e.g.
   : dev x11 -bg white
   Read the macro MACRO.sm and then type 'plot' to see the corresponding plot
   : macro read MACRO.sm
   : plot
- Isoch folder: it contains input parameters that entered the very first step
   of the Isochrone placement and the corresponding output, just for diagnostic
   purpose
 
 
|-----------------------------------------------------------------------------|
| 5) Examples                                                                 |
|-----------------------------------------------------------------------------|
- A filled-in mcmc.dat referring to the analysis of HD 219134 system
- 4 LCs of the two transiting planets HD 219134 b, c
   phot0001.txt; phot0002.txt; phot0003.txt; phot0004.txt
- 1 RV time-series
   rv0001.txt
   

|-----------------------------------------------------------------------------|
| 6) Last updates                                                             |
|-----------------------------------------------------------------------------|
- 06/02/2020. LD coefficients for both quadratic and nonlinear models in TESS
   (TE) and Cheops (Ch) TESS (TE) bandpasses now available. Coefficients have
   been computed from ATLAS models, using the code by Espinoza&Jordan (2015;
   http://arxiv.org/abs/1503.07020)
   Please update quadratic.dat, nonlinear.dat, and mcmcI.f90 files
