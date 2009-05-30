#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <float.h>

#include "metric.h"



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Various physics and model setup parameters that are macros either for performance reasons or since no need to change them at runtime.
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// whether doing super long doubles
#define SUPERLONGDOUBLE 0

#define PRODUCTION 0
// 0: full images, dumps, etc., few more key failure stderr/fail_file messages
// 1: only log density image since too many images (takes alot of time), no utoprim failure messages -- assume debug.out and debug???? will have info if needed

// 0: normal computational zones outputted on diagnostics
// 1: 1 extra boundary zone
// 2: 2 extra boundary zones
// 3: 3 extra boundary zones
#define FULLOUTPUT 0

#define MAILWHENDONE 1
#define MAILFROMREMOTE 0
#define REMOTEHOST "bh.astro.uiuc.edu"
#define EMAILADDRESS "jmckinney@cfa.harvard.edu"

// long double constants
# define M_El           2.7182818284590452353602874713526625L  /* e */
# define M_LOG2El       1.4426950408889634073599246810018922L  /* log_2 e */
# define M_LOG10El      0.4342944819032518276511289189166051L  /* log_10 e */
# define M_LN2l         0.6931471805599453094172321214581766L  /* log_e 2 */
# define M_LN10l        2.3025850929940456840179914546843642L  /* log_e 10 */
# define M_PIl          3.1415926535897932384626433832795029L  /* pi */
# define M_PI_2l        1.5707963267948966192313216916397514L  /* pi/2 */
# define M_PI_4l        0.7853981633974483096156608458198757L  /* pi/4 */
# define M_1_PIl        0.3183098861837906715377675267450287L  /* 1/pi */
# define M_2_PIl        0.6366197723675813430755350534900574L  /* 2/pi */
# define M_2_SQRTPIl    1.1283791670955125738961589031215452L  /* 2/sqrt(pi) */
# define M_SQRT2l       1.4142135623730950488016887242096981L  /* sqrt(2) */
# define M_SQRT1_2l     0.7071067811865475244008443621048490L  /* 1/sqrt(2) */
# define SIXTH          0.1666666666666666666666666666666666L  /* 1/6 */
# define FOURTHIRD      1.3333333333333333333333333333333333L  /* 4/3 */
# define ONE            1.0000000000000000000000000000000000L
# define PTFIVE         0.5L
# define TWO            2.0L
# define ONEPT25        1.25L
# define THREE          3.0L
# define SIX            6.0L
# define EIGHT          8.0L

#define PERFTEST 0
// 0: don't perform performance test
// 1: do

#define DOAVG 0
// 0: don't do time average dumps, so don't allocate memory
// 1: do

// only setup for full Pi theta grid, not Pi/2 or any other theta grid.
// e.g. will fail for defcoord=3 Pi/2 grid.  Causes code to crash
// GODMARK
#define DOJETDIAG 1
// 0: don't do jet diagnostics (ener, not flener that is always done)
// 1: do
#define NUMJETS 2
#define INNERJET 0
#define OUTERJET 1

// -------------> r
// |	     3    
// |	    1-0   
// |	     2    
// v	        
// theta      

#define X1UP	0
#define X1DN	1
#define X2UP	2
#define X2DN	3


#if(DOJETDIAG)
#define NUMENERREGIONS 3
#else
#define NUMENERREGIONS 1
#endif

#define GLOBALENERREGION 0
#define INNERJETREGION 1
#define OUTERJETREGION 2

#define DOAVG2 0 // only make 1 if DOAVG 1 above
// 0: don't split AVG file
// 1: do

#define DODEBUG 1
// 0: don't output debug dump file or ener file(ener is based on dump counts)
// 1: do
#define DOSUPERDEBUG 0
// 0: don't output super debug output
// 1: do

// see diag_source()
#define DOLUMVSR 1

#define DOFIELDLINE 1
// 0: don't output energy@infinity and field line stuff
// 1: do

#define COMPDIM 2


/* size of global arrays */
#define N1	64		/* number of zones */
#define N2	64       	/* number of zones */
#define N3      1


// 3 maximum boundary zones needed if doing Parabolic interpolation#define N1BND 2
#define N1BND 3 
#define N2BND 3

#if(COMPDIM==3)
#define N3BND 2
#else
#define N3BND 0
#endif





#define MAXCPUS 1500
#define MAXFILENAME 200


#define PRECISEINVERSION 1
// whether we use ultimately precise inversion or "workable" inversion precision (see utoprim's for how defined)

#define VEL4 0
#define VEL3 1
#define VELREL4 2


#define WHICHVEL VELREL4
// which velocity to compute in (init can be anything (see init.c's for transformations)
// 0: 4-velocity (leads to ambiguous u^t +- discr part)
// 1: 3-velocity (unambiguous u^t but interpolation is not constrained to be a good 3-velocity)
// 2: relative 4-velocity (unambiguous u^t and any interpolation gives good value)

#define WITHGDET 0
#define WITHNOGDET 1
#define WITHSINSQ 2

#define WHICHEOM WITHGDET
//#define WHICHEOM WITHNOGDET

// if WHICHEOM==WITHNOGDET, then below determines which EOMs get what geometric prefactor.  Notice (as described in phys.c's source_conn() ) that geometry issue applies AFTER additions/subtractions of EOMs (as done by REMOVERESTMASSFROMUU).
// This stays naturally, simply consistent with how code evolves conserved quantities.
#define NOGDETRHO 0
#define NOGDETU0 0
#define NOGDETU1 0
#define NOGDETU2 0
#define NOGDETU3 0
#define NOGDETB1 0
#define NOGDETB2 0
#define NOGDETB3 0
#define NOGDETENTROPY 0



#define REMOVERESTMASSFROMUU 1
// whether to subtract rest-mass from energy equation for divT=0 equation of motion
// 0: use MHD stress tensor with explicit rest-mass included
// 1: as 0, but subtract out rest-mass from conservation and flux terms for evolution
// 2: use MHD stress tensor withOUT rest-mass
// this changes mhd_calc() in phys.c and assumes rest of code uses mhd stress-energy tensor without restmass also!!
// DIAGNOSTICS also without rest-mass for UU terms.



#define RELEOM 0
#define NONRELEOM 1 // NOT FINISHED // NOT RIGHT
// whether relativistic or nonrelativistic EOMs (speed of light limitation)
#define RELTYPE RELEOM

#define EOMGRMHD 0
#define EOMFFDE 1

#define EOMTYPE EOMFFDE
// 0 = GRMHD
// 1 = FF(D)E force-free electrodynamics
// for force-free, must turn off:
// ok now, but effectively setup already the below 2 lines implicitly
// global.h : FIXUPAFTERINIT, FIXUPAFTERRESTART,CHECKSOLUTION,LIMADJUST,FLUXADJUST
// step_ch.h : FIXUPZONES->FIXUPNOZONES


// whether to try other methods for the inversion if primary choices fails
// created because utoprim_2d_final fails for large b^2/rho when other methods (even utoprim_2d) do not fail.
#define UTOPRIMTRYAGAIN 1

// mnenomics
#define DONOENTROPY  0
#define DOEVOLVECOMPAREENTROPY 1
#define DOEVOLVEDIRECTENTROPY 2

// whether to coevolve the entropy to check for shock+reconnection dissipation
#define DOENTROPY DONOENTROPY // normal total energy equation
//#define DOENTROPY DOEVOLVECOMPAREENTROPY // coevolve and compare
//#define DOENTROPY DOEVOLVEDIRECTENTROPY // directly evolve entropy equation instead of total energy equation


#define EVOLVENOENTROPY 0
#define EVOLVESIMPLEENTROPY 1 // should be used with DOENTROPY==DOEVOLVECOMPAREENTOPY
#define EVOLVEFULLENTROPY 2 // should only be used with DOENTROPY==DOEVOLVEDIRECTENTROPY

// generically which type of entropy evolution to do, except not used when doing DOENTROPY==DOEVOLVEDIRECTENTROPY
#define WHICHENTROPYEVOLVE EVOLVESIMPLEENTROPY


// defines how Utoprimgen is used
#define EVOLVEUTOPRIM 0
#define OTHERUTOPRIM 1

// defines data return types for primtoU() and primtoflux()
#define UEVOLVE 0
#define UDIAG 1
#define UNOTHING 2
#define UENTROPY 3 // implicit UNOTHING but cons UU is overwritten by cons entropy (see primtoflux() in phys.c as used by utoprim() in utoprim.c)

// whether to call fixup() after initialization
#define FIXUPAFTERINIT 1

// whether to call fixup() after restart
#define FIXUPAFTERRESTART 1

// whether to check if rho<=0 or u<=0 when restarting.
#if(EOMTYPE==EOMGRMHD)
#define CHECKRHONEGZERORESTART 1
#else
#define CHECKRHONEGZERORESTART 0
#endif

// checks if solution is reasonble and fails a point if not
// only checks if b^2/\rho>>BSQORHOMAX, etc.
#define CHECKSOLUTION 1



// factor by which zone-zone values can be different for \gamma and internal energy, used of CHECKSOLUTION==1
#define GAMMAPERCDIFFMAX 2.0
#define UPERCDIFFMAX 10.0


// whether to limit gamma specially inside ergosphere
#define GAMMAERGOLIMIT 0

#define LIMITERFIXED 0
#define LIMITERBSQORHO 1
#define LIMITERBSQOU 2
#define LIMITERBSQORHOANDU 3

// whether to control the limiter
#define LIMADJUST LIMITERFIXED
//#define LIMADJUST LIMITERBSQORHOANDU
// 0: use fixed limiter
// 1: use limiter based upon b^2/rho
// 2: use limiter based upon b^2/u
// 3: use limiter based upon both b^2/rho or b^2/u

// whether to limit all variables or just hydro variables
#define HYDROLIMADJUSTONLY 1

#define FLUXFIXED 0 // (see get_bsqflags() in fixup.c)
#define FLUXBSQORHO 1
#define FLUXBSQOU 2
#define FLUXBSQORHOANDU 3

// determine how flux is computed per zone
#define FLUXADJUST FLUXFIXED
//#define FLUXADJUST FLUXBSQORHOANDU

// whether to change flux calculation for all variables or just hydro variables
#define HYDROFLUXADJUSTONLY 1


#define UTOPRIMRETURNNOTADJUSTED 0
#define UTOPRIMRETURNADJUSTED 1
// 0: return non-updated quantities
// 1: if u<0 or rho<0, then return updated quantities, else return non-updated quantities
#define UTOPRIMFAILRETURNTYPE UTOPRIMRETURNADJUSTED


// whether to allow negative internal energy in substep
// UTOPRIMFAILRETURNTYPE==UTOPRIMRETURNADJUSTED should be set in global.h
#define STEPOVERNEGU 0
// seems to make things worse (more failures, but also worse failures)


#define UTOPRIMSTATIC 0
#define UTOPRIMAVG 1

//#define UTOPRIMADJUST UTOPRIMSTATIC
#define UTOPRIMADJUST UTOPRIMAVG
// 0=just use static solution
// 1=use average surrounding solution, and if no good surrounding solution use the normal observer velocity with static densities

#define COORDSINGFIX 1
// whether to move polar axis to a bit larger theta
// theta value where singularity is displaced to
//#define SINGSMALL (1E-3)
#define SINGSMALL (1E-20)
// Hawley uses 0.06283 (0.02Pi)

#define VOLUMEDIFF 0
// whether to use volume regularization or not

#define MINDT 1.e-20 // minimum dt

#define JONCHECKS 1 // for vel=3 extra checks
#define JONCHECKS2 1 // for if checks on function calls


#define FLOORDIAGS 1 // whether to compute floor diagnostics

#define ANALYTICCONNECTION 0 // whether to use analytic connection
// only applies to certain metric's and coordinates, see set_grid.c

#define ANALYTICSOURCE 0 // whether to use analytic source function
// very slow with gcc, some extend pgcc, not a bit problem with icc
// only applies to certain metric's and coordaintes, see phys.c

// whicher/which ENO.
#define DOENO 0
// 0: no ENO
// 1: reconstruct U
// 2: reconstruct dU



// 0: original time is on edge and spatial on edge, but spatials are different locations.  old time.
// 1: all centered in space and all time, present time
// 2: like 0, but spatially centered (i.e. old time)
#define WHICHCURRENTCALC 1

#if((WHICHCURRENTCALC==0)||(WHICHCURRENTCALC==2))
#define NUMCURRENTSLOTS 4
#elif(WHICHCURRENTCALC==1)
#define NUMCURRENTSLOTS 5
#endif



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// MNEMONICS or other things that don't often change
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


// setup for various boundary situations
// so doesn't produce differences in irrelevant directions, whether boundary zones or not
// mac(?) macros are for use in definitions of other macros since macro with args needs to be directly a function of what's hardcoded, not some "replacement" since nothing is replaced, not to be used inside code for USING a macro.

#if(N1>1)
#define im1 i-1
#define im1mac(i) i-1
#define ip1 i+1
#define ip1mac(i) i+1
#else

#define im1 i
#define im1mac(i) i
#define ip1 i
#define ip1mac(i) i

#endif

#if(N2>1)
#define jm1 j-1
#define jm1mac(j) j-1
#define jp1 j+1
#define jp1mac(j) j+1
#else
#define jm1 j
#define jm1mac(j) j
#define jp1 j
#define jp1mac(j) j
#endif

#if(N3>1)
#define km1 k-1
#define km1mac(k) k-1
#define kp1 k+1
#define kp1mac(k) k+1
#else
#define km1 k
#define km1mac(k) k
#define kp1 k
#define kp1mac(k) k
#endif

/* NBIG is bigger of N1 and N2 and N3 */
#define NBIG1 ((N1>N2) ? N1 : N2)
#define NBIG  ((NBIG1>N3) ? NBIG1 : N3)

#define NBIGBND1 ((N1BND>N2BND) ? N1BND : N2BND)
#define NBIGBND  ((NBIGBND1>N3BND) ? NBIGBND1 : N3BND)

// N?OFF and N?NOT1 are a bit redundant
#define N1OFF (((N1BND>0)&&(N1>1)) ? 1 : 0)
#define N2OFF (((N2BND>0)&&(N2>1)) ? 1 : 0)
#define N3OFF (((N3BND>0)&&(N3>1)) ? 1 : 0)

#define N1NOT1 ((N1>1) ? 1 : 0)
#define N2NOT1 ((N2>1) ? 1 : 0)
#define N3NOT1 ((N3>1) ? 1 : 0)

/* allocated memory uses this for active zones 0-N1-1 and bc beyond that */
#define N1M (N1+N1BND*2)
#define N2M (N2+N2BND*2)
#define N3M (N3+N3BND*2)

/* NBIGM is bigger of N1M and N2M and N3M */
#define NBIG1M ((N1M>N2M) ? N1M : N2M)
#define NBIGM  ((NBIG1M>N3M) ? NBIG1M : N3M)

// maximal surface of boundary exchange
#if(COMPDIM==3)
#define NBIGS1M ((N1M*N2M>N1M*N3M) ? N1M*N2M : N1M*N3M)
#define NBIGSM ((NBIGS1M>N2M*N3M) ? NBIGS1M : N2M*N3M)
#else
#define NBIGSM (NBIGM)
#endif


// restrict loops only over relevant domain in reduced dimension case
#if(N1>1)
#define INFULL1 -N1BND
#define OUTFULL1 N1+N1BND
#define INHALF1 -N1BND/2
#define OUTHALF1 N1+N1BND/2
#define INP11 -N1BND+1
#define OUTP11 N1+N1BND-1
#else
#define INFULL1 0
#define OUTFULL1 N1
#define INHALF1 0
#define OUTHALF1 N1
#define INP11 0
#define OUTP11 N1
#endif

#if(N2>1)
#define INFULL2 -N2BND
#define OUTFULL2 N2+N2BND
#define INHALF2 -N2BND/2
#define OUTHALF2 N2+N2BND/2
#define INP12 -N2BND+1
#define OUTP12 N2+N2BND-1
#else
#define INFULL2 0
#define OUTFULL2 N2
#define INHALF2 0
#define OUTHALF2 N2
#define INP12 0
#define OUTP12 N2
#endif

#if(N3>1)
#define INFULL3 -N3BND
#define OUTFULL3 N3+N3BND
#define INHALF3 -N3BND/2
#define OUTHALF3 N3+N3BND/2
#define INP13 -N3BND+1
#define OUTP13 N3+N3BND-1
#else
#define INFULL3 0
#define OUTFULL3 N3
#define INHALF3 0
#define OUTHALF3 N3
#define INP13 0
#define OUTP13 N3
#endif










// this is the total stencil half width, where sts is full stencil size as: (sts-1)/2 which represents the one-sided safetey size in number of zones
// 2 even for parabolic, due to magnetic field stencil!
#define SAFESIZE (2)


#define MAX(a,b) ( ((a) > (b)) ? (a) : (b) )
#define MIN(a,b) ( ((a) < (b)) ? (a) : (b) )
#define SIGN(a) ( ((a) <0.) ? -1. : 1. )

#define PROGRADERISCO 0
#define RETROGRADERISCO 1

#define NUMTSCALES 4
// number of times scales to watch failure rates at
#define ALLTS 0 // full cumulative
#define ENERTS 1 // cumulative each dump_ener (over all grid)
#define IMAGETS 2 // cumulative each image dump (full grid)
#define DEBUGTS 3 // debug dump time scale (full grid)

// dump.c's fieldlinedump()
#define NUMFIELDLINEQUANTITIES 11
// rho, u, -hu_t, -T^t_t/U0, u^t, v1,v2,v3,B1,B2,B3

// see failfloorcount counter
#define NUMFAILFLOORFLAGS 9
//  mnemonics
#define COUNTUTOPRIMFAILCONV 0 // if failed to converge
#define COUNTFLOORACT 1 // if floor activated
#define COUNTLIMITGAMMAACT 2 // if Gamma limiter activated
#define COUNTINFLOWACT 3 // if inflow check activated
#define COUNTUTOPRIMFAILRHONEG 4
#define COUNTUTOPRIMFAILUNEG 5
#define COUNTUTOPRIMFAILRHOUNEG 6
#define COUNTGAMMAPERC 7 // see fixup_checksolution()
#define COUNTUPERC 8 // see fixup_checksolution()

// failure codes for utoprim failures
#define UTOPRIMFAILFIXED -1
#define UTOPRIMNOFAIL 0
#define UTOPRIMFAILTYPES 8
#define UTOPRIMFAILCONV 1
#define UTOPRIMFAILCONVUTSQ  2
#define UTOPRIMFAILCONVRET   101
#define UTOPRIMFAILCONVW     3
#define UTOPRIMFAILCONVUTSQ2 4
#define UTOPRIMFAILRHONEG 5
#define UTOPRIMFAILUNEG 6
#define UTOPRIMFAILRHOUNEG 7
#define UTOPRIMFAILGAMMAPERC 8
#define UTOPRIMFAILUPERC 9

/* failure modes */
#define FAIL_UTOPRIM_NEG	1
#define FAILSTR01 "UTOPRIM_NEG"
#define FAIL_UTOPRIM_TEST	2
#define FAILSTR02 "UTOPRIM_TEST"
#define FAIL_VCHAR_DISCR	3
#define FAILSTR03 "VCHAR_DISCR"
#define FAIL_COEFF_NEG		4
#define FAILSTR04 "COEFF_NEG"
#define FAIL_COEFF_SUP		5
#define FAILSTR05 "COEFF_SUP"
#define FAIL_UTCALC_DISCR	6
#define FAILSTR06 "UTCALC_DISCR"
#define FAIL_LDZ	        7
#define FAILSTR07 "FAIL_LDZ"
#define FAIL_BCFIX	        8
#define FAILSTR08 "FAIL_BCFIX"
#define FAIL_VSQ_NEG	        9
#define FAILSTR09 "FAIL_VSQ_NEG"



// flag failures/problems for correction/check in fixup
#define NUMPFLAGS (5)
// the below needs to be bounded since one CPU doesn't know if the other failed, and neighbor failure determines nature of how failure is treated
// also, bounded values at real boundaries need to identify if bad copy
#define FLAGUTOPRIMFAIL 0 // changes behavior of fixup()
// the below flags are done after bound_prim, and can be determined at any time, so just come after bound.
#define FLAGREALLIM 1 // value of limiter to be used
#define FLAGBSQORHO 2 // set when B^2/RHO > BSQORHOLIMIT ; currently changes  behavior of slope_lim
#define FLAGBSQOU 3 // set when B^2/u > BSQOULIMIT
#define FLAGREALFLUX 4 // type of flux to use


#define NUMSOURCES 7
// number of source terms.  Currently includes: 1) geometry, 2) radiative cooling, and 3-7) 3=total neutrino cooling (+ the 4 processes)
#define GEOMSOURCE 0
#define RADSOURCE 1
#define NEUTRINOSOURCE 2 // total of all neutrino processes
#define NEUTRINOANN 3
#define NEUTRINOPLASMA 4
#define NEUTRINOECAP 5
#define NEUTRINOBREM 6



// max number of terms in stress tensor (for term-level flux diagnostic)
#define NUMFLUXTERMS (7)



// allow for pk[0] and pk[1]
#define MAXDTSTAGES 2
// maximum number of allowed temporal integration stages


/* boundary condition mnemonics */
#define OUTFLOW	0
#define SYMM	1
#define ASYMM	2
#define FIXED	3
#define POLARAXIS 4
#define FIXEDOUTFLOW 5 // means fixed inflow but allows outflow -- basically outflow if no inflow, but if inflow then set values to fixed quantities
#define NSSURFACE 6 // whatever in bounds.c for NS surface

/* mnemonics for primitive vars; conserved vars */
#define RHO	0
#define UU	1
#define U1	2
#define U2	3
#define U3	4
#define B1	5
#define B2	6
#define B3	7
#define ENTROPY 8

/* mnemonics for dimensional indices */
#define TT	0
#define RR	1
#define TH	2
#define PH	3

/* mnemonics for centering of grid functions */
#define NUMGRIDPOS 4
#define FACE1	0
#define FACE2	1
#define CORN	2
#define CENT	3

/* mnemonics for diagnostic calls */
#define INIT_OUT	0
#define DUMP_OUT	1
#define IMAGE_OUT	1
#define LOG_OUT		1
#define FINAL_OUT	2

// size of certain dumped tavg quantities

#define NUMNORMDUMP 29 // number of "normal" dump variables
#define NUMFARADAY 6
#define NUMOTHER 1
#define NUMSTRESSTERMS (NUMFLUXTERMS*NDIM*NDIM)

/** GLOBAL ARRAY SECTION **/


/* size of global arrays */
#if(DOENTROPY==DONOENTROPY) // normal total energy equation
#define NPR	8		/* number of primitive variables */
#define NPRDUMP NPR
#define NPRINVERT 5
#define NPRBOUND NPR
#else
#define NPR	9		/* number of primitive variables */
#define NPRDUMP 8
#define NPRINVERT 5
#define NPRBOUND 8
#endif

#define NDIM	4		/* number of total dimensions.  Never
				   changes */
#define NPG	4		/* number of positions on grid for grid 
				   functions */



// size of data type used for all floats
#define FLOATTYPE 0
#define DOUBLETYPE 1
#define LONGDOUBLETYPE 2
#define LONGLONGINTTYPE 3


// define your user type here
// (normal non-sensitive or performance critical datatypes)
#define REALTYPE DOUBLETYPE
 // (non-perf critical or sensitive data types) 
#define SENSITIVE DOUBLETYPE
// WE ASSUME SENSITIVE>=REALTYPE !

#if(REALTYPE>SENSITIVE)
god=deathadflkjasdflkjasdlfkja242424
#endif

#ifndef FLT_EPSILON
#define FLT_EPSILON (1.19209290e-07F)
#endif

#ifndef DBL_EPSILON
#define DBL_EPSILON (2.2204460492503131e-16)
#endif

#ifndef LDBL_EPSILON
#define LDBL_EPSILON (1.08420217248550443401e-19L)
#endif


// need not change below datatype stuff
#if(REALTYPE==FLOATTYPE)
#define NUMEPSILON FLT_EPSILON
#define FTYPE float
#elif(REALTYPE==DOUBLETYPE)
#define NUMEPSILON DBL_EPSILON
#define FTYPE double
#elif(REALTYPE==LONGDOUBLETYPE)
#define NUMEPSILON LDBL_EPSILON
#define FTYPE long double
#endif

#if(SENSITIVE==FLOATTYPE) // for sensitive counters
#define SFTYPE float
#elif(SENSITIVE==DOUBLETYPE)
#define SFTYPE double
#elif(SENSITIVE==LONGDOUBLETYPE)
#define SFTYPE long double
#endif


// used for numerical differencing
#define NUMSQRTEPSILON (sqrt(NUMEPSILON))


// counter (integer) stuff where counts can exceed integer (2 billion)
#define COUNTTYPE DOUBLETYPE // can't make long long int work, so use double
//#define COUNTTYPE LONGLONGINTTYPE // can't make long long int work, so use double

#if(COUNTTYPE==DOUBLETYPE)
#define CTYPE double
#elif(COUNTTYPE==LONGLONGINTTYPE)
#define CTYPE long long int
#endif



#define NUMENERVAR (6+NPR+NPR+3)

/* numerical convenience */
#define BIG (1.e+30)
#define SMALL	(1.e-20)

#define SLEPSILON (1.e-6)


/* size of step in numerical derivative evaluations */
#define HSTEP	1.e-5


#define SURFACETOTAL 0
#define VOLUMETOTAL 1

#define CONSTYPE 0
#define SURFACETYPE 1
#define CUMULATIVETYPE 2
#define CONSJETINNERTYPE 3
#define CONSJETOUTERTYPE 4
#define CUMULATIVETYPE2 5

#define WRITEHEAD 0
#define READHEAD 1

#define TIMESERIESAREAMAP 0
#define FINALTDUMPAREAMAP 1

#define WRITEFILE 0
#define READFILE 1

#define ENERFNAME "ener.out"
#define GENERFNAME "gener.out"

#define NUMDUMPTYPES 8

#define IMAGECOL 0
#define RDUMPCOL 1
#define DUMPCOL 2
#define GDUMPCOL 3
#define AVGCOL 4
#define AVG2COL 5 // used when needing AVG2COL to avoid too large a file size for avgdump
#define DEBUGCOL 6
#define FIELDLINECOL 7

#if(SENSITIVE==LONGDOUBLETYPE)
// assume sensitive>=realtype in precision
#if(REALTYPE==LONGDOUBLETYPE) // was FLOATTYPE==REALTYPE and SENS=DOUBLETYPE
#define HEADERONEIN "%Lf"
#define HEADER2IN "%Lf %Lf"
#define HEADER3IN "%Lf %Lf %Lf"
#define HEADER4IN "%Lf %Lf %Lf %Lf"
#define HEADER5IN "%Lf %Lf %Lf %Lf %Lf"
#define HEADER6IN "%Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADER7IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADER8IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define HEADER9IN "%Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#define RESTARTHEADER "%d %d "\
	      "%Lf %Lf %ld %Lf %Lf %Lf "\
	      "%Lf %Lf %Lf %Lf %ld %Lf %ld %ld %ld %ld %ld "\
	      "%Lf %d %d %d %d %d %d "\
	      "%Lf %Lf %Lf %Lf %d "\
              "%d %d %d %d "\
              "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf %Lf"
#elif(REALTYPE==DOUBLETYPE)
#define HEADERONEIN "%lf"
#define HEADER2IN "%lf %lf"
#define HEADER3IN "%lf %lf %lf"
#define HEADER4IN "%lf %lf %lf %lf"
#define HEADER5IN "%lf %lf %lf %lf %lf"
#define HEADER6IN "%lf %lf %lf %lf %lf %lf"
#define HEADER7IN "%lf %lf %lf %lf %lf %lf %lf"
#define HEADER8IN "%lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER9IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define RESTARTHEADER "%d %d "\
	      "%Lf %Lf %ld %lf %lf %lf "\
	      "%Lf %Lf %Lf %Lf %ld %Lf %ld %ld %ld %ld %ld "\
	      "%Lf %d %d %d %d %d %d "\
	      "%lf %lf %lf %lf %d "\
              "%d %d %d %d "\
              "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#elif(REALTYPE==FLOATTYPE)
#define HEADERONEIN "%f"
#define HEADER2IN "%f %f"
#define HEADER3IN "%f %f %f"
#define HEADER4IN "%f %f %f %f"
#define HEADER5IN "%f %f %f %f %f"
#define HEADER6IN "%f %f %f %f %f %f"
#define HEADER7IN "%f %f %f %f %f %f %f"
#define HEADER8IN "%f %f %f %f %f %f %f %f"
#define HEADER9IN "%f %f %f %f %f %f %f %f %f"
#define RESTARTHEADER "%d %d "\
	      "%lf %lf %ld %f %f %f "\
	      "%lf %lf %lf %lf %ld %lf %ld %ld %ld %ld %ld "\
	      "%lf %d %d %d %d %d %d "\
	      "%f %f %f %f %d "\
              "%d %d %d %d "\
              "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f"
#endif

#elif(SENSITIVE==DOUBLETYPE)
// assume sensitive>=realtype in precision
#if(REALTYPE==DOUBLETYPE)
#define HEADERONEIN "%lf"
#define HEADER2IN "%lf %lf"
#define HEADER3IN "%lf %lf %lf"
#define HEADER4IN "%lf %lf %lf %lf"
#define HEADER5IN "%lf %lf %lf %lf %lf"
#define HEADER6IN "%lf %lf %lf %lf %lf %lf"
#define HEADER7IN "%lf %lf %lf %lf %lf %lf %lf"
#define HEADER8IN "%lf %lf %lf %lf %lf %lf %lf %lf"
#define HEADER9IN "%lf %lf %lf %lf %lf %lf %lf %lf %lf"
#define RESTARTHEADER "%d %d "\
	      "%lf %lf %ld %lf %lf %lf "\
	      "%lf %lf %lf %lf %ld %lf %ld %ld %ld %ld %ld "\
	      "%lf %d %d %d %d %d %d "\
	      "%lf %lf %lf %lf %d "\
              "%d %d %d %d "\
              "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf"
#elif(REALTYPE==FLOATTYPE)
#define HEADERONEIN "%f"
#define HEADER2IN "%f %f"
#define HEADER3IN "%f %f %f"
#define HEADER4IN "%f %f %f %f"
#define HEADER5IN "%f %f %f %f %f"
#define HEADER6IN "%f %f %f %f %f %f"
#define HEADER7IN "%f %f %f %f %f %f %f"
#define HEADER8IN "%f %f %f %f %f %f %f %f"
#define HEADER9IN "%f %f %f %f %f %f %f %f %f"
#define RESTARTHEADER "%d %d "\
	      "%lf %lf %ld %f %f %f "\
	      "%lf %lf %lf %lf %ld %lf %ld %ld %ld %ld %ld "\
	      "%lf %d %d %d %d %d %d "\
	      "%f %f %f %f %d "\
              "%d %d %d %d "\
              "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f"
#endif

#elif(SENSITIVE==FLOATTYPE)

#if(REALTYPE==DOUBLETYPE)
#define RESTARTHEADER "" // dumb, so crash on compile
#elif(REALTYPE==FLOATTYPE)
#define HEADERONEIN "%f"
#define HEADER2IN "%f %f"
#define HEADER3IN "%f %f %f"
#define HEADER4IN "%f %f %f %f"
#define HEADER5IN "%f %f %f %f %f"
#define HEADER6IN "%f %f %f %f %f %f"
#define HEADER7IN "%f %f %f %f %f %f %f"
#define HEADER8IN "%f %f %f %f %f %f %f %f"
#define HEADER9IN "%f %f %f %f %f %f %f %f %f"
#define RESTARTHEADER "%d %d "\
	      "%f %f %ld %f %f %f "\
	      "%f %f %f %f %ld %f %ld %ld %ld %ld %ld "\
	      "%f %d %d %d %d %d %d "\
	      "%f %f %f %f %d "\
              "%d %d %d %d "\
              "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %f %f %f %f %f %f %f %f %f %f %f"
#endif

#endif


// 2
// 6
// 10
// 3
// 5
// 4
// SUM=30
// 23+10=33
// total=63
#define WRITERESTARTHEADER "%d %d " \
		 "%21.15g %21.15g %ld %21.15g %21.15g %21.15g " \
		 "%21.15g %21.15g %21.15g %21.15g %ld %21.15g %ld %ld %ld %ld %ld" \
		 "%21.15g %d %d %d %d %d %d " \
		 "%21.15g %21.15g %21.15g %21.15g %d " \
                 "%d %d %d %d " \
                 "%ld %d %d %d %ld %ld %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g "

#define HEADERONEOUT "%21.15g "

// now that all hashes have been defined, get mpi header
#include "mympi.h"


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// LOOPS
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

// these loops used for general purposes
#define LOOPF3 for(k=INFULL3;k<OUTFULL3;k++)
#define LOOPF2 for(j=INFULL2;j<OUTFULL2;j++)
#define LOOPF1 for(i=INFULL1;i<OUTFULL1;i++)

#define LOOPH3 for(k=INHALF3;k<OUTHALF3;k++)
#define LOOPH2 for(j=INHALF2;j<OUTHALF2;j++)
#define LOOPH1 for(i=INHALF1;i<OUTHALF1;i++)

#define LOOPP13 for(k=INP13;k<OUTP13;k++)
#define LOOPP12 for(j=INP12;j<OUTP12;j++)
#define LOOPP11 for(i=INP11;i<OUTP11;i++)

#define LOOPN3 for(k=0;k<N3;k++)
#define LOOPN2 for(j=0;j<N2;j++)
#define LOOPN1 for(i=0;i<N1;i++)

#define LOOPFMHP3 for(k=INFULL3;k<OUTHALF3;k++)
#define LOOPFMHP2 for(j=INFULL2;j<OUTHALF2;j++)
#define LOOPFMHP1 for(i=INFULL1;i<OUTHALF1;i++)

#define LOOPHMFP3 for(k=INHALF3;k<OUTFULL3;k++)
#define LOOPHMFP2 for(j=INHALF2;j<OUTFULL2;j++)
#define LOOPHMFP1 for(i=INHALF1;i<OUTFULL1;i++)

#define LOOPHP3 for(k=0;k<OUTHALF3;k++)
#define LOOPHP2 for(j=0;j<OUTHALF2;j++)
#define LOOPHP1 for(i=0;i<OUTHALF1;i++)

// below used for initialization and such, not a computational issue
#define LOOPF LOOPF3 LOOPF2 LOOPF1
#define LOOPH LOOPH3 LOOPH2 LOOPH1
#define LOOPP1 LOOPP13 LOOPP12 LOOPP11
#define LOOP LOOPN3 LOOPN2 LOOPN1
#define LOOPFMHP LOOPFMHP3 LOOPFMHP2 LOOPFMHP1
#define LOOPHMFP LOOPHMFP3 LOOPHMFP2 LOOPHMFP1
#define LOOPHP LOOPHP3 LOOPHP2 LOOPHP1

#define LOOPINT3 for(k=intix3;k<intox3;k++)
#define LOOPINT2 for(j=intix2;j<intox2;j++)
#define LOOPINT1 for(i=intix1;i<intox1;i++)

#define LOOPDIAGOUTPUT(bzones3,bzones2,bzones1)\
for(k=((-bzones3<INFULL3) ? INFULL3 : -bzones3);k<((N3+bzones3>OUTFULL3) ? OUTFULL3 : N3+bzones3);k++)\
for(j=((-bzones2<INFULL2) ? INFULL2 : -bzones2);j<((N2+bzones2>OUTFULL2) ? OUTFULL2 : N2+bzones2);j++)\
for(i=((-bzones1<INFULL1) ? INFULL1 : -bzones1);i<((N1+bzones1>OUTFULL1) ? OUTFULL1 : N1+bzones1);i++)

#define LOOPLIMIT(bzones3,bzones2,bzones1) LOOPDIAGOUTPUT(-bzones3,-bzones2,-bzones1)

#define LOOPDIAGOUTPUTFULL(bzones3,bzones2,bzones1)\
for(k=((-bzones3<INFULL3) ? INFULL3 : -bzones3);k<((totalsize[3]+bzones3>OUTFULL3) ? OUTFULL3 : totalsize[3]+bzones3);k++)\
for(j=((-bzones2<INFULL2) ? INFULL2 : -bzones2);j<((totalsize[2]+bzones2>OUTFULL2) ? OUTFULL2 : totalsize[2]+bzones2);j++)\
for(i=((-bzones1<INFULL1) ? INFULL1 : -bzones1);i<((totalsize[1]+bzones1>OUTFULL1) ? OUTFULL1 : totalsize[1]+bzones1);i++)




#define LOOPFC LOOPF
#define LOOPHC LOOPH
#define LOOPFMHPC LOOPFMHP
#define LOOPHMFPC LOOPHMFP
#define LOOPHPC LOOPHP


#define LOOPC3 LOOPN3
#define LOOPC2 LOOPN2
#define LOOPC1 LOOPN1

#define LOOPC LOOPC3 LOOPC2 LOOPC1


/* loop over all active zones */
#define ZLOOP for(i=0;i<N1;i++)for(j=0;j<N2;j++)

/* specialty loop */
#define ZSLOOP(istart,istop,jstart,jstop) \
	for(i=istart;i<=istop;i++)\
	for(j=jstart;j<=jstop;j++)

#define FULLLOOP ZSLOOP(-N1BND, N1 -1 + N1BND, -N2BND, N2 -1 + N2BND)
#define PLUSLOOP ZSLOOP(-1, N1, -1, N2)
// below same as FULLLOOP if NBND=2
#define PLUSPLUSLOOP ZSLOOP(-2, N1+1, -2, N2+1)


#define GENLOOP(i,j,istart,istop,jstart,jstop) \
        for((i)=(istart);(i)<=(istop);(i)++)\
        for((j)=(jstart);(j)<=(jstop);(j)++)


/* specialty loop */
#define LOOPDIVB LOOPP11 LOOPP12
//ZSLOOP(-N1BND+1,N1+1,-1,N2+1)


/* want dump output to be ordered in radius first!! */
#define DUMPLOOP(istart,istop,jstart,jstop) \
	for(j=jstart;j<=jstop;j++)\
	for(i=istart;i<=istop;i++)


#define NUMIMAGEPARMS 2

#define ORIGIN 0
#define LMIN 1



#define IMAGELOOP(istart,istop,jstart,jstop) \
	for(j=jstart;j<=jstop;j++)\
	for(i=istart;i<=istop;i++)

#define OLDIMAGELOOP for(j=N2-1;j>=0;j--) for(i=0;i<N1;i++)	// nasty 
								// to
								// deal 
								// with

/* loop over all Primitive variables */
#define PLOOP for(k=0;k<NPR;k++)
/* loop over all dumped Primitive variables */
#define PDUMPLOOP for(k=0;k<NPRDUMP;k++)
/* loop over all inversion Primitive variables */
#define PINVERTLOOP for(k=0;k<NPRINVERT;k++)
/* loop over all bounding Primitive variables */
#define PBOUNDLOOP for(k=0;k<NPRBOUND;k++)


/* loop over all Dimensions; second rank loop */
#define DLOOP for(j=0;j<NDIM;j++)for(k=0;k<NDIM;k++)
/* loop over all Dimensions; first rank loop */
#define DLOOPA for(j=0;j<NDIM;j++)
/* loop over all Space dimensions; second rank loop */
#define SLOOP for(j=1;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* loop over all Space dimensions; first rank loop */
#define SLOOPA for(j=1;j<NDIM;j++)
/* spatial loop with k */
#define SLOOPAk for(k=1;k<NDIM;k++)
/* loop over all for j and Space for k; second rank loop */
#define DSLOOP for(j=0;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* space-space */
#define SSLOOP for(j=1;j<NDIM;j++)for(k=1;k<NDIM;k++)
/* loop over all for k and Space for j; second rank loop */
#define SDLOOP for(j=1;j<NDIM;j++)for(k=0;k<NDIM;k++)

// loop over directions
#define DIRLOOP for(dir=0;dir<COMPDIM*2;dir++)

// loop over fail flag in boundary code
#define FLOOP for(k=FLAGUTOPRIMFAIL;k<=FLAGUTOPRIMFAIL;k++)

// loop over jet regions
#define JETLOOP for(jetio=0;jetio<NUMJETS;jetio++)

// loop over ener/flux regions
#define ENERREGIONLOOP for(enerregion=0;enerregion<NUMENERREGIONS;enerregion++)

// loop over fair/floor types
#define FLOORLOOP for(floor=0;floor<NUMFAILFLOORFLAGS;floor++)

// loop over debug time scales
#define TSCALELOOP for(tscale=0;tscale<NUMTSCALES;tscale++)


// loop over sources
#define SCLOOP for(sc=0;sc<NUMSOURCES;sc++)

// loop over fluxterms
#define FLLOOP for(fl=0;fl<NUMFLUXTERMS;fl++)


// loop over pflag flags
#define PFLAGLOOP for(pf=0;pf<NUMPFLAGS;pf++)

// for USEMPI&&USEROMIO==1
#define ROMIOCOLLOOP for(romiocoliter=0;romiocoliter<romiocloopend;romiocoliter++)

#define BUFFERINIT nextbuf=0
#define COLINIT nextcol=0

// for mpicombie==0
#define COLLOOP for(coliter=0;coliter<numfiles;coliter++)


#define DTSTAGELOOP for(dtstage=0;dtstage<MAXDTSTAGES;dtstage++)



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Various macros
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



// #define DIVBCONDITION(p,i,j)
// if((i>=-1)&&(j>=-1)&&(startpos[2]+j!=0)&&(startpos[2]+j!=N2TOT))
//#define DIVBCONDITION(p,i,j) if((startpos[1]+i>0)&&(startpos[2]+j>0)&&(startpos[1]+i<totalsize[1])&&(startpos[2]+j<totalsize[2]))
//#define DIVBCONDITION(p,i,j) if((startpos[1]+i!=0)&&(startpos[1]+i!=totalsize[1])&&(startpos[2]+j!=0)&&(startpos[2]+j!=totalsize[2]))
#define DIVBCONDITION(p,i,j) if((startpos[2]+j!=0)&&(startpos[2]+j!=totalsize[2]))

#define DIVBFLUXCT(p,i,j) (MAX(dx[1],dx[2])*(0.5*(\
			      p[i][j][B1]*gdet[i][j][CENT] + p[i][j-1][B1]*gdet[i][j-1][CENT]\
			      - p[i-1][j][B1]*gdet[i-1][j][CENT] - p[i-1][j-1][B1]*gdet[i-1][j-1][CENT]\
			      )/dx[1] +\
		      0.5*(\
			      p[i][j][B2]*gdet[i][j][CENT] + p[i-1][j][B2]*gdet[i-1][j][CENT]\
			      - p[i][j-1][B2]*gdet[i][j-1][CENT] - p[i-1][j-1][B2]*gdet[i-1][j-1][CENT]\
			      )/dx[2])\
                              /(gdet[i][j][CORN]*0.125*fabs(fabs(p[i][j][B1]*gdet[i][j][CENT] + p[i][j-1][B1]*gdet[i][j-1][CENT]\
			      + p[i-1][j][B1]*gdet[i-1][j][CENT] + p[i-1][j-1][B1]*gdet[i-1][j-1][CENT]\
			      + p[i][j][B2]*gdet[i][j][CENT] + p[i-1][j][B2]*gdet[i-1][j][CENT]\
			      + p[i][j-1][B2]*gdet[i][j-1][CENT] + p[i-1][j-1][B2]*gdet[i-1][j-1][CENT]) +SMALL)))
#define DIVBFLUXCD(p,i,j)  (MAX(dx[1],dx[2])*(0.5*(\
				p[i+1][j][B1]*gdet[i+1][j][CENT] - p[i-1][j][B1]*gdet[i-1][j][CENT])/dx[1] +\
		             0.5*(\
				p[i][j+1][B2]*gdet[i][j+1][CENT] - p[i][j-1][B2]*gdet[i][j-1][CENT])/dx[2])\
                                /(gdet[i][j][CORN]*0.25*fabs(fabs(p[i+1][j][B1]*gdet[i+1][j][CENT] + p[i-1][j][B1]*gdet[i-1][j][CENT]\
                                + p[i][j+1][B2]*gdet[i][j+1][CENT] + p[i][j-1][B2]*gdet[i][j-1][CENT])+SMALL)))

// poles defined as divb=0, can't divide due to singularity (could use
// volume regularization)
#define SETFDIVBFLUXCT(divb,p,i,j) {DIVBCONDITION(p,i,j){ divb = fabs(DIVBFLUXCT(p,i,j)) ;} else divb = 0.;}

#define SETFDIVBFLUXCD(divb,p,i,j) {DIVBCONDITION(p,i,j){ divb = fabs(DIVBFLUXCD(p,i,j)) ;} else divb = 0.;}

// for now // GODMARK
#define SETFDIVB(divb,p,i,j) SETFDIVBFLUXCT(divb,p,i,j)




#define MYDMIN(a,b) (mydminarg1=(a),mydminarg2=(b),(mydminarg1) < (mydminarg2) ?\
        (mydminarg1) : (mydminarg2))

#define delta(i,j) ((i == j) ? 1. : 0.)
#define dot(a,b) (a[0]*b[0] + a[1]*b[1] + a[2]*b[2] + a[3]*b[3])

#define mink(I,J) (I != J ? (0.) : (I == 0 ? (-1.) : (1.)))

#define pfixupeach(pr,i,j,which,min) {if(pr[which]<min){ fladd[which]+=dV*gdet[i][j][CENT]*(min-pr[which]); pr[which]=min;}}

#define pfixup(pr,i,j) {pfixupeach(pr,i,j,RHO,RHOMIN); pfixupeach(pr,i,j,UU,UUMIN); }

// #define FAILSTATEMENT(file,function,number) {fprintf(fail_file,"%s
// %d-%s(): failure\n",file,number,function); fflush(fail_file);
// fprintf(fail_file,"rho[i][j]: %21.15g uu[i][j]: %21.15g rho2[i][j]:
// %21.15g uu2[i][j]: %21.15g i: %d j: %d k:
// %d\n",p[i][j][RHO],p[i][j][UU],ph[i][j][RHO],ph[i][j][UU],i,j,k);
// return(1);}

#define FAILSTATEMENT(file,function,number) {if(debugfail>=1){ dualfprintf(fail_file,"%s %d-%s(): failure\n",file,number,function); dualfprintf(fail_file,"i: %d j: %d p: %d\n",icurr,jcurr,pcurr);} return(1);}

#if(JONCHECKS2)
#define MYFUN(fun,one,two,three) if(fun>=1){ FAILSTATEMENT(one,two,three);}
#else
#define MYFUN(fun,one,two,three) {fun;}
#endif



//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// structure definitions
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


struct blink {
  int num;
  struct blink * np;
  // only used by cpu=0
  int cpu; // which cpu
  int i,j,k,col; // starting values for cpu=0
  int ri,rj,rk,rcol; // reference values for first cpu in sequence of nodes for a single buffer
  int end;
};



// structure declarations
/* set global variables that indicate current local metric, etc. */
struct of_geom {
  //FTYPE gcon[NDIM][NDIM];
  //  FTYPE gcov[NDIM][NDIM];
  // bit faster since not all values always used
  FTYPE (*gcov)[NDIM];
  FTYPE (*gcon)[NDIM];
  FTYPE g;
  FTYPE e[NPR]; // eomfunc
  int i,j,p;
};


struct of_state {
  FTYPE rho;
  FTYPE ie;
  FTYPE ucon[NDIM];
  FTYPE ucov[NDIM];
  FTYPE bcon[NDIM];
  FTYPE bcov[NDIM];
};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// function declarations
//
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


extern int main(int argc, char *argv[]);
extern int error_check(int wherefrom);
extern int find_horizon(void);

// initialize DUMP stuff
extern int init_dumps(void);
int setuplinklist(int numcolumns,int which);
extern struct blink * addlink(struct blink * clinkptr);

// ENER file stuff
extern int dump_ener(int doener, int dordump, int call_code);

extern int diag(int call_code);
extern void diag_source(struct of_geom *ptrgeom, FTYPE (*dUcomp)[NPR],SFTYPE Dt);
extern void frdotout(void);

extern void makedirs(void);

extern void appendener(FILE* ener_file,SFTYPE pdot_tot[][NPR],SFTYPE*fladd_tot,SFTYPE*sourceadd_tot);

extern void divbmaxavg(FTYPE p[][N2M][NPR],SFTYPE*ptrdivbmax,SFTYPE*ptrdivbavg);
extern void gettotal(int numvars, SFTYPE* vars[],int*sizes,SFTYPE*vars_tot[]);
extern void gettotali(int numvars, int* vars[],int*sizes,int*vars_tot[]);
extern int constotal(int enerregion, SFTYPE *vars_tot);
extern int integrate(SFTYPE * var,SFTYPE *var_tot,int type, int enerregion);

extern int counttotal(int tscale, int enerregion, CTYPE *vars_tot, int num);
extern int integratel(CTYPE * var,CTYPE *var_tot,int type, int tscale, int enerregion);


// DUMP file stuff
extern int isenoughfreespace(unsigned long long need);

extern int dump_gen(int readwrite, long dump_cnt, int bintxt, int whichdump,MPI_Datatype datatype, char *fileprefix, char *fileformat, char *filesuffix, int (*headerfun) (int bintxt, FILE*headerptr),int (*content) (int i, int j, MPI_Datatype datatype, void*setbuf));

extern int dump(long dump_cnt);
extern int dump_content(int i, int j, MPI_Datatype datatype,void *writebuf);
extern int dump_header(int bintxt, FILE *headerptr);

extern int avgdump(long avg_cnt);
extern int avg_content(int i, int j, MPI_Datatype datatype,void *writebuf);

extern int avgdump2(long avg_cnt);
extern int avg2_content(int i, int j, MPI_Datatype datatype,void *writebuf);

extern int debugdump(long debug_cnt);
extern int debug_content(int i, int j, MPI_Datatype datatype,void *writebuf);

extern int gdump(void);
extern int gdump_content(int i, int j, MPI_Datatype datatype, void *writebuf);

extern int fieldlinedump(long fieldline_cnt);
extern int fieldline_content(int i, int j, MPI_Datatype datatype,void *writebuf);


extern int image_dump(long image_cnt);
extern int imagedefs(int whichk, int scale, int limits, int vartype);
extern int image(long dump_cnt, int whichk, int scale, int limits, int vartype);
extern int image_header(int bintxt, FILE *headerptr);
extern int image_content(int i, int j, MPI_Datatype datatype,void *writebuf);
extern void prminmaxsum(FTYPE p[][N2M][NPR], int start,int nmemb, FTYPE *max, FTYPE*min,FTYPE*sum);

extern int restart_init(int which);

extern int restart_read(long which);
extern int read_restart_header(int bintxt, FILE* headerptr);
extern int restart_read_defs(void);
extern int rdump_read_content(int i, int j, MPI_Datatype datatype,void *writebuf);

extern int restart_write(long dump_cnt);
extern int write_restart_header(int bintxt, FILE* headerptr);
extern int rdump_content(int i, int j, MPI_Datatype datatype,void *writebuf);

extern void myfopen(char*fname, char*fmt, char*message, FILE ** fileptr);
extern void myfclose(FILE ** fileptr,char*message);

extern void myset(MPI_Datatype datatype, void *ptr, int start, size_t nmemb, void*writebuf);
extern void myget(MPI_Datatype datatype, void *ptr, int start, size_t nmemb, void*writebuf);

extern void myfwrite(int bintxt, MPI_Datatype datatype, void *ptr, int start, size_t nmemb, int i, int j, FILE**stream,void*writebuf);

extern void myfread(int bintxt, MPI_Datatype datatype, void *ptr, int start, size_t nmemb, int i, int j, FILE**stream,void*writebuf);


// initialize stuff

extern int init(int argc, char *argv[]);
extern int init_defgrid(void);
extern int init_defglobal(void);
extern int init_defconsts(void);
extern void set_arrays(void), set_grid(void);

extern int init_grid(void);
extern int init_global(void);

extern int init_dsandvels(int *whichvel, int *whichcoord, int i, int j, FTYPE *p);
extern int init_postfield(FTYPE pr[][N2M][NPR]);
extern int init_atmosphere(int *whichvel, int *whichcoord, int i, int j, FTYPE *pr);
extern int normalize_densities(FTYPE p[][N2M][NPR]);
extern int init_vpot(int i, int j, FTYPE *A);
extern int init_vpot2field(FTYPE A[][N2M],FTYPE pr[][N2M][NPR]);
extern int normalize_field(FTYPE p[][N2M][NPR]);
extern int vpot2field(FTYPE A[][N2M],FTYPE p[][N2M][NPR]);

extern int post_init(void);
extern int post_init_specific_init(void);
extern int pre_init_specific_init(void);
extern int pre_init(int argc, char *argv[]);


// some physics

extern int sourcephysics(FTYPE *ph, struct of_geom *geom,
		       struct of_state *q,FTYPE (*dUcomp)[NPR]);

extern int step_ch(void);
extern void postdt(void);
extern int primtoU(int returntype, FTYPE *p, struct of_state *q, struct of_geom *geom,
		   FTYPE *U);
extern int advance(int stage, FTYPE pi[][N2M][NPR], FTYPE ulast[][N2M][NPR], FTYPE pb[][N2M][NPR],
		   FTYPE *CUf, FTYPE pf[][N2M][NPR],FTYPE *Cunew,FTYPE unew[][N2M][NPR],int stagenow, int numstages, FTYPE *ndt);
extern int fluxcalc(int stage, FTYPE pr[][N2M][NPR], FTYPE F[][N2M][NPR],
		    int dir, FTYPE *ndt);
extern void flux_cd(FTYPE F1[][N2M][NPR], FTYPE F2[][N2M][NPR]);

extern int ucon_calc_3vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon);
extern int ucon_calc_rel4vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon);
extern int ucon_calc_4vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon);
#if(RELTYPE==RELEOM)

#if(WHICHVEL==VEL4)
#define ucon_calc ucon_calc_4vel
#define dudp_calc dudp_calc_gen
#elif(WHICHVEL==VEL3)
#define ucon_calc ucon_calc_3vel
#define dudp_calc dudp_calc_3vel
#elif(WHICHVEL==VELREL4)
#define ucon_calc ucon_calc_rel4vel
#define dudp_calc dudp_calc_gen

#elif(RELTYPE==NONRELEOM) // not really right
#define ucon_calc ucon_calc_nonrel
#define dudp_calc dudp_calc_nonrel
#endif

#endif

extern int ucon_calcother(FTYPE *pr, FTYPE *ucon);
extern void ucon_precalc(FTYPE *ucon, FTYPE *AA, FTYPE *BB,
			 FTYPE *CC, FTYPE *discr);



extern FTYPE ranc(int seed);

// fixup stuff

extern int check_pr(FTYPE *pr, FTYPE *prmodel, struct of_geom *geom, int modelpos,int finalstep);
extern int ucon_fix(FTYPE disc, FTYPE AA, FTYPE BB, FTYPE CC,
		    FTYPE *ucon);

/* // dudp stuff */

/* extern void dutdui_calc(FTYPE *ucon, FTYPE *dutdui); */
/* extern void duiduj_calc(FTYPE *ucon, FTYPE *dutdui); */
/* extern void dbtdui_calc(FTYPE *dutdui, FTYPE *pr, FTYPE *dbtdui); */
/* extern void dbiduj_calc(FTYPE *dbtdui, FTYPE *dutdui, FTYPE *ucon, */
/* 			FTYPE *b, FTYPE dbiduj[][NDIM]); */
/* extern void db2dui_calc(FTYPE dbiduj[][NDIM], FTYPE *b, */
/* 			FTYPE *db2dui); */
/* extern void duudud_calc(FTYPE *ucon, FTYPE duudud[][NDIM]); */

/* extern void dbsqdui_calc(FTYPE dbiduj[][NDIM], FTYPE *b, */
/* 			 FTYPE *dbsqdui); */
/* extern void dgdvi_calc(FTYPE *pr,FTYPE *dgdvi); */
/* extern void duidvj_calc(FTYPE *dgdv,FTYPE duidvj[][NDIM]); */
/* extern void dudduu_calc(FTYPE*dutdui, FTYPE dudduu[][NDIM]); */
/* extern void dbdiduj_calc(FTYPE dbiduj[][NDIM],FTYPE dbdiduj[][NDIM]); */
/* extern void ducon_dv3_calc(struct of_state *q,FTYPE ducon_dv[][NDIM]); */
extern int sp_stress_calc(FTYPE *pr, FTYPE tens_matt[][NDIM],
			  FTYPE tens_em[][NDIM], FTYPE *b,
			  FTYPE *ucon);



// log file stuff
extern void myfprintf(FILE* fileptr, char *format, ...);
extern void dualfprintf(FILE* fileptr,char *format, ...);
extern void logsfprintf(char *format, ...);
extern void trifprintf(char *format, ...);

// boundary stuff
extern int bound_prim(int boundstage, FTYPE prim[][N2M][NPR]);
extern int bound_pflag(int boundstage, int primbase[][N2M][NUMPFLAGS]);
extern int inflow_check_4vel(int dir, FTYPE *pr, struct of_geom *ptrgeom, int finalstep);
extern int inflow_check_3vel(int dir, FTYPE *pr, struct of_geom *ptrgeom, int finalstep);
extern int inflow_check_rel4vel(int dir, FTYPE *pr, struct of_geom *ptrgeom, int finalstep);


// transform stuff
extern int bl2met2metp2v(int whichvel, int whichcoord, FTYPE *pr, int ii, int jj);
extern int metp2met2bl(int whichvel, int whichcoord, FTYPE *pr, int ii, int jj);
extern int pr2ucon(int whichvel, FTYPE *pr, struct of_geom *geom, FTYPE*ucon);
extern int coordtrans(int whichcoordin, int whichcoordout, int ii, int jj, FTYPE*ucon);
extern void bltoks(int ii, int jj, FTYPE*ucon);
extern void kstobl(int ii, int jj, FTYPE*ucon);
extern void mettometp(int ii, int jj, FTYPE*ucon);
extern void metptomet(int ii, int jj, FTYPE*ucon);
extern void ucon2pr(int whichvel, FTYPE *ucon, struct of_geom *geom, FTYPE *pr);
extern int vcon2pr(int whichvel, FTYPE *vcon, struct of_geom *geom, FTYPE *pr);



// metric stuff
extern void gset(int getprim, int whichcoord, int i, int j, struct of_geom *geom);
extern FTYPE gdet_func(FTYPE gcov[][NDIM]);
//extern FTYPE bl_gdet_func(FTYPE r, FTYPE th);
//extern void bl_gcov_func(FTYPE r, FTYPE th, FTYPE gcov[][NDIM]);
//extern void bl_gcon_func(FTYPE r, FTYPE th, FTYPE gcon[][NDIM]);
extern void conn_func(int whichcoord, FTYPE *X, struct of_geom *geom,
		      FTYPE lconn[][NDIM][NDIM],FTYPE *conn2);
extern void mks_unitheta_idxvol_func(int i, int j, FTYPE *idxvol);

extern void gcov_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM]);
extern void gcon_func(int getprim, int whichcoord, FTYPE *X, FTYPE gcov[][NDIM], FTYPE gcon[][NDIM]);
extern void matrix_inverse(FTYPE gcov[][NDIM], FTYPE gcon[][NDIM]);

// coordinate stuff
extern void set_coord_parms(void);
extern void write_coord_parms(void);
extern void read_coord_parms(void);
extern void coord(int i, int j, int loc, FTYPE *X);
extern void bl_coord(FTYPE *X, FTYPE *r, FTYPE *th);
extern int setihor(void);
extern FTYPE setRin(int ihor);


// physics stuff
extern FTYPE contract(FTYPE *vcon, FTYPE *wcon);

extern int bsq_calc(FTYPE *pr, struct of_geom *geom, FTYPE *b2);
extern void b_calc(FTYPE *pr, FTYPE *ucon, FTYPE *b);


extern int gamma_calc(FTYPE *pr, struct of_geom *geom,FTYPE *gamma);


extern int dudp_calc_gen(int whichcons,FTYPE *pr, struct of_state *q,
		     struct of_geom *geom, FTYPE **alpha);

extern int dudp_calc_3vel(int whichcons,FTYPE *pr, struct of_state *q,
		     struct of_geom *geom, FTYPE **alpha);


extern void flux2U(int i, int j, FTYPE F1[][N2M][NPR],FTYPE F2[][N2M][NPR],FTYPE *dU, FTYPE *CUf, FTYPE *Cunew, FTYPE *Ui, FTYPE *ulast, FTYPE *Uf, FTYPE *unew);


extern int Utoprimgen(int evolvetype, int inputtype, FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr);
extern int Utoprimloop(FTYPE unew[][N2M][NPR],FTYPE pf[][N2M][NPR]);
extern int primtoUloop(FTYPE pi[][N2M][NPR],FTYPE unew[][N2M][NPR]);

extern int Utoprim(int entropyeom, FTYPE *U, struct of_geom *geom, FTYPE *pr);

extern int Utoprim_ldz(FTYPE *U, struct of_geom *geom, FTYPE *pr);
extern int wvsq_solv_ldz(FTYPE *vsq, FTYPE *W);
extern int nrunsafe(void (*funcd) (FTYPE*, FTYPE*,FTYPE*), FTYPE *guess);
extern int nrunsafeorig(void (*funcd) (FTYPE*, FTYPE*,FTYPE*), FTYPE *guess);
extern void func(FTYPE *x, FTYPE *f, FTYPE *df);
extern FTYPE rtsafe(void (*funcd) (), FTYPE x1, FTYPE x2,
		     FTYPE xacc);

extern int Utoprim_1d(FTYPE *U, struct of_geom *geom, FTYPE *pr);
extern int Utoprim_1d_opt(FTYPE *U, struct of_geom *geom, FTYPE *pr);
extern int Utoprim_2d(FTYPE *U, struct of_geom *geom, FTYPE *pr);
extern int Utoprim_1d_final(FTYPE *U, struct of_geom *geom, FTYPE *pr);
extern int Utoprim_2d_final(FTYPE *U, struct of_geom *geom, FTYPE *pr);
extern int Utoprim_5d2_final(FTYPE *U, struct of_geom *geom, FTYPE *pr);




extern void tetr_func(FTYPE tetr_cov[][NDIM], FTYPE tetr_con[][NDIM]);
extern void get_geometry(int i, int j, int loc, struct of_geom *geom);
extern int get_state(FTYPE *pr, struct of_geom *geom,
		     struct of_state *q);
extern int primtoflux(int returntype, FTYPE *pa, struct of_state *q, int dir,
	       struct of_geom *geom, FTYPE *fl);
extern void mks_source_conn(FTYPE *ph, struct of_geom *ptrgeom,
		     struct of_state *q,FTYPE *dU);
extern int source(FTYPE *pa, struct of_geom *geom,
		  FTYPE (*Uacomp)[NPR], FTYPE *Ua);

extern FTYPE taper_func(FTYPE R,FTYPE rin) ;
extern FTYPE rhor_calc(int which);
extern FTYPE rmso_calc(int which) ;
extern FTYPE uphi_isco_calc(int which,FTYPE r);

extern int set_atmosphere(int whichcond, int whichvel, struct of_geom *geom, FTYPE *pr);
extern int set_density_floors(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *scaler);
extern int fixup(int stage,FTYPE (*var)[N2M][NPR],int finalstep);
extern int fixup1zone(FTYPE *pr,struct of_geom *ptrlgeom, int finalstep);
extern int diag_fixup(FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, int finalstep,int whocalled);

extern int superdebug(FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, int whocalled);

extern void fix_flux(FTYPE (*pb)[N2M][NPR],FTYPE F1[][N2M][NPR], FTYPE F2[][N2M][NPR]) ;
extern int limit_gamma(FTYPE gammamax, FTYPE*pr,struct of_geom *geom, int finalstep);

extern int fixup_checksolution(int stage, FTYPE (*pv)[N2M][NPR],int finalstep);
extern int fixup_utoprim(int stage, FTYPE (*pv)[N2M][NPR],int finalstep);
extern int post_fixup(int stage,FTYPE (*pv)[N2M][NPR],int finalstep);
extern int pre_fixup(int stage, FTYPE (*pv)[N2M][NPR]);
extern int get_bsqflags(int stage, FTYPE (*pv)[N2M][NPR]);


extern void tet_func(FTYPE metr[][NDIM], FTYPE tetr[][NDIM]);
extern int dsyev_(char *jobz, char *uplo, int *n, FTYPE *a, int *lda,
		  FTYPE *w, FTYPE *work, int *lwork, int *iwork,
		  int *liwork, int *info);
extern int fail(int fail_type);
extern void setfailresponse(int restartonfail);
extern int setflux(void);
extern int setjetflux(void);
extern void setrestart(int*appendold);


extern int vchar(FTYPE *pr, struct of_state *q, int dir,
		 struct of_geom *geom, FTYPE *cmax, FTYPE *cmin,int *ignorecourant);
extern FTYPE chk_disp(FTYPE v);
extern void make_co_to_comov(FTYPE *ucon, FTYPE ecov[][NDIM],
			     FTYPE econ[][NDIM]);
extern void transform(FTYPE *vec, FTYPE t[][NDIM]);
extern void coeff_set(FTYPE rho, FTYPE u);
extern void transform(FTYPE *ucon, FTYPE t[][NDIM]);

extern void mhd_calc(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd);

extern int mnewt(int ntrail, int mintrial, FTYPE *p, int n, int startp, FTYPE tolx,
		 FTYPE tolf, FTYPE tolxallowed,
		 FTYPE tolfallowed, FTYPE tolxreport,
		 FTYPE tolfreport);
extern int usrfun(FTYPE *pr, int n, FTYPE *beta, FTYPE **alpha,FTYPE*normf);

extern int area_map(int call_code, int type, int size, int i, int j, FTYPE prim[][N2M][NPR]);
extern void bcon_calc(FTYPE *pr, FTYPE *ucon, FTYPE *ucov,
		      FTYPE *bcon);

extern FTYPE lc4(int updown, FTYPE detg, int mu,int nu,int kappa,int lambda);
extern void faraday_calc(int which, FTYPE *b, FTYPE *u, struct of_geom *geom, FTYPE faraday[][NDIM]);
extern void current_precalc(int which, struct of_geom *geom, struct of_state *q, SFTYPE Dt,FTYPE faraday[][3]);
extern void init_varstavg(void);
extern void final_varstavg(FTYPE IDT);
extern int set_varstavg(FTYPE tfrac);
extern void current_calc(FTYPE cfaraday[][N2M][NUMCURRENTSLOTS][3]);
extern int current_doprecalc(int which, FTYPE p[][N2M][NPR]);
extern int average_calc(int doavg);

// NR STUFF

extern int ludcmp(FTYPE **a, int n, int *indx, FTYPE *d);
extern void lubksb(FTYPE **a, int n, int *indx, FTYPE *d);
//extern FTYPE zbrent(FTYPE (*func) (FTYPE), FTYPE v1, FTYPE v2,
//		     FTYPE tol);


/* NR routines from nrutil.h */
extern int *ivector(long nl, long nh);
extern void free_ivector(int *v, long nl, long nh);
extern FTYPE *dvector(long nl, long nh);
extern void free_dvector(FTYPE *v, long nl, long nh);
extern FTYPE **dmatrix(long nrl, long nrh, long ncl, long nch);
extern void free_dmatrix(FTYPE **m, long nrl, long nrh, long ncl,
			 long nch);
extern FTYPE ***dtensor(long nrl, long nrh, long ncl, long nch,
			 long ndl, long ndh);
extern void free_dtensor(FTYPE ***t, long nrl, long nrh, long ncl,
			 long nch, long ndl, long ndh);
extern void nrerror(char error_text[]);

/* specialty functions */
extern void bondi_solve(FTYPE K, FTYPE gam, FTYPE *Rs, FTYPE *Urs,
			FTYPE *Edot);
extern FTYPE bondi_trace(FTYPE K, FTYPE gam, FTYPE edotf, FTYPE r,
			  FTYPE rs, FTYPE urs);
extern void timestep(FTYPE ndtr, FTYPE ndth);
extern FTYPE dtset(FTYPE ndtr, FTYPE ndth);

extern FTYPE bondi_trace(FTYPE K, FTYPE gam, FTYPE edotf,
			  FTYPE r, FTYPE rs, FTYPE urs);
extern void bondi_solve(FTYPE K, FTYPE gam, FTYPE *Rs,
			FTYPE *Urs, FTYPE *Edot);
extern FTYPE edot_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE dedr_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE dedur_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE d2edr2_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE d2edur2_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);
extern FTYPE d2edrdur_calc(FTYPE r, FTYPE ur, FTYPE g, FTYPE K);

extern void lower(FTYPE *a, struct of_geom *geom, FTYPE *b);
extern void lowerf(FTYPE *a, struct of_geom *geom, FTYPE *b);
extern void raise(FTYPE *v1, struct of_geom *geom, FTYPE *v2);
extern void gaussj(FTYPE **tmp, int n, FTYPE **b, int m);
extern void set_points(void);
// extern FTYPE delta(int j, int k) ;
// extern FTYPE mink(int j, int k) ;
extern void make_tetr(FTYPE *ucon, FTYPE econ[][NDIM]);
extern void dxdxprim(FTYPE *X, FTYPE r, FTYPE th, FTYPE (*dxdxp)[NDIM]);


extern FTYPE sign(FTYPE a);
extern FTYPE sign2(FTYPE a);

extern FTYPE max(FTYPE a, FTYPE b);

extern FTYPE min(FTYPE a, FTYPE b);


#if(SUPERLONGDOUBLE)
#include "mconf.h"
extern long double fabsl ( long double );
extern long double sqrtl ( long double );
extern long double cbrtl ( long double );
extern long double expl ( long double );
extern long double logl ( long double );
extern long double tanl ( long double );
extern long double atanl ( long double );
extern long double sinl ( long double );
extern long double asinl ( long double );
extern long double cosl ( long double );
extern long double acosl ( long double );
extern long double powl ( long double, long double );
extern long double tanhl ( long double );
extern long double atanhl ( long double );
extern long double sinhl ( long double );
extern long double asinhl ( long double );
extern long double coshl ( long double );
extern long double acoshl ( long double );
extern long double exp2l ( long double );
extern long double log2l ( long double );
extern long double exp10l ( long double );
extern long double log10l ( long double );
extern long double gammal ( long double );
extern long double lgaml ( long double );
extern long double jnl ( int, long double );
extern long double ynl ( int, long double );
extern long double ndtrl ( long double );
extern long double ndtril ( long double );
extern long double stdtrl ( int, long double );
extern long double stdtril ( int, long double );
extern long double ellpel ( long double );
extern long double ellpkl ( long double );
long double lgammal(long double);
extern int isfinitel ( long double );
#define finite(arg) isfinitel(arg)
extern int merror;
#else
#include <math.h>
#endif




/* NOTES:

options include:
1) WHICHVEL (which velocity to act as primitive variable) (global.h)
2) FLUXCT/FLUXCD
3) machine precision (alot of internals depend on doubles)
4) bc (just copy, gdet extrap, copy with rescale())
5) ic
6) ldz or old utoprim (UTOPRIMVERSION)
7) conservation equations (any linear sum is fine, if one can detangle solution)
8) slope_lim MC,VANL,MINM, DONOR
9) number of mnewt iterations, dampfactor, dampfactorchange
10) choice of primitive variable?
11) choice of conserved quantity?
12) order of algorithm in time(TIMEORDER)/space(lim)
13) how one splits the equations, how the variables are located
14) how one computes the flux (HLLF LAXF)
15) scale of interpolation used per quantity

want:
1) how to fixup U after calculated OR how to compute U so always good
2) vchar to be computed anisotropically
3) dt to be computed precisely "correct".

fixups include:
1) BSQORHOLIMIT/BSQOULIMIT
(global.h)
(b^2/rho b^2/u limiter determiner)
2) GAMMAMAX
(global.h, limit_gamma or check_pr)
(gamma maximum)
3) UTOPRIMADJUST
(global.h,fixup.c)
(fixes up failures from utoprim())
4) HSTEP : Numerical derivative to get conn
5) SMALL : measure of smallness
6) SLEPSILON: how far to be close to ergosphere where gcon is 0
7) inflow_check:
(bounds.c)
(forces ucon[1] to be 0 for boundary zones if inflow)
7-1) whether to check inflow for flux calculation (step_ch.c) (ZEROOUTFLOWFLUX)
8) rho>rhofl u>ufl floors
9) fix_flux:
(fixup.c)
(forces pole to have F2=0)
9-1) fix up flux (ZEROPOLEFLUX)
10) all sorts of machine precision checks e.g. V^iV^jg_{ij} >0
11) metric avoid exact pole with fabs(th)<SMALL
12) 
*/
