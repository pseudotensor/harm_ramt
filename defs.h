#include "global.h"
#include "mpidefs.h"

FTYPE a_p[N1M][N2M][NPR];	/* space for primitive vars */
FTYPE a_panalytic[N1M][N2M][NPR];	/* space for primitive vars */
FTYPE a_omegafanalytic[N1M][N2M][NPR];

FTYPE a_emf[N1M][N2M];	/* space for emf */
FTYPE a_vconemf[N1M][N2M][NDIM-1];	/* used for Athena EMFs */

FTYPE a_pleft[N1M][N2M][NPR]; /* for parabolic interpolation */
FTYPE a_pright[N1M][N2M][NPR]; /* for parabolic interpolation */

// for higher order RK time stepping integrations
FTYPE a_ulast[N1M][N2M][NPR]; 
FTYPE a_unew[N1M][N2M][NPR];

FTYPE a_dq1[N1M][N2M][NPR];	/* slopes */
FTYPE a_dq2[N1M][N2M][NPR];	/* slopes */
FTYPE a_F1[N1M][N2M][NPR];	/* fluxes */
FTYPE a_F2[N1M][N2M][NPR];	/* fluxes */
FTYPE a_F1CT[N1M][N2M][NPR];	/* fluxes */
FTYPE a_F2CT[N1M][N2M][NPR];	/* fluxes */
FTYPE a_pk[MAXDTSTAGES][N1M][N2M][NPR];	/* next-step primitives */
FTYPE a_prc[N1M][N2M][NPR];	/* rescaled primitives, also used for temporary storage in fixup_utoprim() */

#if(VOLUMEDIFF)
FTYPE a_idxvol[N1M][N2M][NDIM]; // volume regularization 1/dx for each direction
#endif

// current density stuff
FTYPE a_cfaraday[N1M][N2M][NUMCURRENTSLOTS][3]; // 2 for space and 2 for time
FTYPE a_fcon[N1M][N2M][NUMFARADAY];
FTYPE a_jcon[N1M][N2M][NDIM];

#if(DOAVG)
// time average stuff
FTYPE a_normalvarstavg[N1M][N2M][NUMNORMDUMP];
FTYPE a_anormalvarstavg[N1M][N2M][NUMNORMDUMP];

FTYPE a_jcontavg[N1M][N2M][NDIM];
FTYPE a_jcovtavg[N1M][N2M][NDIM];
FTYPE a_ajcontavg[N1M][N2M][NDIM];
FTYPE a_ajcovtavg[N1M][N2M][NDIM];

FTYPE a_massfluxtavg[N1M][N2M][NDIM];
FTYPE a_amassfluxtavg[N1M][N2M][NDIM];

FTYPE a_othertavg[N1M][N2M][NUMOTHER];
FTYPE a_aothertavg[N1M][N2M][NUMOTHER];

FTYPE a_fcontavg[N1M][N2M][NUMFARADAY];
FTYPE a_fcovtavg[N1M][N2M][NUMFARADAY];
FTYPE a_afcontavg[N1M][N2M][NUMFARADAY];
FTYPE a_afcovtavg[N1M][N2M][NUMFARADAY];

FTYPE a_tudtavg[N1M][N2M][NUMSTRESSTERMS];
FTYPE a_atudtavg[N1M][N2M][NUMSTRESSTERMS];
#endif

int a_pflag[N1M][N2M][NUMPFLAGS];	/* space for flag */

/* for debug */
#if(DODEBUG)
CTYPE a_failfloorcount[N1M][N2M][NUMTSCALES][NUMFAILFLOORFLAGS]; // number of failures and floor adjustments for each zone
#endif

/* grid functions */
FTYPE gengcov[NDIM][NDIM]; // gengcov/gengcon used as pointers for get_geometry() in gset() in metric.c .  This allows native coords to be different than requested coords.
FTYPE gengcon[NDIM][NDIM];
FTYPE a_conn[N1M][N2M][NDIM][NDIM][NDIM];
FTYPE a_conn2[N1M][N2M][NDIM];
FTYPE a_gcon[N1M][N2M][NPG][NDIM][NDIM];
FTYPE a_gcov[N1M][N2M][NPG][NDIM][NDIM];
FTYPE a_gdet[N1M][N2M][NPG];
FTYPE a_eomfunc[N1M][N2M][NPG];
//FTYPE a_dxdxp[N1M][N2M][NPG][NDIM][NDIM];

#if(DOENO)
FTYPE a_uenotmp0[N1M][N2M][NPR]; /* for ENO reconstruction of U or dU*/
FTYPE a_uenotmp1[N1M][N2M][NPR]; /* for ENO reconstruction of U or dU*/
FTYPE a_uenotmp2[N1M][N2M][NPR]; /* for ENO reconstruction of U or dU*/
#endif

int (*pflag)[N2M][NUMPFLAGS];
CTYPE (*failfloorcount)[N2M][NUMTSCALES][NUMFAILFLOORFLAGS];
FTYPE (*p)[N2M][NPR];
FTYPE (*panalytic)[N2M][NPR];
FTYPE (*omegafanalytic)[N2M][NPR];


FTYPE (*emf)[N2M];
FTYPE (*vconemf)[N2M][NDIM-1];

FTYPE (*pleft)[N2M][NPR];
FTYPE (*pright)[N2M][NPR];

FTYPE (*unew)[N2M][NPR];
FTYPE (*ulast)[N2M][NPR];

FTYPE (*dq1)[N2M][NPR];
FTYPE (*dq2)[N2M][NPR];

FTYPE (*F1)[N2M][NPR];
FTYPE (*F2)[N2M][NPR];
FTYPE (*F1CT)[N2M][NPR];
FTYPE (*F2CT)[N2M][NPR];
FTYPE (*pk)[N1M][N2M][NPR];
FTYPE (*prc)[N2M][NPR];
FTYPE (*conn)[N2M][NDIM][NDIM][NDIM];
FTYPE (*conn2)[N2M][NDIM];
FTYPE (*gcon)[N2M][NPG][NDIM][NDIM];
FTYPE (*gcov)[N2M][NPG][NDIM][NDIM];
FTYPE (*gdet)[N2M][NPG];
FTYPE (*eomfunc)[N2M][NPG];
//FTYPE (*dxdxp)[N2M][NPG][NDIM];

FTYPE (*idxvol)[N2M][NDIM];

// current density stuff
FTYPE (*cfaraday)[N2M][NUMCURRENTSLOTS][3];
FTYPE (*fcon)[N2M][NUMFARADAY];
FTYPE (*jcon)[N2M][NDIM];

// time average stuff
FTYPE (*normalvarstavg)[N2M][NUMNORMDUMP];
FTYPE (*anormalvarstavg)[N2M][NUMNORMDUMP];

FTYPE (*jcontavg)[N2M][NDIM];
FTYPE (*jcovtavg)[N2M][NDIM];
FTYPE (*ajcontavg)[N2M][NDIM];
FTYPE (*ajcovtavg)[N2M][NDIM];

FTYPE (*massfluxtavg)[N2M][NDIM];
FTYPE (*amassfluxtavg)[N2M][NDIM];

FTYPE (*othertavg)[N2M][NUMOTHER];
FTYPE (*aothertavg)[N2M][NUMOTHER];

FTYPE (*fcontavg)[N2M][NUMFARADAY];
FTYPE (*fcovtavg)[N2M][NUMFARADAY];
FTYPE (*afcontavg)[N2M][NUMFARADAY];
FTYPE (*afcovtavg)[N2M][NUMFARADAY];

FTYPE (*tudtavg)[N2M][NUMSTRESSTERMS];
FTYPE (*atudtavg)[N2M][NUMSTRESSTERMS];

// for image, see image.c
FTYPE (*pimage)[N2M][NPR];

//#if(DOENO)
FTYPE (*uenotmp0)[N2M][NPR]; /* for ENO reconstruction of U or dU*/
FTYPE (*uenotmp1)[N2M][NPR]; /* for ENO reconstruction of U or dU*/
FTYPE (*uenotmp2)[N2M][NPR]; /* for ENO reconstruction of U or dU*/
//#endif



FTYPE imageparms[NUMIMAGEPARMS];
int imagescale,imagewhichk,imagevartype,imagelimits;

SFTYPE *lumvsr,*lumvsr_tot;


/** GLOBAL PARAMETERS SECTION **/

/* physics parameters */
FTYPE a;
FTYPE gam;
FTYPE Bpole,Omegastar; 

/* numerical parameters */
int defcoord;
FTYPE Rin, R0, Rout, hslope;
FTYPE Zin, Zout;
FTYPE Risco,Rhor;
FTYPE cour;
FTYPE dV, dVF, dx[NDIM], startx[NDIM];
SFTYPE dt,t,tf;
FTYPE rcurr, hcurr;
int istart, istop, jstart, jstop;
int isc,iec,jsc,jec;
int isf1,ief1,jsf1,jef1;
int isf2,ief2,jsf2,jef2;
int ise,iee,jse,jee;
int isf1ct,ief1ct,jsf1ct,jef1ct;
int isf2ct,ief2ct,jsf2ct,jef2ct;
int isdq,iedq,jsdq,jedq;
int ispdq,iepdq,jspdq,jepdq;
FTYPE mydminarg1, mydminarg2;
long nstep;

int steppart,numstepparts;


/* output parameters */
SFTYPE DTd;
SFTYPE DTavg;
SFTYPE DTener;
SFTYPE DTi;
SFTYPE DTdebug;
long DTr;
long dump_cnt;
long avg_cnt;
long debug_cnt;
long image_cnt;
long rdump_cnt;
long fieldline_cnt; // assumed to keep track with images (as in diag.c), so no need to include in restart()
int nstroke;

/* global flags */
int failed;
int lim,fluxmethod,FLUXB,UTOPRIMVERSION,TIMEORDER;
FTYPE defcon;

/* diagnostics */
// don't track this separately in other regions except global region
SFTYPE frdot[N1][NPR];
SFTYPE pdottermsjet2[COMPDIM*2][NUMFLUXTERMS][NPR];
CTYPE failfloorcountlocal[NUMTSCALES][NUMFAILFLOORFLAGS]; // don't track this separately in jet
CTYPE failfloorcountlocal_tot[NUMTSCALES][NUMFAILFLOORFLAGS]; // don't track this separately in jet

// general stuff for ener.out file for regions to completely track, including terms within flux
int dofluxreg[NUMENERREGIONS][COMPDIM*2];
int enerposreg[NUMENERREGIONS][COMPDIM*2];
SFTYPE fladdreg[NUMENERREGIONS][NPR];
SFTYPE fladdreg_tot[NUMENERREGIONS][NPR];
SFTYPE fladdtermsreg[NUMENERREGIONS][NUMFAILFLOORFLAGS][NPR];
SFTYPE fladdtermsreg_tot[NUMENERREGIONS][NUMFAILFLOORFLAGS][NPR];
SFTYPE Ureg_init[NUMENERREGIONS][NPR];
SFTYPE Ureg_init_tot[NUMENERREGIONS][NPR];
SFTYPE pcumreg[NUMENERREGIONS][COMPDIM*2][NPR];
SFTYPE pcumreg_tot[NUMENERREGIONS][COMPDIM*2][NPR];
SFTYPE pdotreg[NUMENERREGIONS][COMPDIM*2][NPR];
SFTYPE pdottermsreg[NUMENERREGIONS][COMPDIM*2][NUMFLUXTERMS][NPR];
SFTYPE sourceaddreg[NUMENERREGIONS][NPR];
SFTYPE sourceaddreg_tot[NUMENERREGIONS][NPR];
SFTYPE sourceaddtermsreg[NUMENERREGIONS][NUMSOURCES][NPR];
SFTYPE sourceaddtermsreg_tot[NUMENERREGIONS][NUMSOURCES][NPR];
// used for each region, related to global quantities
// _tot quantities here are global since used in restart.
int *doflux;
int *enerpos;
SFTYPE *fladd;
SFTYPE *fladd_tot;
SFTYPE (*fladdterms)[NPR];
SFTYPE (*fladdterms_tot)[NPR];
SFTYPE *U_init;
SFTYPE *U_init_tot;
SFTYPE (*pcum)[NPR];
SFTYPE (*pcum_tot)[NPR];
SFTYPE (*pdot)[NPR];
SFTYPE (*pdotterms)[NUMFLUXTERMS][NPR];
SFTYPE *sourceadd;
SFTYPE *sourceadd_tot;
SFTYPE (*sourceaddterms)[NPR];
SFTYPE (*sourceaddterms_tot)[NPR];

// end changes after ...

/* current local position */
int icurr, jcurr, pcurr, ihere, jhere, phere;

/* Jon's addition */
int horizoni;
long realnstep;
int partialstep;
int mpicombine;
int mpicombinetype;
int truempicombinetype;
int halftimep;
int whichrestart;
int appendold;
int whocalleducon;
// global flags
long restartsteps[2];
int binaryoutput,sortedoutput;
long steptofaildump,steptofailmap;
int ifail,jfail,dofailmap,dofaildump,restartonfail;
// IC
FTYPE h_over_r;
// BC
FTYPE h_over_r_jet;
int BCtype[COMPDIM*2];
int rescaletype;
int cooling;
int DOENERDIAG,DOGDUMPDIAG,DORDUMPDIAG,DODUMPDIAG,DOAVGDIAG, DOIMAGEDIAG,DOAREAMAPDIAG;
int GAMMIEDUMP,GAMMIEIMAGE,GAMMIEENER,DODIAGS,RESTARTMODE,WHICHFILE,POSDEFMETRIC;
FTYPE RHOMIN,UUMIN,RHOMINLIMIT,UUMINLIMIT;
FTYPE prMAX[NPR];
FTYPE BSQORHOLIMIT,BSQOULIMIT,GAMMAMAX,GAMMADAMP,GAMMAFAIL;
FTYPE SAFE;
int debugfail;
FTYPE uttdiscr; // for check_pr for now
int jonchecks;
int dnumcolumns[NUMDUMPTYPES];
struct blink * blinkptr0[NUMDUMPTYPES];
struct blink * cpulinkptr0[NUMDUMPTYPES];
int DOCOLSPLIT[NUMDUMPTYPES];
int docolsplit; // global var for now
int nextcol;

/* physical consts */
FTYPE msun,lsun,G,H,C,mn,me,kb,arad,sigmasb,sigmamat,mevocsq,ergPmev;
FTYPE mb,M,Mdot,Mdotc;
FTYPE Lunit,Tunit,rho0,Munit,mdotunit,energyunit,edotunit,Pressureunit,Tempunit,Bunit,massunitPmsun;
FTYPE ledd,leddcode;

int NUMBUFFERS;
