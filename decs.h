#include "global.h"
#include "mpidecs.h"

extern FTYPE a_p[N1M][N2M][NPR];	/* space for primitive vars */
extern FTYPE a_panalytic[N1M][N2M][NPR];	/* space for primitive vars */
extern FTYPE a_omegafanalytic[N1M][N2M][NPR];

extern FTYPE a_emf[N1M][N2M];	/* space for emf */
extern FTYPE a_vconemf[N1M][N2M][NDIM-1];	/* used for Athena EMFs */

extern FTYPE a_pleft[N1M][N2M][NPR]; /* for parabolic interpolation */
extern FTYPE a_pright[N1M][N2M][NPR]; /* for parabolic interpolation */

// for higher order RK time stepping integrations
extern FTYPE a_ulast[N1M][N2M][NPR]; 
extern FTYPE a_unew[N1M][N2M][NPR];

extern FTYPE a_dq1[N1M][N2M][NPR];	/* slopes */
extern FTYPE a_dq2[N1M][N2M][NPR];	/* slopes */
extern FTYPE a_F1[N1M][N2M][NPR];	/* fluxes */
extern FTYPE a_F2[N1M][N2M][NPR];	/* fluxes */
extern FTYPE a_F1CT[N1M][N2M][NPR];	/* fluxes */
extern FTYPE a_F2CT[N1M][N2M][NPR];	/* fluxes */
extern FTYPE a_pk[MAXDTSTAGES][N1M][N2M][NPR];	/* next-step primitives */
extern FTYPE a_prc[N1M][N2M][NPR];	/* rescaled primitives, also used for temporary storage in fixup_utoprim() */

#if(VOLUMEDIFF)
extern FTYPE a_idxvol[N1M][N2M][NDIM]; // volume regularization 1/dx for each direction
#endif

// current density stuff
extern FTYPE a_cfaraday[N1M][N2M][NUMCURRENTSLOTS][3]; // 2 for space and 2 for time
extern FTYPE a_fcon[N1M][N2M][NUMFARADAY];
extern FTYPE a_jcon[N1M][N2M][NDIM];

#if(DOAVG)
// time average stuff
extern FTYPE a_normalvarstavg[N1M][N2M][NUMNORMDUMP];
extern FTYPE a_anormalvarstavg[N1M][N2M][NUMNORMDUMP];

extern FTYPE a_jcontavg[N1M][N2M][NDIM];
extern FTYPE a_jcovtavg[N1M][N2M][NDIM];
extern FTYPE a_ajcontavg[N1M][N2M][NDIM];
extern FTYPE a_ajcovtavg[N1M][N2M][NDIM];

extern FTYPE a_massfluxtavg[N1M][N2M][NDIM];
extern FTYPE a_amassfluxtavg[N1M][N2M][NDIM];

extern FTYPE a_othertavg[N1M][N2M][NUMOTHER];
extern FTYPE a_aothertavg[N1M][N2M][NUMOTHER];

extern FTYPE a_fcontavg[N1M][N2M][NUMFARADAY];
extern FTYPE a_fcovtavg[N1M][N2M][NUMFARADAY];
extern FTYPE a_afcontavg[N1M][N2M][NUMFARADAY];
extern FTYPE a_afcovtavg[N1M][N2M][NUMFARADAY];

extern FTYPE a_tudtavg[N1M][N2M][NUMSTRESSTERMS];
extern FTYPE a_atudtavg[N1M][N2M][NUMSTRESSTERMS];
#endif

extern int a_pflag[N1M][N2M][NUMPFLAGS];	/* space for flag */

/* for debug */
#if(DODEBUG)
CTYPE a_failfloorcount[N1M][N2M][NUMTSCALES][NUMFAILFLOORFLAGS]; // number of failures and floor adjustments for each zone
#endif

/* grid functions */
extern FTYPE gengcov[NDIM][NDIM]; // gengcov/gengcon used as pointers for get_geometry() in gset() in metric.c .  This allows native coords to be different than requested coords.
extern FTYPE gengcon[NDIM][NDIM];
extern FTYPE a_conn[N1M][N2M][NDIM][NDIM][NDIM];
extern FTYPE a_conn2[N1M][N2M][NDIM];
extern FTYPE a_gcon[N1M][N2M][NPG][NDIM][NDIM];
extern FTYPE a_gcov[N1M][N2M][NPG][NDIM][NDIM];
extern FTYPE a_gdet[N1M][N2M][NPG];
extern FTYPE a_eomfunc[N1M][N2M][NPG];
//extern FTYPE a_dxdxp[N1M][N2M][NPG][NDIM][NDIM];

#if(DOENO)
extern FTYPE a_uenotmp0[N1M][N2M][NPR]; /* for ENO reconstruction of U or dU*/
extern FTYPE a_uenotmp1[N1M][N2M][NPR]; /* for ENO reconstruction of U or dU*/
extern FTYPE a_uenotmp2[N1M][N2M][NPR]; /* for ENO reconstruction of U or dU*/
#endif

extern int (*pflag)[N2M][NUMPFLAGS];
CTYPE (*failfloorcount)[N2M][NUMTSCALES][NUMFAILFLOORFLAGS];
extern FTYPE (*p)[N2M][NPR];
extern FTYPE (*panalytic)[N2M][NPR];
extern FTYPE (*omegafanalytic)[N2M][NPR];


extern FTYPE (*emf)[N2M];
extern FTYPE (*vconemf)[N2M][NDIM-1];

extern FTYPE (*pleft)[N2M][NPR];
extern FTYPE (*pright)[N2M][NPR];

extern FTYPE (*unew)[N2M][NPR];
extern FTYPE (*ulast)[N2M][NPR];

extern FTYPE (*dq1)[N2M][NPR];
extern FTYPE (*dq2)[N2M][NPR];

extern FTYPE (*F1)[N2M][NPR];
extern FTYPE (*F2)[N2M][NPR];
extern FTYPE (*F1CT)[N2M][NPR];
extern FTYPE (*F2CT)[N2M][NPR];
extern FTYPE (*pk)[N1M][N2M][NPR];
extern FTYPE (*prc)[N2M][NPR];
extern FTYPE (*conn)[N2M][NDIM][NDIM][NDIM];
extern FTYPE (*conn2)[N2M][NDIM];
extern FTYPE (*gcon)[N2M][NPG][NDIM][NDIM];
extern FTYPE (*gcov)[N2M][NPG][NDIM][NDIM];
extern FTYPE (*gdet)[N2M][NPG];
extern FTYPE (*eomfunc)[N2M][NPG];
//extern FTYPE (*dxdxp)[N2M][NPG][NDIM];

extern FTYPE (*idxvol)[N2M][NDIM];

// current density stuff
extern FTYPE (*cfaraday)[N2M][NUMCURRENTSLOTS][3];
extern FTYPE (*fcon)[N2M][NUMFARADAY];
extern FTYPE (*jcon)[N2M][NDIM];

// time average stuff
extern FTYPE (*normalvarstavg)[N2M][NUMNORMDUMP];
extern FTYPE (*anormalvarstavg)[N2M][NUMNORMDUMP];

extern FTYPE (*jcontavg)[N2M][NDIM];
extern FTYPE (*jcovtavg)[N2M][NDIM];
extern FTYPE (*ajcontavg)[N2M][NDIM];
extern FTYPE (*ajcovtavg)[N2M][NDIM];

extern FTYPE (*massfluxtavg)[N2M][NDIM];
extern FTYPE (*amassfluxtavg)[N2M][NDIM];

extern FTYPE (*othertavg)[N2M][NUMOTHER];
extern FTYPE (*aothertavg)[N2M][NUMOTHER];

extern FTYPE (*fcontavg)[N2M][NUMFARADAY];
extern FTYPE (*fcovtavg)[N2M][NUMFARADAY];
extern FTYPE (*afcontavg)[N2M][NUMFARADAY];
extern FTYPE (*afcovtavg)[N2M][NUMFARADAY];

extern FTYPE (*tudtavg)[N2M][NUMSTRESSTERMS];
extern FTYPE (*atudtavg)[N2M][NUMSTRESSTERMS];

// for image, see image.c
extern FTYPE (*pimage)[N2M][NPR];

//#if(DOENO)
extern FTYPE (*uenotmp0)[N2M][NPR]; /* for ENO reconstruction of U or dU*/
extern FTYPE (*uenotmp1)[N2M][NPR]; /* for ENO reconstruction of U or dU*/
extern FTYPE (*uenotmp2)[N2M][NPR]; /* for ENO reconstruction of U or dU*/
//#endif



extern FTYPE imageparms[NUMIMAGEPARMS];
extern int imagescale,imagewhichk,imagevartype,imagelimits;

extern SFTYPE *lumvsr,*lumvsr_tot;


/** GLOBAL PARAMETERS SECTION **/

/* physics parameters */
extern FTYPE a;
extern FTYPE gam;
extern FTYPE Bpole,Omegastar; 

/* numerical parameters */
extern int defcoord;
extern FTYPE Rin, R0, Rout, hslope;
extern FTYPE Zin, Zout;
extern FTYPE Risco,Rhor;
extern FTYPE cour;
extern FTYPE dV, dVF, dx[NDIM], startx[NDIM];
extern SFTYPE dt,t,tf;
extern FTYPE rcurr, hcurr;
extern int istart, istop, jstart, jstop;
extern int isc,iec,jsc,jec;
extern int isf1,ief1,jsf1,jef1;
extern int isf2,ief2,jsf2,jef2;
extern int ise,iee,jse,jee;
extern int isf1ct,ief1ct,jsf1ct,jef1ct;
extern int isf2ct,ief2ct,jsf2ct,jef2ct;
extern int isdq,iedq,jsdq,jedq;
extern int ispdq,iepdq,jspdq,jepdq;
extern FTYPE mydminarg1, mydminarg2;
extern long nstep;

extern int steppart,numstepparts;


/* output parameters */
extern SFTYPE DTd;
extern SFTYPE DTavg;
extern SFTYPE DTener;
extern SFTYPE DTi;
extern SFTYPE DTdebug;
extern long DTr;
extern long dump_cnt;
extern long avg_cnt;
extern long debug_cnt;
extern long image_cnt;
extern long rdump_cnt;
extern long fieldline_cnt; // assumed to keep track with images (as in diag.c), so no need to include in restart()
extern int nstroke;

/* global flags */
extern int failed;
extern int lim,fluxmethod,FLUXB,UTOPRIMVERSION,TIMEORDER;
extern FTYPE defcon;

/* diagnostics */
// don't track this separately in other regions except global region
extern SFTYPE frdot[N1][NPR];
extern SFTYPE pdottermsjet2[COMPDIM*2][NUMFLUXTERMS][NPR];
CTYPE failfloorcountlocal[NUMTSCALES][NUMFAILFLOORFLAGS]; // don't track this separately in jet
CTYPE failfloorcountlocal_tot[NUMTSCALES][NUMFAILFLOORFLAGS]; // don't track this separately in jet

// general stuff for ener.out file for regions to completely track, including terms within flux
extern int dofluxreg[NUMENERREGIONS][COMPDIM*2];
extern int enerposreg[NUMENERREGIONS][COMPDIM*2];
extern SFTYPE fladdreg[NUMENERREGIONS][NPR];
extern SFTYPE fladdreg_tot[NUMENERREGIONS][NPR];
extern SFTYPE fladdtermsreg[NUMENERREGIONS][NUMFAILFLOORFLAGS][NPR];
extern SFTYPE fladdtermsreg_tot[NUMENERREGIONS][NUMFAILFLOORFLAGS][NPR];
extern SFTYPE Ureg_init[NUMENERREGIONS][NPR];
extern SFTYPE Ureg_init_tot[NUMENERREGIONS][NPR];
extern SFTYPE pcumreg[NUMENERREGIONS][COMPDIM*2][NPR];
extern SFTYPE pcumreg_tot[NUMENERREGIONS][COMPDIM*2][NPR];
extern SFTYPE pdotreg[NUMENERREGIONS][COMPDIM*2][NPR];
extern SFTYPE pdottermsreg[NUMENERREGIONS][COMPDIM*2][NUMFLUXTERMS][NPR];
extern SFTYPE sourceaddreg[NUMENERREGIONS][NPR];
extern SFTYPE sourceaddreg_tot[NUMENERREGIONS][NPR];
extern SFTYPE sourceaddtermsreg[NUMENERREGIONS][NUMSOURCES][NPR];
extern SFTYPE sourceaddtermsreg_tot[NUMENERREGIONS][NUMSOURCES][NPR];
// used for each region, related to global quantities
// _tot quantities here are global since used in restart.
extern int *doflux;
extern int *enerpos;
extern SFTYPE *fladd;
extern SFTYPE *fladd_tot;
extern SFTYPE (*fladdterms)[NPR];
extern SFTYPE (*fladdterms_tot)[NPR];
extern SFTYPE *U_init;
extern SFTYPE *U_init_tot;
extern SFTYPE (*pcum)[NPR];
extern SFTYPE (*pcum_tot)[NPR];
extern SFTYPE (*pdot)[NPR];
extern SFTYPE (*pdotterms)[NUMFLUXTERMS][NPR];
extern SFTYPE *sourceadd;
extern SFTYPE *sourceadd_tot;
extern SFTYPE (*sourceaddterms)[NPR];
extern SFTYPE (*sourceaddterms_tot)[NPR];

// end changes after ...

/* current local position */
extern int icurr, jcurr, pcurr, ihere, jhere, phere;

/* Jon's addition */
extern int horizoni;
extern long realnstep;
extern int partialstep;
extern int mpicombine;
extern int mpicombinetype;
extern int truempicombinetype;
extern int halftimep;
extern int whichrestart;
extern int appendold;
extern int whocalleducon;
// global flags
extern long restartsteps[2];
extern int binaryoutput,sortedoutput;
extern long steptofaildump,steptofailmap;
extern int ifail,jfail,dofailmap,dofaildump,restartonfail;
// IC
extern FTYPE h_over_r;
// BC
extern FTYPE h_over_r_jet;
extern int BCtype[COMPDIM*2];
extern int rescaletype;
extern int cooling;
extern int DOENERDIAG,DOGDUMPDIAG,DORDUMPDIAG,DODUMPDIAG,DOAVGDIAG, DOIMAGEDIAG,DOAREAMAPDIAG;
extern int GAMMIEDUMP,GAMMIEIMAGE,GAMMIEENER,DODIAGS,RESTARTMODE,WHICHFILE,POSDEFMETRIC;
extern FTYPE RHOMIN,UUMIN,RHOMINLIMIT,UUMINLIMIT;
extern FTYPE prMAX[NPR];
extern FTYPE BSQORHOLIMIT,BSQOULIMIT,GAMMAMAX,GAMMADAMP,GAMMAFAIL;
extern FTYPE SAFE;
extern int debugfail;
extern FTYPE uttdiscr; // for check_pr for now
extern int jonchecks;
extern int dnumcolumns[NUMDUMPTYPES];
struct blink * blinkptr0[NUMDUMPTYPES];
struct blink * cpulinkptr0[NUMDUMPTYPES];
extern int DOCOLSPLIT[NUMDUMPTYPES];
extern int docolsplit; // global var for now
extern int nextcol;

/* physical consts */
extern FTYPE msun,lsun,G,H,C,mn,me,kb,arad,sigmasb,sigmamat,mevocsq,ergPmev;
extern FTYPE mb,M,Mdot,Mdotc;
extern FTYPE Lunit,Tunit,rho0,Munit,mdotunit,energyunit,edotunit,Pressureunit,Tempunit,Bunit,massunitPmsun;
extern FTYPE ledd,leddcode;

extern int NUMBUFFERS;
