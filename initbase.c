
#include "decs.h"

int pre_init(int argc, char *argv[])
{
  int ii;
  int dir,k,sc,fl,floor,enerregion;
  int tscale;
  int i,j;

  // things initialized whether restarting or init fresh

  ranc(0);  // power up random number generator in case used without init



  /////////////////
  //
  // setup files for writing and reading
  //

  makedirs();



  // below 2 now determined at command line.  See init_mpi() in init_mpi.c (myargs and init_MPI).
  //  RESTARTMODE=0;// whether restarting from rdump or not (0=no, 1=yes)
  //WHICHFILE=0; // see diag.c for dump_cnt and image_cnt starts
  // user defined parameter
  restartonfail=0; // whether we are restarting on failure or not and want special diagnostics

  if(WHICHVEL==VEL3){
    jonchecks=1; // whether to include jon's checks to make sure u^t real and some rho/u limits throughout code
    jonchecks=0;
  }
  else jonchecks=0; // not relevant

  // choice
  periodicx1=0;
  periodicx2=0;

  if(USEMPI&&USEROMIO){
    binaryoutput=MIXEDOUTPUT; // choice: mixed or binary
    sortedoutput=SORTED; // no choice
  }
  else{
    // choice
    binaryoutput=TEXTOUTPUT;
    sortedoutput=SORTED;
  }


  // init MPI (assumes nothing in set_arrays.c used here)
  init_mpi(argc, argv);

  // init arrays
  set_arrays();



  init_dumps();

  // must go here b4 restart if restarting
  ENERREGIONLOOP{
    // used for each region, related to global quantities
    // no need to initialize _tot quantities, they are overwritten during MPI sum in diag.c
    fladd=fladdreg[enerregion];
    fladdterms=fladdtermsreg[enerregion];
    U_init=Ureg_init[enerregion];
    pcum=pcumreg[enerregion];
    pdot=pdotreg[enerregion];
    pdotterms=pdottermsreg[enerregion];
    sourceaddterms=sourceaddtermsreg[enerregion];
    sourceadd=sourceaddreg[enerregion];

    PDUMPLOOP{
      fladd[k] = 0;
      FLOORLOOP fladdterms[floor][k]=0;
      U_init[k] = 0;
      DIRLOOP{
	pcum[dir][k]=0;
	pdot[dir][k]=0;
	FLLOOP pdotterms[dir][fl][k]=0;
	if(enerregion==0) FLLOOP pdottermsjet2[dir][fl][k]=0; // needed for other not-flux cpus!
      }
      sourceadd[k] = 0;
      SCLOOP sourceaddterms[sc][k] = 0;
    }
    if(DOLUMVSR) if(enerregion==0) for(ii=0;ii<ncpux1*N1;ii++) lumvsr[ii]=0;
  }
  
  // start counter
  // full loop since user may choose to count something in boundary zones
  if(DODEBUG) FULLLOOP TSCALELOOP FLOORLOOP failfloorcount[i][j][tscale][floor]=0;

  // zero out jcon since outer boundaries not set ever since j^\mu involves spatial derivatives that don't exist outside a certain point
  for(k=0;k<NDIM;k++) FULLLOOP jcon[i][j][k]=0.0;


  return(0);
}

int init_defgrid(void)
{
  // sets metric
  a=0.0;
  // set coordinates
  defcoord=0;
  // sets parameters of coordinates, default changes
  R0 = 0.0;
  Rin = 0.98 * Rhor;
  Rout = 40.;
  hslope = 0.3;

  return(0);
}


int init_defglobal(void)
{
  int i;

#if(!PRODUCTION)
  debugfail=1;
#else
  debugfail=0; // no messages in production -- assumes all utoprim-like failures need not be debugged
#endif
 // whether to show debug into on failures.  Desirable to turn off if don't care and just want code to burn away using given hacks/kludges
  // 0: no messages
  // 1: critical messages
  // 2: all failure messages

  // included in rdump
  defcon = 1.0;
  /* maximum increase in timestep */
  SAFE=1.3;
  whichrestart = 0;
  restartsteps[0] = 0;
  restartsteps[1] = 0;
  nstep = realnstep = 0;
  failed = 0;
  cour = 0.9;
  //  lim = MC;
#include "step_ch.h" // to get nemonics
#if(EOMTYPE==EOMGRMHD)
  lim = PARA;
  TIMEORDER=4;
  fluxmethod=HLLFLUX;
  FLUXB = FLUXCTTOTH;
  UTOPRIMVERSION=UTOPRIM2DFINAL;
#elif(EOMTYPE==EOMFFDE)
  // PARA and TO=4 and HLL not trustable in FFDE so far
  lim = MC;
  TIMEORDER=2;
  //  fluxmethod=HLLFLUX;
  fluxmethod=LAXFFLUX;
  FLUXB = FLUXCTTOTH;
  UTOPRIMVERSION=UTOPRIM2DFINAL;
#endif
  t = 0.;
  tf = 1.0;
  // initial dt
  dt = 1.0e-5;
  DTd = DTavg = DTdebug = 1.0;
  DTener=1.0;
  DTi=1.0;
  DTr=1;

  GAMMIEDUMP=0;// whether in Gammie output types (sets filename with 3 numbers and doesn't output i,j)
  GAMMIEIMAGE=0; // Gammie filename and single density output
  GAMMIEENER=0; // currently doing gener as well as ener, but this would also output some messages in gammie form

  // DOCOLSPLIT
  //
  // 0: don't ..
  // 1: split dump files into 1 column per file with ending number being column number
  // works in MPI-mode as well.  ROMIO+DOCOLSPLIT is useful for tungsten with low memory and small files to avoid diskspace and memory limits.

  // default
  for(i=0;i<NUMDUMPTYPES;i++){
    DOCOLSPLIT[i]=0;
  }
  // otherwise specify for each dump type

  DODIAGS=1; // whether to do diagnostics
  // specify individual diagnostics to be done
  DOENERDIAG=1;
  DOGDUMPDIAG=1;
  DORDUMPDIAG=1;
  DODUMPDIAG=1;
  if(DOAVG){
    DOAVGDIAG=1; // choice
  }
  else DOAVGDIAG=0; // no choice
  DOIMAGEDIAG=1;
  DOAREAMAPDIAG=1;

  POSDEFMETRIC=0; // see metric.c, bounds.c, and coord.c

  rescaletype=1;
  // 0: normal
  // 1: extended b^2/rho in funnel
  // 2: conserve E and L along field lines
  //   1 or 2 required to conserve E and L along field line

  /** FIXUP PARAMETERS **/
  RHOMIN=1.e-4;
  UUMIN =1.e-6;
  RHOMINLIMIT=1.e-20;
  UUMINLIMIT =1.e-20;

  // limit of B^2/rho if using that flag
  BSQORHOLIMIT=20.0;
  BSQOULIMIT=100.0;
  GAMMADAMP=5.0;

  // GODMARK -- unstable beyond about 25, but can sometimes get away with 1000
  GAMMAMAX=25.0; // when we think gamma is just too high and may cause unstable flow, but solution is probably accurate.
  GAMMAFAIL=100.0*GAMMAMAX; // when we think gamma is rediculous as to mean failure and solution is not accurate.
  prMAX[RHO]=20.0;
  prMAX[UU]=20.0;
  prMAX[U1]=100.0;
  prMAX[U2]=100.0;
  prMAX[U3]=100.0;
  prMAX[B1]=100.0;
  prMAX[B2]=100.0;
  prMAX[B3]=100.0;

  // some physics
  gam=4./3.; // ultrarelativistic gas, assumes pgas<prad and radiation
	     // doesn't escape
  // gam=5/3 for non-relativistic gas, such as neucleons in collapsar model
  cooling=0;
  // cooling: 0: no cooling 1: adhoc thin disk cooling 2: neutrino cooling for collapsar model

  // boundary conditions
  BCtype[X1UP]=OUTFLOW;
  BCtype[X1DN]=OUTFLOW;
  BCtype[X2UP]=POLARAXIS;
  BCtype[X2DN]=POLARAXIS;

  return(0);
}


int init_defconsts(void)
{


  // some constants that define neutrino cooling and spin evolution
  // cgs
  msun=1.989E33;
  lsun=3.89E33;
  G=6.672E-8;
  H=6.6262E-27;
  C=2.99792458E10;
  mn=1.675E-24;
  me=9.10938188E-28;
  kb=1.3807E-16;
  arad=8.*pow(M_PI,5.0)*pow(kb,4.0)/(15*C*C*C*H*H*H);
  sigmasb=arad*C/4.0;
  sigmamat=6.652E-29*100*100;
  mevocsq=1.783E-27;
  ergPmev=1.602E-6;

  // GRB
  mb=mn;
  //  M=3.*msun;
  Mdot=0.1*msun;
  Mdotc=0.35; // approx code units mass accretion rate
  
  // units
  Lunit=G*M/(C*C);
  Tunit=G*M/(C*C*C);
  rho0=Mdot/(Mdotc*Lunit*Lunit*C);
  Munit=rho0*pow(Lunit,3.0);
  mdotunit=Munit/Tunit;
  energyunit=Munit*C*C;
  edotunit=energyunit/Tunit;
  Pressureunit=rho0*C*C;
  Tempunit=Pressureunit*mb/(rho0*kb);
  Bunit=C*sqrt(rho0);
  massunitPmsun=Munit/msun;
  
  // physics stuff
  ledd=4.*M_PI*C*G*M*mb/sigmamat;
  leddcode=ledd/edotunit;

 
  return(0);
}





int init(int argc, char *argv[])
{
  int whichvel, whichcoord;
  int initreturn;
  int i = 0, j = 0, k = 0;
  struct of_geom geom;
  SFTYPE rh,r,th,X[NDIM];
  /* for magnetic field */
  /* for magnetic field */
  //  SFTYPE A[N1 + 1][N2 + 1];
  FTYPE (*A)[N2M];
  extern void filterffde(int i, int j, FTYPE *pr);



  fprintf(stderr,"Start init\n"); fflush(stderr);

  pre_init(argc,argv);
  pre_init_specific_init();

  if(RESTARTMODE==1){
    if (restart_init(WHICHFILE) >= 1) {
      dualfprintf(fail_file, "main:restart_init: failure\n");
      return(1);
    }
  }
  else{
    // rest is for fresh start, and thus the above restart must include all these parameters since not otherwise set!


    // define parameters
    init_defglobal(); // init default global parameters
    init_defconsts(); // init default physical constants
    init_defgrid(); // init default grid
    init_grid(); // request choices for grid/coordinate/metric parameters
    set_coord_parms();
    write_coord_parms(); // output coordinate parameters to file
    init_global(); // request choices for global parameters

    // some calculations, althogh perhaps calculated already, definitely need to make sure computed
    Rhor=rhor_calc(0);
    Risco=rmso_calc(PROGRADERISCO);
  
    // get grid
    set_grid();

    // diagnostic
    // determine nature of inner radial edge (assumes myid==0 is always there)
    if(myid==0){
      coord(-N1BND, 0, FACE1, X);
      bl_coord(X, &r, &th);
      trifprintf("rmin: %21.15g\n", r);
      trifprintf("rmin/rh: %21.15g\n", r / Rhor );
      //    trifprintf("rmin/rsing: %21.15g\n", r / (a+SMALL));
      if(r/Rhor<=1.0){
	trifprintf("inner grid is inside horizon\n");
      }
      else{
	trifprintf("inner grid is outside horizon\n");
      }
    }

    ///////////////////////////////////
    //
    // Assign primitive variables
    //
    ///////////////////////////////////
    trifprintf("Assign primitives\n");

    // assume we start in bl coords and convert to KSprim
    ZSLOOP(-N1BND, N1 - 1+N1BND, -N2BND, N2 - 1+N2BND) {
      initreturn=init_dsandvels(&whichvel, &whichcoord,i,j,p[i][j]); // request densities for all computational centers
      if(initreturn>0) return(1);
      else{
	// transform from whichcoord to MCOORD
	if (bl2met2metp2v(whichvel,whichcoord,p[i][j], i,j) >= 1)
	  FAILSTATEMENT("init.c:init()", "bl2ks2ksp2v()", 1);
      }
    }

    /////////////////////////////
    //
    // normalize density if wanted
    //
    /////////////////////////////// 
    // at this point densities are still standard, so just send "p"
    trifprintf("Normalize densities\n");
    normalize_densities(p);


    /////////////////////////////
    //
    // Define an atmosphere if wanted
    //
    /////////////////////////////// 

#if(EOMTYPE==EOMGRMHD)
    // normalized atmosphere
    trifprintf("Add atmosphere\n");
    ZSLOOP(0, N1 - 1, 0, N2 - 1) {
      initreturn=init_atmosphere(&whichvel, &whichcoord,i,j,p[i][j]);
      if(initreturn>0) return(1);
      else{
	// transform from whichcoord to MCOORD
	if (bl2met2metp2v(whichvel, whichcoord,p[i][j], i,j) >= 1)
	  FAILSTATEMENT("init.c:init()", "bl2ks2ksp2v()", 1);
      }
    }
#endif


    /////////////////////////////
    //
    // Fixup and Bound variables since some primitive quantities may have changed
    // These may be used to define vector potential below
    // Also setup pre_fixup() type quantities
    //
    /////////////////////////////// 
    trifprintf("Fixup and Bound #1\n");

#if(EOMTYPE!=EOMFFDE)
    // assume EOMFFDE doesn't use "density/ie" to set field, so no need to bound, and no field definition is bad for EOMFFDE
#if(FIXUPAFTERINIT)
    if(fixup(STAGEM1,p,0)>=1)
      FAILSTATEMENT("init.c:init()", "fixup()", 1);
#endif

    if (bound_prim(STAGEM1,p) >= 1)
      FAILSTATEMENT("init.c:init()", "bound_prim()", 1);

    if(pre_fixup(STAGEM1,p)>=1)
      FAILSTATEMENT("init.c:init()", "pre_fixup()", 1);
#endif


    /////////////////////////////
    //
    // Initialize field from vector potential
    //
    /////////////////////////////// 
    //  A=(SFTYPE (*)[N2M])(&pk[0][0][0][0]);
    A=emf;


    trifprintf("Initialize field from vector potential\n");
    ZSLOOP(-N1BND, N1-1+N1BND, -N2BND, N2-1+N2BND){
      //    ZSLOOP(0, N1, 0, N2){
      A[i][j] = 0.;
    }
    trifprintf("Initialize field from vector potential assignasdf\n");
    ZSLOOP(-N1BND, N1-1+N1BND, -N2BND, N2-1+N2BND){
      //    ZSLOOP(0, N1, 0, N2){
      init_vpot(i,j,&A[i][j]); // request vector potential for all computational corners
    }
    trifprintf("Initialize field from vector potential assign\n");

    init_vpot2field(A,p);

    // need normalization to account for any modifications
    normalize_field(p);


    // pre-bound dump file
    if (dump(9997) >= 1){
      dualfprintf(fail_file,"unable to print dump file\n");
      return (1);
    }

    // may want to use the field to adjust the solution
    trifprintf("init_postfield\n");
    //    MYFUN(init_postfield(p),"initbase.c:init()","init_postfield()",1);
    init_postfield(p);


    // pre-bound dump file
    if (dump(9998) >= 1){
      dualfprintf(fail_file,"unable to print dump file\n");
      return (1);
    }



    /////////////////////////////
    //
    // Filter to get correct degenerate FFDE solution (last thing that should be done to modify initial solution)
    //
    /////////////////////////////// 

#if(EOMTYPE==EOMFFDE)
    trifprintf("filtered to FFDE\n");
    // filter to get force-free
    ZSLOOP(-N1BND, N1-1+N1BND, -N2BND, N2-1+N2BND){
      filterffde(i,j,p[i][j]);
    }
#endif


    // need normalization to account for any modifications
    normalize_field(p);

    // may want to use the field to adjust the solution
    trifprintf("init_postfield\n");
    init_postfield(p);



    /////////////////////////////
    //
    // Filter to get correct degenerate FFDE solution (last thing that should be done to modify initial solution)
    //
    /////////////////////////////// 

#if(EOMTYPE==EOMFFDE)
    trifprintf("filtered to FFDE\n");
    // filter to get force-free
    ZSLOOP(-N1BND, N1-1+N1BND, -N2BND, N2-1+N2BND){
      filterffde(i,j,p[i][j]);
    }
#endif


    // copy over initial solution as analytic solution
    ZSLOOP(-N1BND, N1-1+N1BND, -N2BND, N2-1+N2BND) PLOOP{
      panalytic[i][j][k]=p[i][j][k];
    }

    if (bound_prim(STAGEM1,panalytic) >= 1)
      FAILSTATEMENT("init.c:init()", "bound_prim()", 3);


    /////////////////////////////
    //
    // Fixup and Bound variables since field may have changed
    // Also setup pre_fixup() type quantities
    //
    /////////////////////////////// 

    trifprintf("Fixup and Bound #2\n");

#if(FIXUPAFTERINIT)
    if(fixup(STAGEM1,p,0)>=1)
      FAILSTATEMENT("init.c:init()", "fixup()", 2);
#endif


    // pre-bound dump file
    if (dump(9999) >= 1){
      dualfprintf(fail_file,"unable to print dump file\n");
      return (1);
    }


    if (bound_prim(STAGEM1,p) >= 1)
      FAILSTATEMENT("init.c:init()", "bound_prim()", 2);

    if(pre_fixup(STAGEM1,p)>=1)
      FAILSTATEMENT("init.c:init()", "pre_fixup()", 2);

  }
  post_init_specific_init();
  post_init();

  trifprintf("end init.c\n");
  return (0);

}

// assumes normal field p
int vpot2field(FTYPE (*A)[N2M],FTYPE p[][N2M][NPR])
{
  int i,j;
  struct of_geom geom;

  ZSLOOP(-N1BND, N1-1+N1BND-1, -N2BND, N2-1+N2BND-1){//  ZLOOP {
    get_geometry(i, j, CENT, &geom);

    /* flux-ct */
    p[i][j][B1] = (-A[i][j]
		   + A[i][j + 1]
		   - A[i + 1][j] 
		   + A[i + 1][j + 1]) 
      / (2. * dx[2] * geom.g);
    p[i][j][B2] = -(-A[i][j] 
		    - A[i][j + 1]
		    + A[i + 1][j]
		    + A[i + 1][j + 1])
      / (2. * dx[1] * geom.g);

    // we assume only a poloidal field is required
    //    p[i][j][B3] = 0.;
  }

  ZSLOOP(-N1BND, N1-1+N1BND, -N2BND, N2-1+N2BND){//  ZLOOP {
    //    if(i==-N1BND){
    //      p[i][j][B1]=p[i+1][j][B1];
    //      p[i][j][B2]=p[i+1][j][B2];
    //    }
    if(i==N1-1+N1BND){
      p[i][j][B1]=p[i-1][j][B1];
      p[i][j][B2]=p[i-1][j][B2];
    }
    else if(j==N2-1+N1BND){
      p[i][j][B1]=p[i][j-1][B1];
      p[i][j][B2]=p[i][j-1][B2];
    }
    //else if(j==-N2BND){
    //      p[i][j][B1]=p[i][j+1][B1];
    //      p[i][j][B2]=p[i][j+1][B2];
    //    }
    //    dualfprintf(fail_file,"%d %d %21.15g %21.15g\n",i,j,p[i][j][B1],p[i][j][B2]);
  }
  p[-N1BND][-N2BND][B1]=p[-N1BND+1][-N2BND][B1];
  p[-N1BND][-N2BND][B2]=p[-N1BND+1][-N2BND][B2];

  p[-N1BND][N2-1+N2BND][B1]=p[-N1BND+1][N2-1+N2BND][B1];
  p[-N1BND][N2-1+N2BND][B2]=p[-N1BND+1][N2-1+N2BND][B2];


  

  return(0);
}


int post_init(void)
{
  int dir,k;

  trifprintf("begin: post_init\n");

  // in synch always here
  if (error_check(3)) {
    dualfprintf(fail_file, "error_check detected failure at main:1\n");
    dualfprintf(fail_file, "Bad initial conditions\n");
    myexit(1);
  }

  // all calculations that do not change for any initial conditions or problem setup or restart conditions

  // these 2 below are used prior, but not initialized otherwised on restart
  // some calculations
  Rhor=rhor_calc(0);
  Risco=rmso_calc(PROGRADERISCO);
  
  find_horizon();
  setflux();
  if(DOJETDIAG)  setjetflux();
 

  if(RESTARTMODE!=0){
    setfailresponse(restartonfail);
  }

#if(1)
  // setup faraday so t=0 dump has diagnostics
  // may want currents for dump=0 (time-derivative terms will be wrong)
  current_doprecalc(0,p);
  current_doprecalc(1,p);
  current_doprecalc(2,p);
  current_doprecalc(3,p);
  // compute current
  current_calc(cfaraday);
#else

 // setup faraday so t=0 dump has diagnostics
  current_doprecalc(3,p);
  if(WHICHCURRENTCALC==1){
    // compute faraday components needed for time centering of J
    current_doprecalc(0,p);
  }



#endif


  trifprintf("end: post_init\n");
  return(0);
}

int find_horizon(void)
{
  int i, j, k;
  FTYPE r1, r2;
  FTYPE X[NDIM];
  FTYPE r, th;
  int horizoncpupos1, gotit;
  FTYPE horizonvalue;
  // called after grid is setup for all cpus


  trifprintf("begin: find_horizon ... ");

  // find cpu column that brackets the horizon and determine the
  // i-offset of horizon

  horizonvalue = Rhor;
  horizoni = -100;
  gotit = 0;
  for (k = numprocs - 1; k >= 0; k--) { // should get done by first row
    if (k == myid) {
      for (i = N1 - 1; i >= 0; i--) {
        j = N2 / 2;             // doesn't matter
        coord(i, j, CENT, X);
        bl_coord(X, &r1, &th);
        coord(i + 1, j, CENT, X);
        bl_coord(X, &r2, &th);
        if (fabs(r1 - horizonvalue) <= (r2 - r1)) {     // find horizon
          horizoni = i;
          horizoncpupos1 = mycpupos[1];
          break;
        }
      }
    }
    if (numprocs > 0) {
#if(USEMPI)
      MPI_Bcast(&horizoni, 1, MPI_INT, k, MPI_COMM_WORLD);
      MPI_Bcast(&horizoncpupos1, 1, MPI_INT, k, MPI_COMM_WORLD);
#endif
    }
    if (horizoni >= 0)
      gotit = 1;                // can stop entire process
    if (mycpupos[1] != horizoncpupos1) {
      horizoni = -100;
    }                           // reset if not right cpu group
    if (gotit)
      break;
  }
  trifprintf("horizoni: %d horizoncpupos1: %d\n", horizoni,
             horizoncpupos1);
  // just a check
  dualfprintf(log_file,"horizoni: %d mycpupos[1]: %d horizoncpupos1: %d\n", horizoni, mycpupos[1], horizoncpupos1);

  trifprintf("end: find_horizon\n");
  return(0);
}


// determine if this cpu is doing what flux
int setflux(void)
{
  int dir;

  // setup pointers
  enerpos=enerposreg[GLOBALENERREGION];
  doflux=dofluxreg[GLOBALENERREGION];

  // all CPUs , total space for global region
  enerpos[X1DN]=0;
  enerpos[X1UP]=N1-1;    
  enerpos[X2DN]=0;
  enerpos[X2UP]=N2-1;


  // only 0 through N-1 mean do flux
  if(mycpupos[1]==0){
    doflux[X1DN]=0; // or horizoni
    trifprintf("proc: %d doing flux X1DN\n",myid);
  }
  else doflux[X1DN]=-100;

  if(mycpupos[1]==ncpux1-1){
    doflux[X1UP]=N1;
    trifprintf("proc: %d doing flux X1UP\n",myid);
  }
  else doflux[X1UP]=-100;

  if(mycpupos[2]==0){
    doflux[X2DN]=0;
    trifprintf("proc: %d doing flux X2DN\n",myid);
  }
  else doflux[X2DN]=-100;

  if(mycpupos[2]==ncpux2-1){
    doflux[X2UP]=N2;
    trifprintf("proc: %d doing flux X2UP\n",myid);
  }
  else doflux[X2UP]=-100;
  // fluxes are on edges of zone, so 0 and N are on edge fluxes

  DIRLOOP trifprintf("proc: %d %d doflux[%d]=%d enerpos[%d]=%d\n",myid,GLOBALENERREGION,dir,doflux[dir],dir,enerpos[dir]);


  return(0);
}

// GODMARK
// the below only works for a grid with full Pi.  Crashes code at runtime otherwise!  should fix.

// determine the flux positions for each CPU for the jet region (jetedge)
// AND the range of volume integration for cons. variables in jet (jetpos).
int setjetflux(void)
{
  FTYPE X[NDIM],r,th;
  int i,j,k;
  int dir;
  FTYPE startth,endth,thetajet;
  int jetedge[NUMJETS];

  if(defcoord==3){
    dualfprintf(fail_file,"setjetflux() not setup to work for non-full-Pi grids\n");
    myexit(1);
  }

  // jet region is assumed to be within a constant theta slice
  // this is theta w.r.t. polar axis
  thetajet=M_PI*0.5-h_over_r_jet;
  // find j for which theta is at our prespecified point

  i=0;j=0;
  coord(i, j, FACE2, X);
  bl_coord(X, &r, &th);
  startth=th;
  i=0;j=N2;
  coord(i, j, FACE2, X);
  bl_coord(X, &r, &th);
  endth=th;

  // assumes 0<thetajet<Pi/2
  if((fabs(startth-thetajet)/thetajet<1E-8)||
     (fabs(endth-thetajet)/thetajet<1E-8)||
     (fabs(startth-(M_PI-thetajet))/thetajet<1E-8)||
     (fabs(endth-(M_PI-thetajet))/thetajet<1E-8)
     ){
    dualfprintf(fail_file,"thetajet is on top of grid, move h_over_r_jet a bit\n");
    myexit(1);
  }


  ////////////////////
  //
  // INNERJET
  //

  // setup pointers
  enerpos=enerposreg[INNERJETREGION];
  doflux=dofluxreg[INNERJETREGION];


  // see if jet edge is related to this CPU
  // assumes increasing j is increasing th
  if((startth<=thetajet)&&(endth<=thetajet)){
    // if cpu entirely within inner theta jet
    enerpos[X1DN]=0;
    enerpos[X1UP]=N1-1;
    enerpos[X2DN]=0;
    enerpos[X2UP]=N2-1;
    jetedge[INNERJET]=-100;
  }
  else if((startth<thetajet)&&(endth>thetajet)){
    // if inner jet edge is on this CPU but not on boundary
    enerpos[X1DN]=0;
    enerpos[X1UP]=N1-1;
    enerpos[X2DN]=0;
    i=0;
    for(j=0;j<=N2;j++){
      coord(i, j, FACE2, X);
      bl_coord(X, &r, &th);
      // look for switch from below to above thetajet at inner theta jet edge
      if(th>thetajet){
	enerpos[X2UP]=j-1;
	jetedge[INNERJET]=j;
	break;
      }
    }
  }
  else if((startth>=thetajet)&&(endth>=thetajet)){
    // if cpu is entirely not contained in inner jet
    enerpos[X1DN]=-100;
    enerpos[X1UP]=-100;
    enerpos[X2DN]=-100;
    enerpos[X2UP]=-100;
    jetedge[INNERJET]=-100;
  }
  else{
    trifprintf("problem with INNERJET setjetflux()\n");
    myexit(1);
  }



  // left edge (any directional condition would do)
  if((enerpos[X1DN]!=-100)&&(mycpupos[1]==0)){
    doflux[X1DN]=0; // or horizoni
    trifprintf("proc: %d doing inner jet flux X1DN\n",myid);
  }
  else doflux[X1DN]=-100;

  // right edge (any directional condition would do)
  if((enerpos[X1DN]!=-100)&&(mycpupos[1]==ncpux1-1)){
    doflux[X1UP]=N1;
    trifprintf("proc: %d doing inner jet flux X1UP\n",myid);
  }
  else doflux[X1UP]=-100;

  // lower theta boundary
  if(mycpupos[2]==0){
    doflux[X2DN]=0;
    trifprintf("proc: %d doing inner jet flux X2DN\n",myid);
  }
  else doflux[X2DN]=-100;
  
  // upper theta boundary
  if(jetedge[INNERJET]!=-100){ // only get flux if CPU has edge
    doflux[X2UP]=jetedge[INNERJET];
    trifprintf("proc: %d doing inner jet flux X2UP\n",myid);
  }
  else doflux[X2UP]=-100;

  DIRLOOP trifprintf("proc: %d %d doflux[%d]=%d enerpos[%d]=%d\n",myid,INNERJETREGION,dir,doflux[dir],dir,enerpos[dir]);


  /////////////////////
  //
  // OUTERJET
  //

  // setup pointers
  enerpos=enerposreg[OUTERJETREGION];
  doflux=dofluxreg[OUTERJETREGION];


  // see if outer jet edge is related to this CPU
  if((startth<=M_PI-thetajet)&&(endth<=M_PI-thetajet)){
    // if cpu entirely not within outer jet region
    enerpos[X1DN]=-100;
    enerpos[X1UP]=-100;
    enerpos[X2DN]=-100;
    enerpos[X2UP]=-100;
    jetedge[OUTERJET]=-100;
  }
  else if((startth<M_PI-thetajet)&&(endth>M_PI-thetajet)){
    enerpos[X1DN]=0;
    enerpos[X1UP]=N1-1;
    // if outer jet edge is on this CPU but not on boundary
    i=0;
    for(j=0;j<=N2;j++){
      coord(i, j, FACE2, X);
      bl_coord(X, &r, &th);
      // look for switch from below to above thetajet at inner theta jet edge
      if(th>M_PI-thetajet){
	enerpos[X2DN]=j-1;
	jetedge[OUTERJET]=j-1;
	break;
      }
    }
    enerpos[X2UP]=N2-1;
  }
  else if((startth>=M_PI-thetajet)&&(endth>=M_PI-thetajet)){
    // if cpu is entirely containe within outer jet
    enerpos[X1DN]=0;
    enerpos[X1UP]=N1-1;
    enerpos[X2DN]=0;
    enerpos[X2UP]=N2-1;
    jetedge[OUTERJET]=-100;
  }
  else{
    trifprintf("problem with OUTERJET setjetflux()\n");
    myexit(1);
  }

  if((enerpos[X1DN]!=-100)&&(mycpupos[1]==0)){
    doflux[X1DN]=0; // or horizoni
    trifprintf("proc: %d doing outer jet flux X1DN\n",myid);
  }
  else doflux[X1DN]=-100;

  if((enerpos[X1DN]!=-100)&&(mycpupos[1]==ncpux1-1)){
    doflux[X1UP]=N1;
    trifprintf("proc: %d doing outer jet flux X1UP\n",myid);
  }
  else doflux[X1UP]=-100;

  if(jetedge[OUTERJET]!=-100){
    doflux[X2DN]=jetedge[OUTERJET];
    trifprintf("proc: %d doing outer jet flux X2DN\n",myid);
  }
  else doflux[X2DN]=-100;

  if(mycpupos[2]==ncpux2-1){
    doflux[X2UP]=N2;
    trifprintf("proc: %d doing outer jet flux X2UP\n",myid);
  }
  else doflux[X2UP]=-100;
  // fluxes are on edges of zone, so 0 and N are on edge fluxes

  DIRLOOP trifprintf("proc: %d %d doflux[%d]=%d enerpos[%d]=%d\n",myid,OUTERJETREGION,dir,doflux[dir],dir,enerpos[dir]);


  return(0);
}


int init_dumps(void)
{
  struct blink * blinkptr;
  struct blink * cpulinkptr;
  int i,numlists,numcells;
  int maxnumcolumns;


  trifprintf("begin: init_dumps\n");


  // now setup the data output/input organization for chunking method for each number of columns
  dnumcolumns[IMAGECOL]=1;
  dnumcolumns[RDUMPCOL]=NPRDUMP;

  // 34 + 3*2*NDIM+2*6 currently
  if(GAMMIEDUMP) dnumcolumns[DUMPCOL]=2*2 + NPRDUMP + 1 + NDIM * 2*2 + 2*2 + 1 // normal 34 terms
		   + 2*NDIM // jcon,jcov 
		   + 2*6;   // fcon,fcov
  else  dnumcolumns[DUMPCOL]=2*3 + NPRDUMP + 1 + NDIM * 2*2 + 2*2 + 1 + 2*NDIM+2*6;    // 36+2*NDIM+2*6 currently

  // 202+4+4*4 currently
  dnumcolumns[GDUMPCOL]=2*3+NDIM*NDIM*NDIM+NPG*NDIM*NDIM*2+NPG+4+4*4;
  // 36+29+8*2+4*2+2+12*2+96*2=339
  dnumcolumns[AVGCOL]=2*3 + 1 + NUMNORMDUMP  // (6+1+29=36)
    + NUMNORMDUMP // |normal terms| (29)
    + NDIM*2 // jcon/jcov (8)
    + NDIM*2 // |jcon|/|jcov| (8)
    + NDIM*2 // massflux/|massflux|
    + NUMOTHER*2 // other stuff and fabs of each
    +6*2 // fcon/fcov (12)
    +6*2 // |fcon|,|fcov| (12)
    +7*16 // Tud all 7 parts, all 4x4 terms (112)
    +7*16 // |Tud| all 7 parts, all 4x4 terms (112)
    ;
  if(DOAVG2){
    dnumcolumns[AVGCOL]-=224;
    dnumcolumns[AVG2COL]=7+224; // otherwise doesn't exist so don't need to set
  }
  else dnumcolumns[AVG2COL]=0;
  
  if(DODEBUG){
    dnumcolumns[DEBUGCOL]=NUMFAILFLOORFLAGS*NUMTSCALES;
  }
  else dnumcolumns[DEBUGCOL]=0;

  if(DOFIELDLINE){
    dnumcolumns[FIELDLINECOL]=NUMFIELDLINEQUANTITIES;
  }
  else dnumcolumns[FIELDLINECOL]=0;

  

  trifprintf("dump number of columns\n\n0=IMAGE, 1=RDUMP, 2=DUMP, 3=GDUMP, 4=AVG, 5=AVG2, 6=DEBUG 7=FIELDLINE (see global.h)\n");
  for(i=0;i<NUMDUMPTYPES;i++) trifprintf("dnumcolumns[%d]=%d\n",i,dnumcolumns[i]);

  // setup number of buffers
  maxnumcolumns=0;
  for(i=0;i<NUMDUMPTYPES;i++){
    if(maxnumcolumns<dnumcolumns[i]) maxnumcolumns=dnumcolumns[i];
  }
  // buffer must at least hold maxcolumns of data, and since buffer is only N1*N2*N3 big, make sure that at least NUMBUFFERS*N1*N2*N3>maxnumcolumns
  if(N1*N2*N3<maxnumcolumns) NUMBUFFERS=(int)ceil((FTYPE)maxnumcolumns/((FTYPE)(N1*N2*N3)));
  else NUMBUFFERS=1;


  for(i=0;i<NUMDUMPTYPES;i++) if(dnumcolumns[i]>0) setuplinklist(dnumcolumns[i],i);


  trifprintf("end setuplinklists: %d\n",NUMDUMPTYPES);


  trifprintf("start per cpu lists\n");
  // check link lists
  for(i=0;i<NUMDUMPTYPES;i++){
    if(dnumcolumns[i]>0){
      fprintf(log_file,"i=%d\n",i); fflush(log_file);
      blinkptr=blinkptr0[i];
      numlists=0;
      numcells=0;
      while(blinkptr!=NULL){
	numcells+=blinkptr->num;
	//      fprintf(log_file,"i=%d num=%d, numtotal=%d\n",i,blinkptr->num,numcells); fflush(log_file);
	numlists++;
	blinkptr=blinkptr->np; // next one
      }
      fprintf(log_file,"i=%d numlists=%d numcells=%d\n",i,numlists,numcells);
      numlists=0;
    }
  }

  // check cpu=0 link list
  if(myid==0){
    trifprintf("start cpu==0 lists\n");
    for(i=0;i<NUMDUMPTYPES;i++){
      if(dnumcolumns[i]>0){
	fprintf(log_file,"i=%d\n",i); fflush(log_file);
	cpulinkptr=cpulinkptr0[i];
	numlists=0;
	numcells=0;
	while(cpulinkptr!=NULL){
	  numcells+=cpulinkptr->num;
	  //	fprintf(log_file,"i=%d num=%d, cpu=%d, li=%d, lj=%d, lk=%d, col=%d, numtotal=%d\n",i,cpulinkptr->num,cpulinkptr->cpu,cpulinkptr->i,cpulinkptr->j,cpulinkptr->k,cpulinkptr->col,numcells); fflush(log_file);
	  numlists++;
	  cpulinkptr=cpulinkptr->np; // next one
	}
	fprintf(log_file,"i=%d numlists=%d numcells=%d\n",i,numlists,numcells);
	numlists=0;
      }
    }
  }


  trifprintf("end: init_dumps\n");


  return(0);
}


int setuplinklist(int numcolumns,int which)
{
  int gcount,lcount,numlinks;
  int i,j,k,col,li,lj,lk,pi,pj,pk,pid,firstlink;
  struct blink * clinkptr0, *clinkptr;
  struct blink * linkptr0for0, *linkptrfor0;
  int *lcountfor0;
  int firstlinkfor0;
  int *firstlijk,*li0,*lj0,*lk0,*lcol0;
  int ri,rj,rk,rcol;
  int *cpulist0;
  int numcpusinlist0,lcpu,itercpu,buffersize;
  int maxnumsize;


  maxnumsize=(int)(ceil(ceil((FTYPE)(N1*N2*NUMBUFFERS)/(FTYPE)numcolumns)*(FTYPE)(numcolumns)));

  if(myid==0){
    // cpulist0's size is maximum possible number of cpus in a list due to buffer size
    //    buffersize=(int)(ceil(ceil((FTYPE)(N1*N2*NUMBUFFERS)/(FTYPE)numcolumns)*(FTYPE)(numcolumns)/(FTYPE)N1));
    buffersize=numprocs;
    fprintf(stderr,"max cpus in a list=%d\n",buffersize); fflush(stderr);
    if((cpulist0=(int*)malloc(sizeof(int)*buffersize))==NULL){
      dualfprintf(fail_file,"can't allocate cpulist0\n");
      myexit(10000);
    }
    if((lcountfor0=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate lcountfor0\n");
      myexit(10000);
    }
    if((firstlijk=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate firstlijk\n");
      myexit(10000);
    }
    if((li0=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate li0\n");
      myexit(10000);
    }
    if((lj0=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate lj0\n");
      myexit(10000);
    }
    if((lk0=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate lk0\n");
      myexit(10000);
    }
    if((lcol0=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate lcol0\n");
      myexit(10000);
    }
    for(i=0;i<buffersize;i++){
      cpulist0[i]=0;
    }
    for(i=0;i<numprocs;i++){
      lcountfor0[i]=firstlijk[i]=li0[i]=lj0[i]=lk0[i]=lcol0[i]=0;
    }
  }



  numcpusinlist0=0;

  clinkptr0=NULL;
  gcount=0;
  lcount=0;
  numlinks=0;
  firstlink=1;
  if(myid==0){
    for(itercpu=0;itercpu<numprocs;itercpu++){  firstlijk[itercpu]=1; }
    linkptr0for0=NULL;
    firstlinkfor0=1;
  }

  /////////////////////////
  // general loop
  for(k=0;k<ncpux3*N3;k++)  for(j=0;j<ncpux2*N2;j++)  for(i=0;i<ncpux1*N1;i++) for(col=0;col<numcolumns;col++){
    // relative local index
    li=i%N1;
    lj=j%N2;
    lk=k%N3;
    // cpu position number
    pi=(int)(i/N1);
    pj=(int)(j/N2);
    pk=(int)(k/N3);
    // cpu id for this data
    pid=pk*ncpux2*ncpux1+pj*ncpux1+pi;
    if(myid==pid) lcount++;
    if(myid==0){
      lcountfor0[pid]++;
      // below is if we have need this cpu's data (pid) and need to mark starting point on full grid
      if(firstlijk[pid]){
	cpulist0[numcpusinlist0++]=pid;
	li0[pid]=i;
	lj0[pid]=j;
	lk0[pid]=k;
	lcol0[pid]=col;
	if(col!=0){
	  dualfprintf(fail_file,"col!=0 col=%d, so chunking bad\n",col);
	  myexit(10000);
	}
	firstlijk[pid]=0;
      }
    }
    gcount++;
    //    if(myid==0){
    //  fprintf(fail_file,"%d %d %d %d\n",numcpusinlist0,gcount,pid,cpulist0[numcpusinlist0]); fflush(fail_file);
    // }
    //    fprintf(log_file,"%d %d %d %d %d %d %d %d\n",li,lj,lk,pi,pj,pk,pid,lcount,gcount); fflush(log_file);
    // 1st below if is to catch every buffer amount, while 2nd if part is needed to account for when the number of buffers is such that the last buffer isn't completely needed
    // this should work for any numcolumns or NUMBUFFERS, even at very last zone no matter what
    // chunk in minimum size of numcolumns
    if((gcount%((int)ceil(N1*N2*N3*NUMBUFFERS/numcolumns)*numcolumns)==0)||(gcount==totalzones*numcolumns)){
      // ok, so numcolumns can't exceed the buffer size, highly unlikely to happen, and checked for!
      if(myid==0){
	// must do in order determined to have data, not numerical order
	for(itercpu=0;itercpu<numcpusinlist0;itercpu++){
	  lcpu=cpulist0[itercpu];
	  if(lcountfor0[lcpu]>0){
	    if(itercpu==0){ // first cpu in list
	      ri=li0[lcpu];
	      rj=lj0[lcpu];
	      rk=lk0[lcpu];
	      rcol=lcol0[lcpu];
	    }
	    if(firstlinkfor0){
	      linkptrfor0=linkptr0for0=addlink(NULL);
	      firstlinkfor0=0;
	    }
	    else{
	      linkptrfor0=addlink(linkptrfor0);
	    }
	    linkptrfor0->cpu=lcpu;
	    linkptrfor0->num=lcountfor0[lcpu];
	    linkptrfor0->i=li0[lcpu];
	    linkptrfor0->j=lj0[lcpu];
	    linkptrfor0->k=lk0[lcpu];
	    linkptrfor0->col=lcol0[lcpu];
	    linkptrfor0->ri=ri;
	    linkptrfor0->rj=rj;
	    linkptrfor0->rk=rk;
	    linkptrfor0->rcol=rcol;
	    linkptrfor0->end=0;
	    
	    lcountfor0[lcpu]=0; // reset counter for this id
	    firstlijk[lcpu]=1; // reset starting value
	  }
	  else{
	    fprintf(fail_file,"wtf: shoudn't be here\n");
	    myexit(10000);
	  }
	}
	// the last link is here identified as the last in the series of cpus to communicate with.  There's at least one new link here!
	linkptrfor0->end=1;
	numcpusinlist0=0; // reset list of cpus for this list
      }
      if(lcount>0){
	fprintf(log_file,"numcolumns=%d lcount=%d\n",numcolumns,lcount); fflush(log_file);
        // initialize another structure
        // set previous structure value to this structure, set this next one to NULL
        if(firstlink){
          clinkptr=clinkptr0=addlink(NULL);
	  clinkptr->num=lcount;
          firstlink=0;
        }
        else{
          clinkptr=addlink(clinkptr);
	  clinkptr->num=lcount;
        }
        lcount=0;
      }
    }// otherwise continue
  }      // now we have a link list for each cpu that determines how long each next buffer is that needs to be sent to cpu=0
  blinkptr0[which]=clinkptr0;
  cpulinkptr0[which]=linkptr0for0;

  return(0);
}

// add link for forward-only link list
struct blink * addlink(struct blink * clinkptr)
{
  struct blink *p;

  p=(struct blink *)malloc(sizeof(struct blink));
  p->np=NULL; // terminate list
  // set last link's pointer to this new structure
  if(clinkptr!=NULL) clinkptr->np=p;

  return(p);
}

