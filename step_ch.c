
/**
 *
 * this contains the generic piece of code for advancing
 * the primitive variables 
 *
**/


#include "decs.h"
#include "step_ch.h"

/** end algorithmic choices **/


int step_ch()
{
  int step_ch_simplempi();
  int step_ch_supermpi();


#if(SIMULBCCALC==-1)
  MYFUN(step_ch_simplempi(),"step_ch.c:step_ch()", "step_ch_simplempi()", 1);
#else
  MYFUN(step_ch_supermpi(),"step_ch.c:step_ch()", "step_ch_supermpi()", 1);
#endif

  /* done! */
  return (0);
}




int step_ch_simplempi()
{
  int boundstage;
  SFTYPE mydt;
  int stage, stagei,stagef;
  FTYPE ndt,lastndt;
  FTYPE (*pi)[N2M][NPR];
  FTYPE (*pb)[N2M][NPR];
  FTYPE (*pf)[N2M][NPR];
  FTYPE (*prevpf)[N2M][NPR];
  FTYPE (*pii[4])[N2M][NPR];
  FTYPE (*pbb[4])[N2M][NPR];
  FTYPE (*pff[4])[N2M][NPR];
  //  FTYPE alphaik[MAXSTAGES][MAXSTAGES],betaik[MAXSTAGES];
  FTYPE Cunew[MAXTIMEORDER][4],CUf[MAXTIMEORDER][3];
  int i, j, k;
  int numtimeorders;
  int timeorder;
  SFTYPE dt4diag;
  int finalstep;
  extern int diag_flux(FTYPE prim[][N2M][NPR], FTYPE F1[][N2M][NPR], FTYPE F2[][N2M][NPR],SFTYPE Dt);
  void setup_rktimestep(int *numtimeorders,
			FTYPE (*p)[N2M][NPR],FTYPE (*pk)[N1M][N2M][NPR],
			FTYPE (*pii[4])[N2M][NPR],FTYPE (*pbb[4])[N2M][NPR],FTYPE (*pff[4])[N2M][NPR],
			FTYPE (*Cunew)[4],FTYPE (*CUf)[3]);


  if((WHICHCURRENTCALC==0)||(WHICHCURRENTCALC==2)){
    // for this method, all J components are located back in time
    current_doprecalc(1,p);
    current_doprecalc(2,p);
    current_doprecalc(3,p);
  }


  // setup time-stepping
  setup_rktimestep(&numtimeorders,p,pk,pii,pbb,pff,Cunew,CUf);

  // global debug tracking var
  numstepparts=numtimeorders;



  // general temporal loop
  for(timeorder=0;timeorder<numtimeorders;timeorder++){
  // global debug tracking var
    steppart=timeorder;
    //    if(timeorder==numtimeorders-1) dt4diag=dt; else dt4diag=-1.0; // set dt for accounting diagnostics
    if(timeorder==numtimeorders-1) finalstep=1; else finalstep=0;

#if(PRODUCTION==0)
    trifprintf("|%ds",timeorder);
#endif

    // force-free doesn't use pre_fixup really, currently.
    //#if(EOMTYPE==EOMGRMHD)
    // prefixup
    MYFUN(pre_fixup(-1, pii[timeorder]),"step_ch.c:advance()", "pre_fixup()", 1);
    //#endif

    // advance
    MYFUN(advance(-1,pii[timeorder], ulast, pbb[timeorder], CUf[timeorder], pff[timeorder], Cunew[timeorder],unew,timeorder,numtimeorders,&ndt),"step_ch.c:step_ch()", "advance()", 1);


    // force-free doesn't use post_fixup
    //#if(EOMTYPE==EOMGRMHD)
    // if CHECKSOLUTION==1, then need values to be bounded right now, since use them to check whether even good solutions are really good.
    // post_fixup() will use previous time step pff boundary values to fixup_utoprim() if this is not called.
    // bound advanced values before post_fixup() so fixup_utoprim() has updated boundary values to base fixup on.
    bound_prim(-1,pff[timeorder]);

    // done when all timeorders are completed, so stencil used doesn't matter

    // postfixup
    // post_fixup: If this modifies values and then uses the modified values for other modified values, then must bound_prim() after this
    // if one doesn't care about MPI being same as non-MPI, then can move bound_prim() above to below error_check() and then remove prc<-pv in fixup_utoprim()
    MYFUN(post_fixup(-1,pff[timeorder],finalstep),"step_ch.c:advance()", "post_fixup()", 1);
#if(PRODUCTION==0)
    trifprintf( "x");
    //#endif


#endif
    // must check before MPI operation (since asymmetries would
    // desynchronize cpus
    MYFUN(error_check(1),"step_ch.c", "error_check", 1);

    // bound final values (comes after post_fixup() since changes made by post_fixup)
    //#if(MPIEQUALNONMPI)
    bound_prim(-1,pff[timeorder]);
    //#endif

    if((WHICHCURRENTCALC==0)||(WHICHCURRENTCALC==2)){
      // puts J at the time center, but hard to know if RK is at mid point in time except for midpoint method
      // compute current_doprecalc if near half-step in time
      if(
	 ((numtimeorders>=3)&&(timeorder==1))
	 ||((numtimeorders<=2)&&(timeorder==0))
	 )
	current_doprecalc(0,pff[timeorder]); // should be called using half-time step data
    }
  }
  
  /////////////////
  //
  // FOR ANY TIME ORDER
  //
  //////////////////

  // compute final faradays
  if(WHICHCURRENTCALC==1){
    // compute faraday components needed for time centering of J
    current_doprecalc(0,p);
    // J is located at this time
    current_doprecalc(1,p);
    current_doprecalc(2,p);
    current_doprecalc(3,p); // for diagnostics
  }
  // compute current
  current_calc(cfaraday);

  // compute flux diagnostics (uses global F1/F2)
  // this doesn't exactly make conservation work -- should have in middle step point using full step.  But if PARA, no middle point that's exact.
  // think about where to put this
  // GODMARK
  diag_flux(p,F1, F2,dt); // should use REAL dt, not within a timeorderd RK integration step

#if(PRODUCTION==0)
  trifprintf( "d");
#endif

  /* check timestep
     if (dt < 1.e-9) {
     trifprintf( "timestep too small\n");
     myexit(11);
     }
  */

  // increment time
  t += dt;
  realnstep++;

  // set next timestep
  // find global minimum value of ndt over all cpus
  mpifmin(&ndt);
  if (ndt > SAFE * dt)    ndt = SAFE * dt;
  dt = ndt;
  // don't step beyond end of run
  if (t + dt > tf) dt = tf - t;


  /* done! */
  return (0);
}




int step_ch_supermpi()
{
  int boundstage;
  SFTYPE mydt;
  int stage, stagei,stagef;
  int timeorder;
  FTYPE ndt,lastndt;
  FTYPE (*pi)[N2M][NPR];
  FTYPE (*pb)[N2M][NPR];
  FTYPE (*pf)[N2M][NPR];
  FTYPE (*prevpf)[N2M][NPR];
  FTYPE (*pii[4])[N2M][NPR];
  FTYPE (*pbb[4])[N2M][NPR];
  FTYPE (*pff[4])[N2M][NPR];
  //  FTYPE alphaik[MAXSTAGES][MAXSTAGES],betaik[MAXSTAGES];
  FTYPE Cunew[MAXTIMEORDER][4],CUf[MAXTIMEORDER][3];
  int i, j, k;
  int numtimeorders;
  int stagenow;
  SFTYPE dt4diag;
  int finalstep;
  extern int diag_flux(FTYPE prim[][N2M][NPR], FTYPE F1[][N2M][NPR], FTYPE F2[][N2M][NPR],SFTYPE Dt);
  void setup_rktimestep(int *numtimeorders,
			FTYPE (*p)[N2M][NPR],FTYPE (*pk)[N1M][N2M][NPR],
			FTYPE (*pii[4])[N2M][NPR],FTYPE (*pbb[4])[N2M][NPR],FTYPE (*pff[4])[N2M][NPR],
			FTYPE (*Cunew)[4],FTYPE (*CUf)[3]);



  // setup time-stepping
  setup_rktimestep(&numtimeorders,p,pk,pii,pbb,pff,Cunew,CUf);


  // SPECIAL BOUNDARY/COMPUTATION MPI METHOD (out of date, and doesn't yet work right even if essentially complete code)
  /* check timestep */
  if (dt < MINDT) {
    trifprintf( "timestep too small\n");
    myexit(11);
  }
  
  lastndt=1E30; // initialize lastndt
  for(timeorder=1;timeorder<=numtimeorders;timeorder++){
#if(PRODUCTION==0)
    trifprintf("-to%d/%d-",timeorder,numtimeorders);
#endif
    if(numtimeorders==2){
      // note that pb (used for flux calc which has a stencil
      // calculation) must be different from pf so new stencils in
      // different stages won't affect stencil calculations -- must
      // use old values, not new from most previous temporary stage
      //
      // pi however can be the same as pf since each pi is replaced 1
      // zone at a time with a 0 stencil.
      if(timeorder==1){
	pi=p;
	pb=p;
	pf=pk[0]; // different already, so good for simulbccalc
	prevpf=p; // previous final true array
	mydt=0.5*dt;
      }
      else if(timeorder==2){
	pi=p;
	pb=pk[0];
	pf=p;
	prevpf=pk[0];
	mydt=dt;
      }
    }
    else if(numtimeorders==1){
      pi=p;
      pb=p;
      if(SIMULBCCALC<=0) pf=p; else pf=pk[0]; // need to be different if doing simulbccalc
      prevpf=p;
      mydt=dt;
    }
    if(SIMULBCCALC<=0){ stagei=STAGEM1; stagef=STAGEM1; }
    else if(SIMULBCCALC==1) { stagei=STAGE0; stagef=STAGE2;}
    else if(SIMULBCCALC==2) { stagei=STAGE0; stagef=STAGE5;}

    
    // initialize bound stage
    if(SIMULBCCALC) boundstage=STAGE0;
    else boundstage=STAGEM1;
    for(stage=stagei;stage<=stagef;stage++){
#if(PRODUCTION==0)
      if(SIMULBCCALC) trifprintf("!s%d!",stage);
#endif
      // setup stage loop
#if(SIMULBCCALC==2)
#if(TYPE2==1)
      // GODMARK: isf1, etc. are NOT defined?!
	STAGECONDITION(0,N1-1,0,N2-1,isc,iec,jsc,jec);
	STAGECONDITION(0,N1,-1,N2,isf1,ief1,jsf1,jef1);
	STAGECONDITION(-1,N1,0,N2,isf2,ief2,jsf2,jef2);
	STAGECONDITION(0,N1,0,N2,ise,iee,jse,jee);
	STAGECONDITION(0,N1,0,N2-1,isf1ct,ief1ct,jsf1ct,jef1ct);
	STAGECONDITION(0,N1-1,0,N2,isf2ct,ief2ct,jsf2ct,jef2ct);
	STAGECONDITION(-1,N1,-1,N2,isdq,iedq,jsdq,jedq);
	STAGECONDITION(-2,N1+1,-2,N2+1,ispdq,iepdq,jspdq,jepdq);
	// GODMARK : probably not right for general boundary condition size
#endif
#endif

      // only bounding if safe zones, unsafe needs bz complete
      if(stage<STAGE2){
	bound_prim(boundstage, prevpf);
	if(stage!=STAGEM1) boundstage++;
      }

      // done here instead of local since pseudo-complicated calculation that might slow the dq calculation if done locally per zone
      MYFUN(pre_fixup(stage, prevpf),"step_ch.c:advance()", "pre_fixup()", 1);

      // go from previous solution to new solution
      partialstep=timeorder;      
      // not right for numtimeorders==4 // GODMARK
      MYFUN(advance(-1,pii[timeorder], ulast, pbb[timeorder], CUf[timeorder], pff[timeorder], Cunew[timeorder],unew,timeorder,numtimeorders,&ndt),"step_ch.c:step_ch()", "advance()", 1);
      //      MYFUN(advance(stage, pi, pb, mydt, pf, 0.0, unew, stage, numtimeorders, &ndt),"step_ch.c:step_ch()", "advance()", 1);

      // must check before MPI operation (since asymmetries would desynchronize cpus)
      if(stage<STAGE2){
	MYFUN(error_check(1),"step_ch.c", "error_check", 1);
      }
      if(stage!=STAGEM1){
	if(stage<STAGE2){
	  bound_prim(boundstage, prevpf);
	  boundstage++;
	}
      }
      if(timeorder==numtimeorders){
	if(ndt>lastndt) ndt=lastndt; // don't change if last was lower
	else lastndt=ndt; // new is lower, keep it
      }
    }
    if(timeorder==numtimeorders){// only do on full step
      // find global minimum value of ndt over all cpus
      mpifmin(&ndt);
    }
    // done when all stages are completed, so stencil used doesn't matter
    MYFUN(post_fixup(-1,pf,dt),"step_ch.c:advance()", "post_fixup()", 1);
  }
  /* evaluate diagnostics based on fluxes on second pass (Dt=dt)*/
  // need to do this every timestep, but only after all stages are complete and on full timestep
  diag_flux(p,F1, F2,dt); // should use REAL dt, not within a timeorderd RK integration step

  // copy the contents to the final working array
  if((numtimeorders==1)&&(SIMULBCCALC)) FULLLOOP PLOOP p[i][j][k]=pf[i][j][k];
  
  
  /* increment time */
  t += dt;
  realnstep++;
  
  // new timestep
  if (ndt > SAFE * dt)    ndt = SAFE * dt;
  dt = ndt;    
  /* but don't step beyond end of run */
  if (t + dt > tf)    dt = tf - t;

  

  /* done! */
  return (0);
}




// for the ith stage:

// Uf^i = ulast^i = CUf^{i0} Ui^i + CUf^{i1} ulast^i + CUf^{i2} dU^i

// unew^i = Cunew^{i0} Ui^i + Cunew^{i1} dU^i + Cunew^{i2} Uf^i
void setup_rktimestep(int *numtimeorders,
		      FTYPE (*p)[N2M][NPR],FTYPE (*pk)[N1M][N2M][NPR],
		      FTYPE (*pii[4])[N2M][NPR],FTYPE (*pbb[4])[N2M][NPR],FTYPE (*pff[4])[N2M][NPR],
		      FTYPE (*Cunew)[4],FTYPE (*CUf)[3])
{


  // to avoid special copying of final pff->p, always use p as final pff
  if(TIMEORDER==4){
    // RK4 stepping
    *numtimeorders=TIMEORDER;

    // Ui ulast dU(pb)
    CUf[0][0]=1.0;  CUf[0][1]=0.0;      CUf[0][2]=0.5;
    CUf[1][0]=1.0;  CUf[1][1]=0.0;      CUf[1][2]=0.5;
    CUf[2][0]=1.0;  CUf[2][1]=0.0;      CUf[2][2]=1.0;
    CUf[3][0]=1.0;  CUf[3][1]=0.0;      CUf[3][2]=1.0;

    // Ui dU(Ub) Uf
    Cunew[0][0]=1.0;  Cunew[0][1]=1.0/6.0;      Cunew[0][2]=0.0;
    Cunew[1][0]=0.0;  Cunew[1][1]=1.0/3.0;      Cunew[1][2]=0.0;
    Cunew[2][0]=0.0;  Cunew[2][1]=1.0/3.0;      Cunew[2][2]=0.0;
    Cunew[3][0]=0.0;  Cunew[3][1]=1.0/6.0;      Cunew[3][2]=0.0;

    pii[0]=p;    pbb[0]=p;       pff[0]=pk[0]; // produces U1
    pii[1]=p;    pbb[1]=pk[0];   pff[1]=pk[1]; // produces U2
    pii[2]=p;    pbb[2]=pk[1];   pff[2]=pk[0]; // produces U3
    pii[3]=p;    pbb[3]=pk[0];   pff[3]=p; // produces U4 (only dU part used)
  }
  else if(TIMEORDER==3){
    // TVD optimal RK3 method as in Shu's report
    *numtimeorders=3;
    
    CUf[0][0]=1.0;      CUf[0][1]=0.0;      CUf[0][2]=1.0;
    CUf[1][0]=3.0/4.0;  CUf[1][1]=1.0/4.0;  CUf[1][2]=1.0/4.0;
    CUf[2][0]=1.0/3.0;  CUf[2][1]=2.0/3.0;  CUf[2][2]=2.0/3.0;
    
    // Ui dU(Ub) Uf
    // unew=U3
    Cunew[0][0]=0.0;   Cunew[0][1]=0.0;      Cunew[0][2]=0.0;
    Cunew[1][0]=0.0;   Cunew[1][1]=0.0;      Cunew[1][2]=0.0;
    Cunew[2][0]=0.0;   Cunew[2][1]=0.0;      Cunew[2][2]=1.0;
    
    pii[0]=p;      pbb[0]=p;       pff[0]=pk[0]; // produces U1
    pii[1]=p;      pbb[1]=pk[0];   pff[1]=pk[1]; // produces U2
    pii[2]=p;      pbb[2]=pk[1];   pff[2]=p; // produces U3
  }
  else if(TIMEORDER==2){
#if(1)
    // midpoint method

    *numtimeorders=TIMEORDER;

    // old unew not used for this method (i.e. [?][1]=0)
    CUf[0][0]=1.0; CUf[0][1]=0.0; CUf[0][2]=0.5;
    CUf[1][0]=1.0; CUf[1][1]=0.0; CUf[1][2]=1.0;

    // Ui dU(Ub) Uf
    // unew=U2
    Cunew[0][0]=0.0;   Cunew[0][1]=0.0;      Cunew[0][2]=0.0;
    Cunew[1][0]=0.0;   Cunew[1][1]=0.0;      Cunew[1][2]=1.0;

    pii[0]=p;    pbb[0]=p;       pff[0]=pk[0];
    pii[1]=p;    pbb[1]=pk[0];   pff[1]=p;

#else
    *numtimeorders=TIMEORDER;
    // TVD RK2 (Chi-Wang Shu 1997 - eq 4.10)
    // actually less robust than generic midpoint method above

    CUf[0][0]=1.0; CUf[0][1]=0.0; CUf[0][2]=1.0;
    CUf[1][0]=0.5; CUf[1][1]=0.5; CUf[1][2]=0.5;

    // Ui dU(Ub) Uf
    // unew=U2
    Cunew[0][0]=0.0;   Cunew[0][1]=0.0;      Cunew[0][2]=0.0;
    Cunew[1][0]=0.0;   Cunew[1][1]=0.0;      Cunew[1][2]=1.0;

    pii[0]=p;    pbb[0]=p;       pff[0]=pk[0];
    pii[1]=p;    pbb[1]=pk[0];   pff[1]=p;
#endif
  }
  else if(TIMEORDER==1){
    // Euler method
    *numtimeorders=TIMEORDER;

    CUf[0][0]=1.0; CUf[0][1]=0.0; CUf[0][2]=1.0;

    // Ui dU(Ub) Uf
    // unew=U1
    Cunew[0][0]=0.0;   Cunew[0][1]=0.0;      Cunew[0][2]=1.0;

    pii[0]=p;    pbb[0]=p;    pff[0]=p;
  }

}







int Utoprimgen(int evolvetype, int inputtype,FTYPE *U,  struct of_geom *ptrgeom, FTYPE *pr)
{
  // debug
  int i, j, k;
  FTYPE Ugeomfree[NPR],Ugeomfree0[NPR];
  FTYPE pr0[NPR];
  FTYPE prother[NPR];
  int whichentropy;
  extern void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);
  extern int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr);
  int Utoprimgen_pick(int which, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, FTYPE *pr);
  int Utoprimgen_compare(int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, FTYPE *pr);
  int Utoprimgen_tryagain(int parameter, FTYPE *Ugeomfree0, FTYPE *Ugeomfree,struct of_geom *ptrgeom, FTYPE *pr0, FTYPE *pr);



  // Notice that reconstruct "standard" geometry-free conserved quantities by first modifying the geometry prefactor and THEN dealing with rest-mass or other subtractions and additions.
  // This is consistent with primtoflux() and how primtoflux() is used in source_conn().
  UtoU(inputtype,UNOTHING,ptrgeom,U,Ugeomfree);

  // backup
  PLOOP{
    Ugeomfree0[k]=Ugeomfree[k];
    pr0[k]=pr[k];
  }

#if(REMOVERESTMASSFROMUU==2)
  if(UTOPRIMVERSION!=UTOPRIM5D1){
    dualfprintf(fail_file,"This code does not handle REMOVERESTMASSFROMUU==2: other Utoprim's besides 5D1\n");
    myexit(1);
  }
#endif

  if(EOMTYPE==EOMGRMHD){

    if(DOENTROPY!=DOEVOLVEDIRECTENTROPY){
    
      if(UTOPRIMVERSION!=UTOPRIMCOMPARE) Utoprimgen_pick(UTOPRIMVERSION, EVOLVENOENTROPY, Ugeomfree, ptrgeom, pr);
      else Utoprimgen_compare(EVOLVENOENTROPY,Ugeomfree,ptrgeom,pr);

      // try other methods (assumes all methods can handle WHICHVEL, etc. used for primary model)
      // right now, all can handle WHICHVEL==VELREL4 and energy equation evolution and REMOVERESTMASSFROMUU=0,1
#if((WHICHVEL==VELREL4)&&(REMOVERESTMASSFROMUU<=1)&&(UTOPRIMTRYAGAIN))
      Utoprimgen_tryagain(EVOLVENOENTROPY, Ugeomfree0, Ugeomfree, ptrgeom, pr0, pr);
#endif
    }


    if(DOENTROPY!=DONOENTROPY){
      // also do separate inversion for entropy version of EOMs
      // only one inversion is setup to handle this
      PLOOP prother[k]=pr[k]; // in any case this is a good guess
      if(DOENTROPY==DOEVOLVECOMPAREENTROPY) whichentropy=WHICHENTROPYEVOLVE;
      else if(DOENTROPY==DOEVOLVEDIRECTENTROPY) whichentropy=EVOLVEFULLENTROPY;

      MYFUN(Utoprim(whichentropy,Ugeomfree, ptrgeom, prother),"step_ch.c:Utoprimgen()", "Utoprim", 1);
      // result now contains internal energy (prother[UU,ENTROPY]) as derived from entropy evolution
      // notice that all other quantities are could also be different (if doentropy==evolvefullentropy), hence the prother variable for temporary storage.

      // at this point, ie version of entropy
      if(DOENTROPY==DOEVOLVEDIRECTENTROPY) PLOOP pr[k]=prother[k]; // full evolution
      else if(DOENTROPY==DOEVOLVECOMPAREENTROPY){
	// just overwrite entropy primitive, leave rest same as from full energy equation
	pr[ENTROPY]=prother[ENTROPY];
	// now compare pr[UU] and pr[ENTROPY] with some kind of diagnostic?
	if(evolvetype==EVOLVEUTOPRIM){
	  // then during evolution and pr[UU]-pr[ENTROPY] is relevant to physical calculation
	  // store difference

	  // use difference to compute something
	}
      }

    }

  }
  else if(EOMTYPE==EOMFFDE){
    MYFUN(Utoprim_ffde(Ugeomfree, ptrgeom, pr),"step_ch.c:Utoprimgen()", "Utoprim_ffde", 1);
    //Utoprim_ffde(Ugeomfree, ptrgeom, pr);
  }

  
  
  return(0);
}

//////////////////////////////
//
// COMPARISON of 2 (currently only 2) methods
int Utoprimgen_compare(int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, FTYPE *pr)
{
  int Utoprimgen_pick(int which, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, FTYPE *pr);
  int k;
  FTYPE Uo1[NPR], Uo2[NPR];
  FTYPE ptest1[NPR],ptest2[NPR];
  FTYPE Uo[NPR], po[NPR];
  int test;


  // backup
  PLOOP{
    ptest1[k]=ptest2[k]=pr[k];
    Uo1[k]=Uo2[k]=Ugeomfree[k];
  }

  // currently comparing 5D1 and LDZ
  Utoprimgen_pick(UTOPRIM5D1, EVOLVENOENTROPY, Uo1, ptrgeom, ptest1);
  Utoprimgen_pick(UTOPRIMLDZ, EVOLVENOENTROPY, Uo2, ptrgeom, ptest2);

  PLOOP{
    if(ptest1[k]!=0.0) test=fabs(ptest1[k]-ptest2[k])/ptest1[k]>1E-11;
    else test=(ptest1[k]-ptest2[k])>1E-11;
    if(test){
      dualfprintf(fail_file,"utoprimdiff: %d %d %d  %21.15g   %21.15g   %21.15g %21.15g\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,k,ptest1[k]-ptest2[k],(ptest1[k]!=0.0) ? fabs(ptest1[k]-ptest2[k])/ptest1[k] : ptest1[k]-ptest2[k],ptest1[k],ptest2[k]);
    }
  }
  // always do (use old utoprim)
  PLOOP pr[k]=ptest1[k];

  return(0);
}


// just picks the algorithm to invert
int Utoprimgen_pick(int which, int parameter, FTYPE *Ugeomfree, struct of_geom *ptrgeom, FTYPE *pr)
{
  if(which==UTOPRIM5D1){      MYFUN(Utoprim(parameter,Ugeomfree, ptrgeom, pr),"step_ch.c:Utoprimgen()", "Utoprim", 1);}
  else if(which==UTOPRIMLDZ){ MYFUN(Utoprim_ldz(Ugeomfree, ptrgeom, pr),"step_ch.c:Utoprimgen()", "Utoprim_ldz", 1);}
  else if(which==UTOPRIM2D){ MYFUN(Utoprim_2d(Ugeomfree, ptrgeom, pr),"step_ch.c:Utoprimgen()", "Utoprim_2d", 1);}
  else if(which==UTOPRIM1D){ MYFUN(Utoprim_1d(Ugeomfree, ptrgeom, pr),"step_ch.c:Utoprimgen()", "Utoprim_1d", 1);}
  else if(which==UTOPRIM1DOPT){ MYFUN(Utoprim_1d_opt(Ugeomfree, ptrgeom, pr),"step_ch.c:Utoprimgen()", "Utoprim_1d_opt", 1);}
  else if(which==UTOPRIM1DFINAL){ MYFUN(Utoprim_1d_final(Ugeomfree, ptrgeom, pr),"step_ch.c:Utoprimgen()", "Utoprim_1d_final", 1);}
  else if(which==UTOPRIM2DFINAL){ MYFUN(Utoprim_2d_final(Ugeomfree, ptrgeom, pr),"step_ch.c:Utoprimgen()", "Utoprim_2d_final", 1);}
  else if(which==UTOPRIM5D2){ MYFUN(Utoprim_5d2_final(Ugeomfree, ptrgeom, pr),"step_ch.c:Utoprimgen()", "Utoprim_5d2_final", 1);}
  
  return(0);
}


// Utoprimgen try again.  Can be whatever sequence or number of inversions
int Utoprimgen_tryagain(int parameter, FTYPE *Ugeomfree0, FTYPE *Ugeomfree,struct of_geom *ptrgeom, FTYPE *pr0, FTYPE *pr)
{
  int Utoprimgen_tryagain_substep(int which, int parameter, FTYPE *Ugeomfree0, FTYPE*Ugeomfree, struct of_geom *ptrgeom, FTYPE *pr0, FTYPE *pr);

  Utoprimgen_tryagain_substep(UTOPRIM5D1, parameter, Ugeomfree0, Ugeomfree, ptrgeom, pr0, pr);
  // LDZ as normally used generates nan's
  //  Utoprimgen_tryagain_substep(UTOPRIMLDZ, parameter, Ugeomfree0, Ugeomfree, ptrgeom, pr0, pr);
  Utoprimgen_tryagain_substep(UTOPRIM2D, parameter, Ugeomfree0, Ugeomfree, ptrgeom, pr0, pr);
  Utoprimgen_tryagain_substep(UTOPRIM1D, parameter, Ugeomfree0, Ugeomfree, ptrgeom, pr0, pr);
  Utoprimgen_tryagain_substep(UTOPRIM1DOPT, parameter, Ugeomfree0, Ugeomfree, ptrgeom, pr0, pr);
  Utoprimgen_tryagain_substep(UTOPRIM1DFINAL, parameter, Ugeomfree0, Ugeomfree, ptrgeom, pr0, pr);
  Utoprimgen_tryagain_substep(UTOPRIM2DFINAL, parameter, Ugeomfree0, Ugeomfree, ptrgeom, pr0, pr);
  Utoprimgen_tryagain_substep(UTOPRIM5D2, parameter, Ugeomfree0, Ugeomfree, ptrgeom, pr0, pr);

  // don't restore in the end, leave to whatever state inversion leaves it in.

  return(0);

}

int Utoprimgen_tryagain_substep(int which, int parameter, FTYPE *Ugeomfree0, FTYPE*Ugeomfree, struct of_geom *ptrgeom, FTYPE *pr0, FTYPE *pr)
{
  int k;

  // if really a bad failure that don't want / can't handle, then try again
  if(
     (pflag[ptrgeom->i][ptrgeom->j][FLAGUTOPRIMFAIL]!=0)
     &&(!((pflag[ptrgeom->i][ptrgeom->j][FLAGUTOPRIMFAIL]==UTOPRIMFAILUNEG)&&(STEPOVERNEGU)))
     ){
    // restore
    PLOOP{
      Ugeomfree[k]=Ugeomfree0[k];
      pr[k]=pr0[k];
    }
    // try again
    MYFUN(Utoprimgen_pick(which,parameter,Ugeomfree, ptrgeom, pr),"step_ch.c:Utoprimgen_tryagain_substep()", "Utoprimgen_pick", 1);
  }

  return(0);
}



// there may be something wrong with this function -- didn't work in TIMEORDER==4, had to do standard method
// could have just been that I wasn't bounding after using this
int Utoprimloop(FTYPE unew[][N2M][NPR],FTYPE pf[][N2M][NPR])
{
  struct of_geom geom;
  int i,j;

  ZLOOP{
    get_geometry(i, j, CENT, &geom);
    // invert True U->p
    MYFUN(Utoprimgen(EVOLVEUTOPRIM,UEVOLVE,unew[i][j], &geom, pf[i][j]),"step_ch.c:advance()", "Utoprimgen", 1);
  }
  return(0);
}


int primtoUloop(FTYPE pi[][N2M][NPR],FTYPE unew[][N2M][NPR])
{
  struct of_geom geom;
  struct of_state q;
  int i,j;

  ZLOOP{
    MYFUN(get_state(pi[i][j], &geom, &q),"step_ch.c:primtoUloop()", "get_state()", 1);
    get_geometry(i, j, CENT, &geom);
    // forward calculate U(p)
    MYFUN(primtoU(UEVOLVE,pi[i][j], &q, &geom, unew[i][j]),"step_ch.c:primtoUloop()", "primtoU()", 1);
  }
  return(0);
}



// choose limiter
int choose_limiter(int i, int j, int k)
{
#if(LIMADJUST==LIMITERFIXED)
  return(lim);
#else

#if(HYDROLIMADJUSTONLY)
  if(k<B1) return(pflag[i][j][FLAGREALLIM]);
  else return(lim);
#else
  return(pflag[i][j][FLAGREALLIM]);
#endif

#endif

}

// choose flux form
void choose_flux(int i, int j, int k,FTYPE *laxffrac,FTYPE *hllfrac)
{
#if(FLUXADJUST)

#if(HYDROFLUXADJUSTONLY)
  if(k<B1){
    if(pflag[i][j][FLAGREALFLUX]==HLLFLUX){ hllfrac[k]=1.0; laxffrac[k]=0.0; }
    else { hllfrac[k]=0.0; laxffrac[k]=1.0; }
  }
  else{
    if(fluxmethod==HLLFLUX){
      hllfrac[k]=1.0;
      laxffrac[k]=0.0;
    }
    else if(fluxmethod==LAXFFLUX){
      hllfrac[k]=0.0;
      laxffrac[k]=1.0;
    }
  }
#else
  if(pflag[i][j][FLAGREALFLUX]==HLLFLUX){ hllfrac[k]=1.0; laxffrac[k]=0.0; }
  else { hllfrac[k]=0.0; laxffrac[k]=1.0; }
#endif

#else
  if(fluxmethod==HLLFLUX){
    hllfrac[k]=1.0;
    laxffrac[k]=0.0;
  }
  else if(fluxmethod==LAXFFLUX){
    hllfrac[k]=0.0;
    laxffrac[k]=1.0;
  }
#endif

}


int fluxcalc(int stage, FTYPE pr[][N2M][NPR],
	     FTYPE F[][N2M][NPR], int dir, FTYPE *ndt)
{
  int choose_limiter(int i, int j, int k);
  void choose_flux(int i, int j, int k,FTYPE *laxffrac,FTYPE *hllfrac);
  int cminmax_calc(FTYPE cmin_l,FTYPE cmin_r,FTYPE cmax_l,FTYPE cmax_r,FTYPE *cmin,FTYPE *cmax,FTYPE *ctop);
  extern void slope_lim(int reallim, int ii, int jj, int kk, FTYPE yll, FTYPE yl, FTYPE yc, FTYPE yr, FTYPE yrr,FTYPE *dq,FTYPE *left,FTYPE *right);
  extern int rescale(int which, int dir, FTYPE *pr, struct of_geom *geom,FTYPE*newvar);
  FTYPE (*dq)[N2M][NPR];
  FTYPE (*p2interp)[N2M][NPR];
  int i, j, k, idel, jdel, face;
  FTYPE p_l[NPR], p_r[NPR], F_l[NPR], F_r[NPR], U_l[NPR], U_r[NPR];
  FTYPE p2interp_l[NPR],p2interp_r[NPR];
  FTYPE cmax_l, cmax_r, cmin_l, cmin_r, cmax, cmin, dtij;
  int ignorecourant;
  FTYPE ctop;
  int failreturn;
  struct of_geom geom;
  FTYPE Xcent[NDIM],Xleft[NDIM];
  FTYPE rleft,rcent,thleft,thcent;
  struct of_state state_l, state_r;
  int reallim;
  int is,js;
  SFTYPE Dt=-1.0;
  FTYPE laxffrac[NPR],hllfrac[NPR];
  int locallim;


  if (dir == 1) {
    is=isf1;
    js=jsf1;
    idel = 1;
    jdel = 0;
    face = FACE1;
    dq=dq1;
  }
  else if (dir == 2) {
    is=isf2;
    js=jsf2;
    idel = 0;
    jdel = 1;
    face = FACE2;
    dq=dq2;
  }
  else {
    myexit(10);
  }


  //////////////////////////
  //
  // rescale before interpolation
  //
  ////////////////////////////
#if(RESCALEINTERP)
  p2interp=prc; // it's different
  PREDQZLOOP{
    // get geometry for center pre-interpolated values
    get_geometry(i, j, CENT, &geom); 
    // assume no need for a guess to p2interp to get pr
    rescale(1,dir,pr[i][j],&geom,p2interp[i][j]);
  }
#else
  p2interp=pr; // it's itself
#endif

 /////////////////////////////////////
 //
 // evaluate slopes of primitive variables
 //
 /////////////////////////////////////
  DQZLOOP PLOOP {
    // find interpolated values
    // if limiter is PARA, returns in pleft, pright.  Other limiters return dq unless LIMADJUST>0, then all limiters return left/right only
    slope_lim(choose_limiter(i,j,k),i,j,k,
	      p2interp[i - 2*idel][j - 2*jdel][k],
	      p2interp[i - idel][j - jdel][k],
	      p2interp[i][j][k],
	      p2interp[i + idel][j + jdel][k],
	      p2interp[i + 2*idel][j + 2*jdel][k],
	      &dq[i][j][k],&pleft[i][j][k],&pright[i][j][k]
	      );
  }




  //////////////////////////////////////
  //
  // flux loop
  //
  ////////////////////////////////////////

  *ndt = 1.e9;
#if((SIMULBCCALC==2)&&(TYPE2==1))
  FZLOOP(is,js)
#else
  FZLOOP(-jdel, -idel)
#endif
  {


    //////////////////////////////////
    //
    // this avoids problems (bad fluxes) on the pole
    //
    ///////////////////////////////////
#if(ZEROPOLEFLUX==1)
    if((ZEROPOLEFLUX==1)&&(dir == 2) && ( ((startpos[2]+j == 0)&&(BCtype[X2DN]==POLARAXIS)) || ((startpos[2]+j == totalsize[2])&&(BCtype[X2UP]==POLARAXIS))  )) {
      // designed to just outright avoid the calculation.
      PLOOP F[i][j][k] = 0. ;
    }
    else {
#endif

      //////////////////////////////////////
      //
      // interpolate primitive using slope
      //
      // |=interface
      // i=zone center of ith zone
      //
      // |              |       dq(i)        |
      // |         pl(i)|pr(i)    i          |
      // |              |pleft(i)   pright(i)|
      //
      //
      //
      //////////////////////////////////////
      
      if((HORIZONSUPERFAST)&&((locallim!=PARA)&&(LIMADJUST==0))&&(dir==1)){
	coord(i-1, j, CENT, Xleft);
	coord(i, j, CENT, Xcent);
	bl_coord(Xleft, &rleft, &thleft);
	bl_coord(Xcent, &rcent, &thcent);

	  PLOOP{
	    locallim=choose_limiter(i,j,k);
	    // get interpolated quantity
	    if((locallim!=PARA)&&(LIMADJUST==0)){
	      if(rleft>Rhor) p2interp_l[k] = p2interp[i - idel][j - jdel][k] + 0.5 * dq[i - idel][j - jdel][k];
	      else p2interp_l[k] = p2interp[i][j][k] - 0.5 * dq[i][j][k];
	      if(rcent>Rhor) p2interp_r[k] = p2interp[i][j][k] - 0.5 * dq[i][j][k];
	      else p2interp_r[k] = p2interp[i+idel][j+jdel][k] - 1.5 * dq[i+idel][j+jdel][k];
	    }
	    else{
	      p2interp_l[k] = pright[i-idel][j-jdel][k];
	      p2interp_r[k] = pleft[i][j][k];
	    }
	  } // end PLOOP
      }// if horizonsuperfast
      else{
	PLOOP{
	  locallim=choose_limiter(i,j,k);
	  // get interpolated quantity
	  if((locallim!=PARA)&&(LIMADJUST==0)){
	    p2interp_l[k] = p2interp[i - idel][j - jdel][k] + 0.5 * dq[i - idel][j - jdel][k];
	    p2interp_r[k] = p2interp[i][j][k] - 0.5 * dq[i][j][k];
	  }
	  else{
	    p2interp_l[k] = pright[i-idel][j-jdel][k];
	    p2interp_r[k] = pleft[i][j][k];
	  }
	}
      }

      ////////////////////////
      //
      // get the geometry for the flux face
      //
      ///////////////////////
      get_geometry(i, j, face, &geom);



      /////////////////////////////////////
      //      
      // unrescale after interpolation
      //
      ///////////////////////////////////
#if(RESCALEINTERP)
      // setup plausible p_l/p_r in case used for some reason (e.g. inversion starting guess)
      // this is sufficient for utoprim_1D...a better guess does no better (see interpU code
      PLOOP{
	p_l[k]=pr[i][j][k];
	p_r[k]=pr[i][j][k];
      }
      rescale(-1,dir,p_l,&geom,p2interp_l);
      rescale(-1,dir,p_r,&geom,p2interp_r);
#else
      PLOOP{
	p_l[k]=p2interp_l[k];
	p_r[k]=p2interp_r[k];
      }
#endif


      //////////////////////////////
      //
      // Must preserve divb in 1D Riemann problem, so B^{dir} must be continuous
      //
      ///////////////////////////////

#if(BDIRCONT)
      if(dir==1) p_l[B1]=p_r[B1]=0.5*(p_l[B1]+p_r[B1]);
      else if(dir==2) p_l[B2]=p_r[B2]=0.5*(p_l[B2]+p_r[B2]);
#endif

      ///////////////////////
      //
      // correct interpolated quantities
      // no fixup accounting for these intermediate quantities
      //
      //////////////////////
#if(EVOLVECHECKS)
#if(WHICHVEL==VEL4)
#if(ZEROOUTFLOWFLUX==1)
      inflow_check_4vel(dir,p_l,&geom,-1.0);
      inflow_check_4vel(dir,p_r,&geom,-1.0);
#endif
#elif(WHICHVEL==VEL3)
#if(ZEROOUTFLOWFLUX==1)
      inflow_check_3vel(dir,p_l,&geom,-1.0);
      inflow_check_3vel(dir,p_r,&geom,-1.0);
#endif
#if(JONCHECKS)
      // must verify if this p makes sense (u^t sense)
      MYFUN(check_pr(p_l,pr[i-idel][j-jdel],&geom,1,-1.0),"step_ch.c:fluxcalc()", "check_pr()", 1);	
      
      MYFUN(check_pr(p_r,pr[i][j],&geom,2,-1.0),"step_ch.c:fluxcalc()", "check_pr()", 2);
#endif
#elif(WHICHVEL==VELREL4)
#if(ZEROOUTFLOWFLUX==1)
      inflow_check_rel4vel(dir,p_l,&geom,-1.0);
      inflow_check_rel4vel(dir,p_r,&geom,-1.0);
#endif
      // need to limit gamma since gamma may be large for interpolated value and would lead to bad fluxes
      MYFUN(limit_gamma(GAMMAMAX,p_l,&geom,-1.0),"step_ch.c:fluxcalc()", "limit_gamma()", 1);
      MYFUN(limit_gamma(GAMMAMAX,p_r,&geom,-1.0),"step_ch.c:fluxcalc()", "limit_gamma()", 2);
#endif// end if WHICHVEL==VEL4REL      
#endif



      //////////////////////
      //
      // setup flux calculation based upon interpolated primitives
      //
      //////////////////////
      MYFUN(get_state(p_l, &geom, &state_l),"step_ch.c:fluxcalc()", "get_state()", 1);
      MYFUN(get_state(p_r, &geom, &state_r),"step_ch.c:fluxcalc()", "get_state()", 2);
      
      MYFUN(primtoflux(UEVOLVE,p_l, &state_l, dir, &geom, F_l),"step_ch.c:fluxcalc()",
		      "primtoflux_calc() dir=1/2 l", 1);
      MYFUN(primtoflux(UEVOLVE,p_r, &state_r, dir, &geom, F_r),"step_ch.c:fluxcalc()",
		      "primtoflux_calc() dir=1/2 r", 1);
      
      MYFUN(primtoflux(UEVOLVE,p_l, &state_l, TT, &geom, U_l),"step_ch.c:fluxcalc()", "primtoflux_calc() dir=l0", 1);
      MYFUN(primtoflux(UEVOLVE,p_r, &state_r, TT, &geom, U_r),"step_ch.c:fluxcalc()", "primtoflux_calc() dir=r0", 1);
      

      // characteristic based upon t^n level for 1/2 step and t^{n+1/2} level for the full step.
      MYFUN(vchar(p_l, &state_l, dir, &geom, &cmax_l, &cmin_l,&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 1);
      MYFUN(vchar(p_r, &state_r, dir, &geom, &cmax_r, &cmin_r,&ignorecourant),"step_ch.c:fluxcalc()", "vchar() dir=1or2", 2);


      //////////////////////////////////
      //
      // get wave speeds for flux calculation
      //
      /////////////////////////////////
      cminmax_calc(cmin_l,cmin_r,cmax_l,cmax_r,&cmin,&cmax,&ctop);
      

      //////////////////////////////////
      //
      // decide which flux formula to use
      //
      /////////////////////////////////
      PLOOP choose_flux(i,j,k,laxffrac,hllfrac);


      //////////////////////////////////
      //
      // Decide if using different flux calculation on boundary zones
      //
      /////////////////////////////////
#if(HLLBOUNDARY)
      if(
	 ((dir == 2) && ( ((startpos[2]+j == 0)&&(BCtype[X2DN]==POLARAXIS)) || ((startpos[2]+j == totalsize[2])&&(BCtype[X2UP]==POLARAXIS))  ))||
	 ((dir == 1) && ( ((startpos[1]+i == 0)&&(BCtype[X1DN]==OUTFLOW)) || ((startpos[1]+i == totalsize[1])&&(BCtype[X1UP]==OUTFLOW))  ))
	 )
	{
	  PLOOP{
	    hllfrac[k]=1.0;
	    laxffrac[k]=0.0;
	  }
	}
#endif

      //////////////////////////////////
      //
      // actually compute the flux
      //
      /////////////////////////////////

      if(cmax+cmin!=0.0){
	PLOOP F[i][j][k] =
	  hllfrac[k] * ((cmax * F_l[k] + cmin * F_r[k] - cmax * cmin * (U_r[k] - U_l[k]) ) / (cmax + cmin) )
	  + laxffrac[k] * (0.5 * (F_l[k] + F_r[k] - ctop * (U_r[k] - U_l[k]) ) );
      }
      else{
	PLOOP F[i][j][k] =(0.5 * (F_l[k] + F_r[k] - ctop * (U_r[k] - U_l[k]) ) );
      }
      
      /////////////////////////////
      //
      // evaluate restriction on timestep
      //
      ///////////////////////////////

      cmax = max(cmax, cmin);
      dtij = cour * dx[dir] / cmax;
      if (dtij < *ndt)
	*ndt = dtij;
#if(ZEROPOLEFLUX==1) // close crazy ifdef
    }
#endif
  }

  // no longer needed since acts on different stored quantity, not pr
  //  if(RESCALEINTERP){
    // rescale before interpolation
    //    PREDQZLOOP{
    //  get_geometry(i, j, CENT, &geom);
    //  rescale(-1,dir,pr[i][j],&geom);
    // }
    // }
  
  
  return (0);
}





// GODMARK:

// really HARM is currently using VERY local lax Friedrich.
// maybe try local lax Friedrich, using max wave speed from zones used to reconstruct the zone (most common?)
// also can try more global wave speed, or even speed of light.

// determine cmin,cmax,ctop
int cminmax_calc(FTYPE cmin_l,FTYPE cmin_r,FTYPE cmax_l,FTYPE cmax_r,FTYPE *cmin,FTYPE *cmax,FTYPE *ctop)
{
  FTYPE lmin,lmax,ltop;

  // need to make use of ROE averages
      
      
  // As per the HLLE method, one must use the wave speed at the interface.  The below is arbitrarily choosing the largest of the interpolated states.

  // Apparently one should use Roe's linearisation to compare
  // against, not the left/right interpolations.  Then compare the
  // left most going and right most going Rho linearisation
  // eigenvalues with the left and right numerical approximations
  // and choose the minimum and maximum values as cmin and cmax,
  // all respectively.
  // Then one defines the new cmin as the minumum of cmin and 0, and new cmax as maximum of cmax and 0.
  // Thus this is slightly different and doesn't avoid negative mass/ie densities like HLLE.


  // LAXF here is really the Rusanov flux.  Lax-Friedrichs would say ctop->dx/dt



  lmin=min(cmin_l,cmin_r);
  lmax=max(cmax_l,cmax_r);


#if(HYPERHLL==0)
  // sharp change at c=0
  lmin = fabs(max(0., -lmin));
  lmax = fabs(max(0., lmax));
#else

#define HYPERFUN(c,x0) (0.5*( (c)+(M_SQRT2l-1.0)*sqrt(2.0*(M_SQRT2l+1.0)*(x0)*(x0)+(3.0+2.0*M_SQRT2l)*(c)*(c)) ) )
  // GODMARK
  // need a good way to set x0 based upon fraction of local light speed
  // or something...

  // function is always positive for given lmin/lmax as original sign
  lmin=HYPERFUN(-lmin,1E-4);
  lmax=HYPERFUN(lmax,1E-4);
#endif

  // returns positive cmin,cmax,ctop
  *cmin=lmin;
  *cmax=lmax;
  *ctop=max(lmax,lmin);


  return(0);
      


}





#define CORNGDETVERSION 1
  // whether to specify gdet at end when setting EMF or to have internal to variables before averaging.
  // point is that gdet at end is probably better, esp. at coordinate singularities.


int flux_ct(int stage, FTYPE pb[][N2M][NPR],FTYPE F1[][N2M][NPR], FTYPE F2[][N2M][NPR])
{
  int i, j, k, l;
  // Below stuff is for Athena 1 and Athena 2
  FTYPE ucon[NDIM];
  struct of_state state;
  int dir;
  struct of_geom geom,geom1,geom2,geom3,geom4;
  FTYPE cmax1,cmin1,cmax2,cmin2,ctop1,ctop2;
  int ignorecourant;
  // Gammie stuff
  FTYPE v1_av,w1,v2_av,w2,w_12 ;
  FTYPE rho_av,sqrtrho,b1_av,b2_av,va1,va2 ;
  FTYPE emfmm,emfpm,emfmp,emfpp,alpha ;
  FTYPE B1pp,B1pm,B1mp,B1mm,B2pp,B2pm,B2mp,B2mm ;
  FTYPE U1pp,U1pm,U1mp,U1mm,U2pp,U2pm,U2mp,U2mm ;
  //  FTYPE cms_func(FTYPE *prim_var) ;
  FTYPE B1d,B1u,B1l,B1r,B2d,B2u,B2l,B2r ;
  FTYPE pbavg[NPR];

  // typically EMF is located at corner (CORN)


  ///////////////////////
  //
  //  COMPUTE PRE-FUNCTIONS used by Athena method
  if((FLUXB==ATHENA1)||(FLUXB==ATHENA2)){
    // compute v^i

    // loop must go over EMFZLOOP's range minus 1 for i and j
    //for(i=-1;i<N1+1;i++)	for(j=-1;j<N2+1;j++) {
    PREEMFZLOOP{
      // use ucon_calc() to get v^i 

      // EMF below is based upon averaging of zone-centered quantities, so use CENT here (i.e. not CORN)
      get_geometry(i, j, CENT, &geom);
      if(ucon_calc(pb[i][j], &geom, ucon) >= 1) return(1);
      // geom.g is \detg for EMF flux (geom.e is EOM factor for flux equation)
#if(CORNGDETVERSION)
      for(l=U1;l<=U3;l++) vconemf[i][j][l]=(ucon[l-U1+1]/ucon[TT]); // put in at end
#else
      for(l=U1;l<=U3;l++) vconemf[i][j][l]=(ucon[l-U1+1]/ucon[TT])*(geom.e[l]);
#endif
    }
  }


  /////////////////////////
  //
  // COMPUTE EMF
  

  if(FLUXB==FLUXCTHLL){
    // do nothing
  }
  else if(FLUXB==FLUXCTTOTH){
    /* calculate EMFs */
    /* Toth approach: just average */

    // Stone & Gardiner point out that not consistent with underlying integration algorithm for plane-parallel, grid-aligned flows.
    // fix is to change 0.25 to 0.5 and use a diffusive term as in ATHENA1

    // flux part is just average of same emf term at 4 different edge locations of (B^2 v^1 - B^1 v^2)
#if(CORNGDETVERSION)
    EMFZLOOP{ 
      get_geometry(i, j,   CORN, &geom);
      get_geometry(i, j,   FACE1, &geom1);
      get_geometry(i, j-1, FACE1, &geom2);
      get_geometry(i, j,   FACE2, &geom3);
      get_geometry(i-1, j, FACE2, &geom4);

      // obviously geom.e[B1] has to be equal to geom.e[B2] for this method
      emf[i][j] =
      0.25 * (geom.e[B2])*(F1[i][j][B2]/(geom1.e[B2]) + F1[i][j - 1][B2]/(geom2.e[B2])
	           - F2[i][j][B1]/(geom3.e[B1]) - F2[i - 1][j][B1]/(geom4.e[B1]));
    }
#else
    EMFZLOOP
      emf[i][j] =
      0.25 * (F1[i][j][B2] + F1[i][j - 1][B2]
	      - F2[i][j][B1] - F2[i - 1][j][B1]);
#endif
  }
  else if(FLUXB==FLUXCD){
    // here emf is actually electric field
    // in Toth 2000, the F1=f^x F2=f^y
    // 0.125 above also takes care of larger differencing used in advance
    EMFZLOOP emf[i][j] =
      0.125 * (-F1[i][j][B2] - F1[i+1][j][B2]
	       + F2[i][j][B1] + F2[i][j+1][B1]);
  }
  else if(FLUXB==ATHENA1){
    /* Stone & Gardiner eq. 39 */
    // Charles Gammie (8/17/05) says this is simple, but repairs HARM defect that 2D does not reduce to 1D when waves are along coordinate lines


    EMFZLOOP{
      // {emf}_i=-\epsilon_{ijk} v^i B^k


      // average of the below results in averaged emf located at corner (CORN)
      // same sign as F's (B^2 v^1 - B^1 v^2) with gdet built into vconemf

      // could remove gdet from vconemf and F1/F2 and put gdet at CORN for final result!  Avoids axis problems?
      // applies to above Toth version as well!
      // GODMARK

      emfmp = pb[i-1][j  ][B2]*vconemf[i-1][j  ][U1] -
	      pb[i-1][j  ][B1]*vconemf[i-1][j  ][U2] ;
      emfmm = pb[i-1][j-1][B2]*vconemf[i-1][j-1][U1] -
	      pb[i-1][j-1][B1]*vconemf[i-1][j-1][U2] ;
      emfpm = pb[i  ][j-1][B2]*vconemf[i  ][j-1][U1] -
	      pb[i  ][j-1][B1]*vconemf[i  ][j-1][U2] ;
      emfpp = pb[i  ][j  ][B2]*vconemf[i  ][j  ][U1] -
	      pb[i  ][j  ][B1]*vconemf[i  ][j  ][U2] ;

      // flux part is just average of same emf term, while here is corrected with another diffusive term.
#if(CORNGDETVERSION)
      get_geometry(i, j,   CORN, &geom);
      get_geometry(i, j,   FACE1, &geom1);
      get_geometry(i, j-1, FACE1, &geom2);
      get_geometry(i, j,   FACE2, &geom3);
      get_geometry(i-1, j, FACE2, &geom4);


      emf[i][j] =
      0.5 * (geom.e[B2])*(F1[i][j][B2]/(geom1.e[B2]) + F1[i][j - 1][B2]/(geom2.e[B2])
	            - F2[i][j][B1]/(geom3.e[B1]) - F2[i - 1][j][B1]/(geom4.e[B1])
                    - 0.25*(emfmp + emfmm + emfpm + emfpp)) ;
#else
      emf[i][j] = 0.5*(F1[i][j][B2] + F1[i][j-1][B2]
    	             - F2[i][j][B1] - F2[i-1][j][B1])
	             - 0.25*(emfmp + emfmm + emfpm + emfpp) ;
#endif
    }
    

  }
  else if(FLUXB==ATHENA2){
    /* Stone & Gardiner eq. 48 */
    // Charles Gammie (8/17/05) says does much better on flux loop advection test than ordinary HARM
    EMFZLOOP{

      // average of the below results in averaged emf located at corner (CORN)
      emfmp = pb[i-1][j  ][B2]*vconemf[i-1][j  ][U1] -
	      pb[i-1][j  ][B1]*vconemf[i-1][j  ][U2] ;
      emfmm = pb[i-1][j-1][B2]*vconemf[i-1][j-1][U1] -
 	      pb[i-1][j-1][B1]*vconemf[i-1][j-1][U2] ;
      emfpm = pb[i  ][j-1][B2]*vconemf[i  ][j-1][U1] -
	      pb[i  ][j-1][B1]*vconemf[i  ][j-1][U2] ;
      emfpp = pb[i  ][j  ][B2]*vconemf[i  ][j  ][U1] -
 	      pb[i  ][j  ][B1]*vconemf[i  ][j  ][U2] ;



      // dq1 and dq2 are well-defined flux fluxcalc() (both directions)
      // simple average of 1-D linear extrapolation here, where could use p2interp data saved from step_ch.c?  still available by this function call?
      B1d = 0.5*(
		 pb[i-1][j-1][B1] + 0.5*dq1[i-1][j-1][B1] +
		 pb[i][j-1][B1] - 0.5*dq1[i][j-1][B1]) ;
      B1u = 0.5*(
		 pb[i-1][j][B1] + 0.5*dq1[i-1][j][B1] +
		 pb[i][j][B1] - 0.5*dq1[i][j][B1]) ;

      B2l = 0.5*(
		 pb[i-1][j-1][B2] + 0.5*dq2[i-1][j-1][B2] +
		 pb[i-1][j][B2] - 0.5*dq2[i-1][j][B2]) ;
      B2r = 0.5*(
		 pb[i][j-1][B2] + 0.5*dq2[i][j-1][B2] +
		 pb[i][j][B2] - 0.5*dq2[i][j][B2]) ;

      B1mm = pb[i-1][j-1][B1] ;
      B1mp = pb[i-1][j][B1] ;
      B1pm = pb[i][j-1][B1] ;
      B1pp = pb[i][j][B1] ;
      B2mm = pb[i-1][j-1][B2] ;
      B2mp = pb[i-1][j][B2] ;
      B2pm = pb[i][j-1][B2] ;
      B2pp = pb[i][j][B2] ;

      // compute characteristic velocity -- only for Athena2 method


      // average pb to CORN for average phase speed there
      PLOOP pbavg[k]=0.25*(pb[i][j][k]+pb[i][j-1][k]+pb[i-1][j][k]+pb[i-1][j-1][k]);
      get_geometry(i, j, CORN, &geom); // used here and below emf's
      MYFUN(get_state(pbavg, &geom, &state),"step_ch.c:flux_ct()", "get_state()", 1);
      dir=1; MYFUN(vchar(pbavg, &state, dir, &geom, &cmax1, &cmin1,&ignorecourant),"step_ch.c:flux_ct()", "vchar() dir=1", 1);
      dir=2; MYFUN(vchar(pbavg, &state, dir, &geom, &cmax2, &cmin2,&ignorecourant),"step_ch.c:flux_ct()", "vchar() dir=2", 2);
      ctop1 = max(fabs(cmax1), fabs(cmin1));
      ctop2 = max(fabs(cmax2), fabs(cmin2));
      //      alpha=0.5*(ctop1+ctop2); // use average?
      alpha=max(ctop1,ctop2); // use maximum?
      // seems alpha can be arbitrary since 0 is ATHENA1

      //      alpha = dx1/dt ;	/* crude approx */



#if(CORNGDETVERSION)
      get_geometry(i, j,   FACE1, &geom1);
      get_geometry(i, j-1, FACE1, &geom2);
      get_geometry(i, j,   FACE2, &geom3);
      get_geometry(i-1, j, FACE2, &geom4);


      emf[i][j] = 0.5 * (geom.e[B2])*(F1[i][j][B2]/(geom1.e[B2]) + F1[i][j - 1][B2]/(geom2.e[B2])
	                        - F2[i][j][B1]/(geom3.e[B1]) - F2[i - 1][j][B1]/(geom4.e[B1])

                     - 0.25*(emfmp + emfmm + emfpm + emfpp)

      	             - 0.125*alpha*(
		         B1d - B1mm - B1u + B1mp
		       + B1d - B1pm - B1u + B1pp
		       + B2r - B2pm - B2l + B2mm
		       + B2r - B2pp - B2l + B2mp)) ;

#else

      emf[i][j] = 0.5*(F1[i][j][B2] + F1[i][j-1][B2]
                     - F2[i][j][B1] - F2[i-1][j][B1]) 

	             - 0.25*(emfmp + emfmm + emfpm + emfpp) 

 	             - 0.125*alpha*geom.e[B1]*(
		         B1d - B1mm - B1u + B1mp
		       + B1d - B1pm - B1u + B1pp
		       + B2r - B2pm - B2l + B2mm
		       + B2r - B2pp - B2l + B2mp) ;
    
#endif

    }// end EMF loop
  }// end if athena2
      

  ////////////////
  //
  // compute flux from EMF


  if(FLUXB==FLUXCTHLL){
    // no change
  }
  else if((FLUXB==FLUXCTTOTH)||(FLUXB==ATHENA2)||(FLUXB==ATHENA1)){
    /* rewrite EMFs as fluxes, after Toth */
    F1CTZLOOP {
      F1CT[i][j][B1]=F1[i][j][B1] = 0.;
      F1CT[i][j][B2]=F1[i][j][B2] = 0.5 * (emf[i][j] + emf[i][j + 1]);
    }
    F2CTZLOOP {
      F2CT[i][j][B1]=F2[i][j][B1] = -0.5 * (emf[i][j] + emf[i + 1][j]);
      F2CT[i][j][B2]=F2[i][j][B2] = 0.;
    }
  }
  else if(FLUXB==FLUXCD){
    F1CTZLOOP {
      F1CT[i][j][B1]=F1[i][j][B1] = 0.0;
      F1CT[i][j][B2]=F1[i][j][B2] = emf[i][j];
    }
    F2CTZLOOP {
      F2CT[i][j][B1]=F2[i][j][B1] = -emf[i][j];
      F2CT[i][j][B2]=F2[i][j][B2] = 0.;
    }
  }


  return(0);
}
