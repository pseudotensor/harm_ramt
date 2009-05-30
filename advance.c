#include "decs.h"
#include "step_ch.h"


// pi: initial values at t=t0 to compute Ui
// pb: values used to compute flux/source
// pf: solution using flux(pb) from pi's Ui -> Uf

// pi, pb, and pf can all be the same since
// 1) pb used first on a stencil, not modified, to compute fluxes
// 2) pf=pi is assigned by value at each zone
// 3) pf is modified using Utoprim at each zone using pb for sources (to correspond to fluxes which used pb)
//
// So in the end only pf is modified at each zone, so the loop changing p at previous (i,j) location doesn't affect the any new location in (i,j)
int advance(int stage,
	    FTYPE pi[][N2M][NPR],
	    FTYPE ulast[][N2M][NPR], 
	    FTYPE pb[][N2M][NPR],
	    FTYPE *CUf, FTYPE pf[][N2M][NPR],
	    FTYPE *Cunew, FTYPE unew[][N2M][NPR],
	    int timeorder, int numtimeorders,
	    FTYPE *ndt)
{
  int advance_standard(int stage, FTYPE pi[][N2M][NPR], FTYPE ulast[][N2M][NPR], FTYPE pb[][N2M][NPR], FTYPE *CUf, FTYPE pf[][N2M][NPR], FTYPE *Cunew, FTYPE unew[][N2M][NPR], int timeorder, int numtimeorders, FTYPE *ndt);
  int advance_eno_u(int stage, FTYPE pi[][N2M][NPR], FTYPE ulast[][N2M][NPR], FTYPE pb[][N2M][NPR], FTYPE *CUf, FTYPE pf[][N2M][NPR], FTYPE *Cunew, FTYPE unew[][N2M][NPR], int timeorder, int numtimeorders, FTYPE *ndt);



#if(DOENO==0)
  MYFUN(advance_standard(stage,pi,ulast,pb,CUf,pf,Cunew,unew,timeorder,numtimeorders,ndt),"step_ch.c:advance()", "advance_standard()", 1);
  return(0);
#elif(DOENO==1)
  MYFUN(advance_eno_u(stage,pi,ulast,pb,CUf,pf,Cunew,unew,timeorder,numtimeorders,ndt),"step_ch.c:advance()", "advance_eno_u()", 1);
  return(0);
#elif(DOENO==2)
  MYFUN(advance_standard(stage,pi,ulast,pb,CUf,pf,Cunew,unew,timeorder,numtimeorders,ndt),"step_ch.c:advance()", "advance_standard()", 1);
  return(0);
#endif


}




int advance_standard(int stage,
	    FTYPE pi[][N2M][NPR],
	    FTYPE ulast[][N2M][NPR], 
	    FTYPE pb[][N2M][NPR],
	    FTYPE *CUf, FTYPE pf[][N2M][NPR],
	    FTYPE *Cunew, FTYPE unew[][N2M][NPR],
	    int timeorder, int numtimeorders,
	    FTYPE *ndt)
{
  int i, j, k, sc;
  FTYPE ndt1, ndt2;
  FTYPE Uf[NPR], Ui[NPR], Ub[NPR];
  FTYPE dUgeom[NPR],dUriemann[NPR],dUriemann1[NPR],dUriemann2[NPR],dUcomp[NUMSOURCES][NPR];
  struct of_geom geom;
  struct of_state q;
  FTYPE dUtot;
  FTYPE idx1,idx2;
  SFTYPE dt4diag;
  int finalstep;
  int flux_ct(int stage, FTYPE pb[][N2M][NPR],FTYPE F1[][N2M][NPR], FTYPE F2[][N2M][NPR]);
  void flux2dUavg(int i, int j, FTYPE F1[][N2M][NPR],FTYPE F2[][N2M][NPR],FTYPE *dUavg1,FTYPE *dUavg2);
  void dUtoU(FTYPE *dUgeom, FTYPE *dUriemann, FTYPE *CUf, FTYPE *Cunew, FTYPE *Ui,  FTYPE *ulast, FTYPE *Uf, FTYPE *unew);
  FTYPE (*dUriemannavg1)[N2M][NPR]; // for ENO point->avg approach
  FTYPE (*dUriemannavg2)[N2M][NPR]; // for ENO point->avg approach
  extern FTYPE a2cij( int dir, int i, int j, FTYPE ua[][N2M][NPR], FTYPE *uc );
  extern FTYPE c2aij( int dir, int i, int j, FTYPE uc[][N2M][NPR], FTYPE *ua );

  
  // tells diagnostics functions if should be accounting or not
  if(timeorder==numtimeorders-1){
    dt4diag=dt;
    finalstep=1;
  }
  else{
    dt4diag=-1.0;
    finalstep=0;
  }

#if(PRODUCTION==0)
  trifprintf( "#0f");
#endif
  // pb used here on a stencil, so if pb=pf or pb=pi in pointers, shouldn't change pi or pf yet -- don't currently
  MYFUN(fluxcalc(stage, pb, F1, 1, &ndt1),"step_ch.c:advance()", "fluxcalc", 1);
  MYFUN(fluxcalc(stage, pb, F2, 2, &ndt2),"step_ch.c:advance()", "fluxcalc", 2);
#if(FIXUPFLUX)
  fix_flux(pb,F1, F2);
#endif

  MYFUN(flux_ct(stage, pb,F1, F2),"step_ch.c:advance()", "flux_ct",1);


#if(PRODUCTION==0)
  trifprintf( "1f");
#endif
  // from here on, pi/pb/pf are only used a zone at a time rather than on a stencil

#if(DOENO==2)
  // allocate memory of dUriemannavg
  dUriemannavg1=uenotmp0;
  dUriemannavg2=uenotmp1;
  //   dUriemannavg=(FTYPE (*[N2M][NPR]))(&ua[0][0][0][0]);
  //   dUriemannavg=(FTYPE *[N2M][NPR])(&ua[2][2][0]);
  // get dU^{n+1}
  CZLOOP flux2dUavg(i,j,F1,F2,dUriemannavg1[i][j],dUriemannavg2[i][j]);
#endif


  /** now update pi to pf **/
  CZLOOP {

    /////////////
    //
    // Utoprim as initial conditions : can't assume want these to be same in the end, so assign
    //
    ////////////
    PLOOP pf[i][j][k] = pi[i][j][k];

    // initialize ulast and unew if very first time here
    if(timeorder==0) PLOOP ulast[i][j][k]=unew[i][j][k]=0.0;


    // set geometry for centered zone to be updated
    get_geometry(i, j, CENT, &geom);

    // find Ui(pi)
    MYFUN(get_state(pi[i][j], &geom, &q),"step_ch.c:advance()", "get_state()", 1);
    MYFUN(primtoU(UEVOLVE,pi[i][j], &q, &geom, Ui),"step_ch.c:advance()", "primtoU()", 1);


    // find dU(pb)
    // source() doesn't actually use CUf[2]=dt right now
    //    MYFUN(source(pb[i][j], &geom, i, j, dUcomp, dU,CUf[2]),"step_ch.c:advance()", "source", 1);
    MYFUN(source(pb[i][j], &geom, dUcomp, dUgeom),"step_ch.c:advance()", "source", 1);
    // assumes final dUcomp is nonzero and representative of source term over this timestep
    diag_source(&geom,dUcomp,dt4diag);

#if(DOENO==0)
    // dUriemann is actually average quantity, but we treat is as a point quantity at the zone center
    flux2dUavg(i,j,F1,F2,dUriemann1,dUriemann2);
    PLOOP dUriemann[k]=dUriemann1[k]+dUriemann2[k];
#elif(DOENO==2)
    //convert average dU to point dU at zone center
#if((N1BND>0)&&(N2BND>0))
    //now ua contains cell averaged conserved quantities, need to convert to cell-centered ones
    //convert to cell centered Uf
    a2cij(1,i,j,dUriemannavg1, dUriemann1 );
    a2cij(2,i,j,dUriemannavg2, dUriemann2 );
    PLOOP dUriemann[k]=dUriemann1[k]+dUriemann2[k];
#elif((N1BND>0)&&(N2BND==0))
    a2cij(1,i,j,dUriemannavg1, dUriemann );
#elif((N1BND==0)&&(N2BND>0))
    a2cij(2,i,j,dUriemannavg2, dUriemann );
#endif
#endif

    // find ulast==Uf and additional terms to unew
    dUtoU(dUgeom, dUriemann, CUf, Cunew, Ui, ulast[i][j], Uf, unew[i][j]);
    //    flux2U(i,j,F1,F2,dU,CUf,Cunew,Ui, ulast[i][j],Uf,unew[i][j]);

    // invert U->p
    if(timeorder==numtimeorders-1){ // last call, so unew is cooked and ready to eat!
      MYFUN(Utoprimgen(EVOLVEUTOPRIM,UEVOLVE,unew[i][j], &geom, pf[i][j]),"step_ch.c:advance()", "Utoprimgen", 1);
    }
    else{ // otherwise still iterating on primitives
      MYFUN(Utoprimgen(EVOLVEUTOPRIM,UEVOLVE,Uf, &geom, pf[i][j]),"step_ch.c:advance()", "Utoprimgen", 1);
    }



    // immediate local (i.e. 1-zone) fix
#if(FIXUPZONES==FIXUP1ZONE)
    if((STEPOVERNEGU==0)||(timeorder==numtimeorders-1)){
      MYFUN(fixup1zone(pf[i][j],&geom,finalstep),"fixup.c:fixup()", "fixup1zone()", 1);
    }
#endif
  }// end CZLOOP

#if(FIXUPZONES==FIXUPALLZONES)
  fixup(stage,pf,finalstep);
#endif  




  *ndt = defcon * 1. / (1. / ndt1 + 1. / ndt2);

#if(PRODUCTION==0)
  trifprintf( "2f");
#endif

  return (0);
}





// get dUavg
void flux2dUavg(int i, int j, FTYPE F1[][N2M][NPR],FTYPE F2[][N2M][NPR],FTYPE *dU1avg,FTYPE *dU2avg)
{
  FTYPE idx1,idx2;
  int k;

#if(VOLUMEDIFF==0)
  idx1=1.0/dx[RR];
  idx2=1.0/dx[TH];
#else
  idx1=idxvol[i][j][RR];
  idx2=idxvol[i][j][TH];
#endif


  if(FLUXB==FLUXCD){ // don't use volume reg. since differencing is large
    for(k=0;k<=U3;k++){
      dU1avg[k]=(
		 -(F1[i + 1][j][k] - F1[i][j][k]) *idx1
		 );
      dU2avg[k]=(
		 - (F2[i][j + 1][k] - F2[i][j][k]) *idx2
		 );
    }
    k=B1;
    dU1avg[k]=(
	       0.0
	       );
    dU2avg[k]=(
	       - (F2[i][j + 1][k] - F2[i][j - 1][k]) *idx2
	       );
    k=B2;
    dU1avg[k]=(
	       - (F1[i+1][j][k] - F1[i-1][j][k]) *idx1
	       );
    dU2avg[k]=(
	       0.0
	       );
    for(k=B3;k<NPR;k++){
      dU1avg[k]=(
		 -(F1[i + 1][j][k] - F1[i][j][k]) *idx1
		 );
      dU2avg[k]=(
		 - (F2[i][j + 1][k] - F2[i][j][k]) *idx2
		 );
    }
  }
  else{
    // other (normal) FLUXB methods
    PLOOP {
      dU1avg[k]=(
		 -(F1[i + 1][j][k] - F1[i][j][k]) *idx1
		 );
      dU2avg[k]=(
		 - (F2[i][j + 1][k] - F2[i][j][k]) *idx2
		 );
    }
  }




}






// convert point versions of U_i^{n} and dU -> U_i^{n+1} and other versions
void dUtoU(FTYPE *dUgeom, FTYPE *dUriemann, FTYPE *CUf, FTYPE *Cunew, FTYPE *Ui,  FTYPE *ulast, FTYPE *Uf, FTYPE *unew)
{
    int k;
    // finally assign new Uf and unew
    // store ulast to avoid recomputing U(pf) used later as pb for advance()
    PLOOP ulast[k]=Uf[k]   = CUf[0]*Ui[k] + CUf[1]*ulast[k] + CUf[2]*dt*(dUriemann[k]+dUgeom[k]);

    // how much of Ui, dU, and Uf to keep for final solution
    // ultimately unew is actual solution used to find final pf
    PLOOP unew[k] += Cunew[0]*Ui[k] + Cunew[1]*dt*(dUriemann[k]+dUgeom[k]) + Cunew[2]*Uf[k];



}











///////////////////////////////////////
//
// Sasha's "ENO" version of advance of conserved quantities
//
// ENO acts on U
//
///////////////////////////////////////

int advance_eno_u(int stage,
	    FTYPE pi[][N2M][NPR],
	    FTYPE ulast[][N2M][NPR], 
	    FTYPE pb[][N2M][NPR],
	    FTYPE *CUf, FTYPE pf[][N2M][NPR],
	    FTYPE *Cunew, FTYPE unew[][N2M][NPR],
	    int timeorder, int numtimeorders,
	    FTYPE *ndt)
{
  int i, j, k, sc;
  FTYPE ndt1, ndt2;
  FTYPE Uf[NPR], Ui[NPR], Ub[NPR];
  FTYPE dUgeom[NPR],dUriemann[NPR],dUriemann1[NPR],dUriemann2[NPR],dUcomp[NUMSOURCES][NPR];
  struct of_geom geom;
  struct of_state q;
  FTYPE dUtot;
  FTYPE idx1,idx2;
  SFTYPE dt4diag;
  int finalstep;
  int flux_ct(int stage, FTYPE pb[][N2M][NPR],FTYPE F1[][N2M][NPR], FTYPE F2[][N2M][NPR]);
  void flux2dUavg(int i, int j, FTYPE F1[][N2M][NPR],FTYPE F2[][N2M][NPR],FTYPE *dUavg1,FTYPE *dUavg2);
  void dUtoU(FTYPE *dUgeom, FTYPE *dUriemann, FTYPE *CUf, FTYPE *Cunew, FTYPE *Ui,  FTYPE *ulast, FTYPE *Uf, FTYPE *unew);
  extern FTYPE a2cij( int dir, int i, int j, FTYPE ua[][N2M][NPR], FTYPE *uc );
  extern FTYPE c2aij( int dir, int i, int j, FTYPE uc[][N2M][NPR], FTYPE *ua );
  FTYPE ucenter[NPR];
  FTYPE (*Uicenter)[N2M][NPR]; // for ENO point->avg approach
  FTYPE (*Uiavgtmp)[N2M][NPR]; // for ENO point->avg approach
  FTYPE (*Uiavg)[N2M][NPR]; // for ENO point->avg approach

  FTYPE (*Ufavg)[N2M][NPR];
  FTYPE (*Ufcentertmp)[N2M][NPR]; // for ENO point->avg approach
  FTYPE (*Ufcenter)[N2M][NPR]; // for ENO point->avg approach


  // tells diagnostics functions if should be accounting or not
  if(timeorder==numtimeorders-1){
    dt4diag=dt;
    finalstep=1;
  }
  else{
    dt4diag=-1.0;
    finalstep=0;
  }


  trifprintf( "#0f");
  // pb used here on a stencil, so if pb=pf or pb=pi in pointers, shouldn't change pi or pf yet -- don't currently
  MYFUN(fluxcalc(stage, pb, F1, 1, &ndt1),"step_ch.c:advance()", "fluxcalc", 1);
  MYFUN(fluxcalc(stage, pb, F2, 2, &ndt2),"step_ch.c:advance()", "fluxcalc", 2);
  fix_flux(pb, F1, F2);
  flux_ct(stage, pb, F1, F2);



  trifprintf( "1f");
  // from here on, pi/pb/pf are only used a zone at a time rather than on a stencil


  /** now update pi to pf **/


  // first get Ui
  Uicenter=uenotmp0;

  CZLOOP {

    /////////////
    //
    // Utoprim as initial conditions : can't assume want these to be same in the end, so assign
    //
    ////////////
    PLOOP pf[i][j][k] = pi[i][j][k];

    // set geometry for centered zone to be updated
    get_geometry(i, j, CENT, &geom);

    // find Ui(pi)
    MYFUN(get_state(pi[i][j], &geom, &q),"step_ch.c:advance()", "get_state()", 1);
    MYFUN(primtoU(UEVOLVE,pi[i][j], &q, &geom, Ui),"step_ch.c:advance()", "primtoU()", 1);
    
    // convert zone centered conserved quantities in Ui to zone averaged ones (atchekho)
    PLOOP Uicenter[i][j][k] = Ui[k];
  }

  // keeping Uicenter (uenotmp0) from above
  Uiavgtmp=uenotmp1;
  Uiavg=uenotmp2;

  // second, convert point Ui to average Ui
#if((N1BND>0)&&(N2BND>0))
  //now ua contains cell averaged conserved quantities, need to convert to cell-centered ones
  //convert to cell centered Uf
  CZLOOP c2aij( 1, i, j, Uicenter, Uiavgtmp[i][j] );
  CZLOOP c2aij( 2, i, j, Uiavgtmp, Uiavg[i][j] );
#elif((N1BND>0)&&(N2BND==0))
  CZLOOP c2aij( 1, i, j, Uicenter, Uiavg[i][j] );
#elif((N1BND==0)&&(N2BND>0))
  CZLOOP c2aij( 2, i, j, Uicenter, Uiavg[i][j] );
#endif
  //  CZLOOP c2aij( 1, i, j, uc, ua[i][j] );
  //  c2a( ua, uc );
  //now ua contains cell averaged values

  // only Uiavg is kept here (using uenotmp2)
  Ufavg=uenotmp0;

  // get average Uf
  CZLOOP {
    // set geometry for centered zone to be updated
    get_geometry(i, j, CENT, &geom); //added by atch

    // find dU(pb) //added by atch
    MYFUN(source(pb[i][j], &geom, dUcomp, dUgeom),"step_ch.c:advance()", "source", 1);
    diag_source(&geom,dUcomp,dt4diag);

    // Ui average from ua
    PLOOP Ui[k] = Uiavg[i][j][k];  //replace with cell averaged values
  
    // initialize ulast and unew if very first time here
    if(timeorder==0) PLOOP ulast[i][j][k]=unew[i][j][k]=0.0;

    // find ulast==Uf and additional terms to unew
    flux2dUavg(i,j,F1,F2,dUriemann1,dUriemann2);
    PLOOP dUriemann[k]=dUriemann1[k]+dUriemann2[k];
    dUtoU(dUgeom, dUriemann, CUf, Cunew, Ui, ulast[i][j], Uf, unew[i][j]);
    //    flux2U(i,j,F1,F2,dU,CUf,Cunew,Ui, ulast[i][j],Uf,unew[i][j]);

    // invert U->p (but do conversion from averages to cell-centered quantities first!)
    if(timeorder==numtimeorders-1){ // last call, so unew is cooked and ready to eat!
      // instead of unew  need to supply cell centered quantities
      PLOOP Ufavg[i][j][k] = unew[i][j][k];
      //MYFUN(Utoprimgen(EVOLVEUTOPRIM,UEVOLVE,unew[i][j], &geom, pf[i][j]),"step_ch.c:advance()", "Utoprimgen", 1);
    }
    else{ // otherwise still iterating on primitives
      // instead of Uf  need to supply cell centered quantities
      PLOOP Ufavg[i][j][k] = Uf[k];
      //MYFUN(Utoprimgen(EVOLVEUTOPRIM,UEVOLVE,Uf, &geom, pf[i][j]),"step_ch.c:advance()", "Utoprimgen", 1);
    }
  }

  // only Ufavg is kept here (using uenotmp0)
  Ufcenter=uenotmp1;
  Ufcentertmp=uenotmp2;

#if((N1BND>0)&&(N2BND>0))
  //now ua contains cell averaged conserved quantities, need to convert to cell-centered ones
  //convert to cell centered Uf
  CZLOOP a2cij( 1, i, j, Ufavg, Ufcentertmp[i][j] );
  CZLOOP a2cij( 2, i, j, Ufcentertmp, Ufcenter[i][j] );
#elif((N1BND>0)&&(N2BND==0))
  CZLOOP a2cij( 1, i, j, Ufavg, Ufcenter[i][j] );
#elif((N1BND==0)&&(N2BND>0))
  CZLOOP a2cij( 2, i, j, Ufavg, Ufcenter[i][j] );
#endif
  
  // now take Uf average (ua) and get Uf at point and invert to point primitive
  CZLOOP {
    // set geometry for centered zone to be updated
    get_geometry(i, j, CENT, &geom); //added by atch

    // invert U->p
    MYFUN(Utoprimgen(EVOLVEUTOPRIM,UEVOLVE,Ufcenter[i][j], &geom, pf[i][j]),"step_ch.c:advance()", "Utoprimgen", 1);
    
    // immediate local (i.e. 1-zone) fix
#if(FIXUPZONES==FIXUP1ZONE)
    if((STEPOVERNEGU==0)||(timeorder==numtimeorders-1)){
      MYFUN(fixup1zone(pf[i][j],&geom,finalstep),"fixup.c:fixup()", "fixup1zone()", 1);
    }
#endif
  }// end CZLOOP

#if(FIXUPZONES==FIXUPALLZONES)
  fixup(stage,pf,finalstep);
#endif  


  *ndt = defcon * 1. / (1. / ndt1 + 1. / ndt2);

  trifprintf( "2f");

  return (0);
}









