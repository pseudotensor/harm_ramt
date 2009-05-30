
#include "decs.h"

/* apply floors to density, internal energy */

// currently called before bound, which assumes bound sets boundary
// values exactly as wanted without any fixing.

#define JONFIXUP 1 // 0=gammie 1=jon's

#if 0
void fixup(int stage,FTYPE (*pv)[N2M][NPR],int finalstep)
{
  int i, j;



  CZLOOP{ pfixup(pv[i][j], i, j);}

}
#endif

// operations that require synch of boundary zones in MPI, or that require use of boundary zones at all
// operations that only need to be done inside computational loop
int pre_fixup(int stage,FTYPE (*pv)[N2M][NPR])
{

  // grab b^2 flags (above fixup may change u or rho, so must do this after)
  get_bsqflags(stage,pv);


  return(0);
}

// operations that require synch of boundary zones in MPI, or that require use of boundary zones at all
// this function actually changes primitives
int post_fixup(int stageit,FTYPE (*pv)[N2M][NPR],int finalstep)
{
  int stage,stagei,stagef;
  int boundstage;

#if(UTOPRIMADJUST!=0)
  ////////////////////////////////////
  //
  // utoprim fixup of primitive solution



  if(SIMULBCCALC<=0){ stagei=STAGEM1; stagef=STAGEM1; }
  else if(SIMULBCCALC==1) { stagei=STAGE0; stagef=STAGE2;}
  else if(SIMULBCCALC==2) { stagei=STAGE0; stagef=STAGE5;}

  if(SIMULBCCALC) boundstage=STAGE0;
  else boundstage=STAGEM1;
  for(stage=stagei;stage<=stagef;stage++){

    // first bound failure flag
    // could optimize bound of pflag since often failures don't occur (just ask if any failures first), although probably negligible performance hit
    if(stage<STAGE2){
      bound_pflag(boundstage, pflag);
      if(stage!=STAGEM1) boundstage++;
    }


    // check for bad solutions
#if(CHECKSOLUTION)
    fixup_checksolution(stage,pv,finalstep);
#endif

#if(EOMTYPE==EOMGRMHD)
    fixup_utoprim(stage,pv,finalstep); // fixup before new solution (has to be here since need previous stage's failure flag)
#endif

    if(stage<STAGE2){
      if(stage!=STAGEM1){
	bound_pflag(boundstage, pflag);
	boundstage++;
      }
    }

  }
#endif

  ////////////////////////////////////
  //
  // standard fixup of floor
  // in this case stageit=-1, so does all stages
  //  fixup(stageit,pv,finalstep);




  return(0);
}

#if(JONFIXUP==1)

int fixup(int stage,FTYPE (*pv)[N2M][NPR], int finalstep)
{
  int i, j;
  struct of_geom geom;


  CZLOOP{
    //    fprintf(fail_file,"i=%d j=%d\n",i,j); fflush(fail_file);
    get_geometry(i,j,CENT,&geom);
    if(fixup1zone(pv[i][j],&geom,finalstep)>=1)
      FAILSTATEMENT("fixup.c:fixup()", "fixup1zone()", 1);
  }
  return(0);
}
#else

// GAMMIE OLD FIXUP (not kept up to date)
int fixup(int stage,FTYPE (*pv)[N2M][NPR],int finalstep)
{
  int i, j, k;
  int ip, jp, im, jm;
  FTYPE bsq, del;
  FTYPE ftempA,ftempB;
  struct of_geom geom ;
  FTYPE prfloor[NPR];

  if(EOMTYPE==EOMFFDE) return(0); // nothing to do


  CZLOOP {
    get_geometry(i,j,CENT,&geom) ;
    // densities
    set_density_floors(&geom,pv[i][j],prfloor);
    
    /* floor on density (momentum *not* conserved) */
    if (pv[i][j][RHO] < prfloor[RHO]) {
#if(FLOORDIAGS)
      fladd[RHO] +=
	  dVF * geom->g * (prfloor[RHO] - pv[i][j][RHO]);
#endif
      pv[i][j][RHO] = prfloor[RHO];
    }

    /* floor on internal energy */
    if (pv[i][j][UU] < prfloor[UU]) {
#if(FLOORDIAGS)
      fladd[UU] +=
	  dVF * geom->g * (prfloor[UU] - pv[i][j][UU]);
#endif
      pv[i][j][UU] = prfloor[UU];
    }

    /* limit gamma wrt normal observer */
#if(WHICHVEL==VELREL4)
      if(limit_gamma(GAMMAMAX,pv[i][j],&geom,-1.0)>=1)  // need general accounting for entire routine.
	FAILSTATEMENT("fixup.c:fixup()", "limit_gamma()", 1);
#endif
    /* 
       if(bsq_calc(pv[i][j],&geom,&bsq)>=1){
       fprintf(fail_file,"1init:bsq_calc: failure\n");
       fflush(fail_file); return(1); } if(bsq >
       MAXBSQOVERRHO*pv[i][j][RHO]){
       fladd[RHO]+=dVF*gdet[i][j][CENT]*(prfloor[RHO]-pv[i][j][RHO]);
       pv[i][j][RHO] = bsq/MAXBSQOVERRHO ; } if(bsq >
       MAXBSQOVERUU*pv[i][j][UU]){
       fladd[UU]+=dVF*gdet[i][j][CENT]*(prfloor[UU]-pv[i][j][UU]);
       pv[i][j][UU] = bsq/MAXBSQOVERUU ; } */
  }

  return(0);}


#endif

// account for changes by tracking conserved quantities
// accounts for both failures and floor recoveries
int diag_fixup(FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, int finalstep,int whocalled)
{
  struct of_state q;
  FTYPE Ui[NPR],Uf[NPR];
  FTYPE ftemp[NPR];
  int failreturn;
  int k,enerregion, tscale;

  //  fprintf(stderr,"before diag_fixup %d %d %d\n",ptrgeom->i,ptrgeom->j,whocalled);

#if(DOSUPERDEBUG)
      superdebug(pr0,pr,ptrgeom,whocalled);
      // collect values for non-failed and failed zones
#endif


  // count every time corrects, not just on conserved quantity tracking time
  if(DODEBUG){
    TSCALELOOP failfloorcount[ptrgeom->i][ptrgeom->j][tscale][whocalled]++;
  }

  if(finalstep){ // only account if on full timestep
    ENERREGIONLOOP{ // could be designed more efficiently, but not called too often
      enerpos=enerposreg[enerregion];
      fladd=fladdreg[enerregion];
      fladdterms=fladdtermsreg[enerregion];

      // only account if within active zones for that region
      if(((ptrgeom->i)>=enerpos[X1DN])&&((ptrgeom->i)<=enerpos[X1UP])&&((ptrgeom->j)>=enerpos[X2DN])&&((ptrgeom->j)<=enerpos[X2UP])){
	// before any changes
	failreturn=get_state(pr0,ptrgeom,&q);
	if(failreturn>=1) dualfprintf(fail_file,"get_state(1) failed in fixup.c, why???\n");
	failreturn=primtoU(UDIAG,pr0,&q,ptrgeom,Ui);
	if(failreturn>=1) dualfprintf(fail_file,"primtoU(1) failed in fixup.c, why???\n");
	
	// after any changes
	failreturn=get_state(pr,ptrgeom,&q);
	if(failreturn>=1) dualfprintf(fail_file,"get_state(2) failed in fixup.c, why???\n");
	failreturn=primtoU(UDIAG,pr,&q,ptrgeom,Uf);
	if(failreturn>=1) dualfprintf(fail_file,"primtoU(2) failed in fixup.c, why???\n");

	PLOOP{
	  ftemp[k]=dVF * (Uf[k]-Ui[k]);
	  fladdterms[whocalled][k] += (SFTYPE)ftemp[k];
	  fladd[k] += ftemp[k];
	}
      }
    }
  }



  //  fprintf(stderr,"after diag_fixup\n"); fflush(stderr);

  return(0);
}


// account for changes by tracking conserved quantities
// accounts for both failures and floor recoveries
int diag_fixup_U(FTYPE *Ui, FTYPE *Uf, struct of_geom *ptrgeom, int finalstep,int whocalled)
{
  struct of_state q;
  FTYPE ftemp[NPR];
  int failreturn;
  int k,enerregion, tscale;


  // count every time corrects, not just on conserved quantity tracking time
  if(DODEBUG){
    TSCALELOOP failfloorcount[ptrgeom->i][ptrgeom->j][tscale][whocalled]++;
  }

  if(1||finalstep){ // only account if on full timestep
    ENERREGIONLOOP{ // could be designed more efficiently, but not called too often
      enerpos=enerposreg[enerregion];
      fladd=fladdreg[enerregion];
      fladdterms=fladdtermsreg[enerregion];

      // only account if within active zones for that region
      if(((ptrgeom->i)>=enerpos[X1DN])&&((ptrgeom->i)<=enerpos[X1UP])&&((ptrgeom->j)>=enerpos[X2DN])&&((ptrgeom->j)<=enerpos[X2UP])){

	PLOOP{
	  ftemp[k]=dVF * (Uf[k]-Ui[k]);
	  fladdterms[whocalled][k] += (SFTYPE)ftemp[k];
	  fladd[k] += ftemp[k];
	}
      }
    }
  }



  //  fprintf(stderr,"after diag_fixup\n"); fflush(stderr);

  return(0);
}


// finalstep==0 is non-accounting, finalstep==1 is accounting
int fixup1zone(FTYPE *pr, struct of_geom *ptrgeom, int finalstep)
{
  int k;
  int ip, jp, im, jm;
  FTYPE bsq, del;
  FTYPE r, th, X[NDIM];
  FTYPE ftempA,ftempB;
  struct of_state q;
  FTYPE prfloor[NPR];
  FTYPE pr0[NPR];
  int checkfl[NPR];
  int failreturn;
  int didchangeprim;
  FTYPE scalemin[NPR];


  //  if(EOMTYPE==EOMFFDE) return(0); // nothing to do


  if(EOMTYPE==EOMGRMHD){
    // assign general floor variables
    // whether to check floor condition
    PLOOP{ checkfl[k]=0;}
    checkfl[RHO]=1;
    checkfl[UU]=1;

    // densities
    set_density_floors(ptrgeom,pr,prfloor);

    scalemin[RHO]=RHOMINLIMIT;
    scalemin[UU]=UUMINLIMIT;


    PLOOP{
      if(checkfl[k]){
	if(prfloor[k]<scalemin[k]) prfloor[k]=scalemin[k];
      }
      pr0[k]=pr[k];
    }

    /////////////////////////////
    //
    // perform changes

    // shouldn't fail since before and after states should be ok, as
    // long as would have changed the value.  Check will occur if
    // simulation continues ok.  Could place check inside if below.
    didchangeprim=0;


    PLOOP{
      if ( checkfl[k]&&(prfloor[k] > pr[k]) ){
	didchangeprim=1;
	//dualfprintf(fail_file,"%d : %d %d : %d : %d : %21.15g - %21.15g\n",k,ptrgeom->i,ptrgeom->j,ptrgeom->p,checkfl[k],prfloor[k],pr[k]); 
	// only add on full step since middle step is not really updating primitive variables
	pr[k]=prfloor[k];
      }
    }

    // since inflow check is on boundary values, no need for inflow check here
  }


#if(WHICHVEL==VELREL4)
  failreturn=limit_gamma(GAMMAMAX,pr,ptrgeom,-1.0);
  if(failreturn>=1) FAILSTATEMENT("fixup.c:fixup()", "limit_gamma()", 1);
  if(failreturn==-1) didchangeprim=1;
#endif// end if WHICHVEL==VEL4REL



  /* limit gamma wrt normal observer */
#if(WHICHVEL==VELREL4)
#endif

  /* 
     if(bsq_calc(pr,ptrgeom,&bsq)>=1){
     fprintf(fail_file,"1init:bsq_calc: failure\n");
     fflush(fail_file); return(1); } if(bsq >
     MAXBSQOVERRHO*pr[RHO]){
     fladd[RHO]+=dVF*gdet[i][j][gpos]*(prfloor[RHO]-pr[RHO]);
     pr[RHO] = bsq/MAXBSQOVERRHO ; } if(bsq >
     MAXBSQOVERUU*pr[UU]){
     fladd[UU]+=dVF*gdet[i][j][gpos]*(prfloor[UU]-pr[UU]);
     pr[UU] = bsq/MAXBSQOVERUU ; }
  */

  //////////////////////////
  // now keep track of modified primitives via conserved quantities

  if(didchangeprim){
    // assume once we go below floor, all hell will break loose unless we calm the storm by shutting down this zone's relative velocity
    // normal observer velocity
    // i.e. consider this a failure
    //pflag[ptrgeom->i][ptrgeom->j][FLAGUTOPRIMFAIL]= 1;
  }

  
  if(didchangeprim&&FLOORDIAGS){// FLOORDIAGS includes fail diags
    diag_fixup(pr0, pr, ptrgeom, finalstep,COUNTFLOORACT);
  }

  return(0);
}


// number of vote checks per zone
#define MAXVOTES 8

// check whether solution seems reasonable
// useful if b^2/rho\gg 1 or try to approach stationary model where variations in space shouldn't be large zone to zone
// checks whether u or gamma is much different than surrounding zones.
// if true, then flag as failure, else reasonable solution and keep
// can't assume failed zones are reasonably set
int fixup_checksolution(int stage, FTYPE (*pv)[N2M][NPR],int finalstep)
{
  int i,j,k,l;
  struct of_geom geom;
  FTYPE gamma,ucon[NDIM];
  static FTYPE gammacheck[N1+2][N2+2];
  FTYPE percdiff[8],percdiff2[8];
  int pp,qq;
  int vote[2];
  int numvotes;

  //  if(EOMTYPE==EOMFFDE) return(0); // nothing to do

  // determine gamma to be used
  DQZLOOP {
    if((pflag[i][j][FLAGBSQORHO]||pflag[i][j][FLAGBSQOU])&&(pflag[i][j][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)){
      get_geometry(i,j,CENT,&geom);
      
#if(WHICHVEL==VELREL4)
      MYFUN(gamma_calc(pv[i][j],&geom,&gammacheck[i+1][j+1]),"fixup_checksolution: gamma calc failed\n","fixup.c",1);
#else
      if (ucon_calc(pv[i][j], &geom, ucon) >= 1)  FAILSTATEMENT("fixup.c:fixup_checksolution()", "ucon_calc()", 1);
      gammacheck[i+1][j+1]=ucon[TT];
#endif
    }
  }

  // see if percent differences (actually factors) are too large
  CZLOOP {
    // doesn't need bound of pflag since don't check pflag of surrounding values (assumes if comparing with failure point, failure point is reasonable afterwards if before)
    if((pflag[i][j][FLAGBSQORHO]||pflag[i][j][FLAGBSQOU])&&(pflag[i][j][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)){// if b^2/{rho,u}\gg 1 and not failure already, check if solution is reasonable
      pp=i+1;
      qq=j+1;
      /*
      percdiff[0]=2.0*(gammacheck[pp][qq+1]-gammacheck[pp][qq])/(gammacheck[pp][qq+1]+gammacheck[pp][qq]);
      percdiff[1]=2.0*(gammacheck[pp][qq-1]-gammacheck[pp][qq])/(gammacheck[pp][qq-1]+gammacheck[pp][qq]);
      percdiff[2]=2.0*(gammacheck[pp+1][qq]-gammacheck[pp][qq])/(gammacheck[pp+1][qq]+gammacheck[pp][qq]);
      percdiff[3]=2.0*(gammacheck[pp-1][qq]-gammacheck[pp][qq])/(gammacheck[pp-1][qq]+gammacheck[pp][qq]);

      percdiff[4]=2.0*(gammacheck[pp+1][qq+1]-gammacheck[pp][qq])/(gammacheck[pp+1][qq+1]+gammacheck[pp][qq]);
      percdiff[5]=2.0*(gammacheck[pp+1][qq-1]-gammacheck[pp][qq])/(gammacheck[pp+1][qq-1]+gammacheck[pp][qq]);
      percdiff[6]=2.0*(gammacheck[pp+1][qq+1]-gammacheck[pp][qq])/(gammacheck[pp+1][qq+1]+gammacheck[pp][qq]);
      percdiff[7]=2.0*(gammacheck[pp-1][qq-1]-gammacheck[pp][qq])/(gammacheck[pp-1][qq-1]+gammacheck[pp][qq]);

      percdiff2[0]=2.0*(pv[i][j+1][UU]-pv[i][j][UU])/(pv[i][j+1][UU]+pv[i][j][UU]);
      percdiff2[1]=2.0*(pv[i][j-1][UU]-pv[i][j][UU])/(pv[i][j-1][UU]+pv[i][j][UU]);
      percdiff2[2]=2.0*(pv[i+1][j][UU]-pv[i][j][UU])/(pv[i+1][j][UU]+pv[i][j][UU]);
      percdiff2[3]=2.0*(pv[i-1][j][UU]-pv[i][j][UU])/(pv[i-1][j][UU]+pv[i][j][UU]);

      percdiff2[4]=2.0*(pv[i+1][j+1][UU]-pv[i][j][UU])/(pv[i+1][j+1][UU]+pv[i][j][UU]);
      percdiff2[5]=2.0*(pv[i+1][j-1][UU]-pv[i][j][UU])/(pv[i+1][j-1][UU]+pv[i][j][UU]);
      percdiff2[6]=2.0*(pv[i+1][j+1][UU]-pv[i][j][UU])/(pv[i+1][j+1][UU]+pv[i][j][UU]);
      percdiff2[7]=2.0*(pv[i-1][j-1][UU]-pv[i][j][UU])/(pv[i-1][j-1][UU]+pv[i][j][UU]);
      */
      percdiff[0]=(pflag[i][j+1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL) ? (gammacheck[pp][qq+1]/gammacheck[pp][qq]) : -1;
      percdiff[1]=(pflag[i][j-1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL) ? (gammacheck[pp][qq-1]/gammacheck[pp][qq]) : -1;
      percdiff[2]=(pflag[i+1][j][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL) ? (gammacheck[pp+1][qq]/gammacheck[pp][qq]) : -1;
      percdiff[3]=(pflag[i-1][j][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL) ? (gammacheck[pp-1][qq]/gammacheck[pp][qq]) : -1;

      percdiff[4]=(pflag[i+1][j+1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL) ? (gammacheck[pp+1][qq+1]/gammacheck[pp][qq]) : -1;
      percdiff[5]=(pflag[i+1][j-1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL) ? (gammacheck[pp+1][qq-1]/gammacheck[pp][qq]) : -1;
      percdiff[6]=(pflag[i+1][j+1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL) ? (gammacheck[pp+1][qq+1]/gammacheck[pp][qq]) : -1;
      percdiff[7]=(pflag[i-1][j-1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL) ? (gammacheck[pp-1][qq-1]/gammacheck[pp][qq]) : -1;

      if(EOMTYPE==EOMGRMHD){
	percdiff2[0]=(pflag[i][j+1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL) ? (pv[i][j+1][UU]/pv[i][j][UU]) : -1;
	percdiff2[1]=(pflag[i][j-1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL) ? (pv[i][j-1][UU]/pv[i][j][UU]) : -1;
	percdiff2[2]=(pflag[i+1][j][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL) ? (pv[i+1][j][UU]/pv[i][j][UU]) : -1;
	percdiff2[3]=(pflag[i-1][j][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL) ? (pv[i-1][j][UU]/pv[i][j][UU]) : -1;

	percdiff2[4]=(pflag[i+1][j+1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL) ? (pv[i+1][j+1][UU]/pv[i][j][UU]) : -1;
	percdiff2[5]=(pflag[i+1][j-1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL) ? (pv[i+1][j-1][UU]/pv[i][j][UU]) : -1;
	percdiff2[6]=(pflag[i+1][j+1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL) ? (pv[i+1][j+1][UU]/pv[i][j][UU]) : -1;
	percdiff2[7]=(pflag[i-1][j-1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL) ? (pv[i-1][j-1][UU]/pv[i][j][UU]) : -1;
      }

      vote[0]=vote[1]=0;
      numvotes=0;
      for(l=0;l<MAXVOTES;l++){
	if(percdiff[0]>=0.0){//assumes failure same for any type of check.  No vote for failed zones
	  if( (fabs(percdiff[l])>GAMMAPERCDIFFMAX)||(fabs(percdiff[l])<1.0/GAMMAPERCDIFFMAX) ) vote[0]++;
	  if((EOMTYPE==EOMGRMHD)&&( (fabs(percdiff2[l])>UPERCDIFFMAX)||(fabs(percdiff2[l])<1.0/UPERCDIFFMAX) ) ) vote[1]++;
	  numvotes++;
	}
      }
      // use majority rule.  This allows shock fronts, but no faster members.
      if(vote[0]>numvotes*0.5){// then majority rules
	//	fprintf(stderr,"caught one-0: %d %d\n",i,j);
	pflag[i][j][FLAGUTOPRIMFAIL]=UTOPRIMFAILGAMMAPERC;
      }
      if((EOMTYPE==EOMGRMHD)&&(vote[1]>numvotes*0.5)){// then majority rules
	// overrides above if both true
	//	fprintf(stderr,"caught one-1: %d %d\n",i,j);
	pflag[i][j][FLAGUTOPRIMFAIL]=UTOPRIMFAILUPERC;
      }

    } // end if b^2/rho too large
  }// end CZLOOP

  return(0);
}






#define AVG4_1(pr,i,j,k) (0.25*(pr[i][j+1][k]+pr[i][j-1][k]+pr[i-1][j][k]+pr[i+1][j][k]))
#define AVG4_2(pr,i,j,k) (0.25*(pr[i+1][j+1][k]+pr[i+1][j-1][k]+pr[i-1][j+1][k]+pr[i-1][j-1][k]))

#define AVG2_1(pr,i,j,k) (0.5*(pr[i][j+1][k]+pr[i][j-1][k]))
#define AVG2_2(pr,i,j,k) (0.5*(pr[i+1][j][k]+pr[i-1][j][k]))
#define AVG2_3(pr,i,j,k) (0.5*(pr[i+1][j+1][k]+pr[i-1][j-1][k]))
#define AVG2_4(pr,i,j,k) (0.5*(pr[i+1][j-1][k]+pr[i-1][j+1][k]))

#if(MPIEQUALNONMPI==1)
#define ORDERINDEPENDENT 1 // no choice
// whether fixup_utoprim should be order-independent or not
#else
#define ORDERINDEPENDENT 1 // choice
#endif

#define NOGOODSTATIC 0
#define NOGOODRESET 1
#define NOGOODAVERAGE 2

// no good values.  Choices: STATIC, RESET to something, and AVERAGE failed (t-dt surrounding) values.
#define TODONOGOOD 2
// 0=STATIC // model-independent, but really bad idea -- causes death of computation for accretion disks in jet region.
// 1=RESET // model-dependent
// 2=AVERAGE // model-independent -- probably best idea -- and then hope for the best.

#define HANDLEUNEG 0
// seems to keep failing with this, so probably treating u<0 like failure is better idea
// 0: treat as full failure
// 1: treat as floor issue

// fix the bad solution as determined by utoprim() and fixup_checksolution()
// needs fail flag over -1..N, but uses p at 0..N-1
int fixup_utoprim(int stage, FTYPE (*pv)[N2M][NPR],int finalstep)
{
  int i,j,k;
  struct of_geom geom;
  FTYPE gamma,alpha,vsq,ucon[NDIM];
  FTYPE pr0[NPR],Uf[NPR],Ui[NPR];
  struct of_state q;
  int failreturn;
  int mypflag,utoprimfailtype;
  FTYPE (*ptoavg)[N2M][NPR];
  FTYPE prguess[NPR];
  int fixed;

  // this average only works if using 4 velocity since only then guaranteed solution is good after interpolation
  if(WHICHVEL==VEL3) return(0); // just stick with static, best can do
  if(EOMTYPE==EOMFFDE) return(0); // nothing to do

  if(ORDERINDEPENDENT){
    ptoavg=prc;
    // prc is used since in step_ch.c not needed after each advance()
    // prc is template for values used to fixup failures.
    // In fixup-loop below, only change/fixup pv, not prc.  Otherwise order-dependent and not MPI friendly
    // prc is only used for surrounding zones used to fixup local zone
    FULLLOOP PLOOP prc[i][j][k]=pv[i][j][k];
  }
  else ptoavg=pv;

  // first check for bad solutions and try to fix based upon average surrounding solution
  if(UTOPRIMADJUST==UTOPRIMAVG){
    CZLOOP {
      mypflag=pflag[i][j][FLAGUTOPRIMFAIL];
      /////////////////
      //
      // see if utoprim() failed
      //
      //////////////////
      
      // set pre-primitive (needed outside failure "if" when superdebugging
#if(DOSUPERDEBUG==1)
      PLOOP    pr0[k]=pv[i][j][k];
      get_geometry(i,j,CENT,&geom);
#endif

      if(mypflag!=0){
	fixed=0; // assume not fixed yet

	// set pre-fixed primitives
	// put back inside "if" when not superdebugging since wasteful of cpu
#if(DOSUPERDEBUG==0)
	PLOOP    pr0[k]=pv[i][j][k];
	get_geometry(i,j,CENT,&geom);
#endif

	////////////////////////////
	//
	// deal with negative internal energy as special case of failure
	//
	if(!fixed){
	  if(mypflag==UTOPRIMFAILUNEG){

	    if((STEPOVERNEGU==0)||(STEPOVERNEGU&&finalstep)){

	      if(HANDLEUNEG==1){
		// set back to floor level
		get_geometry(i,j,CENT,&geom);
		set_density_floors(&geom,pv[i][j],prguess);
		// GODMARK -- maybe too agressive, maybe allow more negative?
		
		if(UTOPRIMFAILRETURNTYPE==UTOPRIMRETURNADJUSTED){
		  // then pv is previous timestep value and can use to make fix
		  if(-pv[i][j][UU]<prguess[UU]){ fixed=1; pv[i][j][UU]=prguess[UU];} // otherwise assume really so bad that failure
		  else fixed=0;
		}
		else{
		  // just treat as floor for all failures since do not know what updated quantity is
		  pv[i][j][UU]=prguess[UU];
		  fixed=1;
		}
	      }// end if handling uneg in special way
	    }// end if not allowing negative u or if allowing but not yet final step
	    else if((STEPOVERNEGU)&&(!finalstep)){
	      fixed=1; // tells rest of routine to leave alone and say ok solution, but don't use it to fix convergence failures for other zones
	    }
	  }// end if u<0
	}// end if not fixed


	//////////////////
	// other kinds of failures
	//
	if(fixed==0){
#if(WHICHVEL==VELREL4)
	  //////////////
	  //
	  // check gamma to so calibrate new gamma to no larger than previous gamma
	  //
	  /////////////
	  MYFUN(gamma_calc(pv[i][j],&geom,&gamma),"fixup_utoprim: gamma calc failed\n","fixup.c",1);
	  if (ucon_calc(pv[i][j], &geom, ucon) >= 1)
	    FAILSTATEMENT("utoprim_ldz.c:Utoprim_ldz()", "ucon_calc()", 1);
	  alpha = 1. / sqrt(-geom.gcon[0][0]);
	  vsq = 1. - 1. / (alpha * alpha * ucon[0] * ucon[0]);
	  if(debugfail>=2) dualfprintf(fail_file,"initial gamma: %21.15g,  max: %21.15g, initial vsq: %21.15g\n",gamma,GAMMAMAX,vsq);
	  if(gamma>GAMMAMAX) gamma=GAMMAMAX;
#endif


	  /////////////////////
	  //
	  // fix using average of surrounding good values, if they exist
	  //

	  /////////////
	  //
	  // 4 VALUES
	  //
	  /////////////
	  if( // but if surrounded by good values
	     (pflag[i][j+1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)&&
	     (pflag[i][j-1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)&&
	     (pflag[i+1][j][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)&&
	     (pflag[i-1][j][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)
	     ){
	    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d : utoprim corrected1\n",t,startpos[1]+i,startpos[2]+j);
	    // then average
	    // don't mess with B field, it's correct already
	    for(k=0;k<=U3;k++){
	      pv[i][j][k]=AVG4_1(ptoavg,i,j,k);
	      if(debugfail>=2) dualfprintf(fail_file,"uc2: k=%d pv=%21.15g\n",k,pv[i][j][k]);
	    }
	  }
	  else if( // but if surrounded by good values
		  (pflag[i+1][j+1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)&&
		  (pflag[i+1][j-1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)&&
		  (pflag[i-1][j+1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)&&
		  (pflag[i-1][j-1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)
		  ){
	    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d : utoprim corrected2\n",t,startpos[1]+i,startpos[2]+j);
	    // then average
	    for(k=0;k<=U3;k++){
	      pv[i][j][k]=AVG4_2(ptoavg,i,j,k);
	      if(debugfail>=2) dualfprintf(fail_file,"uc2: k=%d pv=%21.15g\n",k,pv[i][j][k]);
	    }
	  }
	  /////////////
	  //
	  // 2 VALUES
	  //
	  /////////////
	  else if( // but if "surrounded" by good values
		  (pflag[i][j+1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)&&
		  (pflag[i][j-1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)
		  ){
	    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d : utoprim corrected3\n",t,startpos[1]+i,startpos[2]+j);
	    // then average
	    for(k=0;k<=U3;k++){
	      pv[i][j][k]=AVG2_1(ptoavg,i,j,k);
	      if(debugfail>=2) dualfprintf(fail_file,"uc2: k=%d pv=%21.15g\n",k,pv[i][j][k]);
	    }
	  }
	  else if( // but if "surrounded" by good values
		  (pflag[i+1][j][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)&&
		  (pflag[i-1][j][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)
		  ){
	    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d : utoprim corrected4\n",t,startpos[1]+i,startpos[2]+j);
	    // then average
	    for(k=0;k<=U3;k++){
	      pv[i][j][k]=AVG2_2(ptoavg,i,j,k);
	      if(debugfail>=2) dualfprintf(fail_file,"uc2: k=%d pv=%21.15g\n",k,pv[i][j][k]);
	    }
	  }
	  else if( // but if "surrounded" by good values
		  (pflag[i+1][j+1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)&&
		  (pflag[i-1][j-1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)
		  ){
	    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d : utoprim corrected5\n",t,startpos[1]+i,startpos[2]+j);
	    // then average
	    for(k=0;k<=U3;k++){
	      pv[i][j][k]=AVG2_3(ptoavg,i,j,k);
	      if(debugfail>=2) dualfprintf(fail_file,"uc2: k=%d pv=%21.15g\n",k,pv[i][j][k]);
	    }
	  }
	  else if( // but if "surrounded" by good values
		  (pflag[i+1][j-1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)&&
		  (pflag[i-1][j+1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)
		  ){
	    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d : utoprim corrected6\n",t,startpos[1]+i,startpos[2]+j);
	    // then average
	    for(k=0;k<=U3;k++){
	      pv[i][j][k]=AVG2_4(ptoavg,i,j,k);
	      if(debugfail>=2) dualfprintf(fail_file,"uc2: k=%d pv=%21.15g\n",k,pv[i][j][k]);
	    }
	  }
	  /////////////
	  //
	  // SINGLE VALUES
	  //
	  /////////////
	  else if( // but if "surrounded" by good value
		  (pflag[i+1][j+1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)
		  ){
	    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j);
	    // then ASSIGN
	    for(k=0;k<=U3;k++){
	      pv[i][j][k]=ptoavg[i+1][j+1][k];
	      if(debugfail>=2) dualfprintf(fail_file,"uc2: k=%d pv=%21.15g\n",k,pv[i][j][k]);
	    }
	  }
	  else if( // but if "surrounded" by good value
		  (pflag[i+1][j][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)
		  ){
	    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j);
	    // then ASSIGN
	    for(k=0;k<=U3;k++){
	      pv[i][j][k]=ptoavg[i+1][j][k];
	      if(debugfail>=2) dualfprintf(fail_file,"uc2: k=%d pv=%21.15g\n",k,pv[i][j][k]);
	    }
	  }
	  else if( // but if "surrounded" by good value
		  (pflag[i+1][j-1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)
		  ){
	    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j);
	    // then ASSIGN
	    for(k=0;k<=U3;k++){
	      pv[i][j][k]=ptoavg[i+1][j-1][k];
	      if(debugfail>=2) dualfprintf(fail_file,"uc2: k=%d pv=%21.15g\n",k,pv[i][j][k]);
	    }
	  }
	  else if( // but if "surrounded" by good value
		  (pflag[i][j-1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)
		  ){
	    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j);
	    // then ASSIGN
	    for(k=0;k<=U3;k++){
	      pv[i][j][k]=ptoavg[i][j-1][k];
	      if(debugfail>=2) dualfprintf(fail_file,"uc2: k=%d pv=%21.15g\n",k,pv[i][j][k]);
	    }
	  }
	  else if( // but if "surrounded" by good value
		  (pflag[i-1][j-1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)
		  ){
	    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j);
	    // then ASSIGN
	    for(k=0;k<=U3;k++){
	      pv[i][j][k]=ptoavg[i-1][j-1][k];
	      if(debugfail>=2) dualfprintf(fail_file,"uc2: k=%d pv=%21.15g\n",k,pv[i][j][k]);
	    }
	  }
	  else if( // but if "surrounded" by good value
		  (pflag[i-1][j][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)
		  ){
	    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j);
	    // then ASSIGN
	    for(k=0;k<=U3;k++){
	      pv[i][j][k]=ptoavg[i-1][j][k];
	      if(debugfail>=2) dualfprintf(fail_file,"uc2: k=%d pv=%21.15g\n",k,pv[i][j][k]);
	    }
	  }
	  else if( // but if "surrounded" by good value
		  (pflag[i-1][j+1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)
		  ){
	    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j);
	    // then ASSIGN
	    for(k=0;k<=U3;k++){
	      pv[i][j][k]=ptoavg[i-1][j+1][k];
	      if(debugfail>=2) dualfprintf(fail_file,"uc2: k=%d pv=%21.15g\n",k,pv[i][j][k]);
	    }
	  }
	  else if( // but if "surrounded" by good value
		  (pflag[i][j+1][FLAGUTOPRIMFAIL]==UTOPRIMNOFAIL)
		  ){
	    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d : utoprim corrected7\n",t,startpos[1]+i,startpos[2]+j);
	    // then ASSIGN
	    for(k=0;k<=U3;k++){
	      pv[i][j][k]=ptoavg[i][j+1][k];
	      if(debugfail>=2) dualfprintf(fail_file,"uc2: k=%d pv=%21.15g\n",k,pv[i][j][k]);
	    }
	  }


	  /////////////
	  //
	  // no good values.  Choices: STATIC, RESET to something, and AVERAGE failed (t-dt surrounding) values.
	  //
	  /////////////
	  else if(TODONOGOOD==NOGOODAVERAGE){// if no solution, revert to normal observer and averaged densities
	    // average out densities+velocity
	    for(k=0;k<=U3;k++) pv[i][j][k]=0.5*(AVG4_1(ptoavg,i,j,k)+AVG4_2(ptoavg,i,j,k));
	  }
	  else if(TODONOGOOD==NOGOODSTATIC){// if no solution, revert to normal observer and averaged densities
	    // don't change -- stay at previous timestep value (or whatever utoprim() left it as).
	  }
	  else if(TODONOGOOD==NOGOODRESET){// if no solution, revert to normal observer and averaged densities

	    // model-dependent: assumes failures occur mostly in jet region near where matter is mostly freefalling near black hole where b^2/rho_0>>1

	    if(debugfail>=2) dualfprintf(fail_file,"t=%21.15g : i=%d j=%d : utoprim corrected9\n",t,startpos[1]+i,startpos[2]+j);
	    // otherwise some surrounding solutions are also bad, so just keep static
	    // keeping static can lead to large regions of high gamma that are static and unstable to interactions.  Thus, reset static to normal observer flow for safety.
	    // setting to normal observer for surrounding high gamma flows leads to large shocks and high internal energy.
	    //
	    // thus not sure what to do. (should have gotten correct U in first place! -- that's what you should do)
	    //
	    // limit_gamma?  Could also reset v completely to normal observer
	    // if(limit_gamma(1.0,pv[i][j],&geom,-1.0)>=1)
	    //	  FAILSTATEMENT("fixup.c:fixup_utoprim()", "limit_gamma()", 1);

	    // average out densities
	    for(k=0;k<=UU;k++) pv[i][j][k]=0.5*(AVG4_1(ptoavg,i,j,k)+AVG4_2(ptoavg,i,j,k));
	    // make velocity the normal observer
	    get_geometry(i,j,CENT,&geom);
	    set_atmosphere(0,WHICHVEL,&geom,prguess);
	    for(k=U1;k<=U3;k++) pv[i][j][k]=prguess[k];
	  }



	  /////////////
	  //
	  // check new gamma to make sure smaller than original
	  //
	  /////////////
#if(WHICHVEL==VELREL4)
	  if(limit_gamma(gamma,pv[i][j],&geom,-1.0)>=1)
	    FAILSTATEMENT("fixup.c:fixup()", "limit_gamma()", 2);
	  // check gamma
	  if(debugfail>=2){
	    MYFUN(gamma_calc(pv[i][j],&geom,&gamma),"fixup_utoprim: gamma calc failed\n","fixup.c",2);
	    dualfprintf(fail_file,"final gamma: %21.15g\n",gamma);
	  }
#endif

	}
      }// end if mypflag>0 (i.e. failure)

      /////////////////////////////////
      //
      // ACCOUNTING
      //

      // account for changes by tracking conserved quantities
      // note that if new utoprim solution was impossible, we are using here prior solution as new state, which means the diagnostics won't be acccurate.  There's no way to fix this unless the fluxes always give a well-defined new primitive variable.
      // specific value of the mypflag>0 communicates the type of failure
      if(mypflag==UTOPRIMNOFAIL){
	utoprimfailtype=-1;
      }
      else if(
	      (mypflag==UTOPRIMFAILCONV)|| // only used by 5D method currently
	      (mypflag==UTOPRIMFAILCONVUTSQ)|| // rest are only used by 1D/2D method currently
	      (mypflag>=UTOPRIMFAILCONVRET)||
	      (mypflag==UTOPRIMFAILCONVW)||
	      (mypflag==UTOPRIMFAILCONVUTSQ2)
	      ){
	utoprimfailtype=COUNTUTOPRIMFAILCONV;
      }
      else if(mypflag==UTOPRIMFAILRHONEG){
	utoprimfailtype=COUNTUTOPRIMFAILRHONEG;
      }
      else if(mypflag==UTOPRIMFAILUNEG){
	// whether to count uneg as failure in diagnostic reporting or not
	// should really have a new diagnostic for substep u<0 's.
	if((STEPOVERNEGU)&&(!finalstep)){
	  utoprimfailtype=-1;
	}
	else utoprimfailtype=COUNTUTOPRIMFAILUNEG;
      }
      else if(mypflag==UTOPRIMFAILRHOUNEG){
	utoprimfailtype=COUNTUTOPRIMFAILRHOUNEG;
      }
      else if(mypflag==UTOPRIMFAILGAMMAPERC){
	utoprimfailtype=COUNTGAMMAPERC;
      }
      else if(mypflag==UTOPRIMFAILUPERC){
	utoprimfailtype=COUNTUPERC;
      }
      else{
	dualfprintf(fail_file,"No such pflag failure type: %d\n",mypflag);
	myexit(1);
      }
      if(utoprimfailtype!=-1){
	// diagnostics
	diag_fixup(pr0, pv[i][j], &geom, finalstep,utoprimfailtype);
	////////////////
	//
	// reset true pflag counter to "no" (fixed) failure
	//
	////////////////
	pflag[i][j][FLAGUTOPRIMFAIL]=UTOPRIMFAILFIXED;
      }

#if(DOSUPERDEBUG)
      // only used to keep track of normal non-failure points
      // some non-utoprim-failure points will be picked up by this.  I.E. those with limitgammaact/flooract/inflowact
      if(utoprimfailtype==-1) superdebug(pr0,pv[i][j],&geom,utoprimfailtype);
      // collect values for non-failed and failed zones
#endif
    }// end over CZLOOP loop
  }// end if doing avg-type fixup of failures
  return(0);
}



int superdebug(FTYPE *pr0, FTYPE *pr, struct of_geom *ptrgeom, int whocalled)
{
  int k;
  struct of_state q;
  FTYPE Ui[NPR],Uf[NPR];
  FTYPE X[NDIM],r,th;
  FTYPE ftemp;
  int failreturn;
  static FILE * super_fail_file;
  static int firsttime=1;
  int output;
  static int countoutput=0;


  if(firsttime){
    super_fail_file=fopen("superdebug.out","wt");
    if(super_fail_file==NULL){
      dualfprintf(fail_file,"Cannot open superdebug.out\n");
      myexit(1);
    }
    countoutput=0;
  }

  if(whocalled==-1){
    // then only output if every so often, randomly in step
    if(ranc(0)>0.99) output=1;
    else output=0;
  }
  else output=1;

  if(output){
    // before any changes
    failreturn=get_state(pr0,ptrgeom,&q);
    if(failreturn>=1) dualfprintf(fail_file,"get_state(1) failed in fixup.c, why???\n");
    failreturn=primtoU(UDIAG,pr0,&q,ptrgeom,Ui);
    if(failreturn>=1) dualfprintf(fail_file,"primtoU(1) failed in fixup.c, why???\n");
    
    // after any changes
    failreturn=get_state(pr,ptrgeom,&q);
    if(failreturn>=1) dualfprintf(fail_file,"get_state(2) failed in fixup.c, why???\n");
    failreturn=primtoU(UDIAG,pr,&q,ptrgeom,Uf);
    if(failreturn>=1) dualfprintf(fail_file,"primtoU(2) failed in fixup.c, why???\n");


    coord(ptrgeom->i, ptrgeom->j, ptrgeom->p, X);
    bl_coord(X, &r, &th);
    
    
    // used to determine nature of pre and post failure quantities
    fprintf(super_fail_file,"%21.15g %ld %d %d %d ",t,realnstep,steppart,numstepparts,whocalled);
    fprintf(super_fail_file,"%21.15g %21.15g %d %d ",r,th,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j);
    PDUMPLOOP fprintf(super_fail_file,"%21.15g %21.15g %21.15g %21.15g ",pr0[k],pr[k],Ui[k],Uf[k]);
    fprintf(super_fail_file,"\n");
    if(!(countoutput%1000)) fflush(super_fail_file);
    countoutput++;
  }



  firsttime=0;

  return(0);
}



int set_density_floors(struct of_geom *ptrgeom, FTYPE *pr, FTYPE *prfloor)
{
  struct of_state q;
  FTYPE U[NPR];
  FTYPE Upbinf[NPR];
  FTYPE r,th,X[NDIM];
  FTYPE prcoef[NPR];



  coord(ptrgeom->i, ptrgeom->j, ptrgeom->p, X);
  bl_coord(X, &r, &th);

  if(BCtype[X1UP]==FIXEDOUTFLOW){ // then doing bondi inflow
    // avoids constant floor activation -- trying to be more physical
    prcoef[RHO]=RHOMIN/100.0;
    prcoef[UU]=UUMIN/100.0;
  }
  else{
    prcoef[RHO]=RHOMIN;
    prcoef[UU]=UUMIN;
  }

  ////////////////////
  // scaling functions
  ////
  if(1){ // choice
    if(rescaletype==0){
      // to be like Bondi

      // even if d/dt->0 and d/dphi->0, this doesn't conserve E and L along magnetic flux lines in the funnel region.  This is bad!
      prfloor[RHO] = prcoef[RHO]*pow(r, -1.5); 
      prfloor[UU] = prcoef[UU]*pow(r,-2.5);
    }
    else if(rescaletype==1){
      // to conserve E and L along magnetic lines, have b^2/\rho\sim const.
      prfloor[RHO] = prcoef[RHO]*pow(r, -2.7);
      prfloor[UU] = prcoef[UU]*pow(r,-2.7) ;
    }
    else if(rescaletype==2){
      // only set absolute limits and let density be set later

      // absolute limits determined by what density contrast can be handled by flux calculation.  This is at least limited by machine precision, but more likely an instability is generated for some  ratio of densities that's much smaller than machine precision.

      // we should fix absolute floor so that this point is not too much different than surrounding points.
      // this is difficult because process depends on how done, so just set reasonable minimum

      // Bondi like
      // coefficient should be small enough so that at outer radius the energy per baryon at infinity limit of (say) 100 can be maintained (b^2 \sim 0.01 r^{-2.7} for \rho=1 density maximum)
      // these will do for rout=10^4
      //      prfloor[RHO] = 1E-12*pow(r, -1.5);
      //prfloor[UU] = 1E-14*prfloor[RHO]/r ;
      


      // to best conserve E and L along magnetic field lines
      if (get_state(pr, ptrgeom, &q) >= 1)
	FAILSTATEMENT("fixup.c:set_density_floors()", "get_state() dir=0", 1);

      MYFUN(primtoU(UDIAG,pr, &q, ptrgeom, U),"fixup.c:set_density_floors()", "primtoU()", 1);

      // now have U[UU] and U[PH]
      // note that U[UU]/(gdet*rho*u^t) is conserved along field lines, even in non-ideal MHD, where A_{\phi} is still a good stream function in non-ideal MHD.
      // same in terms of U[PH]
      // this is really a coordinate dependent quantity, but if field lines connect to infinity, then this is the same at infinity.

      // Conserved quantity per baryon (or rest-mass unit)
      //      PLOOP Upbinf=U[k]/(ptrgeom->g * pr[RHO]*(q.ucon[TT]));
      //Upbinf[UU]*=-1.0; // time component has - sign in this -+++ signature code

      // could inject mass or enthalpy, we choose mass for now since rho>h typically (except at very edge of torus)
      prfloor[RHO]=-U[UU]/(ptrgeom->g * GAMMAMAX * (q.ucon[TT]));
      prfloor[UU]=prfloor[RHO]*0.01; // cold injection

    }
  }
  else{
    prfloor[RHO] = RHOMIN*pow(r, -2.0);
    prfloor[UU] = UUMIN*prfloor[RHO] ;
  }






  // some old attempts
    /*
    if(0){ // choice
      ftempA=(RHOMAX+RHOMIN)*0.5/RHOMIN;
      ftempB=(RHOMAX-RHOMIN)*0.5/RHOMIN;
      // choice
      if(1) rhotscal = ftempA+ftempB*cos(2.0*M_PI*X[2]);
      else  rhotscal = ftempA+ftempB*cos(2.0*M_PI*th);
      ftempA=(UUMAX+UUMIN)*0.5/UUMIN;
      // choice (make same as above)
      ftempB=(UUMAX-UUMIN)*0.5/UUMIN;
      if(1) uutscal = ftempA+ftempB*cos(2.0*M_PI*X[2]);
      else  uutscal = ftempA+ftempB*cos(2.0*M_PI*th);
    }

    //    else{
    if(1){ 
      rhotscal = 1.0;
      uutscal = 1.0;
     }


     prfloor[UU] = UUMIN*uurscal*uutscal;
     prfloor[RHO] = RHOMIN*rhorscal*rhotscal;
    */


  return(0);
}


// for each WHICHVEL possibility, set atmosphere state for any coordinate system
// which=0 : initial condition
// which=1 : evolution condition (might also include a specific angular momentum or whatever)
// which==1 assumes pr set to something locally reasonable, and we adjust to that slowly

#define TAUADJUSTATM (10.0) // timescale for boundary to adjust to using preset inflow
int set_atmosphere(int whichcond, int whichvel, struct of_geom *ptrgeom, FTYPE *pr)
{
  FTYPE rho,u,ur,uh,up;
  FTYPE X[NDIM];
  FTYPE r,th;
  FTYPE prlocal[NPR];
  int k;

  // Bondi like initial atmosphere
  //    rho = RHOMIN * 1.E-14;
  //    u = UUMIN * 1.E-14;
  coord(ptrgeom->i, ptrgeom->j, ptrgeom->p, X);
  bl_coord(X, &r, &th);

  // default
  PLOOP prlocal[k]=pr[k];

#if(EOMTYPE!=EOMFFDE)
  // Bondi-like atmosphere
  prlocal[RHO] = RHOMIN*pow(r,-1.5);
  prlocal[UU]  = UUMIN*pow(r,-2.5);
#else
  prlocal[RHO] = 0;
  prlocal[UU]  = 0;
#endif
    
  // bl-normal observer (4-vel components)
  
  // normal observer velocity in atmosphere
  if(whichvel==VEL4){
    prlocal[U1] = -ptrgeom->gcon[0][1]/sqrt(-ptrgeom->gcon[0][0]) ;
    prlocal[U2] = -ptrgeom->gcon[0][2]/sqrt(-ptrgeom->gcon[0][0]) ;
    prlocal[U3] = -ptrgeom->gcon[0][3]/sqrt(-ptrgeom->gcon[0][0]) ;
  }
  else if(whichvel==VEL3){
    prlocal[U1] = ptrgeom->gcon[0][1]/ptrgeom->gcon[0][0] ;
    prlocal[U2] = ptrgeom->gcon[0][2]/ptrgeom->gcon[0][0] ;
    prlocal[U3] = ptrgeom->gcon[0][3]/ptrgeom->gcon[0][0] ;
    // GAMMIE
    //ur = -1./(r*r);
    //uh=up=0.0;
  }
  else if(whichvel==VELREL4){
    prlocal[U1] = 0.0;
    prlocal[U2] = 0.0;
    prlocal[U3] = 0.0;
  }
  
  if(whichcond==1){
    if(100.0*dt>TAUADJUSTATM){
      dualfprintf(fail_file,"dt=%21.15g and TAUADJUSTATM=%21.15g\n",dt,TAUADJUSTATM);
      myexit(1);
    }
    // TAUADJUSTATM must be >> dt always in order for this to make sense (i.e. critical damping to fixed solution)
    PLOOP pr[k] = pr[k]+(prlocal[k]-pr[k])*dt/TAUADJUSTATM;
  }
  else if(whichcond==0){ 
    PLOOP pr[k] = prlocal[k];
    // very specific
    // always outflow field
    //    pr[B1] = pr[B2] = pr[B3] = 0;
  }
  else if(whichcond==-1){ 
    // t=0, just sets as "floor"
    PLOOP pr[k] = prlocal[k];
    pr[B1] = pr[B2] = pr[B3] = 0;
  }


  return(0);

}



#define FLOORDAMPFRAC (0.1)
#define NUMBSQFLAGS 5

int get_bsqflags(int stage, FTYPE (*pv)[N2M][NPR])
{
  int i,j;
  FTYPE bsq;
  struct of_geom geom;
  struct of_state q;
  FTYPE gamma,X[NDIM],r,th;
  FTYPE prfloor[NPR];
  int flags[NUMBSQFLAGS]={0};

  if((CHECKSOLUTION==0)&&(LIMADJUST==0)&&(FLUXADJUST==0)) return(0); // otherwise do it

  DQZLOOP { // over range dq is computed

    get_geometry(i, j, CENT, &geom);


#if(EOMTYPE==EOMGRMHD)
    // b^2
    if (get_state(pv[i][j], &geom, &q) >= 1)
      FAILSTATEMENT("fixup.c:get_bsqflags()", "get_state()", 1);
    bsq = dot(q.bcon, q.bcov);
    // initial
    pflag[i][j][FLAGBSQORHO]=0;

    // now do checks 10/31/2003 : constants fixed based upon accretion
    // models and locations of disk, corona, and funnel (b^2/rho>1)

    // b^2/\rho
    if(bsq/pv[i][j][RHO]>BSQORHOLIMIT) flags[0]=2;
    else if(bsq/pv[i][j][RHO]>0.9*BSQORHOLIMIT) flags[0]=1;
    else flags[0]=0;

    // b^2/u

    if(bsq/pv[i][j][UU]>BSQOULIMIT) flags[3]=2;
    else if(bsq/pv[i][j][UU]>0.9*BSQOULIMIT) flags[3]=1;
    else flags[3]=0;
#endif

#if(WHICHVEL==VELREL4)
    // gamma
    if(gamma_calc(pv[i][j],&geom,&gamma)>=1){
      dualfprintf(fail_file,"gamma calc failed: get_bsqflags\n");
      if (fail(FAIL_UTCALC_DISCR) >= 1)
        return (1);
    }
    // \gamma relative to 4-velocity
    if(gamma>2.0*GAMMADAMP) flags[1]=2;
    else if(gamma>GAMMADAMP) flags[1]=1;
    else flags[1]=0;
#endif
    // these density mod on limiters may be too much diffusivity in otherwise ok regions

    


#if(0)
    // density floors
    set_density_floors(&geom,pv[i][j],prfloor);
    // rho/rho_{fl}    
    // GODMARK, 0.1 here    
    if(prfloor[RHO]/pv[i][j][RHO]>FLOORDAMPFRAC) flags[2]=2;
    else if(prfloor[RHO]/pv[i][j][RHO]>0.1*FLOORDAMPFRAC) flags[2]=1;
    else flags[2]=0;

    // u/u_{fl}

    if(prfloor[UU]/pv[i][j][UU]>FLOORDAMPFRAC) flags[4]=2;
    else if(prfloor[UU]/pv[i][j][UU]>0.1*FLOORDAMPFRAC) flags[4]=1;
    else flags[4]=0;
#endif


    // now check our flag state
    pflag[i][j][FLAGBSQORHO]=flags[0];
    pflag[i][j][FLAGBSQOU]=flags[3];


#if(LIMADJUST>0)
    // now check our flag state
    if((flags[0]==2)||(flags[1]==2)||(flags[2]==2)||(flags[3]==2)||(flags[4]==2)){ pflag[i][j][FLAGREALLIM]=lim-1; }
    else if((flags[0]==1)||(flags[1]==1)||(flags[2]==1)||(flags[3]==1)||(flags[4]==1)){ pflag[i][j][FLAGREALLIM]=lim-1; }
    else{pflag[i][j][FLAGREALLIM]=lim;}
    if(pflag[i][j][FLAGREALLIM]<0) pflag[i][j][FLAGREALLIM]=0;
#endif

#if(FLUXADJUST>0)
    // now check our flag state
    if((flags[0]==2)||(flags[1]==2)||(flags[2]==2)||(flags[3]==2)||(flags[4]==2)){ pflag[i][j][FLAGREALFLUX]=fluxmethod-1;}
    else if((flags[0]==1)||(flags[1]==1)||(flags[2]==1)||(flags[3]==1)||(flags[4]==1)){ pflag[i][j][FLAGREALFLUX]=fluxmethod; }
    else{pflag[i][j][FLAGREALFLUX]=fluxmethod;}
    if(pflag[i][j][FLAGREALFLUX]<0) pflag[i][j][FLAGREALFLUX]=0;
#endif


  }
  return(0);
}


int limit_gamma(FTYPE gammamax, FTYPE*pr,struct of_geom *ptrgeom,int finalstep)
{
  FTYPE f,gamma,pref;
  FTYPE pr0[NPR];
  int k;
  FTYPE r,th,X[NDIM];
  FTYPE realgammamax;

  coord(ptrgeom->i,ptrgeom->j,ptrgeom->p,X);
  bl_coord(X,&r,&th);


#if(GAMMAERGOLIMIT)
  // force flow to not move too fast inside ergosphere
  if(r<2) realgammamax=3;
  else realgammamax=GAMMAMAX;
#else
  realgammamax=gammamax;
  if(realgammamax<=1.0) return(0); // nothing to do
#endif



  //PLOOP{ dualfprintf(fail_file,"i=%d j=%d pr[%d]=%21.15g\n",startpos[1]+i,startpos[2]+j,k,pv[i][j][k]);}
  if(gamma_calc(pr,ptrgeom,&gamma)>=1){
    dualfprintf(fail_file,"limit_gamma: gamma calc failed\n");
    dualfprintf(fail_file,"i=%d j=%d oldgamma=%21.15g\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,gamma);
    if (fail(FAIL_UTCALC_DISCR) >= 1)
      return (1);
  }
    
  if((gamma > realgammamax && (gamma!=1.0))) {    
    // set pre-primitive
    PLOOP    pr0[k]=pr[k];

    if(debugfail>=2) dualfprintf(fail_file,"i=%d j=%d oldgamma=%21.15g\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,gamma);
    /* rescale velocities to reduce gamma 
     * to realgammamax */
    pref=(realgammamax*realgammamax - 1.)/(gamma*gamma - 1.);
    if(pref<0.0){
      dualfprintf(fail_file,"limit_gamma: pref calc failed pref=%21.15g\n",pref);
      dualfprintf(fail_file,"i=%d j=%d oldgamma=%21.15g\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,gamma);
      if (fail(FAIL_UTCALC_DISCR) >= 1)
	return (1);
    }

    f = sqrt(pref);
    pr[U1] *= f ;	
    pr[U2] *= f ;	
    pr[U3] *= f ;
    
    diag_fixup(pr0, pr, ptrgeom, finalstep,COUNTLIMITGAMMAACT);

    return(-1);// indicates did change primitive
  }
  return(0); // indicates didn't change anything
}

// only correct for polar axis at both inner and outer x2/theta grid edge.
void fix_flux(FTYPE (*pb)[N2M][NPR],FTYPE F1[][N2M][NPR], FTYPE F2[][N2M][NPR])
{
        int i,j,k ;
        FTYPE sth ;
        FTYPE X[NDIM],r,th ;

	// this has nothing to deal with MPI-boundaries, so ok as is
	// only applies for polar axis
	if(mycpupos[2]==0){	 
	  if(BCtype[X2DN]==POLARAXIS){
	    for(i=-1;i<N1;i++) {
	    // emf
	      F1[i][-1][B2] = -F1[i][0][B2] ;
	      // all should be 0 except kinetic energy flux
	      PLOOP if(k!=U2) F2[i][0][k] = 0. ;
	    }
	  }
	}
	if(mycpupos[2]==ncpux2-1){
	  if(BCtype[X2UP]==POLARAXIS){
	    for(i=-1;i<N1;i++) {
	      // emf
	      F1[i][N2][B2] = -F1[i][N2-1][B2] ;	    
	      PLOOP if(k!=U2) F2[i][N2][k] = 0. ;
	    }
	  }
	}

	LOOPF2 LOOPF1{
	  if((i>=0)&&(i<N1)&&( (abs(startpos[2]+j-(totalsize[2]-1)/2)<=0) ) ){
	    PLOOP if(k!=U2) F2[i][j][k] = 0. ;
	  }
	}

	// avoid mass flux in wrong direction, so consistent with velocity fix
	// how to treat other fluxes?
	// GODMARK
	// this probably isn't right for general boundary condition size
	if(mycpupos[1]==0){	 
	  if(BCtype[X1DN]==OUTFLOW){
	    for(j=-1;j<N2;j++) {
	      if(F1[0][j][RHO]>0) F1[0][j][RHO]=0;
	    }
	  }
	}
	if(mycpupos[1]==ncpux1-1){
	  if(BCtype[X1UP]==OUTFLOW){
	    for(j=-1;j<N2;j++) {
	      if(F1[N1][j][RHO]<0) F1[N1][j][RHO]=0;
	    }
	  }
	}

}



#define UTLIMIT (50.0) // limit of u^t given no other info and need to fix u
#define UTFIX (utobs) //  if failed u^t, then fix to this level to keep normal, or use local value as desired value
#define GRADIENTFACTOR (.001) // how to set?

// pr is pr to be fixed
// prmodel is pr that has a model ucon[TT], since we want to interpolate based upon ucon[TT], not just a random one when inside step_ch.c in fluxcalc()
// prmodel is just pr if in utoprim.c since we want a real change, not just an interpolation, thus will not be used as basis
int check_pr(FTYPE *pr,FTYPE *prmodel, struct of_geom *ptrgeom,int modelpos, int finalstep)
{
  int failedcheck;
  FTYPE ucon[NPR],uconmodel[NPR];
  FTYPE prold[NPR],probs[NPR];
  FTYPE gradient[NDIM],normsq;
  int trialcount,ntrials;
  int i,j,k;
  FTYPE lastuttdiscr,uttdiscr0,utobs;
  FTYPE dampfactor,dampfactorchange;
  FTYPE newerr,olderr;
  int method;
  struct of_geom modelgeom;
  struct of_geom *ptrmodelgeom;
  int idel,jdel;
  FTYPE realutlimit,realdiscrlimit;
  FTYPE modeldiscr;
  int utinterp;
  FTYPE toldiscr;
  int modeli,modelj;
  int bctype;

  FTYPE pr0[NPR];



  if(modelpos==-100){ // then utoprim called us from usrfun
    modelpos=0;
    ntrials=30.0; // don't try so hard since failure is more likely, leads to static solution
  }
  else{
    ntrials=200.0; // try really hard since using observer as solution is not fun
  }
  toldiscr=1.E-2;
  dampfactorchange=0.5;
  dampfactor=1.0;
  method=1; // 0=GRADIENTFACTOR 1=damp newton's method
  utinterp=0; // whether to fix interpolation (if bad) based upon a model pr
  

  if(method==1){
    // save old pr
    PLOOP probs[k]=prold[k]=pr[k];
  }

  whocalleducon=1; // turn off failures
  if(ucon_calc(pr, ptrgeom, ucon) >= 1) ucon[TT]=1E30; // bad, so keep going
  uttdiscr0=lastuttdiscr=uttdiscr; // ucon_calc() set uttdiscr
  if(ucon[TT]<UTLIMIT) return(0); // good, so just continue calculations


  /////////////////////////////////////////////////////
  //
  // if here, need to fix
  //
  // set pre-primitive
  PLOOP    pr0[k]=pr[k];



  if(modelpos==-2) return(1); // don't try to correct, just fail since ucon[TT] wasn't less than the limit (otherwise wouldn't reach here).  Used to just check, no fix.
  
  if(modelpos==-3){
    modelpos=-1;
    bctype=1; // so bound_prim called us
  }
  else bctype=0; // non-bound call
  

  //////////////////
  //
  // if we are here, original velocities need to be corrected

  // first find normal observer to get relevant u^t
  for(i=1;i<=3;i++){
    probs[U1+i-1] = ptrgeom->gcon[TT][i]/ptrgeom->gcon[TT][TT] ;
  }
  // ucon[TT] must be good for normal observer (?) hard to solve for u^t for normal observer in general
  // change of whocalleducon forces kill level check on normal observer
  whocalleducon=0;
  if(ucon_calc(probs, ptrgeom, ucon) >= 1){
    dualfprintf(fail_file,"Thought normal observer had to have good u^t!\n");
    return(1);
  }
  else utobs=ucon[TT];
  whocalleducon=1;


  if(utinterp&&(modelpos>=0)){ // use modelpos==-1 for no use of model
    // below for no interpolation but use of model
    if(modelpos==0){
      modeli=ptrgeom->i;
      modelj=ptrgeom->j;
      ptrmodelgeom=ptrgeom;
    }
    // modelpos>=1 for interpolations in fluxcalc() (model value on center)
    else if(modelpos==1){
      if(ptrgeom->p==FACE1){ idel=1; jdel=0; }
      else if(ptrgeom->p==FACE2){ idel=0; jdel=1; }
      modeli=(ptrgeom->i) -idel;
      modelj=(ptrgeom->j) -jdel;
      // then i-1,j is center position
      get_geometry(modeli,modelj,CENT,&modelgeom);
      ptrmodelgeom=&modelgeom;
    }
    else if(modelpos==2){
      modeli=ptrgeom->i;
      modelj=ptrgeom->j;
      // then i,j is center position
      get_geometry(modeli,modelj,CENT,&modelgeom);
      ptrmodelgeom=&modelgeom;
    }
    // determine model u^t
    if(ucon_calc(prmodel, ptrmodelgeom, uconmodel) >= 1){
      // model no good
      if(bctype==0){
	dualfprintf(fail_file,"serious failure.  On-grid values and fixed bc values used have u^t imaginary: modeli: %d modelj: %d uttdiscr: %21.15g\n",startpos[1]+modeli,startpos[2]+modelj,uttdiscr);
	whocalleducon=0; // turn on failures
	if (fail(FAIL_UTCALC_DISCR) >= 1)
	  return (1);
      }
      else uconmodel[TT]=realutlimit=1E30;
      // otherwise normal to sometimes encounter failure if using model in bc (which isn't currently)
    }
    else realutlimit=uconmodel[TT]; // model ut
    modeldiscr=uttdiscr;

    // upper limit (if model has large u^t, assume ok to have one)
    if(realutlimit>UTLIMIT) realutlimit=UTLIMIT;
  }
  else realutlimit=UTFIX; // no model, just fix based upon normal observer since no idea what is "ok" to have

  realdiscrlimit=1.0/(realutlimit*realutlimit);
  newerr=olderr=fabs(lastuttdiscr-realdiscrlimit)/realdiscrlimit;


  // LOOP


  // otherwise need to fix
  failedcheck=0;

  trialcount=0;
  while( ((newerr>toldiscr)&&(method==1)) ||((ucon[TT]>realutlimit)&&(method==0)) ){
    // see if we can fix things since outside limits

    // determine gradient
    normsq=0.0;
    for(i=1;i<=3;i++){
      gradient[i]=2.0*(ptrgeom->gcov[0][i]);
      for(j=1;j<=3;j++){
	// note that ucon is the same as pr here since ucon_calc sets spatial terms to pr
	gradient[i]+=2.0*ucon[j]*ptrgeom->gcov[i][j];
      }
      normsq+=gradient[i]*gradient[i];
    }
    // normalize gradient
    for(i=1;i<=3;i++){
      gradient[i]/=sqrt(normsq);
    }
    // save old pr and change new one    
    if(method==0){
      for(i=1;i<=3;i++){
	
	pr[U1+i-1]-=gradient[i]*GRADIENTFACTOR*((FTYPE)(ntrials)+1.0)/((FTYPE)(ntrials)-(FTYPE)(trialcount)+1.0);
	//pr[U1+i-1]-=gradient[i]*GRADIENTFACTOR;
      }
    }
    else if(method==1){
      for(i=1;i<=3;i++){
	prold[U1+i-1]=pr[U1+i-1];
	if(realdiscrlimit-uttdiscr>0)	pr[U1+i-1]-=gradient[i]*dampfactor;
	else 	pr[U1+i-1]+=gradient[i]*dampfactor;
      }
    }
    // get new ucon[TT], is it ok now?
    if(ucon_calc(pr, ptrgeom, ucon) >= 1) ucon[TT]=1E30; // bad bad, keep going
    newerr=fabs(uttdiscr-realdiscrlimit)/realdiscrlimit;
    if((method==1)&&(newerr>=olderr)){
      // then went too far (if going in right direction at all)
      dampfactor*=dampfactorchange;
      if(dampfactor<1E-10){
	if((fabs(ucon[TT]-realutlimit)/realutlimit)<0.5) break; // just be happy you got out alive
	else{
	  failedcheck=1;
	  if(debugfail>=1) dualfprintf(fail_file,"dumpfactor reached min\n");
	  break;
	}
      }
      // revert to old pr and start again
      for(i=1;i<=3;i++){
	pr[U1+i-1]=prold[U1+i-1];
      }
    }
    else{
      trialcount++; // only iterate if good step
      olderr=newerr;
    }
    if(debugfail>=2) {
      if((myid==0)&&(ptrgeom->i==117)&&(ptrgeom->j==-1)){
	dualfprintf(fail_file,"trialcount=%d uttdiscr0=%21.15g uttdiscr=%21.15g newerr: %21.15g dampfactor=%21.15g\n",trialcount,uttdiscr0,uttdiscr,newerr,dampfactor);
      }
    }
    // even if not bad, could still be too large, so check
    if(trialcount==ntrials){
      if((fabs(ucon[TT]-realutlimit)/realutlimit)<0.5) break; // just be happy you got out alive
      else{
	failedcheck=1;
	if(debugfail>=1) dualfprintf(fail_file,"number of trials reached max\n");
	break;
      }
    }
  }
  whocalleducon=0; // turn back on failures

  if(failedcheck){
    if(debugfail>=1) dualfprintf(fail_file,"couldn't fix ucon[TT]=%21.15g newerr=%21.15g uttdiscr=%21.15g\n",ucon[TT],newerr,uttdiscr);
    if(debugfail>=1) dualfprintf(fail_file,"check_pr failure: t=%21.15g , couldn't fix ucon: i=%d j=%d p=%d ucon[TT]=%21.15g realutlimit=%21.15g uconmodel[TT]=%21.15g\n",t,startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,ptrgeom->p,ucon[TT],realutlimit,uconmodel[TT]);
    if(debugfail>=1) dualfprintf(fail_file,"uttdiscr0=%21.15g uttdiscr=%21.15g realdiscrlimit=%21.15g modeldiscr=%21.15g dampfactor=%21.15g\n",uttdiscr0,uttdiscr,realdiscrlimit,modeldiscr,dampfactor);
    // will still run perhaps with UT>realutlimit, could stop it, but won't for now      
    if(debugfail>=1){
      PLOOP{
	dualfprintf(fail_file,"pr[%d]=%21.15g prmodel[%d]=%21.15g\n",k,pr[k],k,prmodel[k]);
      }
      if(debugfail>=1) dualfprintf(fail_file,"need better algorithm: check_pr failure: couldn't fix ucon: i=%d j=%d p=%d ucon[TT]=%21.15g\n",startpos[1]+ptrgeom->i,startpos[2]+ptrgeom->j,ptrgeom->p,ucon[TT]);
    }

    // force to normal observer solution
    PLOOP pr[k]=probs[k];
  }

  // account for changes
  diag_fixup(pr0, pr, ptrgeom, finalstep,COUNTLIMITGAMMAACT);


  return(0);
}


int inflow_check_4vel(int dir, FTYPE *pr, struct of_geom *ptrgeom, int finalstep)
{
  int ii,jj;
  int iin,iout;
  int jjn,jout;
  FTYPE pr0[NPR];
  int k;


  ii=ptrgeom->i;
  jj=ptrgeom->j;

  if(dir==1){
    // for dir==1
    if((ptrgeom->p==CENT)||(ptrgeom->p==FACE2)){ iin=-1; iout=totalsize[1]; }
    else if((ptrgeom->p==FACE1)||(ptrgeom->p==CORN)){ iin=0; iout=totalsize[1]; }
    else{
      dualfprintf(fail_file,"dir=%d no such location: %d\n",dir,ptrgeom->p);
      myexit(1);
    }
    if( 
       ((startpos[1]+ii<=iin)&&(BCtype[X1DN]==OUTFLOW)&&(pr[U1+dir-1] > 0.)) 
       ||((startpos[1]+ii>=iout)&&(BCtype[X1UP]==OUTFLOW)&&(pr[U1+dir-1] < 0.)) 
       ) {
      // set pre-primitive
      PLOOP    pr0[k]=pr[k];
      pr[U1]=0;
      diag_fixup(pr0, pr, ptrgeom, -1.0,COUNTINFLOWACT);
    }
    if( 
       ((startpos[1]+ii<=iin)&&(BCtype[X1DN]==FIXEDOUTFLOW)&&(pr[U1+dir-1] > 0.)) 
       ||((startpos[1]+ii>=iout)&&(BCtype[X1UP]==FIXEDOUTFLOW)&&(pr[U1+dir-1] < 0.)) 
       ) {
      // set pre-primitive
      PLOOP    pr0[k]=pr[k];
      // then inflow according to Bondi-like atmosphere
      set_atmosphere(1,WHICHVEL,ptrgeom,pr);

      // below never really accounted for since on boundary zones
      diag_fixup(pr0, pr, ptrgeom, -1.0,COUNTINFLOWACT);
    }
  }
  else if(dir==2){
    // for dir==2
    if((ptrgeom->p==CENT)||(ptrgeom->p==FACE2)){ jjn=-1; jout=totalsize[2]; }
    else if((ptrgeom->p==FACE1)||(ptrgeom->p==CORN)){ jjn=0; jout=totalsize[2]; }
    else{
      dualfprintf(fail_file,"dir=%d no such location: %d\n",dir,ptrgeom->p);
      myexit(1);
    }
    if( 
       ((startpos[2]+jj<=jjn)&&(BCtype[X2DN]==OUTFLOW)&&(pr[U1+dir-1] > 0.)) 
       ||((startpos[2]+jj>=jout)&&(BCtype[X2UP]==OUTFLOW)&&(pr[U1+dir-1] < 0.)) 
       ) {
      // set pre-primitive
      PLOOP    pr0[k]=pr[k];
      pr[U2]=0;
      diag_fixup(pr0, pr, ptrgeom, -1.0,COUNTINFLOWACT);
    }
  }
  else return(1); // uh

  return(0);
}


int inflow_check_3vel(int dir, FTYPE *pr, struct of_geom *ptrgeom, int finalstep)
{

  return(inflow_check_4vel(dir,pr,ptrgeom,-1.0));

}

int inflow_check_rel4vel(int dir, FTYPE *pr, struct of_geom *ptrgeom, int finalstep)
{
  FTYPE ucon[NDIM] ;
  int ii,jj;
  int j,k ;
  FTYPE alpha,betacon,gamma,vsq ;
  int iin,iout;
  int jjn,jout;
  int dofix=0;
  FTYPE pr0[NPR];


  ii=ptrgeom->i;
  jj=ptrgeom->j;


  ucon_calc(pr, ptrgeom, ucon) ;

  dofix=0;
  if(dir==1){
    // for dir==1
    if((ptrgeom->p==CENT)||(ptrgeom->p==FACE2)){ iin=-1; iout=totalsize[1]; }
    else if((ptrgeom->p==FACE1)||(ptrgeom->p==CORN)){ iin=0; iout=totalsize[1]; }
    else{
      dualfprintf(fail_file,"dir=%d no such location: %d\n",dir,ptrgeom->p);
      myexit(1);
    }
    if( 
       ((startpos[1]+ii<=iin)&&(BCtype[X1DN]==OUTFLOW)&&(ucon[dir] > 0.)) 
       ||((startpos[1]+ii>=iout)&&(BCtype[X1UP]==OUTFLOW)&&(ucon[dir] < 0.)) 
       ) {
      dofix=1;
    }
  }
  else if(dir==2){
    // for dir==2
    if((ptrgeom->p==CENT)||(ptrgeom->p==FACE2)){ jjn=-1; jout=totalsize[2]; }
    else if((ptrgeom->p==FACE1)||(ptrgeom->p==CORN)){ jjn=0; jout=totalsize[2]; }
    else{
      dualfprintf(fail_file,"dir=%d no such location: %d\n",dir,ptrgeom->p);
      myexit(1);
    }
    if( 
       ((startpos[2]+jj<=jjn)&&(BCtype[X2DN]==OUTFLOW)&&(ucon[dir] > 0.)) 
       ||((startpos[2]+jj>=jout)&&(BCtype[X2UP]==OUTFLOW)&&(ucon[dir] < 0.)) 
       ) {
      dofix=2;
    }
  }
  else return(1); // uh

  if(dofix){
    // set pre-primitive
    PLOOP    pr0[k]=pr[k];

    /* find gamma and remove it from primitives */
    if(gamma_calc(pr,ptrgeom,&gamma)>=1){
      dualfprintf(fail_file,"gamma calc failed: inflow_check_rel4vel\n");
      if (fail(FAIL_UTCALC_DISCR) >= 1)
        return (1);
    }
    pr[U1] /= gamma ;
    pr[U2] /= gamma ;
    pr[U3] /= gamma ;
    alpha = 1./sqrt(-ptrgeom->gcon[0][0]) ;
    
    /* reset radial velocity so radial 4-velocity
     * is zero */
    if(dofix==1){
      betacon = ptrgeom->gcon[0][1]*alpha*alpha ;
      pr[U1] = betacon/alpha ; // gives 3-vel contravariant
    }
    else if(dofix==2){
      betacon = ptrgeom->gcon[0][2]*alpha*alpha ;
      pr[U2] = betacon/alpha ; // gives 3-vel contravariant
    }
    /* now find new gamma and put it back in */
    vsq = 0. ;
    SLOOP vsq += ptrgeom->gcov[j][k]*pr[U1+j-1]*pr[U1+k-1] ;

    if(vsq<0.0){
      if(vsq>-NUMEPSILON*100.0) vsq=0.0; // machine precision thing
      else if (fail(FAIL_VSQ_NEG) >= 1){
        trifprintf("vsq=%21.15g\n",vsq);
        return (1);
      }
    }
 
    // it's possible that setting ucon(bc comp)->0 leads to gamma>1.0, so just reduce gamma to GAMMAMAX if that's the case
    if(vsq>=1.0){
      vsq = 1.0-1.0/(GAMMAMAX*GAMMAMAX);
      //dualfprintf(fail_file,"i=%d j=%d inflow limit required change in gamma (after dofix==%d): vsq=%21.15g\n",ii+startpos[1],jj+startpos[2],dofix,vsq);
    }

    gamma = 1./sqrt(1. - vsq) ;
    pr[U1] *= gamma ;
    pr[U2] *= gamma ;
    pr[U3] *= gamma ;

    // only for boundary conditions, not active zones, hence -1.0 instead of finalstep
    diag_fixup(pr0, pr, ptrgeom, -1.0,COUNTINFLOWACT);


    /* done */
    return(-1);
  }
  else return(0);

}
