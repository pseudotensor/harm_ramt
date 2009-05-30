
#include "decs.h"






/* calculate fluxes in direction dir and conserved variable U; these
   are always needed together, so there is no point in calculated the
   stress tensor twice */

// returntype==0 : flux with geometric factor geom->e (used by evolution code)
// returntype==1 : flux with physical geometry factor geom->g (used by diagnostics)
// see UtoU and source_conn()
int primtoflux(int returntype, FTYPE *pr, struct of_state *q, int dir,
	       struct of_geom *geom, FTYPE *flux)
{
  // sizes: NPR,struct of_state, int, struct of_geom, NPR
  int i = 0, j = 0, k = 0;
  FTYPE dualf[NDIM];
  FTYPE fluxinput[NPR];
  int dualfaradayspatial_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *dualf);
  int massflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *massflux);
  int entropyflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *entropyflux);
  void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);


  massflux_calc(pr, dir, q, &fluxinput[RHO]); // fills RHO only


  // GODMARK WTF!  Problems with code (compiling?) with this
  // if(mhd_calc(pr,dir,geom,q,mhd)>=1)
  // FAILSTATEMENT("phys.c:primtoflux()","mhd_calc() dir=1or2",1);

  // MHD stress-energy tensor w/ first index up, second index down.
  mhd_calc(pr, dir, geom, q, &fluxinput[UU]); // fills UU->U3

  dualfaradayspatial_calc(pr,dir,q,&fluxinput[B1]); // fills B1->B3

#if(DOENTROPY!=DONOENTROPY)
  entropyflux_calc(pr, dir, q, &fluxinput[ENTROPY]); // fills ENTROPY only
  // below is special for utoprim() 5D version for full entropy evolution and inversion
  if(returntype==UENTROPY){
    fluxinput[UU]=fluxinput[ENTROPY]; // overwrite for utoprim()
    returntype=UNOTHING; // reset returntype for UtoU
  }
#endif

  

  // notice that geometry comes after subtractions/additions of EOMs
  UtoU(UNOTHING,returntype,geom,fluxinput,flux); // convert from UNOTHING->returntype

  return (0);
}


int massflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *massflux)
{
  /* particle number flux */
  *massflux = pr[RHO] * q->ucon[dir];

  return(0);
}

int entropyflux_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *entropyflux)
{
  FTYPE entropy;
  int entropy_calc(FTYPE *pr, FTYPE *entropy);

  // get entropy
  entropy_calc(pr,&entropy);
  /* entropy per unit rest-mass flux */
  // entropy=entropy per unit volume, where conserved quantity is specific entropy:
  // d/d\tau(entropy/rho)=0
  // -> \nabla_\mu(entropy u^\mu)=0

  *entropyflux = entropy * q->ucon[dir];

  return(0);
}


// ideal gas entropy
int entropy_calc(FTYPE *pr, FTYPE *entropy)
{
  FTYPE rho,ie,pressure;
  FTYPE indexn;

  rho=pr[RHO];
  ie=pr[UU];

  pressure=(gam-1.0)*ie;
  indexn=1.0/(gam-1.0);
  
  *entropy=rho*log(pow(pressure,indexn)/pow(rho,indexn+1.0));

  return(0);
}

int invertentropyflux_calc(FTYPE entropyflux,int dir, struct of_state *q, FTYPE *pr)
{
  FTYPE entropy;
  int ufromentropy_calc(FTYPE entropy, struct of_state *q, FTYPE *pr);

  // get entropy
  entropy=entropyflux/q->ucon[dir];
  ufromentropy_calc(entropy,q,pr);


  return(0);
}

// ideal gas u from ideal gas entropy (uses pr[RHO])
int ufromentropy_calc(FTYPE entropy, struct of_state *q, FTYPE *pr)
{
  FTYPE rho,ie,pressure;
  FTYPE indexn;

  rho=pr[RHO];
  indexn=1.0/(gam-1.0);
  
  // entropy version of ie
  pr[ENTROPY]=pow(pow(rho,indexn+1.0)*exp(entropy/rho),1.0/indexn)/(gam-1.0);

  return(0);
}


// spatial part of dualfaraday
int dualfaradayspatial_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *dualf)
{
  FTYPE dualffull[NDIM];
  int dualfullfaraday_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *dualf);


  dualfullfaraday_calc(pr,dir,q,dualffull);
  dualf[0]=dualffull[1];
  dualf[1]=dualffull[2];
  dualf[2]=dualffull[3];

  return(0);

}


#define GENMAXWELL 0
#define PRIMMAXWELL 1
//#define MAXWELL GENMAXWELL
#define MAXWELL PRIMMAXWELL

// returns \dF^{\mu dir}
// well, actually returns dualffull[dir], so gives columns instead of rows
int dualfullfaraday_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *dualffull)
{

#if(MAXWELL==GENMAXWELL)
  /* dual of Maxwell tensor */
  dualffull[0] = q->bcon[0] * q->ucon[dir] - q->bcon[dir] * q->ucon[0];
  dualffull[1] = q->bcon[1] * q->ucon[dir] - q->bcon[dir] * q->ucon[1];
  dualffull[2] = q->bcon[2] * q->ucon[dir] - q->bcon[dir] * q->ucon[2];
  dualffull[3] = q->bcon[3] * q->ucon[dir] - q->bcon[dir] * q->ucon[3];
#elif(MAXWELL==PRIMMAXWELL)
  if(dir>0){
    /* dual of Maxwell tensor */
    // dir refers to the direction of the derivative of the dualffull
    // B1,B2,B3 refers to LHS of equation dB^i/dt
    // due to antisymmetry, dir==i is 0
    dualffull[0] = - pr[B1+dir-1] ; // dualffull[i]=\dF^{i dir} where \dF^{0 dir} =-B^{dir}
    dualffull[1] = (pr[B1] * q->ucon[dir] - pr[B1+dir-1] * q->ucon[1])/q->ucon[0];
    dualffull[2] = (pr[B2] * q->ucon[dir] - pr[B1+dir-1] * q->ucon[2])/q->ucon[0];
    dualffull[3] = (pr[B3] * q->ucon[dir] - pr[B1+dir-1] * q->ucon[3])/q->ucon[0];
  }
  else{
    dualffull[0] = 0;
    dualffull[1] = pr[B1];
    dualffull[2] = pr[B2];
    dualffull[3] = pr[B3];
  }
#endif

  return(0);

}

// returns entire space-time(NDIM in size) / EOM(NPR in size) matrix
int primtofullflux(int returntype, FTYPE *pr, struct of_state *q,
	       struct of_geom *ptrgeom, FTYPE (*flux)[NPR])
{
  int j;
  
  // j=0,1,2,3 corresponding to U^j_\nu , where \nu corresponds to all EOMs and j to space-time for each
  DLOOPA primtoflux(returntype,pr,q,j,ptrgeom,flux[j]);

  return(0);
}


/* calculate "conserved" quantities */
int primtoU(int returntype, FTYPE *pr, struct of_state *q, struct of_geom *geom,
	    FTYPE *U)
{
  int i = 0, j = 0, k = 0;
  MYFUN(primtoflux(returntype,pr, q, 0, geom, U) ,"phys.c:primtoU()", "primtoflux_calc() dir=0", 1);

  return (0);
}



// standardized U form is geometry free and 
// \rho u^t , T^t_\nu , *F^{it}

// convert one form of U(or component of Flux) to another form
// UtoU controls meaning of todo and REMOVERESTMASSFROMUU.
// present order means start with geometry-free EOMs, add/subtract them, THEN geometry is assigned to that list of new EOMs.
// can choose to change order so that add geometry terms, THEN add/subtract them.  Rest of code shouldn't care (except source_conn()'s first connection)
void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout)
{
  FTYPE Ugeomfree[NPR];
  int k;

  /////////////////////
  //
  // input
  //
  if(inputtype==UEVOLVE){
    PLOOP Ugeomfree[k]=Uin[k]/ptrgeom->e[k];
    if((REMOVERESTMASSFROMUU==1)&&(EOMTYPE!=EOMFFDE)){
      // go back to standard stress-energy tensor form
      Ugeomfree[UU]  +=  - Ugeomfree[RHO] ; // - means adding back rest-mass
    }
  }
  else if(inputtype==UDIAG){
    PLOOP Ugeomfree[k]=Uin[k]/ptrgeom->g;
  }
  else if(inputtype==UNOTHING){
    PLOOP Ugeomfree[k]=Uin[k];
  }

  // at this point, Ugeomfree is geometry-free standard form of conserved quantities

  /////////////////////////
  //
  // output
  //

  if(returntype==UEVOLVE){
    if((REMOVERESTMASSFROMUU==1)&&(EOMTYPE!=EOMFFDE)){ // diagnostics want normal stress-energy tensor
      // "subtract" rest-mass
      Ugeomfree[UU] += Ugeomfree[RHO];
    }
    PLOOP Uout[k]=Ugeomfree[k]*ptrgeom->e[k];
  }
  else if(returntype==UDIAG){
    PLOOP Uout[k]=Ugeomfree[k]*ptrgeom->g;
  }
  else if(returntype==UNOTHING){
    PLOOP Uout[k]=Ugeomfree[k];
  }



}

// standardized primitive form is assumed to be
// \rho , u, \tilde{u}^\mu , *F^{it}=B^i
// where \tilde{u} is relative 4-velocity, as relative to $n_\mu = (-\alpha,0,0,0)$ and $\alpha^2=-1/g^{tt}$.

// For any space-time with no time-like curves this 4-velocity is always single-valued (i.e. unique).  It can also take on any value, so a reasonable primitive quantity.

// convert from one primitive form to another
void PtoP(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *pin, FTYPE *pout)
{
  FTYPE pstandard[NPR];
  int k;



}


/* calculate magnetic field four-vector */
void bcon_calc(FTYPE *pr, FTYPE *ucon, FTYPE *ucov, FTYPE *bcon)
{
  int j;

  bcon[TT] = pr[B1] * ucov[1] + pr[B2] * ucov[2] + pr[B3] * ucov[3];
  for (j = 1; j <= 3; j++)
    bcon[j] = (pr[B1 - 1 + j] + bcon[TT] * ucon[j]) / ucon[TT];

  return;
}

// inverse of bcon_calc()
void Bcon_calc(struct of_state *q, FTYPE*B)
{
  FTYPE uu0,uu1,uu2,uu3;
  FTYPE ud0,ud1,ud2,ud3;
  FTYPE bu1,bu2,bu3;
  FTYPE denom;

  uu0=q->ucon[TT];
  uu1=q->ucon[RR];
  uu2=q->ucon[TH];
  uu3=q->ucon[PH];

  ud0=q->ucov[TT];
  ud1=q->ucov[RR];
  ud2=q->ucov[TH];
  ud3=q->ucov[PH];

  bu1=q->bcon[RR];
  bu2=q->bcon[TH];
  bu3=q->bcon[PH];
 
  denom=1.0/(1.0+ud1*uu1+ud2*uu2+ud3*uu3);
  
  B[1]=uu0*(-(bu2*ud2+bu3*ud3)*uu1+bu1*(1.0+ud2*uu2+ud3*uu3))*denom;
  B[2]=uu0*(-(bu1*ud1+bu3*ud3)*uu2+bu2*(1.0+ud1*uu1+ud3*uu3))*denom;
  B[3]=uu0*(-(bu2*ud2+bu2*ud2)*uu3+bu3*(1.0+ud2*uu2+ud1*uu1))*denom;


}


// convert (e^\mu=0 case) b^\mu and (3-velocity in coordinate lab frame) v^\mu to pr
void vbtopr(FTYPE *vcon,FTYPE *bcon,struct of_geom *geom, FTYPE *pr)
{
  int k;
  struct of_state q;
  void Bcon_calc(struct of_state *q, FTYPE*B);
  FTYPE prim[NPR];
  FTYPE ucon[NDIM];

  // go ahead and get pr velocity
  PLOOP prim[k]=0.0;
  prim[U1]=vcon[1];
  prim[U2]=vcon[2];
  prim[U3]=vcon[3];

  //  vcon2pr(WHICHVEL,vcon,geom,pr); // need u^\mu, so do below instead
  ucon_calc_3vel(prim,geom,q.ucon);
  ucon2pr(WHICHVEL,q.ucon,geom,pr); // fills pr[U1->U3]
  
  //  q.ucon[TT]=ucon[TT];
  //  q.ucon[RR]=ucon[RR];
  //  q.ucon[TH]=ucon[TH];
  //  q.ucon[PH]=ucon[PH];
  
  lower(q.ucon,geom,q.ucov);

  //  q.bcon[TT]=bcon[TT]; // not used below
  q.bcon[RR]=bcon[RR];
  q.bcon[TH]=bcon[TH];
  q.bcon[PH]=bcon[PH];
  
  Bcon_calc(&q,&pr[B1-1]); // &pr[B1-1] since Bcon_calc() fills 1-3



}

/* MHD stress tensor, with first index up, second index down */
// mhd^dir_j
void mhd_calc(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd)
{
  void mhd_calc_0(FTYPE *pr, int dir, struct of_state *q, FTYPE *mhd);
  void mhd_calc_norestmass(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd);

  if(REMOVERESTMASSFROMUU==2) mhd_calc_norestmass(pr, dir, geom, q, mhd);
  else mhd_calc_0(pr, dir, q, mhd);



}


/* MHD stress tensor, with first index up, second index down */
void mhd_calc_0(FTYPE *pr, int dir, struct of_state *q, FTYPE *mhd)
{
  int j;
  FTYPE r, u, P, w, bsq, eta, ptot;

  // below allows other scalars to be advected but not affect the stress-energy equations of motion
#if(EOMTYPE==EOMGRMHD)
  r = pr[RHO];
  u = pr[UU];
  P = (gam - 1.) * u;
  w = P + r + u;
#elif(EOMTYPE==EOMFFDE)
  r=u=P=w=0.0;
#endif
  bsq = dot(q->bcon, q->bcov);
  eta = w + bsq;
  ptot = P + bsq*0.5;

  /* single row of mhd stress tensor, first index up, second index down 
   */
  // mhd^{dir}_{j} =
  // j=0..3
  DLOOPA mhd[j] = eta * q->ucon[dir] * q->ucov[j]
      + ptot * delta(dir, j) - q->bcon[dir] * q->bcov[j];

}


/* MHD stress tensor, with first index up, second index down */
void mhd_calc_norestmass(FTYPE *pr, int dir, struct of_geom *geom, struct of_state *q, FTYPE *mhd)
{
  int j;
  FTYPE rho, u, P, w, bsq, eta, ptot;
  FTYPE plus1ud0;
  void compute_1plusud0(FTYPE *pr,struct of_geom *geom, struct of_state *q, FTYPE *plus1ud0); // plus1ud0=(1+q->ucov[TT])

  // below allows other scalars to be advected but not affect the stress-energy equations of motion
#if(EOMTYPE==EOMGRMHD)
  rho = pr[RHO];
  u = pr[UU];
  P = (gam - 1.) * u;
  w = P + rho + u;
#elif(EOMTYPE==EOMFFDE)
  rho=u=P=w=0.0;
#endif
  bsq = dot(q->bcon, q->bcov);
  eta = w + bsq;
  ptot = P + bsq*0.5;

  compute_1plusud0(pr,geom,q,&plus1ud0); // plus1ud0=(1+q->ucov[TT])

  /* single row of mhd stress tensor, first index up, second index down 
   */
  // mhd^{dir}_{j} =
  // j=0..3
  j=0;
  // eta u^dir u_j + rho u^dir = (p+u+b^2) u^dir u_j + rho u^dir u_j + rho u^dir
  // = (p+u+b^2) u^dir u_j + rho u^dir (u_j + 1)
  mhd[j] = (P+u+bsq) * q->ucon[dir] *q->ucov[j] + rho * q->ucon[dir] * plus1ud0   + ptot * delta(dir, j) - q->bcon[dir] * q->bcov[j];
  SLOOPA mhd[j] = eta * q->ucon[dir] * q->ucov[j]   + ptot * delta(dir, j) - q->bcon[dir] * q->bcov[j];

}


// plus1ud0=(1+q->ucov[TT])
void compute_1plusud0(FTYPE *pr,struct of_geom *geom, struct of_state *q, FTYPE *plus1ud0)
{
  int j,k;
  volatile FTYPE plus1gv00;
  FTYPE AA,BB,alpha;
  FTYPE vcon[NDIM];

  // 3-velocity in coordinate basis
  SLOOPA vcon[j]=q->ucon[j]/q->ucon[TT];

  plus1gv00=1.0+geom->gcov[TT][TT];

  AA=0.0;
  SLOOPA AA+=2.0*geom->gcov[TT][j]*vcon[j];
  SLOOP AA+=geom->gcov[j][k]*vcon[j]*vcon[k];
  //AA/=geom->gcov[TT][TT];
  BB=geom->gcov[TT][TT];

  alpha=0.0;
  SLOOPA alpha+=geom->gcov[j][TT]*q->ucon[j];

  //  *plus1ud0=(plus1gv00+(2.0*alpha+alpha*alpha)*(1.0+AA)+AA)/((1.0+alpha)*(1.0+AA)+sqrt(-geom->gcov[TT][TT]*(1.0+AA)));

  *plus1ud0=(plus1gv00*BB+(2.0*alpha+alpha*alpha)*(BB+AA)+AA)/((1.0+alpha)*(BB+AA)+BB*sqrt(-(BB+AA)));

}


/* add in source terms to equations of motion */
int source(FTYPE *pr, struct of_geom *ptrgeom,
	   FTYPE (*dUcomp)[NPR], FTYPE *dU)
{
  //  double (*)[8]
  int i = 0, j = 0, k = 0,sc=0;
  struct of_state q;
  int source_conn(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q,FTYPE *dU);



  PLOOP{
    SCLOOP dUcomp[sc][k] = 0.;
    dU[k] = 0.;
  }



  MYFUN(get_state(pr, ptrgeom, &q) ,"phys.c:source()", "get_state() dir=0", 1);


  // physics source terms
  sourcephysics(pr, ptrgeom, &q,dUcomp);

  //////////////////
  //
  // Now deal with equation of motion factors
  SCLOOP PLOOP dUcomp[sc][k] *= ptrgeom->e[k];


  // geometry source (already does contain geometry prefactor term)
  source_conn(pr,ptrgeom, &q,dUcomp[GEOMSOURCE]);

  /////////
  //
  // now compute total since that's all the evolution cares about (comp just for diagnostics).
  SCLOOP PLOOP dU[k]+=dUcomp[sc][k];
  


  /* done! */
  return (0);
}


/* add in source terms related to geometry to equations of motion */
int source_conn(FTYPE *pr, struct of_geom *ptrgeom,
	   struct of_state *q,FTYPE *dU)
{
  int i = 0, j = 0, k = 0, l=0;
  FTYPE dUconn[NPR];
  FTYPE todo[NPR];
  FTYPE ftemp;
  FTYPE mhd[NDIM][NDIM];
  FTYPE flux[NDIM][NPR];
  int primtofullflux(int returntype, FTYPE *pr, struct of_state *q, struct of_geom *ptrgeom, FTYPE (*flux)[NPR]);
  void mhd_calc_0(FTYPE *pr, int dir, struct of_state *q, FTYPE *mhd);


  if((ANALYTICSOURCE)&&(defcoord==0)&&(MCOORD==KSCOORDS)){
    // have both WHICHEOM==WITHGDET and WHICHEOM==WITHNOGDET
    mks_source_conn(pr,ptrgeom, q, dU); // returns without geometry prefactor
    PLOOP dU[k]*=ptrgeom->e[k];
  }
  else { // general non-analytic form, which uses an analytic or numerical origin for conn/conn2 (see set_grid.c)


    /////////////////////
    //
    // define first connection
    //
    // notice mhd_calc(), rather than primtoflux(), is used because this is a special connection only operated on the (at most) 4 energy-momentum EOMs.
    //
    // notice that mhd_calc_0 is used, which includes rest-mass since connection coefficient source term has rest-mass.  The rest-mass was only subtracted out of conserved quantity and flux term.
    DLOOPA  mhd_calc_0(pr, j, q, mhd[j]);

    /* contract mhd stress tensor with connection */
    // mhd^{dir}_{comp} = mhd^j_k
    // dU^{l} = T^j_k C^{k}_{l j}
    PLOOP dUconn[k]=0;
    for(l=0;l<NDIM;l++)  DLOOP{
      dUconn[UU+l] += mhd[j][k] * conn[ptrgeom->i][ptrgeom->j][k][l][j];
    }


    ////////////////////////////
    //
    // use first connection
    //
    // true as long as UU->U3 are not added/subtracted from other EOMs
    // Note that UtoU is ignorant of this connection term and additions/subtractions of various EOMs to/from another must account for this connection here.
    PLOOP dU[k]+=ptrgeom->e[k]*dUconn[k]; 

 

    ///////////////////
    //
    // second connection
    //

#if(WHICHEOM!=WITHGDET)
    // deal with second connection.  Now can use general flux as defined by primtofullflux since all EOMs operate similarly with respect to second connection
    primtofullflux(UEVOLVE,pr,q,ptrgeom,flux); // returns with geometry prefactor


    // todo = whether that EOM has the NOGDET form.  If so, then need 2nd connection.  Do this instead of different connection2 for each EOM since each spatial component is actually the same.

    todo[RHO]=(NOGDETRHO>0) ? 1 : 0;
    todo[UU]=(NOGDETU0>0) ? 1 : 0;
    todo[U1]=(NOGDETU1>0) ? 1 : 0;
    todo[U2]=(NOGDETU2>0) ? 1 : 0;
    todo[U3]=(NOGDETU3>0) ? 1 : 0;
    todo[B1]=(NOGDETB1>0) ? 1 : 0;
    todo[B2]=(NOGDETB2>0) ? 1 : 0;
    todo[B3]=(NOGDETB3>0) ? 1 : 0;
    if(DOENTROPY) todo[ENTROPY]=(NOGDETENTROPY>0) ? 1 : 0;

    // conn2 is assumed to take care of sign
    // conn2 has geom->e and normal d(ln(gdet)) factors combined

    if(REMOVERESTMASSFROMUU){
      if(todo[RHO]!=todo[UU]){
	dualfprintf(fail_file,"You are stupid\n");
	myexit(1);
      }
    }

    //////////////////////////////////
    //
    // notice that we assume equations are differenced first, then one manipulates the connections.
    // Thus, the manipulation of connections applies to final form of EOMs afer differencing or adding.
    // This agrees with how code evolves conserved quantities, so simplest choice to make.
    //
    // Alternatively stated, UtoU() controls this ordering issue.

    PLOOP DLOOPA dU[k] += todo[k]*(flux[j][k] * conn2[ptrgeom->i][ptrgeom->j][j]);

#endif

  }



  /* done! */
  return (0);
}

/* returns b^2 (i.e., twice magnetic pressure) */
int bsq_calc(FTYPE *pr, struct of_geom *ptrgeom, FTYPE *bsq)
{
  int i = 0, j = 0, k = 0;
  struct of_state q;

  MYFUN(get_state(pr, ptrgeom, &q) ,"phys.c:bsq_calc()", "get_state() dir=0", 1);
  *bsq = dot(q.bcon, q.bcov);
  return (0);
}


void lower(FTYPE *ucon, struct of_geom *geom, FTYPE *ucov)
{
        ucov[0] = geom->gcov[0][0]*ucon[0]
                + geom->gcov[0][1]*ucon[1]
                + geom->gcov[0][2]*ucon[2]
	  + geom->gcov[0][3]*ucon[3] ;
        ucov[1] = geom->gcov[0][1]*ucon[0]
                + geom->gcov[1][1]*ucon[1]
                + geom->gcov[1][2]*ucon[2]
	  + geom->gcov[1][3]*ucon[3] ;
        ucov[2] = geom->gcov[0][2]*ucon[0]
                + geom->gcov[1][2]*ucon[1]
                + geom->gcov[2][2]*ucon[2]
	  + geom->gcov[2][3]*ucon[3] ;
        ucov[3] = geom->gcov[0][3]*ucon[0]
                + geom->gcov[1][3]*ucon[1]
                + geom->gcov[2][3]*ucon[2]
	  + geom->gcov[3][3]*ucon[3] ;

        return ;
}

void lowerf(FTYPE *fcon, struct of_geom *geom, FTYPE *fcov)
{
  int j,k;
  int jp,kp;
  FTYPE myfcon[NDIM][NDIM],myfcov[NDIM][NDIM];

  myfcon[0][0]=myfcon[1][1]=myfcon[2][2]=myfcon[3][3]=0;
  myfcon[0][1]=fcon[0];
  myfcon[0][2]=fcon[1];
  myfcon[0][3]=fcon[2];
  myfcon[1][2]=fcon[3];
  myfcon[1][3]=fcon[4];
  myfcon[2][3]=fcon[5];
  //
  myfcon[1][0]=-fcon[0];
  myfcon[2][0]=-fcon[1];
  myfcon[3][0]=-fcon[2];
  myfcon[2][1]=-fcon[3];
  myfcon[3][1]=-fcon[4];
  myfcon[3][2]=-fcon[5];
  
  DLOOP{
    myfcov[j][k]=0;
    for(jp=0;jp<NDIM;jp++) for(kp=0;kp<NDIM;kp++){
      myfcov[j][k]+=myfcon[jp][kp]*geom->gcov[j][jp]*geom->gcov[k][kp];
    }
  }
  fcov[0]=myfcov[0][1];
  fcov[1]=myfcov[0][2];
  fcov[2]=myfcov[0][3];
  fcov[3]=myfcov[1][2];
  fcov[4]=myfcov[1][3];
  fcov[5]=myfcov[2][3];

  return ;
}

void raise(FTYPE *ucov, struct of_geom *geom, FTYPE *ucon)
{

        ucon[0] = geom->gcon[0][0]*ucov[0]
                + geom->gcon[0][1]*ucov[1]
                + geom->gcon[0][2]*ucov[2]
	  + geom->gcon[0][3]*ucov[3] ;
        ucon[1] = geom->gcon[0][1]*ucov[0]
                + geom->gcon[1][1]*ucov[1]
                + geom->gcon[1][2]*ucov[2]
	  + geom->gcon[1][3]*ucov[3] ;
        ucon[2] = geom->gcon[0][2]*ucov[0]
                + geom->gcon[1][2]*ucov[1]
                + geom->gcon[2][2]*ucov[2]
	  + geom->gcon[2][3]*ucov[3] ;
        ucon[3] = geom->gcon[0][3]*ucov[0]
                + geom->gcon[1][3]*ucov[1]
                + geom->gcon[2][3]*ucov[2]
	  + geom->gcon[3][3]*ucov[3] ;

        return ;
}

/* find ucon, ucov, bcon, bcov from primitive variables */
int get_state(FTYPE *pr, struct of_geom *ptrgeom, struct of_state *q)
{
  int i = 0, j = 0, k = 0;

  /* get ucon */
  
  MYFUN(ucon_calc(pr, ptrgeom, q->ucon) ,"phys.c:get_state()", "ucon_calc()", 1);
  lower(q->ucon, ptrgeom, q->ucov);
  bcon_calc(pr, q->ucon, q->ucov, q->bcon);
  lower(q->bcon, ptrgeom, q->bcov);
  return (0);
}

/* load local geometry into structure geom */
void get_geometry(int ii, int jj, int pp, struct of_geom *geom)
{
  int j, k;
  void assign_eomfunc(struct of_geom *geom, FTYPE eomfunc);
  

  //  DLOOP geom->gcov[j][k] = gcov[ii][jj][pp][j][k];
  //DLOOP geom->gcon[j][k] = gcon[ii][jj][pp][j][k];
  // let's vectorize it
  /*
  for(j=0;j<=NDIM*NDIM-1;j++){
    geom->gcon[0][j] = gcon[ii][jj][pp][0][j];
    geom->gcov[0][j] = gcov[ii][jj][pp][0][j];
  }
  
  geom->g = gdet[ii][jj][pp];
  */
  // let's pointer it, faster by a bit
  geom->gcov=(FTYPE (*)[NDIM])(&(gcov[ii][jj][pp][0][0]));
  geom->gcon=(FTYPE (*)[NDIM])(&(gcon[ii][jj][pp][0][0]));

  geom->g = gdet[ii][jj][pp];

  geom->i = ii;
  geom->j = jj;
  geom->p = pp;
#if(JONCHECKS)
  icurr = ii;
  jcurr = jj;
  pcurr = pp;
#endif

  // get eomfunc (see eomfunc_func() and lngdet_func_mcoord() in metric.c)
  assign_eomfunc(geom,eomfunc[ii][jj][pp]);

}


//eomfuncarray[ii][jj][pp];
// stationary factor in equation of motion that multiplies T^\mu_\nu
void assign_eomfunc(struct of_geom *geom, FTYPE eomfunc)
{

  // now set each EOM
#if(NOGDETRHO==0)
  geom->e[RHO]=geom->g;
#else
  geom->e[RHO]=eomfunc;
#endif
#if(NOGDETU0==0)
  geom->e[UU]=geom->g;
#else
  geom->e[UU]=eomfunc;
#endif
#if(NOGDETU1==0)
  geom->e[U1]=geom->g;
#else
  geom->e[U1]=eomfunc;
#endif
#if(NOGDETU2==0)
  geom->e[U2]=geom->g;
#else
  geom->e[U2]=eomfunc;
#endif
#if(NOGDETU3==0)
  geom->e[U3]=geom->g;
#else
  geom->e[U3]=eomfunc;
#endif

#if(NOGDETB1==0)
  geom->e[B1]=geom->g;
#else
  geom->e[B1]=eomfunc;
#endif
#if(NOGDETB2==0)
  geom->e[B2]=geom->g;
#else
  geom->e[B2]=eomfunc;
#endif
#if(NOGDETB3==0)
  geom->e[B3]=geom->g;
#else
  geom->e[B3]=eomfunc;
#endif

#if(DOENTROPY)
#if(NOGDETENTROPY==0)
  geom->e[ENTROPY]=geom->g;
#else
  geom->e[ENTROPY]=eomfunc;
#endif
#endif


}



/* find contravariant four-velocity from the relative 4 velocity */
int ucon_calc_rel4vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon)
{
        FTYPE alpha,gamma ;
        FTYPE beta[NDIM] ;
        int j ;

        alpha = 1./sqrt(-geom->gcon[TT][TT]) ;
        SLOOPA beta[j] = geom->gcon[TT][j]*alpha*alpha ;

        MYFUN(gamma_calc(pr,geom,&gamma),"ucon_calc_rel4vel: gamma calc failed\n","phys.c",1);

        ucon[TT] = gamma/alpha ;
        SLOOPA ucon[j] = pr[U1+j-1] - gamma*beta[j]/alpha ;

        return(0) ;
}

/* find gamma-factor wrt normal observer */
int gamma_calc(FTYPE *pr, struct of_geom *geom,FTYPE*gamma)
{
        FTYPE qsq ;
	int j,k;

        qsq =     geom->gcov[1][1]*pr[U1]*pr[U1]
                + geom->gcov[2][2]*pr[U2]*pr[U2]
                + geom->gcov[3][3]*pr[U3]*pr[U3]
            + 2.*(geom->gcov[1][2]*pr[U1]*pr[U2]
                + geom->gcov[1][3]*pr[U1]*pr[U3]
                + geom->gcov[2][3]*pr[U2]*pr[U3]) ;

#if(JONCHECKS2&&0)
        if(qsq<0.0){
	  if(qsq>-NUMEPSILON*100){ // then assume not just machine precision
	    qsq=0.0;
	  }
	  else{
	    dualfprintf(fail_file,"gamma_calc failed: i=%d j=%d qsq=%21.15g\n",geom->i,geom->j,qsq);
	    PLOOP dualfprintf(fail_file,"pr[%d]=%21.15g\n",k,pr[k]);
	    DLOOP dualfprintf(fail_file,"gcov[%d][%d]=%21.15g\n",j,k,geom->gcov[j][k]);
	    if (fail(FAIL_UTCALC_DISCR) >= 1)
	      return (1);
	  }
	}
#endif
        *gamma = sqrt(1. + qsq) ;

        return(0) ;
}





/* find contravariant four-velocity */
//int ucon_calc(FTYPE *pr, struct of_geom *geom, FTYPE *ucon)
int ucon_calc_3vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon)
{
  FTYPE discr;
  // debug stuff
  int j,k;
  FTYPE r,th,X[NDIM];
  FTYPE vcon[NDIM];

  SLOOPA vcon[j] = pr[U1+j-1];

  discr =
    + geom->gcov[0][0]
    + geom->gcov[1][1] * vcon[1] * vcon[1]
    + geom->gcov[2][2] * vcon[2] * vcon[2]
    + geom->gcov[3][3] * vcon[3] * vcon[3]
    + 2. * (geom->gcov[0][1]* vcon[1]
	    + geom->gcov[0][2] * vcon[2]
	    + geom->gcov[0][3] * vcon[3]
	    + geom->gcov[1][2] * vcon[1] * vcon[2]
	    + geom->gcov[1][3] * vcon[1] * vcon[3]
	    + geom->gcov[2][3] * vcon[2] * vcon[3]
	    );

#if(JONCHECKS)  
  uttdiscr=-discr;
#endif

  if (discr > 0.) {
#if(JONCHECKS&&0)
    if(whocalleducon==0){
      // then report on disc
      dualfprintf(fail_file,"disc=%21.15g, should be negative\n",discr);
      for(k=U1;k<=U3;k++){
	dualfprintf(fail_file,"uconfailed on pr[%d]=%21.15g\n",k,pr[k]);
      }
      coord(geom->i,geom->j,geom->p,X);
      bl_coord(X,&r,&th);
      dualfprintf(fail_file,"i=%d j=%d pcurr=%d\nx1=%21.15g x2=%21.15g \nr=%21.15g th=%21.15g \ng=%21.15g\n",startpos[1]+geom->i,startpos[2]+geom->j,geom->p,X[1],X[2],r,th,geom->g);
      dualfprintf(fail_file,"\ngcon\n");
      dualfprintf(fail_file,"{");
      for(j=0;j<NDIM;j++){
	dualfprintf(fail_file,"{");
	for(k=0;k<NDIM;k++){
	  dualfprintf(fail_file,"%21.15g",geom->gcon[j][k]);
	  if(k!=NDIM-1) dualfprintf(fail_file," , ");
	}
	dualfprintf(fail_file,"}");	
	if(j!=NDIM-1) dualfprintf(fail_file," , ");
      }
      dualfprintf(fail_file,"}");
      dualfprintf(fail_file,"\ngcov\n");
      dualfprintf(fail_file,"{");
      for(j=0;j<NDIM;j++){
	dualfprintf(fail_file,"{");
	for(k=0;k<NDIM;k++){
	  dualfprintf(fail_file,"%21.15g",geom->gcov[j][k]);
	  if(k!=NDIM-1) dualfprintf(fail_file," , ");
	}
	dualfprintf(fail_file,"}");	
	if(j!=NDIM-1) dualfprintf(fail_file," , ");
      }
      dualfprintf(fail_file,"}");
    }
#endif
    //    if (fail(FAIL_UTCALC_DISCR) >= 1)
      return (1);
  }

  ucon[TT] = 1. / sqrt(-discr);
  SLOOPA ucon[j]=vcon[j]*ucon[TT];

  return (0);
}



/* find contravariant time component of four-velocity from the 4velocity (3 terms)*/
int ucon_calc_4vel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon)
{
	FTYPE AA,BB,CC ;
	FTYPE discr ;
	FTYPE bsq,X[NDIM] ;
	int i=0,j=0,k=0 ;

        ucon[1] = pr[U1] ;
        ucon[2] = pr[U2] ;
        ucon[3] = pr[U3] ;

	AA = geom->gcov[TT][TT] ;
	BB = 2.*(geom->gcov[TT][1]*ucon[1] +
		 geom->gcov[TT][2]*ucon[2] +
		 geom->gcov[TT][3]*ucon[3]) ;
	CC = 1. +
	  geom->gcov[1][1]*ucon[1]*ucon[1] +
	  geom->gcov[2][2]*ucon[2]*ucon[2] +
	  geom->gcov[3][3]*ucon[3]*ucon[3] +
	  2.*(geom->gcov[1][2]*ucon[1]*ucon[2] +
	      geom->gcov[1][3]*ucon[1]*ucon[3] +
	      geom->gcov[2][3]*ucon[2]*ucon[3]) ;
	
	discr = BB*BB - 4.*AA*CC ;
	if(discr < 0.) {
		/*
		fprintf(fail_file,"failure %d %d\n",icurr,jcurr) ;
        	ucon[TT] = (-BB - sqrt(-discr))/(2.*AA) ;
        	ucon[TT] = -BB/(2.*AA) ;
		*/
#if(0)
		fprintf(fail_file,"failure: spacelike four-velocity %21.15g\n",
			discr) ;
		fprintf(fail_file,"i=%d j=%d %d\n",startpos[1]+geom->i,startpos[2]+geom->j,geom->p) ;
		coord(geom->i,geom->j,geom->p,X);
		fprintf(fail_file,"%21.15g %21.15g\n",X[1],X[2]) ;

		/*
		  if(bsq_calc(&bsq,pr[icurr][jcurr])>=1) FAILSTATEMENT("phys.c:ucon_calc()","bsq_calc()",1);

		fprintf(fail_file,"bsq/rho: %21.15g\n",bsq/pr[icurr][jcurr][0]) ;

		*/
		for(k=0;k<NPR;k++) fprintf(fail_file,"%d %21.15g\n",k,pr[k]) ;
		// GODMARK -- why did we have failed=1?
		//		failed=1;
#endif
		return(1);
	}

	ucon[TT] = (-BB - sqrt(discr))/(2.*AA) ;
	return(0) ;
}

/* find contravariant time component of four-velocity from the 4velocity (3 terms)*/
int ucon_calc_nonrel(FTYPE *pr, struct of_geom *geom, FTYPE *ucon)
{

  // this isn't really right
  // neglects kinetic energy term.  Need to really re-write T^\mu_\nu

        ucon[1] = pr[U1] ;
        ucon[2] = pr[U2] ;
        ucon[3] = pr[U3] ;
	ucon[TT] = 1.0;

	return(0) ;
}


FTYPE taper_func(FTYPE R,FTYPE rin)
{

  if(R <= rin)
    return(0.) ;
  else
    return(1. - sqrt(rin/R)) ;

}

// compute the radius of the inner most stable circular orbit
FTYPE rmso_calc(int which)
{
  FTYPE rmso,Z1,Z2,sign ;


  if(which==PROGRADERISCO) sign=1; else sign=-1;

  Z1 = 1. + pow(1. - a*a,1./3.)*(pow(1. + a,1./3.) +
                                 pow(1. - a, 1./3.)) ;
  Z2 = sqrt(3.*a*a + Z1*Z1) ;
  rmso=3. + Z2-sign*sqrt((3. - Z1)*(3. + Z1 + 2.*Z2)) ;

  return(rmso) ;
}

FTYPE uphi_isco_calc(int which,FTYPE r)
{
  FTYPE uphi;
  FTYPE sign;
  FTYPE Z1,Z2;

  if(which==PROGRADERISCO) sign=1; else sign=-1;

  Z1=r*r-sign*2.*a*sqrt(r)+a*a;
  Z2=r*(r*r-3.*r+sign*2.*a*sqrt(r));

  uphi=sign*Z1/sqrt(Z2);

  return(uphi);

}

FTYPE rhor_calc(int which)
{
  FTYPE sign;
  if(which==0) sign=1; else sign=-1;

  return(1. +sign*sqrt(1. - a * a));
}




// used Mathematica's MinimumChangePermutations and Signature
// updown = 0 : down
// updown = 1 : up
FTYPE lc4(int updown, FTYPE detg, int mu,int nu,int kappa,int lambda)
{
  int i;
  FTYPE lc4sign; // 1,-1,1,-1... for all 24 entires
  int l1[24]={1, 2, 3, 1, 2, 3, 4, 2, 1, 4, 2, 1, 1, 3, 4, 1, 3, 4, 4, 3, 2, 4, 3, 2};
  int l2[24]={2, 1, 1, 3, 3, 2, 2, 4, 4, 1, 1, 2, 3, 1, 1, 4, 4, 3, 3, 4, 4, 2, 2, 3};
  int l3[24]={3, 3, 2, 2, 1, 1, 1, 1, 2, 2, 4, 4, 4, 4, 3, 3, 1, 1, 2, 2, 3, 3, 4, 4};
  int l4[24]={4, 4, 4, 4, 4, 4, 3, 3, 3, 3, 3, 3, 2, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1};

  for(i=0;i<24;i++){
    if((1+mu==l1[i])&&(1+nu==l2[i])&&(1+kappa==l3[i])&&(1+lambda==l4[i])){
      lc4sign=(i%2) ? -1 : 1;
      if(updown==1) return(-1.0/detg*lc4sign); // up
      else if(updown==0) return(detg*lc4sign); // down
    }
  }
  // if didn't get here, then 0
  return(0.0);
}

// below not used currently
void faraday_calc(int which, FTYPE *b, FTYPE *u, struct of_geom *geom, FTYPE faraday[][NDIM])
{
  int nu,mu,kappa,lambda;

  for(nu=0;nu<NDIM;nu++) for(mu=0;mu<NDIM;mu++){
    faraday[mu][nu]=0.0;
    for(kappa=0;kappa<NDIM;kappa++) for(lambda=0;lambda<NDIM;lambda++){
      faraday[mu][nu]+=lc4(which,geom->g,mu,nu,kappa,lambda)*u[kappa]*b[lambda];
    }
  }

}


int current_doprecalc(int which, FTYPE p[][N2M][NPR])
{
  int i,j;
  int idel, jdel;
  struct of_geom geom;
  struct of_state q;
  FTYPE Dt;
  int face;


  if(WHICHCURRENTCALC==0){ // then need to save 2 times
    // which==1,2 should be using face primitives, but probably not
    if(which==0){ face=CENT; idel=0; jdel=0; Dt=dt*0.5; }
    else if(which==1){ face=FACE1; idel=1; jdel=0; Dt=dt; }
    else if(which==2){ face=FACE2; idel=0; jdel=1; Dt=dt; }
    else if(which==3){ face=CENT; idel=0; jdel=0; Dt=dt; }
  }
  else if(WHICHCURRENTCALC==1){
    face=CENT;
    idel=0;
    jdel=0;
    Dt=dt;
  }
  else if(WHICHCURRENTCALC==2){ // then need to save 2 times
    if(which==0){ face=CENT; idel=0; jdel=0; Dt=dt*0.5; }
    else if(which==1){ face=CENT; idel=1; jdel=0; Dt=dt; }
    else if(which==2){ face=CENT; idel=0; jdel=1; Dt=dt; }
    else if(which==3){ face=CENT; idel=0; jdel=0; Dt=dt; }
  }

  // assume no other conditions GODMARK

  //  FZLOOP(-jdel, -idel) {
  FULLLOOP{
    get_geometry(i, j, face, &geom);
    MYFUN(get_state(p[i][j], &geom, &q),"phys.c:current_doprecalc()", "get_state()", 1);
    current_precalc(which,&geom,&q,Dt,cfaraday[i][j]);
  }

  return(0);
}

// which: 0: somewhere at half step
// which: 1: doing flux calculation in x1 direction, full step
// which: 2: doing flux calculation in x2 direction, full step
// faraday[0-3][0-3]=[0=mid-time,1=radial edge final time, 2= theta edge final time, 3=old mid-time][whatever 3 things needed to compute relevant current]
void current_precalc(int which, struct of_geom *geom, struct of_state *q, FTYPE Dt,FTYPE faraday[][3])
{
  // assume outside loop is like flux, from 0..N in r for r, 0..N in h for h.  And must wait till 2nd timestep before computing the current since need time differences
  if(which==0){
    // assume got here when DT==dt/2 and geom and state set at zone center
    // first save old calculation
    if((WHICHCURRENTCALC==0)||(WHICHCURRENTCALC==2)){ // then need to save 2 times
      faraday[3][0]=faraday[0][0];
      faraday[3][1]=faraday[0][1];
      faraday[3][2]=faraday[0][2];
    }
    else if(WHICHCURRENTCALC==1){ // then need to save 3 times and have [3] as 2 times ago
      // 2 times ago
      faraday[3][0]=faraday[4][0];
      faraday[3][1]=faraday[4][1];
      faraday[3][2]=faraday[4][2];

      // 1 time ago
      faraday[4][0]=faraday[0][0];
      faraday[4][1]=faraday[0][1];
      faraday[4][2]=faraday[0][2];
    }
    // now calculate new version
    faraday[0][0]=-1.0/geom->g * (-q->bcov[3]*q->ucov[2]+q->bcov[2]*q->ucov[3]); // f^{rt}
    faraday[0][1]=-1.0/geom->g * (q->bcov[3]*q->ucov[1]-q->bcov[1]*q->ucov[3]); // f^{ht}
    faraday[0][2]=-1.0/geom->g * (-q->bcov[2]*q->ucov[1]+q->bcov[1]*q->ucov[2]); // f^{pt}
  }
  else if(which==1){
    // assume got here with DT=dt and geom and state at radial zone edge
    faraday[1][0]=-1.0/geom->g * (q->bcov[3]*q->ucov[2]-q->bcov[2]*q->ucov[3]); // f^{tr}=-f^{rt}
    faraday[1][1]=-1.0/geom->g * (-q->bcov[3]*q->ucov[0]+q->bcov[0]*q->ucov[3]); // f^{hr}
    faraday[1][2]=-1.0/geom->g * (q->bcov[2]*q->ucov[0]-q->bcov[0]*q->ucov[2]); // f^{pr}
  }
  else if(which==2){
    // assume got here with DT=dt and geom and state at theta zone edge
    faraday[2][0]=-1.0/geom->g * (-q->bcov[3]*q->ucov[1]+q->bcov[1]*q->ucov[3]); // f^{th}=-f^{ht}
    faraday[2][1]=-1.0/geom->g * (q->bcov[3]*q->ucov[0]-q->bcov[0]*q->ucov[3]); // f^{rh}=-f^{hr}
    faraday[2][2]=-1.0/geom->g * (-q->bcov[1]*q->ucov[0]+q->bcov[0]*q->ucov[1]); // f^{ph}
  }
  else if(which==3){
    // DT==dt, but zone center
    fcon[geom->i][geom->j][0]=-1.0/geom->g * (q->bcov[3]*q->ucov[2]-q->bcov[2]*q->ucov[3]); // f^{tr}
    fcon[geom->i][geom->j][1]=-1.0/geom->g * (-q->bcov[3]*q->ucov[1]+q->bcov[1]*q->ucov[3]); // f^{th}
    fcon[geom->i][geom->j][2]=-1.0/geom->g * (q->bcov[2]*q->ucov[1]-q->bcov[1]*q->ucov[2]); // f^{tp}
    fcon[geom->i][geom->j][3]=-1.0/geom->g * (q->bcov[3]*q->ucov[0]-q->bcov[0]*q->ucov[3]); // f^{rh}
    fcon[geom->i][geom->j][4]=-1.0/geom->g * (-q->bcov[2]*q->ucov[0]+q->bcov[0]*q->ucov[2]); // f^{rp}
    fcon[geom->i][geom->j][5]=-1.0/geom->g * (q->bcov[1]*q->ucov[0]-q->bcov[0]*q->ucov[1]); // f^{hp}
  }
}


void current_calc(FTYPE cfaraday[][N2M][NUMCURRENTSLOTS][3])
{
  void current_calc_0(FTYPE cfaraday[][N2M][NUMCURRENTSLOTS][3]);
  void current_calc_1(FTYPE cfaraday[][N2M][NUMCURRENTSLOTS][3]);
  void current_calc_2(FTYPE cfaraday[][N2M][NUMCURRENTSLOTS][3]);

  if(WHICHCURRENTCALC==0){
    current_calc_0(cfaraday);
  }
  else if(WHICHCURRENTCALC==1){
    current_calc_1(cfaraday);
  }
  else if(WHICHCURRENTCALC==2){
    current_calc_2(cfaraday);
  }

}

// the current is calculated to end up at the zone and time edge
// point is to have J^t at same time as rest of J's, although different spacial points
void current_calc_0(FTYPE cfaraday[][N2M][NUMCURRENTSLOTS][3])
{
  int i,j;
  struct of_geom geomt;
  struct of_geom geomr;
  struct of_geom geomh;
  struct of_geom geomrp1;
  struct of_geom geomhp1;
  static FTYPE lastdt;
  static int calls=0;
  FTYPE idtc,idx1,idx2;

  if(calls>0){ // since need 2 times

    idtc=2.0/(lastdt+dt);
    idx1=1.0/dx[1];
    idx2=1.0/dx[2];

    LOOPP11 LOOPP12{ // largest possible loop for this differencing (could isolate directions)
      get_geometry(i,j,CENT,&geomt);
      get_geometry(i,j,FACE1,&geomr);
      get_geometry(i,j,FACE2,&geomh);
      // geomtp1 is same as geomt since d/dt( geometry) -> 0
      get_geometry(i+1,j,FACE1,&geomrp1);
      get_geometry(i,j+1,FACE2,&geomhp1);
      
      // J^t = F^{tr},r + F^{th},h
      jcon[i][j][0]=
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][0]-geomr.g*cfaraday[i][j][1][0])*idx1+ // F^{tr},r
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][0]-geomh.g*cfaraday[i][j][2][0])*idx2; // F^{th},h
      
      // J^r = F^{rt},t + F^{rh},h
      jcon[i][j][1]=
	(cfaraday[i][j][0][0]-cfaraday[i][j][3][0])*idtc+ // F^{rt},t
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][1]-geomh.g*cfaraday[i][j][2][1])*idx2; // F^{rh},h
      
      // J^h = F^{ht},t + F^{hr},r
      jcon[i][j][2]=
	(cfaraday[i][j][0][1]-cfaraday[i][j][3][1])*idtc+ // F^{ht},t
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][1]-geomr.g*cfaraday[i][j][1][1])*idx1; // F^{hr},r
      
      // J^p = F^{pt},t + F^{pr},r + F^{ph},h
      jcon[i][j][3]=
	(cfaraday[i][j][0][2]-cfaraday[i][j][3][2])*idtc+ // F^{pt},t
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][2]-geomr.g*cfaraday[i][j][1][2])*idx1+ // F^{pr},r
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][2]-geomh.g*cfaraday[i][j][2][2])*idx2; // F^{ph},h
    }
  }
  else{

    idtc=2.0/(lastdt+dt);
    idx1=1.0/dx[1];
    idx2=1.0/dx[2];

    LOOPP11 LOOPP12{ // largest possible loop for this differencing (could isolate directions)
      get_geometry(i,j,CENT,&geomt);
      get_geometry(i,j,FACE1,&geomr);
      get_geometry(i,j,FACE2,&geomh);
      // geomtp1 is same as geomt since d/dt( geometry) -> 0
      get_geometry(i+1,j,FACE1,&geomrp1);
      get_geometry(i,j+1,FACE2,&geomhp1);
      
      // J^t = F^{tr},r + F^{th},h
      jcon[i][j][0]=
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][0]-geomr.g*cfaraday[i][j][1][0])*idx1+ // F^{tr},r
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][0]-geomh.g*cfaraday[i][j][2][0])*idx2; // F^{th},h
      
      // J^r = F^{rt},t + F^{rh},h
      jcon[i][j][1]=
	0+
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][1]-geomh.g*cfaraday[i][j][2][1])*idx2; // F^{rh},h
      
      // J^h = F^{ht},t + F^{hr},r
      jcon[i][j][2]=
	0+
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][1]-geomr.g*cfaraday[i][j][1][1])*idx1; // F^{hr},r
      
      // J^p = F^{pt},t + F^{pr},r + F^{ph},h
      jcon[i][j][3]=
	0+
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][2]-geomr.g*cfaraday[i][j][1][2])*idx1+ // F^{pr},r
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][2]-geomh.g*cfaraday[i][j][2][2])*idx2; // F^{ph},h
    }
  }
  calls++;
  lastdt=dt;

}



// the current is calculated to end up at the zone and time edge
// point is to have J^t at same time as rest of J's and at same spatial location
// the temporal value of things is obtained at same time of everything else by fitting to parabola in time
void current_calc_1(FTYPE cfaraday[][N2M][NUMCURRENTSLOTS][3])
{
  int i,j;
  struct of_geom geomt;
  struct of_geom geomr;
  struct of_geom geomh;
  struct of_geom geomrp1;
  struct of_geom geomhp1;
  struct of_geom geomrm1;
  struct of_geom geomhm1;
  static FTYPE lastdt;
  static long long calls=0;
  FTYPE idtc,idx1,idx2;
  FTYPE AA,BB;
  FTYPE dtl,dtr;
  FTYPE fl,fr,f0,derf;

  if(calls>=2){ // since need 3 times to properly time center without having to worry about what RK is doing

    dtl=lastdt;
    dtr=dt;
    idx1=1.0/(2.0*dx[1]);
    idx2=1.0/(2.0*dx[2]);

    LOOPP11 LOOPP12{ // largest possible loop for this differencing (could isolate directions)
      get_geometry(i,j,CENT,&geomt);
      get_geometry(i,j,CENT,&geomr);
      get_geometry(i,j,CENT,&geomh);
      // geomtp1 is same as geomt since d/dt( geometry) -> 0
      get_geometry(i+1,j,CENT,&geomrp1);
      get_geometry(i,j+1,CENT,&geomhp1);
      get_geometry(i-1,j,CENT,&geomrm1);
      get_geometry(i,j-1,CENT,&geomhm1);
      
      // J^t = F^{tr},r + F^{th},h
      jcon[i][j][0]=
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][0]-geomrm1.g*cfaraday[i-1][j][1][0])*idx1+ // F^{tr},r
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][0]-geomhm1.g*cfaraday[i][j-1][2][0])*idx2; // F^{th},h
      
      // J^r = F^{rt},t + F^{rh},h
      fl=cfaraday[i][j][3][0];
      f0=cfaraday[i][j][4][0];
      fr=cfaraday[i][j][0][0];
      AA=(dtr*(fl-f0)+dtl*(fr-f0))/(dtl*dtr*(dtl+dtr));
      BB=(-dtr*dtr*(fl-f0)+dtl*dtl*(fr-f0))/(dtl*dtr*(dtl+dtr));
      derf=2.0*AA*dtr+BB;

      jcon[i][j][1]=
	derf + // F^{rt},t
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][1]-geomhm1.g*cfaraday[i][j-1][2][1])*idx2; // F^{rh},h

      fl=cfaraday[i][j][3][1];
      f0=cfaraday[i][j][4][1];
      fr=cfaraday[i][j][0][1];
      AA=(dtr*(fl-f0)+dtl*(fr-f0))/(dtl*dtr*(dtl+dtr));
      BB=(-dtr*dtr*(fl-f0)+dtl*dtl*(fr-f0))/(dtl*dtr*(dtl+dtr));
      derf=2.0*AA*dtr+BB;

      
      // J^h = F^{ht},t + F^{hr},r
      jcon[i][j][2]=
	derf + // F^{ht},t
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][1]-geomr.g*cfaraday[i-1][j][1][1])*idx1; // F^{hr},r

      fl=cfaraday[i][j][3][2];
      f0=cfaraday[i][j][4][2];
      fr=cfaraday[i][j][0][2];
      AA=(dtr*(fl-f0)+dtl*(fr-f0))/(dtl*dtr*(dtl+dtr));
      BB=(-dtr*dtr*(fl-f0)+dtl*dtl*(fr-f0))/(dtl*dtr*(dtl+dtr));
      derf=2.0*AA*dtr+BB;

      
      // J^p = F^{pt},t + F^{pr},r + F^{ph},h
      jcon[i][j][3]=
	derf + // F^{pt},t
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][2]-geomrm1.g*cfaraday[i-1][j][1][2])*idx1+ // F^{pr},r
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][2]-geomhm1.g*cfaraday[i][j-1][2][2])*idx2; // F^{ph},h
    }
  }
  else{

    dtl=lastdt;
    dtr=dt;
    idx1=1.0/(2.0*dx[1]);
    idx2=1.0/(2.0*dx[2]);

    LOOPP11 LOOPP12{ // largest possible loop for this differencing (could isolate directions)
      get_geometry(i,j,CENT,&geomt);
      get_geometry(i,j,CENT,&geomr);
      get_geometry(i,j,CENT,&geomh);
      // geomtp1 is same as geomt since d/dt( geometry) -> 0
      get_geometry(i+1,j,CENT,&geomrp1);
      get_geometry(i,j+1,CENT,&geomhp1);
      get_geometry(i-1,j,CENT,&geomrm1);
      get_geometry(i,j-1,CENT,&geomhm1);
      
      // J^t = F^{tr},r + F^{th},h
      jcon[i][j][0]=
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][0]-geomrm1.g*cfaraday[i-1][j][1][0])*idx1+ // F^{tr},r
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][0]-geomhm1.g*cfaraday[i][j-1][2][0])*idx2; // F^{th},h
      
      // J^r = F^{rt},t + F^{rh},h

      jcon[i][j][1]=
	0+
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][1]-geomhm1.g*cfaraday[i][j-1][2][1])*idx2; // F^{rh},h

      // J^h = F^{ht},t + F^{hr},r
      jcon[i][j][2]=
	0+
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][1]-geomr.g*cfaraday[i-1][j][1][1])*idx1; // F^{hr},r

      // J^p = F^{pt},t + F^{pr},r + F^{ph},h
      jcon[i][j][3]=
	0+
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][2]-geomrm1.g*cfaraday[i-1][j][1][2])*idx1+ // F^{pr},r
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][2]-geomhm1.g*cfaraday[i][j-1][2][2])*idx2; // F^{ph},h
    }
  }
  calls++;
  lastdt=dt;

}

// the current is calculated to end up at the zone and time edge
// point is to have J^t at same time as rest of J's, and same spatial points.  time for J is one timestep back for all components.
// like current_calc_0 in time and current_calc_1 in space
void current_calc_2(FTYPE cfaraday[][N2M][NUMCURRENTSLOTS][3])
{
  int i,j;
  struct of_geom geomt;
  struct of_geom geomr;
  struct of_geom geomh;
  struct of_geom geomrp1;
  struct of_geom geomhp1;
  struct of_geom geomrm1;
  struct of_geom geomhm1;
  static FTYPE lastdt;
  static int calls=0;
  FTYPE idtc,idx1,idx2;

  if(calls>0){ // since need 2 times

    idtc=2.0/(lastdt+dt);
    idx1=1.0/(2.0*dx[1]);
    idx2=1.0/(2.0*dx[2]);

    LOOPP11 LOOPP12{ // largest possible loop for this differencing (could isolate directions)
      get_geometry(i,j,CENT,&geomt);
      get_geometry(i,j,CENT,&geomr);
      get_geometry(i,j,CENT,&geomh);
      // geomtp1 is same as geomt since d/dt( geometry) -> 0
      get_geometry(i+1,j,CENT,&geomrp1);
      get_geometry(i,j+1,CENT,&geomhp1);
      get_geometry(i-1,j,CENT,&geomrm1);
      get_geometry(i,j-1,CENT,&geomhm1);
      
      // J^t = F^{tr},r + F^{th},h
      jcon[i][j][0]=
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][0]-geomrm1.g*cfaraday[i-1][j][1][0])*idx1+ // F^{tr},r
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][0]-geomhm1.g*cfaraday[i][j-1][2][0])*idx2; // F^{th},h
      
      // J^r = F^{rt},t + F^{rh},h
      jcon[i][j][1]=
	(cfaraday[i][j][0][0]-cfaraday[i][j][3][0])*idtc+ // F^{rt},t
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][1]-geomhm1.g*cfaraday[i][j-1][2][1])*idx2; // F^{rh},h
      
      // J^h = F^{ht},t + F^{hr},r
      jcon[i][j][2]=
	(cfaraday[i][j][0][1]-cfaraday[i][j][3][1])*idtc+ // F^{ht},t
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][1]-geomr.g*cfaraday[i-1][j][1][1])*idx1; // F^{hr},r
      
      // J^p = F^{pt},t + F^{pr},r + F^{ph},h
      jcon[i][j][3]=
	(cfaraday[i][j][0][2]-cfaraday[i][j][3][2])*idtc+ // F^{pt},t
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][2]-geomrm1.g*cfaraday[i-1][j][1][2])*idx1+ // F^{pr},r
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][2]-geomhm1.g*cfaraday[i][j-1][2][2])*idx2; // F^{ph},h
    }
  }
  else{

    idtc=2.0/(lastdt+dt);
    idx1=1.0/(2.0*dx[1]);
    idx2=1.0/(2.0*dx[2]);

    LOOPP11 LOOPP12{ // largest possible loop for this differencing (could isolate directions)
      get_geometry(i,j,CENT,&geomt);
      get_geometry(i,j,CENT,&geomr);
      get_geometry(i,j,CENT,&geomh);
      // geomtp1 is same as geomt since d/dt( geometry) -> 0
      get_geometry(i+1,j,CENT,&geomrp1);
      get_geometry(i,j+1,CENT,&geomhp1);
      get_geometry(i-1,j,CENT,&geomrm1);
      get_geometry(i,j-1,CENT,&geomhm1);
      
      // J^t = F^{tr},r + F^{th},h
      jcon[i][j][0]=
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][0]-geomrm1.g*cfaraday[i-1][j][1][0])*idx1+ // F^{tr},r
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][0]-geomhm1.g*cfaraday[i][j-1][2][0])*idx2; // F^{th},h
      
      // J^r = F^{rt},t + F^{rh},h
      jcon[i][j][1]=
	0+
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][1]-geomhm1.g*cfaraday[i][j-1][2][1])*idx2; // F^{rh},h
      
      // J^h = F^{ht},t + F^{hr},r
      jcon[i][j][2]=
	0+
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][1]-geomr.g*cfaraday[i-1][j][1][1])*idx1; // F^{hr},r
      
      // J^p = F^{pt},t + F^{pr},r + F^{ph},h
      jcon[i][j][3]=
	0+
	1./geomt.g*(geomrp1.g*cfaraday[i+1][j][1][2]-geomrm1.g*cfaraday[i-1][j][1][2])*idx1+ // F^{pr},r
	1./geomt.g*(geomhp1.g*cfaraday[i][j+1][2][2]-geomhm1.g*cfaraday[i][j-1][2][2])*idx2; // F^{ph},h
    }
  }
  calls++;
  lastdt=dt;

}



// assumes below get inlined

// much faster than macro using ? :
FTYPE sign(FTYPE a)
{
  if(a>0.0) return 1.0;
  else if(a<0.0) return -1.0;
  else return 0.0;
}

// not any faster than above
FTYPE sign2(FTYPE a)
{
#if(SUPERLONGDOUBLE)
  return(sign(a));
  // no such function
#else
  return(copysign(1.0,a));
#endif
}

FTYPE max(FTYPE a, FTYPE b)
{
  if(a>b) return a;
  else return b;
  // if equal, then above is fine
}

FTYPE min(FTYPE a, FTYPE b)
{
  if(a>b) return b;
  else return a;
  // if equal, then above is fine
}


// compute speed of light 3-velocity in particular direction assuming other direction velocities fixed
int sol(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmax, FTYPE *vmin)
{
  int i,j,k;
  FTYPE vu[NDIM],BB,CC,vsol1p,vsol1m;
  FTYPE ftemp1,ftemp2,ftemp3;
  FTYPE disc;
  int diro1,diro2;

  /* 
     m%3+1 gives next 1->2,2->3,3->1
     3-(4-m)%3 gives previous 1->3,2->1,3->2
  */
  diro1=dir%3+1;
  diro2=3-(4-dir)%3;

  // 3-velocity in coordinate frame
  SLOOPA vu[j]=q->ucon[j]/q->ucon[TT];

  ftemp1=0;
  ftemp2=0;
  ftemp3=0;
  SLOOPA{
    if(j!=dir){
      ftemp1+=vu[j]*geom->gcov[dir][j];
      ftemp2+=vu[j]*geom->gcov[TT][j];
      ftemp3+=vu[j]*vu[j]*geom->gcov[j][j];
    }
  }

  BB=2.0*(ftemp1+geom->gcov[0][dir])/geom->gcov[dir][dir];
  CC=(geom->gcov[TT][TT] + 2.0*ftemp2 + ftemp3 + 2.0*vu[diro1]*vu[diro2]*geom->gcov[diro1][diro2])/geom->gcov[dir][dir];
  
  disc=BB*BB-4.0*CC;
  if(disc>0){
    *vmax=0.5*(-BB+sqrt(disc));
    *vmin=0.5*(-BB-sqrt(disc));
  }
  else{
    dualfprintf(fail_file,"disc=%21.15g < 0\n",disc);
    return(1);
  }

  return(0);

}


// limit the 3-velocity to a physically valid velocity (i.e. less than c ), assuming all other velocity directions are the same.
int limitv3(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *v)
{
  FTYPE vmax,vmin;
  int sol(FTYPE *pr, struct of_state *q, int dir, struct of_geom *geom, FTYPE *vmax, FTYPE *vmin);
  FTYPE ratv;
  
  // get speed of light 3-velocity
  MYFUN(sol(pr,q,dir,geom,&vmax,&vmin),"phys.c:limitv3()", "sol()", 1);

  // get ratio of given 3-velocity to speed of light 3-velocity for appropriate direction of coordinate-based velocity
  ratv=(*v>0) ? *v/vmax : *v/vmin;

  // limit 3-velocity to speed of light
  if(ratv>1.0){
    if(*v>0.0) *v=vmax;
    else *v=vmin;
  }

  return(0);
}

// take projection of v onto u, both are 4-vectors
// vcon=0 or 1 : whether or not ucon (1=true, 0=ucov)
// same for vresultcon
void projectionvec(int vcon,int vresultcon, struct of_state *q, struct of_geom *geom,FTYPE *v,FTYPE*vresult)
{
  FTYPE proj[NDIM][NDIM];
  int j,k;


  if((vcon)&&(vresultcon)){ // vresult^\mu = P^\mu_\nu v^\nu
    DLOOP proj[j][k]=delta(j,k) + q->ucon[j]*q->ucov[k];
  }
  if((!vcon)&&(!vresultcon)){ // vresult_\mu = P_\mu^\nu v_\nu
    DLOOP proj[j][k]=delta(j,k) + q->ucov[j]*q->ucon[k];
  }
  else if((!vcon)&&(vresultcon)){ // vresult^\mu = P^{\mu\nu} v_\nu
    DLOOP proj[j][k]=geom->gcon[j][k] + q->ucon[j]*q->ucon[k];
  }
  else if((vcon)&&(!vresultcon)){ // vresult_\mu = P_{\mu\nu} v^\nu
    DLOOP proj[j][k]=geom->gcov[j][k] + q->ucov[j]*q->ucov[k];
  }
  DLOOPA vresult[j]=0.0;
  DLOOP vresult[j]+=proj[j][k]*v[k];
  

}
