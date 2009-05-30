
#include "decs.h"


/*

Generally:

\dF^{\alpha\beta} = (1/2) e^{\alpha\beta\mu\nu} F_{\mu\nu}

F^{\alpha\beta} = (-1/2) e^{\alpha\beta\mu\nu} \dF_{\mu\nu}

e^{\alpha\beta\gamma\delta} = (-1/\detg) [\alpha\beta\gamma\delta]

e_{\alpha\beta\gamma\delta} = (\detg) [\alpha\beta\gamma\delta]

[ijkl] = Signature[{i,j,k,l}]

With

F^{\alpha\beta} = A n^\alpha E^\beta + B n^\beta E^\alpha + C B_\gamma n_\delta e^{\alpha\beta\gamma\delta}

and

\dF^{\alpha\beta} = D n^\alpha B^\beta + F n^\beta B^\alpha + G E_\gamma n_\delta e^{\alpha\beta\gamma\delta}

then antisymmetry gives A=-B and D=-F

where the definition of the dual gives:

C=-F=D
A=-G

which then gives 4 independent equations for 6 unknowns.

Let A=1=-C=-1

then:

F^{\alpha\beta} =  n^\alpha E^\beta - n^\beta E^\alpha - B_\gamma n_\delta e^{\alpha\beta\gamma\delta}

and

\dF^{\alpha\beta} = - n^\alpha B^\beta + n^\beta B^\alpha - E_\gamma n_\delta e^{\alpha\beta\gamma\delta}

As in Ian FFDE paper.

// Gammie:

\dF^{\alpha\beta} = b^\alpha u^\beta - b^\beta u^\alpha = -u^\alpha b^\beta + u^\beta b^\alpha

If \eta->u and then A=1.  Let C=-1 also if setting coefficient of e^\alpha terms.







//Jon notation:

// see grmhd-faraday_allforms.nb

F^{\mu\nu} = -1/\detg {{0,-Eu1,-Eu2,-Eu3},{Eu1,0,-Bd3,Bd2},{Eu2,Bd3,0,-Bd1},{Eu3,-Bd2,Bd1,0}}

\dF_{\mu\nu} = {{0,Bd1,Bd2,Bd3},{-Bd1,0,Eu3,-Eu2},{-Bd2,-Eu3,0,Eu1},{-Bd3,Eu2,-Eu1,0}}


F^{ij} = (1/\detg) [ijk] B_k 

F_{ij}/\detg = B^k [ijk]

[ijk] F_{ij}/\detg = 2B^k


\dF^{ij} = [ijk] E^k and \dF^{ij}[ijk] = 2 E_k


\dF_{ij} = [ijk] E^k


E_i = F_{it}/\detg

E^i = F^{ti}\detg // changed sign at some point (now IanPrimitive^i=F^{it} = -E^i/\detg)

B^i = \dF^{it}  (negative of Ian field Primitive -- IanPrim^i=\dF^{it} = -B^i)

B_i = \dF_{ti} // changed sign at some point

This makes E.B=0


b^\mu = P^\mu_\nu B^\nu/u^t

b_\mu = P_\mu^\nu B_\nu/u_t






*/


// whether to check if B^2-E^2>0 is satisfied
#define CHECKBSQMINUSESQ 1

// whether to account for dissipation of field
#define ACCOUNTDISS 1

// whether to account for lack of energy conservation in accounting
#define ACCOUNTNONCONS 1

// whether to avoid current sheet dissipation by setting the cross-field drift to vanish.
// currently very model dependent
// need a contact discontinuity detector
#define AVOIDCS 1

// whether to fixup B3 to satisfy u.B=0 even when AVOIDCS=1
#define FIXB3 0 // doesn't work -- causes u^t problems


// 0: komissarov
// 1: jon#1
// 2: jon#2 (paper version BEST)
// 3: jon#3 (trial Lorentz factor CS limit)
#define AVOIDUTINF 2


// whether to remove a T^t_i component and replace it with T^t_t
#define REMOVET 0

// assumes geometry-free U in standard form
// latest version
int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr)
{
  void TtoEconOalpha(FTYPE *T0, FTYPE *BcovOalpha, FTYPE BsqOalphasq, struct of_geom *geom, FTYPE *EconOalpha);
  int TtoT(FTYPE BsqOalphasq, FTYPE alphasq, FTYPE *BconOalpha, FTYPE *Tin, struct of_geom *geom, FTYPE *Tout);
  int i,j,k;
  FTYPE alphasq;
  FTYPE Fit[NDIM];
  FTYPE EconOalpha[NDIM],EcovOalpha[NDIM],BcovOalpha[NDIM],BconOalpha[NDIM],realBsqOalphasq;
  FTYPE vcon[NDIM],ucon[NDIM];
  FTYPE BsqOalphasq,EsqOalphasq,BsqmEsqOalphasq;
  FTYPE gammasq,Gsq;
  int modified;
  int diss_account(int modified, FTYPE *Fit,FTYPE *U, struct of_geom *geom,FTYPE*pr);
  void EBcovtovcon(FTYPE alphasq,FTYPE *Ecov,FTYPE *Bcov,FTYPE Bsq, struct of_geom *geom, FTYPE*vcon);
  FTYPE prim[NPR];
  int ti,tj;
  FTYPE ucov[NDIM];
  FTYPE Tin[NDIM],Tout[NDIM];
  void U2BBsq(FTYPE *U, struct of_geom *geom, FTYPE *BconOalpha, FTYPE *BcovOalpha, FTYPE *BsqOalphasq);



  ////////////////////////
  //
  // get lapse
  //
  ///////////////////////

  // define \eta_\alpha
  // assume always has 0 value for space components
  alphasq = 1./(- geom->gcon[0][0]);  
  
  /////////////
  //
  // get (B^\mu/\alpha)
  //
  /////////////

  U2BBsq(U, geom, BconOalpha, BcovOalpha, &BsqOalphasq);

  ///////////////
  //
  // get T^t_\mu
  //
  /////////////
  DLOOPA Tin[j] = Tout[j] = U[UU+j] ; /* T^t_\mu terms (t up, \mu down) */

#if(REMOVET)
  TtoT(BsqOalphasq,alphasq, BconOalpha, Tin, geom, Tout);
#endif

  /////////////////
  //
  // get (E^\mu/\alpha)
  //
  /////////////////
  TtoEconOalpha(Tout, BcovOalpha, BsqOalphasq, geom, EconOalpha);


  /////////////////////////
  //
  // modify velocity so B^2-E^2>0 always
  //
  ////////////////////////

  // set no modification to energy
  modified=0;

  realBsqOalphasq=BsqOalphasq; // default
#if(CHECKBSQMINUSESQ)



  lower(EconOalpha,geom,EcovOalpha);

  EsqOalphasq = 0.0;
  SLOOPA EsqOalphasq += EconOalpha[j]*EcovOalpha[j];

  BsqmEsqOalphasq = BsqOalphasq-EsqOalphasq;
  
  //check if effective dissipation is occuring

#if(AVOIDUTINF==0)
  if(BsqmEsqOalphasq<0.0){
    //    dualfprintf(fail_file,"You are going to fail! %21.15g\n",BsqmEsqOalphasq);
    // then modify drift velocity
    // see Komissarov (2005)
    realBsqOalphasq=(max(BsqOalphasq,EsqOalphasq));
    modified=1;
  }
#elif(AVOIDUTINF==1)
  gammasq=BsqOalphasq/BsqmEsqOalphasq;
  if((gammasq>GAMMAMAX*GAMMAMAX)||(gammasq<1.0)){
    // then limit the Lorentz factor
    //  gamma = sqrt(Bsq/fabs(BsqmEsqOalphasq));
    Gsq=1.0/(GAMMAMAX*GAMMAMAX);
    realBsqOalphasq=(max(BsqOalphasq,EsqOalphasq)*sqrt(1.0+Gsq));
    // else assume fine
    modified=1;
  }
#elif(AVOIDUTINF==2)
  gammasq=BsqOalphasq/BsqmEsqOalphasq;
  if((gammasq>GAMMAMAX*GAMMAMAX)||(gammasq<1.0)){
    Gsq=1.0-1.0/(GAMMAMAX*GAMMAMAX);
    realBsqOalphasq=(sqrt(EsqOalphasq*BsqOalphasq/Gsq));
    modified=1;
  }
#elif(AVOIDUTINF==3)
  // limit Lorentz factor to 1 in current sheets
  // doesn't make sense to limit completely like this since directional dependent
  if((tj>totalsize[2]/2-3)&&(tj<totalsize[2]/2+2)){
    // then limit the Lorentz factor
    //  gamma = sqrt(Bsq/fabs(BsqmEsqOalphasq));
    Gsq=1.0/(1.0);
    realBsqOalphasq=(max(BsqOalphasq,EsqOalphasq)*sqrt(1.0+Gsq));
    // else assume fine
    modified=1;
  }
#endif

#endif // end if checking


  /////////////////////
  //
  // Find v^i
  //
  /////////////////////

  // below assumes E_i B^i = F_{it} \dF^{it} = 0, which it is.
  // note that the velocity is defined only up to a boost along field
  
  // actually, we removed 1/alpha^2 from Bsq thing, and removed alpha from Bcov and Ecov each (total alpha^2), so left with alpha^2
  // so alpha cancels between Ecov*Bcov/Bsq, so removed versions and versions with it in use same function
  EBcovtovcon(alphasq, EcovOalpha, BcovOalpha, realBsqOalphasq, geom, vcon);



  ///////////////
  //
  // avoid current sheets
  //
  /////////////
  // tried various ways to avoid dissipation in current sheet (CS).
  // 1) use gamma~1 in CS : kinda worked
  // 2) set B2=0 in CS: didn't work
  // 3) set vcon2=0 in CS: works very well!


  // symmetric around equator for centered zones
  ti=startpos[1]+geom->i;
  tj=startpos[2]+geom->j;
  if(AVOIDCS&&((tj>totalsize[2]/2-3)&&(tj<totalsize[2]/2+2))){
    vcon[TH]=0.0;
    modified=1;
  }

  //////////////////////////
  //
  // convert from v^i -> u^\mu -> primitive velocity
  //
  //////////////////

  //  PLOOP prim[k]=0.0;
  //  prim[U1]=vcon[1];
  //  prim[U2]=vcon[2];
  //  prim[U3]=vcon[3];

  vcon2pr(WHICHVEL,vcon,geom,pr);

  //  ucon_calc_3vel(prim,geom,ucon);
  //  pr2ucon(VEL3,prim,geom,ucon); // if u^t can't be found, then FFDE breaks down (i.e. dF^2>0 is ok) 

  // convert from u^\mu -> true WHICHVEL primitive
  //  ucon2pr(WHICHVEL,ucon,geom,pr); // fills pr[U1->U3]

  //////////////////////////
  //
  // now assign field primitives
  //
  ///////////////////////
  pr[B1]=U[B1];
  pr[B2]=U[B2];
  pr[B3]=U[B3];

  //////////////////
  //
  // rearrange B^\phi so that u.B=0 even after current sheet changes
  //
  /////////////////
  if(FIXB3&&AVOIDCS&&((tj>totalsize[2]/2-3)&&(tj<totalsize[2]/2+2))){
    // doesn't work
    PLOOP prim[k]=0.0;
    prim[U1]=vcon[1];
    prim[U2]=vcon[2];
    prim[U3]=vcon[3];
    ucon_calc_3vel(prim,geom,ucon);
    lower(ucon,geom,ucov);
    pr[B3]=-ucov[RR]/ucov[PH]*pr[B1];
    // this only works in axisymmetry (else divb!=0)
  }
    


  //////////////////
  //
  // Account for any modifications
  //
  /////////////////

#if(ACCOUNTDISS)
  // account for any energy change
  if(1||modified){ // 1|| because if not conserving energy, then really changed

    // used by diagnostics
    //SLOOPA Econ[j]=-Fit[j];
    Fit[TT]=0.0;
    SLOOPA Fit[j]=-EconOalpha[j];

    diss_account(modified,Fit,U,geom,pr);
  }
#endif
  // done.


  return(0);


}




// assumes geometry-free U in standard form
int Utoprim_ffde_oldnew(FTYPE *U, struct of_geom *geom, FTYPE *pr)
{
  int i,j,k;
  void UtoFitMxBsqgeom(FTYPE *U, struct of_geom *geom, FTYPE *Fit, FTYPE*Mx, FTYPE*oBsqgeom);
  FTYPE alphasq,beta[NDIM];
  FTYPE Fit[NDIM],Mx[NDIM],Mt[NDIM],oBsqgeom;
  FTYPE Econ[NDIM],Ecov[NDIM],Bcov[NDIM],Bcon[NDIM],realoBsqgeom;
  FTYPE vcon[NDIM],ucon[NDIM];
  FTYPE Bsq,Esq,BsqmEsq;
  FTYPE gammasq,Gsq;
  int modified;
  int diss_account(int modified, FTYPE *Fit,FTYPE *U, struct of_geom *geom,FTYPE*pr);
  void EBcovtovcon_old(FTYPE *beta, FTYPE alphasq,FTYPE *Ecov,FTYPE *Bcov,FTYPE realoBsqgeom, FTYPE*vcon);
  //  FTYPE prim[NPR];
  int ti,tj;

  ti=startpos[1]+geom->i;
  tj=startpos[2]+geom->j;

  // set no modification to energy
  modified=0;

  // lapse
  // define \eta_\alpha
  // assume always has 0 value for space components
  alphasq = 1./(- geom->gcon[0][0]);  
  
  SLOOPA beta[j] = geom->gcon[0][j]*alphasq ;

  UtoFitMxBsqgeom(U,geom,Fit,Mx,&oBsqgeom);
  // -> B_\mu = -alpha * M^t_\mu
  // Real Bsq = Bsqgeom*alpha^2

  Econ[TT]=0.0;
  // remove alpha since cancels below
  //  SLOOPA Econ[j]=-alpha*Fit[j];
  SLOOPA Econ[j]=-Fit[j];
  lower(Econ,geom,Ecov);

  // Mx=M^t_\mu = g_{\mu\nu} M^{t\nu} = g_{\mu\nu} B^\nu/(-alpha)
  // remove alpha since cancels below
  //SLOOPA Bcov[j]=-alpha*Mx[j];
  SLOOPA Bcov[j]=-Mx[j]; // Mx = \dF^t_j

  // remove alpha^2 since cancels below
  //  realoBsqgeom=oBsqgeom/alphasq;
  realoBsqgeom=oBsqgeom;

#if(CHECKBSQMINUSESQ)
  // check B^2-E^2
  Bcon[TT]=0.0;
  // Bcon without alpha
  SLOOPA Bcon[j] = U[B1+j-1] ;  /* Bcon=-M^{ti} , where U[Bi]=B^i=M^{it}*/

  Esq = 0.0;
  SLOOPA Esq += Econ[j]*Ecov[j];

  Bsq = 0.0;
  SLOOPA Bsq += Bcon[j]*Bcov[j];

  BsqmEsq = Bsq-Esq;

  // tried various ways to avoid dissipation in current sheet (CS).
  // 1) use gamma~1 in CS : kinda worked
  // 2) set B2=0 in CS: didn't work
  // 3) set vcon2=0 in CS: works very well!

  
  //check if effective dissipation is occuring

#if(AVOIDUTINF==0)
  if(BsqmEsq<0.0){
    //    dualfprintf(fail_file,"You are going to fail! %21.15g\n",BsqmEsq);
    // then modify drift velocity
    // see Komissarov (2005)
    realoBsqgeom=1.0/(max(Bsq,Esq)*geom->g);
    modified=1;
  }
#elif(AVOIDUTINF==1)
  gammasq=Bsq/BsqmEsq;
  if((gammasq>GAMMAMAX*GAMMAMAX)||(gammasq<1.0)){
    // then limit the Lorentz factor
    //  gamma = sqrt(Bsq/fabs(BsqmEsq));
    Gsq=1.0/(GAMMAMAX*GAMMAMAX);
    realoBsqgeom=1.0/(max(Bsq,Esq)*sqrt(1.0+Gsq)*geom->g);
    // else assume fine
    modified=1;
  }
#elif(AVOIDUTINF==2)
  gammasq=Bsq/BsqmEsq;
  if((gammasq>GAMMAMAX*GAMMAMAX)||(gammasq<1.0)){
    Gsq=1.0-1.0/(GAMMAMAX*GAMMAMAX);
    realoBsqgeom=1.0/(sqrt(Esq*Bsq/Gsq)*geom->g);
    modified=1;
  }
#elif(AVOIDUTINF==3)
  // limit Lorentz factor to 1 in current sheets
  // doesn't make sense to limit completely like this since directional dependent
  if((tj>totalsize[2]/2-3)&&(tj<totalsize[2]/2+2)){
    // then limit the Lorentz factor
    //  gamma = sqrt(Bsq/fabs(BsqmEsq));
    Gsq=1.0/(1.0);
    realoBsqgeom=1.0/(max(Bsq,Esq)*sqrt(1.0+Gsq)*geom->g);
    // else assume fine
    modified=1;
  }
#endif

#endif // end if checking

  /////////////////////
  //
  // Find v^i
  //
  // below assumes E_i B^i = F_{it} \dF^{it} = 0, which it is.
  // note that the velocity is defined only up to a boost along field
  
  // actually, we removed 1/alpha^2 from Bsq thing, and removed alpha from Bcov and Ecov each (total alpha^2), so left with alpha^2
  //  vcon[1]=-beta[1]+alphasq*(Ecov[2]*Bcov[3]-Ecov[3]*Bcov[2])*realoBsqgeom;
  //  vcon[2]=-beta[2]+alphasq*(Ecov[3]*Bcov[1]-Ecov[1]*Bcov[3])*realoBsqgeom;
  //  vcon[3]=-beta[3]+alphasq*(Ecov[1]*Bcov[2]-Ecov[2]*Bcov[1])*realoBsqgeom;
  //  vcon[1]=-beta[1]+alphasq*(Ecov[2]*Bcov[3]-Ecov[3]*Bcov[2])*realoBsqgeom;
  //  vcon[2]=-beta[2]+alphasq*(Ecov[3]*Bcov[1]-Ecov[1]*Bcov[3])*realoBsqgeom;
  //  vcon[3]=-beta[3]+alphasq*(Ecov[1]*Bcov[2]-Ecov[2]*Bcov[1])*realoBsqgeom;


  EBcovtovcon_old(beta,alphasq,Ecov,Bcov,realoBsqgeom,vcon);

  // symmetric around equator for centered zones
  if(AVOIDCS&&((tj>totalsize[2]/2-3)&&(tj<totalsize[2]/2+2))){
    vcon[2]=0.0;
    modified=1;
  }
  //////////////////////////
  //
  // convert from v^i -> u^\mu
  //  PLOOP prim[k]=0.0;
  //  prim[U1]=vcon[1];
  //  prim[U2]=vcon[2];
  //  prim[U3]=vcon[3];

  vcon2pr(WHICHVEL,vcon,geom,pr);

  //  ucon_calc_3vel(prim,geom,ucon);
  //  pr2ucon(VEL3,prim,geom,ucon); // if u^t can't be found, then FFDE breaks down (i.e. dF^2>0 is ok) 

  // convert from u^\mu -> true WHICHVEL primitive
  //  ucon2pr(WHICHVEL,ucon,geom,pr); // fills pr[U1->U3]

  // now assign field primitives
  pr[B1]=U[B1];
  pr[B2]=U[B2];
  pr[B3]=U[B3];
  


#if(ACCOUNTDISS)
  // account for any energy change
  if(1||modified){ // 1|| because if not conserving energy, then really changed
    diss_account(modified, Fit,U,geom,pr);
  }
#endif
  // done.


  return(0);


}

void EBcovtovcon(FTYPE alphasq,FTYPE *Ecov,FTYPE *Bcov,FTYPE Bsq, struct of_geom *geom, FTYPE*vcon)
{
  FTYPE oBsqgeom;
  
  oBsqgeom=1.0/(Bsq*geom->g);
  
  vcon[1]=alphasq*(-geom->gcon[TT][1] + (Ecov[2]*Bcov[3]-Ecov[3]*Bcov[2])*oBsqgeom);
  vcon[2]=alphasq*(-geom->gcon[TT][2] + (Ecov[3]*Bcov[1]-Ecov[1]*Bcov[3])*oBsqgeom);
  vcon[3]=alphasq*(-geom->gcon[TT][3] + (Ecov[1]*Bcov[2]-Ecov[2]*Bcov[1])*oBsqgeom);
}

void EBcovtovcon_old(FTYPE *beta, FTYPE alphasq,FTYPE *Ecov,FTYPE *Bcov,FTYPE oBsqgeom, FTYPE*vcon)
{
  vcon[1]=-beta[1]+alphasq*(Ecov[2]*Bcov[3]-Ecov[3]*Bcov[2])*oBsqgeom;
  vcon[2]=-beta[2]+alphasq*(Ecov[3]*Bcov[1]-Ecov[1]*Bcov[3])*oBsqgeom;
  vcon[3]=-beta[3]+alphasq*(Ecov[1]*Bcov[2]-Ecov[2]*Bcov[1])*oBsqgeom;
}



// account for dissipation or truncation level lack of conservation
int diss_account(int modified, FTYPE *Fit,FTYPE *U, struct of_geom *geom,FTYPE*pr)
{
  int j;
  void FitUtoT(FTYPE *Fit,FTYPE *U,struct of_geom *geom, FTYPE T[NDIM][NDIM]);
  extern int diag_fixup_U(FTYPE *Ui, FTYPE *Uf, struct of_geom *ptrgeom, int finalstep,int whocalled);
  extern int primtoU(int returntype, FTYPE *pr, struct of_state *q, struct of_geom *geom,FTYPE *U);
  struct of_state q;
  FTYPE T[NDIM][NDIM];
  FTYPE Ui[NPR],Uf[NPR];


  // do this to get initial T consistent with momentum (excludes energy) equations
  FitUtoT(Fit,U,geom,T);

  // initial energy-momentum and field
  Ui[RHO]=0.0;
  DLOOPA Ui[UU+j]=geom->g*T[TT][j];
  SLOOPA Ui[B1+j-1]=geom->g*U[B1+j-1];
  
  // get final energy-momentum and field
  MYFUN(get_state(pr, geom, &q),"phys.ffde.c:Utoprim_ffde()", "get_state()", 1);
  primtoU(UDIAG,pr,&q,geom,Uf);
  Uf[RHO]=0.0;

#if(ACCOUNTNONCONS)
  Uf[UU]-=(geom->g*U[UU]-Ui[UU]); // add back in lost energy from nonconservation of energy method
#endif
  
  // account for change
  // diag_fixup_U is forced to account on every time substep since can't step through problems like with negative internal energy.
  if(modified) diag_fixup_U(Ui,Uf,geom,0,COUNTUTOPRIMFAILCONV); // real dissipation, absolute limit of Lorentz factor
  else diag_fixup_U(Ui,Uf,geom,0,COUNTFLOORACT); // truncation dissipation
   

  return(0);
}  



// input E and B are defined as in GRFFE paper
void EBtopr_2(FTYPE *Econ,FTYPE *Bcon,struct of_geom *geom, FTYPE *pr)
{
  int j;
  FTYPE alphasq,beta[NDIM],alpha;
  FTYPE Bsq,realoBsqgeom;
  FTYPE Ecov[NDIM],Bcov[NDIM];
  FTYPE vcon[NDIM];
  void EBcovtovcon_old(FTYPE *beta, FTYPE alphasq,FTYPE *Ecov,FTYPE *Bcov,FTYPE realoBsqgeom, FTYPE*vcon);

  

  alphasq = 1./(- geom->gcon[0][0]);  
  alpha=sqrt(alphasq);
  SLOOPA beta[j] = geom->gcon[0][j]*alphasq ;

  lower(Econ,geom,Ecov);
  lower(Bcon,geom,Bcov);

  Bsq = 0.0;
  SLOOPA Bsq += Bcon[j]*Bcov[j];
  realoBsqgeom=1.0/(Bsq*geom->g);

  EBcovtovcon_old(beta,alphasq,Ecov,Bcov,realoBsqgeom,vcon);
  vcon2pr(WHICHVEL,vcon,geom,pr);

  pr[B1]=-Bcon[1]/(-alpha); // \dF^{it}
  pr[B2]=-Bcon[2]/(-alpha);
  pr[B3]=-Bcon[3]/(-alpha);

}


// convert E^\mu and B^\mu to WHICHVEL v^i and B^i (code's primitive quantities)
// input E and B are defined as in GRFFE paper
void EBtopr(FTYPE *Econ,FTYPE *Bcon,struct of_geom *geom, FTYPE *pr)
{
  void EBtoT(FTYPE *Econ,FTYPE *Bcon,struct of_geom *geom, FTYPE T[NDIM][NDIM]);
  int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr);
  FTYPE T[NDIM][NDIM];
  FTYPE U[NPR];
  FTYPE alpha;

  // basically convert T and B -> U

  EBtoT(Econ,Bcon,geom,T);

  alpha = sqrt(1./(- geom->gcon[0][0]));  

  U[RHO]=0.0;
  U[UU]=T[TT][TT]; // T^t_\mu
  U[U1]=T[TT][RR];
  U[U2]=T[TT][TH];
  U[U3]=T[TT][PH];
  U[B1]=-Bcon[1]/(-alpha); // \dF^{it}
  U[B2]=-Bcon[2]/(-alpha);
  U[B3]=-Bcon[3]/(-alpha);

  Utoprim_ffde(U, geom, pr);



}

// convert E^\mu and B^\mu to T
// input E and B are defined as in GRFFE paper
void EBtoT(FTYPE *Econ,FTYPE *Bcon,struct of_geom *geom, FTYPE T[NDIM][NDIM])
{
  void FitUtoT(FTYPE *Fit,FTYPE *U,struct of_geom *geom, FTYPE T[NDIM][NDIM]);
  FTYPE Fit[NDIM];
  FTYPE U[NPR];
  FTYPE alpha;
  

  alpha = sqrt(1./(- geom->gcon[0][0]));  

  Fit[TT]=0;
  Fit[1]=Econ[1]/(-alpha); // F^{it}
  Fit[2]=Econ[2]/(-alpha);
  Fit[3]=Econ[3]/(-alpha);
  U[B1]=-Bcon[1]/(-alpha); // \dF^{it}
  U[B2]=-Bcon[2]/(-alpha);
  U[B3]=-Bcon[3]/(-alpha);
  
  
  FitUtoT(Fit,U,geom,T);


}

// v.B = KK B^2/\sqrt{B^2}
// component of velocity along field (assumed to be zero in GRFFE formulation)
void computeKK(FTYPE *pr, struct of_geom *geom, FTYPE *KK)
{
  FTYPE ucon[NDIM],ucov[NDIM],Bcon[NDIM],Bcov[NDIM];
  FTYPE Bsq;
  int j;

  ucon_calc(pr,geom,ucon);
  lower(ucon,geom,ucov);

  SLOOPA Bcon[j]=pr[B1+j-1];
  lower(Bcon,geom,Bcov);

  Bsq=0.0;
  SLOOPA Bsq+=Bcon[j]*Bcov[j];

  *KK=0.0;
  SLOOPA *KK+=(ucov[j]/ucov[TT])*Bcon[j];
  *KK *=sqrt(1.0/Bsq);
  

}


// for E and B in arbitrary frame moving with coodinate 3-velocity veta^i
void EBvetatopr(FTYPE *Econ, FTYPE *Bcon, FTYPE *veta, struct of_geom *geom, FTYPE *pr)
{
  int j,k,l,m;
  FTYPE T[NDIM][NDIM];
  FTYPE Fcon[NDIM][NDIM],Fcov[NDIM][NDIM],Fud[NDIM][NDIM];
  FTYPE Mcon[NDIM][NDIM];
  FTYPE prim[NPR];
  FTYPE Bcov[NDIM];
  FTYPE etacon[NDIM],etacov[NDIM];
  FTYPE Fsq;
  FTYPE alpha,etacovzamo[NDIM],etaconzamo[NDIM];
  FTYPE beta[NDIM];
  FTYPE U[NPR];

  extern FTYPE lc4(int updown, FTYPE detg, int mu,int nu,int kappa,int lambda);
  void lower_A(FTYPE Acon[NDIM][NDIM], struct of_geom *geom, FTYPE Acov[NDIM][NDIM]);
  int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr);
  void FitUtoT(FTYPE *Fit,FTYPE *U,struct of_geom *geom, FTYPE T[NDIM][NDIM]);
  void MtoF(int which, FTYPE Max[NDIM][NDIM],struct of_geom *geom, FTYPE faraday[NDIM][NDIM]);
  //\begin{equation}\label{faraday}
  //F^{\alpha\beta} \equiv f\eta^\alpha E^\beta - f\eta^\beta E^\alpha -
  //h B_\gamma \eta_\delta \epsilon^{\alpha\beta\gamma\delta}
  //\end{equation}
  //where $f$ and $h$ are arbitrary independent constants that we set to
  //be $f=h=1$ consistent with the conventions in MTW (see also, e.g.,

  //\begin{equation}\label{tmunu}
  //T^{\mu\nu} = {F^\mu}_\lambda F^{\nu\lambda} - \frac{1}{4}g^{\mu\nu}
  //F^{\lambda\kappa} F_{\lambda\kappa} ,
  //\end{equation}

  PLOOP prim[k]=0.0;
  prim[U1]=veta[1];
  prim[U2]=veta[2];
  prim[U3]=veta[3];

  // arbitrary frame where Bcon and Econ are measured
  ucon_calc_3vel(prim,geom,etacon);
  lower(etacon,geom,etacov);

  // B_\mu
  Bcon[TT]=0.0; // force
  lower(Bcon,geom,Bcov);

  Econ[TT]=0.0; // force


  // F^{\mu\nu}
  DLOOP Fcon[j][k]=0.0;
  DLOOP for(l=0;l<NDIM;l++) for(m=0;m<NDIM;m++){
    Fcon[j][k]+=-Bcov[l]*etacov[m]*lc4(1,geom->g,j,k,l,m); // \epsilon^{jklm}
  }
  DLOOP Fcon[j][k] += etacon[j]*Econ[k]-etacon[k]*Econ[j];

  // F_mu^\nu
  DLOOPA lower(Fcon[j],geom,Fud[j]);

  // F_{\mu\nu}
  lower_A(Fcon,geom,Fcov);



  // T^\mu_\nu
  DLOOP T[j][k]=0.0;

  Fsq=0.0;
  DLOOP Fsq+=Fcon[j][k]*Fcov[j][k];

  DLOOP{
    T[j][k]+=-0.25*delta(j,k)*Fsq;
    for(l=0;l<NDIM;l++) T[j][k]+=(-Fud[l][j])*(Fud[k][l]);
  }


  //  DLOOP dualfprintf(fail_file,"Fcov[%d][%d]=%15.7g\n",j,k,Fcov[j][k]);

  MtoF(3,Fcov,geom,Mcon);

  U[RHO]=0.0;
  DLOOPA U[UU+j]=T[TT][j]; // T^t_\mu
  SLOOPA U[B1+j-1]=Mcon[j][TT]; // Jon's B^i = \dF^{it}
  
  //  PLOOP dualfprintf(fail_file,"U[%d]=%21.15g\n",k,U[k]);

  Utoprim_ffde(U, geom, pr);

}


// needed in case T^t_t isn't consistent with T^i_t
void FitUtoT(FTYPE *Fit,FTYPE *U,struct of_geom *geom, FTYPE T[NDIM][NDIM])
{
  FTYPE prffde[NPR];
  FTYPE Mcon[NDIM][NDIM] ; /* contravariant Maxwell */
  void ffdestresstensor(FTYPE Mcon[NDIM][NDIM], struct of_geom *geom, FTYPE T[NDIM][NDIM]);
  void Max_con(FTYPE prffde[NPR], struct of_geom *geom, FTYPE Mcon[NDIM][NDIM]);


    // get initial U (only from E and B, not from T^t_t)
    prffde[U1]=Fit[1];
    prffde[U2]=Fit[2];
    prffde[U3]=Fit[3];
    prffde[B1]=-U[B1];  // M^{ti} since  U(geomfree)=B^i = \dF^{it} = M^{it}
    prffde[B2]=-U[B2];
    prffde[B3]=-U[B3];
    Max_con(prffde, geom, Mcon);
    /* T^i_t terms (i up, t down) */
    ffdestresstensor(Mcon, geom, T);
}


// 1=x, 2=y, 3=z
#define DEFAULTKILL 2 // choose T^t_y to kill

// convert T^t_t,T^t_y to T^t_x , where y is 2 spatial dimensions and x is one spatial dimension.  Eliminates one spatial dimension for T^t_t.
int TtoT(FTYPE BsqOalphasq, FTYPE alphasq, FTYPE *BconOalpha, FTYPE *Tin, struct of_geom *geom, FTYPE *Tout)
//int TtoT(FTYPE Bsq,FTYPE alphasq, FTYPE *Bcon, FTYPE *Uin, struct of_geom *geom, FTYPE *Uout)
{
  void getABC(int which, FTYPE BsqOalphasq,FTYPE alphasq, FTYPE *BconOalpha, FTYPE *Tin, struct of_geom *geom, FTYPE *AA, FTYPE*BB, FTYPE*CC);
  int j,k;
  FTYPE subEsq,EsqOBsq;
  FTYPE det,sqrtdet;
  FTYPE AA,BB,CC,qq;
  FTYPE replacea,replaceb;
  int which;


  // only do this if E^2 not close to 0
  subEsq=0.0; SLOOPA subEsq+=Tin[j]*geom->gcon[TT][j];
  EsqOBsq=2.0*(-Tin[TT]/alphasq+subEsq)/BsqOalphasq-1.0;
  // see if not enough information to change T
  // also check if violating causality, in which case dissipation is taken care of elsewhere.
  if( ((EsqOBsq)<1E-13)||(EsqOBsq>1.0) ) return(0); 

  // choose which T^t_{which} to remove 
  which=DEFAULTKILL;

  // set default as out=in
  DLOOPA Tout[j]=Tin[j]; // j=0..3

  // get quadratic coefficients for quadratic equation with solution T^t_{which} = (-BB\pm\sqrt{BB^2-4*AA*CC})/(2*AA) .
  getABC(which, BsqOalphasq,alphasq,BconOalpha,Tin,geom,&AA,&BB,&CC);

  det=BB*BB-4.0*AA*CC;
  if(det<0.0){
    //    dualfprintf(fail_file,"t=%g i=%d j=%d : TtoT det=%g\n",t,geom->i,geom->j,det);
    //    myexit(1);
    //dualfprintf(fail_file,"t=%g i=%d j=%d det=%21.15g<0: Tin[%d]=%21.15g Tout[%d]=%21.15g\n",t,geom->i,geom->j,det,which,Tin[which],which,Tout[which]);
    return(0);
  }
  if(AA==0.0){
    //    dualfprintf(fail_file,"t=%g i=%d j=%d : TtoT AA=%g\n",t,geom->i,geom->j,AA);
    //    myexit(1);
    dualfprintf(fail_file,"t=%g i=%d j=%d AA==%21.15g: Tin[%d]=%21.15g Tout[%d]=%21.15g\n",t,geom->i,geom->j,AA,which,Tin[which],which,Tout[which]);
    return(0);
  }
  if(BB==0.0){
    sqrtdet=sqrt(det);
    replacea=(-BB+sqrtdet)/(2.0*AA);
    replaceb=(-BB-sqrtdet)/(2.0*AA);
  }
  else{
    sqrtdet=sqrt(det);
    qq=-0.5*(BB+sign(BB)*sqrtdet);
    if(BB>0.0){
      replacea=CC/qq;
      replaceb=qq/AA;
    }
    else{
      replacea=CC/qq;
      replaceb=qq/AA;
    }
  }


  //Tout[which]=replaceb;
  //Tout[which]=replacea;
  //  dualfprintf(fail_file,"t=%g i=%d j=%d qq=%21.15g AA=%21.15g BB=%21.15g CC=%21.15g \n det=%g signBB=%g \n Tin[%d]=%21.15g (A)Tout[%d]=%21.15g (B)Tout[%d]=%21.15g\n",t,geom->i,geom->j,qq,AA,BB,CC, det,sign(BB),  which,Tin[which],which,replacea,which,replaceb);
  //  return(0);


  // figure out which one is closest to original to get directionality information
  if(fabs(replacea-Tin[which])<fabs(replaceb-Tin[which])){
    Tout[which]=replacea;
    //    dualfprintf(fail_file,"t=%g i=%d j=%d a: Tin[%d]=%21.15g Tout[%d]=%21.15g\n",t,geom->i,geom->j,which,Tin[which],which,Tout[which]);
  }
  else{
    Tout[which]=replaceb;
    //    dualfprintf(fail_file,"t=%g i=%d j=%d b: Tin[%d]=%21.15g Tout[%d]=%21.15g\n",t,geom->i,geom->j,which,Tin[which],which,Tout[which]);
  }

  return(0);

}

void getABC(int which, FTYPE BsqOalphasq,FTYPE alphasq, FTYPE *BconOalpha, FTYPE *Tin, struct of_geom *geom, FTYPE *AA, FTYPE*BB, FTYPE*CC)
{
  int j,k;
  FTYPE Q[NDIM][NDIM];
  FTYPE subBB,subCC1,subCC2;



  SSLOOP Q[j][k]=(geom->gcon[j][k])+(alphasq*geom->gcon[TT][j]*geom->gcon[TT][k])-BconOalpha[j]*BconOalpha[k]/BsqOalphasq;

  *AA=-1.0/(2.0*BsqOalphasq)*Q[which][which];

  subBB=0; SLOOPA if(j!=which) subBB+=Tin[j]*Q[which][j];
  *BB=alphasq*(geom->gcon[TT][which]) - subBB/BsqOalphasq;

  subCC1=0; SLOOPA if(j!=which) subCC1+=geom->gcon[TT][j]*Tin[j];
  subCC2=0; SSLOOP if((j!=which)&&(k!=which)) subCC2+=Tin[j]*Tin[k]*Q[j][k];
  *CC=-Tin[TT] + alphasq*subCC1 - alphasq*BsqOalphasq*0.5 - subCC2/(2.0*BsqOalphasq);
  
}





// assumes geometry-free U in standard form
int Utoprim_ffde_old(FTYPE *U, struct of_geom *geom, FTYPE *pr)
{
  int i,j,k;
  FTYPE T[NDIM][NDIM];
  FTYPE Ti0[NDIM] ;   /* stress tensor terms T^i_t */
  FTYPE Fit[NDIM],prffde[NPR];
  FTYPE Mcon[NDIM][NDIM] ; /* contravariant Maxwell */
  FTYPE Mcov[NDIM][NDIM] ; /* covariant Maxwell */
  FTYPE ucon[NDIM];
  FTYPE prim[NPR];
  FTYPE Fcov[NDIM][NDIM];
  FTYPE Fcon[NDIM][NDIM];
  FTYPE Bsq;
  FTYPE vcon[NDIM];
  void UtoFit(FTYPE *U, struct of_geom *geom, FTYPE *Fit);
  void Max_con(FTYPE prffde[NPR], struct of_geom *geom, FTYPE Mcon[NDIM][NDIM]);
  void lower_A(FTYPE Acon[NDIM][NDIM], struct of_geom *geom, FTYPE Acov[NDIM][NDIM]);
  void MtoF(int which, FTYPE Max[NDIM][NDIM],struct of_geom *geom, FTYPE faraday[NDIM][NDIM]);
  void ffdestresstensor(FTYPE Mcon[NDIM][NDIM], struct of_geom *geom, FTYPE T[NDIM][NDIM]);
  void raise_A(FTYPE Acov[NDIM][NDIM], struct of_geom *geom, FTYPE Acon[NDIM][NDIM]);
    FTYPE alpha;
      FTYPE betai[NDIM];
      FTYPE etacov[NDIM],etacon[NDIM];
      FTYPE Ecov[NDIM],Bcov[NDIM],Econ[NDIM],Bcon[NDIM];

      // lapse
      // define \eta_\alpha
      // assume always has 0 value for space components
      alpha = 1./sqrt(- geom->gcon[0][0]);  

      etacov[TT]=-alpha; // any constant will work.
      SLOOPA etacov[j]=0.0; // must be 0

      // shift
      /* beta */
      //      SLOOPA beta[j] = geom->gcon[0][j]*alpha*alpha ;
      // define \eta^\beta
      raise(etacov,geom,etacon);

      SLOOPA betai[j]=-etacon[j]*alpha;


  // get F^{it}
  UtoFit(U,geom,Fit);



  // get M^{\mu\nu}
  // correct order, sign, etc. for Max_con
  prffde[U1]=Fit[1];
  prffde[U2]=Fit[2];
  prffde[U3]=Fit[3];
  prffde[B1]=-U[B1];  // M^{ti} since  U(geomfree)=B^i = \dF^{it} = M^{it}
  prffde[B2]=-U[B2];
  prffde[B3]=-U[B3];

  Max_con(prffde, geom, Mcon);

  // for alternative formulation
  /* T^i_t terms (i up, t down) */
  ffdestresstensor(Mcon, geom, T);
  SLOOPA Ti0[j] = T[j][TT] ; // T^i_t

  /* get covariant Maxwell from contravariant */
  lower_A(Mcon, geom, Mcov) ;

  // get F_{\mu\nu}
  MtoF(0,Mcon,geom,Fcov);

  /* get covariant Maxwell from contravariant */
  raise_A(Fcov, geom, Fcon) ;


  // now can get velocity
  Bsq=0.0;
  DLOOPA Bsq+=-Mcon[j][TT]*Mcov[j][TT]; // sign of Mcov_{jt} consistent with used by Mcov in vcon below.  That's all that matters.

  // v^i = E_j B_k [ijk]/(B^l B_l) in Jon notation

  // S^i = -[ijk] E_j B_k (in Jon notation)
  // S_i = [ijk] E^j B^k = [ijk] \detg F^{tj} \dF^{kt} (Jon notation)
  // T^t_i = - sqrt(-g) * [ijk] * M^tj * F^kt = S_i

  // so v^i = -S^i/(B^l B_l) where S^i = T^i_t (can construct from T^{\mu\nu}) and compare with E/B form


  // below assumes E_i B^i = F_{it} \dF^{it} = 0, which it is.
  // note that the velocity is defined only up to a boost along field
  // that is, the below results in \dF_{it} v^i = 0 and F_{it} v^i = 0
  // notice that for v^i that \eta_t constant!=0 doesn't matter (cancels out)
#if(1)
  DLOOPA Ecov[j]=0;
  DLOOP Ecov[j] += etacon[k]*Fcov[j][k];

  DLOOPA Bcov[j]=0;
  DLOOP Bcov[j] += etacon[k]*Mcov[k][j];

  Bsq=0.0;
  SLOOPA Bsq+=-alpha*prffde[B1+j-1]*Bcov[j];

  vcon[1]=-betai[1]+alpha*alpha*(Ecov[2]*Bcov[3]-Ecov[3]*Bcov[2])/(geom->g*Bsq);
  vcon[2]=-betai[2]+alpha*alpha*(Ecov[3]*Bcov[1]-Ecov[1]*Bcov[3])/(geom->g*Bsq);
  vcon[3]=-betai[3]+alpha*alpha*(Ecov[1]*Bcov[2]-Ecov[2]*Bcov[1])/(geom->g*Bsq);
#elif(0)
  DLOOPA Econ[j]=0;
  DLOOP Econ[j] += etacov[k]*Fcon[j][k];
  lower(Econ,geom,Ecov);

  DLOOPA Bcon[j]=0;
  DLOOP Bcon[j] += etacov[k]*Mcon[k][j];
  lower(Bcon,geom,Bcov);

  Bsq=0.0;
  DLOOPA Bsq+=Bcon[j]*Bcov[j]; // sign of Mcov_{jt} consistent with used by Mcov in vcon below.  That's all that matters.

  vcon[1]=-betai[1]+alpha*alpha*(Ecov[2]*Bcov[3]-Ecov[3]*Bcov[2])/(geom->g*Bsq);
  vcon[2]=-betai[2]+alpha*alpha*(Ecov[3]*Bcov[1]-Ecov[1]*Bcov[3])/(geom->g*Bsq);
  vcon[3]=-betai[3]+alpha*alpha*(Ecov[1]*Bcov[2]-Ecov[2]*Bcov[1])/(geom->g*Bsq);
#elif(0)
  vcon[1]=(Fcov[TT][2]*Mcov[3][TT]-Fcov[TT][3]*Mcov[2][TT])/(geom->g*Bsq);
  vcon[2]=(Fcov[TT][3]*Mcov[1][TT]-Fcov[TT][1]*Mcov[3][TT])/(geom->g*Bsq);
  vcon[3]=(Fcov[TT][1]*Mcov[2][TT]-Fcov[TT][2]*Mcov[1][TT])/(geom->g*Bsq);
#elif(0) 

  // seems to kinda work

  Bsq=0.0;
  SLOOPA Bsq+=-prffde[B1+j-1]*prffde[B1+j-1];

  vcon[1]=(Fcov[TT][2]*Mcon[3][TT]-Fcov[TT][3]*Mcon[2][TT])/(geom->g*Bsq);
  vcon[2]=(Fcov[TT][3]*Mcon[1][TT]-Fcov[TT][1]*Mcon[3][TT])/(geom->g*Bsq);
  vcon[3]=(Fcov[TT][1]*Mcon[2][TT]-Fcov[TT][2]*Mcon[1][TT])/(geom->g*Bsq);
#elif(0)
  SLOOPA vcon[j]=0;
  SLOOP vcon[j]+=((-Fcov[TT][k]/geom->g)*geom->g*Fcon[j][k])/(Bsq);
#elif(0)
  SLOOPA vcon[j]=-Ti0[j]/Bsq;
#endif

  // convert from v^i -> u^\mu
  prim[U1]=vcon[1];
  prim[U2]=vcon[2];
  prim[U3]=vcon[3];

  //ucon_calc_3vel(prim,geom,ucon);
  pr2ucon(VEL3,prim,geom,ucon); // if u^t can't be found, then FFDE breaks down (i.e. dF^2>0 is ok) 

  // convert from u^\mu -> true WHICHVEL primitive
  ucon2pr(WHICHVEL,ucon,geom,pr); // fills pr[U1->U3]

  // now assign field primitives
  pr[B1]=U[B1];
  pr[B2]=U[B2];
  pr[B3]=U[B3];
  

  // done.

  return(0);




}

/*
 * Convert conserved quantities back to primitives
 */

// convert T^t_i & M^{it} -> F^{it}
void UtoFit(FTYPE *U, struct of_geom *geom, FTYPE *Fit)
{
      int k ;
      int j;

      FTYPE Mt[NDIM] ;
      FTYPE Bsq ;
      FTYPE T0[NDIM] ;   /* stress tensor terms T^t_i */
      FTYPE Mx[NDIM] ;   /* mixed Maxwell components M^t_i */


      Mt[TT]=0.0;
      SLOOPA Mt[j] = -U[B1+j-1] ;  /* M^{ti} , where U[Bi]=B^i=M^{it}*/
      // T0's time term not needed
      SLOOPA T0[j] = U[U1+j-1] ; /* T^t_i terms (t up, i down) */

      /* lower second index of Maxwell to get M^t_i */
      lower(Mt,geom,Mx);

      // Bsq = \dF^t_\lambda \dF^{t\lambda} = Mx[1]*Mt1 + Mx[2]*Mt2 + Mx[3]*Mt3 ;
      Bsq=0.0;
      DLOOPA Bsq+=Mx[j]*Mt[j];
      
      /* calculate F^it
       *
       * F^it = Mx[j] T[k] / (Mx[j] Mtj), where i,j,k run
       * cyclic 1,2,3.
       *
       *
       * if Q^{\mu\nu} = \ep^{\mu\nu\lambda\delta} \dF^t_\lambda T^t_\delta
       *
       * then F^{it} = Q^{it} / (\dF^t_\lambda \dF^{t\lambda} ) = Q^{it} / Bsq
       *
       * where Q^{it} = \ep^{i t j k} \dF^t_j T^t_k = -(1/\detg)[itjk] \dF^t_j T^t_k
       *              = (1/\detg) [ijk] \dF^t_j T^t_k
       *       Q^{1t} = (1/\detg) (\dF^t_2 T^t_3 - \dF^t_3 T^t_2)
       *       Q^{2t} = (1/\detg) (\dF^t_3 T^t_1 - \dF^t_1 T^t_3)
       *       Q^{3t} = (1/\detg) (\dF^t_1 T^t_2 - \dF^t_2 T^t_1)
       */
      
      /* F^{it} */
      Fit[0] = 0;
      Fit[1] = (Mx[2]*T0[3] - Mx[3]*T0[2])/(Bsq*geom->g) ;
      Fit[2] = (Mx[3]*T0[1] - Mx[1]*T0[3])/(Bsq*geom->g) ;
      Fit[3] = (Mx[1]*T0[2] - Mx[2]*T0[1])/(Bsq*geom->g) ;


}      


//  realBsqOalphasq=BsqOalphasq; // default
void U2BBsq(FTYPE *U, struct of_geom *geom, FTYPE *BconOalpha, FTYPE *BcovOalpha, FTYPE *BsqOalphasq)
{
  int j;
  // Mx=M^t_\mu = g_{\mu\nu} M^{t\nu} = g_{\mu\nu} B^\nu/(-alpha)
  // remove alpha since cancels below
  //SLOOPA Bcov[j]=-alpha*Mx[j];
  //  SLOOPA BcovOalpha[j]=-Mx[j]; // Mx = \dF^t_j

  // check B^2-E^2
  BconOalpha[TT]=0.0;
  // Bcon without alpha
  SLOOPA BconOalpha[j] = U[B1+j-1] ;  /* Bcon=-M^{ti} , where U[Bi]=B^i=M^{it}*/

  lower(BconOalpha,geom,BcovOalpha);

  // remove alpha^2 since cancels below
  //  realoBsqgeom=oBsqgeom/alphasq;
  *BsqOalphasq=0.0;
  SLOOPA *BsqOalphasq+=BcovOalpha[j]*BconOalpha[j];



}



void TtoEconOalpha(FTYPE *T0, FTYPE *BcovOalpha, FTYPE BsqOalphasq, struct of_geom *geom, FTYPE *EconOalpha)
{
      FTYPE oBsqOalphasqgeom;

      oBsqOalphasqgeom=1.0/(BsqOalphasq*geom->g); // 1/(\detg B^2/\alpha^2)
      
      /* calculate F^it
       *
       * F^it = Mx[j] T[k] / (Mx[j] Mtj), where i,j,k run
       * cyclic 1,2,3.
       *
       *
       * if Q^{\mu\nu} = \ep^{\mu\nu\lambda\delta} \dF^t_\lambda T^t_\delta
       *
       * then F^{it} = Q^{it} / (\dF^t_\lambda \dF^{t\lambda} ) = Q^{it} / Bsq
       *
       * where Q^{it} = \ep^{i t j k} \dF^t_j T^t_k = -(1/\detg)[itjk] \dF^t_j T^t_k
       *              = (1/\detg) [ijk] \dF^t_j T^t_k
       *       Q^{1t} = (1/\detg) (\dF^t_2 T^t_3 - \dF^t_3 T^t_2)
       *       Q^{2t} = (1/\detg) (\dF^t_3 T^t_1 - \dF^t_1 T^t_3)
       *       Q^{3t} = (1/\detg) (\dF^t_1 T^t_2 - \dF^t_2 T^t_1)
       */
      
      // E^i/\alpha = (\alpha^2/(\detg B^2)) ([ijk] (B_j/\alpha) T^t_k)
      EconOalpha[0] = 0;
      EconOalpha[1] = (BcovOalpha[2]*T0[3] - BcovOalpha[3]*T0[2])*(oBsqOalphasqgeom) ;
      EconOalpha[2] = (BcovOalpha[3]*T0[1] - BcovOalpha[1]*T0[3])*(oBsqOalphasqgeom) ;
      EconOalpha[3] = (BcovOalpha[1]*T0[2] - BcovOalpha[2]*T0[1])*(oBsqOalphasqgeom) ;


}      



void UtoFitMxBsqgeom(FTYPE *U, struct of_geom *geom, FTYPE *Fit, FTYPE*Mx, FTYPE*oBsqgeom)
{
      int k ;
      int j;

      FTYPE Mt[NDIM] ;
      FTYPE T0[NDIM] ;   /* stress tensor terms T^t_i */
      //      FTYPE Mx[NDIM] ;   /* mixed Maxwell components M^t_i */
      // Mx=M^t_\mu = g_{\mu\nu} M^{t\nu} = g_{\mu\nu} B^\nu/(-alpha)
      // -> B_\mu = -alpha * M^t_\mu
      // Real Bsq = Bsqgeom*alpha^2
      FTYPE Bsqgeom;

      Mt[TT]=0.0;
      SLOOPA Mt[j] = -U[B1+j-1] ;  /* M^{ti} , where U[Bi]=B^i=M^{it}*/
      // T0's time term not needed
      SLOOPA T0[j] = U[U1+j-1] ; /* T^t_i terms (t up, i down) */

      /* lower second index of Maxwell to get M^t_i */
      lower(Mt,geom,Mx);

      // Bsq = \dF^t_\lambda \dF^{t\lambda} = Mx[1]*Mt1 + Mx[2]*Mt2 + Mx[3]*Mt3 ;
      // RealBsq = Bsq*alpha^2
      Bsqgeom=0.0;
      SLOOPA Bsqgeom+=Mx[j]*Mt[j];
      Bsqgeom*=geom->g;
      *oBsqgeom=1.0/Bsqgeom;
      
      /* calculate F^it
       *
       * F^it = Mx[j] T[k] / (Mx[j] Mtj), where i,j,k run
       * cyclic 1,2,3.
       *
       *
       * if Q^{\mu\nu} = \ep^{\mu\nu\lambda\delta} \dF^t_\lambda T^t_\delta
       *
       * then F^{it} = Q^{it} / (\dF^t_\lambda \dF^{t\lambda} ) = Q^{it} / Bsq
       *
       * where Q^{it} = \ep^{i t j k} \dF^t_j T^t_k = -(1/\detg)[itjk] \dF^t_j T^t_k
       *              = (1/\detg) [ijk] \dF^t_j T^t_k
       *       Q^{1t} = (1/\detg) (\dF^t_2 T^t_3 - \dF^t_3 T^t_2)
       *       Q^{2t} = (1/\detg) (\dF^t_3 T^t_1 - \dF^t_1 T^t_3)
       *       Q^{3t} = (1/\detg) (\dF^t_1 T^t_2 - \dF^t_2 T^t_1)
       */
      
      /* F^{it} */
      Fit[0] = 0;
      Fit[1] = (Mx[2]*T0[3] - Mx[3]*T0[2])*(*oBsqgeom) ;
      Fit[2] = (Mx[3]*T0[1] - Mx[1]*T0[3])*(*oBsqgeom) ;
      Fit[3] = (Mx[1]*T0[2] - Mx[2]*T0[1])*(*oBsqgeom) ;


}      

/* 
 * convert primitives to conserved quantities 
 */

void primtoU_ffde(FTYPE pr[NPR], struct of_geom *geom, FTYPE U[NPR])
{

      int k ;
      FTYPE F1t,F2t,F3t,Mt1,Mt2,Mt3 ;

      /* load prims */
      F1t = pr[U1] ;   /* F^it -- E components */
      F2t = pr[U2] ;
      F3t = pr[U3] ;
      Mt1 = pr[B1] ;   /* M^ti -- B components */
      Mt2 = pr[B2] ;
      Mt3 = pr[B3] ;

      /* M^ti -- just copy!*/
      U[0] = Mt1 ;
      U[1] = Mt2 ;
      U[2] = Mt3 ;

      /* T^t_i terms( t up i down):
       *
       * T^t_i = - sqrt(-g) * [ijk] * M^tj * F^kt
       *
       */

      U[3] = - geom->g * (Mt2*F3t - Mt3*F2t) ;
      U[4] = - geom->g * (Mt3*F1t - Mt1*F3t) ;
      U[5] = - geom->g * (Mt1*F2t - Mt2*F1t) ;

      /* mult by geometric factor */
      PLOOP U[k] *= geom->g ;
}


/*
 * load stress tensor and contravariant Maxwell in structure state
 */

void get_state_ffde(FTYPE pr[NPR], struct of_geom *geom, struct of_state *q)
{
      /* get contravariant Maxwell and stress tensor */
  //      Max_con(pr, geom ,q->Mcon) ;
  //      stresstensor(q->Mcon, geom, q->T) ;
}

// should be able to compare this with mhd_calc
// i.e.
// get_geometry()
// get_state
// mhd_calc(pr,dir,geom,q,mhd) where T^{dir}_j is returned
//
void ffdestresstensor(FTYPE Mcon[NDIM][NDIM], struct of_geom *geom, FTYPE T[NDIM][NDIM])
{
      int i,j,k ;
      void lower_A(FTYPE Acon[NDIM][NDIM], struct of_geom *geom, FTYPE Acov[NDIM][NDIM]);


      FTYPE Mcov[NDIM][NDIM] ; /* covariant Maxwell */
      FTYPE Msq ;

      /* get covariant Maxwell from contravariant */
      lower_A(Mcon, geom, Mcov) ;

      /* out with the old */
      DLOOP  T[j][k] = 0. ;

      /* in with the new:
       *
       * a, b, c, d run 0 to 4:
       *
       * T^a_b = M^ac * M_bc - 1/4 KroneckerDelta[a,b] (M_cd M^cd)
       *
       */

      Msq = 0. ;
      DLOOP Msq += Mcov[j][k]*Mcon[j][k] ;

      DLOOP {
	for(i = 0; i < NDIM; i++) {
	  T[j][k] += Mcon[j][i]*Mcov[k][i] ;
	}
      }
      DLOOPA T[j][j] -= 0.25 * Msq ;
}

/*
 * lower both indices of an anti-symmetric  tensor
 *  -- go from contravariant to covariant
 */

/*
 * calculate the contravariant Maxwell
 */
void Max_con(FTYPE prffde[NPR], struct of_geom *geom, FTYPE Mcon[NDIM][NDIM])
{
      int i,j,k ;

      FTYPE Fit[NDIM];
      FTYPE Bcon[NDIM];
      FTYPE Econ[NDIM], Ecov[NDIM] ;          /* E four-vector */
      FTYPE etacon[NDIM],etacov[NDIM];
      FTYPE alpha;

      // lapse
      // define \eta_\alpha
      // assume always has 0 value for space components
      alpha = 1./sqrt(- geom->gcon[0][0]);  

      etacov[TT]=-alpha; // any constant will work.
      SLOOPA etacov[j]=0.0; // must be 0

      // shift
      /* beta */
      //      SLOOPA beta[j] = geom->gcon[0][j]*alpha*alpha ;
      // define \eta^\beta
      raise(etacov,geom,etacon);

      Fit[0] = 0.0; // Ftt=0
      Fit[1] = prffde[U1] ; /* F^{it} = -E^i/\detg in Jon notation */
      Fit[2] = prffde[U2] ;
      Fit[3] = prffde[U3] ;

      // B^\nu = \eta_\mu \dF^{\mu\nu}
      Bcon[0] = 0.0; // assumes \eta_\alpha={CONST,0,0,0}
      Bcon[1] = etacov[TT]*prffde[B1] ; /* M^{ti} =-B^i in Jon notation */
      Bcon[2] = etacov[TT]*prffde[B2] ;
      Bcon[3] = etacov[TT]*prffde[B3] ;

      /* create E^{\mu} = \eta_\beta F^{\alpha\beta} vector */

      DLOOPA Econ[j]=etacov[TT]*Fit[j];
      lower(Econ, geom, Ecov) ;

      /* diagonal */
      Mcon[0][0] = Mcon[1][1] = Mcon[2][2] = Mcon[3][3] = 0. ;

      /* out with the old */
      //DLOOP  Mcon[j][k] = 0. ;

      /* M^ij terms:
       *
       * M^ij = -beta^i M^tj + beta^j M^ti
       * + alpha * (1/sqrt(-g)) * Ecov[k]
       *
       * use primitive variables for M's
       *
       */

      // again assume \eta_\delta={CONST,0,0,0}
      //
      // \dF^{\alpha\beta} = -eta^\alpha B^\beta + \eta^\beta B^\alpha - E_\gamma \eta_\delta e^{\alpha\beta\gamma\delta}
      //                   = -eta^\alpha B^\beta + \eta^\beta B^\alpha + E_i \eta_t e^{t \alpha\beta i} // 3 sign switches from t
      // \dF^{jk}          = -eta^j B^k + \eta^k B^j - (1/\detg) E_i \eta_t [jki] // 1 sign switch from e^
      // \dF^{jk}          = -eta^j B^k + \eta^k B^j - (\eta_t/\detg) (E_i [jki])
      j=1;k=2;i=3; Mcon[j][k] = -etacon[j]*Bcon[k] + etacon[k]*Bcon[j] - (etacov[TT]/geom->g)*Ecov[i];
      j=2;k=3;i=1; Mcon[j][k] = -etacon[j]*Bcon[k] + etacon[k]*Bcon[j] - (etacov[TT]/geom->g)*Ecov[i];
      j=3;k=1;i=2; Mcon[j][k] = -etacon[j]*Bcon[k] + etacon[k]*Bcon[j] - (etacov[TT]/geom->g)*Ecov[i];


 
      /* copy remaining spacial terms */
   
      Mcon[2][1] = -Mcon[1][2] ;
      Mcon[3][2] = -Mcon[2][3] ;
      Mcon[1][3] = -Mcon[3][1] ;

      /* time terms - easy */
      SLOOPA{
	Mcon[TT][j]=Bcon[j]/etacov[TT];
	Mcon[j][TT]=-Mcon[TT][j];
      }
}


/*
 * calculate the contravariant Maxwell
 */
void Max_con_old(FTYPE prffde[NPR], struct of_geom *geom, FTYPE Mcon[NDIM][NDIM])
{
      int j ;

      FTYPE F1t, F2t, F3t, Mt1, Mt2, Mt3 ;    /* prims */
      FTYPE Econ[NDIM], Ecov[NDIM] ;          /* E four-vector */
      FTYPE beta[NDIM] ;                      /* shift  3-vector  */
      FTYPE alpha;


      /* lapse */
      alpha = 1./sqrt(- geom->gcon[0][0]);  

      /* beta */ 
      // GODMARK: wrongly used below.  Should contract Mt with \eta^\mu = 1/\alpha (1,-\beta^i), so missing 1/\alpha
      // but then left off -alpha on Mt terms in Mcon below, so Mcon below ends up being correct (including signs)
      SLOOPA beta[j] = geom->gcon[0][j]*alpha*alpha ;

      F1t = prffde[U1] ; /* F^{it} = -E^i/\detg in Jon notation */
      F2t = prffde[U2] ;
      F3t = prffde[U3] ;
      Mt1 = prffde[B1] ; /* M^{ti} =-B^i in Jon notation */
      Mt2 = prffde[B2] ;
      Mt3 = prffde[B3] ;

      /* create E^{\mu} = \eta_\beta F^{\alpha\beta} vector */
      Econ[0] = 0. ;
      Econ[1] = - alpha * F1t ;
      Econ[2] = - alpha * F2t ;
      Econ[3] = - alpha * F3t ;

      lower(Econ, geom, Ecov) ;

      /* diagonal */
      Mcon[0][0] = Mcon[1][1] = Mcon[2][2] = Mcon[3][3] = 0. ;

      /* out with the old */
      //DLOOP  Mcon[j][k] = 0. ;

      /* M^ij terms:
       *
       * M^ij = -beta^i M^tj + beta^j M^ti
       * + alpha * (1/sqrt(-g)) * Ecov[k]
       *
       * use primitive variables for M's
       *
       */

      // different sign! for E part (GODMARK) (see above comments -- just crazy mixing of terms leads to same results)
      Mcon[1][2] = (-beta[1] * Mt2 + beta[2] * Mt1 + alpha/(geom->g) * Ecov[3]) ;
      Mcon[2][3] = (-beta[2] * Mt3 + beta[3] * Mt2 + alpha/(geom->g) * Ecov[1]) ;
      Mcon[3][1] = (-beta[3] * Mt1 + beta[1] * Mt3 + alpha/(geom->g) * Ecov[2]) ;
	 
      /* copy remaining spacial terms */
   
      Mcon[2][1] = -Mcon[1][2] ;
      Mcon[3][2] = -Mcon[2][3] ;
      Mcon[1][3] = -Mcon[3][1] ;

      /* time terms - easy */
      Mcon[0][1] = Mt1 ;
      Mcon[0][2] = Mt2 ;
      Mcon[0][3] = Mt3 ;
      Mcon[1][0] = - Mt1 ;
      Mcon[2][0] = - Mt2 ;
      Mcon[3][0] = - Mt3 ;
}

void raise_A(FTYPE Acov[NDIM][NDIM], struct of_geom *geom, FTYPE Acon[NDIM][NDIM])
{
        int j,k ;

	/* out with the old */
	DLOOP  Acon[j][k] = 0. ;

	/* in with the new */
	DLOOP {
       	       Acon[0][1] += (geom->gcon[0][j])*(geom->gcon[1][k])*Acov[j][k] ;
       	       Acon[0][2] += (geom->gcon[0][j])*(geom->gcon[2][k])*Acov[j][k] ;
       	       Acon[0][3] += (geom->gcon[0][j])*(geom->gcon[3][k])*Acov[j][k] ;
	       Acon[1][2] += (geom->gcon[1][j])*(geom->gcon[2][k])*Acov[j][k] ;
	       Acon[2][3] += (geom->gcon[2][j])*(geom->gcon[3][k])*Acov[j][k] ;
	       Acon[3][1] += (geom->gcon[3][j])*(geom->gcon[1][k])*Acov[j][k] ;
	      }

	/*
	 * diagonal is already set to zero,
	 * just copy the rest
	 */

	Acon[1][0] = - Acon[0][1] ;
	Acon[2][0] = - Acon[0][2] ;
	Acon[3][0] = - Acon[0][3] ;
	Acon[2][1] = - Acon[1][2] ;
	Acon[3][2] = - Acon[2][3] ;
	Acon[1][3] = - Acon[3][1] ;
}

void lower_A(FTYPE Acon[NDIM][NDIM], struct of_geom *geom, FTYPE Acov[NDIM][NDIM])
{
        int j,k ;

	/* out with the old */
	DLOOP  Acov[j][k] = 0. ;

	/* in with the new */
	DLOOP {
       	       Acov[0][1] += (geom->gcov[0][j])*(geom->gcov[1][k])*Acon[j][k] ;
       	       Acov[0][2] += (geom->gcov[0][j])*(geom->gcov[2][k])*Acon[j][k] ;
       	       Acov[0][3] += (geom->gcov[0][j])*(geom->gcov[3][k])*Acon[j][k] ;
	       Acov[1][2] += (geom->gcov[1][j])*(geom->gcov[2][k])*Acon[j][k] ;
	       Acov[2][3] += (geom->gcov[2][j])*(geom->gcov[3][k])*Acon[j][k] ;
	       Acov[3][1] += (geom->gcov[3][j])*(geom->gcov[1][k])*Acon[j][k] ;
	      }

	/*
	 * diagonal is already set to zero,
	 * just copy the rest
	 */

	Acov[1][0] = - Acov[0][1] ;
	Acov[2][0] = - Acov[0][2] ;
	Acov[3][0] = - Acov[0][3] ;
	Acov[2][1] = - Acov[1][2] ;
	Acov[3][2] = - Acov[2][3] ;
	Acov[1][3] = - Acov[3][1] ;
}



// Maxwell to Faraday
// which=0 : Mcon -> Fcov (for clean Mcon, Fcov has \detg)
// which=1 : Mcov -> Fcon (for clean Mcov)

// which=2 : Fcon -> Mcov
// which=3 : Fcov -> Mcon

// copies faraday_calc() in phys.c
void MtoF(int which, FTYPE IN[NDIM][NDIM],struct of_geom *geom, FTYPE OUT[NDIM][NDIM])
{
  int nu,mu,kappa,lambda;
  FTYPE prefactor;
  int whichlc;

  if((which==0)||(which==1)){
    prefactor=-0.5;
    if(which==0) whichlc=0;
    if(which==1) whichlc=1;
  }
  if((which==2)||(which==3)){
    prefactor=0.5;
    if(which==2) whichlc=0;
    if(which==3) whichlc=1;
  }

  for(nu=0;nu<NDIM;nu++) for(mu=0;mu<NDIM;mu++){
    OUT[mu][nu]=0.0;
    for(kappa=0;kappa<NDIM;kappa++) for(lambda=0;lambda<NDIM;lambda++){
      OUT[mu][nu]+=prefactor*lc4(whichlc,geom->g,mu,nu,kappa,lambda)*IN[kappa][lambda];
    }
  }


}


// see also phys.ffde.debug.c

// assumes geometry exists
void testffdeinversion(void)
{
  int i,j,k;
  struct of_geom geom;
  FTYPE prout[NPR],prin[NPR];
  FTYPE prout2[NPR];
  FTYPE U[NPR];
  struct of_state q;
  int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr)  ;
  FTYPE Bu[NDIM],uu[NDIM];
  int itest,jtest;
  FTYPE diff;
  FTYPE largestu[NDIM],largestB[NDIM];
  FTYPE largesterror;
  int qq;

#if(0)
  itest=50;
  jtest=32;
  PLOOP{
    prin[k]=p[itest][jtest][k];
    prin[U1]=prin[U3]*0.1;
    prin[U2]=prin[U2]*0.1;
    prin[B3]=prin[B1];
  }

  get_geometry(itest,jtest,CENT,&geom);
  get_state(prin,&geom,&q);

  primtoU(UNOTHING,prin,&q,&geom,U);
  Utoprim_ffde(U,&geom,prout); // no need for initial guess since analytic inversion

  PLOOP prin[k]=prout[k];

  primtoU(UNOTHING,prin,&q,&geom,U);
  Utoprim_ffde(U,&geom,prout); // no need for initial guess since analytic inversion

  // just compare pr in and pr out.
  PLOOP dualfprintf(fail_file,"prold[%d]=%21.15g  prnew[%d]=%21.15g :: %21.15g\n",k,prin[k],k,prout[k],(prin[k]-prout[k])/prin[k]);


  i=48;
  j=32;
  PLOOP{
    prin[k]=p[i][j][k];
    prout[k]=0.0;
  }

  get_geometry(i,j,CENT,&geom);
  get_state(prin,&geom,&q);
  primtoU(UNOTHING,prin,&q,&geom,U);
  Utoprim_ffde(U,&geom,prout); // no need for initial guess since analytic inversion

  // just compare pr in and pr out.
  PLOOP dualfprintf(fail_file,"prold[%d]=%21.15g  prnew[%d]=%21.15g :: %21.15g\n",k,prin[k],k,prout[k],(prin[k]-prout[k])/prin[k]);



  // loop through to test consistency
  PLOOP{
    prin[k]=prout[k];
    prout[k]=0.0;
  }

  get_geometry(i,j,CENT,&geom);
  get_state(prin,&geom,&q);
  primtoU(UNOTHING,prin,&q,&geom,U);

  PLOOP dualfprintf(fail_file,"U[%d]=%21.15g\n",k,U[k]);
  Utoprim_ffde(U,&geom,prout); // no need for initial guess since analytic inversion

  // just compare pr in and pr out.
  PLOOP dualfprintf(fail_file,"prold[%d]=%21.15g  prnew[%d]=%21.15g :: %21.15g\n",k,prin[k],k,prout[k],(prin[k]-prout[k])/prin[k]);

#endif
  

  // now do loop over range of relative 4-velocity and field for \rho=u=0

  dualfprintf(fail_file,"realtest\n");




#if(0)
  largestu[NDIM]=largestB[NDIM]=0.0;
  largesterror=0;
  i=3;
  j=32;
  for(uu[1]=-100.0;uu[1]<=100.0;uu[1]+=10)
    for(uu[2]=-100.0;uu[2]<=100.0;uu[2]+=10)
      for(uu[3]=-100.0;uu[3]<=010.0;uu[3]+=10)

	for(Bu[1]=-10.0;Bu[1]<=10.0;Bu[1]+=1)
	  for(Bu[2]=-10.0;Bu[2]<=10.0;Bu[2]+=1)
	    for(Bu[3]=-10.0;Bu[3]<=10.0;Bu[3]+=1){

	      if((Bu[1]==0)&&(Bu[2]==0)&&(Bu[3]==0)) continue;

	      prin[RHO]=prin[UU]=0;

	      prin[U1]=uu[1];
	      prin[U2]=uu[2];
	      prin[U3]=uu[3];

	      prin[B1]=Bu[1];
	      prin[B2]=Bu[2];
	      prin[B3]=Bu[3];

	      get_geometry(i,j,CENT,&geom);

	      if(get_state(prin,&geom,&q)>=1) dualfprintf(fail_file,"getstate failure in realtest\n");
	      if(primtoU(UNOTHING,prin,&q,&geom,U)>=1) dualfprintf(fail_file,"primtoU failure in realtest\n");
	      
	      Utoprim_ffde(U,&geom,prout); // no need for initial guess since analytic inversion

	      PLOOP prin[k]=prout[k];

	      // only clean solution to test
	      get_state(prin,&geom,&q);
	      primtoU(UNOTHING,prin,&q,&geom,U);
	      
	      Utoprim_ffde(U,&geom,prout); // no need for initial guess since analytic inversion

	      // just compare pr in and pr out.
	      for(k=U1;k<=B3;k++){
		diff=fabs(prin[k]-prout[k])/(fabs(prin[k])+fabs(prout[k]));
		if((fabs(prin[k])>1E-10)&&(diff>1E-5)){
		  dualfprintf(fail_file,"new prold[%d]=%21.15g  prnew[%d]=%21.15g :: %21.15g\n",k,prin[k],k,prout[k],diff); 
		  SLOOPA dualfprintf(fail_file,"new realtest uu[%d]=%21.15g\n",j,uu[j]);
		  SLOOPA dualfprintf(fail_file,"new realtest Bu[%d]=%21.15g\n",j,Bu[j]);
		  if(diff>largesterror){
		    largesterror=diff;
		    for(qq=1;qq<=3;qq++){
		      largestu[qq]=uu[qq];
		      largestB[qq]=Bu[qq];
		    }
		  }
		}
	      }

	      Utoprim_ffde_oldnew(U,&geom,prout); // no need for initial guess since analytic inversion

	      // just compare pr in and pr out.
	      for(k=U1;k<=B3;k++){
		diff=fabs(prin[k]-prout[k])/(fabs(prin[k])+fabs(prout[k]));
		if((fabs(prin[k])>1E-10)&&(diff>1E-5)){
		  dualfprintf(fail_file,"old prold[%d]=%21.15g  prnew[%d]=%21.15g :: %21.15g\n",k,prin[k],k,prout[k],diff); 
		  SLOOPA dualfprintf(fail_file,"old realtest uu[%d]=%21.15g\n",j,uu[j]);
		  SLOOPA dualfprintf(fail_file,"old realtest Bu[%d]=%21.15g\n",j,Bu[j]);
		}
	      }


	      //	      myexit(0);
	    }
  dualfprintf(fail_file,"largesterror=%g\n",largesterror);
  for(qq=1;qq<=3;qq++) dualfprintf(fail_file,"largestu[%d]=%21.15g largestB[%d]=%21.15g\n",qq,largestu[qq],qq,largestB[qq]);


#endif


#if(1)
  i=3;
  j=32;
  //  i=0;
  //  j=40;
  //  uu[1]=uu[2]=uu[3]=-1.9999E3; //yes
  //  uu[1]=uu[2]=uu[3]=-2.4E3; // no
  //  uu[1]=uu[2]=uu[3]=-2.28E3; // top? no?
  //uu[1]=uu[2]=uu[3]=-2.3E3; // works! (with old formulation)
  //  uu[1]=uu[2]=uu[3]=-2.30000099E3; // odd exact machine precision answer
  //uu[1]=uu[2]=uu[3]=-2.300000999E3; // doesn't work.

  // new formulation
  //  uu[1]=uu[2]=uu[3]=-3E4; //yes
  /*
  // I get 5E-5 error for the below with uu2=-1E5
  // 8E-3 error for uu2=-1E6
  // -0.35 error for uu2=-1E7
  uu[1]=0;
  uu[2]=-1E5;
  uu[3]=0;

  Bu[1]=0;
  Bu[2]=0;
  Bu[3]=-1;

  // even uu2=-1E10 works, but uu1=-1E10 doesn't.
  */
  uu[1]=-2270.477;
  uu[2]=-5.407722;
  uu[3]=0;

  Bu[1]=.5;
  Bu[2]=-1;
  Bu[3]=-.25;

  uu[1]=-2.477;
  uu[2]=-1.407722;
  uu[3]=0;

  Bu[1]=15;
  Bu[2]=-11;
  Bu[3]=-13;

  // largest error version
  uu[1]=0;
  uu[2]=0;
  uu[3]=-40;

  Bu[1]=0;
  Bu[2]=-6;
  Bu[3]=0;

  // largest error version
  uu[1]=-3.17661921216517e-08;
  uu[2]=2.07598472203926e-12;
  uu[3]=-13.4160167167174;

  Bu[1]=0;
  Bu[2]=-6;
  Bu[3]=0;
	      prin[RHO]=prin[UU]=0;

	      prin[U1]=uu[1];
	      prin[U2]=uu[2];
	      prin[U3]=uu[3];

	      prin[B1]=Bu[1];
	      prin[B2]=Bu[2];
	      prin[B3]=Bu[3];

	      SLOOPA dualfprintf(fail_file,"realtest uu[%d]=%21.15g\n",j,uu[j]);
	      SLOOPA dualfprintf(fail_file,"realtest Bu[%d]=%21.15g\n",j,Bu[j]);

	      get_geometry(i,j,CENT,&geom);

	      if(get_state(prin,&geom,&q)>=1) dualfprintf(fail_file,"getstate failure in realtest\n");
	      for(k=0;k<NDIM;k++) fprintf(fail_file,"1 uu[%d]=%21.15g\n",k,q.ucon[k]);
	      if(primtoU(UNOTHING,prin,&q,&geom,U)>=1) dualfprintf(fail_file,"primtoU failure in realtest\n");
	      
	      Utoprim_ffde(U,&geom,prout); // no need for initial guess since analytic inversion
	      for(k=U1;k<=B3;k++) fprintf(fail_file,"new: prold[%d]=%21.15g  prnew[%d]=%21.15g :: %21.15g\n",k,prin[k],k,prout[k],(prin[k]-prout[k])/prin[k]); 

	      PLOOP prin[k]=prout[k];

	      Utoprim_ffde_oldnew(U,&geom,prout); // no need for initial guess since analytic inversion
	      for(k=U1;k<=B3;k++) fprintf(fail_file,"old: prold[%d]=%21.15g  prnew[%d]=%21.15g :: %21.15g\n",k,prin[k],k,prout[k],(prin[k]-prout[k])/prin[k]); 


	      // only clean solution to test
	      get_state(prin,&geom,&q);
	      for(k=0;k<NDIM;k++) fprintf(fail_file,"2 uu[%d]=%21.15g\n",k,q.ucon[k]);
	      primtoU(UNOTHING,prin,&q,&geom,U);
	      
	      Utoprim_ffde(U,&geom,prout); // no need for initial guess since analytic inversion

	      // just compare pr in and pr out.
	      for(k=U1;k<=B3;k++) fprintf(fail_file,"new: prold[%d]=%21.15g  prnew[%d]=%21.15g :: %21.15g\n",k,prin[k],k,prout[k],(prin[k]-prout[k])/prin[k]); 
	      fflush(fail_file);

	      Utoprim_ffde_oldnew(U,&geom,prout); // no need for initial guess since analytic inversion
	      for(k=U1;k<=B3;k++) fprintf(fail_file,"old: prold[%d]=%21.15g  prnew[%d]=%21.15g :: %21.15g\n",k,prin[k],k,prout[k],(prin[k]-prout[k])/prin[k]); 


	      //	      myexit(0);
#endif


  myexit(0);




}

// filter out velocity along field line
void filterffde(int i, int j, FTYPE *pr)
{
  int k;
  struct of_geom geom;
  FTYPE prout[NPR],prin[NPR];
  FTYPE U[NPR];
  struct of_state q;
  int Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr)  ;


  // now filter velocity
  get_geometry(i,j,CENT,&geom);
  get_state(pr,&geom,&q);
  primtoU(UNOTHING,pr,&q,&geom,U);

  Utoprim_ffde(U,&geom,prout); // no need for initial guess since analytic inversion

  PLOOP pr[k]=prout[k];
  // kill densities
  pr[RHO]=pr[UU]=0.0;



}
