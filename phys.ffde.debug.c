
#include "decs.h"


/*

Jon notation:

// see grmhd-faraday_allforms.nb

F^{\mu\nu} = -1/\detg {{0,-Eu1,-Eu2,-Eu3},{Eu1,0,-Bd3,Bd2},{Eu2,Bd3,0,-Bd1},{Eu3,-Bd2,Bd1,0}}

\dF_{\mu\nu} = {{0,Bd1,Bd2,Bd3},{-Bd1,0,Eu3,-Eu2},{-Bd2,-Eu3,0,Eu1},{-Bd3,Eu2,-Eu1,0}}


F^{ij} = (1/\detg) [ijk] B_k 

F_{ij}/\detg = B^k [ijk]


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



// assumes geometry-free U in standard form
void Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr)
{
  int i,j,k;
  FTYPE T[NDIM][NDIM];
  FTYPE Ti0[NDIM] ;   /* stress tensor terms T^t_i */
  FTYPE Fit[NDIM],prffde[NPR];
  FTYPE Mcon[NDIM][NDIM] ; /* contravariant Maxwell */
  FTYPE Mcov[NDIM][NDIM] ; /* covariant Maxwell */
  FTYPE ucon[NDIM];
  FTYPE prim[NPR];
  FTYPE Fcov[NDIM][NDIM];
  FTYPE Fcon[NDIM][NDIM];
  FTYPE Bsq;
  FTYPE vcon[NDIM];
  FTYPE Ed[NDIM];
  void UtoFit(FTYPE *U, struct of_geom *geom, FTYPE *Fit);
  void Max_con(FTYPE prffde[NPR], struct of_geom *geom, FTYPE Mcon[NDIM][NDIM]);
  void Max_con_old(FTYPE prffde[NPR], struct of_geom *geom, FTYPE Mcon[NDIM][NDIM]);
  void lower_A(FTYPE Acon[NDIM][NDIM], struct of_geom *geom, FTYPE Acov[NDIM][NDIM]);
  void MtoF(int which, FTYPE Max[NDIM][NDIM],struct of_geom *geom, FTYPE faraday[NDIM][NDIM]);
  void raise_A(FTYPE Acov[NDIM][NDIM], struct of_geom *geom, FTYPE Acon[NDIM][NDIM]);
  void ffdestresstensor(FTYPE Mcon[NDIM][NDIM], struct of_geom *geom, FTYPE T[NDIM][NDIM]);
  FTYPE ftemp;
  FTYPE Mt[NDIM] ;
  FTYPE Mx[NDIM] ;   /* mixed Maxwell components M^t_i */
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

  DLOOPA dualfprintf(fail_file,"Fit[%d]=%21.15g\n",j,Fit[j]);

  // get M^{\mu\nu}
  // correct order, sign, etc. for Max_con
  prffde[U1]=Fit[1]; 
  prffde[U2]=Fit[2];
  prffde[U3]=Fit[3];
  prffde[B1]=-U[B1]; // M^{ti} since  U(geomfree)=B^i = \dF^{it} = M^{it}
  prffde[B2]=-U[B2];
  prffde[B3]=-U[B3];
  
  // gcon
  dualfprintf(fail_file,"gcon={");
  for(j=0;j<NDIM;j++){
    dualfprintf(fail_file,"{");
    for(k=0;k<NDIM;k++){
      dualfprintf(fail_file,"%21.15g``25",geom->gcon[j][k]);
      if(k!=NDIM-1) dualfprintf(fail_file,",");
    }
    if(j!=NDIM-1) dualfprintf(fail_file,"},");
    else dualfprintf(fail_file,"}");
  }
  dualfprintf(fail_file,"}");


  Max_con(prffde, geom, Mcon);
  //  Max_con_old(prffde,geom,Mcon);


  DLOOP dualfprintf(fail_file,"Mcon[%d][%d]=%21.15g\n",j,k,Mcon[j][k]);

  // Mcon
  dualfprintf(fail_file,"Mconcode={");
  for(j=0;j<NDIM;j++){
    dualfprintf(fail_file,"{");
    for(k=0;k<NDIM;k++){
      dualfprintf(fail_file,"%21.15g``25",Mcon[j][k]);
      if(k!=NDIM-1) dualfprintf(fail_file,",");
    }
    if(j!=NDIM-1) dualfprintf(fail_file,"},");
    else dualfprintf(fail_file,"}");
  }
  dualfprintf(fail_file,"}");

   // for alternative formulation
  /* T^i_t terms (i up, t down) */
  ffdestresstensor(Mcon, geom, T);
  SLOOPA Ti0[j] = T[j][TT] ; // T^i_t

  DLOOP dualfprintf(fail_file,"Tud[%d][%d]=%21.15g\n",j,k,T[j][k]);

  // Tud
  dualfprintf(fail_file,"Tudcode={");
  for(j=0;j<NDIM;j++){
    dualfprintf(fail_file,"{");
    for(k=0;k<NDIM;k++){
      dualfprintf(fail_file,"%21.15g``25",T[j][k]);
      if(k!=NDIM-1) dualfprintf(fail_file,",");
    }
    if(j!=NDIM-1) dualfprintf(fail_file,"},");
    else dualfprintf(fail_file,"}");
  }
  dualfprintf(fail_file,"}");


  /* get covariant Maxwell from contravariant */
  lower_A(Mcon, geom, Mcov) ;

  DLOOP dualfprintf(fail_file,"Mcov[%d][%d]=%21.15g\n",j,k,Mcov[j][k]);

  // get F_{\mu\nu}
  MtoF(0,Mcon,geom,Fcov);

  // Fcov
  DLOOP dualfprintf(fail_file,"Fcov[%d][%d]=%21.15g\n",j,k,Fcov[j][k]);

  dualfprintf(fail_file,"Fcovcode={");
  for(j=0;j<NDIM;j++){
    dualfprintf(fail_file,"{");
    for(k=0;k<NDIM;k++){
      dualfprintf(fail_file,"%21.15g``25",Fcov[j][k]);
      if(k!=NDIM-1) dualfprintf(fail_file,",");
    }
    if(j!=NDIM-1) dualfprintf(fail_file,"},");
    else dualfprintf(fail_file,"}");
  }
  dualfprintf(fail_file,"}");


  raise_A(Fcov,geom,Fcon);

  DLOOP dualfprintf(fail_file,"Fcon[%d][%d]=%21.15g\n",j,k,Fcon[j][k]);

  dualfprintf(fail_file,"Fconcode={");
  for(j=0;j<NDIM;j++){
    dualfprintf(fail_file,"{");
    for(k=0;k<NDIM;k++){
      dualfprintf(fail_file,"%21.15g``25",Fcon[j][k]);
      if(k!=NDIM-1) dualfprintf(fail_file,",");
    }
    if(j!=NDIM-1) dualfprintf(fail_file,"},");
    else dualfprintf(fail_file,"}");
  }
  dualfprintf(fail_file,"}");


  dualfprintf(fail_file,"geom->g=%21.15g\n",geom->g);

  //  myexit(0);

  ///////////////////////////////////////////////////////////
  // now can get velocity
  //
  //
  Bsq=0.0;
  DLOOPA Bsq+=-Mcon[j][TT]*Mcov[j][TT]; // sign of Mcov_{jt} consistent with used by Mcov in vcon below.  That's all that matters.

  dualfprintf(fail_file,"Bsq=%21.15g\n",Bsq);

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
  DLOOPA Econ[j]=0;
  DLOOP Econ[j] += etacov[k]*Fcon[j][k];
  lower(Econ,geom,Ecov);

  DLOOPA Bcon[j]=0;
  DLOOP Bcon[j] += etacov[k]*Mcon[k][j];
  lower(Bcon,geom,Bcov);

  Bsq=0.0;
  DLOOPA Bsq+=Bcon[j]*Bcov[j]; // sign of Mcov_{jt} consistent with used by Mcov in vcon below.  That's all that matters.
  SLOOPA dualfprintf(fail_file,"beta[%d]=%21.15g\n",j,betai[j]);

  vcon[1]=-betai[1]+alpha*alpha*(Ecov[2]*Bcov[3]-Ecov[3]*Bcov[2])/(geom->g*Bsq);
  vcon[2]=-betai[2]+alpha*alpha*(Ecov[3]*Bcov[1]-Ecov[1]*Bcov[3])/(geom->g*Bsq);
  vcon[3]=-betai[3]+alpha*alpha*(Ecov[1]*Bcov[2]-Ecov[2]*Bcov[1])/(geom->g*Bsq);
#elif(0)
  vcon[1]=(Fcov[TT][2]*Mcov[3][TT]-Fcov[TT][3]*Mcov[2][TT])/(geom->g*Bsq);
  vcon[2]=(Fcov[TT][3]*Mcov[1][TT]-Fcov[TT][1]*Mcov[3][TT])/(geom->g*Bsq);
  vcon[3]=(Fcov[TT][1]*Mcov[2][TT]-Fcov[TT][2]*Mcov[1][TT])/(geom->g*Bsq);
#elif(0)
  //  Mt[TT]=0.0;
  //  SLOOPA Mt[j] = -U[B1+j-1] ;  /* M^{ti} , where U[Bi]=B^i=M^{it}*/
  /* lower second index of Maxwell to get Mx=M^t_i */
  //  lower(Mt,geom,Mx);
  
  // Bsq = \dF^t_\lambda \dF^{t\lambda} = Mx[1]*Mt1 + Mx[2]*Mt2 + Mx[3]*Mt3 ;
  //  Bsq=0.0;
  //  DLOOPA Bsq+=-Mt[j]*Mt[j];

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

  SLOOPA dualfprintf(fail_file,"vcon[%d]=%21.15g\n",j,vcon[j]);

  // convert from v^i -> u^\mu
  prim[U1]=vcon[1];
  prim[U2]=vcon[2];
  prim[U3]=vcon[3];

  ucon_calc_3vel(prim,geom,ucon);
  DLOOPA dualfprintf(fail_file,"ucon_calc[%d]=%21.15g\n",j,ucon[j]);

  pr2ucon(VEL3,prim,geom,ucon); // if u^t can't be found, then FFDE breaks down (i.e. dF^2>0 is ok) 

  DLOOPA dualfprintf(fail_file,"pr2ucon[%d]=%21.15g\n",j,ucon[j]);

  // convert from u^\mu -> true WHICHVEL primitive
  ucon2pr(WHICHVEL,ucon,geom,pr); // fills pr[U1->U3]

  PLOOP dualfprintf(fail_file,"pr[%d]=%21.15g\n",k,pr[k]);

  // now assign field primitives
  pr[B1]=U[B1];
  pr[B2]=U[B2];
  pr[B3]=U[B3];

  pr[RHO]=0.0;
  pr[UU]=0.0;

  // done.


  // see if E_i = -v^j B^k [ijk]

  
  SLOOPA Ed[j]=0;
  Ed[1] = -(vcon[2]*Mcon[3][TT]-vcon[3]*Mcon[2][TT]);
  Ed[2] = -(vcon[3]*Mcon[1][TT]-vcon[1]*Mcon[3][TT]);
  Ed[3] = -(vcon[1]*Mcon[2][TT]-vcon[2]*Mcon[1][TT]);

  SLOOPA dualfprintf(fail_file,"Ed[%d]=%21.15g Fcov[%d][TT]/detg=%21.15g\n",j,Ed[j],j,Fcov[j][TT]/geom->g);

  ftemp=0.0;
  DLOOPA ftemp+=Ed[j]*Mcon[j][TT];
  dualfprintf(fail_file,"Ed.Bu=%21.15g\n",ftemp);// zero
  ftemp=0.0;
  DLOOPA ftemp+=Fcon[TT][j]*Mcov[j][TT];
  dualfprintf(fail_file,"Eu.Bd=%21.15g\n",ftemp);// zero
  ftemp=0.0;
  DLOOPA ftemp+=Fcov[TT][j]*Mcov[j][TT];
  dualfprintf(fail_file,"Ed.Bd=%21.15g\n",ftemp); // NONzero

  ftemp=0.0;
  DLOOPA ftemp+=vcon[j]*Mcov[j][TT];
  dualfprintf(fail_file,"vu.Bd=%21.15g\n",ftemp); // zero

  ftemp=0.0;
  DLOOPA ftemp+=vcon[j]*Fcov[j][TT];
  dualfprintf(fail_file,"vu.Ed=%21.15g\n",ftemp); // zero


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

      /* lower second index of Maxwell to get Mx=M^t_i */
      lower(Mt,geom,Mx);

      // Bsq = \dF^t_\lambda \dF^{t\lambda} = Mx[1]*Mt1 + Mx[2]*Mt2 + Mx[3]*Mt3 ;
      Bsq=0.0;
      DLOOPA Bsq+=Mx[j]*Mt[j];

      DLOOPA dualfprintf(fail_file,"Mt[%d]=%21.15g T0[%d]=%21.15g Mx[%d]=%21.15g Bsq=%21.15g\n",j,Mt[j],j,T0[j],j,Mx[j],Bsq);
      
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
      if(Bsq!=0.0){
	Fit[0] = 0;
	Fit[1] = (Mx[2]*T0[3] - Mx[3]*T0[2])/(Bsq*geom->g) ;
	Fit[2] = (Mx[3]*T0[1] - Mx[1]*T0[3])/(Bsq*geom->g) ;
	Fit[3] = (Mx[1]*T0[2] - Mx[2]*T0[1])/(Bsq*geom->g) ;
      }
      else DLOOPA Fit[j]=0.0;

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
// mhd_calc(pr,dir,q,mhd) where T^{dir}_j is returned
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

// Maxwell to Faraday
// which=0 : Mcon and Fcov (for clean Mcon, Fcov has \detg)
// which=1 : Mcov and Fcon (for clean Mcov
// copies faraday_calc() in phys.c
void MtoF(int which, FTYPE Max[NDIM][NDIM],struct of_geom *geom, FTYPE faraday[NDIM][NDIM])
{
  int nu,mu,kappa,lambda;
  FTYPE prefactor;
  FTYPE testfactor;


  if((which==0)||(which==1)) prefactor=-0.5;

  for(nu=0;nu<NDIM;nu++) for(mu=0;mu<NDIM;mu++){
    faraday[mu][nu]=0.0;
    for(kappa=0;kappa<NDIM;kappa++) for(lambda=0;lambda<NDIM;lambda++){
      testfactor=prefactor*lc4(which,geom->g,mu,nu,kappa,lambda)*Max[kappa][lambda];
      faraday[mu][nu]+=testfactor;
      //      dualfprintf(fail_file,":: %d %d %d %d : %21.15g %21.15g %21.15g %21.15g\n",mu,nu,kappa,lambda,prefactor,testfactor,faraday[mu][nu],Max[kappa][lambda]);
    }
  }


}

// assumes geometry exists
void testffdeinversion(void)
{
  int i,j,k;
  struct of_geom geom;
  FTYPE prout[NPR],prin[NPR];
  FTYPE U[NPR];
  struct of_state q;
  void Utoprim_ffde(FTYPE *U, struct of_geom *geom, FTYPE *pr)  ;
  int dualfullfaraday_calc(FTYPE *pr, int dir, struct of_state *q, FTYPE *dualffull);
  void lower_A(FTYPE Acon[NDIM][NDIM], struct of_geom *geom, FTYPE Acov[NDIM][NDIM]);
  void MtoF(int which, FTYPE Max[NDIM][NDIM],struct of_geom *geom, FTYPE faraday[NDIM][NDIM]);
  FTYPE dualf[NDIM][NDIM];
  FTYPE Fcon[NDIM][NDIM];
  FTYPE dualfdd[NDIM][NDIM];
  FTYPE mhd[NDIM][NDIM];
  FTYPE Fcov[NDIM][NDIM];
  FTYPE Ed[NDIM];
  int itest,jtest;


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


  // mhd^{dir}_{comp} = mhd^j_k
  DLOOPA  mhd_calc(prin, j, &q, mhd[j]);

  DLOOP dualfprintf(fail_file,"MHD Tud[%d][%d]=%21.15g\n",j,k,mhd[j][k]);

  // MHD Tud
  dualfprintf(fail_file,"Tudanal={");
  for(j=0;j<NDIM;j++){
    dualfprintf(fail_file,"{");
    for(k=0;k<NDIM;k++){
      dualfprintf(fail_file,"%21.15g``25",mhd[j][k]);
      if(k!=NDIM-1) dualfprintf(fail_file,",");
    }
    if(j!=NDIM-1) dualfprintf(fail_file,"},");
    else dualfprintf(fail_file,"}");
  }
  dualfprintf(fail_file,"}");


  // compare F(B,v) with F(E,B)
  dualfullfaraday_calc(prin, 0, &q, dualf[0]); // dualf[dir]=\dF^{dir i}
  dualfullfaraday_calc(prin, 1, &q, dualf[1]);
  dualfullfaraday_calc(prin, 2, &q, dualf[2]);
  dualfullfaraday_calc(prin, 3, &q, dualf[3]);

  DLOOP dualfprintf(fail_file,"dualf[%d][%d]=%21.15g\n",j,k,-dualf[j][k]); // switch sign since want \dF^{i dir}

  // Mcon
  dualfprintf(fail_file,"dualfcode={");
  for(j=0;j<NDIM;j++){
    dualfprintf(fail_file,"{");
    for(k=0;k<NDIM;k++){
      dualfprintf(fail_file,"%21.15g``25",-dualf[j][k]);
      if(k!=NDIM-1) dualfprintf(fail_file,",");
    }
    if(j!=NDIM-1) dualfprintf(fail_file,"},");
    else dualfprintf(fail_file,"}");
  }
  dualfprintf(fail_file,"}");


  lower_A(dualf,&geom,dualfdd); // \dF_{\mu\nu}
  DLOOP dualfdd[j][k]=-dualfdd[j][k]; // switch sign as noted above
  MtoF(1,dualfdd,&geom,Fcon);
  
  DLOOP dualfprintf(fail_file,"Fconin[%d][%d]=%21.15g\n",j,k,Fcon[j][k]);

  dualfprintf(fail_file,"Fconin={");
  for(j=0;j<NDIM;j++){
    dualfprintf(fail_file,"{");
    for(k=0;k<NDIM;k++){
      dualfprintf(fail_file,"%21.15g``25",Fcon[j][k]);
      if(k!=NDIM-1) dualfprintf(fail_file,",");
    }
    if(j!=NDIM-1) dualfprintf(fail_file,"},");
    else dualfprintf(fail_file,"}");
  }
  dualfprintf(fail_file,"}");


  lower_A(Fcon,&geom,Fcov);

  DLOOP dualfprintf(fail_file,"Fcovin[%d][%d]=%21.15g\n",j,k,Fcov[j][k]);


  Ed[1]=-( (q.ucon[2]/q.ucon[TT])*(-dualf[3][TT])-(q.ucon[3]/q.ucon[TT])*(-dualf[2][TT]));
  Ed[2]=-( (q.ucon[3]/q.ucon[TT])*(-dualf[1][TT])-(q.ucon[1]/q.ucon[TT])*(-dualf[3][TT]));
  Ed[3]=-( (q.ucon[1]/q.ucon[TT])*(-dualf[2][TT])-(q.ucon[2]/q.ucon[TT])*(-dualf[1][TT]));

  SLOOPA dualfprintf(fail_file,"\\detg * Ed[%d]=%21.15g\n",j,geom.g*Ed[j]);





  DLOOPA dualfprintf(fail_file,"prin uu[%d]=%21.15g\n",j,q.ucon[j]);
  SLOOPA dualfprintf(fail_file,"prin vu[%d]=%21.15g\n",j,q.ucon[j]/q.ucon[TT]);

  primtoU(UNOTHING,prin,&q,&geom,U);

  PLOOP dualfprintf(fail_file,"U[%d]=%21.15g\n",k,U[k]);
  Utoprim_ffde(U,&geom,prout); // no need for initial guess since analytic inversion

  // just compare pr in and pr out.
  PLOOP dualfprintf(fail_file,"prold[%d]=%21.15g  prnew[%d]=%21.15g\n",k,prin[k],k,prout[k]);








  PLOOP{
    prin[k]=p[itest][jtest][k];
  }




  // loop through to test consistency, keep same densities in case matters small amount
  // but change field and velocities
  for(k=UU;k<=B3;k++){
    prin[k]=prout[k];
  }

  PLOOP{
    prout[k]=0.0;
  }


  get_geometry(itest,jtest,CENT,&geom);
  get_state(prin,&geom,&q);
  DLOOPA dualfprintf(fail_file,"prin2 uu[%d]=%21.15g\n",j,q.ucon[j]);

  primtoU(UNOTHING,prin,&q,&geom,U);

  PLOOP dualfprintf(fail_file,"U[%d]=%21.15g\n",k,U[k]);
  Utoprim_ffde(U,&geom,prout); // no need for initial guess since analytic inversion

  // just compare pr in and pr out.
  PLOOP dualfprintf(fail_file,"prold[%d]=%21.15g  prnew[%d]=%21.15g :: %21.15g\n",k,prin[k],k,prout[k],(prin[k]-prout[k])/prin[k]);






  myexit(0);



  // compare T(B,v) with T(E,B)



}
