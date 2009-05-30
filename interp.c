
#include "decs.h"
#include "interp.h"



// left means L(i) and right means R(i-1) , where left and right refer to direction of interpolation from center toward edge
void slope_lim(int reallim, int ii, int jj, int kk, FTYPE yll, FTYPE yl, FTYPE yc, FTYPE yr, FTYPE yrr,FTYPE *dq,FTYPE *left,FTYPE *right)
{
  void para(FTYPE x1, FTYPE x2, FTYPE x3, FTYPE x4, FTYPE x5, FTYPE *lout, FTYPE *rout);
  void para2(FTYPE x1, FTYPE x2, FTYPE x3, FTYPE x4, FTYPE x5, FTYPE *lout, FTYPE *rout);
  void para3(FTYPE x1, FTYPE x2, FTYPE x3, FTYPE x4, FTYPE x5, FTYPE *lout, FTYPE *rout);
  void para4(FTYPE x1, FTYPE x2, FTYPE x3, FTYPE x4, FTYPE x5, FTYPE *lout, FTYPE *rout);
  void slope_lim_3points(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq);
  void csslope(int kk, FTYPE x1, FTYPE x2, FTYPE x3, FTYPE x4, FTYPE x5, FTYPE *dq);



  // parabolic reconstruction
  if (reallim == PARA){
    // para
#if(WHICHPARA==PARA1)
    para(yll,yl,yc,yr,yrr,left,right);
#elif(WHICHPARA==PARA2)
    para2(yll,yl,yc,yr,yrr,left,right);
#elif(WHICHPARA==PARA3)
    para3(yll,yl,yc,yr,yrr,left,right);
#elif(WHICHPARA==PARA4)
    para4(yll,yl,yc,yr,yrr,left,right);
#endif
  }
  else if(reallim == CSSLOPE){
    csslope(kk, yll,yl,yc,yr,yrr,dq);
  }
  else{
    slope_lim_3points(reallim, yl, yc, yr,dq);
  }

  if((reallim!=PARA)&&(LIMADJUST!=0)){
    *left =yc - 0.5* (*dq);
    *right=yc + 0.5* (*dq);
  }

}




void slope_lim_3points_old(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq)
{
  FTYPE Dqm, Dqp, Dqc, s;

  if (reallim == MC) {
    Dqm = 2.0 * (yc - yl);
    Dqp = 2.0 * (yr - yc);
    Dqc = 0.5 * (yr - yl);
    s = Dqm * Dqp;
    if (s <= 0.)  *dq= 0.;
    else{
      if (fabs(Dqm) < fabs(Dqp) && fabs(Dqm) < fabs(Dqc))
	*dq= (Dqm);
      else if (fabs(Dqp) < fabs(Dqc))
	*dq= (Dqp);
      else
	*dq= (Dqc);
    }
  }
  /* van leer slope limiter */
  else if (reallim == VANL) {
    Dqm = (yc - yl);
    Dqp = (yr - yc);
    s = Dqm * Dqp;
    if (s <= 0.)
      *dq= 0.;
    else
      *dq= (2.0 * s / (Dqm + Dqp));
  }
  /* minmod slope limiter (crude but robust) */
  else if (reallim == MINM) {
    Dqm = (yc - yl);
    Dqp = (yr - yc);
    s = Dqm * Dqp;
    if (s <= 0.) *dq= 0.;
    else{
      if (fabs(Dqm) < fabs(Dqp)) *dq= Dqm;
      else *dq= Dqp;
    }
  }
  else if (reallim == NLIM) {
    Dqc = 0.5 * (yr - yl);
    *dq= (Dqc);
  }
  else if (reallim == DONOR) {
    *dq=(0.0);
  }
  else {
    dualfprintf(fail_file, "unknown slope limiter: %d\n",reallim);
    myexit(10);
  }



}


void slope_lim_3points(int reallim, FTYPE yl, FTYPE yc, FTYPE yr,FTYPE *dq)
{
  FTYPE Dqm, Dqp, Dqc, s;

  if (reallim == MC) { // monotonized central (Woodward) (Barth-Jespersen)
    Dqm = 2.0 * (yc - yl);
    Dqp = 2.0 * (yr - yc);
    Dqc = 0.5 * (yr - yl);
    *dq=MINMOD(Dqc,MINMOD(Dqm,Dqp));
  }
  /* van leer slope limiter */
  else if (reallim == VANL) {
    Dqm = (yc - yl);
    Dqp = (yr - yc);
    s = Dqm * Dqp;
    if (s <= 0.)
      *dq= 0.;
    else
      *dq= (2.0 * s / (Dqm + Dqp));
  }
  /* minmod slope limiter (crude but robust) */
  else if (reallim == MINM) {
    Dqm = (yc - yl);
    Dqp = (yr - yc);
    *dq = MINMOD(Dqm,Dqp);
  }
  else if (reallim == NLIMCENT) { // Centered slope (Fromm)
    *dq = 0.5 * (yr - yl);
  }
  else if (reallim == NLIMUP) { // Upwind slope (Beam-Warming)
    *dq = (yc - yl);
  }
  else if (reallim == NLIMDOWN) { // Downwind slope (Lax-Wendroff)
    *dq = (yr - yc);
  }
  else if (reallim == DONOR) { // no slope
    *dq=(0.0);
  }
  else {
    dualfprintf(fail_file, "unknown slope limiter: %d\n",reallim);
    myexit(10);
  }



}


// see Mignone & Bodo (2005) astro-ph/0506414 equations 27-31
// see Colella (1985), Saltzman (1994)
// second order with 4th order steepeners
void csslope(int kk, FTYPE x1, FTYPE x2, FTYPE x3, FTYPE x4, FTYPE x5, FTYPE *dq)
{
  FTYPE s, sm, sp;
  FTYPE Dql, Dqlm, Dqlp;
  FTYPE Dqp, Dqpp;
  FTYPE Dqm, Dqmm;
  FTYPE Dqc, Dqcm, Dqcp;
  FTYPE Dqbp,Dqbm;
  FTYPE alpha;

  Dqp=(x4-x3);
  Dqpp=(x5-x4);

  Dqm=(x3-x2);
  Dqmm=(x2-x1);

  Dqc=0.5*(x4-x2);
  Dqcm=0.5*(x3-x1);
  Dqcp=0.5*(x5-x3);

  // Dqmm  Dqm Dqp  Dqpp
  //    Dqcm Dqc Dqcp
  //     sm   q   sp
  //    Dqlm Dql Dqlp
  //    Dqbm     Dqbp

  s=0.5*(sign(Dqp)+sign(Dqm));
  sp=0.5*(sign(Dqpp)+sign(Dqp));
  sm=0.5*(sign(Dqm)+sign(Dqmm));

  // MB05 use alpha=2 for 1-D and alpha=2,1.25,1 for rho,v,p (respectively) for 2D.
  // alpha=[1-2].  alpha=2 more compressive, alpha=1 less compressive.
  if(kk==RHO) alpha=2.0;
  else if((kk>=U1)&&(kk<=U3)) alpha=1.25;
  else if(kk==UU) alpha=1.0;
  else if((kk>=B1)&&(kk<=B3)) alpha=1.25;
  else alpha=2.0;

  Dql=alpha*min(fabs(Dqp),fabs(Dqm));
  Dqlp=alpha*min(fabs(Dqpp),fabs(Dqp));
  Dqlm=alpha*min(fabs(Dqm),fabs(Dqmm));

  Dqbp=sp*min(Dqlp,fabs(Dqcp));
  Dqbm=sm*min(Dqlm,fabs(Dqcm));

  *dq=s*min(fabs(FOURTHIRD*Dqc-SIXTH*(Dqbp+Dqbm)),Dql);
  

}



/*
 * parabolic interpolation subroutin  
 * ref. Colella && Woodward's paper
 * Colella, P., & Woodward, P. R. 1984, J. Comput. Phys., 54, 174-201
 *
 * using zone-centered value of 5 continuous zones 
 * to get left and right value of the middle zone.
 *  
 * 
 */

// lout/rout is left and right sides of cell
// note how used in step_ch.c to get appropriate interface value

// given by Xiaoyue Guan to Scott Noble on Nov 9, 2004, given to me Jan 7, 2005
void para(FTYPE x1, FTYPE x2, FTYPE x3, FTYPE x4, FTYPE x5, FTYPE *lout, FTYPE *rout)
{
  int i ;
  FTYPE y[5], dq[5];
  FTYPE Dqm, Dqc, Dqp, aDqm,aDqp,aDqc,s,l,r,qa, qd, qe;

  y[0]=x1;
  y[1]=x2;
  y[2]=x3;
  y[3]=x4;
  y[4]=x5;

  /*CW1.7 */
  for(i=1 ; i<4 ; i++) {
    Dqm = 2.0 *(y[i]-y[i-1]);
    Dqp = 2.0 *(y[i+1]-y[i]);
    Dqc = 0.5 *(y[i+1]-y[i-1]);
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    s = Dqm*Dqp;

    if (s <=0.) dq[i]=0.;       //CW1.8
    else dq[i]=min(aDqc,min(aDqm,aDqp))*sign(Dqc);
  }

  /* CW1.6 */

  l=0.5*(y[2]+y[1])-(dq[2]-dq[1])/6.0;
  r=0.5*(y[3]+y[2])-(dq[3]-dq[2])/6.0;

  qa=(r-y[2])*(y[2]-l);
  qd=(r-l);
  qe=6.0*(y[2]-0.5*(l+r));


  if (qa <=0. ) {
    l=y[2];
    r=y[2];
  }

  if (qd*(qd-qe)<0.0) l=3.0*y[2]-2.0*r;
  else if (qd*(qd+qe)<0.0) r=3.0*y[2]-2.0*l;
	 

  *lout=l;   //a_L,j
  *rout=r;
  //*dw=r-l;                      //CW1.5
  //*w6=6.0*(y[2]-0.5*(l+r));
}



// given by Xiaoyue Guan on Jan 9, 2005
void para2(FTYPE x1, FTYPE x2, FTYPE x3, FTYPE x4, FTYPE x5, FTYPE *lout, FTYPE *rout)
{
  int i ;
  FTYPE y[5], dq[5];
  FTYPE Dqm, Dqc, Dqp, Dqvanl,aDqm,aDqp,aDqc,aDqvanl,s,l,r,qa, qd, qe;

  y[0]=x1;
  y[1]=x2;
  y[2]=x3;
  y[3]=x4;
  y[4]=x5;

  /*CW1.7 */
  for(i=1 ; i<4 ; i++) {
    Dqm = 2.0 *(y[i]-y[i-1]);
    Dqp = 2.0 *(y[i+1]-y[i]);
    Dqc = 0.5 *(y[i+1]-y[i-1]);
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    s = Dqm*Dqp;
	       
#if(PARA2LIM == VANL) 
    Dqvanl=2.0*Dqm*Dqp/(Dqm+Dqp);
    aDqvanl=fabs(Dqvanl);
    if (s <=0.) dq[i]=0.;       //CW1.8
    else dq[i]=min(min(aDqc,aDqvanl),min(aDqm,aDqp))*sign(Dqc);
#elif(PARA2LIM == PMC)
    if (s <=0.) dq[i]=0.;       //CW1.8
    else dq[i]=min(aDqc,min(aDqm,aDqp))*sign(Dqc);
#elif(PARA2LIM == MC)
    dq[i] =Dqc;
#endif
  }
  /* CW1.6 */

  l=0.5*(y[2]+y[1])-(dq[2]-dq[1])/6.0;
  r=0.5*(y[3]+y[2])-(dq[3]-dq[2])/6.0;
	 
  /*
    l=max(min(y[2],y[1]),l);
    l=min(max(y[2],y[1]),l);
    r=max(min(y[2],y[3]),r);
    r=min(max(y[2],y[3]),r);
  */
	 
  qa=(r-y[2])*(y[2]-l);
  qd=(r-l);
  qe=6.0*(y[2]-0.5*(l+r));


  if (qa <=0. ) {
    l=y[2];
    r=y[2];
  }

  else if (qd*(qd-qe)<0.0) l=3.0*y[2]-2.0*r;
  else if (qd*(qd+qe)<0.0) r=3.0*y[2]-2.0*l;
	 

  *lout=l;   //a_L,j
  *rout=r;
  //*dw=r-l;                      //CW1.5
  //*w6=6.0*(y[2]-0.5*(l+r));
}




// 3rd para from Xiaoyue that she bundled with a new TVD-optimal RK3
// given on 02/17/2005
void para3(FTYPE x1, FTYPE x2, FTYPE x3, FTYPE x4, FTYPE x5, FTYPE *lout, FTYPE *rout)
{
  int i ;
  FTYPE y[5], dq[5];
  FTYPE Dqm, Dqc, Dqp, Dqvanl,aDqm,aDqp,aDqc,aDqvanl,s,l,r,qa, qd, qe;

  y[0]=x1;
  y[1]=x2;
  y[2]=x3;
  y[3]=x4;
  y[4]=x5;

  /*CW1.7 */
  for(i=1 ; i<4 ; i++) {
    Dqm = 2.0 *(y[i]-y[i-1]);
    Dqp = 2.0 *(y[i+1]-y[i]);
    Dqc = 0.5 *(y[i+1]-y[i-1]);
    aDqm = fabs(Dqm) ;
    aDqp = fabs(Dqp) ;
    aDqc = fabs(Dqc) ;
    s = Dqm*Dqp;
	       
#if(PARA2LIM == VANL) 
    Dqvanl=2.0*Dqm*Dqp/(Dqm+Dqp);
    aDqvanl=fabs(Dqvanl);

    if (s <=0.) dq[i]=0.;
    else dq[i] = -aDqvanl*sign(Dqc);
    //else dq[i]=min(min(aDqc,aDqvanl),min(aDqm,aDqp))*sign(Dqc);
#elif(PARA2LIM == MC)
    if (s <=0.) dq[i]=0.;       //CW1.8
    else dq[i]=-min(aDqc,min(aDqm,aDqp))*sign(Dqc);
#elif(PARA2LIM == MINM)
    if (s<=0.) dq[i] = 0.;
    else if (aDqm<aDqp) dq[i] = -aDqm*sign(Dqc);
    else dq[i]=-aDqp*sign(Dqc);
#elif(PARA2LIM == NLIM) //w/o slope limiter
    //if(s<=0.) dq[i] = 0.; // DONOR
    dq[i] = Dqc;
#endif
  }

  /* CW1.6 */

  l=0.5*(y[2]+y[1])-(dq[2]-dq[1])/6.0;
  r=0.5*(y[3]+y[2])-(dq[3]-dq[2])/6.0;
	 
	 
  l=max(min(y[2],y[1]),l);
  l=min(max(y[2],y[1]),l);
  r=max(min(y[2],y[3]),r);
  r=min(max(y[2],y[3]),r);
	 
	 
  qa=(r-y[2])*(y[2]-l);
  qd=(r-l);
  qe=6.0*(y[2]-0.5*(l+r));
	 
  /*
    if (qa <=0. ) {
    l=y[2];
    r=y[2];
    }

    else if (qd*(qd-qe)<0.0) l=3.0*y[2]-2.0*r;
    else if (qd*(qd+qe)<0.0) r=3.0*y[2]-2.0*l;
	 

    *lout=l;   //a_L,j
    *rout=r;
    */

  if (qa <=0. ) {
    *lout=y[2];
    *rout=y[2];
  }  
  else {
    *lout = l;
    *rout = r;
  }
		  
  //2.0 at top/bottom of a steep gradient 
  if (qd*(qd-qe)<0.0) *lout=3.0*y[2]-2.0*r;
  else *lout = l;
 
  if (qd*(qd+qe)<0.0) *rout=3.0*y[2]-2.0*l;
  else *rout = r;
  //*dw=r-l;                      //CW1.5
  //*w6=6.0*(y[2]-0.5*(l+r));
}






// Xiaoyue given on 03/25/05
// she realized sign error in 1st der's in para3()
// noted Matt's paper astro-ph/0503420 suggested CW1.6 uses 1/8 rather than 1/6
// I noted Matt uses MC for field variables and PPM+ for hydro variables
void para4(FTYPE x1, FTYPE x2, FTYPE x3, FTYPE x4, FTYPE x5, FTYPE *lout, FTYPE *rout)
{
  int i ;
  FTYPE y[5], dq[5];
  FTYPE Dqm, Dqc, Dqp, Dqvanl,aDqm,aDqp,aDqc,aDqvanl,s,l,r,qa, qd, qe;

  y[0]=x1;
  y[1]=x2;
  y[2]=x3;
  y[3]=x4;
  y[4]=x5;

  /*CW1.7 */
  for(i=1 ; i<4 ; i++) 
    {
      Dqm = 2.0 *(y[i]-y[i-1]);
      Dqp = 2.0 *(y[i+1]-y[i]);
      Dqc = 0.5 *(y[i+1]-y[i-1]);
      aDqm = fabs(Dqm) ;
      aDqp = fabs(Dqp) ;
      aDqc = fabs(Dqc) ;
      s = Dqm*Dqp;


#if(PARA2LIM == VANL) 
      Dqvanl=2.0*Dqm*Dqp/(Dqm+Dqp);
      aDqvanl=fabs(Dqvanl);
	     
      if (s <=0.) dq[i]=0.;
      //else dq[i] = aDqvanl*sign(Dqc);
      else dq[i]=min(min(aDqc,aDqvanl),min(aDqm,aDqp))*sign(Dqc);
	     
#elif(PARA2LIM == MC)
	     
      if (s <=0.) dq[i]=0.;       //CW1.8
      else dq[i]= min(aDqc,min(aDqm,aDqp))*sign(Dqc);
	     
#elif(PARA2LIM == MINM)
	     
      if (s<=0.) dq[i] = 0.;
      else if (aDqm<aDqp) dq[i] = aDqm*sign(Dqc);
      else dq[i]=aDqp*sign(Dqc);
	     
#elif(PARA2LIM == NLIM) //w/o slope limiter
	     
      dq[i] = Dqc;
#endif
    }
         
	 
  /* CW1.6 */

  // modified as per Matt's paper
  l=0.5*(y[2]+y[1])-(dq[2]-dq[1])/8.0;
  r=0.5*(y[3]+y[2])-(dq[3]-dq[2])/8.0;
	 
	 
  l=max(min(y[2],y[1]),l);
  l=min(max(y[2],y[1]),l);
  r=max(min(y[2],y[3]),r);
  r=min(max(y[2],y[3]),r);
	 
	 
  // modified as per Matt's paper
  qa=(r-y[2])*(y[2]-l);
  qd=(r-l);
  qe=8.0*(y[2]-0.5*(l+r));
	 
	 
  if (qa <=0. ) {
    l=y[2];
    r=y[2];
  }
  else{

    if (qd*(qd-qe)<0.0) 
      l=3.0*y[2]-2.0*r;
     
	 
    if (qd*(qd+qe)<0.0) 
      r=3.0*y[2]-2.0*l;
  }
	 

  *lout=l;   //a_L,j
  *rout=r;
	 


}















// notice that input order is always:
// 2nd arg: true primitive
// last arg: scaled primitive

// this also allows any fancy remapping, such as characteristic interpolation -- rest of code is setup to allow any remapping as long as you have an inversion
int rescale(int which, int dir, FTYPE *pr, struct of_geom *ptrgeom,FTYPE *p2interp)
{
  FTYPE scale[NPR],r,th,X[NDIM];
  int ii,jj,k;
  struct of_state q;




  ii=ptrgeom->i;
  jj=ptrgeom->j;
  coord(ii,jj,ptrgeom->p,X);
  bl_coord(X,&r,&th);

#if(VARTOINTERP==PRIMTOINTERP)
  if(dir==1){
    // optimized for pole
    if(rescaletype==0){
      scale[RHO]=pow(r,-1.5);
    }
    else if(rescaletype==1){
      scale[RHO]=pow(r,-2.7);
    }
    scale[UU]=scale[RHO]/r;
    scale[U1]=scale[RHO];
    scale[U2]=1.0;
    scale[U3]=1.0/(r*r);
    if(rescaletype==0){
      scale[B1]=scale[U3];
    }
    else if(rescaletype==1){
      scale[B1]=pow(r,-2.4);
    }
    //    if(statpos[2]+jj < 0 || startpos[2]+jj >= totalsize[2]) scale[B1] *= -1. ;
    scale[B2]=scale[B1];
    scale[B3]=scale[B1];
    
    if(DOENTROPY) scale[ENTROPY]=1.0;
  }
  else if(dir==2){
    scale[RHO]=1.0;
    scale[UU]=1.0;
    scale[U1]=1.0;
    scale[U2]=1.0;
    scale[U3]=1.0;
    scale[B1]=1.0;
    scale[B2]=1.0;
    scale[B3]=1.0;
    if(DOENTROPY) scale[ENTROPY]=1.0;
  }
  else{
    dualfprintf(fail_file,"rescale(): no such direction! dir=%d\n",dir);
    myexit(100);
  }

  if(which==1){ // rescale before interpolation
    PLOOP p2interp[k]=pr[k]/scale[k];
  }
  else if(which==-1){ // unrescale after interpolation
    PLOOP pr[k]=p2interp[k]*scale[k];
  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }




#elif(VARTOINTERP==CONSTOINTERP)
  // this doesn't work at all, even if no bug.
  // doesn't work even if setup better guess, as in interpU code.

  
  if(which==1){ // rescale before interpolation
    MYFUN(get_state(pr, ptrgeom, &q),"interp.c:rescale()", "get_state()", 1);
    MYFUN(primtoU(UDIAG,pr, &q, ptrgeom, p2interp),"interp.c:rescale()", "primtoU()", 1);
  }
  else if(which==-1){ // unrescale after interpolation
    MYFUN(Utoprimgen(OTHERUTOPRIM,UDIAG,p2interp, ptrgeom, pr),"interp.c:rescale()", "Utoprimgen", 1);
  }
  else{
    dualfprintf(fail_file,"rescale(): no such rescale type! which=%d\n",which);
    myexit(100);
  }



#endif
  return(0);
}
