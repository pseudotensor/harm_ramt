
#include "decs.h"

/** 
 *
 * this file contains all the coordinate dependent
 * parts of the code, except the initial and boundary
 * conditions 
 *
 **/

// static variables with global scope to this file
// could make any or all of these true global if want to change them in, say, init.c

// for defcoord==4
static FTYPE der0=9;
//static FTYPE R=8.0; // too unresolved for 256x128
static FTYPE R=20.0;

// for defcoord==5
static FTYPE x2trans;
static FTYPE m2,d2,c2,m3,b3,myx2,flip1,flip2,thetatores;

// for defcoord==7
static FTYPE h0,hf,rh0,myrout,myhslope1,dmyhslope1dr,myhslope1,myhslope2,myhslope,dmyhslope2dx1,dmyhslopedx1,x1in,x1out;
static FTYPE npow;

// for defcoord=8
static FTYPE r1jet,njet,rpjet,dmyhslopedr;

// for defcoord=9
static FTYPE r0jet,rsjet,Qjet;

// for defcoord=10
static FTYPE hinner,houter;


void set_coord_parms(void)
{
  // assumes R0, Rin, Rout, and hslope are so general that are set in init.c
  if (defcoord == 0) {
  }
  else if (defcoord == 1) {
  }
  else if(defcoord==2){
  }
  else if (defcoord == 3) {
  }
  else if(defcoord == 4) {
  }
  else if(defcoord == 5) {
    x2trans=0.1; // user settable, must be same as below in dxdxp
    thetatores=2.5*h_over_r;

    // fixed coefficients
    m2=(3.*(-2.*thetatores + M_PI))/(2.*x2trans) + (4.*thetatores)/(-1. + 2.*x2trans);
    d2=(2.*thetatores - M_PI + 2.*M_PI*x2trans)/(-2.*pow(x2trans,3.) + 4.*pow(x2trans,4.));
    c2=(6.*thetatores - 3.*M_PI + 6.*M_PI*x2trans)/(2.*pow(x2trans,2.) - 4.*pow(x2trans, 3.));
    m3=(2.*thetatores)/(1. - 2.*x2trans);
    b3=M_PI/2. + thetatores/(-1. + 2.*x2trans);
  }
  else if (defcoord == 6) { // uniform theta and log in radius
  }
  else if (defcoord == 7) {
    // optimal is npow=10 R0=-3
    npow=1.0;
    //R0=0.0;

    // must be same as in dxdxp()
    h0=hslope;
    hf=2.0-0.22;
    rh0=40.0;
    myrout=Rout;
    dmyhslope1dr = (hf-h0)/(myrout-rh0);
    dmyhslope2dx1=(hf-h0)/(x1out-x1in);
    x1in=log(Rin-R0);
    x1out=log(Rout-R0);

  }
  else if (defcoord == 8) {
    npow=1.0;

    // must be same as in dxdxp()
    if(1){
      r1jet=16.0;
      njet=0.3;
      rpjet=0.9;
    }
    else{
      r1jet=9.0;
      njet=0.3;
      rpjet=.9;
    }
  }
  else if (defcoord == 9) {
    npow=10.0;

    // must be same as in dxdxp()
    if(0){ // first attempt
      r1jet=2.8;
      njet=0.3;
      r0jet=7.0;
      rsjet=21.0;
      Qjet=1.7;
    }
    else if(0){ // chosen to resolve disk then resolve jet
      r1jet=2.8;
      njet=0.3;
      r0jet=20.0;
      rsjet=80.0;
      Qjet=1.8;
    }
    else if(1){
      r1jet=2.8;
      njet=0.3;
      r0jet=20.0;
      rsjet=80.0;
      Qjet=1.3; // chosen to help keep jet resolved even within disk region
    }
  }
  else if (defcoord == 10) {
    // pulsar_grid.nb for theta part and for the radial part:
    // see pulsar_gridnew.nb
    // for Rout=10^6 and R0=0.786*Rin Rin=4.84, npow=10 gives same dr/r as npow=1 R0=0.9*Rin at r=Rin
    npow=1.0;

    // must be same as in dxdxp()
    hinner=hslope; // hslope specifies inner hslope
    houter=hslope*0.05; // reduce by some arbitrary factor (currently 1/20)
    r0jet=5.0; // spread in radius over which hslope changes
    rsjet=18.0; // location of current sheet beginning for NS pulsar

  }
  else if (defcoord == 666) {
    npow=10.0; // exponential rate
  }
  else{
    dualfprintf(fail_file,"Shouldn't reach end of set_coord_parms\n");
    myexit(1);
  }

}


void write_coord_parms(void)
{
  FILE *out;

  if(myid==0){
    if((out=fopen("coordparms.dat","wt"))==NULL){
      dualfprintf(fail_file,"Couldn't write coordparms.dat file\n");
      myexit(1);
    }
    else{

      // same for all coords (notice no carraige return)
      fprintf(out,"%21.15g %21.15g %21.15g %21.15g ",R0,Rin,Rout,hslope);

      if (defcoord == 0) {
      }
      else if (defcoord == 1) {
      }
      else if(defcoord==2){
      }
      else if (defcoord == 3) {
      }
      else if(defcoord == 4) {
      }
      else if(defcoord == 5) {
	fprintf(out,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",x2trans,thetatores,m2,d2,c2,m3,b3,h_over_r);
      }
      else if (defcoord == 6) { // uniform theta and log in radius
      }
      else if (defcoord == 7) {
	fprintf(out,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",npow,h0,hf,rh0,myrout,dmyhslope1dr,dmyhslope2dx1,x1in,x1out);
      }
      else if (defcoord == 8) {
	fprintf(out,"%21.15g %21.15g %21.15g %21.15g\n",npow,r1jet,njet,rpjet);
      }
      else if (defcoord == 9) {
	fprintf(out,"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",npow,r1jet,njet,r0jet,rsjet,Qjet);
      }
      else if (defcoord == 10) {
	fprintf(out,"%21.15g %21.15g %21.15g %21.15g %21.15g\n",npow,hinner,houter,r0jet,rsjet);
      }
      else if (defcoord == 666) {
	fprintf(out,"%21.15g\n",npow);
      }
      else{
	dualfprintf(fail_file,"Shouldn't reach end of write_coord_parms\n");
	myexit(1);
      }

      fclose(out);
    }
  }
}



void read_coord_parms(void)
{
  FILE *in;
  FTYPE ftemp;
  
  if(myid==0){
    in=fopen("coordparms.dat","rt");
    if(in==NULL){
      dualfprintf(fail_file,"Couldn't read coordparms.dat file.  I'll assume coded coordinates and let restart header overwrite any global restart parameters\n");
      set_coord_parms();
    }
    else{
      // don't want to overwrite since restart file sets this
      //      fscanf(in,HEADER4IN,&R0,&Rin,&Rout,&hslope);
      fscanf(in,HEADER4IN,&ftemp,&ftemp,&ftemp,&ftemp);

      if (defcoord == 0) {
      }
      else if (defcoord == 1) {
      }
      else if(defcoord==2){
      }
      else if (defcoord == 3) {
      }
      else if(defcoord == 4) {
      }
      else if(defcoord == 5) {
	fscanf(in,HEADER8IN,&x2trans,&thetatores,&m2,&d2,&c2,&m3,&b3,&h_over_r);
      }
      else if (defcoord == 6) { // uniform theta and log in radius
      }
      else if (defcoord == 7) {
	fscanf(in,HEADER9IN,&npow,&h0,&hf,&rh0,&myrout,&dmyhslope1dr,&dmyhslope2dx1,&x1in,&x1out);
      }
      else if (defcoord == 8) {
	fscanf(in,HEADER4IN,&npow,&r1jet,&njet,&rpjet);
      }
      else if (defcoord == 9) {
	fscanf(in,HEADER6IN,&npow,&r1jet,&njet,&r0jet,&rsjet,&Qjet);
      }
      else if (defcoord == 10) {
	fscanf(in,HEADER5IN,&npow,&hinner,&houter,&r0jet,&rsjet);
      }
      else if (defcoord == 666) {
	fscanf(in,HEADERONEIN,&npow);
      }
      else{
	dualfprintf(fail_file,"Shouldn't reach end of read_coord_parms\n");
	myexit(1);
      }

      fclose(in);
    }
  }

#if(USEMPI)
  // broadcast
  MPI_Bcast(&R0, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Rin, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&Rout, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&hslope, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);

  if (defcoord == 0) {
  }
  else if (defcoord == 1) {
  }
  else if(defcoord==2){
  }
  else if (defcoord == 3) {
  }
  else if(defcoord == 4) {
  }
  else if(defcoord == 5) {
    MPI_Bcast(&x2trans, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&thetatores, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m2, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&d2, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&c2, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&m3, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&b3, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    //  MPI_Bcast(&h_over_r, 1, MPI_FTYPE, 0, MPI_COMM_WORLD); // set by pre_init_specific_init() in init.c
  }
  else if (defcoord == 6) { // uniform theta and log in radius
  }
  else if (defcoord == 7) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&h0, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&hf, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rh0, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&myrout, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dmyhslope1dr, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&dmyhslope2dx1, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&x1in, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&x1out, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  }
  else if (defcoord == 8) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&r1jet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&njet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rpjet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  }
  else if (defcoord == 9) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&r1jet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&njet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&r0jet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rsjet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&Qjet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  }
  else if (defcoord == 10) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&hinner, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&houter, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&r0jet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
    MPI_Bcast(&rsjet, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  }
  else if (defcoord == 666) {
    MPI_Bcast(&npow, 1, MPI_FTYPE, 0, MPI_COMM_WORLD);
  }
  else{
    dualfprintf(fail_file,"Shouldn't reach end of read_coord_parms\n");
    myexit(1);
  }
  
#endif

}





/* Returns boyer-lindquist coordinte of point */
void bl_coord(FTYPE *X, FTYPE *r, FTYPE *th)
{


  if (defcoord == 0) {
    //    *r = Rin*exp(X[1]) ;
    *r = R0+exp(X[1]) ;
    //*r = Rin * exp(X[1]);
    *th = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
  } else if (defcoord == 1) {
    *r = R0+Rin * exp(X[1] * log(Rout / Rin));
    *th =
	((-49. * hslope + 60. * M_PI) * X[2]) / 12. +
	((247. * hslope - 240. * M_PI) * pow(X[2],2)) / 12. +
	((-83. * hslope + 80. * M_PI) * pow(X[2],3)) / 2. -
	(5. * (-25. * hslope + 24. * M_PI) * pow(X[2], 4)) / 3. +
	(2. * (-25. * hslope + 24. * M_PI) * pow(X[2], 5)) / 3.;
  }
  else if(defcoord==2){
    *r=X[1];
    *th = M_PI* X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
  }
  else if (defcoord == 3) {
    // MIRROR at equator, equator is outer theta edge
    *r = R0+exp(X[1]) ;
    *th = M_PI * X[2] + ((1. - hslope) / 2.) * sin(2. * M_PI * X[2]);
  }
  else if(defcoord == 4) {
    *r = R0+exp(X[1]) ;

    *th = (der0*X[2]*(-32.*pow(-1. + X[2],3.)*pow(X[2],2.)*(-1. + 2.*X[2]) - 
		   R*(-1. + X[2])*pow(-1. + 2.*X[2],3.)*
		   (-1. + 7.*(-1. + X[2])*X[2])) + 
	   M_PI*R*pow(X[2],3.)*(70. + 
			     3.*X[2]*(-105. + 2.*X[2]*(91. + 10.*X[2]*(-7. + 2.*X[2])))))/R;

  }
  else if(defcoord == 5) {

    *r = R0+exp(X[1]) ;

    // now assign values
    if(X[2]<0.5){ myx2=X[2]; flip1=0.0; flip2=1.0;}
    else{ myx2=1.0-X[2]; flip1=M_PI; flip2=-1.0;}

    if(myx2<=x2trans){
      *th = flip1+flip2*(d2*pow(myx2,3.0)+c2*pow(myx2,2.0)+m2*myx2);
    }
    else{
      *th = flip1+flip2*(m3*myx2+b3);
    }

  }
  else if (defcoord == 6) { // uniform theta and log in radius
    *r = R0+exp(X[1]) ;
    *th = M_PI * X[2] ;
  }
  else if (defcoord == 7) {
    *r = R0+exp(pow(X[1],npow)) ;

    myhslope1 = h0+dmyhslope1dr*((*r)-rh0);
    myhslope2 = h0+(hf-h0)*(X[1]-x1in)/(x1out-x1in);

    myhslope=(myhslope1+myhslope2)*0.5;
    *th = M_PI * X[2] + ((1. - myhslope) / 2.) * sin(2. * M_PI * X[2]);
  }
  else if (defcoord == 8) {
    *r = R0+exp(pow(X[1],npow)) ;

    myhslope=2.0-pow(*r/r1jet,njet*(-1+exp(1)/exp(*r+rpjet)));

    *th = M_PI * X[2] + ((1. - myhslope) / 2.) * sin(2. * M_PI * X[2]);
  }
  else if (defcoord == 9) {
    *r = R0+exp(pow(X[1],npow)) ;

    myhslope=2.0-Qjet*pow(*r/r1jet,-njet*(0.5+1.0/M_PI*atan(*r/r0jet-rsjet/r0jet)));

    *th = M_PI * X[2] + ((1. - myhslope) / 2.) * sin(2. * M_PI * X[2]);
  }
  else if (defcoord == 10) {
    *r = R0+exp(pow(X[1],npow)) ;

    myhslope=(0.5+1.0/M_PI*atan((*r-rsjet)/r0jet))*(houter-hinner)+hinner;

    *th = M_PI * X[2] + ((1. - myhslope) / 2.) * sin(2. * M_PI * X[2]);
  }
  else if (defcoord == 666) {
    // R : cylindrical radius, assumes X[1]=0..1
    // exponential grid
    *r = ((Rout-Rin)*exp(npow*X[1])+Rin*exp(npow)-Rout)/(exp(npow)-1.0);

    // z : cylindrical height, assumes X[2]=-1..1
    // bi-exponential grid
    // here the grid goes from Zin to Zout in a bi-log way, and X[2]=0 is Z=0
    if(X[2]>0.0) *th = ((Zout-0)*exp(npow*fabs(X[2])) + 0*exp(npow)-Zout)/(exp(npow)-1.0);
    else *th = ((Zin-0)*exp(npow*fabs(X[2])) + 0*exp(npow)-Zin)/(exp(npow)-1.0);
    
  }  else{
    dualfprintf(fail_file,"Shouldn't reach end of bl_coord\n");
    myexit(1);
  }

  // don't allow to be smaller to avoid singularity
  // noted this caused problems with jon_interp in calculating jacobian
  if(POSDEFMETRIC){
    if(*th<0) *th = -*th;
    if(*th>M_PI) *th=2.0*M_PI-*th;
  }

  if(COORDSINGFIX){
    if (fabs(*th) < SINGSMALL){
      if(*th>=0) *th=SINGSMALL;
      if(*th<0) *th=-SINGSMALL;
    }
  }



}



// Jacobian for dx uniform per dx nonuniform (dx/dr / dx/dr')
// i.e. Just take d(bl-coord)/d(ksp uniform coord)
// e.g. dr/dx1 d\theta/dx2

// take note of the ordering of indicies
// dxdxp[j][k]=dxdxp[mu][nu]=(dx^\mu_{BL}/dx^\nu_{KSP uni})

// should make this numerical like connection, then to conserve CPU, would need all over grid
void dxdxprim(FTYPE *X, FTYPE r, FTYPE th, FTYPE (*dxdxp)[NDIM])
{
  void dxdxp_numerical(FTYPE *X, FTYPE (*dxdxp)[NDIM]);
  void dxdxp_analytic(FTYPE *X, FTYPE r, FTYPE th, FTYPE (*dxdxp)[NDIM]);

  if(defcoord<=10){ // then have analytic dxdxp
    dxdxp_analytic(X,r,th,dxdxp);
  }
  else{
    dxdxp_numerical(X,dxdxp);
  }

}

// should make this numerical like connection, then to conserve CPU, would need all over grid
void dxdxp_analytic(FTYPE *X, FTYPE r, FTYPE th, FTYPE (*dxdxp)[NDIM])
{
  int j,k;

  // default identity transformation
  DLOOP dxdxp[j][k]=0.0;
  DLOOPA dxdxp[j][j]=1.0;

  if (defcoord == 0) {
    dxdxp[1][1] = r-R0;
    dxdxp[2][2] = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);
  } else if (defcoord == 1) {
    dxdxp[1][1] = (r-R0) * log(Rout / Rin);

    dxdxp[2][2] = (-49. * hslope + 60. * M_PI) / 12. +
	((247. * hslope - 240. * M_PI) * X[2]) / 6. +
	(3. * (-83. * hslope + 80. * M_PI) * pow(X[2], 2)) / 2. -
	(20. * (-25. * hslope + 24. * M_PI) * pow(X[2], 3)) / 3. +
	(10. * (-25. * hslope + 24. * M_PI) * pow(X[2], 4)) / 3.;

  } else if(defcoord==2){
    dxdxp[2][2] = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);
  }
  else if (defcoord == 3) {
    dxdxp[1][1] = r-R0;
    dxdxp[2][2] = M_PI + (1. - hslope) * M_PI * cos(2. * M_PI * X[2]);
  } 
  else if (defcoord == 4) {
    dxdxp[1][1] = r-R0;
    dxdxp[2][2] = (210.*M_PI*R*pow(1. - 2.*X[2],2.)*pow(-1. + X[2],2.)*
		pow(X[2],2.) + der0*
		(-32.*pow(-1. + X[2],2.)*pow(X[2],2.)*
		 (3. + 14.*(-1. + X[2])*X[2]) - 
		 R*pow(1. - 2.*X[2],2.)*
		 (-1. + 2.*(-1. + X[2])*X[2]*(2. + 49.*(-1. + X[2])*X[2]))))/R;
  } 
  else if(defcoord == 5) {
    dxdxp[1][1] = r-R0;

    // now assign values
    if(X[2]<0.5){ myx2=X[2]; flip1=0.0; flip2=1.0;}
    else{ myx2=1.0-X[2]; flip1=M_PI; flip2=-1.0;}

    if(myx2<=x2trans){
      dxdxp[2][2] = (3.0*d2*pow(myx2,2.0)+2.0*c2*pow(myx2,1.0)+m2);
    }
    else{
      dxdxp[2][2] = (m3);
    }




  }
  else if (defcoord == 6) {
    dxdxp[1][1] = r-R0;
    dxdxp[2][2] = M_PI;
  }
  else if (defcoord == 7) {

    //dxdxp[1][1] = npow*(r-R0)*pow(log(r-R0),(npow-1.0)/npow);
    dxdxp[1][1] = npow*(r-R0)*pow(X[1],npow-1.0);

    myhslope1 = h0+dmyhslope1dr*(r-rh0);
    myhslope2 = h0+dmyhslope2dx1*(X[1]-x1in);

    dmyhslopedx1=0.5*(dmyhslope1dr*dxdxp[1][1]+dmyhslope2dx1);
    myhslope=0.5*(myhslope1+myhslope2);

    dxdxp[2][2] = M_PI + (1. - myhslope) * M_PI * cos(2. * M_PI * X[2]);
    // d\theta/dx1 not 0
    // d\theta/dx1 = (d\theta/dr)*(dr/dx1)  (is this generally true or just when r(x1)?
    dxdxp[2][1] = -0.5*dmyhslopedx1* sin(2. * M_PI * X[2]);
  }
  else if (defcoord == 8) {
    // drdx1
    dxdxp[1][1] = npow*(r-R0)*pow(X[1],npow-1.0);


    myhslope=2.0-pow(r/r1jet,njet*(-1.0+exp(1.0)/exp(r+rpjet)));

    dmyhslopedr=-((pow(exp(1.0),-r - rpjet)*myhslope*njet*(-exp(1.0) + pow(exp(1.0),r + rpjet) + exp(1.0)*r*log(r/r1jet)))/r);
    dmyhslopedx1=dmyhslopedr*dxdxp[1][1];

    dxdxp[2][2] = M_PI + (1. - myhslope) * M_PI * cos(2. * M_PI * X[2]);
    dxdxp[2][1] = -0.5*dmyhslopedx1* sin(2. * M_PI * X[2]);
  }
  else if(defcoord==9){
    // drdx1
    dxdxp[1][1] = npow*(r-R0)*pow(X[1],npow-1.0);


    myhslope=2.0-Qjet*pow(r/r1jet,-njet*(0.5+1.0/M_PI*atan(r/r0jet-rsjet/r0jet)));

    dmyhslopedr=-((Qjet*(-((njet*(0.5 + atan(r/r0jet - rsjet/r0jet)/M_PI))/r) - (njet*r0jet*log(r/r1jet))/(M_PI*(pow(r0jet,2) + pow(r - rsjet,2)))))/pow(r/r1jet,njet*(0.5 + atan(r/r0jet - rsjet/r0jet)/M_PI)));

    if(!finite(dmyhslopedr)){
      dualfprintf(fail_file,"Problem with dmyhslopedr=%g\n",dmyhslopedr);
      dualfprintf(fail_file,"Qjet=%g njet=%g r=%g rsjet=%g r0jet=%g r1jet=%g\n",Qjet,njet,r,rsjet,r0jet,r1jet);
      myexit(1);
    }

    dmyhslopedx1=dmyhslopedr*dxdxp[1][1];

    dxdxp[2][2] = M_PI + (1. - myhslope) * M_PI * cos(2. * M_PI * X[2]);
    dxdxp[2][1] = -0.5*dmyhslopedx1* sin(2. * M_PI * X[2]);



    

  }
  else if(defcoord==10){
    // drdx1
    dxdxp[1][1] = npow*(r-R0)*pow(X[1],npow-1.0);

    myhslope=(0.5+1.0/M_PI*atan((r-rsjet)/r0jet))*(houter-hinner)+hinner;
    dmyhslopedr=(houter-hinner)*r0jet/(M_PI*(r0jet*r0jet+(r-rsjet)*(r-rsjet)));
    dmyhslopedx1=dmyhslopedr*dxdxp[1][1];

    dxdxp[2][2] = M_PI + (1. - myhslope) * M_PI * cos(2. * M_PI * X[2]);
    dxdxp[2][1] = -0.5*dmyhslopedx1* sin(2. * M_PI * X[2]);






    

  }



  

  else{
    dualfprintf(fail_file,"Shouldn't reach end of dxdxp\n");
    myexit(1);
  }




}



#define GAMMIEDERIVATIVE 0
#define NUMREC 1  // at present this is too slow since dxdxp is called many times.  Could setup permenant memory space, but kinda stupid

// Note that can't use NUMREC for connection if using numerical dxdxp.  Any other form of connection and then any dxdxp can be used.

// For example, if both connection and dxdxp are computed using GAMMIEDERIVATIVE, seems to work fine as a decent numerical approximation

// Also, one can use NUMREC for dxdxp if using analytic connection or GAMMIEDERIVATIVE connection.

//#define DXDERTYPE NUMREC 
#define DXDERTYPE GAMMIEDERIVATIVE

// see conn_func() for notes
#if((REALTYPE==DOUBLETYPE)||(REALTYPE==FLOATTYPE))
#define DXDELTA 1E-5
#elif(REALTYPE==LONGDOUBLETYPE)
#define DXDELTA 1E-6
#endif

void dxdxp_numerical(FTYPE *X, FTYPE (*dxdxp)[NDIM])
{
  int j,k,l;
  FTYPE Xh[NDIM], Xl[NDIM];
  FTYPE coordh[NDIM],coordl[NDIM];
  FTYPE blcoordsimple(FTYPE*X, int i, int j);
  extern FTYPE dfridr(FTYPE (*func)(FTYPE*,int,int), FTYPE *X,int ii, int jj, int kk);

  if(DXDERTYPE==GAMMIEDERIVATIVE){

    for(k=0;k<NDIM;k++){
      for(j=0;j<NDIM;j++){

	if((j==TT)||(k==TT)){
	  if(j!=k) dxdxp[j][k]=0.0;
	  else dxdxp[j][k]=1.0;
	}
	else if((j==PH)||(k==PH)){
	  if(j!=k) dxdxp[j][k]=0.0;
	  else dxdxp[j][k]=1.0;
	}
	else{
	  for(l=0;l<NDIM;l++) Xl[l]=Xh[l]=X[l]; // location of derivative
	  Xh[k]+=DXDELTA; // shift up
	  Xl[k]-=DXDELTA; // shift down
	  // below 2 lines redundant because gets both coordinates, but ok
	  bl_coord(Xh, &coordh[RR], &coordh[TH]);
	  bl_coord(Xl, &coordl[RR], &coordl[TH]);
	  dxdxp[j][k] = (coordh[j] - coordl[j]) / (Xh[k] - Xl[k]);
	}
      }
    }

  }
  else if(DXDERTYPE==NUMREC){

    for(k=0;k<NDIM;k++) for(j=0;j<NDIM;j++){
      dxdxp[j][k]=dfridr(blcoordsimple,X,0,j,k);
    }

  }
}

#undef GAMMIEDERIVATIVE
#undef NUMREC
#undef DXDERTYPE
#undef DXDELTA

FTYPE blcoordsimple(FTYPE*X, int i, int j) // i not used
{
  FTYPE coord[NDIM];
  
  if((j==TT)||(j==PH)) return(X[j]); // dummy linear relationship
  else{
    bl_coord(X, &coord[RR], &coord[TH]);
    return(coord[j]);
  }

}



// /////////////////////////////////////////////////////////////////
// 
// Below set X uniform grid -- usually doesn't change.
// Can usually force startx[1]=startx[2]=0. and dx[1]=1./N1 dx[2]=1./N2
// 
// /////////////////////////////////////////////////////////////////


/* some grid location, dxs */
// could find this by root finding.  Needed if no obvious bounds
// alternatively, could always define grid so x1=0..1 and x2=0..1 (likely more reasonable!)
void set_points()
{

  if (defcoord == 0) {
    startx[1] = log(Rin-R0);
    startx[2] = 0.;
    dx[1] = log((Rout-R0)/(Rin-R0)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 2.0*M_PI;
  } else if (defcoord == 1) {
    startx[1] = 0.;
    startx[2] = 0.;
    dx[1] = 1. / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 2.0*M_PI;
  } else if (defcoord == 2) {
    startx[1] = Rin;
    startx[2] = 0.324;
    dx[1] = (Rout-Rin) / totalsize[1];
    dx[2] = (1.0-2*.324) / totalsize[2];
    dx[3] = 2.0*M_PI;
  }
  else if (defcoord == 3) {
    startx[1] = log(Rin-R0);
    startx[2] = 0.;
    dx[1] = log((Rout-R0)/(Rin-R0)) / totalsize[1];
    dx[2] = 0.5 / totalsize[2];
    dx[3] = 2.0*M_PI;
  }
  else if (defcoord == 4) {
    startx[1] = log(Rin-R0);
    startx[2] = 0.;
    dx[1] = log((Rout-R0)/(Rin-R0)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 2.0*M_PI;
  }
  else if (defcoord == 5) {
    startx[1] = log(Rin-R0);
    startx[2] = 0.;
    dx[1] = log((Rout-R0)/(Rin-R0)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 2.0*M_PI;
  }
  else if (defcoord == 6) {
    startx[1] = log(Rin-R0);
    startx[2] = 0.;
    dx[1] = log((Rout-R0)/(Rin-R0)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 2.0*M_PI;
  } 
  else if (defcoord == 7) {
    startx[1] = pow(log(Rin-R0),1.0/npow);
    startx[2] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 2.0*M_PI;
  } 
  else if (defcoord == 8) {
    startx[1] = pow(log(Rin-R0),1.0/npow);
    startx[2] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 2.0*M_PI;
  } 
  else if (defcoord == 9) {
    startx[1] = pow(log(Rin-R0),1.0/npow);
    startx[2] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 2.0*M_PI;
  } 
  else if (defcoord == 10) {
    startx[1] = pow(log(Rin-R0),1.0/npow);
    startx[2] = 0.;
    dx[1] = (pow(log(Rout-R0),1.0/npow)-pow(log(Rin-R0),1.0/npow)) / totalsize[1];
    dx[2] = 1. / totalsize[2];
    dx[3] = 2.0*M_PI;
  } 
  else if (defcoord == 666) {
    startx[1] = 0.;
    startx[2] = -1.0;
    dx[1] = 1.0 / totalsize[1];
    dx[2] = 2.0 / totalsize[2];
    dx[3] = 2.0*M_PI;
  }
  else{
    dualfprintf(fail_file,"Shouldn't reach end of set_points\n");
    myexit(1);
  }


}

#define MAXIHOR 10
#define FRACN1 (0.1)
#define ADJUSTFRACT (0.25)

int setihor(void)
{
  // set to smaller of either totalsize[1]*0.1 or MAXIHOR
  if(totalsize[1]*FRACN1>MAXIHOR) return((int)MAXIHOR);
  else return((int)((FTYPE)totalsize[1]*(FTYPE)FRACN1));
}


// there's probably a way to do this in general
// probably can root find to get this

// set Rin so horizon exactly on FACE1 at i=ihor
FTYPE setRin(int ihor)
{
 
  FTYPE ftemp;
  FTYPE ihoradjust;


  ihoradjust=((FTYPE)ihor)+ADJUSTFRACT; // can't have grid edge exactly on horizon due to ucon_calc()

  //  fprintf(stderr,"ihoradjust = %21.15g\n",ihoradjust);

  if(defcoord==0){
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord==1){ // even though form appears different in X1, same Rin results
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord==2){ // uniform
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return((Rhor-ftemp*Rout)/(1.0-ftemp));
  }
  else if(defcoord==3){ // as defcoord=0
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord==4){ // as above
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord==5){ // as above
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord==6){ // as above
    ftemp=ihoradjust/(FTYPE)totalsize[1];
    return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
  }
  else if(defcoord==7){
    if(npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else{
      // SUPERGODMARK : need to change for npow>1.0
      return(1.1);
    }
  }
  else if(defcoord==8){
    if(npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else{
      // SUPERGODMARK : need to change for npow>1.0
      return(1.1);
    }
  }
  else if(defcoord==9){
    if(npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else{
      // SUPERGODMARK : need to change for npow>1.0
      return(1.1);
    }
  }
  else if(defcoord==10){
    if(npow==1.0){
      ftemp=ihoradjust/(FTYPE)totalsize[1];
      return(R0+pow((Rhor-R0)/pow(Rout-R0,ftemp),1.0/(1.0-ftemp)));
    }
    else{
      // SUPERGODMARK : need to change for npow>1.0
      return(1.1);
    }
  }
  else{
    dualfprintf(fail_file,"Shouldn't reach end of setRin\n");
    myexit(1);
    return(-1);
  }

}





void coord(int i, int j, int loc, FTYPE *X)
{
  if (loc == FACE1) {
    X[1] = startx[1] + (i + startpos[1]) * dx[1];
    X[2] = startx[2] + ((j + startpos[2]) + 0.5) * dx[2];
  } else if (loc == FACE2) {
    X[1] = startx[1] + ((i + startpos[1]) + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2]) * dx[2];
  } else if (loc == CENT) {
    X[1] = startx[1] + ((i + startpos[1]) + 0.5) * dx[1];
    X[2] = startx[2] + ((j + startpos[2]) + 0.5) * dx[2];
  } else {
    X[1] = startx[1] + (i + startpos[1]) * dx[1];
    X[2] = startx[2] + (j + startpos[2]) * dx[2];
  }

  X[TT]=X[PH]=0.0;

  return;
}


void coordf(FTYPE i, FTYPE j, int loc, FTYPE *X)
{
  if (loc == FACE1) {
    X[1] = startx[1] + (i + startpos[1]) * dx[1];
    X[2] = startx[2] + ((j + startpos[2]) + 0.5) * dx[2];
  } else if (loc == FACE2) {
    X[1] = startx[1] + ((i + startpos[1]) + 0.5) * dx[1];
    X[2] = startx[2] + (j + startpos[2]) * dx[2];
  } else if (loc == CENT) {
    X[1] = startx[1] + ((i + startpos[1]) + 0.5) * dx[1];
    X[2] = startx[2] + ((j + startpos[2]) + 0.5) * dx[2];
  } else {
    X[1] = startx[1] + (i + startpos[1]) * dx[1];
    X[2] = startx[2] + (j + startpos[2]) * dx[2];
  }

  X[TT]=X[PH]=0.0;

  return;
}

void icoord(FTYPE *X,int loc, int *i, int *j)
{
  if(loc == CENT){
    *i = (int)((X[1]-startx[1])/dx[1] - 0.5) ;
    *j = (int)((X[2]-startx[2])/dx[2] - 0.5) ;
  }

  if(startpos[1]+ *i>=totalsize[1]+N1BND) *i=totalsize[1]-1+N1BND;
  if(startpos[2]+ *j>=totalsize[2]+N2BND) *j=totalsize[2]-1+N2BND;
  
  if(startpos[1]+ *i<-N1BND) *i=-N1BND;
  if(startpos[2]+ *j<-N2BND) *j=-N2BND;

}
