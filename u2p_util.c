
//#include "u2p_defs.h"


/* void primtoU_g( FTYPE prim[], FTYPE gcov[][4], FTYPE gcon[][4],  FTYPE U[] ); */
/* void ucon_calc_g(FTYPE prim[],FTYPE gcov[][4],FTYPE gcon[][4],FTYPE ucon[]); */
/* void raise_g(FTYPE vcov[], FTYPE gcon[][4], FTYPE vcon[]); */
/* void lower_g(FTYPE vcon[], FTYPE gcov[][4], FTYPE vcov[]); */
/* void ncov_calc(FTYPE gcon[][4],FTYPE ncov[]) ; */
/* void bcon_calc_g(FTYPE prim[],FTYPE ucon[],FTYPE ucov[],FTYPE ncov[],FTYPE bcon[]);  */
/* FTYPE pressure_rho0_u(FTYPE rho0, FTYPE u); */
/* FTYPE pressure_rho0_w(FTYPE rho0, FTYPE w); */
/* void dualfprintf(FILE* fileptr, char *format, ...); */

/* 
 * reference implementation of transformation from
 * primitive to conserved variables.
 *
 * cfg 7-6-04
 *
 * input: 
 * 
 * primitive variables in 8 element array 
 * metric in contravariant and covariant form
 *
 * output:
 * 
 * conserved variables in 8 element array 
 * 
 */

static void primtoU_g(
	FTYPE prim[8],       /* primitive variables */
	FTYPE gcov[4][4],    /* covariant (index dn) form of metric */
	FTYPE gcon[4][4],    /* contravariant (index up) form of metric */
	FTYPE U[8]           /* matrix of derivatives */
) {
	int i,j ;
	FTYPE rho0 ;
	FTYPE ucon[4],ucov[4],bcon[4],bcov[4],ncov[4] ;
	FTYPE gamma,n_dot_b,bsq,u,p,w ;
	static void ucon_calc_g(FTYPE prim[],FTYPE gcov[][4],FTYPE gcon[][4],FTYPE ucon[]);
	static void raise_g(FTYPE vcov[], FTYPE gcon[][4], FTYPE vcon[]);
	static void lower_g(FTYPE vcon[], FTYPE gcov[][4], FTYPE vcov[]);
	static void ncov_calc(FTYPE gcon[][4],FTYPE ncov[]) ;
	static FTYPE pressure_rho0_u(FTYPE rho0, FTYPE u);
	static FTYPE pressure_rho0_w(FTYPE rho0, FTYPE w);
	static void bcon_calc_g(FTYPE prim[],FTYPE ucon[],FTYPE ucov[],FTYPE ncov[],FTYPE bcon[]); 

	/* preliminaries */
	ucon_calc_g(prim,gcov,gcon,ucon) ;
	lower_g(ucon,gcov,ucov) ;
	ncov_calc(gcon,ncov) ;

	gamma = -ncov[0]*ucon[0] ;

	bcon_calc_g(prim,ucon,ucov,ncov,bcon) ;
	lower_g(bcon,gcov,bcov) ;

   	n_dot_b = 0. ;
	for(i=0;i<4;i++) n_dot_b += ncov[i]*bcon[i] ;
	bsq = 0. ;
	for(i=0;i<4;i++) bsq += bcov[i]*bcon[i] ;

	rho0 = prim[RHO] ;
	u = prim[UU] ;
	p = pressure_rho0_u(rho0,u) ;
	w = rho0 + u + p ;

	U[RHO] = gamma*rho0 ;

	for(i=0;i<4;i++) 
		U[QCOV0+i] = gamma*(w + bsq)*ucov[i] 
			- (p + bsq/2.)*ncov[i] 
			+ n_dot_b*bcov[i] ;

	U[BCON1] = prim[BCON1] ;
	U[BCON2] = prim[BCON2] ;
	U[BCON3] = prim[BCON3] ;

	return ;
}

/* find the contravariant fluid four-velocity from primitive 
   variables plus the metric */
static void ucon_calc_g(FTYPE prim[8],FTYPE gcov[4][4],FTYPE gcon[4][4],FTYPE ucon[4])
{
	FTYPE u_tilde_con[4] ;
	FTYPE u_tilde_sq ;
	FTYPE gamma,lapse ;
	int i,j ;
	
	u_tilde_con[0] = 0. ;
	u_tilde_con[1] = prim[UTCON1] ;
	u_tilde_con[2] = prim[UTCON2] ;
	u_tilde_con[3] = prim[UTCON3] ;

	u_tilde_sq = 0. ;
	for(i=0;i<4;i++)
	for(j=0;j<4;j++)
		u_tilde_sq += gcov[i][j]*u_tilde_con[i]*u_tilde_con[j] ;
	u_tilde_sq = fabs(u_tilde_sq) ;

	gamma = sqrt(1. + u_tilde_sq) ;

	lapse = sqrt(-1./gcon[0][0]) ;

	for(i=0;i<4;i++) ucon[i] = u_tilde_con[i] - lapse*gamma*gcon[0][i] ;

	return ;
}

/* raise covariant vector vcov using gcon, place result in vcon */
static void raise_g(FTYPE vcov[4], FTYPE gcon[4][4], FTYPE vcon[4])
{
	int i,j;

	for(i=0;i<4;i++) {
		vcon[i] = 0. ;
		for(j=0;j<4;j++) 
			vcon[i] += gcon[i][j]*vcov[j] ;
	}

	return ;
}
/* lower contravariant vector vcon using gcov, place result in vcov */
static void lower_g(FTYPE vcon[4], FTYPE gcov[4][4], FTYPE vcov[4])
{
	int i,j;

	for(i=0;i<4;i++) {
		vcov[i] = 0. ;
		for(j=0;j<4;j++) 
			vcov[i] += gcov[i][j]*vcon[j] ;
	}

	return ;
}

/* set covariant normal observer four-velocity */
static void ncov_calc(FTYPE gcon[4][4],FTYPE ncov[4]) 
{
	FTYPE lapse ;

	lapse = sqrt(-1./gcon[0][0]) ;

	ncov[0] = -lapse ;
	ncov[1] = 0. ;
	ncov[2] = 0. ;
	ncov[3] = 0. ;

	return ;
}

/* calculate contravariant magnetic field four-vector b */
static void bcon_calc_g(FTYPE prim[8],FTYPE ucon[4],FTYPE ucov[4],FTYPE ncov[4],FTYPE bcon[4]) 
{
	FTYPE Bcon[4] ;
	FTYPE u_dot_B ;
	FTYPE gamma ;
	int i ;

	Bcon[0] = 0. ;
	for(i=1;i<4;i++) Bcon[i] = prim[BCON1+i-1] ;

	u_dot_B = 0. ;
	for(i=0;i<4;i++) u_dot_B += ucov[i]*Bcon[i] ;

	gamma = -ucon[0]*ncov[0] ;
	for(i=0;i<4;i++) bcon[i] = (Bcon[i] + ucon[i]*u_dot_B)/gamma ;
}


/* 

pressure as a function of rho0 and u 

this is used by primtoU and Utoprim_?D

*/
static FTYPE pressure_rho0_u(FTYPE rho0, FTYPE u)
{
	return((GAMMA - 1.)*u) ;
}


  
/* 

pressure as a function of rho0 and w = rho0 + u + p 

this is used by primtoU and Utoprim_1D

*/
static FTYPE pressure_rho0_w(FTYPE rho0, FTYPE w)
{
	return((GAMMA-1.)*(w - rho0)/GAMMA) ;
}


/* void dualfprintf(FILE* fileptr, char *format, ...) */
/* { */
/*   va_list arglist; */
 
/*   va_start (arglist, format); */
 
/*   vfprintf (stderr, format, arglist); */
/*   fflush(stderr); */

/*   va_end (arglist); */
/* } */
