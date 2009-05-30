
#include "decs.h"

/* bound array containing entire set of primitive variables */

#define EXTRAP 2
// 0: just copy
// 1: gdet or other extrapolation
// 2: copy (with rescale())

int bound_prim(int boundstage, FTYPE prim[][N2M][NPR])
{
  int i, j, k;
  struct of_geom geom,rgeom;
#if(WHICHVEL==VEL3)
  int failreturn;
#endif
  int ri, rj; // reference i,j



  // real boundary zones
  if((boundstage==STAGE0)||(boundstage==STAGEM1)){

    if (mycpupos[1] == 0) {
      if(BCtype[X1DN]==OUTFLOW){
	/* inner r boundary condition: u, just copy */
	for (j = -2; j <= N2+1; j++) {
#if(EXTRAP==0)
	  ri=0;
	  rj=j;
	  for(i=-1;i>=-2;i--)	  PBOUNDLOOP	  prim[i][j][k] = prim[ri][rj][k];
#elif(EXTRAP==1)
	  prim[-1][j][RHO] = prim[0][j][RHO] *
	    gdet[0][j][CENT]/gdet[-1][j][CENT] ;
	  prim[-2][j][RHO] = prim[0][j][RHO] *
	    gdet[0][j][CENT]/gdet[-2][j][CENT] ;
	
	  prim[-1][j][UU] = prim[0][j][UU] *
	    gdet[0][j][CENT]/gdet[-1][j][CENT] ;
	  prim[-2][j][UU] = prim[0][j][UU] *
	    gdet[0][j][CENT]/gdet[-2][j][CENT] ;
	
	  prim[-1][j][U1] = prim[0][j][U1] * (1. + 1.*dx[1]) ;
	  prim[-2][j][U1] = prim[0][j][U1] * (1. + 2.*dx[1]) ;
	  /*
	    prim[-1][j][U1] = prim[0][j][U1] ;
	    prim[-2][j][U1] = prim[0][j][U1] ;
	  */
	
	  prim[-1][j][U2] = prim[0][j][U2] * (1. - 1.*dx[1]) ;
	  prim[-2][j][U2] = prim[0][j][U2] * (1. - 2.*dx[1]) ;
	
	  prim[-1][j][U3] = prim[0][j][U3] * (1. - 1.*dx[1]) ;
	  prim[-2][j][U3] = prim[0][j][U3] * (1. - 2.*dx[1]) ;
	
	  prim[-1][j][B1] = prim[0][j][B1] *
	    gdet[0][j][CENT]/gdet[-1][j][CENT] ;
	  prim[-2][j][B1] = prim[0][j][B1] *
	    gdet[0][j][CENT]/gdet[-2][j][CENT] ;
	
	  prim[-1][j][B2] = prim[0][j][B2] * (1. - 1.*dx[1]) ;
	  prim[-2][j][B2] = prim[0][j][B2] * (1. - 2.*dx[1]) ;
	
	  prim[-1][j][B3] = prim[0][j][B3] * (1. - 1.*dx[1]) ;
	  prim[-2][j][B3] = prim[0][j][B3] * (1. - 2.*dx[1]) ;
#elif(EXTRAP==2)
	  ri=0;
	  rj=j;
	  get_geometry(ri, rj, CENT, &rgeom);
	  rescale(1,1,prim[ri][rj],&rgeom,prim[ri][rj]);
	  for(i=-1;i>=-2;i--){
	    PBOUNDLOOP	  prim[i][j][k] = prim[ri][rj][k];
	    get_geometry(i, j, CENT, &geom);
	    rescale(-1,1,prim[i][j],&geom,prim[i][j]);
	  }
	  rescale(-1,1,prim[ri][rj],&rgeom,prim[ri][rj]);	
#endif
	  for(i=-1;i>=-2;i--){
#if(WHICHVEL==VEL4)
	    get_geometry(i, j, CENT, &geom);
	    inflow_check_4vel(1,prim[i][j],&geom) ;
#elif(WHICHVEL==VEL3)
	    get_geometry(i, j, CENT, &geom);
	    inflow_check_3vel(1,prim[i][j],&geom) ;
	    // projection may not preserve u^t to be real and rho>rhoscal u>uuscal
	    if(jonchecks){
	      //fixup1zone(prim[i][j],&geom,0);
	      failreturn=check_pr(prim[i][j],prim[i][j],&geom,-3);
	      if(failreturn){
		dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d\n",startpos[1]+i,startpos[2]+j);
		if (fail(FAIL_BCFIX) >= 1) return (1);
	      }
	    }
#elif(WHICHVEL==VELREL4)
	    get_geometry(i,j,CENT,&geom) ;
	    inflow_check_rel4vel(1,prim[i][j],&geom) ;
	    if(limit_gamma(GAMMAMAX,prim[i][j],&geom)>=1)
	      FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 1);
#endif	
	  }
	}
      }
    }

    // outer r BC:
    if (mycpupos[1] == ncpux1 - 1) {
      if(BCtype[X1UP]==OUTFLOW){
	/* outer r BC: outflow */
      
	for (j = -2; j <= N2+1; j++) {
#if(EXTRAP==0)
	  ri=N1-1;
	  rj=j;
	  for(i=N1;i<=N1+1;i++)	  PBOUNDLOOP prim[i][j][k] = prim[ri][rj][k];
#elif(EXTRAP==1)
	  prim[N1][j][RHO] = prim[N1-1][j][RHO] *
	    gdet[N1-1][j][CENT]/gdet[N1][j][CENT] ;
	  prim[N1+1][j][RHO] = prim[N1-1][j][RHO] *
	    gdet[N1-1][j][CENT]/gdet[N1+1][j][CENT] ;

	  prim[N1][j][UU] = prim[N1-1][j][UU] *
	    gdet[N1-1][j][CENT]/gdet[N1][j][CENT] ;
	  prim[N1+1][j][UU] = prim[N1-1][j][UU] *
	    gdet[N1-1][j][CENT]/gdet[N1+1][j][CENT] ;

	  prim[N1][j][U1] = prim[N1-1][j][U1] * (1. - 2.*dx[1]) ;
	  prim[N1+1][j][U1] = prim[N1-1][j][U1] * (1. - 4.*dx[1]) ;

	  prim[N1][j][U2] = prim[N1-1][j][U2] * (1. - 1.*dx[1]) ;
	  prim[N1+1][j][U2] = prim[N1-1][j][U2] * (1. - 2.*dx[1]) ;

	  prim[N1][j][U3] = prim[N1-1][j][U3] * (1. - 1.*dx[1]) ;
	  prim[N1+1][j][U3] = prim[N1-1][j][U3] * (1. - 2.*dx[1]) ;

	  prim[N1][j][B1] = prim[N1-1][j][B1] *
	    gdet[N1-1][j][CENT]/gdet[N1][j][CENT] ;
	  prim[N1+1][j][B1] = prim[N1-1][j][B1] *
	    gdet[N1-1][j][CENT]/gdet[N1+1][j][CENT] ;

	  prim[N1][j][B2] = prim[N1-1][j][B2] * (1. - 1.*dx[1]) ;
	  prim[N1+1][j][B2] = prim[N1-1][j][B2] * (1. - 2.*dx[1]) ;

	  prim[N1][j][B3] = prim[N1-1][j][B3] * (1. - 1.*dx[1]) ;
	  prim[N1+1][j][B3] = prim[N1-1][j][B3] * (1. - 2.*dx[1]) ;
#elif(EXTRAP==2)
	  ri=N1-1;
	  rj=j;
	  get_geometry(ri, rj, CENT, &rgeom);
	  rescale(1,1,prim[ri][rj],&rgeom,prim[ri][rj]);
	  for(i=N1;i<=N1+1;i++){
	    PBOUNDLOOP prim[i][j][k] = prim[ri][rj][k];
	    get_geometry(i, j, CENT, &geom);
	    rescale(-1,1,prim[i][j],&geom,prim[i][j]);
	  }
	  rescale(-1,1,prim[ri][rj],&rgeom,prim[ri][rj]);
#endif

	  for(i=N1;i<=N1+1;i++){
#if(WHICHVEL==VEL4)
	    get_geometry(i, j, CENT, &geom);
	    inflow_check_4vel(1,prim[i][j],&geom) ;
#elif(WHICHVEL==VEL3)
	    get_geometry(i, j, CENT, &geom);
	    inflow_check_3vel(1,prim[i][j],&geom) ;
	    // projection may not preserve u^t to be real and rho>rhoscal u>uuscal
	    if(jonchecks){
	      //fixup1zone(prim[i][j],&geom,0);
	      failreturn=check_pr(prim[i][j],prim[i][j],&geom,-3);
	      if(failreturn){
		dualfprintf(fail_file,"Bad boundary zone, couldn't fix: i=%d j=%d\n",startpos[1]+i,startpos[2]+j);
		if (fail(FAIL_BCFIX) >= 1) return (1);
	      }
	    }
#elif(WHICHVEL==VELREL4)
	    get_geometry(i,j,CENT,&geom) ;
	    inflow_check_rel4vel(1,prim[i][j],&geom) ;
	    if(limit_gamma(GAMMAMAX,prim[i][j],&geom)>=1)
	      FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 2);
#endif	
	  }
	}
      }
      /* if fixed BC: do nothing */
    }

    /* inner polar BC (preserves u^t rho and u) */
    if (mycpupos[2] == 0) {
      for (i = -2; i <= N1 + 1; i++){
	PBOUNDLOOP {
	  prim[i][-1][k] = prim[i][0][k];
	  prim[i][-2][k] = prim[i][1][k];
	}
      }
    }

    /* outer polar BC  (preserves u^t rho and u) */
    if (mycpupos[2] == ncpux2 - 1) {
      for (i = -2; i <= N1 + 1; i++){
	PBOUNDLOOP {
	  prim[i][N2][k] = prim[i][N2 - 1][k];
	  prim[i][N2 + 1][k] = prim[i][N2 - 2][k];
	}
      }
    }
    /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
    /* inner pole */
    if (mycpupos[2] == 0) {
      for (i = -2; i < N1 + 2; i++) {
	for (j = -2; j < 0; j++) {
	  if(POSDEFMETRIC==0){
	    // u^t must be symmetric across pole, which is functions of u2 and u3 as well as their squares and othe products.  u2 in KS happens to be independent of sign, but in general is could be for some other metric.
	    // for now, assume KS-like metric where u2 is antisymmetric and u^t dep only on u2^2, not u2
	    prim[i][j][U2] *= -1.;
	    prim[i][j][U3] *= 1.;
	    prim[i][j][B2] *= -1.;
	    prim[i][j][B3] *= 1.;
	  }
	  else{
	    prim[i][j][U2] *= -1.;
	    prim[i][j][U3] *= 1.;
	    prim[i][j][B2] *= -1.;
	    prim[i][j][B3] *= 1.;
	  }
	}
      }
    }
    /* outer pole */
    if (mycpupos[2] == ncpux2 - 1) {
      for (i = -2; i < N1 + 2; i++) {
	for (j = N2; j < N2 + 2; j++) {
	  if(POSDEFMETRIC==0){
	    prim[i][j][U2] *= -1.;
	    prim[i][j][U3] *= 1.;
	    prim[i][j][B2] *= -1.;
	    prim[i][j][B3] *= 1.;
	  }
	  else{
	    prim[i][j][U2] *= -1.;
	    prim[i][j][U3] *= 1.;
	    prim[i][j][B2] *= -1.;
	    prim[i][j][B3] *= 1.;
	  }
	}
      }
    }
  }

  if (USEMPI) bound_mpi(boundstage, prim);


  return (0);
}
