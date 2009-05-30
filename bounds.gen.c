
#include "decs.h"

/* bound array containing entire set of primitive variables */

#define EXTRAP 2
// 0: just copy
// 1: gdet or other extrapolation (no longer supported here, use old bounds.c and boundsint.c)
// 2: copy (with rescale())

int bound_prim(int boundstage, FTYPE prim[][N2M][NPR])
{
  int i, j, k;
  int is,ie,js,je;
  struct of_geom geom,rgeom;
#if(WHICHVEL==VEL3)
  int failreturn;
#endif
  int ri, rj; // reference i,j
  int dir;
  int ref;
  
  // real boundary zones
  if((boundstage==STAGE0)||(boundstage==STAGEM1)){

    DIRLOOP{
      if((mycpupos[1]==0)&&(dir==X1DN)){
	is=-2;
	ie=-1;
	js=-2;
	je=N2+1;
	ref=1;
      }
      else if((mycpupos[1]==ncpux1-1)&&(dir==X1UP)){
	is=N1;
	ie=N1+1;
	js=-2;
	je=N2+1;
	ref=1;
      }
      else if((mycpupos[2]==0)&&(dir==X2DN)){
	is=-2;
	ie=N1+1;
	js=-2;
	je=-1;
	ref=2;
      }
      else if((mycpupos[2]==ncpux2-1)&&(dir==X2UP)){
	is=-2;
	ie=N1+1;
	js=N2;
	je=N2+1;
	ref=2;
      }
      else continue;

      if(BCtype[dir]==OUTFLOW){
	/* inner r boundary condition: u, just copy */
	ZSLOOP(is,ie,js,je){
	  if(dir==X1DN){
	    ri=0;
	    rj=j;
	  } else if(dir==X1UP){
	    ri=N1-1;
	    rj=j;
	  } else if(dir==X2DN){
	    ri=i;
	    rj=0;
	  } else if(dir==X2UP){
	    ri=i;
	    rj=N2-1;
	  }
#if(EXTRAP==0)
	  PLOOP	  prim[i][j][k] = prim[ri][rj][k];
#elif(EXTRAP==2)
	  get_geometry(ri, rj, CENT, &rgeom);
	  rescale(1,ref,prim[ri][rj],&rgeom,prim[ri][rj]);
	  PLOOP	  prim[i][j][k] = prim[ri][rj][k];
	  get_geometry(i, j, CENT, &geom);
	  rescale(-1,ref,prim[i][j],&geom,prim[i][j]);
	  rescale(-1,ref,prim[ri][rj],&rgeom,prim[ri][rj]);	
#endif
#if(WHICHVEL==VEL4)
	  get_geometry(i, j, CENT, &geom);
	  inflow_check_4vel(ref,prim[i][j],&geom) ;
#elif(WHICHVEL==VEL3)
	  get_geometry(i, j, CENT, &geom);
	  inflow_check_3vel(ref,prim[i][j],&geom) ;
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
	  inflow_check_rel4vel(ref,prim[i][j],&geom) ;
	  if(limit_gamma(GAMMAMAX,prim[i][j],&geom)>=1)
	    FAILSTATEMENT("bounds.c:bound_prim()", "limit_gamma()", 1);
#endif	
	}
      }
      else if(BCtype[dir]==POLARAXIS){
	ZSLOOP(is,ie,js,je){
	  if(dir==X1DN){
	    if(j==-1){
	      ri=0;
	      rj=j;
	    }
	    else if(j==-2){
	      ri=1;
	      rj=j;
	    }
	  }
	  else if(dir==X1UP){
	    if(j==N2){
	      ri=N1-1;
	      rj=j;
	    }
	    else if(j==N2+1){
	      ri=N1-2;
	      rj=j;
	    }
	  } else if(dir==X2DN){
	    if(i==-1){
	      ri=i;
	      rj=0;
	    }
	    else if(i==-2){
	      ri=i;
	      rj=1;
	    }
	  } else if(dir==X2UP){
	    if(i==N1){
	      ri=i;
	      rj=N2-1;
	    }
	    else if(i==N1+1){
	      ri=i;
	      rj=N2-2;
	    }
	  }
	
	  /* inner polar BC (preserves u^t rho and u) */
	  PLOOP prim[i][j][k] = prim[ri][rj][k];
	  if(POSDEFMETRIC==0){
	    /* make sure b and u are antisymmetric at the poles   (preserves u^t rho and u) */
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
  }

  if (USEMPI) bound_mpi(boundstage, prim);


  return (0);
}
