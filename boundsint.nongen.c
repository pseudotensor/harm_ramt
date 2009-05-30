
#include "decs.h"

int bound_pflag(int boundstage, int prim[][N2M][NUMPFLAGS])
{
  int i, j, k;
  int failreturn;
  int ri, rj; // reference i,j



 if((boundstage==STAGE0)||(boundstage==STAGEM1)){
    if (mycpupos[1] == 0) {
      if(BCtype[X1DN]==OUTFLOW){
	/* inner r boundary condition: u, just copy */
	for (j = -2; j <= N2+1; j++) {
	  ri=0;
	  rj=j;
	  for(i=-1;i>=-2;i--) FLOOP	  	  prim[i][j][k] = prim[ri][rj][k];
	}
      }
    }

    // outer r BC:
    if (mycpupos[1] == ncpux1 - 1) {
      if(BCtype[X1UP]==OUTFLOW){
	/* outer r BC: outflow */
      
	for (j = -2; j <= N2+1; j++) {
	  ri=N1-1;
	  rj=j;
	  for(i=N1;i<=N1+1;i++)	   FLOOP prim[i][j][k] = prim[ri][rj][k];
	}
      }
      /* if fixed BC: do nothing */
    }

    /* inner polar BC (preserves u^t rho and u) */
    if (mycpupos[2] == 0) {
      for (i = -2; i <= N1 + 1; i++){
	FLOOP{
	  prim[i][-1][k] = prim[i][0][k];
	  prim[i][-2][k] = prim[i][1][k];
	}
      }
    }

    /* outer polar BC  (preserves u^t rho and u) */
    if (mycpupos[2] == ncpux2 - 1) {
      for (i = -2; i <= N1 + 1; i++){
	FLOOP{
	  prim[i][N2][k] = prim[i][N2 - 1][k];
	  prim[i][N2 + 1][k] = prim[i][N2 - 2][k];
	}
      }
    }
  }

  if (USEMPI) bound_mpi_int(boundstage, prim);

  return (0);
}
