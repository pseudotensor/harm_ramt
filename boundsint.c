
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
	for (j = -N2BND; j <=N2-1+N2BND; j++) {
	  ri=0;
	  rj=j;
	  for(i=-1;i>=-N1BND;i--) FLOOP	  	  prim[i][j][k] = prim[ri][rj][k];
	}
      }
    }

    // outer r BC:
    if (mycpupos[1] == ncpux1 - 1) {
      if(BCtype[X1UP]==OUTFLOW){
	/* outer r BC: outflow */
      
	for (j = -N2BND; j <=N2-1+N2BND; j++) {
	  ri=N1-1;
	  rj=j;
	  for(i=N1;i<=N1-1+N1BND;i++)	   FLOOP prim[i][j][k] = prim[ri][rj][k];
	}
      }
      /* if fixed BC: do nothing */
    }

    /* inner polar BC (preserves u^t rho and u) */
    if (mycpupos[2] == 0) {
      for (i = -N1BND; i <=N1-1+N1BND; i++){
	ri=i;
	rj=0;
	for(j=-N2BND;j<=-1;j++)	FLOOP  prim[i][j][k] = prim[ri][rj+(rj-j-1)][k];
      }
    }

    /* outer polar BC  (preserves u^t rho and u) */
    if (mycpupos[2] == ncpux2 - 1) {
      for (i = -N1BND; i <=N1-1+N1BND; i++){
	ri=i;
	rj=N2-1;
	for(j=N2;j<=N2-1+N2BND;j++) FLOOP  prim[i][j][k] = prim[ri][rj+(rj-j+1)][k];
      }
    }
  }

  if (USEMPI) bound_mpi_int(boundstage, prim);

  return (0);
}
