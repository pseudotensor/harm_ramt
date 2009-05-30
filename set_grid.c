
#include "decs.h"


/* set up all grid functions */
void set_grid()
{
  int i, j, k, l, m;
  int ii, jj;
  FTYPE X[NDIM];
  struct of_geom geom;
  int loc;
  FTYPE r,th;
  extern void eomfunc_func(int getprim, int whichcoord, FTYPE *X, FTYPE *eomfunc);

  /* set up boundaries, steps in coordinate grid */
  set_points();
  dV = dx[1] * dx[2]; // computational volume
  dVF = dV * dx[3] ; // full 3d volume (used for diagnostics only)

  DLOOPA X[j] = 0.;

  FULLLOOP {

    // over grid locations needing these quantities
    for (loc = NUMGRIDPOS - 1; loc >= 0; loc--) {

      coord(i, j, loc, X);

      /////////////////
      //
      // (1,MCOORD) here actually means PRIMCOORDS since the "1" means convert MCOORD to PRIMCOORDS.

      //      dxdxp_func(X,dxdxp[i][j][loc]); // future numerical version
      gcov_func(1,MCOORD,X, gcov[i][j][loc]);
      gdet[i][j][loc] = gdet_func(gcov[i][j][loc]);
      gcon_func(1,MCOORD,X,gcov[i][j][loc],gcon[i][j][loc]);
      eomfunc_func(1,MCOORD,X,&eomfunc[i][j][loc]);

      // check if near static limit since can't divide by the below in ucon_calc
      // GODMARK
      if (fabs(gcon[i][j][loc][TT][TT]) < SMALL) {
	bl_coord(X,&r,&th);
	dualfprintf(fail_file, "grid location too near g_{tt}==0: %d %d : r=%21.15g th=%21.15g : Rin=%21.15g %21.15g\n", i, j,r,th,Rin,gcon[i][j][loc][TT][TT]);
	myexit(1);
      }
      if (fabs(gcon[i][j][loc][RR][RR]) < SMALL) {
	bl_coord(X,&r,&th);
	dualfprintf(fail_file, "grid location too near g^{rr}==0: %d %d : r=%21.15g th=%21.15g :  Rin=%21.15g %21.15g\n", i, j,r,th,Rin,gcon[i][j][loc][RR][RR]);
	myexit(1);
      }
      if (fabs(gcon[i][j][loc][TH][TH]) < SMALL) {
	bl_coord(X,&r,&th);
	dualfprintf(fail_file,"grid location too near g^{\\theta\\theta}==0: %d %d : r=%21.15g th=%21.15g :  Rin=%21.15g %21.15g\n", i, j,r,th,Rin,gcon[i][j][loc][TH][TH]);
	myexit(1);
      }
      if (fabs(gcon[i][j][loc][PH][PH]) < SMALL) {
	bl_coord(X,&r,&th);
	dualfprintf(fail_file,"grid location too near g^{\\phi\\phi}==0: %d %d : r=%21.15g th=%21.15g :  Rin=%21.15g %21.15g\n", i, j,r,th,Rin,gcon[i][j][loc][PH][PH]);
	myexit(1);
      }
      // what about g_{tt}==0? Do I ever divide by g_{tt}?
      // Yes, for ucon[TT] for 4 velocity, which is done if WHICHVEL==VEL4 or init.c
      // what about g^{rr}==0? Do I ever divide by g^{rr}?
      // Doesn't appear so
      // g^{pp} inf taken care of in metric.c by avoiding theta=0,Pi
      if (loc == CENT) {
	get_geometry(i, j, loc, &geom);
	conn_func(MCOORD,X, &geom, conn[i][j],conn2[i][j]);
      }
    }
  }
  if(VOLUMEDIFF){
    ZLOOP{
      // only at 1 location, centered, using surrounding edge values
      if((defcoord==6)&&(MCOORD==KSCOORDS)){
	mks_unitheta_idxvol_func(i,j,idxvol[i][j]);
      }
      else{
	idxvol[i][j][TT]=1.0; // really 1/dt, but changes in time      
	idxvol[i][j][RR]=1.0/dx[1];
	idxvol[i][j][TH]=1.0/dx[2];
	idxvol[i][j][PH]=1.0/dx[3];
      }
    }
  }

  /* done! */
}
