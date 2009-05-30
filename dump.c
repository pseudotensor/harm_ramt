#include "decs.h"

// mpi.h has following datatypes corresponding to the C types
// pick one per dump file. (no per column types yet)
// same for image.c
// same for restart.c
/*
#define MPI_CHAR           ((MPI_Datatype)1)
#define MPI_UNSIGNED_CHAR  ((MPI_Datatype)2)
#define MPI_BYTE           ((MPI_Datatype)3)
#define MPI_SHORT          ((MPI_Datatype)4)
#define MPI_UNSIGNED_SHORT ((MPI_Datatype)5)
#define MPI_INT            ((MPI_Datatype)6)
#define MPI_UNSIGNED       ((MPI_Datatype)7)
#define MPI_LONG           ((MPI_Datatype)8)
#define MPI_UNSIGNED_LONG  ((MPI_Datatype)9)
#define MPI_FLOAT          ((MPI_Datatype)10)
#define MPI_DOUBLE         ((MPI_Datatype)11)
#define MPI_LONG_DOUBLE    ((MPI_Datatype)12)
#define MPI_LONG_LONG_INT  ((MPI_Datatype)13)
*/

int dump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];


  trifprintf("begin dumping dump# %ld ... ",dump_cnt);

  whichdump=DUMPCOL;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/dump");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,dump_content)>=1) return(1);

  trifprintf("end dumping dump# %ld ... ",dump_cnt);


  return(0);
  
}

int dump_header(int bintxt, FILE *headerptr)
{
  // 15 elements total
  if(bintxt==BINARYOUTPUT){
    fwrite(&t,sizeof(FTYPE),1,headerptr);
    fwrite(&totalsize[1],sizeof(int),1,headerptr);
    fwrite(&totalsize[2],sizeof(int),1,headerptr);
    fwrite(&startx[1],sizeof(FTYPE),1,headerptr);
    fwrite(&startx[2],sizeof(FTYPE),1,headerptr);
    fwrite(&dx[1],sizeof(FTYPE),1,headerptr);
    fwrite(&dx[2],sizeof(FTYPE),1,headerptr);
    fwrite(&realnstep,sizeof(long),1,headerptr);
    fwrite(&gam,sizeof(FTYPE),1,headerptr);
    fwrite(&a,sizeof(FTYPE),1,headerptr);
    fwrite(&R0,sizeof(FTYPE),1,headerptr);
    fwrite(&Rin,sizeof(FTYPE),1,headerptr);
    fwrite(&Rout,sizeof(FTYPE),1,headerptr);
    fwrite(&hslope,sizeof(FTYPE),1,headerptr);
    fwrite(&dt,sizeof(FTYPE),1,headerptr);
    fwrite(&defcoord,sizeof(int),1,headerptr);
  }
  else{
#if(REALTYPE==DOUBLETYPE)
    //    fprintf(headerptr, "%21.15g %d %d %21.15g %21.15g %21.15g %21.15g %ld %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %d\n", t, totalsize[1]+2*N1BND, totalsize[2]+2*N2BND, startx[1], startx[2], dx[1], dx[2],realnstep,gam,a,R0,Rin,Rout,hslope,dt,defcoord);
    fprintf(headerptr, "%21.15g %d %d %21.15g %21.15g %21.15g %21.15g %ld %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %d\n", t, totalsize[1], totalsize[2], startx[1], startx[2], dx[1], dx[2],realnstep,gam,a,R0,Rin,Rout,hslope,dt,defcoord);
#elif(REALTYPE==LONGDOUBLETYPE)
    fprintf(headerptr, "%31.25Lg %d %d %31.25Lg %31.25Lg %31.25Lg %31.25Lg %ld %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %31.25Lg %d\n", t, totalsize[1], totalsize[2], startx[1], startx[2], dx[1], dx[2],realnstep,gam,a,R0,Rin,Rout,hslope,dt,defcoord);
#endif
  }
  fflush(headerptr);
  return(0);
}	

int dump_content(int i, int j, MPI_Datatype datatype,void *writebuf)
{
  int k = 0;
  FTYPE r, th, vmin1, vmax1, vmin2, vmax2;
  int ignorecourant;
  struct of_geom geom;
  struct of_state q;
  FTYPE X[NDIM];
  FTYPE divb;
  FTYPE b[NDIM],ucon[NDIM];
  FTYPE U[NPR];
  FTYPE ftemp;
  FTYPE jcov[NDIM];
  FTYPE fcov[NUMFARADAY];


  //////////////
  //
  // some calculations
  //

  coord(i, j, CENT, X);
  bl_coord(X, &r, &th);
  // if failed, then data output for below invalid, but columns still must exist    
  get_geometry(i, j, CENT, &geom);
  if (!failed) {
    if (get_state(p[i][j], &geom, &q) >= 1)
      FAILSTATEMENT("dump.c:dump()", "get_state() dir=0", 1);
    if (vchar(p[i][j], &q, 1, &geom, &vmax1, &vmin1,&ignorecourant) >= 1)
      FAILSTATEMENT("dump.c:dump()", "vchar() dir=1or2", 1);
    if (vchar(p[i][j], &q, 2, &geom, &vmax2, &vmin2,&ignorecourant) >= 1)
      FAILSTATEMENT("dump.c:dump()", "vchar() dir=1or2", 2);
  }
  else {// do a per zone check, otherwise set to 0
    whocalleducon=1; // force no failure mode, just return like failure, and don't return if failure, just set to 0 and continue
    if (get_state(p[i][j], &geom, &q) >= 1){
      for (k = 0; k < NDIM; k++)
	q.ucon[k]=0;
      for (k = 0; k < NDIM; k++)
	q.ucov[k]=0;
      for (k = 0; k < NDIM; k++)
	q.bcon[k]=0;
      for (k = 0; k < NDIM; k++)
	q.bcov[k]=0;
    }
    if (vchar(p[i][j], &q, 1, &geom, &vmax1, &vmin1,&ignorecourant) >= 1){
      vmax1=vmin1=0;
    }
    
    if (vchar(p[i][j], &q, 2, &geom, &vmax2, &vmin2,&ignorecourant) >= 1){
      vmax2=vmin2=0;
    }
    whocalleducon=0; // return to normal state
    
  }
  SETFDIVB(divb, p, i, j);


  //////////////////////////
  //
  // do the assignments
  //
  // if you change # of outputted vars, remember to change numcolumns


  //static
  if(!GAMMIEDUMP){
    ftemp=(FTYPE)(i+startpos[1]);
    myset(datatype,&ftemp,0,1,writebuf);
    ftemp=(FTYPE)(j+startpos[2]);
    myset(datatype,&ftemp,0,1,writebuf);
  }
  myset(datatype,X,1,2,writebuf);
  myset(datatype,&r,0,1,writebuf);
  myset(datatype,&th,0,1,writebuf); // 6

  // rest dynamic
  myset(datatype,p[i][j],0,NPRDUMP,writebuf);
  myset(datatype,&divb,0,1,writebuf); // 15

  for (k = 0; k < NDIM; k++)
    myset(datatype,&(q.ucon[k]),0,1,writebuf);
  for (k = 0; k < NDIM; k++)
    myset(datatype,&(q.ucov[k]),0,1,writebuf);
  for (k = 0; k < NDIM; k++)
    myset(datatype,&(q.bcon[k]),0,1,writebuf);
  for (k = 0; k < NDIM; k++)
    myset(datatype,&(q.bcov[k]),0,1,writebuf);
    
  myset(datatype,&vmin1,0,1,writebuf);
  myset(datatype,&vmax1,0,1,writebuf);
  myset(datatype,&vmin2,0,1,writebuf);
  myset(datatype,&vmax2,0,1,writebuf);

  // one static term
  myset(datatype,&geom.g,0,1,writebuf);

  // updated 11/16/2003
  // new 10/23/2003
  // current density (2*NDIM)
  lower(jcon[i][j],&geom,jcov);
  myset(datatype,jcon[i][j],0,NDIM,writebuf);
  myset(datatype,jcov,0,NDIM,writebuf);
  // faraday (2*6)
  lowerf(fcon[i][j],&geom,fcov);
  myset(datatype,fcon[i][j],0,NUMFARADAY,writebuf);
  myset(datatype,fcov,0,NUMFARADAY,writebuf);

  return (0);
}






int debugdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];
  int i,j, floor;


  trifprintf("begin dumping debug dump# %ld ... ",dump_cnt);

  whichdump=DEBUGCOL;
  datatype=MPI_CTYPE;
  strcpy(fileprefix,"dumps/debug");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");
  
  // same header as dump
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,debug_content)>=1) return(1);

  trifprintf("end dumping debug# %ld ... ",dump_cnt);

  return(0);

}




int debug_content(int i, int j, MPI_Datatype datatype,void *writebuf)
{
    // could also make everything FTYPE and convert like for normal i,j dump file
    myset(datatype,failfloorcount[i][j][0],0,NUMTSCALES*NUMFAILFLOORFLAGS,writebuf);
    
    return(0);
}


int avgdump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];
	


  trifprintf("begin dumping avgdump# %ld ... ",dump_cnt);

  whichdump=AVGCOL;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/avg");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");
  
  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,avg_content)>=1) return(1);

  trifprintf("end dumping avgdump# %ld ... ",dump_cnt);


  return(0);

}

int avg_content(int i, int j, MPI_Datatype datatype,void *writebuf)
{
  int k = 0, l = 0, col = 0;
  FTYPE r, th, vmin1, vmax1, vmin2, vmax2;
  struct of_geom geom;
  struct of_state q;
  FTYPE X[NDIM];
  FTYPE divb;
  FTYPE b[NDIM],ucon[NDIM];
  FTYPE U[NPR];
  FTYPE ftemp;


  coord(i, j, CENT, X);
  bl_coord(X, &r, &th);
  get_geometry(i, j, CENT, &geom);

  if(!GAMMIEDUMP){
    ftemp=(FTYPE)(i+startpos[1]);
    myset(datatype,&ftemp,0,1,writebuf);
    ftemp=(FTYPE)(j+startpos[2]);
    myset(datatype,&ftemp,0,1,writebuf);
  }
  myset(datatype,X,1,2,writebuf);
  myset(datatype,&r,0,1,writebuf);
  myset(datatype,&th,0,1,writebuf);

  myset(datatype,&geom.g,0,1,writebuf);

  // now do time average stuff
  myset(datatype,normalvarstavg[i][j],0,NUMNORMDUMP,writebuf);
  myset(datatype,anormalvarstavg[i][j],0,NUMNORMDUMP,writebuf);

  myset(datatype,jcontavg[i][j],0,NDIM,writebuf);
  myset(datatype,jcovtavg[i][j],0,NDIM,writebuf);
  myset(datatype,ajcontavg[i][j],0,NDIM,writebuf);
  myset(datatype,ajcovtavg[i][j],0,NDIM,writebuf);

  myset(datatype,massfluxtavg[i][j],0,NDIM,writebuf);
  myset(datatype,amassfluxtavg[i][j],0,NDIM,writebuf);

  myset(datatype,othertavg[i][j],0,NUMOTHER,writebuf);
  myset(datatype,aothertavg[i][j],0,NUMOTHER,writebuf);

  myset(datatype,fcontavg[i][j],0,NUMFARADAY,writebuf);
  myset(datatype,fcovtavg[i][j],0,NUMFARADAY,writebuf);
  myset(datatype,afcontavg[i][j],0,NUMFARADAY,writebuf);
  myset(datatype,afcovtavg[i][j],0,NUMFARADAY,writebuf);

#if(DOAVG2==0)
  myset(datatype,tudtavg[i][j],0,NUMSTRESSTERMS,writebuf);
  myset(datatype,atudtavg[i][j],0,NUMSTRESSTERMS,writebuf);
#endif

  return(0);

}


int avg2dump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];



  trifprintf("begin dumping avg2dump# %ld ... ",dump_cnt);

  whichdump=AVG2COL;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/avg2");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");

  if(dump_gen(WRITEFILE,dump_cnt,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,avg2_content)>=1) return(1);

  trifprintf("end dumping avg2dump# %ld ... ",dump_cnt);


  return(0);

}


int avg2_content(int i, int j, MPI_Datatype datatype,void *writebuf)
{
  int k = 0, l = 0, col = 0;
  FTYPE r, th, vmin1, vmax1, vmin2, vmax2;
  struct of_geom geom;
  struct of_state q;
  FTYPE X[NDIM];
  FTYPE divb;
  FTYPE b[NDIM],ucon[NDIM];
  FTYPE U[NPR];
  FTYPE ftemp;


  coord(i, j, CENT, X);
  bl_coord(X, &r, &th);
  get_geometry(i, j, CENT, &geom);
  // if you change # of outputted vars, remember to change numcolumns above

  if(!GAMMIEDUMP){
    ftemp=(FTYPE)(i+startpos[1]);
    myset(datatype,&ftemp,0,1,writebuf);
    ftemp=(FTYPE)(j+startpos[2]);
    myset(datatype,&ftemp,0,1,writebuf);
  }
  myset(datatype,X,1,2,writebuf);
  myset(datatype,&r,0,1,writebuf);
  myset(datatype,&th,0,1,writebuf);

  myset(datatype,&geom.g,0,1,writebuf);
  // 7

  myset(datatype,tudtavg[i][j],0,NUMSTRESSTERMS,writebuf);
  myset(datatype,atudtavg[i][j],0,NUMSTRESSTERMS,writebuf);
  // 112*2

  // total=7+112*2=231

    return(0);
}


int gdump(void)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];
	

  trifprintf("begin dumping gdump# %ld ... ",dump_cnt);

  whichdump=GDUMPCOL;
  datatype=MPI_FTYPE;
  strcpy(fileprefix,"dumps/gdump");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");

  if(dump_gen(WRITEFILE,-1,binaryoutput,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,gdump_content)>=1) return(1);

  trifprintf("end dumping gdump# %ld ... ",dump_cnt);

  return(0);
}




int gdump_content(int i, int j, MPI_Datatype datatype, void *writebuf)
{
  int k = 0, l = 0, m = 0, n = 0, col = 0;
  FTYPE r, th;
  FTYPE X[NDIM];
  FTYPE ftemp;
  FTYPE *ptrftemp;
  FTYPE dxdxp[NDIM][NDIM];

  coord(i, j, CENT, X);
  bl_coord(X, &r, &th);
  dxdxprim(X, r, th, dxdxp);


  ftemp=(FTYPE)(i+startpos[1]);
  myset(datatype,&ftemp,0,1,writebuf);
  ftemp=(FTYPE)(j+startpos[2]);
  myset(datatype,&ftemp,0,1,writebuf);
  // 2
  myset(datatype,X,1,2,writebuf);
  myset(datatype,&r,0,1,writebuf);
  myset(datatype,&th,0,1,writebuf);
  // 2+4
  ptrftemp=(FTYPE*)(&conn[i][j][0][0][0]);
  myset(datatype,ptrftemp,0,NDIM*NDIM*NDIM,writebuf);
  ptrftemp=(FTYPE*)(&gcon[i][j][0][0][0]);
  myset(datatype,ptrftemp,0,NPG*NDIM*NDIM,writebuf);
  ptrftemp=(FTYPE*)(&gcov[i][j][0][0][0]);
  myset(datatype,ptrftemp,0,NPG*NDIM*NDIM,writebuf);
  ptrftemp=(FTYPE*)(&gdet[i][j][0]);
  myset(datatype,ptrftemp,0,NPG,writebuf);
  ptrftemp=(FTYPE*)(&conn2[i][j][0]);
  myset(datatype,ptrftemp,0,NDIM,writebuf);
  // 4*4
  ptrftemp=(FTYPE*)(&dxdxp[0][0]);
  myset(datatype,ptrftemp,0,NDIM*NDIM,writebuf);


  return(0);

}



int fieldlinedump(long dump_cnt)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];
	


  trifprintf("begin dumping fieldlinedump# %ld ... ",dump_cnt);

  whichdump=FIELDLINECOL;
  datatype=MPI_FLOAT; // don't need good precision
  strcpy(fileprefix,"dumps/fieldline");
  strcpy(fileformat,"%04ld");
  strcpy(filesuffix,"");
  
  // MIXEDOUTPUT means text header and forced binary data
  if(dump_gen(WRITEFILE,dump_cnt,MIXEDOUTPUT,whichdump,datatype,fileprefix,fileformat,filesuffix,dump_header,fieldline_content)>=1) return(1);

  trifprintf("end dumping fieldlinedump# %ld ... ",dump_cnt);


  return(0);

}

int fieldline_content(int i, int j, MPI_Datatype datatype,void *writebuf)
{
  int k = 0, l = 0, col = 0;
  struct of_geom geom;
  struct of_state q;
  //FTYPE U[NPR];
  FTYPE FL[NPR];
  // must be same precision as written content
  float ftemp;

  //////////////
  //
  // some calculations
  //

  // if failed, then data output for below invalid, but columns still must exist    
  get_geometry(i, j, CENT, &geom);
  if (!failed) {
    if (get_state(p[i][j], &geom, &q) >= 1)
      FAILSTATEMENT("dump.c:dump()", "get_state() dir=0", 1);
  }
  else {// do a per zone check, otherwise set to 0
    whocalleducon=1; // force no failure mode, just return like failure, and don't return if failure, just set to 0 and continue
    if (get_state(p[i][j], &geom, &q) >= 1){
      for (k = 0; k < NDIM; k++)
	q.ucon[k]=0;
      for (k = 0; k < NDIM; k++)
	q.ucov[k]=0;
      for (k = 0; k < NDIM; k++)
	q.bcon[k]=0;
      for (k = 0; k < NDIM; k++)
	q.bcov[k]=0;
    }
    whocalleducon=0; // return to normal state
    
  }

  MYFUN(primtoflux(UDIAG,p[i][j], &q, RR, &geom, FL),"step_ch.c:fluxcalc()", "primtoflux_calc() dir=1/2 l", RR);


  //////////////////////////
  //
  // do the assignments
  //
  // if you change # of outputted vars, remember to change numcolumns



  ////////////////////
  //
  // 2 various things

  // rho (for various things)
  ftemp=p[i][j][RHO];
  myset(datatype,&ftemp,0,1,writebuf);

  // u (for various things)
  ftemp=p[i][j][UU];
  myset(datatype,&ftemp,0,1,writebuf);


  //////////////////////
  //
  // 2 things for jet/energy per baryon at infinity

  // -u_t (-hu_t can be found from this and rho/u/p above)
  ftemp=-q.ucov[0];
  myset(datatype,&ftemp,0,1,writebuf);

  // -T^t_t/(gdet rho u^t)
  //  ftemp=-U[UU]/(geom.g * p[i][j][RHO]*q.ucon[TT]);
  //myset(datatype,&ftemp,0,1,writebuf);

  // -T^r_t/(rho u^r)
  if(q.ucon[RR]!=0.0){
    ftemp=-FL[UU]/(geom.g * p[i][j][RHO]*q.ucon[RR]);
  }
  else ftemp=0.0;
  myset(datatype,&ftemp,0,1,writebuf);


  // 1 extra thing

  // u^t
  ftemp=q.ucon[0];
  myset(datatype,&ftemp,0,1,writebuf);


  ///////////////////////////
  //
  // 6 things for the field line stuff

  // v^r [ in grid frame]
  ftemp=q.ucon[1]/q.ucon[0];
  myset(datatype,&ftemp,0,1,writebuf);

  // v^\theta
  ftemp=q.ucon[2]/q.ucon[0];
  myset(datatype,&ftemp,0,1,writebuf);

  // v^\phi
  ftemp=q.ucon[3]/q.ucon[0];
  myset(datatype,&ftemp,0,1,writebuf);

  // B^r
  ftemp=p[i][j][B1];
  myset(datatype,&ftemp,0,1,writebuf);

  // B^\theta
  ftemp=p[i][j][B2];
  myset(datatype,&ftemp,0,1,writebuf);

  // B^\phi
  ftemp=p[i][j][B3];
  myset(datatype,&ftemp,0,1,writebuf);

  // see grmhd-dualfcon2omegaf.nb
  // below can be obtained from above set of v and B
  // \Omega_F_1
  //  ftemp=v3-B3*v2/(B2+SMALL);
  //myset(datatype,&ftemp,0,1,writebuf);

  // \Omega_F_2
  //ftemp=v3-B3*v1/(B1+SMALL);
  // myset(datatype,&ftemp,0,1,writebuf);

  return(0);

}













// I handle DOCOLSPLIT==1 simply with non-ROMIO MPI transfers.  With
// ROMIO, to obtain a small memory foot (i.e. smaller than just with
// DOCOLSPLIT==0), I loop over the dump process "numcolumns" times.
// This uses more CPU in favor of lower memory.  MINMEM version does
// this automatically.  original MPI version does DOCOLSPLIT, but not
// memory optimal.


int dump_gen(int readwrite, long dump_cnt, int bintxt, int whichdump, MPI_Datatype datatype, char *fileprefix, char *fileformat, char *filesuffix, int (*headerfun) (int bintxt, FILE*headerptr),int (*setgetcontent) (int i, int j, MPI_Datatype datatype, void*setbuf))
{
  int i = 0, j = 0, k = 0, l = 0, col = 0;
  FILE **fpp;
  char dfnam[MAXFILENAME];
  char dfnamreal[MAXFILENAME];
  char localfileformat[MAXFILENAME];
  void *jonio;
  void *writebuf;
  char truemyidtxt[MAXFILENAME];
  char filerw[MAXFILENAME];
  FILE *headerptr;
  int numfiles,coliter;
  void *setbuf;
  int sizeofdatatype;
  int romiocloopend;
  int headerbintxt;
  int dumpbintxt;
  fpos_t headerendpos;
  long headerbytesize;
  int binextension;


  ////////////
  //
  // setup file format for header and dump
  //
  ////////////

  if(readwrite==READFILE){
    if(bintxt==BINARYOUTPUT) strcpy(filerw,"r");
    else strcpy(filerw,"rt");
  }
  else if(readwrite==WRITEFILE){
    if(bintxt==BINARYOUTPUT) strcpy(filerw,"w");
    else strcpy(filerw,"wt");
  }

  if(bintxt==BINARYOUTPUT) headerbintxt=dumpbintxt=BINARYOUTPUT;
  else if(bintxt==TEXTOUTPUT) headerbintxt=dumpbintxt=TEXTOUTPUT;
  else if(bintxt==MIXEDOUTPUT){
    headerbintxt=TEXTOUTPUT;
    dumpbintxt=BINARYOUTPUT;
  }



  numcolumns=dnumcolumns[whichdump];
  docolsplit=DOCOLSPLIT[whichdump]; // docolsplit global var for now

  ////////////////////
  //
  // See if enough HD space
  //
  ////////////////////

  if(mpicombine){
    if(dumpbintxt==BINARYOUTPUT){
      if(myid==0) isenoughfreespace((unsigned long long)(sizeof(FTYPE))*(unsigned long long)numcolumns*(unsigned long long)(totalsize[1])*(unsigned long long)(totalsize[2]));
      else isenoughfreespace(0);
    }
    else{// text
      if(myid==0) isenoughfreespace((unsigned long long)(22)*(unsigned long long)numcolumns*(unsigned long long)(totalsize[1])*(unsigned long long)(totalsize[2]));
      else isenoughfreespace(0);
    }
  }
  else{
    if(dumpbintxt==BINARYOUTPUT){
      isenoughfreespace((unsigned long long)(sizeof(FTYPE))*(unsigned long long)numcolumns*(unsigned long long)(N1)*(unsigned long long)(N2));
    }
    else{// text
      isenoughfreespace((unsigned long long)(21)*(unsigned long long)numcolumns*(unsigned long long)(N1)*(unsigned long long)(N2));
    }
  }

  /////////////////////
  //
  // Allocate memory for setting up setbuf
  //
  //////////////////////

  sizeofdatatype=getsizeofdatatype(datatype);
  setbuf=malloc(numcolumns*sizeofdatatype);
  if(setbuf==NULL){
    dualfprintf(fail_file,"cannot allocate memory to setbuf in %s %s with numcolumns=%d and sizeofdatatype=%d\n",fileprefix,filesuffix,numcolumns,sizeofdatatype);
    myexit(1);
  }


  //////////////////////////////////
  //
  //  Set up DOCOLSPLIT for normal and ROMIO loop
  //
  ///////////////////////////////////



  if(docolsplit){
    numfiles=numcolumns;
    if(mpicombine&&USEMPI&&USEROMIO) romiocloopend=numfiles;
    else romiocloopend=1;
  }
  else{
    numfiles=1;
    romiocloopend=1;
  }


  
  //////////////////////////////////
  //
  //  Define file output and open it
  //
  ///////////////////////////////////

  // say whether .bin is allowed or not if binary
  if(fileprefix[0]=='i') binextension=0; // images don't need binary extension
  else binextension=1;


  // sometimes all CPUs need to know filename (e.g. ROMIO)
  // setup file suffix
  if((dumpbintxt==BINARYOUTPUT)&&(binextension)){
    if(USEMPI&&(mpicombine==0)&&(numprocs>1)) sprintf(truemyidtxt,".bin.%04d",myid);
    else strcpy(truemyidtxt,".bin");
  }
  else{
    if(USEMPI&&(mpicombine==0)&&(numprocs>1)) sprintf(truemyidtxt,".%04d",myid);
    else strcpy(truemyidtxt,"");
  }

  // setup filename
  if(dump_cnt>=0){
    strcpy(localfileformat,"%s");
    strcat(localfileformat,fileformat);
    strcat(localfileformat,"%s");
    strcat(localfileformat,"%s");
    sprintf(dfnam, localfileformat, fileprefix, dump_cnt, filesuffix, truemyidtxt);
  }
  else{ // then no file number wanted (i.e. for gdump())
    sprintf(dfnam, "%s%s%s", fileprefix, filesuffix, truemyidtxt);
  }



  ////////////////
  //
  // open files, or open files for header if mpicombine==1, which for mpicombine==1 gets reopened later by MPI routines
  //
  ///////////////

  if((USEMPI&&(myid==0)&&(mpicombine==1))||(mpicombine==0)){// for mpicombine==1 even with ROMIO, real filename and header not needed

    // only one CPU does header if mpicombine==1, header+dump done in all CPUs if mpicombine==0
    // create files for each column, or each column's header if mpicombine==1
    if((fpp=(FILE**)malloc(sizeof(FILE*)*numfiles))==NULL){
      dualfprintf(fail_file,"couldn't open fpp in dump()\n");
      myexit(2);
    }// now fpp[i] indexes a list of file pointers


    // setup each file corresponding to each column
    COLLOOP{
      if(docolsplit&&(numfiles>1)){
	sprintf(dfnamreal,"%s-col%04d",dfnam,coliter);
      }
      else strcpy(dfnamreal,dfnam);
      
      if ((fpp[coliter] = fopen(dfnamreal, filerw)) == NULL) {
	dualfprintf(fail_file, "error opening %s : %s %s file\n",dfnamreal,fileprefix,filesuffix);
	myexit(2);
      }
      //////////////////////////////////
      //
      //  read or write header: header is read/written in whatever style chosen to the top of each dump file created
      //
      ///////////////////////////////////
      headerfun(headerbintxt,fpp[coliter]); // outputs header to each column file (or just one file, or all CPU files, etc.)
      headerbytesize=ftell(fpp[coliter]);
    }
    
    // don't close if mpicombine==0, since used in a moment, else mpicombine==1 it's reopened by MPI routines
    if(USEMPI&&(myid==0)&&(mpicombine==1)) COLLOOP fclose(fpp[coliter]); // will get reopened later by MPI routines


  }
  // need to broadcast the header size to other CPUs for ROMIO
#if(USEMPI&&USEROMIO)
  MPI_Bcast(&headerbytesize,1,MPI_LONG,0,MPI_COMM_WORLD);
#endif

    
  ///////////////////////////////////////////////////////////
  //
  // loop over columns for per-column buffer ROMIO dump
  //
  //
  ////////////////////////////////////////////////////////////

  ROMIOCOLLOOP{ // only loop if mpicombine&&USEMPI&&USEROMIO&&docolsplit==1
    if(romiocloopend>1) trifprintf("romiocoliter=%d of romiocloopend=%d\n",romiocoliter,romiocloopend);


    
    // setup MPI buffer if mpicombine==1
    if( mpicombine == 0 ) { // then one file per CPU if USEMPI or just normal file writing on 1CPU
      writebuf=NULL;
    }
    else mpiio_init(dumpbintxt,sortedoutput, fpp, headerbytesize, readwrite, dfnam, numcolumns, datatype, &jonio, &writebuf);
    // if USEROMIO==1 then numcolumns interpreted properly for docolsplit


    if(readwrite==READFILE){
      //////////////////////////////////
      //
      // read DUMP 
      //
      //////////////////////////////////
      
      if (mpicombine == 1) {
#if(USEMPI)
	mpiio_seperate(binaryoutput,sortedoutput, STAGE1, numcolumns, datatype, fpp, jonio, writebuf);
#endif
      }
    }




    //////////////////
    //
    // DUMP LOOP
    //
    //////////////////



    if(readwrite==READFILE){
      BUFFERINIT0;
      //      DUMPLOOP(-N1BND, N1 - 1+N1BND, -N2BND, N2 - 1+N2BND) {
      DUMPLOOP(0, N1 - 1, 0, N2 - 1) {
	// buffer init starts the parallel index
	BUFFERINIT;
	// initialize to 0th column
	COLINIT;
	///////////////////////////////////////
	//
	// READFILE
	//
	//////////////////////
	if((mpicombine)&&(truempicombinetype==MPICOMBINEMINMEM)) mpiio_minmem(READFILE,whichdump,i,j,dumpbintxt,sortedoutput,numcolumns,datatype, fpp,jonio,writebuf);
	
	// read all at once
	myfread(dumpbintxt,datatype,setbuf,0,numcolumns,i,j,fpp,writebuf);
	
	// check
	if(nextbuf!=numcolumns){
	  dualfprintf(fail_file,"Number of columns (numcolumns=%d) isn't equal to number of columns/buffers attempted (nextbuf=%d)\n",numcolumns,nextbuf);
	  myexit(1);
	}
	
	// get the content of 1 row
	setgetcontent(i,j,datatype,setbuf);
	
	// check
	if(nextcol!=numcolumns){
	  dualfprintf(fail_file,"Number of columns (numcolumns=%d) isn't equal to number of columns attempted (nextcol=%d)\n",numcolumns,nextcol);
	  myexit(1);
	}
      }// end DUMPLOOP
    }// end readwrite==READFILE
    else if(readwrite==WRITEFILE){
      BUFFERINIT0;
      //      DUMPLOOP(-N1BND, N1-1+N1BND, -N2BND, N2 - 1 + N2BND) {
      DUMPLOOP(0, N1-1, 0, N2 - 1) {
	// buffer init starts the parallel index
	BUFFERINIT;
	// initialize to 0th column
	COLINIT;
	///////////////////////////////////////
	//
	// WRITEFILE
	//
	//////////////////////

	// set the content of 1 row
	setgetcontent(i,j,datatype,setbuf);
	
	// check
	if(nextcol!=numcolumns){
	  dualfprintf(fail_file,"Number of columns (numcolumns=%d) isn't equal to number of columns attempted (nextcol=%d)\n",numcolumns,nextcol);
	  myexit(1);
	}
	
	// write all at once
	myfwrite(dumpbintxt,datatype,setbuf,0,numcolumns,i,j,fpp,writebuf);
	
	// check
	if(nextbuf!=numcolumns){
	  dualfprintf(fail_file,"Number of columns (numcolumns=%d) isn't equal to number of columns/buffers attempted (nextbuf=%d)\n",numcolumns,nextbuf);
	  myexit(1);
	}
	
	// finish up this row
	if((mpicombine==0)&&(dumpbintxt==TEXTOUTPUT)) COLLOOP fprintf(fpp[coliter],"\n");
	if((mpicombine)&&(truempicombinetype==MPICOMBINEMINMEM)) mpiio_minmem(WRITEFILE,whichdump,i,j,dumpbintxt,sortedoutput,numcolumns,datatype, fpp, jonio,writebuf);
      }// end DUMPLOOP
    }//end readwrite==WRITEFILE
  



    //////////////////
    //
    // Close dump file and write/close file if mpicombine==1
    //
    //////////////////

    if (mpicombine == 0){
      COLLOOP if (fpp[coliter] != NULL) fclose(fpp[coliter]);
    }
    else{
#if(USEMPI)
      if(readwrite==WRITEFILE) mpiio_combine(dumpbintxt, sortedoutput, numcolumns, datatype, fpp, jonio, writebuf);
      else if(readwrite==READFILE) mpiio_seperate(binaryoutput,sortedoutput, STAGE2, numcolumns, datatype, fpp, jonio, writebuf);
#endif
    }

  }// end column loop for ROMIO&&docolsplit

  // free the set/get buffer
  if(setbuf!=NULL) free(setbuf);

  return (0);
}


