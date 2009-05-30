/* 
   produces an "r8" file. */

#include "decs.h"

#define GAMMIESTARTI (0)
#define GAMMIEENDI (0)

#define JONSTARTI (0)
#define JONENDI (NPRDUMP-1)

#define NORMAL (0)
#define ZOOM (1)
// revert to only normal for now
//#define NUMLIMITS (2)
#define NUMLIMITS (1)

#define LOG (0)
#define LINEAR (1)
#define NUMSCALE (2)

#define DOCONS (1)
// whether to do (0=no, 1=yes) conserved quantities in image

#define MINVECTOR (1E-10)



int image_dump(long dump_cnt)
{
  int startk,endk,starts,ends,startl,endl,startv,endv;
  int whichk,limits,scale,vartype;
  int i,j, floor;
  ////////////////////////////
  //
  // Image Loop
  //
  ////////////////////////////

  // STARTI=0 is normal
  // ENDI=NPRDUMP is normal, but can be NPRDUMP+2 to get linear RHO/UU
#if(PRODUCTION==0)
  if(GAMMIEIMAGE){
    startk=GAMMIESTARTI;
    endk=GAMMIEENDI;
    starts=0;
    ends=1;
    startl=0;
    endl=0;
    startv=0;
    endv=0;
  }
  else{
    startk=JONSTARTI;
    endk=JONENDI;
    starts=0;
    ends=NUMSCALE-1;
    startl=0;
    endl=NUMLIMITS-1;
    startv=0;
    //    endv=1;
    endv=2; // including failures now
  }
#else
  // no special images for production mode, just basic log density
  if(GAMMIEIMAGE){
#if(EOMTYPE==EOMGRMHD)
    startk=0;
    endk=0;
#elif(EOMTYPE==EOMFFDE)
    startk=B1;
    endk=B1;
#endif
    starts=0;
    ends=0;
    startl=0;
    endl=0;
    startv=0;
    endv=0;
  }
  else{
#if(EOMTYPE==EOMGRMHD)
    startk=0;
    endk=0;
#elif(EOMTYPE==EOMFFDE)
    startk=B1;
    endk=B1;
#endif
    starts=0;
    ends=0;
    startl=0;
    endl=0;
    startv=0;
    endv=0;
  }
#endif

  for(vartype=startv;vartype<=endv;vartype++){
    for(limits=startl;limits<=endl;limits++){
      for(scale=starts;scale<=ends;scale++){
	for (whichk = startk; whichk <= endk; whichk++) {
	  if((vartype<=1)||((vartype==2)&&(DODEBUG)&&(whichk<=NUMFAILFLOORFLAGS-1))){
	    if(image(dump_cnt,whichk,scale,limits,vartype)>=1) return(1);
	  }
	}
      }
    }
  }

  return(0);
}


int imagedefs(int whichk, int scale, int limits, int vartype)
{
  int i = 0, j = 0, l = 0, col = 0, floor;
  FILE *fp;
  // whichk : whichk primitive variable
  FTYPE pr,iq, liq, lmax;
  unsigned char liqb;
  FTYPE min,max,sum;
  FTYPE minptr[NPR], maxptr[NPR], sumptr[NPR];
  char truemyidtxt[MAXFILENAME];
  FTYPE U[NPR];
  struct of_geom geom;
  struct of_state q;
  FTYPE X[NDIM],r,th;
  FTYPE lmin,aa;


  ////////////////////////////
  //
  // Image output setup/definition
  //
  // Purpose is to set pimage to correct variable type (primitive or conservative), limits, scale, and which k.
  // Then image() outputs that one thing to file
  //
  ////////////////////////////

  pimage=pk[0]; // assume pk[0] isn't needed here since should be done with timestep
  if(limits==ZOOM){ // zoom in on dynamic range of values to see fine details
    ZLOOP{
      if(vartype==0){
	if(whichk<=1){
	  coord(i,j,CENT,X);
	  bl_coord(X,&r,&th);
	  if(whichk==0) pimage[i][j][whichk]=p[i][j][whichk]/(RHOMIN*pow(r,-1.5));
	  if(whichk==1) pimage[i][j][whichk]=p[i][j][whichk]/(UUMIN*pow(r,-2.5));
	}
	else{
	  if(scale==LINEAR) pimage[i][j][whichk]=p[i][j][whichk];
	  else if(scale==LOG) pimage[i][j][whichk]=fabs(p[i][j][whichk])+MINVECTOR;
	}
      }
      else if(vartype==1){// conserved quantity
	// computes too much (all conserved quantites every time)
	get_geometry(i,j,CENT,&geom) ;
	if(!failed){
	  if(get_state(p[i][j],&geom,&q)>=1) return(1);
	  if(primtoU(UDIAG,p[i][j],&q,&geom,U)>=1) return(1);
	}
	if(scale==LINEAR) pimage[i][j][whichk]=U[whichk];
	else if(scale==LOG) pimage[i][j][whichk]=fabs(U[whichk]/geom.g)+MINVECTOR;
      }
      else if(vartype==2){ // failure quantity (no diff from below right now -- could zoom in on single failure regions)
	if(whichk<NUMFAILFLOORFLAGS){
	  floor=whichk;
	  if(scale==LINEAR) pimage[i][j][whichk]=(FTYPE)failfloorcount[i][j][IMAGETS][floor];
	  else if(scale==LOG) pimage[i][j][whichk]=fabs((FTYPE)failfloorcount[i][j][IMAGETS][floor]+1);
	}
      }
    }
  }
  else{
    ZLOOP{
      if(vartype==0){
	if(whichk<=1) pimage[i][j][whichk]=p[i][j][whichk];
	else{
	  if(scale==LINEAR) pimage[i][j][whichk]=p[i][j][whichk];
	  else if(scale==LOG) pimage[i][j][whichk]=fabs(p[i][j][whichk])+MINVECTOR;
	}
      }
      else if(vartype==1){// conserved quantity
	// computes too much (all conserved quantites every time)
	get_geometry(i,j,CENT,&geom) ;
	if(!failed){
	  if(get_state(p[i][j],&geom,&q)>=1) return(1);
	  if(primtoU(UDIAG,p[i][j],&q,&geom,U)>=1) return(1);
	}
	if(scale==LINEAR) pimage[i][j][whichk]=U[whichk];
	else if(scale==LOG) pimage[i][j][whichk]=fabs(U[whichk]/geom.g)+MINVECTOR;
      }
      else if(vartype==2){ // failure quantity
	if(whichk<NUMFAILFLOORFLAGS){
	  floor=whichk;
	  if(scale==LINEAR) pimage[i][j][whichk]=(FTYPE)failfloorcount[i][j][IMAGETS][floor];
	  else if(scale==LOG) pimage[i][j][whichk]=fabs((FTYPE)failfloorcount[i][j][IMAGETS][floor]+1);
	}
      }
    }
  }



  ////////////////////////////
  //
  // Image FILE open/initialize
  //
  ////////////////////////////



  ////////////////////
  //
  // Image paramters setup (whole purpose currently is to find lmin and aa)
  // 
  /////////////////////

  /* density mapping is logarithmic, in 255 steps between e^lmax and
     e^lmin */

#define ZOOMFACTOR (10000)

  prminmaxsum(pimage,whichk,1,maxptr,minptr,sumptr);
  if(limits==NORMAL){
    max=maxptr[whichk];
    min=minptr[whichk];
    sum=sumptr[whichk];
  }
  else{
    if(whichk<=1){
      max=maxptr[whichk]/ZOOMFACTOR;
      min=minptr[whichk];
    }
    else{
      if(scale==LINEAR){
	max=maxptr[whichk]/ZOOMFACTOR;
	min=minptr[whichk]/ZOOMFACTOR;
      }
      else{
	max=maxptr[whichk]/ZOOMFACTOR;
	min=minptr[whichk];
      }
    }
  }
  sum=sumptr[whichk];
  logsfprintf("whichk: %d scale: %d limits: %d : min,max,avg: %21.15g %21.15g %21.15g\n",whichk,scale,limits,min,max,sum/totalzones);

  if(scale==LOG){
    lmax = log(max);
    lmin = log(min);
  } else if(scale==LINEAR) {
    lmax = max;
    lmin = min;
  }
  else{
    dualfprintf(fail_file,"no such scale=%d\n",scale);
    myexit(1);
  }

  if (lmax != lmin)
    aa = 256. / (lmax - lmin);
  else
    aa = 0;



  // set the only paramters needed to dump image
  imageparms[ORIGIN]=aa;
  imageparms[LMIN]=lmin;

  return(0);
}




int image(long dump_cnt, int whichk, int scale, int limits, int vartype)
{
  MPI_Datatype datatype;
  int whichdump;
  char fileprefix[MAXFILENAME];
  char filesuffix[MAXFILENAME];
  char fileformat[MAXFILENAME];

  trifprintf("begin dumping image# %ld whichk: %d scale: %d limits: %d vartype: %d\n ",dump_cnt,whichk,scale,limits,vartype);

  // global vars so can avoid passing through general functions
  imagescale=scale;
  imagevartype=vartype;
  imagelimits=limits;
  imagewhichk=whichk;



  whichdump=IMAGECOL;
  datatype=MPI_UNSIGNED_CHAR;

  // actual prefix
  if(GAMMIEIMAGE&&(GAMMIESTARTI==GAMMIEENDI)&&(GAMMIESTARTI==0)){
    sprintf(fileprefix, "images/im");
  }
  else{
    if(vartype==0) sprintf(fileprefix, "images/im%1dp%1ds%1dl", whichk, scale, limits);
    else if(vartype==1)  sprintf(fileprefix, "images/im%1dc%1ds%1dl", whichk, scale, limits);
    else if(vartype==2)  sprintf(fileprefix, "images/im%1df%1ds%1dl", whichk, scale, limits);
  }
  strcpy(fileformat,"%04ld"); // actual format
  strcpy(filesuffix,".r8"); // actual suffix

  // setup the image definitions (min,max, what data, etc.)
  if(imagedefs(whichk,scale,limits,vartype)>=1) return(1);


  // MIXEDOUTPUT tells dump_gen() to treat as .r8 file, with text header and binary dump
  if(dump_gen(WRITEFILE,dump_cnt,MIXEDOUTPUT,whichdump,datatype,fileprefix,fileformat,filesuffix,image_header,image_content)>=1) return(1);

  trifprintf("end dumping image# %ld whichk: %d scale: %d limits: %d vartype: %d\n ",dump_cnt,whichk,scale,limits,vartype);

  return(0);

}


int image_header(int bintxt, FILE *headerptr)
{
  ////////////////////////////
  //
  // HEADER file open/initialize
  //
  ////////////////////////////

  // write header
  if(bintxt==TEXTOUTPUT){// should always be TEXTOUTPUT
    fprintf(headerptr, "RAW\n# t=%21.15g vartype=%2d  whichk=%2d scale=%2d limits=%2d \n%i %i\n255\n", t, imagevartype, imagewhichk, imagescale, imagelimits, totalsize[1],totalsize[2]);
    //fprintf(headerptr, "RAW\n# t=%21.15g vartype=%2d  whichk=%2d scale=%2d limits=%2d \n%i %i\n255\n", t, imagevartype, imagewhichk, imagescale, imagelimits, totalsize[1]+2*N1BND,totalsize[2]+2*N2BND);
  }
  else{
    dualfprintf(fail_file,"Shouldn't be trying to write binary header to image file\n");
    myexit(1);
  }
  fflush(headerptr);

  return(0);
}


// uses global vars to get aa and lmin
int image_content(int i, int j, MPI_Datatype datatype,void *writebuf)
{
  unsigned char liqb;
  FTYPE pr,iq,liq;
  FTYPE aa;
  FTYPE lmin;

  aa=imageparms[ORIGIN];
  lmin=imageparms[LMIN];

  pr=pimage[i][j][imagewhichk];
  if (imagescale==LOG) iq = log(pr);
  else if(imagescale==LINEAR) iq = pr;
  
  liq = aa * (iq - lmin);
  if (liq > 255.)
    liq = 255.;
  if (liq < 0.)
    liq = 0.;
  
  liqb=(unsigned char)(liq);
  
  myset(datatype,&liqb,0,1,writebuf);

  return(0);
}
