#include "decs.h"

int dump_ener(int doener, int dordump, int call_code)
{

  static FILE *flenerreg_file[NUMENERREGIONS];
  static FILE *flener_file;
  char FLENERREGIONNAME[NUMENERREGIONS][MAXFILENAME];
  char *FLENERNAME;

  static FILE *enerreg_file[NUMENERREGIONS];
  static FILE *ener_file;
  char ENERREGIONNAME[NUMENERREGIONS][MAXFILENAME];
  char *ENERNAME;

  static FILE *generreg_file[NUMENERREGIONS];
  static FILE *gener_file;
  char GENERREGIONNAME[NUMENERREGIONS][MAXFILENAME];
  char *GENERNAME;

  // per cpu file
  static FILE *probe_file;

  // single region files
  static FILE *debug_file;
  static FILE *lumener_file;

  char dfnam[MAXFILENAME], ifnam[MAXFILENAME];
  int i, j, k, l, dir,sc,fl,floor,tscale;
  //  FILE *imagecnt_file, *dumpcnt_file,*avgcnt_file,*debugcnt_file,*lumvsrcnt_file;

  // full grid only
  SFTYPE divb, divbmax = 0, divbavg =0;
  // special jet2
  SFTYPE pdottermsjet2_tot[COMPDIM*2][NUMFLUXTERMS][NPR];

  // all regions, local quantities
  SFTYPE Ureg_tot[NUMENERREGIONS][NPR];
  SFTYPE Ureg_final[NUMENERREGIONS][NPR];
  SFTYPE pdotreg_tot[NUMENERREGIONS][COMPDIM*2][NPR];
  SFTYPE pdottermsreg_tot[NUMENERREGIONS][COMPDIM*2][NUMFLUXTERMS][NPR];
  // each region, local quantities
  SFTYPE *U_tot;
  SFTYPE *U_final;
  SFTYPE (*pdot_tot)[NPR];
  SFTYPE (*pdotterms_tot)[NUMFLUXTERMS][NPR];

  int enerregion;

  static int firsttime = 1;
  SFTYPE ftemp0,ftemp1;

  int ii;


  //////////////////////////////////
  //
  // Open/append some files
  //

  if ((call_code == INIT_OUT) || (firsttime == 1)) {
    // PER CPU FILE
    sprintf(dfnam,"probe.dat%s",myidtxt);
    if((probe_file=fopen(dfnam,"at"))==NULL){
      dualfprintf(fail_file,"Can't open probe file\n");
      myexit(1);
    }
    
    
    if(myid==0) if(DODEBUG){
      // CPU=0 only
      sprintf(dfnam,"debug.out");
      myfopen(dfnam,"a+","error opening debug output file\n",&debug_file);
    }

    if(myid==0) if(DOLUMVSR){
      // CPU=0 only
      sprintf(dfnam,"lumvsr.out");
      myfopen(dfnam,"a+","error opening lumvsr output file\n",&lumener_file);
    }
      
    // CPU=0 only
    if(myid==0) ENERREGIONLOOP{
      // set some pointers
      FLENERNAME=FLENERREGIONNAME[enerregion];

      // when setting file pointers, need pointer to pointer or just set pointers directly (we do latter for naming reasons)
      /* 	flener_file=flenerreg_file[enerregion]; */
      /* 	ener_file=&enerreg_file[enerregion]; */
      /* 	gener_file=generreg_file[enerregion]; */
      ENERNAME=ENERREGIONNAME[enerregion];

      GENERNAME=GENERREGIONNAME[enerregion];

      pcum_tot=pcumreg_tot[enerregion];
      fladd_tot=fladdreg_tot[enerregion];
      sourceadd_tot=sourceaddreg_tot[enerregion];

      // flener
      if(enerregion==0) sprintf(FLENERNAME,"flener.out");
      else sprintf(FLENERNAME,"flenerjet%d.out",enerregion-1); // specific naming convention
      myfopen(FLENERNAME,"a+","error opening FLenergy output file\n",&flenerreg_file[enerregion]);

      // ener
      if(enerregion==0) sprintf(ENERNAME,ENERFNAME);
      else sprintf(ENERNAME,"enerjet%d.out",enerregion-1); // specific naming convention

      if (appendold == 0) {
	// a+ in general, or w for overwrites, but let's append for now
	myfopen(ENERNAME,"a+","error opening  energy output file\n",&enerreg_file[enerregion]);
      }
      else {			// if appendold==1
	trifprintf("Start setup of %s file append\n", ENERNAME);
	myfopen(ENERNAME,"a+","error opening  energy output file for append\n",&enerreg_file[enerregion]);
	appendener(ener_file,pcum_tot,fladd_tot,sourceadd_tot);
      }

      // gener
      if(enerregion==0) sprintf(GENERNAME,GENERFNAME);
      else sprintf(GENERNAME,"generjet%d.out",enerregion-1); // specific naming convention
      myfopen(GENERNAME,"a+","error opening  energy output file\n",&generreg_file[enerregion]);
    }
  }


  /////////////////////
  //
  // compute divB and output to logs
  //
  divbmaxavg(p,&divbmax,&divbavg);



  /////////////////////////////////////////
  //
  // ENER/FLENER/FRDOT/DEBUG output integrated diagnostics
  // for any number of predetermined regions
  //
  ////////////////////////////////////////

  ENERREGIONLOOP{
    trifprintf(".BE.%d",enerregion);

    //////////////////////////////
    //
    // setup some pointers
    //
    //
    // each region, local quantities
    U_tot=Ureg_tot[enerregion];
    U_final=Ureg_final[enerregion];
    pdot_tot=pdotreg_tot[enerregion];
    pdotterms_tot=pdottermsreg_tot[enerregion];
    // used for each region, related to global quantities
    fladd=fladdreg[enerregion];
    fladd_tot=fladdreg_tot[enerregion];
    fladdterms=fladdtermsreg[enerregion];
    fladdterms_tot=fladdtermsreg_tot[enerregion];
    U_init=Ureg_init[enerregion];
    U_init_tot=Ureg_init_tot[enerregion];
    pcum=pcumreg[enerregion];
    pcum_tot=pcumreg_tot[enerregion];
    pdot=pdotreg[enerregion];
    pdotterms=pdottermsreg[enerregion];
    sourceaddterms=sourceaddtermsreg[enerregion];
    sourceaddterms_tot=sourceaddtermsreg_tot[enerregion];
    sourceadd=sourceaddreg[enerregion];
    sourceadd_tot=sourceaddreg_tot[enerregion];

    if(myid==0){ // only cpu=0 writes to these files, although doesn't matter
      ener_file=enerreg_file[enerregion];
      gener_file=generreg_file[enerregion];
      flener_file=flenerreg_file[enerregion];
    }

    /////////////////////////////////
    //
    // do some integrations
    //
    if(dordump||doener){

      trifprintf("BI%d",enerregion);

      // compute total conserved quantity
      if(integrate(U_tot,U_tot,CONSTYPE,enerregion)>=1) return(1);

      DIRLOOP if(integrate(pdot[dir],pdot_tot[dir],SURFACETYPE,enerregion)>=1) return(1);
      // above should be sum of below, approximately
      DIRLOOP FLLOOP if(integrate(pdotterms[dir][fl],pdotterms_tot[dir][fl],SURFACETYPE,enerregion)>=1) return(1);
      DIRLOOP if(integrate(pcum[dir],pcum_tot[dir],SURFACETYPE,enerregion)>=1) return(1);
      
      FLOORLOOP if(integrate(fladdterms[floor],fladdterms_tot[floor],CUMULATIVETYPE,enerregion)>=1) return(1);
      if(integrate(fladd,fladd_tot,CUMULATIVETYPE,enerregion)>=1) return(1);

      SCLOOP if(integrate(sourceaddterms[sc],sourceaddterms_tot[sc],CUMULATIVETYPE,enerregion)>=1) return(1);
      if(integrate(sourceadd,sourceadd_tot,CUMULATIVETYPE,enerregion)>=1) return(1);

      if(enerregion==0){
	// CUMULATIVETYPE2 for 1 data object per ii
	if(DOLUMVSR) for(ii=0;ii<ncpux1*N1;ii++) if(integrate(&lumvsr[ii],&lumvsr_tot[ii],CUMULATIVETYPE2,enerregion)>=1) return(1);

	// special jet2 accounting
	DIRLOOP FLLOOP if(integrate(pdottermsjet2[dir][fl],pdottermsjet2_tot[dir][fl],SURFACETYPE,enerregion)>=1) return(1);
	
	// debug
	if(DODEBUG){
	  TSCALELOOP if(integratel(failfloorcountlocal[tscale],failfloorcountlocal_tot[tscale],CONSTYPE,tscale, enerregion)>=1) return(1);
	}
      }
      trifprintf("AI%d",enerregion);
    }

    /////////////////////////////
    //
    // initialize the total conserved counters
    //
    if (call_code == INIT_OUT || firsttime) {
      if(RESTARTMODE==0)    PLOOP U_init_tot[k] = U_init[k] = U_tot[k];
      // otherwise read from restart file
    }

    /////////////////////////////
    //
    // output some interesting diagnostics to stderr/log files.
    //
    if (call_code == FINAL_OUT) {
      PLOOP U_final[k] = U_tot[k];

      if(GAMMIEENER){
	trifprintf("\n\nEnergy: ini,fin,del: %21.15g %21.15g %21.15g\n",
		   U_init[UU], U_final[UU], (U_final[UU] - U_init[UU]) / U_init[UU]);

	trifprintf("\n\nMass: ini,fin,del: %21.15g %21.15g %21.15g\n",
		   U_init[RHO], U_final[RHO], (U_final[RHO] - U_init[RHO]) / U_init[RHO]);
      }
      else{
	trifprintf("\n");
	PDUMPLOOP{
	  ftemp0=(U_final[k] - U_init[k]);
	  ftemp1=(U_final[k]-fladd_tot[k]-sourceadd_tot[k]-(pcum_tot[X1DN][k]+pcum_tot[X2DN][k])+(pcum_tot[X1UP][k]+pcum_tot[X2UP][k]) - U_init[k]);
	  trifprintf("U[%d]: ini,fin,fdf,del,tfdf,tdel: %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",
		     k,
		     U_init[k],
		     U_final[k],
		     U_init[k]!=0.0 ? ftemp0/U_init[k] : 0.0,
		     ftemp0,
		     U_init[k]!=0.0 ? ftemp1/U_init[k] : 0.0,
		     ftemp1 );
	}	
      }
    }


    /////////////////////
    //
    //  output some ener diagnostics to ener, flener, debug, and probe(only cpu specific one) files
    //
    if (doener){
      if(1||GAMMIEENER){
	// 6+3=9 terms
	myfprintf(gener_file, "%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g ", t, U_tot[RHO], U_tot[U3], U_tot[UU], p[N1 / 2][N2 / 2][UU] * pow(p[N1 / 2][N2 / 2][RHO], -gam), p[N1 / 2][N2 / 2][UU]);
	myfprintf(gener_file, "%21.15g %21.15g %21.15g ", pdot_tot[X1DN][RHO], pdot_tot[X1DN][RHO]-pdot_tot[X1DN][UU], pdot_tot[X1DN][U3]);
      }
      if(1){
	// SM use gammie.m macro, jrdpener, gammieener's
	//////////////////////////
	// ENER FILE (only dir and k)
	
	// 1+NPR+COMPDIM*2*NPRDUMP+NPRDUMP+(2) terms
	myfprintf(ener_file,"%21.15g %ld ",t ,realnstep);
	PDUMPLOOP myfprintf(ener_file, "%21.15g ", U_tot[k]);
	DIRLOOP PDUMPLOOP myfprintf(ener_file, "%21.15g ", pdot_tot[dir][k]);
	// COMPDIM*2*NPRDUMP more terms
	DIRLOOP PDUMPLOOP myfprintf(ener_file, "%21.15g ", pcum_tot[dir][k]);
	PDUMPLOOP myfprintf(ener_file, "%21.15g ", fladd_tot[k]);
	PDUMPLOOP myfprintf(ener_file, "%21.15g ", sourceadd_tot[k]);
	// no need for divb separately in other regions
	if(enerregion==0) myfprintf(ener_file, "%21.15g %21.15g ", divbmax,divbavg);
      }
      ////////////////////////////////
      // FLENER FILE (dir, k, and linear summed terms for flux, floor, and source quantities)
      // (note the change in ordering for external file read in macros/functions)
      myfprintf(flener_file,"%21.15g %ld ",t,realnstep);
      DIRLOOP PDUMPLOOP FLLOOP myfprintf(flener_file, "%21.15g ", pdotterms_tot[dir][fl][k]);
      PDUMPLOOP FLOORLOOP myfprintf(flener_file, "%21.15g ", fladdterms_tot[floor][k]);
      PDUMPLOOP SCLOOP myfprintf(flener_file, "%21.15g ", sourceaddterms_tot[sc][k]);
      
      if(enerregion==0){ // only for total region for now
	// only care about inner and outer radial part
	for(dir=0;dir<=1;dir++) PDUMPLOOP FLLOOP myfprintf(flener_file, "%21.15g ", pdottermsjet2_tot[dir][fl][k]); // jet/pole values only

	if(DOLUMVSR){
	  // luminosity vs radius
	  myfprintf(lumener_file,"%21.15g %ld ",t,realnstep);
	  for(ii=0;ii<ncpux1*N1;ii++) myfprintf(lumener_file, "%21.15g ", lumvsr_tot[ii]);
	}
	
	/////////////////////////
	// PROBE FILE
	// !!per CPU!! probe file NPRDUMP*3 terms
#define ITER1 MAX(N1/4,1)
#define ITER2 MAX(N2/4,1)
	for(i=0;i<N1;i+=ITER1) for(j=0;j<N2;j+=ITER2){
	  PDUMPLOOP fprintf(probe_file, "%d %d %d %ld %21.15g %21.15g\n",startpos[1]+i,startpos[2]+j,k,realnstep, t, p[i][j][k]);
	}
	fflush(probe_file);
	if(DODEBUG){
	  /////////////////////////
	  // DEBUG FILE
	  myfprintf(debug_file,"%21.15g %ld ",t,realnstep);
#if(COUNTTYPE==LONGLONGINTTYPE)
	    TSCALELOOP FLOORLOOP myfprintf(debug_file, "%lld ", failfloorcountlocal_tot[tscale][floor]);
#elif(COUNTTYPE==DOUBLETYPE)
	    TSCALELOOP FLOORLOOP myfprintf(debug_file, "%21.15g ", failfloorcountlocal_tot[tscale][floor]);
#endif
	}
      }

      myfprintf(flener_file,"\n");
      myfprintf(ener_file,"\n");
      myfprintf(gener_file,"\n");
      if((DOLUMVSR)&&(enerregion==0)) myfprintf(lumener_file,"\n");
      if((DODEBUG)&&(enerregion==0)) myfprintf(debug_file,"\n");

      trifprintf("W%d",enerregion);

    }// end if doener // done if writing to file

    ///////////////
    //
    // close the files at end of simulation
    //
    if(call_code==FINAL_OUT){
      myfclose(&flener_file,"Couldn't close flener_file\n");
      if(1) myfclose(&ener_file,"Couldn't close ener_file\n");
      if(1||GAMMIEENER) myfclose(&gener_file,"Couldn't close gener_file\n");
      if(enerregion==0){
	if(DOLUMVSR) myfclose(&lumener_file,"Couldn't close lumener_file\n");
	if(DODEBUG) myfclose(&debug_file,"Couldn't close debug_file\n");
	fclose(probe_file);
      }
    }
  }// end enerregion loop



  trifprintf("E");
  firsttime = 0;
  return (0);




}



void appendener(FILE* ener_file,SFTYPE pcum_tot[][NPR],SFTYPE*fladd_tot,SFTYPE *sourceadd_tot)
{
  int gotit;
  SFTYPE tcheck;
  int l,k,dir;
  long fpos0;
  FILE *ener_file_temp;
  char dfnam[MAXFILENAME],dfnamback[MAXFILENAME], dfnamtemp[MAXFILENAME];

  // only CPU=0 does anything here
  if(myid==0){
    
    rewind(ener_file);	// go to start
    
    gotit = 0;
    while ((!feof(ener_file)) && (gotit == 0)) {
      
      fscanf(ener_file, "%lf", &tcheck);
      
      if (fabs(tcheck - t) < 0.5 * DTener) {
	gotit = 1;
	for (l = 1; l <= NUMENERVAR; l++) {
	  if ((l > 3+NPR+COMPDIM*2*NPR) && (l < 3+NPR+2*COMPDIM*2*NPR+NPR)) {
	    DIRLOOP PDUMPLOOP {
	      fscanf(ener_file, "%lf", &pcum_tot[dir][k]);
	      l++;
	    }
	    PDUMPLOOP {
	      fscanf(ener_file, "%lf", &fladd_tot[k]);
	      l++;
	    }
	    PDUMPLOOP {
	      fscanf(ener_file, "%lf", &sourceadd_tot[k]);
	      l++;
	    }
	  }
	}
      } else {
	// skip this bad line 
	while ((fgetc(ener_file) != '\n') && (!feof(ener_file)));
      }
      // continue after successful get since successful get is good 
      // data and should keep since corresponds to dump one is
      // keeping
      fpos0 = ftell(ener_file);	// position to continue
      // writting at if successful get
    }
    if (gotit == 0) {
      dualfprintf(fail_file,
		  "Never found right time in loss file when appending: looking for t=%21.15g lastt=%21.15g\n",
		  t, tcheck);
      myexit(1);
    } else {
      dualfprintf(logfull_file,
		  "found goodtime t=%21.15g (wanted %21.15g) to restart ener file\n",
		  tcheck, t);
      sprintf(dfnamtemp, "%s0_ener%s.temp", DATADIR, ".out");
      sprintf(dfnam, "%sener%s", DATADIR, ".out");
      sprintf(dfnamback, "%s0_ener%s.back", DATADIR, ".out");
      
      // now that done, fix up file
      if ((ener_file_temp = fopen(dfnamtemp, "wt")) == NULL) {
	dualfprintf(fail_file,
		    "Cannot open temp ener file for appending: %s\n",
		    dfnamtemp);
	myexit(1);
      } else {
	rewind(ener_file);
	while (ftell(ener_file) < fpos0 + 1) {	// +1 is for
	  // '\n' at end
	  // of line
	  fputc(fgetc(ener_file), ener_file_temp);
	}
	fclose(ener_file_temp);
	fclose(ener_file);
	rename(dfnam, dfnamback);	// move old to backup location
	rename(dfnamtemp, dfnam);	// move new to old name(normal
	// name)
	// reopen loss_file (now normal name)

	if ((ener_file = fopen(dfnam, "at")) == NULL) {
	  dualfprintf(fail_file,
		      "2: error opening ener output file %s\n", dfnam);
	  myexit(1);
	}
	trifprintf("End setup of ener file append\n");
      }			// end else if can open temp file
    }			// end else if gotit==1
  }
}

void divbmaxavg(FTYPE p[][N2M][NPR],SFTYPE*ptrdivbmax,SFTYPE*ptrdivbavg)
{
  int i,j,k;
  int imax=0,jmax=0;
  SFTYPE divb;
  SFTYPE divbmax=0,divbavg=0;
  SFTYPE divbmaxsend,divbavgsend;

  LOOPDIVB {
    // doesn't need geom, just use global gdet
    SETFDIVB(divb, p, i, j);
    if (divb > divbmax) {
      imax = i;
      jmax = j;
      divbmax = divb;
    }
    divbavg += divb;
  }

  // PER CPU
  fprintf(log_file,"  proc: %04d : divbmax: %d %d %21.15g divbavg: %21.15g\n",
	  myid, imax, jmax, divbmax, divbavg / ((SFTYPE) (N1*N2))); fflush(log_file);

#if(USEMPI)			// give CPU=0 total
  divbmaxsend = divbmax;
  divbavgsend = divbavg;
  MPI_Reduce(&divbmaxsend, &divbmax, 1, MPI_SFTYPE, MPI_MAX, 0,
	     MPI_COMM_WORLD);
  MPI_Reduce(&divbavgsend, &divbavg, 1, MPI_SFTYPE, MPI_SUM, 0,
	     MPI_COMM_WORLD);
#endif
  divbavg /= (SFTYPE) (totalzones);

  // Total over all CPUs
  myfprintf(logfull_file,"  divbmax: %21.15g divbavg: %21.15g\n",divbmax, divbavg);
  *ptrdivbmax=divbmax;
  *ptrdivbavg=divbavg;
}





/* gettotal accepts an arbitrary pointer set each of different sizes
 * i.e. one could do:
 * numptrs=2+NPR;
 * totalptrs[0]=pdot;  totalsizes[0]=NPR; totaloptrs[0]=pdot;
 * totalptrs[1]=fladd; totalsizes[1]=NPR; totaloptrs[1]=fladd_tot;
 * PLOOP{ totalptrs[2+k]=&U_tot[k]; totalsizes[k]=1;}
 * gettotal(numptrs,totalptrs,totaloptrs,totalsizes);
 *
 */

void gettotal(int numvars, SFTYPE* vars[],int*sizes,SFTYPE*vars_tot[])
{
  int j,k;
  SFTYPE send;
  
  // for 1 CPU
  if (numprocs == 1) {
    // must use _tot since can't overwrite normal global integrator 
    // variable
    for(k=0;k<numvars;k++){
      for(j=0;j<sizes[k];j++){
	vars_tot[k][j] = vars[k][j];
      }
    }
  } else {
    // give CPU=0 the totals
#if(USEMPI)
    for(k=0;k<numvars;k++){
      for(j=0;j<sizes[k];j++){
	send=vars[k][j];
	// send and receive can't be same address, hence "send" variable
	MPI_Reduce(&send, &vars_tot[k][j], 1, MPI_SFTYPE, MPI_SUM, 0,
		   MPI_COMM_WORLD);
      }
    }
#endif
  }
}

void gettotall(int numvars, CTYPE* vars[],int*sizes,CTYPE *vars_tot[])
{
  int j,k;
  CTYPE send;
  
  // for 1 CPU
  if (numprocs == 1) {
    // must use _tot since can't overwrite normal global integrator 
    // variable
    for(k=0;k<numvars;k++){
      for(j=0;j<sizes[k];j++){
	vars_tot[k][j] = vars[k][j];
      }
    }
  } else {
    // give CPU=0 the totals
#if(USEMPI)
    for(k=0;k<numvars;k++){
      for(j=0;j<sizes[k];j++){
	send=vars[k][j];
	// send and receive can't be same address, hence "send" variable
	MPI_Reduce(&send, &vars_tot[k][j], 1, MPI_CTYPE, MPI_SUM, 0,
		   MPI_COMM_WORLD);
      }
    }
#endif
  }
}


// each CPU does constotal
int constotal(int enerregion, SFTYPE *vars)
{
  int i,j,k;
  FTYPE U[NPR];
  struct of_geom geom;
  struct of_state q;
  FTYPE ftemp[NPR];

  PLOOP vars[k]= 0.0;

  enerpos=enerposreg[enerregion];

  ZLOOP{
    if((i>=enerpos[X1DN])&&(i<=enerpos[X1UP])&&(j>=enerpos[X2DN])&&(j<=enerpos[X2UP])){

      get_geometry(i,j,CENT,&geom) ;
      if(!failed){
	if(get_state(p[i][j],&geom,&q)>=1) return(1);
	if(primtoU(UDIAG,p[i][j],&q,&geom,U)>=1) return(1);
      }
      PLOOP{
	ftemp[k]=U[k]*dVF;
	vars[k] += ftemp[k];
      }
    }
  }
  return(0);
}




// each CPU does constotal
int counttotal(int tscale, int enerregion, CTYPE *vars, int num)
{
  int i,j,k;

  for(k=0;k<num;k++) vars[k]= 0;

  enerpos=enerposreg[enerregion];


  ZLOOP {
    if((i>=enerpos[X1DN])&&(i<=enerpos[X1UP])&&(j>=enerpos[X2DN])&&(j<=enerpos[X2UP])){
      for(k=0;k<num;k++){
	vars[k] += failfloorcount[i][j][tscale][k];
      }
    }
  }
  return(0);
}



#define MAXPTRS 10

int integrate(SFTYPE * var,SFTYPE *var_tot,int type, int enerregion)
{
  SFTYPE *totalptrs[MAXPTRS],*totaloptrs[MAXPTRS];
  int totalsizes[MAXPTRS],numptrs;

  switch(type){
  case CONSTYPE:
    if(constotal(enerregion,var)>=1) return(1);
    totalsizes[0]=NPR;
    gettotal(1,&var,totalsizes,&var_tot); // if only 1 object, then just send address
    break;
  case SURFACETYPE:
  case CUMULATIVETYPE:
    totalsizes[0]=NPR;
    gettotal(1,&var,totalsizes,&var_tot);
    break;
  case CUMULATIVETYPE2:
    totalsizes[0]=1;
    gettotal(1,&var,totalsizes,&var_tot);
    break;
  default:
    dualfprintf(fail_file,"No defined type=%d in integrate\n",type);
    myexit(1);
  }
  return(0);
}

int integratel(CTYPE * var, CTYPE *var_tot,int type, int tscale, int enerregion)
{
  CTYPE *totalptrs[MAXPTRS],*totaloptrs[MAXPTRS];
  int totalsizes[MAXPTRS],numptrs;

  switch(type){
  case CONSTYPE:
    if(counttotal(tscale, enerregion,var,NUMFAILFLOORFLAGS)>=1) return(1);
    totalsizes[0]=NUMFAILFLOORFLAGS;
    gettotall(1,&var,totalsizes,&var_tot);
    break;
  default:
    dualfprintf(fail_file,"No defined type=%d in integratel\n",type);
    myexit(1);
  }
  return(0);
}


void setrestart(int*appendold)
{
  *appendold = 0;
  // 0: deal with ener.out manually
  // 1: append automatically (no longer used)
}
