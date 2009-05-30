
#include "decs.h"

/* diagnostics subroutine */

#define DIAGREPORT {trifprintf("t=%21.15g to do: tener=%21.15g (dt=%21.15g): dump_cnt=%ld @ t=%21.15g (dt=%21.15g) : avg_cnt=%ld @ t=%21.15g (dt=%21.15g) : debug_cnt=%ld @ t=%21.15g (dt=%21.15g) : image_cnt=%ld @ t=%21.15g (dt=%21.15g): restart=%d @ nstep=%ld (dt=%ld)\n",t,tener,DTener,dump_cnt,tdump,DTd,avg_cnt,tavg,DTavg,debug_cnt,tdebug,DTdebug,image_cnt,timage,DTi,whichrestart,nrestart,DTr);}


int diag(int call_code)
{
  //////////////////
  //
  // BEGIN DIAG VARS
  //
  static int firsttime = 1;

  FILE *imagecnt_file, *dumpcnt_file,*avgcnt_file,*debugcnt_file,*fieldlinecnt_file;

  static SFTYPE tlastener,tlastdump,tlastimage,tlastareamap,tlastavg,tlastdebug,tlastfieldline;
  static int dumpc, avgc, imagec, restartc, enerc,debugc,fieldlinec;
  static SFTYPE tdump, tavg, timage, tener,tdebug,tfieldline;
  static long nlastrestart,nrestart;
  int dogdump,doavg, dordump,dodump,doener,doimagedump,doareamap,dodebug,dofieldlinedump;

  int floor,i,j;


  ///////////////////////////
  //
  // setup timing of writes to files
  //

  if ((call_code == INIT_OUT) || (firsttime == 1)) {
    
    // no need to repeat if restarting
    tlastener  = DTener*(SFTYPE)((int)(t/DTener))-SMALL;
    nlastrestart = DTr*(SFTYPE)((realnstep/DTr))-1;
    tlastdump  = DTd*(SFTYPE)((int)(t/DTd))-SMALL;
    tlastavg  = DTavg*(SFTYPE)((int)(t/DTavg))-SMALL;
    tlastimage = DTi*(SFTYPE)((int)(t/DTi))-SMALL;
    tlastdebug  = DTdebug*(SFTYPE)((int)(t/DTdebug))-SMALL;


    tlastareamap = t-SMALL;
    
    dumpc = avgc = imagec = restartc = enerc = debugc = 0;

    if (RESTARTMODE == 0) {
      tdump = timage = tener = t;
      tavg=t+DTavg; // do next time
      nrestart = nstep;

      dump_cnt = 0;
      image_cnt = 0;
      rdump_cnt = 0;
      avg_cnt = 0;
      debug_cnt = 0;

      appendold = 0;
    } else {
      setrestart(&appendold);

      // assuming started at t=0 and nstep=0 for original run
      // time to dump NEXT
      tener  = DTener*(SFTYPE)((int)(t/DTener)+1);
      nrestart = DTr*(SFTYPE)((realnstep/DTr)+1);
      tdump  = DTd*(SFTYPE)((int)(t/DTd)+1);
      timage = DTi*(SFTYPE)((int)(t/DTi)+1);
      tavg  = DTavg*(SFTYPE)((int)(t/DTavg)+1);
      tdebug  = DTdebug*(SFTYPE)((int)(t/DTdebug)+1);
      DIAGREPORT;
    }
  }



  ///////////////////////
  //
  // setup what we will dump this call to diag()
  //

  // output grid (probaly want both fullgrid (to make sure ok) and compute grid to compare with data dumps
  if((DOGDUMPDIAG)&&(!GAMMIEDUMP)&&(firsttime&&(RESTARTMODE==0))){
    dogdump=1;
  }
  else dogdump=0;


  
  if((DODUMPDIAG)&&((t!=tlastdump)&&(t >= tdump || (RESTARTMODE&&dofaildump&&(nstep>=steptofaildump))))){
    dodump=1;
  }
  else dodump=0;

  if((DODEBUG)&&((t!=tlastdebug)&&(t >= tdebug))){
    dodebug=1;
  }
  else dodebug=0;

  
  if((DOAVGDIAG)&&((t!=tlastavg)&&(t >= tavg))){
    doavg=1;
  }
  else doavg=0;
    
  if((DORDUMPDIAG)&&( ((nlastrestart!=nrestart)&&(failed == 0) && (realnstep >= nrestart ))||(call_code==FINAL_OUT) ) ){
    dordump=1;
  }
  else dordump=0;

  // t!=tlast avoids duplicate entries
  if (DOENERDIAG&&(t!=tlastener)&&( (t >= tener)||(call_code==INIT_OUT)||(call_code==FINAL_OUT)||firsttime)){
    doener=1;
  }
  else doener=0;

  if((DOAREAMAPDIAG)&&(t != tlastareamap)&&(dofailmap)&&(nstep>=steptofailmap)){
    doareamap=1;
  }
  else doareamap=0;

  /* image dump at regular intervals */
  if((DOIMAGEDIAG)&&((t!=tlastimage)&&(t >= timage))){
    doimagedump=1;
  }
  else doimagedump=0;

  /* fieldline dump at regular intervals */
  // use image time period
  if((DOFIELDLINE)&&((t!=tlastfieldline)&&(t >= tfieldline))){
    dofieldlinedump=1;
    fieldline_cnt=image_cnt; // force to go with image dump (needed also so restart() knows the count)
  }
  else dofieldlinedump=0;

  //////////////////////////////////////////////////////////////////////
  //
  ////////////////// NOW WRITE TO FILES
  //
  //////////////////////////////////////////////////////////////////////





  //////////////////////
  //
  // Grid dump
  //
  /////////////////////

  if(dogdump) gdump();



  // if doing simulbccalc type loop in step_ch.c then need to bound when doing diagnostics since not done yet
  if(SIMULBCCALC>=1){
    if(dodump||dordump||doener){
      // for dump, rdump, and divb in ener
      bound_prim(-1,p);
    }
  }


  //////////////////////
  //
  // ener dump (integrated quantities: integrate and possibly  dump them too (if doener==1))
  //
  /////////////////////

  if(doener||dordump||(call_code==FINAL_OUT)||(call_code==INIT_OUT)){ // need integratd quantities for restart dump, but don't write them to ener file.

    dump_ener(doener,dordump,call_code);
    if(doener||(call_code==FINAL_OUT)||(call_code==INIT_OUT)){
      frdotout(); // need to include all terms and theta fluxes on horizon/outer edge at some point GODMARK
    }
    while (t >= tener) {
      enerc++;
      tener = enerc * DTener;
    }
    tlastener=t;
  }




  ///////////////////////
  //
  // RESTART DUMP
  //
  ///////////////////////
  

  if(dordump){
    DIAGREPORT;
    trifprintf("dumping: restart: %d nstep: %ld nlastrestart: %ld nrestart: %ld restartc: %d\n", whichrestart,realnstep,nlastrestart,nrestart,restartc);
    restart_write((long)whichrestart);	// 0 1 0 1 0 1 ...
    restartsteps[whichrestart] = realnstep;
    whichrestart = !whichrestart;
    while (realnstep >= nrestart) {
      restartc++;
      nrestart = restartc * DTr;
    }
    nlastrestart=realnstep;
  }
      
  
  ///////////////////////
  //
  // DUMP
  //
  ///////////////////////
  
  if(dodump){
    DIAGREPORT;
    trifprintf("dumping: dump_cnt=%ld : t=%21.15g tlastdump=%21.15g tdump=%21.15g dumpc=%d\n", dump_cnt,t,tlastdump,tdump,dumpc);
    
    /* make regular dump file */
    if (dump(dump_cnt) >= 1){
      dualfprintf(fail_file,"unable to print dump file\n");
      return (1);
    }
    restart_write(-(long)dump_cnt-1);// so can restart at a dump without reconstructing the rdump from a dump.

    // iterate counter
    dump_cnt++;
    while (t >= tdump) {
      dumpc++;
      tdump = dumpc * DTd;
    }
    // output number of dumps
    myfopen("dumps/0_numdumps.dat","w","error opening dump count file\n",&dumpcnt_file);      
    myfprintf(dumpcnt_file, "# Number of dumps\n%ld\n", dump_cnt);
    myfclose(&dumpcnt_file,"Couldn't close dumpcnt_file");

    tlastdump=t;
  }

  ///////////////////////
  //
  // DEBUG DUMP
  //
  ///////////////////////
  
  if(dodebug){
    DIAGREPORT;
    trifprintf("debug dumping: debug_cnt=%ld : t=%21.15g tlastdebug=%21.15g tdebug=%21.15g debugc=%d\n", debug_cnt,t,tlastdebug,tdebug,debugc);
    
    /* make regular dump file */
    if (debugdump(debug_cnt) >= 1){
      dualfprintf(fail_file,"unable to print debug dump file\n");
      return (1);
    }
    // iterate counter
    debug_cnt++;
    while (t >= tdebug) {
      debugc++;
      tdebug = debugc * DTdebug;
    }
    // output number of dumps
    myfopen("dumps/0_numdebug.dat","w","error opening debug dump count file\n",&debugcnt_file);      
    myfprintf(debugcnt_file, "# Number of debug dumps\n%ld\n", debug_cnt);
    myfclose(&debugcnt_file,"Couldn't close debugcnt_file");

    tlastdebug=t;
  }


  ///////////////////////
  //
  // AVG
  //
  ///////////////////////

  if(DOAVGDIAG){
    // do every time step
    // assume can't fail, but can apparently
    if(average_calc(doavg)>=1) return(1);
  }
  
  if(doavg){
    DIAGREPORT
    trifprintf("avging dump: avg_cnt=%ld : t=%21.15g tlastavg=%21.15g tavg=%21.15g avgc=%d\n", avg_cnt,t,tlastavg,tavg,avgc);
    
    /* make avg dump file */
    if (avgdump(avg_cnt) >= 1){
      dualfprintf(fail_file,"unable to print avg file\n");
      return (1);
    }
#if(DOAVG2)
    /* make avg dump file */
    if (avgdump2(avg_cnt) >= 1){
      dualfprintf(fail_file,"unable to print avg2 file\n");
      return (1);
    }
#endif

    // iterate counter
    avg_cnt++;
    while (t >= tavg) {
      avgc++;
      tavg = avgc * DTavg;
    }
    // output number of avgs
    myfopen("dumps/0_numavgs.dat","w","error opening avg count file\n",&avgcnt_file);      
    myfprintf(avgcnt_file, "# Number of avgs\n%ld\n", avg_cnt);
    myfclose(&avgcnt_file,"Couldn't close avgcnt_file");

    tlastavg=t;
  }

  ///////////////////////
  //
  // AREA MAP
  //
  ///////////////////////
  if(doareamap){
    if(area_map(call_code, TIMESERIESAREAMAP, 20, ifail, jfail, p)>=1) return(1);
    tlastareamap=t;
  }

  ///////////////////////
  //
  // IMAGE
  //
  ///////////////////////

  if(doimagedump){
    DIAGREPORT;
    trifprintf("image dump %ld : t=%21.15g tlastimage=%21.15g timage=%21.15g imagec=%d\n", image_cnt, t,tlastimage,timage,imagec);
    
    
    /* make regular image file */
    if(image_dump(image_cnt)>=1) return(1);

    // iterate counter
    image_cnt++;
    while (t >= timage) {
      imagec++;
      timage = imagec * DTi;
    }

    // output number of images
    myfopen("images/0_numimages.dat","w","error opening image count file\n",&imagecnt_file);      
    myfprintf(imagecnt_file, "# Number of images\n%ld\n", image_cnt);
    myfclose(&imagecnt_file,"Couldn't close imagecnt_file");
    
    tlastimage=t;
  }


  ///////////////////////
  //
  // FIELDLINE
  //
  ///////////////////////


  if(dofieldlinedump){
    DIAGREPORT;
    trifprintf("fieldline dump %ld : t=%21.15g tlastfieldline=%21.15g tfieldline=%21.15g fieldlinec=%d\n", fieldline_cnt, t,tlastfieldline,tfieldline,fieldlinec);
    
    // (after processing) equivalent to image in interest
    if(fieldlinedump(fieldline_cnt)>=1) return(1);

    // iterate counter
    fieldline_cnt++;
    while (t >= tfieldline) {
      fieldlinec++;
      tfieldline = fieldlinec * DTi;
    }

    // output number of fieldlines
    myfopen("dumps/0_numfieldlines.dat","w","error opening fieldline count file\n",&fieldlinecnt_file);      
    myfprintf(fieldlinecnt_file, "# Number of fieldlines\n%ld\n", fieldline_cnt);
    myfclose(&fieldlinecnt_file,"Couldn't close fieldlinecnt_file");
    
    tlastfieldline=t;
  }


  ///////////////////////
  //
  // DIAG clearings
  //
  // finally some clearing based upon what called
  //
  ////////////////////////

  if(DODEBUG){ // shouldn't clear these till after ener and debug dump done so both have all timescale data.
    if(doener){
      // cleanse the ener time scale for the failure diag
      ZLOOP FLOORLOOP failfloorcount[i][j][ENERTS][floor]=0;
    }
    if(dodebug){
      // clense failure diag
      ZLOOP FLOORLOOP failfloorcount[i][j][DEBUGTS][floor]=0;
    }
    if(doimagedump){
      // clense the failure counts for this time scale
      ZLOOP FLOORLOOP failfloorcount[i][j][IMAGETS][floor]=0;
    }
  }



  firsttime = 0;
  return (0);
}

/** some diagnostic routines **/

/* map out region around failure point */
int area_map(int call_code, int type, int size, int i, int j, FTYPE prim[][N2M][NPR])
{
  int k;
  int l,m,ll,mm;
  FTYPE r, th, vmin1, vmax1, vmin2, vmax2;
  int ignorecourant;
  struct of_geom geom;
  struct of_state q;
  FTYPE X[NDIM];
  FTYPE divb;
  FTYPE tens_em[NDIM][NDIM], tens_matt[NDIM][NDIM], b[NDIM],
    ucon[NDIM];
  FTYPE U[NPR];
  int lowersizex1,uppersizex1;
  int lowersizex2,uppersizex2;

  static FILE* fileptr;
  static int firsttime=1;
  static int domap=0;
  static int doclose=0;

  trifprintf("\nStart area_map function ... ");

  if(i-(-N1BND)<size/2) lowersizex1=i-(-N1BND);
  else lowersizex1=size/2;
  if((N1-1+N1BND)-i<size/2) uppersizex1=(N1-1+N1BND)-i;
  else uppersizex1=size/2;
  
  if(j-(-N2BND)<size/2) lowersizex2=j-(-N2BND);
  else lowersizex2=size/2;
  if((N2-1+N2BND)-j<size/2) uppersizex2=(N2-1+N2BND)-j;
  else uppersizex2=size/2;


  if(firsttime){
    if((type==TIMESERIESAREAMAP)&&(dofailmap)){
      if((fileptr=fopen("areamap","wt"))==NULL){
	dualfprintf(fail_file,"Cannot open ./areamap on proc=%d\n",myid);
	domap=0;
      }
      else domap=1;
    }
  }

  if((type==TIMESERIESAREAMAP)&&domap&&(call_code==2)){
    doclose=1;
  }
  else doclose=0;

  if(type==FINALTDUMPAREAMAP){
    dualfprintf(fail_file, "area map\n");
    dualfprintf(fail_file, "failure at: i=%d j=%d\n",i+startpos[1],j+startpos[2]);
    coord(i,j,CENT,X);
    dualfprintf(fail_file, "failure at: i=%d j=%d\n",i+startpos[1],j+startpos[2]);



    PLOOP {// all vars
      
      dualfprintf(fail_file, "variable %d \n", k);
      
      dualfprintf(fail_file, "i = \t ");
      for(l=i-lowersizex1;l<=i+uppersizex1;l++){
	ll=l+startpos[1];
	dualfprintf(fail_file, "%21d", ll);
      }
      dualfprintf(fail_file, "\n");
      for(m=j-lowersizex2;m<=j+uppersizex2;m++){
	mm=m+startpos[2];
	dualfprintf(fail_file, "j = %d \t ",mm);
	for(l=i-lowersizex1;l<=i+lowersizex1;l++){
	  ll=l+startpos[1];
	  dualfprintf(fail_file, "%21.15g ",prim[l][m][k]);
	}
	dualfprintf(fail_file, "\n");
      }
    }
  }
  else if((type==TIMESERIESAREAMAP)&&(domap)){
    if(firsttime){
      fprintf(fileptr,"%21.15g %d %d %21.15g %21.15g %21.15g %21.15g %d %d %d %d %d %d %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g\n",
	      t,totalsize[1],totalsize[2],startx[1],startx[2],dx[1],dx[2],lowersizex1,uppersizex1,lowersizex2,uppersizex2,startpos[1]+i,startpos[2]+j,gam,a,R0,Rin,Rout,hslope);
      fflush(fileptr);
    }
    for(m=j-size/2;m<=j+size/2;m++){
      if((m<-N2BND)||(m>N2-1+N2BND)) continue;
      mm=m+startpos[2];
      for(l=i-size/2;l<=i+size/2;l++){
	if((l<-N1BND)||(l>N1-1+N1BND)) continue;
	ll=l+startpos[1];
	
	
	coord(l, m, CENT, X);
	bl_coord(X, &r, &th); 
	get_geometry(l, m, CENT, &geom);
	if (!failed) {
	  if (get_state(p[l][m], &geom, &q) >= 1)
	    FAILSTATEMENT("diag.c:areamap()", "get_state() dir=0", 1);
	  if (vchar(p[l][m], &q, 1, &geom, &vmax1, &vmin1,&ignorecourant) >= 1)
	    FAILSTATEMENT("diag.c:areamap()", "vchar() dir=1or2", 1);
	  if (vchar(p[l][m], &q, 2, &geom, &vmax2, &vmin2,&ignorecourant) >= 1)
	    FAILSTATEMENT("diag.c:areamap()", "vchar() dir=1or2", 2);
	}
	if((l>=-1)&&(l<=N1+1)&&(m>=-1)&&(m<=N2+1)){ SETFDIVB(divb, p, l, m);}
	else divb=0.0;

	// same order as dump.c for first columns (easy sm read)
	fprintf(fileptr,
		"%d %d "
		"%21.15g %21.15g "
		"%21.15g %21.15g "
		"%21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g %21.15g "
		"%21.15g "
		"%21.15g %21.15g %21.15g %21.15g "
		"%21.15g %21.15g %21.15g %21.15g "
		"%21.15g %21.15g %21.15g %21.15g "
		"%21.15g %21.15g %21.15g %21.15g "
		"%21.15g %21.15g %21.15g %21.15g "
		"%21.15g "
		"%21.15g %ld\n",
		ll,mm,
		X[1],X[2],
		r,th,
		p[l][m][0],
		p[l][m][1],
		p[l][m][2],
		p[l][m][3],
		p[l][m][4],
		p[l][m][5],
		p[l][m][6],
		p[l][m][7],
		divb,
		q.ucon[0],q.ucon[1],q.ucon[2],q.ucon[3],
		q.ucov[0],q.ucov[1],q.ucov[2],q.ucov[3],
		q.bcon[0],q.bcon[1],q.bcon[2],q.bcon[3],
		q.bcov[0],q.bcov[1],q.bcov[2],q.bcov[3],
		vmin1,vmax1,vmin2,vmax2,
		geom.g,
		t,realnstep);
      }
    }
    fflush(fileptr);
  }

  if(doclose) if(fileptr!=NULL) fclose(fileptr);


  /* print out other diagnostics here */

  firsttime=0;
  trifprintf("end area_map function.\n");  
  return(0);
}



/* evaluate fluxed based diagnostics; put results in global variables */

#define JETBSQORHO (3.162)


// notice that F1 and F2 have arbitrary eomfunc that must be divided out to get real flux
int diag_flux(FTYPE prim[][N2M][NPR], FTYPE F1[][N2M][NPR], FTYPE F2[][N2M][NPR],SFTYPE Dt)
{
  int fluxdir;
  int i, j, k, l, dir,fl,enerregion;
  SFTYPE surface,surf2;
  SFTYPE surgdet;
  int *iter1,*iter2;
  FTYPE (*flux)[N2M][NPR];
  int start1,start2,stop1,stop2;
  int gpos;
  int ii;
  FTYPE ftemp;
  FTYPE ftemp0,ftemp1,ftemp2,ftemp3,ftemp4,ftemp5,ftemp6;
  FTYPE pgas,bsq,bsqorho;
  FTYPE r, th;
  struct of_geom geom;
  struct of_state q;
  FTYPE X[NDIM];
  FTYPE b[NDIM],ucon[NDIM];
  FTYPE U[NPR];
  int condjet2;
  FTYPE Ftemp[NPR],Ftempdiag[NPR];
  int firstinloop;
  FTYPE pr[NPR];
  extern void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);

  // initialize
  ENERREGIONLOOP{
    doflux=dofluxreg[enerregion];
    enerpos=enerposreg[enerregion];
    pcum=pcumreg[enerregion];
    pdot=pdotreg[enerregion];
    pdotterms=pdottermsreg[enerregion];
    

    DIRLOOP PLOOP{
      pdot[dir][k]=0;
      FLLOOP pdotterms[dir][fl][k]=0;
      if(enerregion==0) FLLOOP pdottermsjet2[dir][fl][k]=0;
    }
    // true outer boundary surface fluxes (per direction, per conserved variable)
    DIRLOOP {
      if (doflux[dir] >= 0) {// then this cpu is involved in flux BOUNDARY calculation
	// otherwise don't add to it
	if((dir==X1UP)||(dir==X1DN)){
	  surface = dx[2]*dx[3];
	  start1=stop1=doflux[dir];
	  start2=enerpos[X2DN];
	  stop2=enerpos[X2UP];
	  flux=F1;
	  gpos=FACE1;
	  fluxdir=1;
	}
	else if((dir==X2UP)||(dir==X2DN)){
	  surface = dx[1]*dx[3];
	  start2=stop2=doflux[dir];
	  start1=enerpos[X1DN];
	  stop1=enerpos[X1UP];
	  flux=F2;
	  gpos=FACE2;
	  fluxdir=2;
	}
	else{
	  dualfprintf(fail_file,"no such direction: %d\n",dir);
	  myexit(1);
	}

	// zero out summation quantities
	PLOOP{
	  pdot[dir][k]=0;
	  FLLOOP pdotterms[dir][fl][k]=0;
	  if(enerregion==0) FLLOOP pdottermsjet2[dir][fl][k]=0;
	}
	GENLOOP(i,j,start1,stop1,start2,stop2){
	  // now add up all zones for each conserved quantity (PLOOP)

	  ////////////
	  //
	  // do standard flux for accounting
	  //
	  ////////////

	  get_geometry(i, j, gpos, &geom);
	  PLOOP Ftemp[k]=flux[i][j][k]*surface; // in UEVOLVE form
	  UtoU(UEVOLVE,UDIAG,&geom,Ftemp,Ftempdiag); // convert to diag form
	  PLOOP pdot[dir][k]  += Ftempdiag[k];




	
	  ////////////
	  //
	  // do term-by-term flux for physics accounting
	  // somewhat like dumps and somewhat like avg of terms
	  //
	  ////////////
	  get_geometry(i, j, CENT, &geom);
	  PLOOP pr[k]=prim[i][j][k];

	  //coord(i, j, CENT, X);
	  //bl_coord(X, &r, &th);
	    // if failed, then data output for below invalid, but columns still must exist    
	    // CENT since p at center
	  if (!failed) {
	    if (get_state(pr, &geom, &q) >= 1)
	      FAILSTATEMENT("diag.c:diag_flux()", "get_state() dir=0", 1);
	  }
	  else {// do a per zone check, otherwise set to 0
	    whocalleducon=1; // force no failure mode, just return like failure, and don't return if failure, just set to 0 and continue
	    if (get_state(pr, &geom, &q) >= 1){
	      for (l = 0; l < NDIM; l++)
		q.ucon[l]=0;
	      for (l = 0; l < NDIM; l++)
		q.ucov[l]=0;
	      for (l = 0; l < NDIM; l++)
		q.bcon[l]=0;
	      for (l = 0; l < NDIM; l++)
		q.bcov[l]=0;
	    }
	    whocalleducon=0; // return to normal state
	  }
	  
	  // somewhat like function which averages stress terms
	  pgas=(gam-1.0)*pr[UU];
	  bsq=0; for(l=0;l<NDIM;l++) bsq+=(q.bcon[l])*(q.bcov[l]);
	  if(enerregion==0){
	    bsqorho=bsq/pr[RHO]; // b^2/\rho
	    // we assume user will check if this condition makes sense for a particular simulation.
	    // condition answer for recording for jet2 region
	    //
	    // a real analysis would trace a field line from the horizon to the outer edge in the funnel-type region and only include the region inside, where eventually we have at the outer edge (or will have) unbound/outbound flow.
	    if(dir==X1DN){
	      condjet2=(bsqorho>JETBSQORHO); // assumes that plunging region never develops such large values.  Can occur, but generally not so.  Can raise if problems.
	    }
	    else if(dir==X1UP){ // assumes in jet2 region, but non-jet2 region could have this property.
	      condjet2=((q.ucon[1]>0.0)&&(-q.ucov[0]-1.0>0.0)); // outgoing and unbound at outer edge
	    }
	    else{
	      condjet2=0;
	    }
	  }
	  else condjet2=0; // never touches pdottermsjet2 then
	  
	  
	  surgdet=(geom.g)*surface;

	  // loop and if's since some calculations are redundantly simple for similar k
	  PLOOP{
	    if(k==RHO){
	      ftemp0=pr[k]*(q.ucon[fluxdir])*surgdet;
	      pdotterms[dir][0][k]+=ftemp0; // only one part
	      if(condjet2) pdottermsjet2[dir][0][k]+=ftemp0; // only one part
	      //	  pdot[dir][k]  += ftemp0 * surgdet;
	    }
	    // part0-6
	    else if((k>=UU)&&(k<=U3)){
	      l=k-UU;
	      // we currently DO NOT add flux[RHO] to flux[UU], just assume reader knows this is from the native stress
	      ftemp0=pgas*(q.ucon[fluxdir])*(q.ucov[l])*surgdet;
	      pdotterms[dir][0][k]+=ftemp0;
	      if(condjet2)	    pdottermsjet2[dir][0][k]+=ftemp0;

	      ftemp1=p[i][j][RHO]*(q.ucon[fluxdir])*(q.ucov[l])*surgdet;
	      pdotterms[dir][1][k]+=ftemp1;
	      if(condjet2)	    pdottermsjet2[dir][1][k]+=ftemp1;


	      ftemp2=p[i][j][UU]*(q.ucon[fluxdir])*(q.ucov[l])*surgdet;
	      pdotterms[dir][2][k]+=ftemp2;
	      if(condjet2)	    pdottermsjet2[dir][2][k]+=ftemp2;


	      ftemp3=bsq*(q.ucon[fluxdir])*(q.ucov[l])*surgdet;
	      pdotterms[dir][3][k]+=ftemp3;
	      if(condjet2)	    pdottermsjet2[dir][3][k]+=ftemp3;


	      ftemp4=pgas*delta(fluxdir,k-UU)*surgdet;
	      pdotterms[dir][4][k]+=ftemp4;
	      if(condjet2)	    pdottermsjet2[dir][4][k]+=ftemp4;


	      ftemp5=0.5*bsq*delta(fluxdir,k-UU)*surgdet;
	      pdotterms[dir][5][k]+=ftemp5;
	      if(condjet2)	    pdottermsjet2[dir][5][k]+=ftemp5;


	      ftemp6=-(q.bcon[fluxdir])*(q.bcov[l])*surgdet;
	      pdotterms[dir][6][k]+=ftemp6;
	      if(condjet2)	    pdottermsjet2[dir][6][k]+=ftemp6;

	    }
	    else if(k==B1){
	      ftemp0=(q.bcon[1])*(q.ucon[fluxdir])*surgdet; // flux_b1 term1
	      pdotterms[dir][0][k]+=ftemp0;
	      if(condjet2)	    pdottermsjet2[dir][0][k]+=ftemp0;


	      ftemp1=-(q.bcon[fluxdir])*(q.ucon[1])*surgdet; // flux_b1 term2
	      pdotterms[dir][1][k]+=ftemp1;
	      if(condjet2)	    pdottermsjet2[dir][1][k]+=ftemp1;

	    }
	    else if(k==B2){
	      ftemp0=(q.bcon[2])*(q.ucon[fluxdir])*surgdet; // flux_b2 term1
	      pdotterms[dir][0][k]+=ftemp0;
	      if(condjet2)	    pdottermsjet2[dir][0][k]+=ftemp0;

	    
	      ftemp1=-(q.bcon[fluxdir])*(q.ucon[2])*surgdet; // flux_b2 term2
	      pdotterms[dir][1][k]+=ftemp1;
	      if(condjet2)	    pdottermsjet2[dir][1][k]+=ftemp1;

	    }
	    else if(k==B3){
	      ftemp0=(q.bcon[3])*(q.ucon[fluxdir])*surgdet; // flux_b3 term1
	      pdotterms[dir][0][k]+=ftemp0;
	      if(condjet2)	    pdottermsjet2[dir][0][k]+=ftemp0;

	      ftemp1=-(q.bcon[fluxdir])*(q.ucon[3])*surgdet; // flux_b3 term2
	      pdotterms[dir][1][k]+=ftemp1;
	      if(condjet2)	    pdottermsjet2[dir][1][k]+=ftemp1;

	    }

	  }// end PLOOP over term-by-term fluxes

	}// end GENLOOP

	// cumulative only based upon pdot
	// pdot contains entire sum over grid of relevant surface integral value for this direction and ener-region.
	PLOOP pcum[dir][k]+=pdot[dir][k]*Dt;

      }// end if doflux
    }// end DIRLOOP
  }// end ENERloop
#if(1)
  // radial flux vs. radius
  flux=F1;
  surface = dx[2]*dx[3];
  for(i=0;i<N1;i++){
    PLOOP frdot[i][k]=0;
    for(j=0;j<N2;j++){
      get_geometry(i, j, FACE1, &geom);
      PLOOP Ftemp[k]=flux[i][j][k]*surface; // UEVOLVE form
      UtoU(UEVOLVE,UDIAG,&geom,Ftemp,Ftempdiag); // convert to diag form
      PLOOP frdot[i][k]+=Ftempdiag[k];
    }
  }
#endif
  // GODMARK
  // want all fluxes vs theta on horizon


  return(0);

}


#define DEBUGFRLOOP 1

// write the flux vs. radius
void frdotout(void)
{
  int i,j,k,l;
  SFTYPE ftemp;
#if(USEMPI)
  MPI_Request rrequest;
  MPI_Request srequest;
#endif
  SFTYPE frdottemp[N1][NPR];
  SFTYPE *frtot;
  int ospos1;
  FILE*frout;
  static int firsttime=1;

  if(numprocs==1){
    frtot=(SFTYPE (*))(&frdot[0][0]);
  }
  else{
#if(USEMPI)
    if(myid==0){
      frtot=(SFTYPE*) malloc(sizeof(SFTYPE)*totalsize[1]*NPR);
      if(frtot==NULL){
	dualfprintf(fail_file,"Cannot get frtot memory\n");
	myexit(1);
      }
      else{
	for(i=0;i<totalsize[1];i++) PLOOP{
	  frtot[i*NPR+k]=0;
	}
      }
      for(l=0;l<numprocs;l++){ // just go over all cpus and assume only added to frdot per cpu for correct cpus.
	ospos1=(l%ncpux1)*N1;
	if(l==0){ // assumes cpu=0 is main cpu and is on horizon
	  for(i=0;i<N1;i++) PLOOP{
	    frdottemp[i][k]=frdot[i][k];
	  }
	}
	else{
	  MPI_Irecv(frdottemp,N1*NPR,MPI_SFTYPE,l,l,MPI_COMM_WORLD,&rrequest);
	  MPI_Wait(&rrequest,&mpichstatus);
	}
	for(i=0;i<N1;i++) PLOOP{
	  frtot[(ospos1+i)*NPR+k]+=frdottemp[i][k];
#if(DEBUGFRLOOP)
	  if((ospos1+i)*NPR+k>=totalsize[1]*NPR){
	    dualfprintf(fail_file,"outside bounds: %d\n",(ospos1+i)*NPR+k);
	    myexit(1);
	  }
#endif
	}
      }
    }
    else{
      MPI_Isend(frdot,N1*NPR,MPI_SFTYPE,0,myid,MPI_COMM_WORLD,&srequest);
      MPI_Wait(&srequest,&mpichstatus);
    }
#endif
  }
  // now we have frtot with full fluxes vs. radius (totalsize[1]), so output

  if(myid==0){
    frout=fopen("frdot.out","at");
    if(frout==NULL){
      dualfprintf(fail_file,"Cannot open frdot.out\n");
      myexit(1);
    }
    if(firsttime){
      fprintf(frout,"%21.15g %ld %d %d %21.15g %21.15g %21.15g %21.15g\n",
	      t,realnstep,totalsize[1],totalsize[2],startx[1],startx[2],dx[1],dx[2]);
      fflush(frout);
    }

    for(i=0;i<totalsize[1];i++){
      fprintf(frout,"%21.15g %d ",t,i);
      PDUMPLOOP{// dump only dump prims
	fprintf(frout,"%21.15g ",frtot[i*NPR+k]);
      }
      fprintf(frout,"\n");
    }
    
    if(frout!=NULL) fclose(frout);
#if(USEMPI)
    free(frtot);
#endif
  }
  firsttime=0;
}


void makedirs(void)
{

  if ((USEMPI == 0) || (USEMPI && (!USEGM))) {
    if ((mpicombine && (myid == 0)) || (mpicombine == 0)) {
      system("mkdir dumps");
      system("mkdir images");
    }
#if(USEMPI)
    MPI_Barrier(MPI_COMM_WORLD);	// all cpus wait for directory
    // to be created
#endif
  }
}





void init_varstavg(void)
{
  int i,j,ii;

  ZLOOP{
    for(ii=0;ii<NUMNORMDUMP;ii++){
      normalvarstavg[i][j][ii]=0.0;
      anormalvarstavg[i][j][ii]=0.0;
    }      
    for(ii=0;ii<NDIM;ii++){
      jcontavg[i][j][ii]=0.0;
      jcovtavg[i][j][ii]=0.0;
      ajcontavg[i][j][ii]=0.0;
      ajcovtavg[i][j][ii]=0.0;
      massfluxtavg[i][j][ii]=0.0;
      amassfluxtavg[i][j][ii]=0.0;
    }
    for(ii=0;ii<NUMOTHER;ii++){
      othertavg[i][j][ii]=0.0;      
      aothertavg[i][j][ii]=0.0;      
    }
    for(ii=0;ii<NUMFARADAY;ii++){
      fcontavg[i][j][ii]=0.0;
      fcovtavg[i][j][ii]=0.0;
      afcontavg[i][j][ii]=0.0;
      afcovtavg[i][j][ii]=0.0;
    }
    for(ii=0;ii<NUMSTRESSTERMS;ii++){
      tudtavg[i][j][ii]=0.0;
      atudtavg[i][j][ii]=0.0;
    }      
    
  }
}

void final_varstavg(FTYPE IDT)
{
  int i,j,ii;

  ZLOOP{
    for(ii=0;ii<NUMNORMDUMP;ii++){
      normalvarstavg[i][j][ii]=normalvarstavg[i][j][ii]*IDT;
      anormalvarstavg[i][j][ii]=anormalvarstavg[i][j][ii]*IDT;
    }      
    for(ii=0;ii<NDIM;ii++){
      jcontavg[i][j][ii]=jcontavg[i][j][ii]*IDT;
      jcovtavg[i][j][ii]=jcovtavg[i][j][ii]*IDT;
      ajcontavg[i][j][ii]=ajcontavg[i][j][ii]*IDT;
      ajcovtavg[i][j][ii]=ajcovtavg[i][j][ii]*IDT;
      massfluxtavg[i][j][ii]*=IDT;
      amassfluxtavg[i][j][ii]*=IDT;
    }
    for(ii=0;ii<NUMOTHER;ii++){
      othertavg[i][j][ii]*=IDT;
      aothertavg[i][j][ii]*=IDT;
    }
    for(ii=0;ii<NUMFARADAY;ii++){
      fcontavg[i][j][ii]=fcontavg[i][j][ii]*IDT;
      fcovtavg[i][j][ii]=fcovtavg[i][j][ii]*IDT;
      afcontavg[i][j][ii]=afcontavg[i][j][ii]*IDT;
      afcovtavg[i][j][ii]=afcovtavg[i][j][ii]*IDT;
    }
    for(ii=0;ii<NUMSTRESSTERMS;ii++){
      tudtavg[i][j][ii]=tudtavg[i][j][ii]*IDT;
      atudtavg[i][j][ii]=atudtavg[i][j][ii]*IDT;
    }      
  }
}

int set_varstavg(FTYPE tfrac)
{
  int i,j,k,l,ii,aii;
  FTYPE ftemp;
  FTYPE ftemp0,ftemp1,ftemp2,ftemp3,ftemp4,ftemp5,ftemp6;
  FTYPE pgas,bsq;
  FTYPE jcov[NDIM];
  FTYPE fcov[NUMFARADAY];
  FTYPE r, th, vmin1, vmax1, vmin2, vmax2;
  int ignorecourant;
  struct of_geom geom;
  struct of_state q;
  FTYPE X[NDIM];
  FTYPE divb;
  FTYPE b[NDIM],ucon[NDIM];
  FTYPE U[NPR];


  ZLOOP{

    // just like dumps
    coord(i, j, CENT, X);
    bl_coord(X, &r, &th);
    // if failed, then data output for below invalid, but columns still must exist    
    get_geometry(i, j, CENT, &geom);
    if (!failed) {
      if (get_state(p[i][j], &geom, &q) >= 1)
	FAILSTATEMENT("diag.c:set_varstavg()", "get_state() dir=0", 1);
      if (vchar(p[i][j], &q, 1, &geom, &vmax1, &vmin1,&ignorecourant) >= 1)
	FAILSTATEMENT("diag.c:set_varstavg()", "vchar() dir=1or2", 1);
      if (vchar(p[i][j], &q, 2, &geom, &vmax2, &vmin2,&ignorecourant) >= 1)
	FAILSTATEMENT("diag.c:set_varstavg()", "vchar() dir=1or2", 2);
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


    ii=0;
    aii=0;

    for(k=0;k<NPR;k++){
      normalvarstavg[i][j][ii++]+=p[i][j][k]*tfrac;
      anormalvarstavg[i][j][aii++]+=fabs(p[i][j][k])*tfrac;
    }
    normalvarstavg[i][j][ii++]+=divb*tfrac;
    anormalvarstavg[i][j][aii++]+=fabs(divb)*tfrac;

    for (k = 0; k < NDIM; k++) normalvarstavg[i][j][ii++]+=q.ucon[k]*tfrac;
    for (k = 0; k < NDIM; k++) anormalvarstavg[i][j][aii++]+=fabs(q.ucon[k])*tfrac;

    for (k = 0; k < NDIM; k++) normalvarstavg[i][j][ii++]+=q.ucov[k]*tfrac;
    for (k = 0; k < NDIM; k++) anormalvarstavg[i][j][aii++]+=fabs(q.ucov[k])*tfrac;

    for (k = 0; k < NDIM; k++) normalvarstavg[i][j][ii++]+=q.bcon[k]*tfrac;
    for (k = 0; k < NDIM; k++) anormalvarstavg[i][j][aii++]+=fabs(q.bcon[k])*tfrac;

    for (k = 0; k < NDIM; k++) normalvarstavg[i][j][ii++]+=q.bcov[k]*tfrac;
    for (k = 0; k < NDIM; k++) anormalvarstavg[i][j][aii++]+=fabs(q.bcov[k])*tfrac;

    normalvarstavg[i][j][ii++]+=vmin1*tfrac;
    anormalvarstavg[i][j][aii++]+=fabs(vmin1)*tfrac;

    normalvarstavg[i][j][ii++]+=vmax1*tfrac;
    anormalvarstavg[i][j][aii++]+=fabs(vmax1)*tfrac;

    normalvarstavg[i][j][ii++]+=vmin2*tfrac;
    anormalvarstavg[i][j][aii++]+=fabs(vmin2)*tfrac;

    normalvarstavg[i][j][ii++]+=vmax2*tfrac;
    anormalvarstavg[i][j][aii++]+=fabs(vmax2)*tfrac;

    lower(jcon[i][j],&geom,jcov);
    for(ii=0;ii<NDIM;ii++){
      jcontavg[i][j][ii]+=jcon[i][j][ii]*tfrac;
      jcovtavg[i][j][ii]+=jcov[ii]*tfrac;
      ajcontavg[i][j][ii]+=fabs(jcon[i][j][ii])*tfrac;
      ajcovtavg[i][j][ii]+=fabs(jcov[ii])*tfrac;
    }
    
    for(ii=0;ii<NDIM;ii++){
      ftemp=(geom.g)*p[i][j][RHO]*(q.ucon[ii]);
      massfluxtavg[i][j][ii]+=ftemp*tfrac;
      amassfluxtavg[i][j][ii]+=fabs(ftemp)*tfrac;
    }

    ii=0;
    aii=0;
    ftemp=(q.ucon[3])/(q.ucon[0]);
    othertavg[i][j][ii++]=ftemp*tfrac;
    aothertavg[i][j][aii++]=fabs(ftemp)*tfrac;

    lowerf(fcon[i][j],&geom,fcov);
    for(ii=0;ii<NUMFARADAY;ii++){
      fcontavg[i][j][ii]+=fcon[i][j][ii]*tfrac;
      fcovtavg[i][j][ii]+=fcov[ii]*tfrac;
      afcontavg[i][j][ii]+=fabs(fcon[i][j][ii])*tfrac;
      afcovtavg[i][j][ii]+=fabs(fcov[ii])*tfrac;
    }

    pgas=(gam-1.0)*p[i][j][UU];
    bsq=0; for(k=0;k<NDIM;k++) bsq+=(q.bcon[k])*(q.bcov[k]);

    // part0
    ii=0;
    for(k=0;k<NDIM;k++) for(l=0;l<NDIM;l++){
      ftemp0=pgas*(q.ucon[k])*(q.ucov[l]);
      tudtavg[i][j][ii]+=ftemp0*tfrac;
      atudtavg[i][j][ii]+=fabs(ftemp0)*tfrac;
      ii++;
    }
    // part1
    for(k=0;k<NDIM;k++) for(l=0;l<NDIM;l++){
      ftemp1=p[i][j][RHO]*(q.ucon[k])*(q.ucov[l]);
      tudtavg[i][j][ii]+=ftemp1*tfrac;
      atudtavg[i][j][ii]+=fabs(ftemp1)*tfrac;
      ii++;
    }
    // part2
    for(k=0;k<NDIM;k++) for(l=0;l<NDIM;l++){
      ftemp2=p[i][j][UU]*(q.ucon[k])*(q.ucov[l]);
      tudtavg[i][j][ii]+=ftemp2*tfrac;
      atudtavg[i][j][ii]+=fabs(ftemp2)*tfrac;
      ii++;
    }
    // part3
    for(k=0;k<NDIM;k++) for(l=0;l<NDIM;l++){
      ftemp3=bsq*(q.ucon[k])*(q.ucov[l]);
      tudtavg[i][j][ii]+=ftemp3*tfrac;
      atudtavg[i][j][ii]+=fabs(ftemp3)*tfrac;
      ii++;
    }
    // part4
    for(k=0;k<NDIM;k++) for(l=0;l<NDIM;l++){
      ftemp4=pgas*delta(k,l);
      tudtavg[i][j][ii]+=ftemp4*tfrac;
      atudtavg[i][j][ii]+=fabs(ftemp4)*tfrac;
      ii++;
    }
    // part5
    for(k=0;k<NDIM;k++) for(l=0;l<NDIM;l++){
      ftemp5=0.5*bsq*delta(k,l);
      tudtavg[i][j][ii]+=ftemp5*tfrac;
      atudtavg[i][j][ii]+=fabs(ftemp5)*tfrac;
      ii++;
    }
    // part6
    for(k=0;k<NDIM;k++) for(l=0;l<NDIM;l++){
      ftemp6=-(q.bcon[k])*(q.bcov[l]);
      tudtavg[i][j][ii]+=ftemp6*tfrac;
      atudtavg[i][j][ii]+=fabs(ftemp6)*tfrac;
      ii++;
    }

  }

  return(0);

}


// if doavg==1, then assume this is call before dumping
int average_calc(int doavg)
{
  static FTYPE lastdt;
  static int calls=0;
  static FTYPE tavgi,tavgf;
  static int tavgflag=1;
 
  if(calls>0){ // since need 2 times

    if(tavgflag){
      // gets reached on next call after dump call or first time
      init_varstavg();
      tavgflag=0;
      tavgi=t;
    }

    // always do
    if(set_varstavg(0.5*(lastdt+dt))>=1) return(1);

    if(doavg==1){
      tavgflag=1;
      tavgf=t;
      final_varstavg(1.0/(tavgf-tavgi));
      // expect to dump after this function ends and before next call to this function
    }
  }
  calls++;
  lastdt=dt;

  return(0);
}


void diag_source(struct of_geom *ptrgeom, FTYPE (*dUcomp)[NPR],SFTYPE Dt)
{
  int sc,k,enerregion;
  FTYPE ftemp[NPR];
  FTYPE ftempdiag[NPR];
  extern void UtoU(int inputtype, int returntype,struct of_geom *ptrgeom,FTYPE *Uin, FTYPE *Uout);

  // does not matter what stage
  if(Dt==dt){

    ENERREGIONLOOP{
      enerpos=enerposreg[enerregion];
      sourceaddterms=sourceaddtermsreg[enerregion];
      sourceadd=sourceaddreg[enerregion];

      if((ptrgeom->i>=enerpos[X1DN])&&(ptrgeom->i<=enerpos[X1UP])&&(ptrgeom->j>=enerpos[X2DN])&&(ptrgeom->j<=enerpos[X2UP])){
	SCLOOP{
	  PLOOP ftemp[k]=Dt*dVF*dUcomp[sc][k]; // in UEVOLVE form
	  UtoU(UEVOLVE,UDIAG,ptrgeom,ftemp,ftempdiag); // convert to diag form

	  // now assign diagnostic form of source
	  PLOOP{
	    sourceaddterms[sc][k]+=ftempdiag[k];
	    sourceadd[k]+=ftempdiag[k];
#if(DOLUMVSR)
	    // only correct for diagonal coordinate Jacobian in which each i is same radius for all j
	    if(k==UU) if(enerregion==0) lumvsr[startpos[1]+ptrgeom->i]+=ftempdiag[k];
#endif
	  } // end PLOOP on diag
	} // end SCLOOP
      }
    }
  }

}
