#include "decs.h"

/* 
   modifications to gammie code: 1) add mympi.h to decs.h 2) add
   mpidecs.h to decs.h 3) add mpidefs.h to defs.h 4) add init_mpi() to
   main.c:main() 5) add bound_mpi() to bounds.c:bound_prim() 6) modify
   set_grid.c for mpi (base *position* on global geometry instead of
   local) 7) modify diag.c for output of files to separate names for
   each cpu 8) "" for rest of fopen (restart.c, postmort.c, etc.) 9)
   modify makefile for MPIability 10) modify step_ch.c for timestep and 
   flux */

void init_mpi(int argc, char *argv[])
{
  int size;

#if(USEMPI)
  // initialize MPI
  workbc =
      (FTYPE(*)[COMPDIM * 2][NPRMPI * NBIGBND * NBIGSM]) (&(workbca[-1][0]
							 [0]));
  workbc_int =
      (int(*)[COMPDIM * 2][NBIGBND * NBIGSM]) (&(workbc_inta[-1][0]
							 [0]));
  init_MPI(argc, argv);
#else
  ncpux1 = 1;
  ncpux2 = 1;
  ncpux3 = 1;

  if(argc==1){}
  else if(argc==3){
    RESTARTMODE=atoi(argv[1]);
    WHICHFILE=atoi(argv[2]);    
  }
  else{
    if(myid==0){
      fprintf(stderr,"<progname>\n");
      fprintf(stderr,"OR\n");
      fprintf(stderr,"<progname> RESTARTMODE WHICHFILE\n");
    }
    exit(1);
  }


  myid = 0;			// defines single process run
  sprintf(myidtxt, "");
  numprocs = 1;
#endif
  if(USEMPI){
    ////////////////////
    //
    // choose to combine files or not
    //
    mpicombine = 1;    // choice
    //mpicombine=0;

    // 
    if(mpicombine){
      if(USEROMIO==0){
	// choice
	if(sortedoutput==SORTED) mpicombinetype=MPICOMBINEMINMEM;
	else if(sortedoutput==UNSORTED) mpicombinetype=MPICOMBINESIMPLE; //forced to happen since no unsorted method for the advanced combine technique
	//mpicombinetype=MPICOMBINESIMPLE; // forced for testing
      }
      else truempicombinetype=mpicombinetype=MPICOMBINEROMIO;
    }
  }
  else{
    // no choice
    mpicombine = 0;
  }
  // always done
  init_genfiles(0);
  init_placeongrid();

#if(USEMPI)
  
  if(sizeof(CTYPE)==sizeof(long long int)){
    MPI_Type_size(MPI_LONG_LONG_INT,&size);
    if(size!=8){
      dualfprintf(fail_file,"size of the long long int in MPI=%d, should be 8\n",size);
      myexit(1000);
    }
  }
  
  if((sizeof(REALTYPE)==sizeof(long double))||(sizeof(SENSITIVE)==sizeof(long double))){
    MPI_Type_size(MPI_LONG_DOUBLE,&size);
    if(size!=16){
      dualfprintf(fail_file,"size of the long double in MPI=%d, should be 16\n",size);
      myexit(1000);
    }
  }
#endif

  trifprintf("done with init_mpi()\n");  fflush(log_file);

}

void myargs(int argc, char *argv[])
{
  if(argc<COMPDIM+1){
    if(myid==0){
      fprintf(stderr,"proc: %04d : Incorrect command line: argc: %d needed=%d, please specify:\n",myid,argc,COMPDIM+1);
      fprintf(stderr,"proc: %04d : mpirun <mpirunoptions> <progname> ncpux1 ncpux2\n",myid);
      fprintf(stderr,"proc: %04d : OR\n",myid);
      fprintf(stderr,"proc: %04d : mpirun <mpirunoptions> <progname> ncpux1 ncpux2 RESTARTMODE WHICHFILE\n",myid);
    }
    exit(1);
  }
  ncpux1=atoi(argv[1]);
  ncpux2=atoi(argv[2]);
  ncpux3=1;
  // default unless user adds below
  RESTARTMODE=0;
  WHICHFILE=0;

  if(argc==COMPDIM+3){
    RESTARTMODE=atoi(argv[3]);
    WHICHFILE=atoi(argv[4]);
  }// no failure stuff since have to recompile anyways for failure tracking
  
}


void init_genfiles(int gopp)
{
  char temps[MAXFILENAME];
  char extension[MAXFILENAME];

  fprintf(stderr, "begin: init_genfiles ... ");
  fflush(stderr);

  if (gopp == 1) {
    strcpy(extension, PPEXT);
  } else if (gopp == 0) {
    strcpy(extension, OUTEXT);
  }
  // always have fail and general log open

  sprintf(temps, "%s0_fail%s%s", DATADIR, extension, myidtxt);


  if ((fail_file = fopen(temps, "at")) == NULL) {
    fprintf(stderr, "fail: Cannot open: %s\n", temps);
    exit(1);
  }
  fprintf(stderr, "opened: %s\n", temps);
  sprintf(temps, "%s0_log%s%s", DATADIR, extension, myidtxt);

  if ((log_file = fopen(temps, "at")) == NULL) {
    fprintf(stderr, "log: Cannot open: %s\n", temps);
    exit(1);
  }
  fprintf(stderr, "opened: %s\n", temps);
  fprintf(log_file, "fail_file: %d log_file: %d\n", (int)fail_file,
	  (int)log_file);
  fflush(log_file);
  if (myid == 0) {
    sprintf(temps, "%s0_logfull%s", DATADIR, extension);

    if ((logfull_file = fopen(temps, "at")) == NULL) {
      fprintf(stderr, "logfull: Cannot open: %s\n", temps);
      exit(1);
    }
    fprintf(stderr, "opened: %s\n", temps);
    fprintf(logfull_file, "logfull_file: %d \n", (int)logfull_file);
    fflush(logfull_file);
  }


  // ok now
  trifprintf("end: init_genfiles\n");
}



#if(USEMPI)

int init_MPI(int argc, char *argv[])
{

  fprintf(stderr, "begin: init_MPI\n");
  fflush(stderr);

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);
  sprintf(myidtxt, CPUTXT, myid);
  MPI_Get_processor_name(processor_name, &procnamelen);

  // currently INIT provides args to rest of processes
  myargs(argc,argv);

  if (MAXCPUS < numprocs) {
    fprintf(stderr,
	    "Must increase MAXCPUS in global.h, %d is too many\n",
	    numprocs);
    myexit(1);
  }

  myfprintf(stderr,
	    "numprocs=%d ncpux1=%d ncpux2=%d percpusize: N1=%d N2=%d\n",
	    numprocs, ncpux1, ncpux2, N1, N2);

  fprintf(stderr, "proc: %s on %s\n", myidtxt, processor_name);

  fprintf(stderr, "end: init_MPI\n");
  fflush(stderr);
  return (0);
}


#endif


void init_placeongrid(void)
{
  int stage,stagei,stagef;
  int i, j, m, l;
  int N[COMPDIM + 1];
  int numbercpu[COMPDIM + 1];
  int dir;
  int opp[COMPDIM*2];

  trifprintf("begin: init_placeongrid ... ");

  N[1] = N1;
  N[2] = N2;
  //  N[3] = N3;

  numbercpu[1] = ncpux1;
  numbercpu[2] = ncpux2;
  //  numbercpu[3] = ncpux3;

  mycpupos[1] = myid % ncpux1;
  mycpupos[2] = (int) ((myid % (numprocs)) / ncpux1);


  for (m = 1; m <= COMPDIM; m++) {
    startpos[m] = mycpupos[m] * N[m];
    endpos[m] = (mycpupos[m] + 1) * N[m] - 1;

    // add up sizes for total size of grid
    totalsize[m] = 0;
    itotalsize[m] = 0;
    for (i = 0; i < numbercpu[m]; i++) {
      totalsize[m] += N[m];
      itotalsize[m] += N[m];
    }
  }

  for(m=1;m<=COMPDIM;m++){
    if((mycpupos0[m]=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate mycpupos0[%d]\n",m);
      myexit(1);
    }
    if((startpos0[m]=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate startpos0[%d]\n",m);
      myexit(1);
    }
    if((endpos0[m]=(int*)malloc(sizeof(int)*numprocs))==NULL){
      dualfprintf(fail_file,"can't allocate endpos0[%d]\n",m);
      myexit(1);
    }
  }
  // for cpu=0 as master to rest, needs this info
  for(i=0;i<numprocs;i++){
    mycpupos0[1][i] = i % ncpux1;
    mycpupos0[2][i] = (int) ((i % (numprocs)) / ncpux1);
    
    for (m = 1; m <= COMPDIM; m++) {
      startpos0[m][i] = mycpupos0[m][i] * N[m];
      endpos0[m][i] = (mycpupos0[m][i] + 1) * N[m] - 1;
    }
  }

  realtotalzones = totalzones =
    totalsize[1] * totalsize[2];
  itotalzones = itotalsize[1] * itotalsize[2];

  /////////////// standard interior MPI data transfer setup
  //
  for(dir=0;dir<COMPDIM*2;dir++) for(j=0;j<DIRNUMVARS;j++){
    dirset[dir][j]=0;
  }
  // see where this cpu needs to send/recv

  // figure out left/right send/recv
  if (mycpupos[1] > 0) {
    dirset[X1DN][DIRIF] = 1;		// do -x1 dir
  }
  if (mycpupos[1] < ncpux1 - 1) {
    dirset[X1UP][DIRIF] = 1;		// do +x1 dir
  }
  // figure out up/down send/recv
  if (mycpupos[2] > 0) {
    dirset[X2DN][DIRIF] = 1;		// -x2 dir
  }
  if (mycpupos[2] < ncpux2 - 1) {
    dirset[X2UP][DIRIF] = 1;		// towards and from +x2 dir
  }

  // only do periodic mpi if 
  if(periodicx1&&(ncpux1>1)){
    if(mycpupos[1]==0) dirset[X1DN][DIRIF]=1;
    else if(mycpupos[1]==ncpux1-1) dirset[X1UP][DIRIF]=1;
  }

  if(periodicx2&&(ncpux2>1)){
    if(mycpupos[2]==0) dirset[X2DN][DIRIF]=1;
    else if(mycpupos[2]==ncpux2-1) dirset[X2UP][DIRIF]=1;
  }

  opp[X1DN]=X1UP;
  opp[X1UP]=X1DN;
  opp[X2DN]=X2UP;
  opp[X2UP]=X2DN;


  
  // tags are defined by sender's ID and direction sent
  // tag=(myid*COMPDIM*2)+{0,1,2,3,4,5}
  // 0=right, 1=up,2=left,3=down,4=out,5=in
  // works/v bc[1=output/2=input][0,1,2,3,4,5]
  // so sends are like: (sendtoid,myid*COMPDIM*2+?) and recv's are
  // like: (fromid,otherid*COMPDIM*2+*) where ? and * are
  // opposites(i.e. 0 and 2, 1 and 3, 4 and 5)

  for(dir=0;dir<COMPDIM*2;dir++){

    // zones to copy from (packing)
    if(dirset[dir][DIRIF]){
      if(dir==X1UP){ // right
	dirset[dir][DIRPSTART1]=N1-N1BND;
	dirset[dir][DIRPSTOP1]=N1-1;
      }
      else if(dir==X1DN){ // left
	dirset[dir][DIRPSTART1]=0;
	dirset[dir][DIRPSTOP1]=N1BND-1;
      }
      if((dir==X1UP)||(dir==X1DN)){
	dirset[dir][DIRPSTART2]=-N2BND;
	dirset[dir][DIRPSTOP2]=N2-1+N2BND;
      }

      if(dir==X2UP){ // up
	dirset[dir][DIRPSTART2]=N2-N2BND;
	dirset[dir][DIRPSTOP2]=N2-1;
      }
      else if(dir==X2DN){ // down
	dirset[dir][DIRPSTART2]=0;
	dirset[dir][DIRPSTOP2]=N2BND-1;
      }
      if((dir==X2UP)||(dir==X2DN)){
	dirset[dir][DIRPSTART1]=-N1BND;
	dirset[dir][DIRPSTOP1]=N1-1+N1BND;
      }

      // sets opposite direction
      dirset[dir][DIROPP]=opp[dir];

      // sets size of transfer
      if((dir==X1UP)||(dir==X1DN)) dirset[dir][DIRSIZE]=N1BND*N2M*NPRMPI;
      else if((dir==X2UP)||(dir==X2DN)) dirset[dir][DIRSIZE]=N2BND*N1M*NPRMPI;

      // matching CPU to transfer to/from
      if((dir==X1UP)||(dir==X1DN)){
	if((periodicx1==0)||((mycpupos[1]>0)&&(mycpupos[1]<ncpux1-1))){
	  if(dir==X1UP) dirset[dir][DIROTHER]=myid+1;
	  if(dir==X1DN) dirset[dir][DIROTHER]=myid-1;
	}
	else if(periodicx1){
	  if(mycpupos[1]==0) dirset[dir][DIROTHER]=myid+(ncpux1-1);
	  else if(mycpupos[1]==ncpux1-1) dirset[dir][DIROTHER]=myid-(ncpux1-1);
	}
      }
      if((dir==X2UP)||(dir==X2DN)){
	if((periodicx2==0)||((mycpupos[2]>0)&&(mycpupos[2]<ncpux2-1))){
	  if(dir==X2UP) dirset[dir][DIROTHER]=myid+ncpux1;
	  if(dir==X2DN) dirset[dir][DIROTHER]=myid-ncpux1;
	}
	else if(periodicx2){
	  if(mycpupos[2]==0) dirset[dir][DIROTHER]=myid+(ncpux2-1)*ncpux1;
	  else if(mycpupos[2]==ncpux2-1) dirset[dir][DIROTHER]=myid-(ncpux2-1)*ncpux1;
	}
      }

      // MPI tags that label transfer, must be unique while doing multiple transfers
      dirset[dir][DIRTAGS]= myid     * COMPDIM * 2 + dir;
      dirset[dir][DIRTAGR]= dirset[dir][DIROTHER] * COMPDIM * 2 + dirset[dir][DIROPP];
    
      // zones to copy into (unpacking)
      if(dir==X1UP){ // right
	dirset[dir][DIRUSTART1]=N1;
	dirset[dir][DIRUSTOP1]=N1-1+N1BND;
      }
      else if(dir==X1DN){ // left
	dirset[dir][DIRUSTART1]=-N1BND;
	dirset[dir][DIRUSTOP1]=-1;
      }
      if((dir==X1UP)||(dir==X1DN)){
	dirset[dir][DIRUSTART2]=-N2BND;
	dirset[dir][DIRUSTOP2]=N2-1+N2BND;
      }
      if(dir==X2UP){ // up
	dirset[dir][DIRUSTART2]=N2;
	dirset[dir][DIRUSTOP2]=N2-1+N2BND;
      }
      else if(dir==X2DN){ // down
	dirset[dir][DIRUSTART2]=-N2BND;
	dirset[dir][DIRUSTOP2]=-1;
      }
      if((dir==X2UP)||(dir==X2DN)){
	dirset[dir][DIRUSTART1]=-N1BND;
	dirset[dir][DIRUSTOP1]=N1-1+N1BND;
      }
    }
  }


  fprintf(log_file,"per: %d %d\n", periodicx1, periodicx2);
  for (m = 1; m <= COMPDIM; m++) {
    fprintf(log_file,"mycpupos[%d]: %d\n", m, mycpupos[m]);
    fprintf(log_file, "startpos[%d]: %d\n", m, startpos[m]);
    fprintf(log_file, "endpos[%d]: %d\n", m, endpos[m]);
    fprintf(log_file, "totalsize[%d]: %d\n", m, totalsize[m]);
  }
  for (m = 0; m < COMPDIM*2; m++) {
    for(l = 0 ; l < DIRNUMVARS ; l++) {
      fprintf(log_file, "dirset[%d][%d]: %d\n", m, l, dirset[m][l]);
    }
  }
  trifprintf("totalzones: %d\n", totalzones);

#if(SIMULBCCALC!=-1)

  if(SIMULBCCALC==2){
    if(SIMULBCCALC<=0){ stagei=STAGEM1; stagef=STAGEM1; }
    else if(SIMULBCCALC==1) { stagei=STAGE0; stagef=STAGE2;}
    else if(SIMULBCCALC==2) { stagei=STAGE0; stagef=STAGE5;}
    
    if(SIMULBCCALC){
      for(stage=stagei;stage<=stagef;stage++){
	STAGECONDITION(0,N1-1,0,N2-1,isc,iec,jsc,jec);
	fprintf(log_file,"CZLOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
	STAGECONDITION(0,N1,-1,N2,isc,iec,jsc,jec);
	fprintf(log_file,"F1LOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
	STAGECONDITION(-1,N1,0,N2,isc,iec,jsc,jec);
	fprintf(log_file,"F2LOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
	STAGECONDITION(0,N1,0,N2,isc,iec,jsc,jec);
	fprintf(log_file,"EMFLOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
	STAGECONDITION(0,N1,0,N2-1,isc,iec,jsc,jec);
	fprintf(log_file,"F1CTLOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
	STAGECONDITION(0,N1-1,0,N2,isc,iec,jsc,jec);
	fprintf(log_file,"F2CTLOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
	STAGECONDITION(-1,N1,-1,N2,isc,iec,jsc,jec);
	fprintf(log_file,"DQLOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
	STAGECONDITION(-2,N1+1,-2,N2+1,isc,iec,jsc,jec);
	fprintf(log_file,"PREDQLOOP: stage=%d : %d %d %d %d\n",stage,isc,iec,jsc,jec);
	fprintf(log_file,"\n");
      }
    }
  }
#endif

  fflush(log_file);




  trifprintf("end: init_placeongrid\n");
}


int myexit(int call_code)
{
  int i, j, k, l;
  int cleanfinish,dofaildump;
  FILE *faildump;
  char mysys[MAXFILENAME];

  trifprintf("proc: %s : Exiting cc: %d nstep: %ld\n", myidtxt,
	  call_code, nstep);



#if(MAILWHENDONE)
    if(myid==0){
      system("echo \"done with `pwd`\" > done.txt");
      if(MAILFROMREMOTE){
	sprintf(mysys,"scp done.txt %s ; ssh %s \"mail %s < done.txt\"",REMOTEHOST,REMOTEHOST,EMAILADDRESS);
	system(mysys);
      }
      else{
	sprintf(mysys,"mail %s < done.txt",EMAILADDRESS);
	system(mysys);
      }
    }
#endif

  if (call_code >= 0) {
    if (fail_file)
      fclose(fail_file);
    if (log_file)
      fclose(log_file);
    myfclose(&logfull_file,"Can't close logfull_file\n");
  }
  dofaildump=0;
  if (call_code > 0) {
    fprintf(stderr,
	    "proc: %s : Failure.  Please check failure file: cc: %d\n",
	    myidtxt, call_code);

    if(call_code<1000) cleanfinish = 1;
    else cleanfinish=0; // assume this means dump procedure failed, so don't get into infinite failure loop
    // should never have non-clean finish, but sometimes do have them in code, but not marked right now
    if(cleanfinish) dofaildump=1;
    if(!cleanfinish){
#if(USEMPI)
      // must abort since no clear to communicate to other cpus now
      MPI_Abort(MPI_COMM_WORLD, 1);
#endif
    }
  }
  else{
    dofaildump=0;
    cleanfinish=1;
  }

  if (dofaildump) {
    fprintf(stderr, "proc: %s : dumping failure dump with callcode=2\n",
	    myidtxt);
      
    // assume want previous timestep data, not bad just-computed
    // data\n");
    // now diag should not fail if last timestep was non-fail type
    if (DODIAGS)
      diag(2);
  }
  if(cleanfinish){
    fprintf(stderr,
	    "Ending Computation on proc: %s, holding for other cpus\n",
	    myidtxt);

#if(USEMPI)
    // finish up MPI
    MPI_Barrier(MPI_COMM_WORLD);	// required!
    MPI_Finalize();
#endif
    
    myfprintf(stderr, "Ended Computation on all processors\n");
  }    
  fprintf(stderr, "END\n");
  fflush(stderr);
  exit(0);
  return (0);
}

// note, this may be called in different locations of the code by
// different CPUs
int error_check(int wherefrom)
{
  int i, j, k;
  int errorsend = 0;
  // check if error exists and exit if so

  if (failed > 0) {
    dualfprintf(fail_file,
	    "Detected failure on proc: %d failed: %d nstep: %ld realnstep: %ld t: %21.15g wherefrom = %d\n",
	    myid, failed, nstep, realnstep, t,wherefrom);
  }

  if (numprocs > 1) {
    errorsend = failed;
#if(USEMPI)
    // fprintf(fail_file,"wtf: %d %d\n",errorsend,failed);
    // fflush(fail_file);
    MPI_Allreduce(&errorsend, &failed, 1, MPI_INT, MPI_MAX,
		  MPI_COMM_WORLD);
    // fprintf(fail_file,"wtf: %d %d\n",errorsend,failed);
    // fflush(fail_file);
#endif
  }
  if (failed > 0) {
    dualfprintf(fail_file,
	    "Result: Detected failure on proc: %d failed: %d nstep: %ld realnstep: %ld t: %21.15g\n",
	    myid, failed, nstep, realnstep, t);
    // control behavior of failure here (i.e. could return(1) and
    // continue or something)
    // if(failed==1) myexit(1);
    // if(failed==2) myexit(1);
    // if(failed==3) myexit(1);
    myexit(1);
    return (1);
  }
  return (0);
}


void mpiio_init(int bintxt, int sorted, FILE ** fpptr, long headerbytesize, int which, char *filename, int numcolumns,
		MPI_Datatype datatype, void **jonioptr, void **writebufptr)
{

  logsfprintf("\nmpiio_init begin\n");

  // this check covers combine and seperate
  if(!sorted){
    if(bintxt==BINARYOUTPUT){
      dualfprintf(fail_file,"No such thing as binary unsorted output\n");
      myexit(1);
    }
  }

  if(mpicombinetype==MPICOMBINEROMIO){
    // "true" value same
    // doesn't use jonioptr or bintxt or sorted or fpptr
    // writebuf not set yet for combine or seperate, have to allocate at address writebufptr first.
    mpiioromio_init_combine(INITROMIO, which, headerbytesize, filename,numcolumns, datatype,writebufptr,(void*)0x0);
  }
  else{
    mpiios_init(bintxt, sorted, fpptr, which, headerbytesize, filename, numcolumns, datatype,
		jonioptr, writebufptr);
  }

  logsfprintf("mpiio_init end\n");

}


void mpiio_combine(int bintxt, int sorted,
		   int numcolumns, MPI_Datatype datatype,
		   FILE ** fpptr, void *jonio, void *writebuf)
{
#if(USEMPI)

  logsfprintf("mpiio start combine\n");

  if(USEROMIO){
    // doesn't use jonioptr or bintxt or sorted or fpptr
    // address not used, just value of writebuf
    // headerbytesize no longer needed
    mpiioromio_init_combine(WRITECLOSEROMIO, WRITEFILE,  0, "", numcolumns, datatype,&writebuf,writebuf);
  }
  else{
    if(truempicombinetype==MPICOMBINESIMPLE){
      
      if (sorted) {
	mpiios_combine(bintxt, datatype, numcolumns,fpptr, jonio, writebuf);
      }
      else{
	mpiiotu_combine(datatype, numcolumns, fpptr, writebuf);
      }
  }
    else if(truempicombinetype==MPICOMBINEMINMEM){
      mpiiomin_final(numcolumns,fpptr, jonio, writebuf);
    }
  }
#endif
  logsfprintf("mpiio end combine\n");

}

#define DEBUGMINMEM 0


// this can be made to work without mpi for super-efficient single cpu mode, but not important.  Can still do 1 cpu with MPI and use this!
void mpiio_minmem(int readwrite, int whichdump, int i, int j, int bintxt, int sorted,
		   int numcolumns, MPI_Datatype datatype,
		   FILE ** fp, void *jonio, void *writebuf)
{
  static struct blink * blinkptr;
  static struct blink * cpulinkptr;
  static int datasofar0,datasofarc0;
  static int done0=0;
  static int datagot;
  int bfi;
  int sii,uii;// sorted and unsorted index sometimes
#if(USEMPI)
  MPI_Request request;
  MPI_Request request0;
#endif
  int joniooffset;
  void *joniosubmit;

  unsigned char *jonio1;
  float *jonio4;
  double *jonio8;
  long double *jonio16;
  int *jonio4i;
  long long int *jonio8i;

  int mapjoniosorted,mapjoniounsorted;
  int datatoget0,datatogetc0,datasent0;
  int doset;
  int dofull;
  int dolocal;
  int doing0;
  int gi,gj;//global i,j starting from first cpu reference and adding the sorted index
  int lastnum;
  static int masterdone;
  static int thisslavedone;

  int numfiles;
  int sizeofdatatype;

#if(USEMPI)

  sizeofdatatype=getsizeofdatatype(datatype);

#if(DEBUGMINMEM)
  fprintf(fail_file,"got here 0: i=%d j=%d\n",i,j); fflush(fail_file);
#endif

  if(sorted==UNSORTED){
    dualfprintf(fail_file,"Not setup to do unsorted with this method\n");
    myexit(10000);
  }

  // very quickly determine whether at required point to begin read or write of data, or still busy with buffer.
  if((i==0)&&(j==0)){
    // then first time here for this cpu for this dump
    blinkptr=blinkptr0[whichdump];
    datagot=0;
    if(myid==0) masterdone=0;
    thisslavedone=0;
  }
  if(masterdone&&(myid==0)) return;
  if(thisslavedone&&(myid!=0)) return;


  if(readwrite==WRITEFILE){
    // i+1 since called inside loop before i++, but really has done i++ things
    bfi=j*N1*numcolumns+(i+1)*numcolumns; // recall that we are always at col=0 since we always complete a column
  }
  else if(readwrite==READFILE){
    bfi=j*N1*numcolumns+i*numcolumns;
  }


  if(blinkptr==NULL) dolocal=0; // end of local list
  else{
    // see if completed a single node in the list yet
    if(readwrite==WRITEFILE){
      if(bfi==(datagot+blinkptr->num)) dolocal=1; // at a completion point, need to send to cpu=0
      else   dolocal=0; // not yet completed
    }
    else if(readwrite==READFILE){
      if(bfi==datagot) dolocal=1; // at a completion point, need to recv from cpu=0
      else   dolocal=0; // not yet completed
    }
  }

  if(dolocal){
#if(DEBUGMINMEM)
    fprintf(fail_file,"got here 0.5: %d %d %d %d\n",bfi,datagot+blinkptr->num,blinkptr->num,dolocal); fflush(fail_file);
#endif
    // We are ready to read or write from/to cpu=0
#if(DEBUGMINMEM)
    fprintf(fail_file,"got here 0.6\n"); fflush(fail_file);
#endif
    
    if(readwrite==WRITEFILE){
      // means we have put what we want into writebuf, now send to cpu=0
#if(DEBUGMINMEM)
      fprintf(fail_file,"got here 0.65 : %d\n",writebuf); fflush(fail_file);
#endif
      // must keep myid>0 cpus stalled until cpu=0 needs their data
      // that is, Wait below continues if data copied out of buffer, but we want pseudo-blocking call here
      if(myid>0) MPI_Issend(writebuf,blinkptr->num,datatype,0,myid,MPI_COMM_WORLD,&request);
      // can't stall cpu=0 since only below has recv, but ok since cpu=0 stuck below until cpu=0 data needed again
      else MPI_Isend(writebuf,blinkptr->num,datatype,0,myid,MPI_COMM_WORLD,&request);
#if(DEBUGMINMEM)
      fprintf(fail_file,"got here 0.66 : %d\n",writebuf); fflush(fail_file);
#endif
      // adjust offset to keep standard writebuf BUFFERMAP format with just an offset (offset keeps true memory area always as written part)
      bufferoffset-=blinkptr->num;// or =-datagot+blinkptr->num
    }
    else if(readwrite==READFILE){
      // means we need to fill writebuf with what cpu=0 will send us
      // recv's will Wait below until data is coming.
#if(DEBUGMINMEM)
      fprintf(fail_file,"got here 0.65 : %d\n",writebuf); fflush(fail_file);
#endif
      MPI_Irecv(writebuf,blinkptr->num,datatype,0,myid,MPI_COMM_WORLD,&request);      
#if(DEBUGMINMEM)
      fprintf(fail_file,"got here 0.66 : %d %d\n",writebuf,blinkptr); fflush(fail_file);
#endif
      bufferoffset=-datagot;
    }
    datagot+=(blinkptr->num);
    lastnum=blinkptr->num;
    // now iterate to next node in list
    blinkptr=(blinkptr->np);
#if(DEBUGMINMEM)
    fprintf(fail_file,"got here 0.67 : %d\n",blinkptr); fflush(fail_file);
#endif
    if(blinkptr==NULL){
      thisslavedone=1;
      // then done with list, check to see if everything adds up
      if(
	 ((readwrite==WRITEFILE)&&(-bufferoffset!=N1*N2*numcolumns))
	||((readwrite==READFILE)&&(-bufferoffset!=N1*N2*numcolumns-lastnum))
	 )
	{
	dualfprintf(fail_file,"local total doesn't equal expected value\n got=%d demand=%d\n",-bufferoffset,N1*N2*numcolumns);
	myexit(10000);
      }
    }
#if(DEBUGMINMEM)
    fprintf(fail_file,"got here 0.7\n"); fflush(fail_file);
#endif

  }

#if(DEBUGMINMEM)
  fprintf(fail_file,"got here .8\n"); fflush(fail_file);
#endif

  ///////////////////////
  // CPU==0 stuff
  
  // for cpu=0, don't wait on local send/recv since cpu=0 needs to setup it's receives first (since cpu=0 sends to itself)
  // only do this if cpu=0 is just done with local data send/recv(i.e. dolocal==1) or cpu=0 is done (done0==1)
  // no way to decouple cpu=0 local from master, so master waits for local
  if(myid==0){
    // first time hit, get initial node in list and reset counters
    if((i==0)&&(j==0)){
      cpulinkptr=cpulinkptr0[whichdump];
      datasofar0=0;
      datasofarc0=0;
    }
  }


  // check to see if done with cpu=0 data (we use cpu=0 to continue grabbing data as needed)
  if(readwrite==WRITEFILE){
    // wait until all data is grabbed from cpu=0 to send to disk before holding up here.
    if((myid==0)&&(bfi==N1*N2*numcolumns)){
      done0=1;
    }
    else done0=0;
  }
  else if(readwrite==READFILE){
    // wait until last read done.  Still can't use writebuf since will do last local cpu=0 assignment using writebuf after done with distribution of rest of data using cpu=0 to other cpus
    if((myid==0)&&(thisslavedone)){
      done0=1; // then will be done after final read and grab for cpu=0 (then cpu=0's writebuf will be written to after all cpus done)
    }
    else done0=0;
  }

  if((myid==0)&&((dolocal==1)||(done0==1))){
    // We are ready to read or write from/to other cpus

    if (datatype == MPI_UNSIGNED_CHAR) jonio1 = (unsigned char *) jonio;
    else if (datatype == MPI_FLOAT) jonio4 = (float *) jonio;
    else if (datatype == MPI_DOUBLE) jonio8 = (double *) jonio;
    else if (datatype == MPI_LONG_DOUBLE) jonio16 = (long double *) jonio;
    else if (datatype == MPI_INT) jonio4i = (int *) jonio;
    else if (datatype == MPI_LONG_LONG_INT) jonio8i = (long long int *) jonio;



    // we need to loop over list until hit a marked end-cpu.  This will give us a chunk of data to sort for direct output
    // this is the offset to the sorted part of jonio
    joniooffset=joniosize/2;

#if(DEBUGMINMEM)
    fprintf(fail_file,"got here 2\n"); fflush(fail_file);
#endif

    if(readwrite==WRITEFILE){
      dofull=1;
      while(dofull){
#if(DEBUGMINMEM)
	fprintf(fail_file,"got here 3\n"); fflush(fail_file);
#endif
	/////////
	// determine amount of data to get from cpu group
	
	// limit the amount of data once doing last grab
	if(datasofar0+joniosize/2>totalsize[1]*totalsize[2]*numcolumns) datatogetc0=-datasofar0+totalsize[1]*totalsize[2]*numcolumns;
	else datatogetc0=joniosize/2;
	// new amount of data (that will be read total)
	datasofarc0+=datatogetc0;      


	doset=1;
	doing0=0;
	// every receieved dataset is kept in first half of jonio from start of jonio, sorted, then next data set received.
	datatoget0=0;
	while(doset){
#if(DEBUGMINMEM)
	  fprintf(fail_file,"got here 3.4 : cpu=%d\n",cpulinkptr->cpu); fflush(fail_file);
#endif
	  datatoget0+=cpulinkptr->num;
	  if(cpulinkptr->cpu==0) doing0=1;
#if(DEBUGMINMEM)
	  fprintf(fail_file,"got here 3.5: %d %d %d\n",jonio,cpulinkptr->num,cpulinkptr->cpu); fflush(fail_file);
	  fprintf(fail_file,"got here 3.51: %d %d\n",datatogetc0,totalsize[1]); fflush(fail_file);
#endif
	  MPI_Irecv(jonio,cpulinkptr->num,datatype,cpulinkptr->cpu,cpulinkptr->cpu,MPI_COMM_WORLD,&request0);
	  // have myid wait before continuing to make sure receive is complete
	  MPI_Wait(&request0,&mpichstatus);
	  // easiest algorithm/mapping is done using loop over full sorted size, reduced per cpu by if statement and checked
	  //	  for(sii=0;sii<cpulinkptr->num;sii++){
	  uii=0;
	  for(sii=0;sii<datatogetc0;sii++){
	    gi=(sii/numcolumns+cpulinkptr->ri)%totalsize[1];
	    gj=(sii/numcolumns+(cpulinkptr->rj)*totalsize[1]+(cpulinkptr->ri))/totalsize[1];
#if(DEBUGMINMEM)
	    fprintf(fail_file,"got here 3.55: %d %d %d %d %d\n",sii, gi,gj,cpulinkptr->ri,cpulinkptr->rj); fflush(fail_file);
#endif

	    if(
	       (gi>=startpos0[1][cpulinkptr->cpu])&&
	       (gi<=endpos0[1][cpulinkptr->cpu])&&
	       (gj>=startpos0[2][cpulinkptr->cpu])&&
	       (gj<=endpos0[2][cpulinkptr->cpu])
	       ){
	      if (datatype == MPI_UNSIGNED_CHAR) jonio1[sii+joniooffset]=jonio1[uii++];
	      else if (datatype == MPI_FLOAT) jonio4[sii+joniooffset]=jonio4[uii++];
	      else if (datatype == MPI_DOUBLE) jonio8[sii+joniooffset]=jonio8[uii++];
	      else if (datatype == MPI_LONG_DOUBLE) jonio16[sii+joniooffset]=jonio16[uii++];
	      else if (datatype == MPI_INT) jonio4i[sii+joniooffset]=jonio4i[uii++];
	      else if (datatype == MPI_LONG_LONG_INT) jonio8i[sii+joniooffset]=jonio8i[uii++];
	    }
	  }
#if(DEBUGMINMEM)
	  fprintf(fail_file,"got here 3.6\n"); fflush(fail_file);
#endif
	  // check!
	  // see if local data total is as expected
	  if(uii!=cpulinkptr->num){
	    dualfprintf(fail_file,"uii and num for this cpu not same, algorithm failure: uii=%d num=%d datatogetc0=%d\n",uii,cpulinkptr->num,datatogetc0);
	    myexit(10000);
	  }
	  if(cpulinkptr->end) doset=0;
	  cpulinkptr=cpulinkptr->np;
	}
	// see if total data is as expected
	datasofar0+=datatoget0; // diagnostic
	if(datasofar0!=datasofarc0){
	  dualfprintf(fail_file,"cumulative data received via MPI and expected data is different: got=%d expected=%d\n",datasofar0,datasofarc0);
	  myexit(10000);
	}

	// now take sorted part and write it to disk
	// now write out collected data using CPU=0
	/////////////
	// send data to file
	if (datatype == MPI_UNSIGNED_CHAR) joniosubmit=(void*) (jonio1+joniooffset);
	else if (datatype == MPI_FLOAT) joniosubmit=(void*) (jonio4+joniooffset);
	else if (datatype == MPI_DOUBLE) joniosubmit=(void*) (jonio8+joniooffset);
	else if (datatype == MPI_LONG_DOUBLE) joniosubmit=(void*) (jonio16+joniooffset);
	else if (datatype == MPI_INT) joniosubmit=(void*) (jonio4i+joniooffset);
	else if (datatype == MPI_LONG_LONG_INT) joniosubmit=(void*) (jonio8i+joniooffset);

	if(docolsplit){
	  numfiles=numcolumns;
	}
	else numfiles=1;

	if(bintxt==BINARYOUTPUT){
	  if(numfiles==1) fwrite(joniosubmit, sizeofdatatype,datatoget0, *fp);
	  else{
	    for(sii=0;sii<datatoget0;sii++){
	      if (datatype == MPI_UNSIGNED_CHAR) fwrite(&jonio1[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
	      else if (datatype == MPI_FLOAT) fwrite(&jonio4[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
	      else if (datatype == MPI_DOUBLE) fwrite(&jonio8[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
	      else if (datatype == MPI_LONG_DOUBLE) fwrite(&jonio16[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
	      else if (datatype == MPI_INT) fwrite(&jonio4i[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
	      else if (datatype == MPI_LONG_LONG_INT) fwrite(&jonio8i[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
	    }
	  }
	}
	else if(bintxt==TEXTOUTPUT){ // properly ordered, so just dump it
	  for(sii=0;sii<datatoget0;sii++){
	    if (datatype == MPI_UNSIGNED_CHAR) fprintf(fp[sii%numfiles],"%c",jonio1[sii+joniooffset]);
	    else if (datatype == MPI_FLOAT) fprintf(fp[sii%numfiles],"%15.7g",jonio4[sii+joniooffset]);
	    else if (datatype == MPI_DOUBLE) fprintf(fp[sii%numfiles],"%21.15g",jonio8[sii+joniooffset]);
	    else if (datatype == MPI_LONG_DOUBLE) fprintf(fp[sii%numfiles],"%31.25Lg",jonio16[sii+joniooffset]);
	    else if (datatype == MPI_INT) fprintf(fp[sii%numfiles],"%d",jonio4i[sii+joniooffset]);
	    else if (datatype == MPI_LONG_LONG_INT) fprintf(fp[sii%numfiles],"%lld",jonio8i[sii+joniooffset]);
	    if(numfiles==1){
	      if((sii+1)%numcolumns) fprintf(fp[sii%numfiles]," ");
	      else fprintf(fp[sii%numfiles],"\n");
	    }
	    else fprintf(fp[sii%numfiles],"\n");
	  }
	}
	if((done0==0)&&(doing0==1)) dofull=0; // cpu==0 still needs to deal with it's own data
	// if wasn't locally doing cpu=0 by master, then can't exit yet, continue till cpu=0 local data needed by master
	if(cpulinkptr==NULL){
	  dofull=0; // end of list
	  masterdone=1;
	  if(datasofar0!=totalsize[1]*totalsize[2]*numcolumns){
	    dualfprintf(fail_file,"write: total data written not equal to expected amount: wrote=%d expected=%d\n",datasofar0,totalsize[1]*totalsize[2]*numcolumns);
	    myexit(10000);
	  }
	}
	// otherwise continue
      }
    }
    else if(readwrite==READFILE){
      dofull=1;
      while(dofull){
#if(DEBUGMINMEM)
	fprintf(fail_file,"got here 4\n"); fflush(fail_file);
#endif
	// we read data into 2nd half of jonio (sorted part), then desort for one cpu, then continue for next cpu 
	
	/////////
	// determine amount of data to get from file


	// limit the amount of data once doing last grab
	if(datasofar0+joniosize/2>totalsize[1]*totalsize[2]*numcolumns) datatogetc0=-datasofar0+totalsize[1]*totalsize[2]*numcolumns;
	else datatogetc0=joniosize/2;
	// new amount of data (that will be read total)
	datasofarc0+=datatogetc0;      

#if(DEBUGMINMEM)
	fprintf(fail_file,"got here1 : %d %d\n",datatogetc0,datasofarc0); fflush(fail_file);
#endif

	/////////////
	// get data from file
	if (datatype == MPI_UNSIGNED_CHAR) joniosubmit=(void*) (jonio1+joniooffset);
	else if (datatype == MPI_FLOAT) joniosubmit=(void*) (jonio4+joniooffset);
	else if (datatype == MPI_DOUBLE) joniosubmit=(void*) (jonio8+joniooffset);
	else if (datatype == MPI_LONG_DOUBLE) joniosubmit=(void*) (jonio16+joniooffset);
	else if (datatype == MPI_INT) joniosubmit=(void*) (jonio4i+joniooffset);
	else if (datatype == MPI_LONG_LONG_INT) joniosubmit=(void*) (jonio8i+joniooffset);
	

	if(docolsplit){
	  numfiles=numcolumns;
	}
	else numfiles=1;

	if(bintxt==BINARYOUTPUT){
	  // first let cpu=0 read data
	  if(numfiles==1) fread(joniosubmit, sizeofdatatype,datatogetc0,*fp);
	  else{
	    for(sii=0;sii<datatoget0;sii++){
	    if (datatype == MPI_UNSIGNED_CHAR) fread(&jonio1[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
	    else if (datatype == MPI_FLOAT) fread(&jonio4[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
	    else if (datatype == MPI_DOUBLE) fread(&jonio8[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
	    else if (datatype == MPI_LONG_DOUBLE) fread(&jonio16[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
	    else if (datatype == MPI_INT) fread(&jonio4i[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
	    else if (datatype == MPI_LONG_LONG_INT) fread(&jonio8i[sii+joniooffset], sizeofdatatype,1, fp[sii%numfiles]);
	    }
	  }
	}
	else if(bintxt==TEXTOUTPUT){ // properly ordered, so just grab it
	  for(sii=0;sii<datatogetc0;sii++){
	    if (datatype == MPI_UNSIGNED_CHAR) fscanf(fp[sii%numfiles],"%uc",&jonio1[sii+joniooffset]);
	    else if (datatype == MPI_FLOAT) fscanf(fp[sii%numfiles],"%f",&jonio4[sii+joniooffset]);
	    else if (datatype == MPI_DOUBLE) fscanf(fp[sii%numfiles],"%lf",&jonio8[sii+joniooffset]);
	    else if (datatype == MPI_LONG_DOUBLE) fscanf(fp[sii%numfiles],"%Lf",&jonio16[sii+joniooffset]);
	    else if (datatype == MPI_INT) fscanf(fp[sii%numfiles],"%d",&jonio4i[sii+joniooffset]);
	    else if (datatype == MPI_LONG_LONG_INT) fscanf(fp[sii%numfiles],"%lld",&jonio8i[sii+joniooffset]);
	  }
	}
#if(DEBUGMINMEM)
	fprintf(fail_file,"got here2\n"); fflush(fail_file);
#endif
	
	/////////////////
	// send data to other cpus
	
	// now that data is read into 2nd half of jonio, need to desort data into 1st half per cpu in a loop
	datasent0=0;
	datatoget0=0;
	doset=1;
	doing0=0;

	while(doset){
#if(DEBUGMINMEM)
	  fprintf(fail_file,"got here 2.4\n"); fflush(fail_file);
#endif

	  if(cpulinkptr->cpu==0) doing0=1;
	  datatoget0+=cpulinkptr->num;

	  // copy from 2nd half of jonio to first half the node cpu's data
	  uii=0;
	  for(sii=0;sii<datatogetc0;sii++){
#if(DEBUGMINMEM)
	    fprintf(fail_file,"got here2.5: %d\n",sii); fflush(fail_file);
#endif
	    gi=(sii/numcolumns+cpulinkptr->ri)%totalsize[1];
	    gj=(sii/numcolumns+(cpulinkptr->rj)*totalsize[1]+(cpulinkptr->ri))/totalsize[1];
#if(DEBUGMINMEM)
	    fprintf(fail_file,"got here2.6: %d %d\n",gi,gj); fflush(fail_file);
#endif
	    if(
	       (gi>=startpos0[1][cpulinkptr->cpu])&&
	       (gi<=endpos0[1][cpulinkptr->cpu])&&
	       (gj>=startpos0[2][cpulinkptr->cpu])&&
	       (gj<=endpos0[2][cpulinkptr->cpu])
	       ){
#if(DEBUGMINMEM)
	      fprintf(fail_file,"got here2.7 %d\n",cpulinkptr->cpu); fflush(fail_file);
#endif

	      if (datatype == MPI_UNSIGNED_CHAR) jonio1[uii++]=jonio1[sii+joniooffset];
	      else if (datatype == MPI_FLOAT) jonio4[uii++]=jonio4[sii+joniooffset];
	      else if (datatype == MPI_DOUBLE) jonio8[uii++]=jonio8[sii+joniooffset];
	      else if (datatype == MPI_LONG_DOUBLE) jonio16[uii++]=jonio16[sii+joniooffset];
	      else if (datatype == MPI_INT) jonio4i[uii++]=jonio16[sii+joniooffset];
	      else if (datatype == MPI_LONG_LONG_INT) jonio8i[uii++]=jonio16[sii+joniooffset];
	    }
	  }
#if(DEBUGMINMEM)
	  fprintf(fail_file,"got here3\n"); fflush(fail_file);
#endif

	  // check!
	  if(uii!=cpulinkptr->num){
	    dualfprintf(fail_file,"uii and num for this cpu not same, algorithm failure: uii=%d num=%d datatogetc0=%d\n",uii,cpulinkptr->num,datatogetc0);
	    myexit(10000);
	  }
#if(DEBUGMINMEM)
	  fprintf(fail_file,"got here4\n"); fflush(fail_file);	  
#endif
	  // jonio is the unsorted bit here starting at index=0 (for all cpus)
	  MPI_Isend(jonio,cpulinkptr->num,datatype,cpulinkptr->cpu,cpulinkptr->cpu,MPI_COMM_WORLD,&request0);
	  // have myid wait before continuing to make sure receive is complete
	  // alternative is to have many sends, but then need to desort all at once which isn't easy since cycling through link list one at a time only once (could store starter pointer, etc.)
	  MPI_Wait(&request0,&mpichstatus);
	  
	  datasent0+=cpulinkptr->num; // diagnostic

	  if(cpulinkptr->end) doset=0;
	  cpulinkptr=cpulinkptr->np;
#if(DEBUGMINMEM)
	  fprintf(fail_file,"got here5\n"); fflush(fail_file);	  
#endif

	}
	// check to see if total data is correct
#if(DEBUGMINMEM)
	fprintf(fail_file,"got here6\n"); fflush(fail_file);	  
#endif

	datasofar0+=datatoget0; // diagnostic
	if(datasofar0!=datasofarc0){
	  dualfprintf(fail_file,"cumulative data received via MPI and expected data is different: got=%d expected=%d\n",datasofar0,datasofarc0);
	  myexit(10000);
	}
	if(datasent0!=datatoget0){
	  dualfprintf(fail_file,"data sent doesn't match data read\n");
	  myexit(10000);
	}
	if((done0==0)&&(doing0==1)) dofull=0; // cpu==0 still needs to deal with more reads for it's own data
	if(cpulinkptr==NULL){
	  masterdone=1; // i.e. don't come back
	  dofull=0; // end of list
	  if(datasofar0!=totalsize[1]*totalsize[2]*numcolumns){
	    dualfprintf(fail_file,"read: total data written not equal to expected amount: wrote=%d expected=%d\n",datasofar0,totalsize[1]*totalsize[2]*numcolumns);
	    myexit(10000);
	  }
	}
	// otherwise continue
#if(DEBUGMINMEM)
	fprintf(fail_file,"got here7\n"); fflush(fail_file);	  
#endif

      }
#if(DEBUGMINMEM)
      fprintf(fail_file,"got here8\n"); fflush(fail_file);	  
#endif

    }
  }

#if(DEBUGMINMEM)
  fprintf(fail_file,"got here 6\n"); fflush(fail_file);
#endif

  // have myid wait before continuing so buffer can be released for more writing to
  if(dolocal) MPI_Wait(&request,&mpichstatus); // myid==0 handled specially since also master, and can't wait, even if master done first in this function, since master may not require cpu=0 at all for its current node set 

  //  logsfprintf("end mpiminmem_read\n");
#if(DEBUGMINMEM)
  fprintf(fail_file,"got here9\n"); fflush(fail_file);	  
#endif
#endif
}



void mpiio_seperate(int bintxt, int sorted, int stage,
		    int numcolumns, MPI_Datatype datatype,
		    FILE ** fpptr, void *jonio, void *writebuf)
{
#if(USEMPI)

  logsfprintf("mpiio begin seperate\n");

  if(truempicombinetype==MPICOMBINEROMIO){
    // doesn't use jonioptr or bintxt or sorted
    // headerbytesize no longer needed
    if(stage==STAGE1) mpiioromio_init_combine(READROMIO, READFILE,  0, "", numcolumns, datatype,&writebuf,writebuf);
    if(stage==STAGE2) mpiioromio_init_combine(READFREEROMIO, READFILE,  0, "", numcolumns, datatype,&writebuf,writebuf);
  }
  else{

    if(truempicombinetype==MPICOMBINESIMPLE){
      
      if(sorted){
	mpiios_seperate(bintxt, stage, datatype, numcolumns, fpptr, jonio, writebuf);
      }
      else{
	mpiiotu_seperate(stage, datatype, numcolumns, fpptr, writebuf);
      }
      
      
    }
    else if(truempicombinetype==MPICOMBINEMINMEM){
      // do nothing on STAGE1
      if(stage==STAGE2) mpiiomin_final(numcolumns,fpptr, jonio, writebuf);
    }
  }

#endif
  logsfprintf("mpiio end seperate\n");

}

// filename only needed for INITROMIO
// operationtype=INITROMIO or WRITECLOSEROMIO
// mpi io using ROMIO
// which=WRITEFILE or READFILE
void mpiioromio_init_combine(int operationtype, int which,  long headerbytesize, char *filename, int numcolumns,MPI_Datatype datatype, void **writebufptr,void *writebuf)
{
  int i;

  static int sizeofmemory;
  static int sizeofdatatype;

  static long double **writebuf16;
  static double **writebuf8;
  static float **writebuf4;
  static unsigned char **writebuf1;
  static int **writebuf4i;
  static long long int **writebuf8i;

  static int numfiles;
  static int romiocolumns;
  static char realfilename[MAXFILENAME];
#if(USEMPI&&USEROMIO)
  static MPI_Datatype newtype;
  static MPI_File fh;
  static MPI_Status status;
  static MPI_Request request;
#endif
  static int ndims, array_of_gsizes[4], array_of_distribs[4];
  static int order, len;
  static int array_of_dargs[4], array_of_psizes[4];
  static int bufcount, array_size;






  if(operationtype==INITROMIO){

    sizeofdatatype=getsizeofdatatype(datatype);

    logsfprintf("mpiioromio_init begin\n");


    //////////////////////////
    //
    // set dimensionality

    // see if splitting into individual columns
    if(docolsplit){
      numfiles=numcolumns;
      romiocolumns=1;
    }
    else{
      numfiles=1;
      romiocolumns=numcolumns;
    }


#if(USEMPI&&USEROMIO)
    if((COMPDIM==3)&&(romiocolumns>1)){
      //create the distributed array filetype
      ndims = 4;
      order = MPI_ORDER_C;
    
      array_of_gsizes[3] = romiocolumns;
      array_of_gsizes[2] = totalsize[1];
      array_of_gsizes[1] = totalsize[2];
      array_of_gsizes[0] = totalsize[3];
    
      array_of_distribs[3] = MPI_DISTRIBUTE_BLOCK;
      array_of_distribs[2] = MPI_DISTRIBUTE_BLOCK;
      array_of_distribs[1] = MPI_DISTRIBUTE_BLOCK;
      array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
    
      array_of_dargs[3] = MPI_DISTRIBUTE_DFLT_DARG;
      array_of_dargs[2] = MPI_DISTRIBUTE_DFLT_DARG;
      array_of_dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
      array_of_dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
    
      array_of_psizes[3]=1;
      array_of_psizes[2]=ncpux1;
      array_of_psizes[1]=ncpux2;
      array_of_psizes[0]=ncpux3;
    }
    else if((COMPDIM==2)&&(romiocolumns>1)){
      //create the distributed array filetype
      ndims = 3;
      order = MPI_ORDER_C;
    
      array_of_gsizes[2] = romiocolumns;
      array_of_gsizes[1] = totalsize[1];
      array_of_gsizes[0] = totalsize[2];
    
      array_of_distribs[2] = MPI_DISTRIBUTE_BLOCK;
      array_of_distribs[1] = MPI_DISTRIBUTE_BLOCK;
      array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
    
      array_of_dargs[2] = MPI_DISTRIBUTE_DFLT_DARG;
      array_of_dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
      array_of_dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
    
      array_of_psizes[2]=1;
      array_of_psizes[1]=ncpux1;
      array_of_psizes[0]=ncpux2;
    }
    else if((COMPDIM==3)&&(romiocolumns==1)){
      //create the distributed array filetype
      ndims=3;
      order = MPI_ORDER_C;
      
      array_of_gsizes[2] = totalsize[1];
      array_of_gsizes[1] = totalsize[2];
      array_of_gsizes[0] = totalsize[3];
      
      array_of_distribs[2] = MPI_DISTRIBUTE_BLOCK;
      array_of_distribs[1] = MPI_DISTRIBUTE_BLOCK;
      array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
      
      array_of_dargs[2] = MPI_DISTRIBUTE_DFLT_DARG;
      array_of_dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
      array_of_dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
      
      array_of_psizes[2]=ncpux1;
      array_of_psizes[1]=ncpux2;
      array_of_psizes[0]=ncpux3;
    }
    else  if((COMPDIM==2)&&(romiocolumns==1)){
      //create the distributed array filetype
      ndims=2;
      order = MPI_ORDER_C;
      
      array_of_gsizes[1] = totalsize[1];
      array_of_gsizes[0] = totalsize[2];
      
      array_of_distribs[1] = MPI_DISTRIBUTE_BLOCK;
      array_of_distribs[0] = MPI_DISTRIBUTE_BLOCK;
      
      array_of_dargs[1] = MPI_DISTRIBUTE_DFLT_DARG;
      array_of_dargs[0] = MPI_DISTRIBUTE_DFLT_DARG;
      
      array_of_psizes[1]=ncpux1;
      array_of_psizes[0]=ncpux2;
    }
    else{
      dualfprintf(fail_file,"Shouldn't reach to end of ROMIO selection\n");
      myexit(1);
    }
#endif
      
    // initialize filename
    // must be same filename as in dump_gen() when opening files here, so header is included at top of the dump
    if(docolsplit&&(numcolumns>1)){
      sprintf(realfilename,"%s-col%04d",filename,romiocoliter);
    }
    else strcpy(realfilename,filename);
     


    //////////////////////////
    //
    // initialize grid and file handler for ROMIO

#if(USEMPI&&USEROMIO)
    MPI_Type_create_darray(numprocs, myid, ndims, array_of_gsizes, 
			   array_of_distribs, array_of_dargs,
			   array_of_psizes, order, datatype, &newtype);
    MPI_Type_commit(&newtype);
    MPI_Type_size(newtype, &bufcount);
    sizeofmemory=bufcount; //includes type
    bufcount = bufcount/sizeofdatatype; // number of elements of type

    // setup file handler
    // setup to append, in case wanted binary header at front
    if(which==WRITEFILE){
      MPI_File_open(MPI_COMM_WORLD, realfilename, MPI_MODE_APPEND | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
    }
    else if(which==READFILE){
      MPI_File_open(MPI_COMM_WORLD, realfilename, MPI_MODE_RDONLY , MPI_INFO_NULL, &fh);
    }
    // this sets the distributed nature of the file writing
    MPI_File_set_view(fh, headerbytesize, datatype, newtype, "native", MPI_INFO_NULL);
    // all that needs to be done now is initialize and later fill writebuf with the data
#endif

    //////////////////////////
    //
    // allocate memory for writebuf
   

    
    if(datatype==MPI_UNSIGNED_CHAR){ writebuf1=(unsigned char **)writebufptr; }
    else if(datatype==MPI_FLOAT){  writebuf4=(float **)writebufptr; }
    else if(datatype==MPI_DOUBLE){  writebuf8=(double**)writebufptr; }
    else if(datatype==MPI_LONG_DOUBLE){  writebuf16=(long double**)writebufptr; }
    else if(datatype==MPI_INT){  writebuf4i=(int **)writebufptr; }
    else if(datatype==MPI_LONG_LONG_INT){  writebuf8i=(long long int **)writebufptr; }


    if(datatype==MPI_UNSIGNED_CHAR) *writebuf1=(unsigned char*)malloc(sizeofmemory);
    else if(datatype==MPI_FLOAT) *writebuf4=(float*)malloc(sizeofmemory);
    else if(datatype==MPI_DOUBLE) *writebuf8 =(double*)malloc(sizeofmemory);
    else if(datatype==MPI_LONG_DOUBLE) *writebuf16 =(long double*)malloc(sizeofmemory);
    else if(datatype==MPI_INT) *writebuf4i=(int*)malloc(sizeofmemory);
    else if(datatype==MPI_LONG_LONG_INT) *writebuf8i=(long long int*)malloc(sizeofmemory);
    if(
       ((datatype==MPI_UNSIGNED_CHAR)&&(*writebuf1 == NULL)) ||
       ((datatype==MPI_FLOAT)&&(*writebuf4 == NULL)) ||
       ((datatype==MPI_DOUBLE)&&(*writebuf8 == NULL)) ||
       ((datatype==MPI_LONG_DOUBLE)&&(*writebuf16 == NULL)) ||
       ((datatype==MPI_INT)&&(*writebuf4i == NULL)) ||
       ((datatype==MPI_LONG_LONG_INT)&&(*writebuf8i == NULL)) 
       ){
      dualfprintf(fail_file,"Can't initialize writebuf memory for mpiioromio\n");
      myexit(10000);
    }

    logsfprintf("mpiioromio_init end\n");

  }

#if(USEMPI&&USEROMIO)
  if(operationtype==READROMIO){
    logsfprintf("mpiioromio_seperate begin\n");
    MPI_File_read_all(fh, writebuf, bufcount, datatype, &status);
    MPI_File_close(&fh);
    logsfprintf("mpiioromio_seperate end\n");
  }
  else if(operationtype==READFREEROMIO){
    logsfprintf("mpiioromio_seperate_free begin\n");
    free(writebuf);
    MPI_Type_free(&newtype);
    logsfprintf("mpiioromio_seperate_free end\n");
  }
  else  if(operationtype==WRITECLOSEROMIO){
    logsfprintf("mpiioromio_combine begin\n");

    // now write the buffer:
    MPI_File_write_all(fh, writebuf, bufcount, datatype, &status);
    MPI_File_close(&fh);
    // free buffer and type
    free(writebuf);
    MPI_Type_free(&newtype);
    logsfprintf("mpiioromio_combine end\n");
  
  }
#endif

}


#define DEBUGSINIT 0


// initializes memory buffers for MPI combining for all MPI combine types(simple, minmem, romio)
// mpi io sorted or not
void mpiios_init(int bintxt, int sorted, FILE ** fp, int which, int headerbytesize, char *filename, int numcolumns,
		MPI_Datatype datatype, void **jonioptr, void **writebufptr)
{
  int i;

  int sizeofmemory;
  int sizeofdatatype;

  long double **jonio16;
  double **jonio8;
  float **jonio4;
  unsigned char **jonio1;
  int **jonio4i;
  long long int **jonio8i;

  long double **writebuf16;
  double **writebuf8;
  float **writebuf4;
  unsigned char **writebuf1;
  int **writebuf4i;
  long long int **writebuf8i;

  int numfiles;
  char newfilename[200];

#if(USEMPI)

  sizeofdatatype=getsizeofdatatype(datatype);

  if(datatype==MPI_UNSIGNED_CHAR){ jonio1=(unsigned char **)jonioptr; writebuf1=(unsigned char **)writebufptr; }
  else if(datatype==MPI_FLOAT){ jonio4=(float **)jonioptr; writebuf4=(float **)writebufptr; }
  else if(datatype==MPI_DOUBLE){ jonio8=(double**)jonioptr; writebuf8=(double**)writebufptr; }
  else if(datatype==MPI_LONG_DOUBLE){ jonio16=(long double**)jonioptr; writebuf16=(long double**)writebufptr; }
  else if(datatype==MPI_INT){ jonio4i=(int **)jonioptr; writebuf4i=(int **)writebufptr; }
  else if(datatype==MPI_LONG_LONG_INT){ jonio8i=(long long int **)jonioptr; writebuf8i=(long long int **)writebufptr; }
  
  logsfprintf("mpiios_init begin\n");


  /////////////////////
  //
  //  open files for non-ROMIO writing
  //
  /////////////////////
  
  if(myid==0){		// total on CPU=0, always, since mpicombine=1
    if(docolsplit){
      numfiles=numcolumns;
    }
    else numfiles=1;

    for(i=0;i<numfiles;i++){
      if(docolsplit&&(numfiles>1)){
	sprintf(newfilename,"%s-col%04d",filename,i); // must be same name as in dump_gen()
      }
      else{
	sprintf(newfilename,"%s",filename);
      }
      if (which == WRITEFILE){
	if(bintxt==BINARYOUTPUT)         fp[i] = fopen(newfilename, "a");
	else if(bintxt==TEXTOUTPUT)      fp[i] = fopen(newfilename, "at");
	// file pointer set correctly upon append
      }
      else if (which == READFILE){
	if(bintxt==BINARYOUTPUT)     fp[i] = fopen(newfilename, "rb");
	else if(bintxt==TEXTOUTPUT)  fp[i] = fopen(newfilename, "rt");
	// must set file pointer to after header
	fseek(fp[i],headerbytesize,SEEK_SET);
      }
      if (fp[i] == NULL) {
	fprintf(fail_file, "error opening file: %s\n", newfilename);
	myexit(2);
      }
    }
  }
  
  if ( (sorted==SORTED)&&(myid == 0) ){		// total on CPU=0 is stored in memory somehow for sorting before output

    /////////////
    //
    // determine the memory needed for the mpicombinetype for cpu=0
    //
    ///////////////

    truempicombinetype=mpicombinetype; // initial action
    if(truempicombinetype==MPICOMBINESIMPLE){
      sizeofmemory = sizeofdatatype * totalsize[1] * totalsize[2] * numcolumns;
    }
    else if(truempicombinetype==MPICOMBINEMINMEM){
      // 2 needed since need to read sequentially, then order it into the other buffer for writing
      // check this calculation against setuplinklist()'s cpulist0 array size!
      sizeofmemory = sizeofdatatype * (int)ceil(N1 * N2 * NUMBUFFERS/numcolumns)*numcolumns * 2 ; // default
      if(sizeofmemory>sizeofdatatype * totalsize[1] * totalsize[2] * numcolumns){
	sizeofmemory = sizeofdatatype * totalsize[1] * totalsize[2] * numcolumns;
	truempicombinetype=MPICOMBINESIMPLE; // then switch to simple method
      }
      // need memory to be at least larger than number of columns (X2 for 2 buffers)
      // don't want to work with chunks smaller than # of columns, and all chunks should come in # of column chunks times some integer multiple
      if(sizeofmemory<sizeofdatatype*numcolumns*2) sizeofmemory=sizeofdatatype*numcolumns*2;
      if(sizeofmemory<2*numcolumns){
	dualfprintf(fail_file,"problem, sizeofmemory=%d < %d=2*numcolumns\n",sizeofmemory,2*numcolumns);
	myexit(10000);
      }

    }
    joniosize=sizeofmemory/sizeofdatatype;
    if(datatype==MPI_UNSIGNED_CHAR) *jonio1=(unsigned char*)malloc(sizeofmemory);
    else if(datatype==MPI_FLOAT) *jonio4=(float*)malloc(sizeofmemory);
    else if(datatype==MPI_DOUBLE) *jonio8 =(double*)malloc(sizeofmemory);
    else if(datatype==MPI_LONG_DOUBLE) *jonio16 =(long double*)malloc(sizeofmemory);
    else if(datatype==MPI_INT) *jonio4i=(int*)malloc(sizeofmemory);
    else if(datatype==MPI_LONG_LONG_INT) *jonio8i=(long long int*)malloc(sizeofmemory);
    if(
       (datatype==MPI_UNSIGNED_CHAR)&&(*jonio1 == NULL) ||
       (datatype==MPI_FLOAT)&&(*jonio4 == NULL) ||
       (datatype==MPI_DOUBLE)&&(*jonio8 == NULL) ||
       (datatype==MPI_LONG_DOUBLE)&&(*jonio16 == NULL) ||
       (datatype==MPI_INT)&&(*jonio4i == NULL) ||
       (datatype==MPI_LONG_LONG_INT)&&(*jonio8i == NULL) 
       ){
      dualfprintf(fail_file, "Can't initialize jonio memory for mpiios_init\n");
      myexit(10000);
    }
#if(DEBUGSINIT)
    /*
    for(i=0;i<sizeofmemory/datatype;i++){
    if(datatype==sizeof(long double)) (*jonio16)[i]=-10000.0000;
    if(datatype==sizeof(double)) (*jonio8)[i]=-10000.0000;
    if(datatype==sizeof(float)) (*jonio4)[i]=-10000.0000;
    if(datatype==sizeof(unsigned char)) (*jonio1)[i]=100;
    }
    dualfprintf(fail_file,"got here1\n");
    */
#endif
  }


  //////////////////////
  //
  // need to tell all CPUS the status of changed global vars
  //
  ///////////////////////

  MPI_Bcast(&truempicombinetype,1,MPI_INT,0,MPI_COMM_WORLD);


  ///////////////////////////////
  //
  // open per CPU writebuf
  //
  /////////////////////////////////

  if(truempicombinetype==MPICOMBINESIMPLE){
    sizeofmemory = sizeofdatatype * N1 * N2 * numcolumns ;
  }
  else if(truempicombinetype==MPICOMBINEMINMEM){
    // maximum cpu=0 could require under any case
    sizeofmemory = sizeofdatatype * (int)ceil(N1 * N2 * NUMBUFFERS/numcolumns)*numcolumns ; // default
    if(sizeofmemory>sizeofdatatype*N1*N2*numcolumns) sizeofmemory=sizeofdatatype*N1*N2*numcolumns; // can't ask for more!
    if(sizeofmemory<sizeofdatatype*numcolumns) sizeofmemory=sizeofdatatype*numcolumns; // minimum, chunk at minimum of number of columns
    if(sizeofmemory<numcolumns){
      dualfprintf(fail_file,"problem, sizeofmemory=%d < %d=numcolumns\n",sizeofmemory,numcolumns);
      myexit(10000);
    }
  }
  writebufsize=sizeofmemory/sizeofdatatype; // never used currently
  if(datatype==MPI_UNSIGNED_CHAR) *writebuf1=(unsigned char*)malloc(sizeofmemory);
  else if(datatype==MPI_FLOAT) *writebuf4=(float*)malloc(sizeofmemory);
  else if(datatype==MPI_DOUBLE) *writebuf8 =(double*)malloc(sizeofmemory);
  else if(datatype==MPI_LONG_DOUBLE) *writebuf16 =(long double*)malloc(sizeofmemory);
  else if(datatype==MPI_INT) *writebuf4i=(int*)malloc(sizeofmemory);
  else if(datatype==MPI_LONG_LONG_INT) *writebuf8i=(long long int*)malloc(sizeofmemory);
  if(
     (datatype==MPI_UNSIGNED_CHAR)&&(*writebuf1 == NULL) ||
     (datatype==MPI_FLOAT)&&(*writebuf4 == NULL) ||
     (datatype==MPI_DOUBLE)&&(*writebuf8 == NULL) ||
     (datatype==MPI_LONG_DOUBLE)&&(*writebuf16 == NULL) ||
     (datatype==MPI_INT)&&(*writebuf4i == NULL) ||
     (datatype==MPI_LONG_LONG_INT)&&(*writebuf8i == NULL) 
     ){
    dualfprintf(fail_file, "Can't initialize writebuf memory for mpiios_init: datatype=%d sizeofmemory=%d\n",datatype,sizeofmemory);
    myexit(10000);
  }
#if(DEBUGSINIT)
  /*
  for(i=0;i<sizeofmemory/sizeofdatatype;i++){
  if(datatype==sizeof(long double)) (*writebuf16)[i]=-10000.0000;
    if(datatype==sizeof(double)) (*writebuf8)[i]=-10000.0000;
    if(datatype==sizeof(float)) (*writebuf4)[i]=-10000.0000;
    if(datatype==sizeof(unsigned char)) (*writebuf1)[i]=100;
  }
  dualfprintf(fail_file,"got here2\n");
  */
#endif
#endif
  logsfprintf("mpiios_init end\n");

}

void mpiiomin_final(int numcolumns,FILE **fp, void *jonio, void *writebuf)
{
  int i;
  int numfiles;

  free(writebuf);		// writebuf used by each CPU
  if(myid==0){
    free(jonio);		// used by CPU=0
    if(docolsplit){
      numfiles=numcolumns;
    }
    else numfiles=1;
    for(i=0;i<numfiles;i++){
      fclose(fp[i]);
      fp[i] = NULL;
    }
  }
}

#if(USEMPI)


void mpiios_combine(int bintxt, MPI_Datatype datatype, int numcolumns,
		   FILE ** fp, void *jonio, void *writebuf)
{
  int i, j, k, l, col, mapvaluejonio, mapvaluetempbuf;
#if(USEMPI)
  MPI_Request rrequest;
  MPI_Request srequest;
#endif
  int othercpupos[COMPDIM + 1];

  unsigned char *jonio1;
  float *jonio4;
  double *jonio8;
  long double *jonio16;
  int *jonio4i;
  long long int *jonio8i;

  unsigned char *writebuf1;
  float *writebuf4;
  double *writebuf8;
  long double *writebuf16;
  int *writebuf4i;
  long long int *writebuf8i;
  
  int numfiles;
  int sizeofdatatype;
  
  logsfprintf("mpiios begin combine\n");

  sizeofdatatype=getsizeofdatatype(datatype);

  if (datatype == MPI_UNSIGNED_CHAR) jonio1 = (unsigned char *) jonio;
  else if (datatype == MPI_FLOAT) jonio4 = (float *) jonio;
  else if (datatype == MPI_DOUBLE) jonio8 = (double *) jonio;
  else if (datatype == MPI_LONG_DOUBLE) jonio16 = (long double *) jonio;
  else if (datatype == MPI_INT) jonio4i = (int *) jonio;
  else if (datatype == MPI_LONG_LONG_INT) jonio8i = (long long int *) jonio;
  
  if (datatype == MPI_UNSIGNED_CHAR) writebuf1 = (unsigned char *) writebuf;
  else if (datatype == MPI_FLOAT) writebuf4 = (float *) writebuf;
  else if (datatype == MPI_DOUBLE) writebuf8 = (double *) writebuf;
  else if (datatype == MPI_LONG_DOUBLE) writebuf16 = (long double *) writebuf;
  else if (datatype == MPI_INT) writebuf4i = (int *) writebuf;
  else if (datatype == MPI_LONG_LONG_INT) writebuf8i = (long long int *) writebuf;


#if(USEMPI)
  // no need for tempbuf, works since first write to jonio is CPU=0's writebuf
  if(myid!=0) MPI_Isend(writebuf, N1 * N2 * numcolumns, datatype, 0, myid,
			MPI_COMM_WORLD, &srequest);
  if (myid == 0) {
    for (l = 0; l < numprocs; l++) {
      if(l!=0){
	MPI_Irecv(writebuf, N1 * N2 * numcolumns, datatype, l, l, MPI_COMM_WORLD, &rrequest);
	MPI_Wait(&rrequest, &mpichstatus);
      }
      othercpupos[1] = l % ncpux1;
      othercpupos[2] = (int) ((l % (numprocs)) / ncpux1);
      // now fill jonio with proper sequence (i.e. tiled mapping)
      for (j = 0; j < N2; j++)
	for (i = 0; i < N1; i++)
	  for (col = 0; col < numcolumns; col++) {
	    mapvaluejonio =
	      + ncpux1 * N1 * numcolumns * (j + othercpupos[2] * N2)
	      + numcolumns * (i + othercpupos[1] * N1)
	      + col;
	    mapvaluetempbuf =
	      + j * N1 * numcolumns
	      +  i * numcolumns + col;
	    
	    if (datatype == MPI_UNSIGNED_CHAR) jonio1[mapvaluejonio] = writebuf1[mapvaluetempbuf];
	    else if (datatype == MPI_FLOAT) jonio4[mapvaluejonio] = writebuf4[mapvaluetempbuf];
	    else if (datatype == MPI_DOUBLE) jonio8[mapvaluejonio] = writebuf8[mapvaluetempbuf];
	    else if (datatype == MPI_LONG_DOUBLE) jonio16[mapvaluejonio] = writebuf16[mapvaluetempbuf];
	    else if (datatype == MPI_INT) jonio4i[mapvaluejonio] = writebuf4i[mapvaluetempbuf];
	    else if (datatype == MPI_LONG_LONG_INT) jonio8i[mapvaluejonio] = writebuf8i[mapvaluetempbuf];
	  }
    }
  }
  if(myid!=0) MPI_Wait(&srequest, &mpichstatus);
  free(writebuf);		// writebuf used by each CPU

  if (myid == 0) {
    // now write out collected data using CPU=0
    if(docolsplit){
      numfiles=numcolumns;
    }
    else numfiles=1;

    if(bintxt==BINARYOUTPUT){
      if(numfiles==1) fwrite(jonio, sizeofdatatype,totalsize[1] * totalsize[2] * numcolumns, *fp);
      else{
	for(i=0;i<totalsize[1]*totalsize[2]*numcolumns;i++){
	  if (datatype == MPI_UNSIGNED_CHAR) fwrite(&jonio1[i], sizeofdatatype,1, fp[i%numfiles]);
	  else if (datatype == MPI_FLOAT) fwrite(&jonio4[i], sizeofdatatype,1, fp[i%numfiles]);
	  else if (datatype == MPI_DOUBLE) fwrite(&jonio8[i], sizeofdatatype,1, fp[i%numfiles]);
	  else if (datatype == MPI_LONG_DOUBLE) fwrite(&jonio16[i], sizeofdatatype,1, fp[i%numfiles]);
	  else if (datatype == MPI_INT) fwrite(&jonio4i[i], sizeofdatatype,1, fp[i%numfiles]);
	  else if (datatype == MPI_LONG_LONG_INT) fwrite(&jonio8i[i], sizeofdatatype,1, fp[i%numfiles]);
	}
      }
    }
    else if(bintxt==TEXTOUTPUT){ // properly ordered, so just dump it
      for(i=0;i<totalsize[1]*totalsize[2]*numcolumns;i++){
	if (datatype == MPI_UNSIGNED_CHAR) fprintf(fp[i%numfiles],"%c",jonio1[i]);
	else if (datatype == MPI_FLOAT) fprintf(fp[i%numfiles],"%15.7g",jonio4[i]);
	else if (datatype == MPI_DOUBLE) fprintf(fp[i%numfiles],"%21.15g",jonio8[i]);
	else if (datatype == MPI_LONG_DOUBLE) fprintf(fp[i%numfiles],"%31.25Lg",jonio16[i]);
	else if (datatype == MPI_INT) fprintf(fp[i%numfiles],"%d",jonio4i[i]);
	else if (datatype == MPI_LONG_LONG_INT) fprintf(fp[i%numfiles],"%lld",jonio8i[i]);
	if(numfiles==1){
	  if((i+1)%numcolumns) fprintf(*fp," ");
	  else fprintf(*fp,"\n");
	}
	else  fprintf(fp[i%numfiles],"\n");
      }
    }
    free(jonio);		// used by CPU=0
    for(i=0;i<numfiles;i++){
      fclose(fp[i]);
      fp[i] = NULL;
    }
  }
#endif

  logsfprintf("mpiios end combine\n");

}

void mpiios_seperate(int bintxt, int stage, MPI_Datatype datatype, int numcolumns,
		    FILE ** fp, void *jonio,
		    void *writebuf)
{
  int i, j, k, l, col, mapvaluejonio, mapvaluetempbuf;
#if(USEMPI)
  MPI_Request rrequest;
  MPI_Request srequest;
#endif
  int othercpupos[COMPDIM + 1];

  int sizeofdatatype;

  unsigned char *jonio1;
  float *jonio4;
  double *jonio8;
  long double *jonio16;
  int *jonio4i;
  long long int *jonio8i;

  unsigned char *writebuf1;
  float *writebuf4;
  double *writebuf8;
  long double *writebuf16;
  int *writebuf4i;
  long long int *writebuf8i;


  int numfiles;


  logsfprintf("mpiios begin seperate\n");

  sizeofdatatype=getsizeofdatatype(datatype);

  if (datatype == MPI_UNSIGNED_CHAR) jonio1 = (unsigned char *) jonio;
  else if (datatype == MPI_FLOAT) jonio4 = (float *) jonio;
  else if (datatype == MPI_DOUBLE) jonio8 = (double *) jonio;
  else if (datatype == MPI_LONG_DOUBLE) jonio16 = (long double *) jonio;
  else if (datatype == MPI_INT) jonio4i = (int *) jonio;
  else if (datatype == MPI_LONG_LONG_INT) jonio8i = (long long int *) jonio;


  if (datatype == MPI_UNSIGNED_CHAR) writebuf1 = (unsigned char *) writebuf;
  else if (datatype == MPI_FLOAT) writebuf4 = (float *) writebuf;
  else if (datatype == MPI_DOUBLE) writebuf8 = (double *) writebuf;
  else if (datatype == MPI_LONG_DOUBLE) writebuf16 = (long double *) writebuf;
  else if (datatype == MPI_INT) writebuf4i = (int *) writebuf;
  else if (datatype == MPI_LONG_LONG_INT) writebuf8i = (long long int *) writebuf;


#if(USEMPI)

  if (stage == STAGE1) {
    if (myid == 0) {
      if(docolsplit){
	numfiles=numcolumns;
      }
      else numfiles=1;
      
      if(bintxt==BINARYOUTPUT){
	// first let cpu=0 read data
	if(numfiles==1) fread(jonio, sizeofdatatype,totalsize[1] * totalsize[2] * numcolumns, *fp);
	else{
	  for(i=0;i<totalsize[1]*totalsize[2]*numcolumns;i++){
	    if (datatype == MPI_UNSIGNED_CHAR) fread(&jonio1[i], sizeofdatatype,1, fp[i%numfiles]);
	    else if (datatype == MPI_FLOAT) fread(&jonio4[i], sizeofdatatype,1, fp[i%numfiles]);
	    else if (datatype == MPI_DOUBLE)  fread(&jonio8[i], sizeofdatatype,1, fp[i%numfiles]);
	    else if (datatype == MPI_LONG_DOUBLE) fread(&jonio16[i], sizeofdatatype,1, fp[i%numfiles]);
	    else if (datatype == MPI_INT) fread(&jonio4i[i], sizeofdatatype,1, fp[i%numfiles]);
	    else if (datatype == MPI_LONG_LONG_INT) fread(&jonio8i[i], sizeofdatatype,1, fp[i%numfiles]);
	  }
	}
      }
      else if(bintxt==TEXTOUTPUT){ // properly ordered, so just dump it
	for(i=0;i<totalsize[1]*totalsize[2]*numcolumns;i++){
	  if (datatype == MPI_UNSIGNED_CHAR) fscanf(fp[i%numfiles],"%uc",&jonio1[i]);
	  else if (datatype == MPI_FLOAT) fscanf(fp[i%numfiles],"%f",&jonio4[i]);
	  else if (datatype == MPI_DOUBLE) fscanf(fp[i%numfiles],"%lf",&jonio8[i]);
	  else if (datatype == MPI_LONG_DOUBLE) fscanf(fp[i%numfiles],"%Lf",&jonio16[i]);
	  else if (datatype == MPI_INT) fscanf(fp[i%numfiles],"%d",&jonio4i[i]);
	  else if (datatype == MPI_LONG_LONG_INT) fscanf(fp[i%numfiles],"%lld",&jonio8i[i]);
	}
      }
    }
    // writebuf is CPU=0's tempbuf for each CPU, including CPU=0, which is done last
    if (myid == 0) {
      for (l = numprocs-1 ; l >=0; l--) {
	othercpupos[1] = l % ncpux1;
	othercpupos[2] = (int) ((l % (numprocs)) / ncpux1);
	// now unfill jonio with proper sequence (i.e. tiled mapping)
	for (j = 0; j < N2; j++) {
	  for (i = 0; i < N1; i++) {
	    for (col = 0; col < numcolumns; col++) {
	      mapvaluejonio =
		+ ncpux1 * N1 * numcolumns * (j +
					      othercpupos[2] * N2)
		+ numcolumns * (i + othercpupos[1] * N1)
		+ col;
	      mapvaluetempbuf =
		+ j * N1 * numcolumns
		+ i * numcolumns + col;
	      // debug check
	      if((mapvaluetempbuf<0)||(mapvaluetempbuf>=N1*N2*numcolumns)){
		dualfprintf(fail_file,"mapvaluetempbuf out of range: %d\n",mapvaluetempbuf);
		myexit(10000);
	      }
	      if((mapvaluejonio<0)||(mapvaluejonio>=totalsize[1]*totalsize[2]*numcolumns)){
		dualfprintf(fail_file,"mapvaluejonio out of range: %d\n",mapvaluejonio);
		myexit(10000);
	      }
	      if (datatype == MPI_UNSIGNED_CHAR) writebuf1[mapvaluetempbuf] = jonio1[mapvaluejonio];
	      if (datatype == MPI_FLOAT) writebuf4[mapvaluetempbuf] = jonio4[mapvaluejonio];
	      if (datatype == MPI_DOUBLE) writebuf8[mapvaluetempbuf] = jonio8[mapvaluejonio];
	      if (datatype == MPI_LONG_DOUBLE) writebuf16[mapvaluetempbuf] = jonio16[mapvaluejonio];
	      if (datatype == MPI_INT) writebuf4i[mapvaluetempbuf] = jonio4i[mapvaluejonio];
	      if (datatype == MPI_LONG_LONG_INT) writebuf8i[mapvaluetempbuf] = jonio8i[mapvaluejonio];
	    }
	  }
	}
	if(l!=0){
	  
	  MPI_Isend(writebuf, N1 * N2 * numcolumns, datatype, l, l,
	      MPI_COMM_WORLD, &srequest);
	  MPI_Wait(&srequest, &mpichstatus);
	}
      }
      free(jonio); // done with jonio after loop
    }
    else{
      // chosen CPU to receive data from CPU=0
      
            MPI_Irecv(writebuf, N1 * N2 * numcolumns, datatype, 0, myid,
      		MPI_COMM_WORLD, &rrequest);
      MPI_Wait(&rrequest, &mpichstatus);	// writebuf used until      
    }
  } else if (stage == STAGE2) {
    free(writebuf);
    if (myid == 0) {
      fclose(*fp);
      *fp = NULL;
    }
  }
#endif

 logsfprintf("mpiios end seperate\n");



}


void mpiiotu_combine(MPI_Datatype datatype, int numcolumns, FILE ** fp, void *writebuf)
{
  int i, j, k, l, col;
#if(USEMPI)
  MPI_Request rrequest;
  MPI_Request srequest;
#endif
  int othercpupos[COMPDIM + 1];

  int sizeofdatatype;

  unsigned char *writebuf1;
  float *writebuf4;
  double *writebuf8;
  long double *writebuf16;
  int *writebuf4i;
  long long int *writebuf8i;

  int numfiles;

  sizeofdatatype=getsizeofdatatype(datatype);

  if (datatype == MPI_UNSIGNED_CHAR) writebuf1 = (unsigned char *) writebuf;
  else if (datatype == MPI_FLOAT) writebuf4 = (float *) writebuf;
  else if (datatype == MPI_DOUBLE) writebuf8 = (double *) writebuf;
  else if (datatype == MPI_LONG_DOUBLE) writebuf16 = (long double *) writebuf;
  else if (datatype == MPI_INT) writebuf4i = (int *) writebuf;
  else if (datatype == MPI_LONG_LONG_INT) writebuf8i = (long long int *) writebuf;
  

#if(USEMPI)
  if(myid!=0) MPI_Isend(writebuf, N1 * N2 * numcolumns, datatype, 0, myid,
			MPI_COMM_WORLD, &srequest);
  if (myid == 0) {   
    if(docolsplit){
      numfiles=numcolumns;
    }
    else numfiles=1;

    // done in forward order, no need to use tempbuf since CPU=0's writebuf is first out
    for (l = 0; l <numprocs; l++) {
      if(l!=0){
	MPI_Irecv(writebuf, N1 * N2 * numcolumns, datatype, l, l,
		  MPI_COMM_WORLD, &rrequest);
	MPI_Wait(&rrequest, &mpichstatus);
      }
      // now write writebuf
      DUMPLOOP(0,N1-1,0,N2-1){
	for(col=0;col<numcolumns;col++){
	  if(datatype== MPI_UNSIGNED_CHAR){
	    fprintf(fp[(col+numcolumns*(i+N1*j))%numfiles],"%c ",writebuf1[col+numcolumns*(i+N1*j)]);
	  }
	  else if(datatype==MPI_FLOAT){
	    fprintf(fp[(col+numcolumns*(i+N1*j))%numfiles],"%15.7g ",writebuf4[col+numcolumns*(i+N1*j)]);
	  }
	  else if(datatype==MPI_DOUBLE){
	    fprintf(fp[(col+numcolumns*(i+N1*j))%numfiles],"%21.15g ",writebuf8[col+numcolumns*(i+N1*j)]);
	  }
	  else if(datatype==MPI_LONG_DOUBLE){
	    fprintf(fp[(col+numcolumns*(i+N1*j))%numfiles],"%31.25Lg ",writebuf16[col+numcolumns*(i+N1*j)]);
	  }
	  else if(datatype==MPI_INT){
	    fprintf(fp[(col+numcolumns*(i+N1*j))%numfiles],"%d ",writebuf4i[col+numcolumns*(i+N1*j)]);
	  }
	  else if(datatype==MPI_LONG_LONG_INT){
	    fprintf(fp[(col+numcolumns*(i+N1*j))%numfiles],"%lld ",writebuf8i[col+numcolumns*(i+N1*j)]);
	  }
	  if(numfiles>1) fprintf(fp[(col+numcolumns*(i+N1*j))%numfiles],"\n");
	}
	if(numfiles==1) fprintf(*fp,"\n");
      }
    }    
  }
  if(myid!=0) MPI_Wait(&srequest, &mpichstatus);
  free(writebuf);		// writebuf used by each CPU

  if (myid == 0) {
    for(i=0;i<numfiles;i++){
      fclose(fp[i]);
      fp[i] = NULL;
    }
  }
#endif

}

// fill writebuf with each cpu's data set,using CPU=0 to process the file
void mpiiotu_seperate(int stage, MPI_Datatype datatype, int numcolumns,
		      FILE ** fp,void *writebuf)
{
  int i, j, k, l, col;
#if(USEMPI)
  MPI_Request rrequest;
  MPI_Request srequest;
#endif
  int othercpupos[COMPDIM + 1];

  void *tempbuf;
  void *sendbuf;

  int numfiles;

  int sizeofdatatype;

  unsigned char *writebuf1;
  float *writebuf4;
  double *writebuf8;
  long double *writebuf16;
  int *writebuf4i;
  long long int *writebuf8i;

  unsigned char *tempbuf1;
  float *tempbuf4;
  double *tempbuf8;
  long double *tempbuf16;
  int *tempbuf4i;
  long long int *tempbuf8i;

  unsigned char *sendbuf1;
  float *sendbuf4;
  double *sendbuf8;
  long double *sendbuf16;
  int *sendbuf4i;
  long long int *sendbuf8i;

  sizeofdatatype=getsizeofdatatype(datatype);

  if (datatype == MPI_UNSIGNED_CHAR) writebuf1 = (unsigned char *) writebuf;
  else if (datatype == MPI_FLOAT) writebuf4 = (float *) writebuf;
  else if (datatype == MPI_DOUBLE) writebuf8 = (double *) writebuf;
  else if (datatype == MPI_LONG_DOUBLE) writebuf16 = (long double *) writebuf;
  else if (datatype == MPI_INT) writebuf4i = (int *) writebuf;
  else if (datatype == MPI_LONG_LONG_INT) writebuf8i = (long long int *) writebuf;

  if(myid==0){
    if((tempbuf=malloc(sizeofdatatype*N1*N2*numcolumns))==NULL){
      dualfprintf(fail_file,"Can't open tempbuf in gammieio_sep\n");
      myexit(10000);
    }

    if (datatype == MPI_UNSIGNED_CHAR) tempbuf1 = (unsigned char *) tempbuf;
    else if (datatype == MPI_FLOAT) tempbuf4 = (float *) tempbuf;
    else if (datatype == MPI_DOUBLE) tempbuf8 = (double *) tempbuf;
    else if (datatype == MPI_LONG_DOUBLE) tempbuf16 = (long double *) tempbuf;
    else if (datatype == MPI_INT) tempbuf4i = (int *) tempbuf;
    else if (datatype == MPI_LONG_LONG_INT) tempbuf8i = (long long int *) tempbuf;

  }

#if(USEMPI)

  if (stage == 1) {
    if (myid == 0) {
      if(docolsplit){
	numfiles=numcolumns;
      }
      else numfiles=1;

      for (l = 0; l < numprocs; l++) {
	if(l==0){
	  sendbuf=writebuf;
	  sendbuf1=writebuf1;
	  sendbuf4=writebuf4;
	  sendbuf8=writebuf8;
	  sendbuf4i=writebuf4i;
	  sendbuf8i=writebuf8i;
	}
	else{
	  sendbuf=tempbuf;
	  sendbuf1=tempbuf1;
	  sendbuf4=tempbuf4;
	  sendbuf8=tempbuf8;
	  sendbuf4i=tempbuf4i;
	  sendbuf8i=tempbuf8i;
	}

	DUMPLOOP(0,N1-1,0,N2-1) for (col = 0; col < numcolumns; col++) {
	  if(datatype==MPI_UNSIGNED_CHAR) fscanf(fp[(col+numcolumns*(i+N2*j))%numfiles],"%uc",&sendbuf1[col+numcolumns*(i+N2*j)]);
	  else if(datatype==MPI_FLOAT) fscanf(fp[(col+numcolumns*(i+N2*j))%numfiles],"%f",&sendbuf4[col+numcolumns*(i+N2*j)]);
	  else if(datatype==MPI_DOUBLE) fscanf(fp[(col+numcolumns*(i+N2*j))%numfiles],"%lf",&sendbuf8[col+numcolumns*(i+N2*j)]);
	  else if(datatype==MPI_LONG_DOUBLE) fscanf(fp[(col+numcolumns*(i+N2*j))%numfiles],"%Lf",&sendbuf16[col+numcolumns*(i+N2*j)]);
	  else if(datatype==MPI_INT) fscanf(fp[(col+numcolumns*(i+N2*j))%numfiles],"%d",&sendbuf4i[col+numcolumns*(i+N2*j)]);
	  else if(datatype==MPI_LONG_LONG_INT) fscanf(fp[(col+numcolumns*(i+N2*j))%numfiles],"%lld",&sendbuf8i[col+numcolumns*(i+N2*j)]);
	}
	if(l!=0){
	  MPI_Isend(sendbuf, N1 * N2 * numcolumns, datatype, l, l,
		    MPI_COMM_WORLD, &srequest);
	  // have to wait before filling sendbuf buffer again for next CPU
	  MPI_Wait(&srequest, &mpichstatus);
	}
      }
      free(tempbuf);
    }
    else{
      MPI_Irecv(writebuf, N1 * N2 * numcolumns, datatype, 0, myid, MPI_COMM_WORLD, &rrequest);
      MPI_Wait(&rrequest, &mpichstatus);	// writebuf used until
    }
  } else if (stage == 2) {
    free(writebuf);
    if (myid == 0) {
      for(i=0;i<numfiles;i++){
	fclose(fp[i]);
	fp[i] = NULL;
      }
    }
  }
#endif



}




#endif


int getsizeofdatatype(MPI_Datatype datatype)
{
  int sizeofdatatype;

  if (datatype == MPI_UNSIGNED_CHAR) sizeofdatatype=sizeof(unsigned char);
  else if (datatype == MPI_FLOAT) sizeofdatatype=sizeof(float);
  else if (datatype == MPI_DOUBLE) sizeofdatatype=sizeof(double);
  else if (datatype == MPI_LONG_DOUBLE) sizeofdatatype=sizeof(long double);
  else if (datatype == MPI_INT) sizeofdatatype=sizeof(int);
  else if (datatype == MPI_LONG_LONG_INT) sizeofdatatype=sizeof(long long int);
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(1);
  }

  return(sizeofdatatype);
}


// a simple max, assumes local cpu max already found
// sends results back to all cpus
void mpimax(SFTYPE*maxptr)
{
  SFTYPE send;
  
#if(USEMPI)
  send = *maxptr;
  MPI_Allreduce(&send, maxptr, 1, MPI_SFTYPE, MPI_MAX,MPI_COMM_WORLD);
#endif
}

// a simple max, assumes local cpu max already found
// sends results back to all cpus
void mpiisum(int*sumptr)
{
  int send;
  
#if(USEMPI)
  send = *sumptr;
  MPI_Allreduce(&send, sumptr, 1, MPI_INT, MPI_SUM,MPI_COMM_WORLD);
#endif
}
// as above, but only myid=recvid gets result
void mpiisum0(int*sumptr, int recvid)
{
  int send;
  
#if(USEMPI)
  send = *sumptr;
  MPI_Reduce(&send, sumptr, 1, MPI_INT, MPI_SUM,recvid,MPI_COMM_WORLD);
#endif
}

// as above, but only myid=recvid gets result
void mpildsum0(long int*sumptr, int recvid)
{
  long int send;
  
#if(USEMPI)
  send = *sumptr;
  MPI_Reduce(&send, sumptr, 1, MPI_LONG, MPI_SUM,recvid,MPI_COMM_WORLD);
#endif
}

// a simple max, assumes local cpu max already found
// sends results back to all cpus
void mpimin(SFTYPE*minptr)
{
  SFTYPE send;
  
#if(USEMPI)
  send = *minptr;
  MPI_Allreduce(&send, minptr, 1, MPI_SFTYPE, MPI_MIN,MPI_COMM_WORLD);
#endif
}

// a simple max, assumes local cpu max already found
// sends results back to all cpus
void mpifmin(FTYPE*minptr)
{
  FTYPE send;
  
#if(USEMPI)
  send = *minptr;
  MPI_Allreduce(&send, minptr, 1, MPI_FTYPE, MPI_MIN,MPI_COMM_WORLD);
#endif
}


void prminmaxsum(FTYPE p[][N2M][NPR], int start,int nmemb, FTYPE *maxptr, FTYPE*minptr,FTYPE*sumptr)
{
  int i,j,k;
  FTYPE maxsend,minsend,sumsend;
  int domin,domax,dosum;

  if(maxptr==NULL) domax=0; else domax=1;
  if(minptr==NULL) domin=0; else domin=1;
  if(sumptr==NULL) dosum=0; else dosum=1;
  
  for(k=start;k<start+nmemb;k++){
    
    if(domin) minptr[k]=1E30;
    if(domax) maxptr[k]=-1E30;
    if(dosum) sumptr[k]=0;
  }
  ZLOOP {
    for(k=start;k<start+nmemb;k++){
      if(domax) if (p[i][j][k] > maxptr[k]) maxptr[k] = p[i][j][k];
      if(domin) if (p[i][j][k] < minptr[k]) minptr[k] = p[i][j][k];
      if(dosum) sumptr[k]+=p[i][j][k];
    }
  }
#if(USEMPI)
  for(k=start;k<start+nmemb;k++){    
    if(domax){
      maxsend = maxptr[k];
      MPI_Allreduce(&maxsend, &maxptr[k], 1, MPI_FTYPE, MPI_MAX, MPI_COMM_WORLD);
    }
    if(domin){
      minsend = minptr[k];
      MPI_Allreduce(&minsend, &minptr[k], 1, MPI_FTYPE, MPI_MIN, MPI_COMM_WORLD);
    }
    if(dosum){
      sumsend = sumptr[k];
      MPI_Allreduce(&sumsend, &sumptr[k], 1, MPI_FTYPE, MPI_SUM, MPI_COMM_WORLD);
    }
  }
#endif
}



void myfprintf(FILE* fileptr, char *format, ...)
{
  va_list arglist;
  if  (myid==0) {
    va_start (arglist, format);

    if(fileptr==NULL){
      fprintf(stderr,"tried to print to null file pointer: %s\n",format);
      fflush(stderr);
    }
    else{
      vfprintf (fileptr, format, arglist);
      fflush(fileptr);
    }
    va_end (arglist);
  }
}

// prints to stderr(only cpu=0) AND file pointer of choice (all cpus)
void dualfprintf(FILE* fileptr, char *format, ...)
{
  va_list arglist;

  va_start (arglist, format);

  if(fileptr==NULL){
    fprintf(stderr,"tried to print to null file pointer: %s\n",format);
    fflush(stderr);
  }
  else{
    vfprintf (fileptr, format, arglist);
    fflush(fileptr);
  }
  if(myid==0){
    vfprintf (stderr, format, arglist);
    fflush(stderr);
  }
  va_end (arglist);
}

// prints to both logfull_file(cpu=0) and log_file(all cpus)
void logsfprintf(char *format, ...)
{
  va_list arglist;
  
  va_start (arglist, format);

  if  ((myid==0)&&(logfull_file)){
    vfprintf (logfull_file, format, arglist);    
    fflush(logfull_file);
  }
  if(log_file){
    vfprintf (log_file, format, arglist);    
    fflush(log_file);
  }
  va_end (arglist);
}

// prints to logfull_file, log_file, and stderr (but only using cpu=0)
void trifprintf(char *format, ...)
{
  va_list arglist;
  
  va_start (arglist, format);
  if  ((myid==0)&&(logfull_file)){
    vfprintf (logfull_file, format, arglist);    
    fflush(logfull_file);
  }
  if(log_file){
    vfprintf (log_file, format, arglist);    
    fflush(log_file);
  }
  if(myid==0){
    vfprintf (stderr, format, arglist);    
    fflush(stderr);
  }
  va_end (arglist);
}


// cpu==0 opens a file
void myfopen(char*fname, char*fmt,char*message,FILE**fileptrptr)
{
  if(myid==0){
    *fileptrptr = fopen(fname, fmt);
    if (*fileptrptr == NULL) {
      dualfprintf(fail_file, message);
      myexit(10000);
    }
  }
}

// cpu==0 closes a file
void myfclose(FILE ** fileptrptr,char*message)
{
  int reterror;

  if(myid==0){
    if(*fileptrptr!=NULL){
      reterror = fclose(*fileptrptr);
      if (reterror == EOF) {
	dualfprintf(fail_file, message);
	myexit(10000);
      }
    }
    else{
      dualfprintf(fail_file,"file already closed: %s\n",message);
      myexit(10000);
    }
  }
}


// write to file or MPI buffer
void myfwrite(int bintxt, MPI_Datatype datatype, void *ptr, int start, size_t nmemb, int i, int j, FILE**stream,void*writebuf)
{
  int k;
  void *voidbuf;
  int sizeofdatatype;

  static long double *writebuf16;
  static double *writebuf8;
  static float *writebuf4;
  static unsigned char *writebuf1;
  static int *writebuf4i;
  static long long int *writebuf8i;

  static long double *ptr16;
  static double *ptr8;
  static float *ptr4;
  static unsigned char *ptr1;
  static int *ptr4i;
  static long long int *ptr8i;

  int streamnum;



  sizeofdatatype=getsizeofdatatype(datatype);

  if (datatype == MPI_UNSIGNED_CHAR) writebuf1 = (unsigned char *) writebuf;
  else if (datatype == MPI_FLOAT) writebuf4 = (float *) writebuf;
  else if (datatype == MPI_DOUBLE) writebuf8 = (double *) writebuf;
  else if (datatype == MPI_LONG_DOUBLE) writebuf16 = (long double *) writebuf;
  else if (datatype == MPI_INT) writebuf4i = (int *) writebuf;
  else if (datatype == MPI_LONG_LONG_INT) writebuf8i = (long long int *) writebuf;
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(1);
  }

  if (datatype == MPI_UNSIGNED_CHAR) ptr1 = (unsigned char *) ptr;
  else if (datatype == MPI_FLOAT) ptr4 = (float *) ptr;
  else if (datatype == MPI_DOUBLE) ptr8 = (double *) ptr;
  else if (datatype == MPI_LONG_DOUBLE) ptr16 = (long double *) ptr;
  else if (datatype == MPI_INT) ptr4i = (int *) ptr;
  else if (datatype == MPI_LONG_LONG_INT) ptr8i = (long long int *) ptr;
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(1);
  }

  if(mpicombine==0){
    if(bintxt==BINARYOUTPUT){ // binaryoutput
      if(!docolsplit){
	fwrite(ptr+(sizeofdatatype*start), sizeofdatatype, nmemb, stream[0]);
	nextbuf+=nmemb;
      }
      else{
	for(k=start;k<start+nmemb;k++){
	  streamnum=nextbuf;
	  fwrite(ptr+(sizeofdatatype*k), sizeofdatatype, 1, stream[streamnum]);
	  nextbuf++;
	}
      }
    }
    else{ // text output
      for(k=start;k<start+nmemb;k++){
	if(docolsplit) streamnum=nextbuf;
	else streamnum=0;
	if (datatype == MPI_UNSIGNED_CHAR) fprintf(stream[streamnum],"%c ",ptr1[k]);
	else if (datatype == MPI_FLOAT) fprintf(stream[streamnum],"%15.7g ",ptr4[k]);
	else if (datatype == MPI_DOUBLE) fprintf(stream[streamnum],"%21.15g ",ptr8[k]);
	else if (datatype == MPI_LONG_DOUBLE) fprintf(stream[streamnum],"%31.25Lg ",ptr16[k]);
	else if (datatype == MPI_INT) fprintf(stream[streamnum],"%d ",ptr4i[k]);
	else if (datatype == MPI_LONG_LONG_INT) fprintf(stream[streamnum],"%lld ",ptr8i[k]);
	nextbuf++;
      }
    }
  }
  else{ // mpicombine==1
    if(docolsplit&&USEROMIO){ // column splitting with ROMIO
      for(k=start;k<start+nmemb;k++){
	if(nextbuf==romiocoliter){ // only write if doing that column
	  // BUFFERMAP2 only depends on i,j, not column number
	  if (datatype == MPI_UNSIGNED_CHAR) writebuf1[BUFFERMAP2] = ptr1[k];
	  if (datatype == MPI_FLOAT) writebuf4[BUFFERMAP2] = ptr4[k];
	  if (datatype == MPI_DOUBLE) writebuf8[BUFFERMAP2] = ptr8[k];
	  if (datatype == MPI_LONG_DOUBLE) writebuf16[BUFFERMAP2] = ptr16[k];
	  if (datatype == MPI_INT) writebuf4i[BUFFERMAP2] = ptr4i[k];
	  if (datatype == MPI_LONG_LONG_INT) writebuf8i[BUFFERMAP2] = ptr8i[k];
	}
	nextbuf++;
      }
    }
    else{ // no ROMIO column splitting, just normal MPI buffer writing
      for(k=start;k<start+nmemb;k++){
	if (datatype == MPI_UNSIGNED_CHAR) writebuf1[BUFFERMAP] = ptr1[k];
	if (datatype == MPI_FLOAT) writebuf4[BUFFERMAP] = ptr4[k];
	if (datatype == MPI_DOUBLE) writebuf8[BUFFERMAP] = ptr8[k];
	if (datatype == MPI_LONG_DOUBLE) writebuf16[BUFFERMAP] = ptr16[k];
	if (datatype == MPI_INT) writebuf4i[BUFFERMAP] = ptr4i[k];
	if (datatype == MPI_LONG_LONG_INT) writebuf8i[BUFFERMAP] = ptr8i[k];
      }
    }
  }
}

// same kind of process as myfwrite, see comments there
void myfread(int bintxt, MPI_Datatype datatype, void *ptr, int start, size_t nmemb, int i, int j, FILE**stream,void*writebuf)
{
  int k;
  void *voidbuf;
  int sizeofdatatype;

  static long double *writebuf16;
  static double *writebuf8;
  static float *writebuf4;
  static unsigned char *writebuf1;
  static int *writebuf4i;
  static long long int *writebuf8i;

  static long double *ptr16;
  static double *ptr8;
  static float *ptr4;
  static unsigned char *ptr1;
  static int *ptr4i;
  static long long int *ptr8i;

  int streamnum;



  sizeofdatatype=getsizeofdatatype(datatype);

  if (datatype == MPI_UNSIGNED_CHAR) writebuf1 = (unsigned char *) writebuf;
  else if (datatype == MPI_FLOAT) writebuf4 = (float *) writebuf;
  else if (datatype == MPI_DOUBLE) writebuf8 = (double *) writebuf;
  else if (datatype == MPI_LONG_DOUBLE) writebuf16 = (long double *) writebuf;
  else if (datatype == MPI_INT) writebuf4i = (int *) writebuf;
  else if (datatype == MPI_LONG_LONG_INT) writebuf8i = (long long int *) writebuf;
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(1);
  }

  if (datatype == MPI_UNSIGNED_CHAR) ptr1 = (unsigned char *) ptr;
  else if (datatype == MPI_FLOAT) ptr4 = (float *) ptr;
  else if (datatype == MPI_DOUBLE) ptr8 = (double *) ptr;
  else if (datatype == MPI_LONG_DOUBLE) ptr16 = (long double *) ptr;
  else if (datatype == MPI_INT) ptr4i = (int *) ptr;
  else if (datatype == MPI_LONG_LONG_INT) ptr8i = (long long int *) ptr;
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(1);
  }

  if(mpicombine==0){
    if(bintxt==BINARYOUTPUT){
      if(!docolsplit){
	fread(ptr+(sizeofdatatype*start), sizeofdatatype, nmemb, stream[0]);
	nextbuf+=nmemb;
      }
      else{
	for(k=start;k<start+nmemb;k++){
	  streamnum=nextbuf;
	  fread(ptr+(sizeofdatatype*k), sizeofdatatype, 1, stream[streamnum]);
	  nextbuf++;
	}
      }
    }
    else{
      for(k=start;k<start+nmemb;k++){
	if(docolsplit) streamnum=nextbuf;
	else streamnum=0;
	if (datatype == MPI_UNSIGNED_CHAR) fscanf(stream[streamnum],"%uc",&ptr1[k]);
	else if (datatype == MPI_FLOAT) fscanf(stream[streamnum],"%f",&ptr4[k]);
	else if (datatype == MPI_DOUBLE) fscanf(stream[streamnum],"%lf",&ptr8[k]);
	else if (datatype == MPI_LONG_DOUBLE) fscanf(stream[streamnum],"%Lf",&ptr16[k]);
	else if (datatype == MPI_INT) fscanf(stream[streamnum],"%d",&ptr4i[k]);
	else if (datatype == MPI_LONG_LONG_INT) fscanf(stream[streamnum],"%lld",&ptr8i[k]);
	nextbuf++;
      }
    }
  }
  else{
    if(docolsplit&&USEROMIO){
      for(k=start;k<start+nmemb;k++){
	if(nextbuf==romiocoliter){
	  if (datatype == MPI_UNSIGNED_CHAR) ptr1[k]=writebuf1[BUFFERMAP2]; 
	  if (datatype == MPI_FLOAT) ptr4[k]=writebuf4[BUFFERMAP2]; 
	  if (datatype == MPI_DOUBLE) ptr8[k]=writebuf8[BUFFERMAP2]; 
	  if (datatype == MPI_LONG_DOUBLE) ptr16[k]=writebuf16[BUFFERMAP2]; 
	  if (datatype == MPI_INT) ptr4i[k]=writebuf4i[BUFFERMAP2]; 
	  if (datatype == MPI_LONG_LONG_INT) ptr8i[k]=writebuf8i[BUFFERMAP2]; 
	}
	nextbuf++;
      }
    }
    else for(k=start;k<start+nmemb;k++){
      if (datatype == MPI_UNSIGNED_CHAR) ptr1[k]=writebuf1[BUFFERMAP]; 
      if (datatype == MPI_FLOAT) ptr4[k]=writebuf4[BUFFERMAP]; 
      if (datatype == MPI_DOUBLE) ptr8[k]=writebuf8[BUFFERMAP]; 
      if (datatype == MPI_LONG_DOUBLE) ptr16[k]=writebuf16[BUFFERMAP]; 
      if (datatype == MPI_INT) ptr4i[k]=writebuf4i[BUFFERMAP]; 
      if (datatype == MPI_LONG_LONG_INT) ptr8i[k]=writebuf8i[BUFFERMAP]; 
    }
  }
}


// sets values between 2 pointers, typically to cumulate values into writebuf array for later use.
void myset(MPI_Datatype datatype, void *ptr, int start, size_t nmemb, void*writebuf)
{
  int k;

  static long double *writebuf16;
  static double *writebuf8;
  static float *writebuf4;
  static unsigned char *writebuf1;
  static int *writebuf4i;
  static long long int *writebuf8i;

  static long double *ptr16;
  static double *ptr8;
  static float *ptr4;
  static unsigned char *ptr1;
  static int *ptr4i;
  static long long int *ptr8i;

  if (datatype == MPI_UNSIGNED_CHAR) writebuf1 = (unsigned char *) writebuf;
  else if (datatype == MPI_FLOAT) writebuf4 = (float *) writebuf;
  else if (datatype == MPI_DOUBLE) writebuf8 = (double *) writebuf;
  else if (datatype == MPI_LONG_DOUBLE) writebuf16 = (long double *) writebuf;
  else if (datatype == MPI_INT) writebuf4i = (int *) writebuf;
  else if (datatype == MPI_LONG_LONG_INT) writebuf8i = (long long int *) writebuf;
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(1);
  }

  if (datatype == MPI_UNSIGNED_CHAR) ptr1 = (unsigned char *) ptr;
  else if (datatype == MPI_FLOAT) ptr4 = (float *) ptr;
  else if (datatype == MPI_DOUBLE) ptr8 = (double *) ptr;
  else if (datatype == MPI_LONG_DOUBLE) ptr16 = (long double *) ptr;
  else if (datatype == MPI_INT) ptr4i = (int *) ptr;
  else if (datatype == MPI_LONG_LONG_INT) ptr8i = (long long int *) ptr;
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(1);
  }

  for(k=start;k<start+nmemb;k++){
    // nextcol is the iterator here, a global variable currently
    if (datatype == MPI_UNSIGNED_CHAR) writebuf1[nextcol++] = ptr1[k];
    else if (datatype == MPI_FLOAT) writebuf4[nextcol++] = ptr4[k];
    else if (datatype == MPI_DOUBLE) writebuf8[nextcol++] = ptr8[k];
    else if (datatype == MPI_LONG_DOUBLE) writebuf16[nextcol++] = ptr16[k];
    else if (datatype == MPI_INT) writebuf4i[nextcol++] = ptr4i[k];
    else if (datatype == MPI_LONG_LONG_INT) writebuf8i[nextcol++] = ptr8i[k];
  }
}



// very similar to myset, just switched assignments
void myget(MPI_Datatype datatype, void *ptr, int start, size_t nmemb, void*writebuf)
{
  int k;

  static long double *writebuf16;
  static double *writebuf8;
  static float *writebuf4;
  static unsigned char *writebuf1;
  static int *writebuf4i;
  static long long int *writebuf8i;

  static long double *ptr16;
  static double *ptr8;
  static float *ptr4;
  static unsigned char *ptr1;
  static int *ptr4i;
  static long long int *ptr8i;

  if (datatype == MPI_UNSIGNED_CHAR) writebuf1 = (unsigned char *) writebuf;
  else if (datatype == MPI_FLOAT) writebuf4 = (float *) writebuf;
  else if (datatype == MPI_DOUBLE) writebuf8 = (double *) writebuf;
  else if (datatype == MPI_LONG_DOUBLE) writebuf16 = (long double *) writebuf;
  else if (datatype == MPI_INT) writebuf4i = (int *) writebuf;
  else if (datatype == MPI_LONG_LONG_INT) writebuf8i = (long long int *) writebuf;
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(1);
  }

  if (datatype == MPI_UNSIGNED_CHAR) ptr1 = (unsigned char *) ptr;
  else if (datatype == MPI_FLOAT) ptr4 = (float *) ptr;
  else if (datatype == MPI_DOUBLE) ptr8 = (double *) ptr;
  else if (datatype == MPI_LONG_DOUBLE) ptr16 = (long double *) ptr;
  else if (datatype == MPI_INT) ptr4i = (int *) ptr;
  else if (datatype == MPI_LONG_LONG_INT) ptr8i = (long long int *) ptr;
  else{
    dualfprintf(fail_file,"No such datatype=%d\n",datatype);
    myexit(1);
  }

  for(k=start;k<start+nmemb;k++){
    // nextcol is the iterator here, a global variable currently
    if (datatype == MPI_UNSIGNED_CHAR) ptr1[k] = writebuf1[nextcol++];
    else if (datatype == MPI_FLOAT) ptr4[k] = writebuf4[nextcol++];
    else if (datatype == MPI_DOUBLE) ptr8[k] = writebuf8[nextcol++];
    else if (datatype == MPI_LONG_DOUBLE) ptr16[k] = writebuf16[nextcol++];
    else if (datatype == MPI_INT) ptr4i[k] = writebuf4i[nextcol++];
    else if (datatype == MPI_LONG_LONG_INT) ptr8i[k] = writebuf8i[nextcol++];
  }
}

