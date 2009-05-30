#include "decs.h"


void bound_mpi_int(int boundstage, int prim[][N2M][NUMPFLAGS])
{
  int dir;

#if(USEMPI)
  /* These arrays contain designations that identify 
   * each recv and send */
  static MPI_Request requests[COMPDIM * 2 * 2];
  // format of map for requests[dir*2+recv/send(0/1)]
#endif

#if(USEMPI)

  /*
   *
   * 1. do outflow/inflow/etc. boundary conditions (in bounds)
   * 2. pack data into workbc arrays and do send/recv's
   * 3. for each transfer do a wait/unpack
   * 4. do all send waits
   *
   * NOTE: order is important so corner zones are assigned correctly
   *
   * workbc[PACK][][] is data that is being sent
   * workbc[UNPACK][][] is data that is being recvd
   *
   */

  /*
   *  -> higher values of x1
   * 2 -> lower values of x1
   *
   */

  
  /* bounds has already done non-MPI boundary conditions; now do MPI
   * boundary conditions */

  // must go in this order (0,1) then (2,3) or visa versa, a single
  // directions order doesn't matter.  This method of ordering is as
  // opposed to directly transfering the corner zones to the corner
  // CPU.

  // This may be faster since all transfers can proceed at once,
  // although this may be slower since no transfers can occur until
  // packing is completed.  This way packing and transfering occur
  // simultaneously.

  // Although l/r are packed together, since in the end we have to
  // wait for both l/r to complete, so equal time completion is
  // favored //over asynch completion.

  // Also, transfering corner zones with small message sizes increases
  // the importance of latency.

  // I choose left-right N1M/N2M first, then up/down N1M/N2M.  Could
  // just do N1/N2 for interior for L/R, but true boundary needs full
  // N1M/N2M exchanged since cpu sets boundary using normal bc code
  // which needs to get transfered to neight(i.e. currently if corner
  // has bctype 99/?  then doesn't do corner)

  
  if((boundstage==STAGE0)||(boundstage==STAGEM1)){
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[dir][DIRIF]) pack_int(dir,prim,workbc_int);
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[dir][DIRIF]) sendrecv_int(dir,workbc_int,requests);
  }
  if((boundstage==STAGE1)||(boundstage==STAGEM1)){
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[dir][DIRIF]) recvwait(dir,requests);
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[dir][DIRIF]) unpack_int(dir,workbc_int,prim);
    for(dir=X1UP;dir<=X1DN;dir++) if(dirset[dir][DIRIF]) sendwait(dir,requests);
  }
  if((boundstage==STAGE2)||(boundstage==STAGEM1)){
    // now dir=0,1(X1UP,X1DN) is done, so can start 2,3(X2UP,X2DN)
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[dir][DIRIF]) pack_int(dir,prim,workbc_int);
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[dir][DIRIF]) sendrecv_int(dir,workbc_int,requests);
  }
  if((boundstage==STAGE3)||(boundstage==STAGEM1)){
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[dir][DIRIF]) recvwait(dir,requests);
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[dir][DIRIF]) unpack_int(dir,workbc_int,prim);
    for(dir=X2UP;dir<=X2DN;dir++) if(dirset[dir][DIRIF]) sendwait(dir,requests);
  }



  // end if mpi
#endif

}	
// end function


// PACKLOOP allows one to alter which i or j is faster iterated
#define PACKLOOP_INT(i,j,istart,istop,jstart,jstop) GENLOOP(i,j,istart,istop,jstart,jstop) FLOOP

// packs data for shipment
void pack_int(int dir,int prim[][N2M][NUMPFLAGS],int workbc_int[][COMPDIM * 2][NBIGBND * NBIGSM])
{
  // dir=direction sending
  int i,j,k;
  int bci;
  int start1,start2,stop1,stop2;

  bci=0;
  PACKLOOP_INT(i,j
	   ,dirset[dir][DIRPSTART1]
	   ,dirset[dir][DIRPSTOP1]
	   ,dirset[dir][DIRPSTART2]
	   ,dirset[dir][DIRPSTOP2]){
    /*
    if(bci>=dirset[dir][DIRSIZE]){
      dualfprintf(fail_file,"pack memory leak: bci: %d dirset[%d][DIRSIZE]: %d\n",bci,dirset[dir][DIRSIZE]);
      myexit(10);
    }
    */
    workbc_int[PACK][dir][bci++] = prim[i][j][k];
  }
}

#if(USEMPI)
void sendrecv_int(int dir,int workbc_int[][COMPDIM * 2][NBIGBND * NBIGSM],MPI_Request *requests)
{
  MPI_Irecv(workbc_int[UNPACK][dir],
	    dirset[dir][DIRSIZE]/NPR,
	    MPI_INT,
	    dirset[dir][DIROTHER],
	    dirset[dir][DIRTAGR],
	    MPI_COMM_WORLD,
	    &requests[dir*2+REQRECV]);

  MPI_Isend(workbc_int[PACK][dir],
	    dirset[dir][DIRSIZE]/NPR,
	    MPI_INT,
	    dirset[dir][DIROTHER],
	    dirset[dir][DIRTAGS],
	    MPI_COMM_WORLD,
	    &requests[dir*2+REQSEND]);

}
#endif


void unpack_int(int dir,int workbc_int[][COMPDIM * 2][NBIGBND * NBIGSM],int prim[][N2M][NUMPFLAGS])
{
  // dir is direction receiving from
  int i,j,k;
  int bci;

  bci=0;
  PACKLOOP_INT(i,j
	   ,dirset[dir][DIRUSTART1]
	   ,dirset[dir][DIRUSTOP1]
	   ,dirset[dir][DIRUSTART2]
	   ,dirset[dir][DIRUSTOP2]){
    prim[i][j][k]=workbc_int[UNPACK][dir][bci++];
  }
}

