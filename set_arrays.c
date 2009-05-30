
#include "decs.h"

void set_arrays()
{
  int i, j, k, l;
  int ii;
  int floor,pf, tscale,dtstage;
  FTYPE valueinit;

  p = (FTYPE (*)[N2M][NPR]) (&(a_p[N1BND][N2BND][0]));
  panalytic = (FTYPE (*)[N2M][NPR]) (&(a_panalytic[N1BND][N2BND][0]));
  omegafanalytic = (FTYPE (*)[N2M][NPR]) (&(a_omegafanalytic[N1BND][N2BND][0]));

  emf = (FTYPE (*)[N2M]) (&(a_emf[N1BND][N2BND]));
  vconemf = (FTYPE (*)[N2M][NDIM-1]) (&(a_vconemf[N1BND][N2BND][-U1]));

  pleft = (FTYPE (*)[N2M][NPR]) (&(a_pleft[N1BND][N2BND][0]));
  pright = (FTYPE (*)[N2M][NPR]) (&(a_pright[N1BND][N2BND][0]));

  unew = (FTYPE (*)[N2M][NPR]) (&(a_unew[N1BND][N2BND][0]));
  ulast = (FTYPE (*)[N2M][NPR]) (&(a_ulast[N1BND][N2BND][0]));

  dq1 = (FTYPE (*)[N2M][NPR]) (&(a_dq1[N1BND][N2BND][0]));
  dq2 = (FTYPE (*)[N2M][NPR]) (&(a_dq2[N1BND][N2BND][0]));

  

  F1 = (FTYPE (*)[N2M][NPR]) (&(a_F1[N1BND][N2BND][0]));
  F2 = (FTYPE (*)[N2M][NPR]) (&(a_F2[N1BND][N2BND][0]));
  F1CT = (FTYPE (*)[N2M][NPR]) (&(a_F1CT[N1BND][N2BND][0]));
  F2CT = (FTYPE (*)[N2M][NPR]) (&(a_F2CT[N1BND][N2BND][0]));
  pk = (FTYPE (*)[N1M][N2M][NPR]) (&(a_pk[0][N1BND][N2BND][0]));
  prc = (FTYPE (*)[N2M][NPR]) (&(a_prc[N1BND][N2BND][0]));

  pflag = (int (*)[N2M][NUMPFLAGS]) (&(a_pflag[N1BND][N2BND][0]));

#if(DODEBUG)
    failfloorcount = (CTYPE (*)[N2M][NUMTSCALES][NUMFAILFLOORFLAGS]) (&(a_failfloorcount[N1BND][N2BND][0]));
#endif
  // this faraday needed for current calculation
  cfaraday =  (FTYPE (*)[N2M][NUMCURRENTSLOTS][3]) (&(a_cfaraday[N1BND][N2BND][0][0]));
  fcon =  (FTYPE (*)[N2M][NUMFARADAY]) (&(a_fcon[N1BND][N2BND][0]));
  jcon = (FTYPE (*)[N2M][NDIM]) (&(a_jcon[N1BND][N2BND][0]));
  

  // assume time average stuff gets zeroed in avg routine
#if(DOAVG)
  normalvarstavg =  (FTYPE (*)[N2M][NUMNORMDUMP]) (&(a_normalvarstavg[N1BND][N2BND][0]));
  anormalvarstavg =  (FTYPE (*)[N2M][NUMNORMDUMP]) (&(a_anormalvarstavg[N1BND][N2BND][0]));

  fcontavg =  (FTYPE (*)[N2M][NUMFARADAY]) (&(a_fcontavg[N1BND][N2BND][0]));
  fcovtavg =  (FTYPE (*)[N2M][NUMFARADAY]) (&(a_fcovtavg[N1BND][N2BND][0]));

  afcontavg =  (FTYPE (*)[N2M][NUMFARADAY]) (&(a_afcontavg[N1BND][N2BND][0]));
  afcovtavg =  (FTYPE (*)[N2M][NUMFARADAY]) (&(a_afcovtavg[N1BND][N2BND][0]));

  massfluxtavg =  (FTYPE (*)[N2M][NDIM]) (&(a_massfluxtavg[N1BND][N2BND][0]));
  amassfluxtavg =  (FTYPE (*)[N2M][NDIM]) (&(a_amassfluxtavg[N1BND][N2BND][0]));

  othertavg =  (FTYPE (*)[N2M][NUMOTHER]) (&(a_othertavg[N1BND][N2BND][0]));
  aothertavg =  (FTYPE (*)[N2M][NUMOTHER]) (&(a_aothertavg[N1BND][N2BND][0]));

  jcontavg = (FTYPE (*)[N2M][NDIM]) (&(a_jcontavg[N1BND][N2BND][0]));
  jcovtavg = (FTYPE (*)[N2M][NDIM]) (&(a_jcovtavg[N1BND][N2BND][0]));

  ajcontavg = (FTYPE (*)[N2M][NDIM]) (&(a_ajcontavg[N1BND][N2BND][0]));
  ajcovtavg = (FTYPE (*)[N2M][NDIM]) (&(a_ajcovtavg[N1BND][N2BND][0]));

  tudtavg = (FTYPE (*)[N2M][NUMSTRESSTERMS]) (&(a_tudtavg[N1BND][N2BND][0]));
  atudtavg = (FTYPE (*)[N2M][NUMSTRESSTERMS]) (&(a_atudtavg[N1BND][N2BND][0]));
#endif  

#if(DOENO)
  uenotmp0 = (FTYPE (*)[N2M][NPR]) (&(a_uenotmp0[N1BND][N2BND][0]));
  uenotmp1 = (FTYPE (*)[N2M][NPR]) (&(a_uenotmp1[N1BND][N2BND][0]));
  uenotmp2 = (FTYPE (*)[N2M][NPR]) (&(a_uenotmp2[N1BND][N2BND][0]));
#endif
  

  // initialize things to NAN in order to (hopefully) trigger memory leaks to be noticed
  //  valueinit=sqrt(-1.0);
  valueinit=0.0;
  FULLLOOP {
    if(DODEBUG) TSCALELOOP FLOORLOOP failfloorcount[i][j][tscale][floor]=valueinit;
    PFLAGLOOP pflag[i][j][pf] = -100;
    PLOOP {
      p[i][j][k] = valueinit;
      panalytic[i][j][k] = valueinit;
      omegafanalytic[i][j][k] = valueinit;
      pleft[i][j][k] = valueinit;
      pright[i][j][k] = valueinit;
      unew[i][j][k] = valueinit;
      DTSTAGELOOP pk[dtstage][i][j][k] = valueinit;
      prc[i][j][k] = valueinit;
      dq1[i][j][k] = valueinit;
      dq2[i][j][k] = valueinit;
      F1[i][j][k] = valueinit;
      F2[i][j][k] = valueinit;
      F1CT[i][j][k] = valueinit;
      F2CT[i][j][k] = valueinit;
    }
    for(k=0;k<NUMCURRENTSLOTS;k++) for(l=0;l<3;l++){
      cfaraday[i][j][k][l]=valueinit;
    }
    for(k=0;k<NUMFARADAY;k++){
      fcon[i][j][k]=valueinit;
    }
    for(k=0;k<NDIM;k++){
      jcon[i][j][k]=valueinit;
    }
  }


  /* grid functions */
  conn = (FTYPE (*)[N2M][NDIM][NDIM][NDIM])
      (&(a_conn[N1BND][N2BND][0][0][0]));
  conn2 = (FTYPE (*)[N2M][NDIM])
      (&(a_conn2[N1BND][N2BND][0]));
  gcon = (FTYPE (*)[N2M][NPG][NDIM][NDIM])
      (&(a_gcon[N1BND][N2BND][0][0][0]));
  gcov = (FTYPE (*)[N2M][NPG][NDIM][NDIM])
      (&(a_gcov[N1BND][N2BND][0][0][0]));
  gdet = (FTYPE (*)[N2M][NPG])
      (&(a_gdet[N1BND][N2BND][0]));
  eomfunc = (FTYPE (*)[N2M][NPG])
      (&(a_eomfunc[N1BND][N2BND][0]));
  /*
  dxdxp = (FTYPE (*)[N2M][NPG][NDIM][NDIM])
      (&(a_dxdxp[N1BND][N2BND][0][0][0]));
  */
#if(VOLUMEDIFF)
    idxvol = (FTYPE (*)[N2M][NDIM])(&(a_idxvol[N1BND][N2BND][0]));
#endif


#if(DOLUMVSR)
    // yes, for each cpu
    lumvsr=(SFTYPE*)malloc(ncpux1*N1*sizeof(SFTYPE));
    if(lumvsr==NULL){
      dualfprintf(fail_file,"Couldn't open lumvsr memory\n");
      myexit(1);
    }
  
    lumvsr_tot=(SFTYPE*)malloc(ncpux1*N1*sizeof(SFTYPE));
    if(lumvsr_tot==NULL){
      dualfprintf(fail_file,"Couldn't open lumvsr_tot memory\n");
      myexit(1);
    }
#endif
    //for(ii=0;ii<ncpux1*N1;ii++) lumvsr[ii]=0;
    //for(ii=0;ii<ncpux1*N1;ii++) lumvsr_tot[ii]=0;
}
