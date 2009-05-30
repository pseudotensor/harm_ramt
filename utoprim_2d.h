#include "utoprim_1d2d.h"

#define NEWT_DIM 2                      

#define MAX_NEWT_RETRIES 0    /*Max. # of retries of N-R procedure, while increasing tol. by *10 after each time*/
#define CYCLE_BREAK_PERIOD 10 /* change newton step by random factor every this number of newton iterations*/  
#define CHECK_FOR_STALL 1     /* check for stationary newton stepping */                                       

#define SCALEMAX      1.0e2    /* Max. value of the factor used to scale the Newton step */
#define TOL_LINE_STEP NEWT_TOL /* Minimum value of Max(dx/x) in line search algorithm */
#define GRADMIN       1.0e-10  /* Magnitude of gradient below which we say we are at a local min. */
#define NEWT_FUNC_TOL 1.0e-5  /* Max. ratio of the final and initial resid magnitudes to be considered converged */



/* functions used in grmhd */
static void func_vsq(FTYPE x[2], FTYPE dx[2], FTYPE *f, FTYPE *df, int n);
static void bcon_calc_g(FTYPE prim[8],FTYPE ucon[4],FTYPE ucov[4],FTYPE ncov[4],FTYPE bcon[4])  ;
static void lower_g(FTYPE vcon[4], FTYPE gcov[4][4], FTYPE vcov[4]) ;
static void ncov_calc(FTYPE gcon[4][4],FTYPE ncov[4]) ;
static void raise_g(FTYPE vcov[4], FTYPE gcon[4][4], FTYPE vcon[4]) ;
static int Utoprim_new_body(FTYPE U[NPR], struct of_geom *ptrgeom,  FTYPE prim[NPR]);
static void ucon_calc_g(FTYPE prim[8],FTYPE gcov[4][4],FTYPE gcon[4][4],FTYPE ucon[4]) ;

static FTYPE dpdvsq_calc(FTYPE W, FTYPE vsq);
static FTYPE dpdW_calc_vsq(FTYPE W, FTYPE vsq);
static FTYPE pressure_W_vsq(FTYPE W, FTYPE vsq) ;
static FTYPE pressure_rho0_u(FTYPE rho0, FTYPE u) ;
static FTYPE pressure_rho0_w(FTYPE rho0, FTYPE w) ;
static FTYPE utsq_calc(FTYPE W) ;

static int find_root_2D_gen(FTYPE *W) ;

static int general_newton_raphson( FTYPE x[], int n, void (*funcd) (FTYPE [], FTYPE [], FTYPE *, FTYPE *, int) );  
static void validate_x(FTYPE x[2], FTYPE x0[2], int vartype ) ;
static FTYPE x1_of_x0(FTYPE x0, int vartype ) ;

/* "Debugging" routines, that may be useful in the end product */
//int    is_nan_inf( FTYPE x );

static void primtoU_g( FTYPE prim[8], FTYPE gcov[4][4], FTYPE gcon[4][4],FTYPE U[8] );


