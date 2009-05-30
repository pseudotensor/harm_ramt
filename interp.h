
/////////////////////////
//
// for interp.c




#include "step_ch.h" // for nemonics


// definition of minmod operator
#define MINMOD(a,b) ( ((a)*(b)<=0) ? 0.0 : MINMODB(a,b) )
#define MINMODB(a,b) ( (fabs(a)<fabs(b)) ? (a) : (b) )

///////////////////////////////
//
// parabolic interpolation stuff
//
////////////////////////////////

#define PARA1 0 // old
#define PARA2 1 // works
#define PARA3 2 // broken
#define PARA4 3 // latest

#define WHICHPARA PARA4

//
// below for para2() and para3()
//

#define PMC 100
// PMC: originally used in PPM paper -- Xiaoyue says fails on nonlinear tests since has no limiter -- but they used it, right?

#if(WHICHPARA==PARA2)
#define PARA2LIM VANL
// jon tests show for PARA2:
// PARA2LIM = MC completely fails
// PARA2LIM = VANL best ok
// PARA2LIM = PMC kinda ok
#elif(WHICHPARA==PARA3)
// MC recommended by Xiaoyue
#define PARA2LIM MC
#elif(WHICHPARA==PARA4)
// MC recommended by Xiaoyue
#define PARA2LIM MC
#endif




//////////////////////////////////
//
// which variable to interpolate
//
/////////////////////////////////

#define PRIMTOINTERP 0
#define CONSTOINTERP 1

#define VARTOINTERP PRIMTOINTERP
