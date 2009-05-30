#define EVOLVECHECKS 1
// whether to check boundary conditions and limit gamma during advance()

#define FIXUP1ZONE 0
#define FIXUPALLZONES 1
#define FIXUPNOZONES 2

#define FIXUPZONES FIXUP1ZONE

/** algorithmic choices **/

#define HLLBOUNDARY 0
// use HLL on boundary fluxes

/* mnemonics for flux method (Riemann solver) */
// ordered from most diffusive to least diffusive, so can back track
// 0 should be reasonable most diffusive
#define LAXFFLUX 0
#define HLLFLUX 1


// DIVB constraint method
#define FLUXCTHLL 0
#define FLUXCTTOTH 1
#define FLUXCD 2
#define ATHENA1 3
#define ATHENA2 4
/* these are different ways of calculating the EMFs */
//#define FLUXB FLUXCTTOTH
// 0: HLL
// 1: FLUXCT TOTH version (toth 2000 eq. 25)
// 2: FLUXCD TOTH version (toth 2000 eq. 31)
// 3: Athena type eq 39
// 4: Athena type eq 48


//#define UTOPRIMVERSION 6
// 0: original gammie 5D method
#define UTOPRIM5D1 0
// 1: ldz method
#define UTOPRIMLDZ 1
// 2: SCN 2D method
#define UTOPRIM2D 2
// 3: SCN 1D method
#define UTOPRIM1D 3
// 4: SCN 1D OPTIMIZED method -- not sure if identical to 3 otherwise
#define UTOPRIM1DOPT 4
// 5: SCN 1D final and optimized
#define UTOPRIM1DFINAL 5
// 6: SCN 2D final and optimized and recommended by Scott
#define UTOPRIM2DFINAL 6
// 7: SCN 5D final -- bit less accurate compared to 1D and 2D
#define UTOPRIM5D2 7
// 100: use 5D, but compare with ldz in runtime
#define UTOPRIMCOMPARE 100


/* mnemonics for slope limiter */
// ordered from most diffusive to least diffusive, so can start high and go down if needed
// 0 should be reasonble most diffusive
#define DONOR	0
#define VANL	1
#define MINM	2
#define MC      3
#define PARA    4
#define CSSLOPE      90 // not tested/compared against others
#define NLIM    100 // no limiter
#define NLIMCENT    101 // no limiter
#define NLIMUP    102 // no limiter
#define NLIMDOWN    103 // no limiter


#define MAXTIMEORDER 4

//#define TIMEORDER 3
// order of algorithm in time from 1 to 4.
// TIMEORDER: 1 : single step (Euler method -- error term is 2nd order for smooth flows)
// TIMEORDER: 2 : 2 steps in halfs (midpoint method -- error term is 3rd order for smooth flows)
// TIMEORDER: 3 : 4 steps (classic RK3 method -- error term is 4th order for smooth flows)
// TIMEORDER: 4 : 4 steps (classic RK4 method -- error term is 5th order for smooth flows)

#define FIXUPFLUX 1
// fix up the flux using fix_flux() in fixup.c
// caution...

#define ZEROOUTFLOWFLUX 0
// 0: don't do anything special, let boundary conditions control fux
// 1: zero the x1-velocity at the outflow boundary if inflow in flux routine so interpolation is overridden
// seems to cause problems at outer edge in internal energy

// seems to cause problems with EOMFFDE at boudaries.
// GODMARK



// assumes inner and outer theta boundaries are at 0 and Pi.
#define ZEROPOLEFLUX 0
// 0: don't do anything special
// 1: zero polar theta flux since should be 0

// seems to cause problems at 0 and pi boundaries with FFDE, or even just one boundary???
// GODMARK

// REASON: Cannot just set F=0.  The flux of the kinetic energy still has the pressure term even if the velocity into the wall is 0!



// whether to rescale interpolation
#define RESCALEINTERP 0
// 0: don't rescale
// 1: do rescale

#define BDIRCONT 1
// 0: don't
// 1: make field along flux direction continuous

// whether to use hyperbolic function for cmin and cmax instead of sharp max(0,l,r)
#define HYPERHLL 0

// whether to shift stencil inside horizon
#define HORIZONSUPERFAST 0
