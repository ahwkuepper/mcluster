//Constants
#define G 0.0043009211 //in pc*km^2*s^-2*Msun
#define Pi 3.14159265
#define PI 3.14159265
#define TWOPI 6.2831853   /* 2PI */
#define GNBODY 1.0
#define MNBODY 1.0
#define PARSEC  3.08568E13      /* KM pro PC */
#define GBIN    1.327126E11     /* G in (km/sec)^3/Msun */
#define RSUN    6.96265E5       /* Solar radius in km */

//Functions
#define max(a,b)         (a < b) ?  (b) : (a)
#define min(a,b)         (a < b) ?  (a) : (b)
#define Lifetime(Mstar)  1.13E4*pow(Mstar,-3)+0.6E2*pow(Mstar,-0.75)+1.2 //Myr	[Prantzos 2007]


//MW potential consisting of bulge, disk & NFW halo - constants:
//based on Law, Majewski & Johnston (2009)
#define b1_LMJ 0.7000;       //[kpc]
#define M1_LMJ 3.4e10;       //[solar masses]
#define a2_LMJ 6.5000;       //[kpc]
#define b2_LMJ 0.2600;        //[kpc]
#define M2_LMJ 1.0e11;       //[solar masses]
#define R_NFW 42.9647;       //[kpc]
#define M_NFW 9.832538e+11;  //[solar masses]
#define q_halo 1.0;          //[dimensionless]

//Allen & Santillan (1991) MW potential - constants:
#define b1allen 0.3873            //kpc
#define M1allen 606.0*2.32e07    //solar masses
#define a2allen 5.3178
#define b2allen 0.2500
#define M2allen 3690.0*2.32e07
#define a3allen 12.0000
#define M3allen 4615.0*2.32e07

//Additional parameters for Sverre's MW potential
#define VCIRC 220.0			//circular velocity at RCIRC
#define RCIRC 8.500			//kpc

//new input parameters for bulge potential
#define GMB 0.0				//Central bulge mass (Msun).
#define AR 1.0				//Scale radius in gamma/eta model (kpc, Dehnen 1993).
#define GAM 1.0				//gamma (in the range 0 =< gamma < 3).

//Point-mass potential - constants:
#define M1pointmass 9.565439E+10    //solar masses

//Constants for user-defined profile
#define EPS 1.E-05  //precision of integrator
#define JMAX 5      //maximum number of integration steps pow(3, JMAX)
#define BIGNUMBER 1.0e30

//Mass function variables
#define MAX_AN  10
#define MAX_MN  11

//Plumix
#define SWAP(a,b) { temp = star1[a].mass; star1[a].mass = star1[b].mass; star1[b].mass = temp; temp = star1[a].m0; star1[a].m0 = star1[b].m0; star1[b].m0 = temp; temp = star1[a].kw; star1[a].kw = star1[b].kw; star1[b].kw = temp; temp = star1[a].epoch; star1[a].epoch = star1[b].epoch; star1[b].epoch = temp; temp = star1[a].spin; star1[a].spin = star1[b].spin; star1[b].spin = temp; temp = star1[a].rstar; star1[a].rstar = star1[b].rstar; star1[b].rstar = temp; temp = star1[a].lum; star1[a].lum = star1[b].lum; star1[b].lum = temp; temp = star1[a].epochstar; star1[a].epochstar = star1[b].epochstar; star1[b].epochstar = temp; temp = star1[a].zstar; star1[a].zstar = star1[b].zstar; star1[b].zstar = temp;};
#define _G      6.673e-8
#define _pc     3.0856e18
#define _M_sun  1.989e33
#define _R_sun  6.96265e10
#define _AU     1.49597870e13
#define _GMsun  1.3272597e26
#define RSCALE  0.58904862254809  // = 3pi / 16 = 3pi / (16 * |U_tot|)
#define BUFF_STEP 1024





struct t_star1
{
	double mass, U_tmp, Ui;
	double r[3];
	double m0, kw, epoch, spin, rstar, lum, epochstar, zstar;
};
struct t_star2
{
	double v[3];
	double U, U_sub, E;
	double UT_add, UT;
	double rad, M_sub;
	long   ntry;
};
struct t_star1 *star1;
struct t_star2 *star2;

//functions and COMMON blocks for SSE (Hurley, Pols & Tout, 2000, MNRAS, 315, 543) and BSE (Hurley, Tout & Pols, 2002, MNRAS, 329, 897)
#ifdef SSE
extern void zcnsts_(double *z, double *zpars);
extern void evolv1_(int *kw, double *mass, double *mt, double *r, double *lum, double *mc, double *rc, double *menv, double *renv, double *ospin, double *epoch, double *tms, double *tphys, double *tphysf, double *dtp, double *z, double *zpars, double *vkick);
extern void evolv2_(int *kw, double *mass, double *mt, double *r, double *lum, double *mc, double *rc, double *menv, double *renv, double *ospin, double *epoch, double *tms, double *tphys, double *tphysf, double *dtp, double *z, double *zpars, double *tb, double *ecc, double *vkick);
#else
void zcnsts_(double *z, double *zpars) {/*DUMMY*/};
void evolv1_(int *kw, double *mass, double *mt, double *r, double *lum, double *mc, double *rc, double *menv, double *renv, double *ospin, double *epoch, double *tms, double *tphys, double *tphysf, double *dtp, double *z, double *zpars, double *vkick) {/*DUMMY*/};
void evolv2_(int *kw, double *mass, double *mt, double *r, double *lum, double *mc, double *rc, double *menv, double *renv, double *ospin, double *epoch, double *tms, double *tphys, double *tphysf, double *dtp, double *z, double *zpars, double *tb, double *ecc, double *vkick) {/*DUMMY*/};
#endif

struct{
	double neta;
	double bwind;
	double hewind;
	double mxns;
} value1_;

struct{
	double pts1;
	double pts2;
	double pts3;
} points_;

struct{
	double sigma;
	int bhflag;
} value4_;

struct{
	int idum;
} value3_;

struct{
	int ceflag;
	int tflag;
	int ifflag;
	int nsflag;
	int wdflag;
} flags_;

struct{
	double beta;
	double xi;
	double acc2;
	double epsnov;
	double eddfac;
	double gamma;
} value5_;

struct{
	double alpha1;
	double lambda;
} value2_;
