//Constants
#define G 0.00449850214 //in pc*km^2*s^-2*Msun
#define Pi 3.14159265
#define PI 3.14159265
#define TWOPI 6.2831853   /* 2PI */
#define GNBODY 1.0
#define MNBODY 1.0
#define PARSEC  3.08568E13      /* KM pro PC */
#define GBIN    1.327126E11     /* G in (km/sec)^3/Msun */
#define RSUN    6.96265E5       /* Solar radius in km */
#define RSUN2PC_MC 2.25643942E-08 /* conversion solar radius in pc*/

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

//external value from mocca.ini
extern struct
{
   int npop[10], initmodel[10], imfg[10],\
       pairing[10], adis[10],\
       eigen[10];
} mclusterarri_;  

extern struct
{
   double fracb[10], w0[10], conc_pop[10], Seg[10],\
          fractal[10], equalmass[10], mlow[10],\
          mup[10], alpha_L3[10], beta_L3[10],\
          mu_L3[10], amin[10],\
          amax[10], epoch_pop[10], zini_pop[10];
} mclusterarrd_; 

extern struct
{
   int potential_energy, tf, mclusteron,\
       seedmc, outputf, check_en;
} mclusteri_;

extern struct
{
   char alphaimfchar[10000],\
	  mlimimfchar[10000];

} mclusterchar_;
extern struct
{
   double qvir, rbar, rh_mcl;
} mclusterd_; 

int generate_m1(int *N, double **star, double mlow, double mup, double *M, double *mmean, double MMAX, double Mcl, double epoch, double Z, double Rh, int remnant, int N2);
int generate_m2(int an, double *mlim, double *alpha, double Mcl, double M_tmp, double *subcount, int *N, double *mmean, double *M, double **star, double MMAX, double epoch, double Z, double Rh, int remnant, int N2);
int generate_m3(int *N, double **star, double mlow, double mup, double *M, double *mmean, double MMAX, double Mcl, int N2);
double subint(double min, double max, double alpha);
double mlow_func(double mhigh, double alpha, double norma, double delta);
int generate_m4(int *N, double **star, double alpha, double beta, double mu,  double mlow, double mup, double *M, double *mmean, double MMAX, double Mcl, double epoch, double Z, double Rh, int remnant, int N2);
double alogam(double x, int *ifault);
double betain(double x, double p, double q, double beta, int *ifault);
double r8_abs(double x);
int generate_plummer(int N, double **star, double rtide, double rvir, double D, int symmetry, double Q, int N2, double **rho_dens, double conc, double M, double mtot);
int generate_king(int N, double W0, double **star, double *rvir, double *rh, double *rking, double D, int symmetry, int N2, double **rho_dens, double cc, double Mass, double mtot); 
double densty(double z);
int odeint(double ystart0, double ystart1, double x1, double x2, double den, int *kount, double *xp, double **yp, int M, int KMAX);	 
int derivs(double x, double *y, double *dydx, double den);
int rkqc(double *y,double *dydx, double *x, double *h,double den, double *yscal, double *hdid, double *hnext, double TOL);	
int rk4(double x, double *y, double *dydx, double h, double *yout, double den);
int generate_subr(int N, double S, double **star, double rtide, double rvir, int N2);
void quick(int start, int stop);
void position(int id, double U_sub, double UT_sub, double S);
void find_beta(double *avg, double *beta);
void isorand(double *r);
int cmpmy(double *x1, double *x2);
int cmpmy_reverse(double *x1, double *x2);
double generate_profile (int N, double **star, double Rmax, double Mtot, double *p, double *Rh, double D, int symmetry, int N2);
double dfridr(double (*func)(double, double*), double x, double h, double *err, double *p);
double midexp(double (*func)(double, double*), double aa, double bb, int n, double *p);
double midsql(double (*func)(double, double*), double aa, double bb, int n, double *p);
double midsqu(double (*func)(double, double*), double aa, double bb, int n, double *p);
double midinf(double (*func)(double, double*), double aa, double bb, int n, double *p);
double midpnt(double (*func)(double, double*), double a, double b, int n, double *p);
double kernel (double x, double *p);
double rho(double r, double *p);
double rhoR (double x, double *p);
double rho_kernel (double x, double *p);
double M_func(double r, double *p);
double sigma_kernel (double x, double *p);
double sigma_func(double r, double *p);
double get_gauss(void);
double fractalize(double D, int N, double **star, int radial, int symmetry, int N2);
int standalone_rzamsf(double m, double *radius);
int get_binaries(int nbin, double **star, double M, double rvir, int pairing, int *N, int adis, double amin, double amax, double Rh, int Ntot, int eigen, int BSE, double epoch, double Z, int remnant, int OBperiods, double msort, int N2, int N3, double *eccbinaries, double *abinaries, double **cmb);
void shellsort_reverse_1d(double *array, int N);
void shellsort_1d(double *array, int N);	
void shellsort(double **array, int N, int k);
void shellsort_reverse(double **array, int N, int k);
int order(double **star, int N, double M, double msort, int pairing, int N2, double fracb);
int segregate(double **star, int N, double S, int N2);
int energy_order(double **star, int N, int Nstars, int N2);
int radial_distance_order(double **star, int N);
double coepot(double **star, int N);
int randomize(double **star, int N, int N2);
double rtnewt (double ecc, double ma);
int eigenevolution_old(double *m1, double *m2, double *ecc, double *abin, double M, double rvir);
int eigenevolution_new(double *m1, double *m2, double *ecc, double *abin, double M, double rvir);
int cmd(double **star, int l, double Rgal, double *abvmag, double *vmag, double *BV, double *Teff, double *dvmag, double *dBV);
int multiplearraytodouble(char str[], double **array);
int arraytodouble(char str[], double *array);
void info(char *output, int N, double Mcl, int profile, double W0, double S, double D, double Q, double Rh, double gamma[], double a, double Rmax, double tcrit, int tf, double RG[], double VG[], int mfunc, double single_mass, double mlow, double mup, double alpha[], double mlim[], double alpha_L3, double beta_L3, double mu_L3, int weidner, int mloss, int remnant, double epoch, double Z, int prantzos, int nbin, double fbin, int pairing, double msort, int adis, double amin, double amax, int eigen, int BSE, double extmass, double extrad, double extdecay, double extstart, int code, int seed, double dtadj, double dtout, double dtplot, int gpu, int regupdate, int etaupdate, int esc, int units, int match, int symmetry, int OBperiods);
double interpl_density(int N, int Ni, int i, int j_star, double **inputJE_vect, double **star_temp, double ***rho_dens);
double interpolation(double **x, int n, double a);
