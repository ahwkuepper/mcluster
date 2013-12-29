/***************************************************************************
 *   Copyright (C) 2009 by Andreas H.W. Kuepper                            *
 *   akuepper@astro.uni-bonn.de                                            *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
/***************************************************************************
 *   Compile using the command: cc -lm -o mcluster mcluster.c              *
 *   or use the Makefile, i.e. type: make mcluster or make mcluster_sse    *
 ***************************************************************************/
/***************************************************************************
 * Main array with stellar parameters:                                     *
 *      star[][0] = mass(t = epoch)                                        *
 *      star[][1-3] = {x, y, z}                                            *
 *      star[][4-6] = {vx, vy, vz}                                         *
 *      star[][7] = mass(t = 0)                                            *
 *      star[][8] = kstar, stellar type in (SSE/BSE)                       *
 *      star[][9] = epoch1, age within evolutionary phase (SSE/BSE)        *
 *      star[][10] = ospin, spin of star (SSE/BSE)                         *
 *      star[][11] = rstar, radius of star (SSE/BSE)                       *
 *      star[][12] = lstar, luminosity (SSE/BSE)                           *
 *      star[][13] = epochstar, age of star (SSE/BSE)                      *
 *      star[][14] = Zstar, metallicity of star                            *
 ***************************************************************************/


#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include<sys/stat.h>
#include<getopt.h>
#include "main.h"

//use OpenMP if not specified otherwise
#ifndef NOOMP
#include<omp.h>
#endif



int main (int argv, char **argc) {

	/*******************
	 * Input variables *
	 *******************/
	
	//Basic physical parameters
	int N = 0;					//Number of stars, Mcl will be set to 0 if specified!
	double Mcl = 1000.0;			//Total mass of the cluster, only used when N is set to 0, necessary for usage of maximum stellar mass relation of Weidner & Kroupa (2006)
	int profile = 0;				//Density profile; =0 Plummer model, =1 King model (based on king0.f by D.C. Heggie), =2 mass segregated Subr model (Subr et al. 2007), =3 2-dimensional EFF template (Elson, Fall & Freeman 1987) or Nuker template (Lauer et al. 1995), =-1 no density gradient
	double W0 = 5.0;				//King's W0 paramter [0.3-12.0]
	double S = 0.0;			    	//Fraction of mass segregation for profile; =0.0 unsegregated, =1.0 completely segregated (take maximally S=0.99 for profile=2)
	double D = 3.0;			    	//Fractal dimension; =3.0 unfractal, =2.6 2/8 fractal, =2.0 4/8 fractal, =1.6 6/8 fractal, (2^D children per parent, following Goodwin & Whitworth 2004);
	double Q = 0.5;					//Initial virial ratio; =0.5 virial equilibrium, <0.5 collapsing, >0.5 expanding
	double Rh = 0.8;				//Half-mass radius [pc], ignored if profile = 3, set =-1 for using Marks & Kroupa (2012) Mcl-Rh relation
	double gamma[3] = {2.0, 0.0, 2.0}; //Power-law slopes of EFF/Nuker templates (outer slope, inner slope, transition); set gamma[1] = 0.0 and gamma[2] = 2.0 for EFF (profile = 3)  
	double a = 1.0;					//Scale radius of EFF/Nuker template (profile = 3) [pc]
	double Rmax = 100.0;			//Cut-off radius for EFF/Nuker template (profile = 3) [pc]
	double tcrit = 100.0;			//Simulation time [N-body units (Myr in Nbody6 custom)]
	int tf = 3;						//Tidal field: =0 no tidal field, =1 Near-field approximation, =2 point-mass galaxy, =3 Central point mass, disk & log-halo
	double RG[3] = {8500.0,0.0,0.0}; //Initial Galactic coordinates of the cluster [pc]
	double VG[3] = {0.0,220.0,0.0};  //Initial velocity of the cluster [km/s]
	
	//Mass function parameters
	int mfunc = 1;					//0 = single mass stars; 1 = use Kroupa (2001) mass function; 2 = use multi power law (based on mufu.c by L.Subr)
	double single_mass = 1.0;		//Stellar mass in case of single-mass cluster 
	double mlow = 0.08;				//Lower mass limit for mfunc = 1 & mfunc = 4
	double mup = 100.0;				//Upper mass limit for mfunc = 1 & mfunc = 4
	double alpha[MAX_AN] = {-1.35, -2.35, -2.7, 0.0, 0.0};		//alpha slopes for mfunc = 2
	double mlim[MAX_MN] = {0.08, 0.5, 4.0, 100.0, 0.0, 0.0};	//mass limits for mfunc = 2
	double alpha_L3 = 2.3;			//alpha slope for mfunc = 4 (L3 mass function, Maschberger 2012)
	double beta_L3 = 1.4;			//beta slope for mfunc = 4
	double mu_L3 = 0.2;				//mu parameter for mfunc = 4
	int weidner = 0;				//Usage of Weidner & Kroupa (2006) relation for most massive star; =0 off, =1 on
	int mloss = 3;					//Stellar evolution; 0 = off, 3 = Eggleton, Tout & Hurley [KZ19]
	int remnant = 1;				//Use random kick velocity and present-day escape velocity to determine retention of compact remnants & evolved binary components (only for SSE/BSE version); =0 off, =1 on
	double epoch = 0.0;			    //Age of the cluster, i.e. star burst has been ... Myr before [e.g. 1000.0, default = 0.0] [needs special compiling and SSE by Hurley, Pols & Tout]
	double Z = 0.02;				//Metallicity [0.0001-0.03, 0.02 = solar]
	double FeH = -1.41;				//Metallicity [Fe/H], only used when Z is set to 0
	int prantzos = 0;				//Usage of Prantzos (2007) relation for the life-times of stars. Set upper mass limit to Lifetime(mup) >= epoch
	
	//Binary parameters
	int nbin = 0;				    //Number of primordial binary systems
	double fbin = 0.0;				//Primordial binary fraction, number of binary systems = 0.5*N*fbin, only used when nbin is set to 0 
	int pairing = 2;				//Pairing of binary components; 0= random pairing, 1= ordered pairing for components with masses M>msort, 2= random but separate pairing for components with masses m>Msort
	double msort = 5.0;				//Stars with masses > msort will be sorted and preferentially paired into binaries if pairing = 1
	int adis = 1;					//Semi-major axis distribution; 0= flat ranging from amin to amax, 1= based on Kroupa (1995) period distribution, 2= based on Duquennoy & Mayor (1991) period distribution
	int OBperiods = 1;				//Use period distribution for massive binaries with M_primary > msort from Sana & Evans (2011) if OBperiods = 1
	double amin = 0.0001;			//Minimum semi-major axis for adis = 0 [pc]
	double amax = 0.01;				//Maximum semi-major axis for adis = 0 [pc]
#ifdef SSE
	int eigen = 0;					//Use Kroupa (1995) eigenevolution for pre-main sequence short-period binaries; =0 off, =1 on [use either eigenevolution or BSE; BSE recommended when using SSE]
	int BSE = 1;					//Apply binary star evolution using BSE (Hurley, Tout & Pols 2002) =0 off, =1 on [use either eigenevolution or BSE; BSE recommended when using SSE]
#else
	int eigen = 1;					//Use Kroupa (1995) eigenevolution for pre-main sequence short-period binaries; =0 off, =1 on [use either eigenevolution or BSE; BSE recommended when using SSE]
	int BSE = 0;					//Apply binary star evolution using BSE (Hurley, Tout & Pols 2002) [needs special compiling and BSE]; =0 off, =1 on [use either eigenevolution or BSE; BSE recommended when using SSE]
#endif
	
	//Gas parameters (only used for Nbody6 input)
	double extmass = 0.0;			//external Plummer (gas) sphere mass [Msun]
	double extrad = 0.0;			//external Plummer (gas) sphere scale factor [pc]
	double extdecay = 0.0;			//decay time for gas expulsion (set 0 for no decay) [Myr] 
	double extstart = 0.0;			//delay time for start of gas expulsion [Myr]
	
	//Code parameters
	int code = 3;					//Nbody version: =0 Nbody6, =1 Nbody4, =2 Nbody6 custom, =3 only create output list of stars, =4 Nbody7 (not yet fully functional)
	unsigned int seed = 0;			//Number seed for random number generator; =0 for randomization by local time
	char *output = "test";   		//Name of output files
	double dtadj = 1.0;				//DTADJ [N-body units (Myr in Nbody6 custom)], energy-check time step
	double dtout = 1.0;				//DELTAT [N-body units (Myr in Nbody6 custom)], output interval, must be multiple of DTADJ
	double dtplot = 100.0;			//DTPLOT [N-body units (Myr in Nbody6 custom)], output of HRdiagnostics, should be multiple of DTOUT, set to zero if output not desired
	int gpu = 0;					//Use of GPU, 0= off, 1= on
	int regupdate = 1;				//Update of regularization parameters during computation; 0 = off, 0 > on
	int etaupdate = 0;				//Update of ETAI & ETAR during computation; 0 = off, 0 > on
	int esc = 2;					//Removal of escapers; 0 = no removal, 1 = regular removal at 2*R_tide; 2 = removal and output in ESC
	int units = 1;				    //Units of McLuster output; 0= Nbody-Units, 1= astrophysical units
	
	//McLuster internal parameters
	int match = 1;					//Make cluster half-mass radius exactly match the expected/desired value; =0 off, =1 on (recommended)
	int symmetry = 1;				//Force spherical symmetry for fractal clusters; =0 off, =1 on (recommended)
	int check = 0;					//Make energy check at end of McLuster; =0 off, =1 on
	int create_radial_profile = 1;	//Creates a radial density profile and prints it to the screen; =0 off, =1 on
	int create_cumulative_profile = 1;	//Creates a radial cumulative profile and prints it to the screen; =0 off, =1 on
	double Rgal = 10000.0;			//Distance of cluster from sun for artificial CMD with observational errors [pc] 
	double Zsun = 0.02;				//Solar metallicity
	int NMAX = 1500000;	     		//Maximum number of stars & orbits allowed in McLuster
	int NNBMAX_NBODY6 = 500;		//Maximum number of neighbours allowed in NBODY6
	double upper_IMF_limit = 150.0; //Maximum stellar mass allowed in McLuster [Msun]
	int an = 0;						//Counter for number of alpha slopes for mfunc = 2
	int mn = 0;						//Counter for number of mass limits for mfunc = 1, 2 & 4
	int xn = 0;						//Counter for components of galactocentric radius vector
	int vn = 0;						//Counter for components of cluster velocity vector 
	int xx = 0;						//Counter for external potential input parameters 
	int gn = 0;						//Counter for EFF/Nuker profile parameters
	double extgas[4];				//Input array for external potential parameters
	
	
	//SSE internal parameters (see Hurley, Pols & Tout 2000) 
	value1_.neta = 0.5;			//Reimers mass-loss coefficent (neta*4x10^-13; 0.5 normally)
	value1_.bwind = 0.0;		//Binary enhanced mass loss parameter (inactive for single)
	value1_.hewind = 1.0;		//Helium star mass loss factor (1.0 normally)
	value1_.mxns = 3.0;			//Maximum NS mass (1.8, nsflag=0; 3.0, nsflag=1)
	points_.pts1 = 0.05;		//Time-step parameter in evolution phase: MS (0.05)
	points_.pts2 = 0.01;		//Time-step parameter in evolution phase: GB, CHeB, AGB, HeGB (0.01)
	points_.pts3 = 0.02;		//Time-step parameter in evolution phase: HG, HeMS (0.02)
	value4_.sigma = 190.0;		//Kick velocities
	value4_.bhflag = 1;			//bhflag > 0 allows velocity kick at BH formation

	//BSE internal parameters (see Hurley, Pols & Tout 2002) 
	flags_.ceflag = 0;			//ceflag > 0 activates spin-energy correction in common-envelope (0) #defunct# 
	flags_.tflag = 1;			//tflag > 0 activates tidal circularisation (1)
	flags_.ifflag = 0;			//ifflag > 0 uses WD IFMR of HPE, 1995, MNRAS, 272, 800 (0)
	flags_.nsflag = 1;			//nsflag > 0 takes NS/BH mass from Belczynski et al. 2002, ApJ, 572, 407 (1)
	flags_.wdflag = 1;			//wdflag > 0 uses modified-Mestel cooling for WDs (0)
	
	value5_.beta = 1.0/8.0;		//beta is wind velocity factor: proportional to vwind**2 (1/8)
	value5_.xi = 1.0;			//xi is the wind accretion efficiency factor (1.0)
	value5_.acc2 = 3.0/2.0;		//acc2 is the Bondi-Hoyle wind accretion factor (3/2)
	value5_.epsnov = 0.001;		//epsnov is the fraction of accreted matter retained in nova eruption (0.001)
	value5_.eddfac = 1.0;		//eddfac is Eddington limit factor for mass transfer (1.0)
	value5_.gamma = -1.0;		//gamma is the angular momentum factor for mass lost during Roche (-1.0)

	value2_.alpha1 = 1.0;		//alpha1 is the common-envelope efficiency parameter (1.0)
	value2_.lambda = 0.5;		//lambda is the binding energy factor for common envelope evolution (0.5)


	//Command line input
	int option;
	while ((option = getopt(argv, argc, "N:M:P:W:R:r:c:g:S:D:T:Q:C:A:O:G:o:f:a:m:B:b:p:s:t:e:Z:X:x:V:u:h:?")) != -1) switch (option)
	{
		case 'N': N = atoi(optarg); Mcl = 0.0; break;
		case 'M': Mcl = atof(optarg); N = 0; break;
		case 'P': profile = atoi(optarg); break;
		case 'W': W0 = atof(optarg); break;
		case 'R': Rh = atof(optarg); break;
		case 'r': a = atof(optarg); break;
		case 'c': Rmax = atof(optarg); break;
		case 'g': 
			if (gn < 3) { gamma[gn] = atof(optarg); gn++; break; }
		case 'S': S = atof(optarg); break;
		case 'D': D = atof(optarg); break;
		case 'T': tcrit = atof(optarg); break;
		case 'Q': Q = atof(optarg); break;
		case 'C': code = atoi(optarg); break;
		case 'A': dtadj = atof(optarg); break;
		case 'O': dtout = atof(optarg); break;
		case 'G': gpu = atoi(optarg); break;
		case 'o': output = optarg; break;
		case 'f': mfunc = atoi(optarg); break;
		case 'a' :
			if (an < MAX_AN) { 
				alpha[an] = atof(optarg);
				if (an == 0) alpha_L3 = atof(optarg);
				if (an == 1) beta_L3 = atof(optarg);
				if (an == 2) mu_L3 = atof(optarg);
				an++; 
				break; 
			} else { printf("\nError: Number of alphas exceeded maximum limit of %d\n", MAX_AN); return 1; }
		case 'm' :
			if (mn < MAX_MN) { 
				mlim[mn] = atof(optarg);
				if (mn == 0) mlow = atof(optarg);
				if (mn == 1) mup = atof(optarg);
				mn++; 
				break;
			} else { printf("\nError: Number of mass params exceded maximum limit of %d\n", MAX_MN); return 1; }
		case 'B': nbin = atoi(optarg); break;
		case 'b': fbin = atof(optarg); break;
		case 'p': pairing = atoi(optarg); break;
		case 's': seed = atoi(optarg); break;
		case 't': tf = atoi(optarg); break;
		case 'e': epoch = atof(optarg); break;
		case 'Z': Z = atof(optarg); break;
		case 'X' :
			if (xn < 3) { RG[xn] = atof(optarg); xn++; break; }
		case 'V' :
			if (vn < 3) { VG[vn] = atof(optarg); vn++; break; }
		case 'x' : 
			if (xx < 4) { extgas[xx] = atof(optarg); xx++; break; }
		case 'u': units = atoi(optarg); break;
		case ':':
		case 'h':	help(msort); return 1;
		case '?':	help(msort); return 1;
	};
	
	//print summary of input parameters to .info file
	info(output, N, Mcl, profile, W0, S, D, Q, Rh, gamma, a, Rmax, tcrit, tf, RG, VG, mfunc, single_mass, mlow, mup, alpha, mlim, alpha_L3, beta_L3, mu_L3, weidner, mloss, remnant, epoch, Z, prantzos, nbin, fbin, pairing, msort, adis, amin, amax, eigen, BSE, extmass, extrad, extdecay, extstart, code, seed, dtadj, dtout, dtplot, gpu, regupdate, etaupdate, esc, units, match, symmetry, OBperiods);
	
	
	/*********
	 * Start *
	 *********/
	
	printf("\n-----START----         \n"); 

#ifdef NOOMP
	clock_t t1, t2;
	t1 = clock();							//start stop-watch
#else
	double t1, t2;
#pragma omp parallel
	{
		if (omp_get_thread_num() == 0) 
			if (omp_get_num_threads() > 1) printf("\n\nUsing OpenMP with %d threads\n", omp_get_num_threads());

		t1 = omp_get_wtime(); //start stop-watch
	}
#endif
  
	if (seed) srand48(seed);				//initialize random number generator by seed
	else {
		seed = (unsigned) time(NULL);	//initialize random number generator by local time
		seed %= 100000;
		srand48(seed);
	}
	printf ("\n\nRandom seed = %i\n\n\n", seed);
	if (seed) value3_.idum = seed;			//idum is the random number seed used in the kick routine. 
	else value3_.idum = 10000000.0*drand48();
	int i,j;
	double M;								//Total mass [M_sun]
	double mmean;							//Mean stellar mass [M_sun]
	int NNBMAX;								//Maximum neighbour number (Nbody6 only)
	double RS0;								//Initial radius of neighbour sphere [pc], Nbody6 only
	double rtide;							//Tidal radius [pc]
	double omega;							//Angular velocity of cluster around the galaxy
	double rvir;							//Virial radius [pc]
	double cmr[7];							//For CoM correction
	double rking, rplummer;					//King-, Plummer radius
	double MMAX;							//most massive star
	double tscale;							//time-scale factor
	double ekin = 0.0;						//kinetic energy
	double epot = 0.0;						//potential energy
	double sigma = 0.0;						//velocity dispersion
	int bin;								//KZ(22) parameter (binaries yes/no) 
	int sse;								//(evolved stellar population yes/no)
	double submass[MAX_AN], subcount[MAX_AN], norm[MAX_AN], N_tmp, M_tmp; //mass function parameters for mfunc = 2	
	double Rh2D, Rh3D;						//actual 2D/3D half-mass radius of the model
	
	if (profile == 3) Rh = a; //set half-mass radius temporarily to scale radius for computation of escape velocity

	if ((Mcl) && (N)) {
		printf("\nWARNING: specify either Mcl (-M) or N (-N)!\n\n");
		exit (1);
	} else if ((!Mcl) && (mfunc == 3)) {
		printf("\nWARNING: specify Mcl (-M) when using optimal sampling (-f 3)!\n\n");
		exit (1);
	}
    
	if (xx > 0) {
		if ((xx == 4) && (extgas[2])) {
			extmass = extgas[0];
			extrad = extgas[1];
			extdecay = extgas[2];
			extstart = extgas[3];
		} else {
			printf("\nWARNING: Insufficient or incorrect parameters specified for external Plummer potential!\n\n");
			exit (1);
		}
	}
		
	
	/***********************
	 * Generate star array *
	 ***********************/

	int columns = 15;
	double **star;
	star = (double **)calloc(NMAX,sizeof(double *));
	for (j=0;j<NMAX;j++){
		star[j] = (double *)calloc(columns,sizeof(double));
		if (star[j] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}

	
	
	/*******************************************	
	 * Evaluate Z from [Fe/H] if Z is set to 0 *
	 *******************************************/
	
	 if (!Z) {
		Z = pow(10.0, 0.977*FeH)*Zsun; //Bertelli, Bressan, Chiosi, Fagotto, Nasi, 1994, A&AS, 106, 275
		printf("\nUsing Bertelli et al. (1994) relation to convert FeH = %.3f into Z = %.3f\n", FeH, Z);
	}
	
	
	
	/**********************************
	 * Calculate maximum stellar mass *
	 **********************************/

	 if (mfunc == 3) weidner = 1; //always use Weidner relation when using optimal sampling!

	if (!N && weidner && mfunc) {
		mup = upper_IMF_limit;  
	
		printf("\nUsing maximum stellar mass-cluster mass relation for upper stellar mass limit\n");
		printf("\n(Weidner & Kroupa 2006, Pflamm-Altenburg & Kroupa 2007)\n");

		//analytic fit to the observational data from Pflamm-Altenburg & Kroupa (2007), implementation by M. Kruckow
		MMAX = pow(10,2.56*log10(Mcl)*pow(pow(3.82,9.17)+pow(log10(Mcl),9.17),-1/9.17)-0.38);

		/*	
		//three-part power-law fit to the observational data
		if (Mcl < 1000.0) {
			MMAX = (log10(Mcl)*0.540563175) - 0.14120167;
			MMAX = pow(10.0,MMAX);
		} else if (Mcl < 3300.0) {
			MMAX = (log10(Mcl)*0.19186051) + 0.9058611;
			MMAX = pow(10.0,MMAX);
		} else {
			MMAX = (log10(Mcl)*0.360268003) + 0.313342031;
			MMAX = pow(10.0,MMAX);
		}*/
	} else {
		MMAX = mup;
	}	
	
	if (mfunc && epoch && prantzos) {
		printf("\nUsing Prantzos (2007) relation to reduce upper mass limit to Lifetime(mup) > epoch\n");
		while (Lifetime(MMAX) < sqrt(pow(epoch,2))) {
			MMAX -= 0.01;
		}
	}
	
	
	/*******************
	 * Generate masses *
	 *******************/
	
	printf("\n\n-----GENERATE MASSES-----     \n"); 

	//This loop has to be called multiple times with different N (or Mcl), Z and epoch in order to create multiple stellar populations
	//make sure that epoch and Z are stored together with the stars, then concatenate all arrays
	//for now, all Z and epochs are the same, i.e. a single stellar population
	if (mfunc == 1) {
		printf("\nMaximum stellar mass set to: %.2f\n", MMAX);
		generate_m1(&N, star, mlow, mup, &M, &mmean, MMAX, Mcl, epoch, Z, Rh, remnant);
	} else if (mfunc == 2) {
		if (mn) {
			for (i = mn+1; i < MAX_MN; i++) mlim[i] = 0.0;
		} else {
			for (i=0; i<MAX_MN; i++) {
				if (mlim[i]) mn++;
			}
		}
		if (an) {
			for (i = an+1; i < MAX_AN; i++) alpha[i] = 0.0;		
		} else {
			for (i=0; i<MAX_AN; i++) {
				if (alpha[i]) an++;
			}
		}
		if (an >= mn) an = mn - 1;
		mn = an + 1;
		if (!mn){
			printf("\nError: at least one mass limit has to be specified\n");
			return 1;
		} 	else if (mn == 1) {
			single_mass = mlim[0];
			printf("\nSetting stellar masses to %g solar mass\n",single_mass);
			if (!N) N = Mcl/single_mass;
			for (j=0;j<N;j++) star[j][0] = 1.0/N;
			mmean = single_mass;
			M = N*mmean;
			printf("\nM = %g\n", M);
		}	else {
			printf("\nMaximum stellar mass set to: %.2f\n",MMAX);
			norm[an-1] = 1.; //normalization factor of integral
			N_tmp = subcount[an-1] = subint(mlim[an-1], mlim[an], alpha[an-1] + 1.); //integrated number of stars in interval [mlim[an-1]:mlim[an]]
			M_tmp = submass[an-1] = subint(mlim[an-1], mlim[an], alpha[an-1] + 2.); //integrated mass of stars in interval [mlim[an-1]:mlim[an]]
			for (i = an - 2; i >= 0; i--) {
				norm[i] = norm[i+1] * pow(mlim[i+1], alpha[i+1] - alpha[i]);
				subcount[i] = norm[i] * subint(mlim[i], mlim[i+1], alpha[i] + 1.);
				N_tmp += subcount[i];
				submass[i] = norm[i] * subint(mlim[i], mlim[i+1], alpha[i] + 2.);
				M_tmp += submass[i];
			}
			generate_m2(an, mlim, alpha, Mcl, M_tmp, subcount, &N, &mmean, &M, star, MMAX, epoch, Z, Rh, remnant);
		}
 	} else if (mfunc == 3) {
		printf("\nMaximum stellar mass set to: %.2f\n",MMAX);
		generate_m3(&N, star, mlow, mup, &M, &mmean, MMAX, Mcl);
		randomize(star, N);
	} else if (mfunc == 4) {
		printf("\nMaximum stellar mass set to: %.2f\n", MMAX);
		printf("\nUsing L3 IMF (Maschberger 2012)\n");
		generate_m4(&N, star, alpha_L3, beta_L3, mu_L3, mlow, mup, &M, &mmean, MMAX, Mcl, epoch, Z, Rh, remnant);		
	} else {
		printf("\nSetting stellar masses to %.1f solar mass\n", single_mass);
		if (!N) N = Mcl/single_mass;
		for (j=0;j<N;j++) {
			star[j][0] = single_mass;
			star[j][7] = single_mass;
			star[j][8] = 0;
			star[j][9] = 0.0;
			star[j][10] = 0.0;
			star[j][11] = 0.0;
			star[j][12] = 0.0;
			star[j][13] = 0.0;
			star[j][14] = 0.0;
		}
		mmean = single_mass;
		M = N*mmean;
		printf("\nM = %g\n", M);
		mloss = 0;
	}
	
	//set all stars to the same metallicity and age for now
	double epochstar, zstar;
	epochstar = 0.0; //age compared to the oldest stars in the cluster [Myr]
	zstar = Z;
	for (j=0;j<N;j++) {
		star[j][13] = epochstar;
		star[j][14] = zstar;
	}
	
	//Pair binary masses and convert to centre-of-mass particles
	int Nstars;
	if (!nbin) nbin = 0.5*N*fbin;

	double **mbin;	//component mass & stellar evol parameter array
	mbin = (double **)calloc(nbin,sizeof(double *));
	for (j=0;j<nbin;j++){
		mbin[j] = (double *)calloc(20,sizeof(double));
		if (mbin[j] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}	
	
	Nstars = N;

	if (nbin) {
		printf("\nPreparing binary components.\n");
		//Specify component pairing
		if (pairing) {
			order(star, N, M, msort, pairing);
			if (pairing == 1) printf("\nApplying ordered pairing for stars with masses > %.1f Msun.\n",msort);
			else if (pairing == 2) printf("\nApplying random pairing for stars with masses > %.1f Msun.\n",msort);
		} else {
			randomize(star, N);
			printf("\nApplying random pairing.\n");
		}
		N -= nbin;
		for (j=0;j<nbin;j++) {
			mbin[j][0] = star[2*j][0]+star[2*j+1][0]; //system mass
			mbin[j][1] = star[2*j][0]; //primary mass
			mbin[j][2] = star[2*j+1][0]; //secondary mass
			mbin[j][3] = star[2*j][7]; //primary m0
			mbin[j][4] = star[2*j][8]; //primary kw
			mbin[j][5] = star[2*j][9]; //primary epoch1
			mbin[j][6] = star[2*j][10]; //primary spin
			mbin[j][7] = star[2*j][11]; //primary r
			mbin[j][8] = star[2*j][12]; //primary lum
			mbin[j][9] = star[2*j+1][7]; //secondary m0
			mbin[j][10] = star[2*j+1][8]; //secondary kw
			mbin[j][11] = star[2*j+1][9]; //secondary epoch1
			mbin[j][12] = star[2*j+1][10]; //secondary spin
			mbin[j][13] = star[2*j+1][11]; //secondary r
			mbin[j][14] = star[2*j+1][12]; //secondary lum
			mbin[j][15] = 1000+j; //identifier
			mbin[j][16] = star[2*j][13]; //primary epochstar
			mbin[j][17] = star[2*j+1][13]; //secondary epochstar
			mbin[j][18] = star[2*j][14]; //primary zstar
			mbin[j][19] = star[2*j+1][14]; //secondary zstar

			star[2*j][0] += star[2*j+1][0]; //system mass
			star[2*j+1][0] = 0.0;
			star[2*j][7] = 1000+j; //identifier
			star[2*j+1][7] = 0.0; //identifier
			star[2*j][12] += star[2*j+1][12]; //system luminosity
			star[2*j+1][12] = 0.0;
		}
		order(star, Nstars, M, 0.0, 0);
		randomize(star, N);
	}


	//prepare mass segregation
	double mlowest = MMAX; //search lowest mass star
	double mhighest = 0; //search highest mass star
	double mmeancom = 0.0;
	for (i=0;i<N;i++) {
		star[i][0] /= M; //scale masses to Nbody units
		mmeancom += star[i][0];
		if (star[i][0] < mlowest) mlowest = star[i][0];
		if (star[i][0] > mhighest) mhighest = star[i][0];
	}
	mmeancom /= N;
	int Nseg = ceil(N*mmeancom/mlowest);  //number of necessary pos & vel pairs for Baumgardt et al. (2008) mass segregation routine
	int Nunseg = N;

	double *Mcum;
	Mcum = (double *)calloc(N,sizeof(double));
		
	if ((S) && !(profile == 2)) {//sort masses when mass segregation parameter > 0
		printf("\nApplying mass segregation with S = %f\n",S);
		segregate(star, N, S);
		for (i=0;i<N;i++) {//calculate cumulative mass function Mcum
			Mcum[i] = 0.0;
			for (j=0;j<=i;j++) Mcum[i] = Mcum[i] + star[j][0];
		}
		N = Nseg;
	}
	
	
	/*************************************
	 * Generate positions and velocities *
	 *************************************/
	
	printf("\n\n-----GENERATE POSITIONS & VELOCITIES-----   \n"); 

	//calculate half-mass radius according to Marks & Kroupa 2012 if Rh is set to -1
	if ((Rh == -1) && (Mcl)) {
		Rh = 0.1*pow(Mcl,0.13); //Marks & Kroupa (2012), implementation by M. Kruckow
		printf("\nUsing Marks & Kroupa (2012) relation to derive half-mass radius from cluster mass: %g (pc)\n", Rh);
	}
	
	
	//evaluate approximate tidal radius assuming circular orbit
	if (tf == 3) {
		//in the case of log-halo potential, assume kappa = 1.4omega (eq. 9 in Kuepper et al. 2010)
		omega = sqrt(VG[0]*VG[0]+VG[1]*VG[1]+VG[2]*VG[2])/sqrt(RG[0]*RG[0]+RG[1]*RG[1]+RG[2]*RG[2]);
		rtide = pow(G*M/(2.0*omega*omega),1.0/3.0);
	} else if (!tf) {
		rtide = 1.0E5;
	} else if ((tf == 1) && (code == 0)) {
		//in case of Sverre's Nbody6 standard tidal field
		rtide = pow(1.0*M/(3.0*M1pointmass),1.0/3.0)*8000.0;
	} else { 
		//in the case of a point mass potential or near field approximation
		rtide = pow(1.0*M/(3.0*M1pointmass),1.0/3.0)*sqrt(RG[0]*RG[0]+RG[1]*RG[1]+RG[2]*RG[2]);
	}
	printf("\nApproximate tidal radius: %g (pc)\n", rtide);

	
	//generate scaled pos & vel, postpone scaling for Plummer and King in case of mass segregation
	double rhtemp, rvirtemp;
	if (profile == 1) {
		printf("\nGenerating King model with parameters: N = %i\t W0 = %g\t Rh = %.3f\t D = %.2f\n",N, W0, Rh, D);
		generate_king(N, W0, star, &rvirtemp, &rhtemp, &rking, D, symmetry);
	} else if (profile == 2) {
		N = Nunseg;
		printf("\nGenerating segregated Subr model with parameters: N = %i\t S = %g\t Rh = %.3f\n",N, S, Rh);
		rvir = Rh/0.76857063065978; //(value provided by L. Subr) 
		generate_subr(N, S, star, rtide, rvir);
		printf ("\nrvir = %.5f\t rh = %.5f\t rtide = %.5f (pc)\n", rvir, Rh, rtide);
	} else if (profile == 3) {
		if (gamma[1] == 0.0 && gamma[2] == 2.0) printf("\nGenerating EFF model with parameters: N = %i\t a = %.1f\t gamma = %.3f\t Rmax = %.1f\t D = %.2f\n",N, a, gamma[0], Rmax, D);
		else printf("\nGenerating Nuker model with parameters: N = %i\t a = %.1f\t gamma (outer) = %.3f\t gamma (inner) = %.3f\t transition = %.3f\t Rmax = %.1f\t D = %.2f\n",N, a, gamma[0],gamma[1],gamma[2], Rmax, D);
		double p[6];
		p[1] = 1.0; //rho0, will be scaled according to Rmax and Mtot
		p[2] = a;//scale radius
		p[3] = gamma[0];//outer power-law slope (>0.5)
		p[4] = gamma[1];//inner power-law slope (0.0 for EFF template)
		p[5] = gamma[2];//transition parameter (2.0 for EFF template)
		generate_profile(N, star, Rmax, M, p, &Rh, D, symmetry);
		printf("\nRh = %.1f pc\n", Rh);
	} else if (profile == -1) {
		printf("\nGenerating fractal distribution with parameters: N = %i\t Rh = %.3f\t D = %.2f\n", N, Rh, D);
		fractalize(D, N, star, 0, symmetry);
		rvir = Rh;
	} else {
		printf("\nGenerating Plummer model with parameters: N = %i\t Rh = %.3f\t D = %.2f\n", N, Rh, D);
		rvir = Rh/0.772764;
		rplummer = Rh/1.305;
		generate_plummer(N, star, rtide, rvir, D, symmetry);
	}
	

	//Apply Baumgardt et al. (2008) mass segregation
	if (!(profile == 2) && (S)) {
		double **m_temp;
		m_temp = (double **)calloc(Nunseg,sizeof(double *));
		for (j=0;j<Nunseg;j++){
			m_temp[j] = (double *)calloc(9,sizeof(double));
			if (m_temp[j] == NULL) {
				printf("\nMemory allocation failed!\n");
				return 0;
			}
		}	
				
		for (i=0;i<Nunseg;i++) {//temporarily store the masses & stellar evol parameters
			m_temp[i][0] = star[i][0];
			m_temp[i][1] = star[i][7];
			m_temp[i][2] = star[i][8];
			m_temp[i][3] = star[i][9];
			m_temp[i][4] = star[i][10];
			m_temp[i][5] = star[i][11];
			m_temp[i][6] = star[i][12];
			m_temp[i][7] = star[i][13];
			m_temp[i][8] = star[i][14];
		}
		
		printf("\nOrdering orbits by energy.\n");
		energy_order(star, N, Nstars);

		int nlow, nhigh, nrandom;
		for (i=0;i<Nunseg;i++) {
			nhigh = Nseg*Mcum[i];
			if (i) {
				nlow = Nseg*Mcum[i-1];
			} else {
				nlow = 0;
			}
			nrandom = (nhigh-nlow)*drand48()+nlow;
			star[i][0] = m_temp[i][0];
			star[i][1] = star[nrandom][1];
			star[i][2] = star[nrandom][2];
			star[i][3] = star[nrandom][3];
			star[i][4] = star[nrandom][4];
			star[i][5] = star[nrandom][5];
			star[i][6] = star[nrandom][6];
			star[i][7] = m_temp[i][1];
			star[i][8] = m_temp[i][2];
			star[i][9] = m_temp[i][3];
			star[i][10] = m_temp[i][4];
			star[i][11] = m_temp[i][5];
			star[i][12] = m_temp[i][6];
			star[i][13] = m_temp[i][7];
			star[i][14] = m_temp[i][8];
		}

		
		for (j=0;j<Nunseg;j++) free (m_temp[j]);
		free(m_temp);		
		
		N = Nunseg;
	}
		
	
	//CoM correction
	printf("\nApplying centre-of-mass correction.\n");
	for (j=0; j<7; j++) cmr[j] = 0.0;
	
	for (j=0; j<N; j++) {
		for (i=1;i<7;i++) 
			cmr[i] += star[j][0]*star[j][i]; 
	} 
	
	for (j=0; j<N; j++) {
		for (i=1;i<7;i++)
			star[j][i] -= cmr[i];
	}
	
	
	//apply scaling to Nbody-units
	if (profile == 0) {
		printf("\nRe-scaling of orbits (dt ~ N^2!)\n");
		double ke = 0.0;
		double pe = 0.0;
		double sx, sv, r2;			
#ifndef NOOMP	
#pragma omp parallel shared(N, star)  private(i, j, r2)
		{
#pragma omp for reduction(+: pe, ke) schedule(dynamic)
#endif
		for (i=0;i<N;i++) {
			if (i) {
				for (j=0;j<i-1;j++) {
					r2 = (star[i][1]-star[j][1])*(star[i][1]-star[j][1]) + (star[i][2]-star[j][2])*(star[i][2]-star[j][2]) +(star[i][3]-star[j][3])*(star[i][3]-star[j][3]) ;
					pe -=  star[i][0]*star[j][0]/sqrt(r2);
				}
			}
			ke += star[i][0]*(pow(star[i][4],2)+pow(star[i][5],2)+pow(star[i][6],2));
		}
#ifndef NOOMP
		}
#endif
		rvir = -GNBODY*pow(MNBODY,2)/(2.0*pe);
		sx = 1.0/rvir;
		ke *= 0.5;
		sv = sqrt(4.0*ke);
		for (i=0;i<N;i++) {
			star[i][1] *= sx;
			star[i][2] *= sx;
			star[i][3] *= sx;
			star[i][4] /= sv;
			star[i][5] /= sv;
			star[i][6] /= sv;
		}

		ke = 0;
		for (i=0;i<N;i++) {
			ke += M*star[i][0]*(pow(star[i][4],2)+pow(star[i][5],2)+pow(star[i][6],2));
		}
		ke /= 0.5*N;
		printf("Dynamical temperature of centre-of-mass particles kT = %lf\n\n",ke);
		
		//make half-mass radius of the system match the desired one
		radial_profile(star, N, rvir, M, 0, 0, code, &NNBMAX, &RS0, &Rh2D, &Rh3D, NNBMAX_NBODY6);
		if (match) {
			printf("\nmeasuring half-mass radius: %.7f \t %.7f (should/is)\nand correcting for this factor\n",Rh, Rh3D);
			rvir = rvir *Rh/Rh3D;
		}
		printf ("\nrvir = %.5f\t rh = %.5f\t rplummer = %.5f\t rtide = %.5f (pc)\n", rvir, Rh, rplummer, rtide);
	} else if ((profile == 1) || (profile == 3) || (profile == -1)) {
		printf("\nRe-scaling of orbits (dt ~ N^2!)\n");
		double pe = 0.0;
		double ke = 0.0;
		double r2, vscale;
#ifndef NOOMP
#pragma omp parallel shared(N, star)  private(i, j, r2)
		{
#pragma omp for reduction(+: pe, ke) schedule(dynamic)
#endif
		for (i=0;i<N;i++) {
			for (j=0;j<i-1;j++) {
				r2 = (star[i][1]-star[j][1])*(star[i][1]-star[j][1]) + (star[i][2]-star[j][2])*(star[i][2]-star[j][2]) +(star[i][3]-star[j][3])*(star[i][3]-star[j][3]) ;
				pe -=  star[i][0]*star[j][0]/sqrt(r2);
			}
			ke += star[i][0]*(pow(star[i][4],2)+pow(star[i][5],2)+pow(star[i][6],2));
		}
#ifndef NOOMP
		}		
#endif
		ke *= 0.5;
		rvir = -GNBODY*pow(MNBODY,2)/(2.0*pe);
		vscale = sqrt(4.0*ke);
			
		ke = 0;
		for (i=0;i<N;i++) {
			for (j=0;j<3;j++) {
				star[i][j+1] /= rvir;
				star[i][j+4] /= vscale;
			}
			ke += M*star[i][0]*(pow(star[i][4],2)+pow(star[i][5],2)+pow(star[i][6],2));
		}
		ke /= 0.5*N;
		printf("Dynamical temperature of centre-of-mass particles kT = %lf\n\n",ke);
		
		//make half-mass radius of the system match the desired one
		radial_profile(star, N, rvir, M, 0, 0, code, &NNBMAX, &RS0, &Rh2D, &Rh3D, NNBMAX_NBODY6);
		if (match) {
			printf("\nmeasuring half-mass radius: %.7f \t %.7f (should/is)\nand correcting for this factor\n",Rh, Rh3D);
			rvir = rvir *Rh/Rh3D;
		}
		//printf ("\nrvir = %.5f\t rh = %.5f\t rtide = %.5f (pc)\n", rvir, Rh, rtide);
		//printf("Edge radius (King units) = %g\t(Nbody units) = %g\n", rking, rking/rvirtemp);
		//printf("Core radius (King units) = %g\t(Nbody units) = %g\n\n", 1.0, 1.0/rvirtemp);
		//printf("Concentration = %g\n", log10(rking));		
	} else if (profile == 2) {
		double ke = 0;
#ifndef NOOMP
#pragma omp parallel shared(N, star)  private(i)
		{
#pragma omp for reduction(+: ke) schedule(dynamic)
#endif
		for (i=0;i<N;i++) {
			ke += M*star[i][0]*(pow(star[i][4],2)+pow(star[i][5],2)+pow(star[i][6],2));
		}
#ifndef NOOMP
		}
#endif
		ke /= 0.5*N;
		printf("Dynamical temperature of centre-of-mass particles kT = %lf\n\n",ke);
		//make half-mass radius of the system match the desired one
		radial_profile(star, N, rvir, M, 0, 0, code, &NNBMAX, &RS0, &Rh2D, &Rh3D, NNBMAX_NBODY6);
		printf("\nmeasuring half-mass radius: %.7f \t %.7f (should/is)\nand correcting for this factor\n",Rh, Rh3D);
		rvir = rvir *Rh/Rh3D;		
	}
		
	
	
	//Calculate radial density profile, estimate NNBMAX and RS0 (important for Nbody6 only)
	radial_profile(star, N, rvir, M, create_radial_profile, create_cumulative_profile, code, &NNBMAX, &RS0, &Rh2D, &Rh3D, NNBMAX_NBODY6);
	printf("\nActual half-mass radius of the cluster = (%.4f / %.4f) pc (3D / 2D)\n", Rh3D, Rh2D);

	//scale RS0 to nbody units for Nbody6
	RS0 /= 1.0*rvir;
	
	
	
	/*********************
	 * Generate Binaries *
	 *********************/
	
	printf("\n\n-----GENERATE BINARIES-----   \n"); 
	
	if ((!fbin) && (!nbin)) {
		printf("\nNo primordial binaries!\n");
	} else {
		//re-create original array with Nstars (original N) entries
		columns = 15;
		double **star_temp;
		star_temp = (double **)calloc(Nstars,sizeof(double *));
		for (j=0;j<Nstars;j++){
			star_temp[j] = (double *)calloc(columns,sizeof(double));
			if (star_temp[j] == NULL) {
				printf("\nMemory allocation failed!\n");
				return 0;
			}
		}

		double **mbin_index; //sort mbin by identifier
		mbin_index = (double **)calloc(nbin,sizeof(double *));
		for (j=0;j<nbin;j++){
			mbin_index[j] = (double *)calloc(2,sizeof(double));
			if (mbin_index[j] == NULL) {
				printf("\nMemory allocation failed!\n");
				return 0;
			}
		}
		for (j=0;j<nbin;j++) {
			mbin_index[j][0] = mbin[j][15];
			mbin_index[j][1] = j;
		}
		shellsort(mbin_index,nbin,2);
		
		double **star_index;//sort star by identifier
		star_index = (double **)calloc(N,sizeof(double *));
		for (j=0;j<N;j++){
			star_index[j] = (double *)calloc(2,sizeof(double));
			if (star_index[j] == NULL) {
				printf("\nMemory allocation failed!\n");
				return 0;
			}
		}
		for (j=0;j<N;j++) {
			star_index[j][0] = star[j][7];
			star_index[j][1] = j;
		}
		shellsort(star_index,N,2);

		int p;
		for (j=0;j<nbin;j++) {
			if (mbin[(int) mbin_index[j][1]][15] == star[(int) star_index[j][1]][7]) {
				star_temp[2*j][0] = mbin[(int) mbin_index[j][1]][1]/(1.0*M);//primary mass
				star_temp[2*j+1][0] = mbin[(int) mbin_index[j][1]][2]/(1.0*M);//secondary mass
				star_temp[2*j][7] = mbin[(int) mbin_index[j][1]][3];//primary m0
				star_temp[2*j+1][7] = mbin[(int) mbin_index[j][1]][9];//secondary m0
				star_temp[2*j][8] = mbin[(int) mbin_index[j][1]][4];//primary kw
				star_temp[2*j+1][8] = mbin[(int) mbin_index[j][1]][10];//secondary kw
				star_temp[2*j][9] = mbin[(int) mbin_index[j][1]][5];//primary epoch
				star_temp[2*j+1][9] = mbin[(int) mbin_index[j][1]][11];//secondary epoch
				star_temp[2*j][10] = mbin[(int) mbin_index[j][1]][6];//primary spin
				star_temp[2*j+1][10] = mbin[(int) mbin_index[j][1]][12];//secondary spin
				star_temp[2*j][11] = mbin[(int) mbin_index[j][1]][7];//primary r
				star_temp[2*j+1][11] = mbin[(int) mbin_index[j][1]][13];//secondary r
				star_temp[2*j][12] = mbin[(int) mbin_index[j][1]][8];//primary lum
				star_temp[2*j+1][12] = mbin[(int) mbin_index[j][1]][14];//secondary lum				
				star_temp[2*j][13] = mbin[(int) mbin_index[j][1]][16];//primary epochstar
				star_temp[2*j+1][13] = mbin[(int) mbin_index[j][1]][17];//secondary epochstar				
				star_temp[2*j][14] = mbin[(int) mbin_index[j][1]][18];//primary Zstar
				star_temp[2*j+1][14] = mbin[(int) mbin_index[j][1]][19];//secondary Zstar				
				
				for (p=1;p<7;p++) {
					star_temp[2*j][p] = star[(int) star_index[j][1]][p];
					star_temp[2*j+1][p] = star[(int) star_index[j][1]][p];
				}
			}
		}
		
		for (j=nbin;j<N;j++) {			
			for (p=0;p<columns;p++) {
				star_temp[j+nbin][p] = star[(int) star_index[j][1]][p];
			}
		}

		N += nbin;

		for (j=0;j<N;j++) 
			for (p=0;p<columns;p++) star[j][p] = star_temp[j][p]; //copy temporary array back to original
		
		for (j=0;j<Nstars;j++) free (star_temp[j]);
		free(star_temp);

		printf("\nCreating %i primordial binary systems, fraction: %6.2f percent.\n", nbin, 2.0*nbin/N*100.0);
		if (seed) srand48(seed);
		get_binaries(nbin, star, M, rvir, pairing, &N, adis, amin, amax, Rh, eigen, BSE, epoch, Z, remnant, OBperiods, msort);

	} 
	
	//Specify KZ(22) & the sse parameter
#ifdef SSE
	if (epoch) sse = 1; //If feeding an evolved stellar population to Nbody6, KZ(12) has to be =2 in order to read-in fort.12
	else sse = 0;
#else
	sse = 0;
#endif
	bin = 4; //KZ(22)
	
	
	/***********
	 * Scaling * 
	 ***********/

	printf("\n\n-----SCALING-----      \n"); 
		
	//scale masses, pos & vel to astrophysical units or Nbody units
	tscale = sqrt(rvir*rvir*rvir/(G*M));

	if (units) {		
		printf("\nScaling to astrophysical units.\n");
		for (j=0; j<N; j++) star[j][0] *= M;

		for (j=0; j<N; j++) {
			for (i=1;i<4;i++)
				star[j][i] *= rvir;
		}
	
		for (j=0; j<N; j++) {
			for (i=4;i<7;i++)
				star[j][i] *= rvir/tscale;
		}
		bin = -1; //KZ(22)
	} else {
		printf("\nScaling to Nbody units.\n");
	}

	//scale mass, radius and decay time of external (gas) potential to Nbody units
	if (extmass) extmass /= M;
	if (extrad) extrad /= rvir;
	if (extdecay) extdecay = 1.0/(extdecay/tscale);
	if (extstart) extstart = extstart/tscale;

	
	
	/**********
	 * Output * 
	 **********/
	
	printf("\n\n-----OUTPUT-----      \n"); 

	if (code == 0) 
		output0(output, N, NNBMAX, RS0, dtadj, dtout, tcrit, rvir, mmean, tf, regupdate, etaupdate, mloss, bin, esc, M, mlow, mup, MMAX, epoch, dtplot, Z, nbin, Q, RG, VG, rtide, gpu, star, sse, seed, extmass, extrad, extdecay, extstart);
	else if (code == 1)
		output1(output, N, dtadj, dtout, tcrit, rvir, mmean, tf, regupdate, etaupdate, mloss, bin, esc, M, mlow, mup, MMAX, epoch, Z, nbin, Q, RG, VG, rtide, gpu, star);
	else if (code == 2)
		output2(output, N, NNBMAX, RS0, dtadj, dtout, tcrit, rvir, mmean, tf, regupdate, etaupdate, mloss, bin, esc, M, mlow, mup, MMAX, epoch, dtplot, Z, nbin, Q, RG, VG, rtide, gpu, star, sse, seed, extmass, extrad, extdecay, extstart);
	else if (code == 3)
		output3(output, N, rvir, Rh, mmean, M, epoch, Z, RG, VG, rtide, star, Rgal);
	else if (code == 4) 
		output4(output, N, NNBMAX, RS0, dtadj, dtout, tcrit, rvir, mmean, tf, regupdate, etaupdate, mloss, bin, esc, M, mlow, mup, MMAX, epoch, dtplot, Z, nbin, Q, RG, VG, rtide, gpu, star, sse, seed, extmass, extrad, extdecay, extstart);
	
	
	
	
	
	/**********************
	 * Final energy check *
	 **********************/
	
	printf("\n\n-----FINISH-----  \n"); 

	if (check) {
		printf("\nMaking final energy check... (may take a while but can be aborted by pressing CTRL+c)\n");

#ifndef NOOMP
#pragma omp parallel shared(N, star)  private(i, j)
		{
#pragma omp for reduction(+: ekin, epot, sigma) schedule(dynamic)
#endif
			for (j=0; j<N; j++) {
			ekin += star[j][0]*((star[j][4]*star[j][4])+(star[j][5]*star[j][5])+(star[j][6]*star[j][6]));
			if (j) {
				for (i=0;i<j-1;i++) 
					epot -= star[i][0]*star[j][0]/sqrt((star[i][1]-star[j][1])*(star[i][1]-star[j][1])+(star[i][2]-star[j][2])*(star[i][2]-star[j][2])+(star[i][3]-star[j][3])*(star[i][3]-star[j][3])); 
			}
			sigma += star[j][4]*star[j][4]+star[j][5]*star[j][5]+star[j][6]*star[j][6];
		} 
#ifndef NOOMP
		}
#endif
		if (units) epot *= G;
		ekin *= 0.5;

		sigma = sqrt(sigma/N);
		tscale = sqrt(rvir*rvir*rvir/(G*M));
		
		printf("\nEkin = %g\t Epot = %g\t Etot = %g \t kT = %g", ekin, epot, ekin+epot, ekin/(N-nbin));
		printf("\nVel.Disp. = %g\tCross.Time = %g \n", sigma, 2.0/sigma);

		if (units) printf("Vel.Disp. = %g\tCross.Time = %g (Nbody units)\n", sigma/rvir*tscale, 2.0/sigma/tscale);
		else printf("Vel.Disp. = %g\tCross.Time = %g (physical units, km/s, Myr)\n", sigma*rvir/tscale, 2.0/sigma*tscale);
	}
	
	free(Mcum);
	
	for (j=0;j<nbin;j++) free (mbin[j]);
	free(mbin);
	
	for (j=0;j<NMAX;j++) free (star[j]);
	free(star);
	
#ifdef NOOMP
	t2 = clock();														//stop stop-watch
	printf("\nElapsed time: %g sec\n",(double)(t2-t1)/CLOCKS_PER_SEC);	//print stopped time	
#else
#pragma omp parallel
	{
		t2 = omp_get_wtime();//stop stop-watch
	}
	printf("\nElapsed time: %g sec\n",t2-t1);	//print stopped time
#endif
	
	return 0;
}




/*************
 * Functions *
 *************/

int generate_m1(int *N, double **star, double mlow, double mup, double *M, double *mmean, double MMAX, double Mcl, double epoch, double Z, double Rh, int remnant) {
	int ty, i;
	double alpha1, alpha2, c1, c2, k1, k2, xx, mth;

	//set up parameters and variables for SSE (Hurley, Pols & Tout 2002)
	int kw;          //stellar type
	double mass;  //initial mass
	double mt;    //actual mass
	double r = 0.0;     //radius
	double lum = 0.0;   //luminosity
	double mc = 0.0;    //core mass
	double rc = 0.0;    //core radius
	double menv = 0.0;  //envelope mass
	double renv = 0.0;  //envelope radius
	double ospin;		//spin
	double epoch1;    //time spent in current evolutionary state
	double tms = 0.0;   //main-sequence lifetime
	double tphys;       //initial age
	double tphysf = epoch;//final age
	double dtp = epoch+1; //data store value, if dtp>tphys no data will be stored
	double z = Z;      //metallicity
	double zpars[20];     //metallicity parameters
	double vkick;	 //kick velocity for compact remnants
	double vesc;	//escape velocity of cluster
	int lostremnants = 0; //number of ejected compact remnants
	double lostremnantsmass = 0.0; //mass of ejected compact remnants
	if (Mcl && remnant && (Rh != -1)) {
		vesc = sqrt(2.0*G*Mcl/Rh);
		printf("Escape velocity of cluster = %.4f km/s\n", vesc);
	} else if ((!remnant) || (Rh == -1)) {
		vesc = 1.0E10;
		printf("Keeping all compact remnants\n");
	} else {
		vesc = sqrt(2.0*0.4**N/Rh);
		printf("Estimated escape velocity of cluster assuming mean stellar mass of 0.4 Msun = %.4f km/s\n", vesc);
	}
	for (i=0; i<20; i++) zpars[i] = 0;
	zcnsts_(&z,zpars);  //get metallicity parameters

	printf("\nSetting up stellar population with Z = %.4f.\n",Z);
	
	if (epoch) printf("\nEvolving stellar population for %.1f Myr.\n",epoch);
	
	//set up mass function parameters
	ty = 2;   
	alpha1 = 1.3;
	alpha2 = 2.3;
	
	c1 = 1.0-alpha1;  
	c2 = 1.0-alpha2; 
	
	k1 = 2.0/c1*(pow(0.5,c1)-pow(mlow,c1)); 
	if (mlow>0.5) {
        k1 = 0;
        k2 = 1.0/c2*(pow(mup,c2)-pow(mlow,c2)); 
	} else 
        k2 = k1 + 1.0/c2*(pow(mup,c2)-pow(0.5,c2));
	if (mup<0.5) {
		k1 = 2.0/c1*(pow(mup,c1)-pow(mlow,c1));
		k2 = k1;
	}
	

	//determine theoretical mean mass from mass function
	c1 = 2.0-alpha1;
	c2 = 2.0-alpha2;
	
	if (mlow != mup) {
		if (mlow>0.5) {
			mth = (1.0/c2*(pow(mup,c2)-pow(mlow,c2)))/k2;
		} else if (mup<0.5) {
			mth = (2.0/c1*(pow(mup,c1)-pow(mlow,c1)))/k2;
		} else
			mth = (2.0/c1*(pow(0.5,c1)-pow(mlow,c1))+1.0/c2*(pow(mup,c2)-pow(0.5,c2)))/k2;
	} else {
		mth = mlow;
	} 
	
	if (!*N) {
		*N = max(floor((Mcl-MMAX)/mth), 1);
		if (!epoch) printf("Estimated number of necessary stars: %i\n", *N);
		*N = 1;
	}
	
	
	c1 = 1.0-alpha1;  
	c2 = 1.0-alpha2; 
	*mmean = 0.0;				
	*M = 0.0;
	double mostmassive = 0.0;
	
	for (i=0; i<*N; i++) {
		do{
			do {
				xx = drand48();		
				if (xx<k1/k2)   
					star[i][0] = pow(0.5*c1*xx*k2+pow(mlow,c1),1.0/c1);
				else 
					star[i][0] = pow(c2*(xx*k2-k1)+pow(max(0.5,mlow),c2),1.0/c2);
			} while (star[i][0] > MMAX);

			//evolve star for deltat = epoch with SSE (Hurley, Pols & Tout 2002)
			tphys = 0.0;
			kw = 1;
			mass = star[i][0];  //initial mass
			mt = mass;    //actual mass
			ospin = 0.0;
			tphysf = epoch;
			epoch1 = 0.0;
			vkick = 0.0;
			//printf("MASS %.2f", mass);
			star[i][7] = mass;
			evolv1_(&kw, &mass, &mt, &r, &lum, &mc, &rc, &menv, &renv, &ospin, &epoch1, &tms, &tphys, &tphysf, &dtp, &z, zpars, &vkick);
			//printf("-> %.2f\n", mass);
			//if (vkick) printf("KICK: %.5f\n", vkick);
			lostremnants++;
			lostremnantsmass += mt;
		} while (vkick > vesc);
		lostremnants--;
		lostremnantsmass -= mt;
		star[i][0] = mt;
		star[i][8] = kw;
		star[i][9] = epoch1;
		star[i][10] = ospin;
		star[i][11] = r;
		star[i][12] = lum;
		if (star[i][0] > mostmassive) mostmassive = star[i][0];
		*M += star[i][0];
		if ((i==*N-1) && (*M<Mcl)) *N += 1;
	}
	if (lostremnants) printf("Number of ejected compact remnants: %i (%.1f Msun)\n", lostremnants, lostremnantsmass);

	printf("Total mass: %g\t(%i stars)\n",*M,*N);
	*mmean = *M/ *N;
	printf("Most massive star: %g\n",mostmassive);
	if (!epoch) printf("Mean masses theoretical/data: %f %f\n",mth,*mmean);
	
	return 0;
}

int generate_m2(int an, double *mlim, double *alpha, double Mcl, double M_tmp, double *subcount, int *N, double *mmean, double *M, double **star, double MMAX, double epoch, double Z, double Rh, int remnant) {
	int i, j;

	//set up parameters and variables for SSE (Hurley, Pols & Tout 2002)
	int kw;          //stellar type
	double mass;  //initial mass
	double mt;    //actual mass
	double r = 0.0;   //radius
	double lum = 0.0; //luminosity
	double mc = 0.0;  //core mass
	double rc = 0.0;  //core radius
	double menv = 0.0;//envelope mass
	double renv = 0.0;//envelope radius
	double ospin;     //spin
	double epoch1;  //time of birth?
	double tms = 0.0; //main-sequence lifetime
	double tphys;     //initial age
	double tphysf = epoch;//final age
	double dtp = epoch+1; //data store value, if dtp>tphys no data will be stored
	double z = Z;      //metallicity
	double zpars[20];     //metallicity parameters
	double vkick;	 //kick velocity for compact remnants
	double vesc;	//escape velocity of cluster
	int lostremnants = 0; //number of ejected compact remnants
	double lostremnantsmass = 0.0; //mass of ejected compact remnants
	if (Mcl && remnant && (Rh != -1)) {
		vesc = sqrt(2.0*G*Mcl/Rh);
		printf("Escape velocity of cluster = %.4f km/s\n", vesc);
	} else if ((!remnant) || (Rh == -1)) {
		vesc = 1.0E10;
		printf("Keeping all compact remnants\n");
	} else {
		vesc = sqrt(2.0*0.4**N/Rh);
		printf("Estimated escape velocity of cluster assuming mean stellar mass of 0.4 Msun = %.4f km/s\n", vesc);
	}
	
	for (i=0; i<20; i++) zpars[i] = 0;
	zcnsts_(&z,zpars);  //get metallicity parameters
		
	printf("\nSetting up stellar population with Z = %.4f.\n",Z);
	
	if (epoch) printf("\nEvolving stellar population for %.1f Myr.\n",epoch);
	
	double tmp, ml, mup;
	double mostmassive = 0.0;
	*mmean = 0.0;				
	*M = 0.0;
	if (!*N) *N = 1;
		
	for (i = 0; i < an; i++) {
			printf("# <%.2f , %.2f> .. %.2f\n", mlim[i], mlim[i+1], alpha[i]);
	}
	
	for (i = 1; i < an; i++) subcount[i] += subcount[i-1];
	
	for (i = 0; i < *N; i++) {
		do {
			do {
				tmp = drand48() * subcount[an-1];
				for (j = 0; (j < an) && (subcount[j] < tmp); j++);
				if (alpha[j] != -1.) {
					ml = pow(mlim[j], 1. + alpha[j]);
					mup = pow(mlim[j+1], 1. + alpha[j]);
				} else {
					ml = log(mlim[j]);
					mup = log(mlim[j+1]);
				}
				tmp = ml + drand48() * (mup - ml);
				if (alpha[j] != -1.) star[i][0] = pow(tmp, 1. / (1. + alpha[j]));
				else star[i][0] = exp(tmp);
			} while (star[i][0] > MMAX);
			//printf("%8.4f\n", star[i][0]);
		
			//evolve star for deltat = epoch with SSE (Hurley, Pols & Tout 2002)
			tphys = 0.0;
			kw = 1;
			mass = star[i][0];  //initial mass
			mt = mass;    //actual mass
			ospin = 0.0;
			vkick = 0.0;
			tphysf = epoch;
			epoch1 = 0.0;
			//printf("MASS %.2f", mass);
			star[i][7] = mass;
			if (epoch) evolv1_(&kw, &mass, &mt, &r, &lum, &mc, &rc, &menv, &renv, &ospin, &epoch1, &tms, &tphys, &tphysf, &dtp, &z, zpars, &vkick);
			//printf("-> %.2f\n", mass);
			//printf("KICK: %.5f\n", vkick);
			lostremnants++;
			lostremnantsmass += mt;
		} while (vkick > vesc);
		lostremnants--;
		lostremnantsmass -= mt;
		star[i][0] = mt;
		star[i][8] = kw;
		star[i][9] = epoch1;
		star[i][10] = ospin;
		star[i][11] = r;
		star[i][12] = lum;
		
		if (star[i][0] > mostmassive) mostmassive = star[i][0];
		*M += star[i][0];
		if ((i==*N-1) && (*M<Mcl)) *N += 1;
	}

	if (lostremnants) printf("Number of ejected compact remnants: %i (%.1f Msun)\n", lostremnants, lostremnantsmass);
	printf("Total mass: %g\t(%i stars)\n",*M,*N);
	*mmean = *M/ *N;
	printf("Most massive star: %g\n",mostmassive);
	printf("Mean mass: %f\n",*mmean);
	
	return 0;

}

int generate_m3(int *N, double **star, double mlow, double mup, double *M, double *mmean, double MMAX, double Mcl) {
	int i;
	double a1, a2, k1, k2, mb, m;

	//set up mass function parameter
	a1 = 1.3;
	a2 = 2.3;
	mb = 0.5;
	
	k2 = Mcl/(pow(mb, a1-a2)/(2.0-a1)*(pow(mb,2.0-a1) - pow(mlow,2.0-a1)) + 1.0/(2.0-a2)*(pow(MMAX,2.0-a2)-pow(mb,2.0-a2)));
	k1 = k2*pow(mb, a1-a2);

	printf( "Mcl = %f\tMMAX = %f\tk1 = %f\tk2 = %f\n\n",Mcl,MMAX, k1, k2);
	
	*mmean = 0.0;				
	*M = 0.0;
	*N = 0;

	i = 0;
	
	star[i][0] = MMAX; //set first star to MMAX
	star[i][7] = MMAX;
	star[i][8] = 0;
	star[i][9] = 0.0;
	star[i][10] = 0.0;
	star[i][11] = 0.0;
	star[i][12] = 0.0;
	star[i][13] = 0.0;
	star[i][14] = 0.0;

	*M += star[i][0];
	
	do{	
		i++;
		
		if (star[i-1][0] > mb) {
			m = pow(pow(star[i-1][0], 2.0-a2) - star[i-1][0]/k2*(2.0-a2), 1.0/(2.0-a2));
			if (m < mb) {
				m = pow(pow(mb, 2.0-a1) - (2.0-a1)/pow(mb, a1-a2)*(star[i-1][0]/k2 + pow(mb, 2.0-a2)/(2.0-a2) - 1.0/(2.0-a2)*(pow(star[i-1][0], 2.0-a2))), 1.0/(2.0-a1));
			}
		} else {
			m = pow((a1-2.0)/k1*star[i-1][0] + pow(star[i-1][0], 2.0-a1), 1.0/(2.0-a1));
		}		
		if (m<mlow) {
			i--;
			break;
		}
		star[i][0] = m;
		star[i][7] = m;
		star[i][8] = 0;
		star[i][9] = 0.0;
		star[i][10] = 0.0;
		star[i][11] = 0.0;
		star[i][12] = 0.0;
		star[i][13] = 0.0;
		star[i][14] = 0.0;

		*M += star[i][0];
		printf("-------%f\n",*M);
		
	} while ( *M < Mcl);

	*N = i+1;
	
	printf("Total mass: %g\t(%i stars)\n",*M,*N);
	*mmean = *M/ *N;
	printf("Most massive star: %g\n",MMAX);

	return 0;
}

double subint(double min, double max, double alpha) {
	if (alpha == 0.) return log(max / min);
	else return (pow(max, alpha) - pow(min, alpha)) / alpha;
}

double mlow(double mhigh, double alpha, double norma, double delta) {
	if (alpha == 0.) return mhigh / exp(delta / norma);
	else return pow(pow(mhigh, alpha) - delta * alpha / norma, 1. / alpha);
}

int generate_m4(int *N, double **star, double alpha, double beta, double mu,  double mlow, double mup, double *M, double *mmean, double MMAX, double Mcl, double epoch, double Z, double Rh, int remnant) {
	
	// L_3 IMF
	// needs functions
	//          double alogam ( double x, int *ifault );
	//          double betain ( double x, double p, double q, double beta, int *ifault );
	//          double r8_abs ( double x );
	// taken from the Applied Statistics Algorithm 63; ASA063 is a C library which evaluates the incomplete Beta function, by KL Majumder and G Bhattacharjee.

	//set up parameters and variables for SSE (Hurley, Pols & Tout 2002)
	int kw;          //stellar type
	double mass;  //initial mass
	double mt;    //actual mass
	double r = 0.0;     //radius
	double lum = 0.0;   //luminosity
	double mc = 0.0;    //core mass
	double rc = 0.0;    //core radius
	double menv = 0.0;  //envelope mass
	double renv = 0.0;  //envelope radius
	double ospin;		//spin
	double epoch1;    //time spent in current evolutionary state
	double tms = 0.0;   //main-sequence lifetime
	double tphys;       //initial age
	double tphysf = epoch;//final age
	double dtp = epoch+1; //data store value, if dtp>tphys no data will be stored
	double z = Z;      //metallicity
	double zpars[20];     //metallicity parameters
	double vkick;	 //kick velocity for compact remnants
	double vesc;	//escape velocity of cluster
	int lostremnants = 0; //number of ejected compact remnants
	double lostremnantsmass = 0.0; //mass of ejected compact remnants
	double mostmassive = 0.0;	
	

	double G_l;
	double G_u;
	double u;
	double x;
	int i;
	double mth;
	
	double gamma;
	double mpeak;
	// variables for the mean
	double B_l;
	double B_u;
	double a;
	double b;
	double beta_log;
	int ifault;
	
	printf("\nL3 Parameters:\n");
	printf("alpha = %1.2f beta = %1.2f mu = %1.2f m_low = %3.3f m_up = %3.3f\n",alpha,beta,mu,mlow,mup);

	mpeak = mu*pow((beta-1.0),1.0/(alpha-1.0));
	gamma = alpha + beta*(1.0-alpha);

	printf("'Peak' mass = %3.3f   Effective low mass exponent gamma = %1.2f\n\n",mpeak,gamma);
	
	// normalisation
	G_l = pow(1.0 + pow(mlow/mu, 1.0-alpha) , 1.0-beta);
	G_u = pow(1.0 + pow(mup/mu, 1.0-alpha) , 1.0-beta);
	
	
	// Calculate the mean mass
	a = (2.0-alpha)/(1.0-alpha);
	b = beta - a;
	
	// logarithm of the beta function necessary to calculate the incomplete beta fn
	beta_log = alogam( a, &ifault ) + alogam( b, &ifault ) - alogam( a + b, &ifault );
	
	x = pow(mlow/mu,1.0-alpha);
	x = x/(1.0+x);
	B_l = betain(x, a, b, beta_log, &ifault)*exp(beta_log);
	
	x = pow(mup/mu,1.0-alpha);
	x = x/(1.0+x);
	B_u = betain(x, a, b, beta_log, &ifault)*exp(beta_log);
	mth = mu*(1.0-beta)*(B_u - B_l)/(G_u - G_l) ;
	
	if (Mcl && remnant && (Rh != -1)) {
		vesc = sqrt(2.0*G*Mcl/Rh);
		printf("Escape velocity of cluster = %.4f km/s\n", vesc);
	} else if ((!remnant) || (Rh == -1)) {
		vesc = 1.0E10;
		printf("Keeping all compact remnants\n");
	} else {
		vesc = sqrt(2.0*0.4**N/Rh);
		printf("Estimated escape velocity of cluster assuming mean stellar mass of 0.4 Msun = %.4f km/s\n", vesc);
	}
	for (i=0; i<20; i++) zpars[i] = 0;
	zcnsts_(&z,zpars);  //get metallicity parameters
	
	printf("\nSetting up stellar population with Z = %.4f.\n",Z);
	
	if (epoch) printf("\nEvolving stellar population for %.1f Myr.\n",epoch);
	
	if (!*N) {
		*N = max(floor((Mcl-mup)/mth), 1);
		if (!epoch) printf("Estimated number of necessary stars: %i\n", *N);
		*N = 1;
	}	

	// generate random masses
	for (i=0; i<*N; i++) {
		do {
			u = drand48();
			x = u*(G_u-G_l)+G_l;
			x = pow(x,1/(1-beta)) - 1.0;
			star[i][0] = mu*pow(x,1.0/(1.0-alpha));

			tphys = 0.0;
			kw = 1;
			mass = star[i][0];  //initial mass
			mt = mass;    //actual mass
			ospin = 0.0;
			tphysf = epoch;
			epoch1 = 0.0;
			vkick = 0.0;
			//printf("MASS %.2f", mass);
			star[i][7] = mass;
			evolv1_(&kw, &mass, &mt, &r, &lum, &mc, &rc, &menv, &renv, &ospin, &epoch1, &tms, &tphys, &tphysf, &dtp, &z, zpars, &vkick);
			//printf("-> %.2f\n", mass);
			//if (vkick) printf("KICK: %.5f\n", vkick);
			lostremnants++;
			lostremnantsmass += mt;
		} while (vkick > vesc);
		lostremnants--;
		lostremnantsmass -= mt;
		star[i][0] = mt;
		star[i][8] = kw;
		star[i][9] = epoch1;
		star[i][10] = ospin;
		star[i][11] = r;
		star[i][12] = lum;
		if (star[i][0] > mostmassive) mostmassive = star[i][0];
		
		*M += star[i][0];

		if ((i==*N-1) && (*M<Mcl)) *N += 1;
	}
	if (lostremnants) printf("Number of ejected compact remnants: %i (%.1f Msun)\n", lostremnants, lostremnantsmass);

	*mmean = *M/(1.0**N);

	printf("Total mass: %g\t(%i stars)\n",*M,*N);
	printf("Most massive star: %g\n",mostmassive);
	if (!epoch) printf("Mean masses theoretical/data: %f %f\n",mth,*mmean);
	
	return 0;
	
}

double alogam(double x, int *ifault) {
	/*
	 Purpose:
	 
	 ALOGAM computes the logarithm of the Gamma function.
	 
	 Licensing:
	 
	 This code is distributed under the GNU LGPL license. 
	 
	 Modified:
	 
	 20 October 2010
	 
	 Author:
	 
	 Original FORTRAN77 version by Malcolm Pike, David Hill.
	 C version by John Burkardt.
	 
	 Reference:
	 
	 Malcolm Pike, David Hill,
	 Algorithm 291:
	 Logarithm of Gamma Function,
	 Communications of the ACM,
	 Volume 9, Number 9, September 1966, page 684.
	 
	 Parameters:
	 
	 Input, double X, the argument of the Gamma function.
	 X should be greater than 0.
	 
	 Output, int *IFAULT, error flag.
	 0, no error.
	 1, X <= 0.
	 
	 Output, double ALOGAM, the logarithm of the Gamma
	 function of X.
	 */
	
	double f;
	double value;
	double y;
	double z;
	
	if ( x <= 0.0 )
	{
		*ifault = 1;
		value = 0.0;
		return value;
	}
	
	*ifault = 0;
	y = x;
	
	if ( x < 7.0 )
	{
		f = 1.0;
		z = y;
		
		while ( z < 7.0 )
		{
			f = f * z;
			z = z + 1.0;
		}
		y = z;
		f = - log ( f );
	}
	else
	{
		f = 0.0;
	}
	
	z = 1.0 / y / y;
	
	value = f + ( y - 0.5 ) * log ( y ) - y 
    + 0.918938533204673 + 
    ((( 
	   - 0.000595238095238   * z 
	   + 0.000793650793651 ) * z 
	  - 0.002777777777778 ) * z 
	 + 0.083333333333333 ) / y;
	
	return value;
}

double betain(double x, double p, double q, double beta, int *ifault) {
	/*
	 Purpose:
	 
	 BETAIN computes the incomplete Beta function ratio.
	 
	 Licensing:
	 
	 This code is distributed under the GNU LGPL license. 
	 
	 Modified:
	 
	 31 October 2010
	 
	 Author:
	 
	 Original FORTRAN77 version by KL Majumder, GP Bhattacharjee.
	 C version by John Burkardt.
	 
	 Reference:
	 
	 KL Majumder, GP Bhattacharjee,
	 Algorithm AS 63:
	 The incomplete Beta Integral,
	 Applied Statistics,
	 Volume 22, Number 3, 1973, pages 409-411.
	 
	 Parameters:
	 
	 Input, double X, the argument, between 0 and 1.
	 
	 Input, double P, Q, the parameters, which
	 must be positive.
	 
	 Input, double BETA, the logarithm of the complete
	 beta function.
	 
	 Output, int *IFAULT, error flag.
	 0, no error.
	 nonzero, an error occurred.
	 
	 Output, double BETAIN, the value of the incomplete
	 Beta function ratio.
	 */
	
	
	double acu = 0.1E-14;
	double ai;
	double cx;
	int indx;
	int ns;
	double pp;
	double psq;
	double qq;
	double rx;
	double temp;
	double term;
	double value;
	double xx;
	
	value = x;
	*ifault = 0;
	/*
	 Check the input arguments.
	 */
	if ( p <= 0.0 || q <= 0.0 )
	{
		*ifault = 1;
		return value;
	}
	
	if ( x < 0.0 || 1.0 < x )
	{
		*ifault = 2;
		return value;
	}
	/*
	 Special cases.
	 */
	if ( x == 0.0 || x == 1.0 )
	{
		return value;
	}
	/*
	 Change tail if necessary and determine S.
	 */
	psq = p + q;
	cx = 1.0 - x;
	
	if ( p < psq * x )
	{
		xx = cx;
		cx = x;
		pp = q;
		qq = p;
		indx = 1;
	}
	else
	{
		xx = x;
		pp = p;
		qq = q;
		indx = 0;
	}
	
	term = 1.0;
	ai = 1.0;
	value = 1.0;
	ns = ( int ) ( qq + cx * psq );
	/*
	 Use the Soper reduction formula.
	 */
	rx = xx / cx;
	temp = qq - ai;
	if ( ns == 0 )
	{
		rx = xx;
	}
	
	for ( ; ; )
	{
		term = term * temp * rx / ( pp + ai );
		value = value + term;;
		temp = r8_abs ( term );
		
		if ( temp <= acu && temp <= acu * value )
		{
			value = value * exp ( pp * log ( xx ) 
								 + ( qq - 1.0 ) * log ( cx ) - beta ) / pp;
			
			if ( indx )
			{
				value = 1.0 - value;
			}
			break;
		}
		
		ai = ai + 1.0;
		ns = ns - 1;
		
		if ( 0 <= ns )
		{
			temp = qq - ai;
			if ( ns == 0 )
			{
				rx = xx;
			}
		}
		else
		{
			temp = psq;
			psq = psq + 1.0;
		}
	}
	
	return value;
}

double r8_abs(double x) {
	/*
	 Purpose:
	 
	 R8_ABS returns the absolute value of an R8.
	 
	 Licensing:
	 
	 This code is distributed under the GNU LGPL license. 
	 
	 Modified:
	 
	 07 May 2006
	 
	 Author:
	 
	 John Burkardt
	 
	 Parameters:
	 
	 Input, double X, the quantity whose absolute value is desired.
	 
	 Output, double R8_ABS, the absolute value of X.
	 */
	
	double value;
	
	if ( 0.0 <= x )
	{
		value = + x;
	} 
	else
	{
		value = - x;
	}
	return value;
}

int generate_plummer(int N, double **star, double rtide, double rvir, double D, int symmetry){
	int i, h;
	double a[9], ri, sx, sv, rcut;
	double r_norm, v_norm;

	//Scale length
	sx = 1.5*TWOPI/16.0;
	sv = sqrt(1.0/sx);
	
	printf("Setting cut-off radius of Plummer sphere to approximate tidal radius\n");
	rcut = rtide/(sx*rvir);		//cut-off radius for Plummer sphere = tidal radius in scaled length

	printf("\nGenerating Orbits:\n");	

	if (D>=3.0) {
		
		for (i=0;i<N;i++) {
		
		if ((i/1000)*1000 == i) printf("Generating orbit #%i\n", i);
			
		//Positions
		do {
			do { 
				a[1] = drand48();
			} while (a[1]<1.0E-10);
			ri = 1.0/sqrt(pow(a[1],-2.0/3.0) - 1.0);
			
			a[2] = drand48(); 
			a[3] = drand48();
		
			star[i][3] = (1.0 - 2.0*a[2])*ri;
			star[i][1] = sqrt(ri*ri-pow(star[i][3],2))*cos(TWOPI*a[3]);
			star[i][2] = sqrt(ri*ri-pow(star[i][3],2))*sin(TWOPI*a[3]); 
		} while (sqrt(pow(star[i][1],2)+pow(star[i][2],2)+pow(star[i][3],2))>rcut); //reject particles beyond tidal radius
		
		//velocities
		do {
			a[4] = drand48(); 
			a[5] = drand48(); 
			a[6] = pow(a[4],2)*pow(1.0 - pow(a[4],2),3.5);
		} while (0.1*a[5]>a[6]);
		
		a[8] = a[4]*sqrt(2.0)/pow(1.0 + ri*ri,0.25);
		a[6] = drand48(); 
		a[7] = drand48(); 
			
		star[i][6] = (1.0 - 2.0*a[6])*a[8];
		star[i][4] = sqrt(a[8]*a[8] - pow(star[i][6],2))*cos(TWOPI*a[7]);
		star[i][5] = sqrt(a[8]*a[8] - pow(star[i][6],2))*sin(TWOPI*a[7]);
		}
	} else {
		double xcut;
		xcut = 0.000;
		while (1.0/sqrt(pow(xcut,-2.0/3.0) - 1.0)<=rcut) xcut+=0.00001;

		fractalize(D, N, star, 1, symmetry);
		for (i=0;i<N;i++) {
			if ((i/1000)*1000 == i) printf("Generating orbit #%i\n", i);
			ri = sqrt(pow(star[i][1],2)+pow(star[i][2],2)+pow(star[i][3],2))*xcut;
			r_norm = 1.0/sqrt(pow(ri,-2.0/3.0) - 1.0);
			v_norm = sqrt(2.0)/pow(1.0 + ri*ri,0.25);
			for (h=1;h<4;h++) star[i][h] *= r_norm/ri;
			for (h=4;h<7;h++) star[i][h] *= v_norm;
		}
		
	}

	return 0;
}

int generate_king(int N, double W0, double **star, double *rvir, double *rh, double *rking, double D, int symmetry){
	
	//ODE variables
	int M = 10001;				//Number of interpolation points
	int KMAX = 10000;			//Maximum number of output steps of the integrator
	
	int i,j,k;
	double h;
	double den;
	double xstart, ystart0, ystart1;
	double x1, x2;
	double xp[KMAX], x[M];	
	int kount = 0;
	double yking[M][2], mass[M];
	double rmin = 0.2;
	
	double **yp;
	yp = (double **)calloc(KMAX,sizeof(double *));
	for (i=0;i<KMAX;i++){
		yp[i] = (double *)calloc(2,sizeof(double));
		if (yp[i] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}
	
	double pot = 0.0;
	double totmas;
	double zh;
	
	double ve, cg1;
	double fmass, r2;
	double costh, sinth, phi, r, xstar, ystar, zstar;
	double w1, w, wj1, wj;
	double vmax, speed, fstar, ustar, vstar, wstar, mstar;
	double ri;
	
	double **coord;
	coord = (double **)calloc(N,sizeof(double *));
	for (j=0;j<N;j++){
		coord[j] = (double *)calloc(6,sizeof(double));
		if (coord[j] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}

	
	
	if (W0>12.0) {
		printf("W0 too large\n");
		return 0;
	} else if (W0 < 0.2) {
		printf("W0 too small\n");
		return 0;
	}
	
	
	
	
	/***************************
	 * INTERPOLATE KING VALUES *
	 ***************************/
	
	h = pow(W0-rmin,2)/(1.0*M-1.0);	//step size for interpolation
	den = densty(W0);			//central density
	
	//x is King's W, y[0] is z**2, y[1] is 2*z*dz/dW, where z is scaled radius
	//so x(y[0]) is energy as function of radius^2 and x(y[1]) is the derivative
	
	xstart = 0.000001;
	ystart0 = 2.0*xstart/3.0;
	ystart1 = -2.0/3.0;
	
	x1 = W0 - xstart;
	x2 = 0.0;
	
	//integrate Poisson's eqn
	printf("Integrating Poisson's equation\n");

	odeint(ystart0, ystart1, x1, x2, den, &kount, xp, yp, M, KMAX);
	
	printf("No of integration steps = %i\n",kount);
	
	
	
	//interpolate yking 
	for (k=0;k<M;k++) {
		x[k] = W0-rmin-sqrt(h*k); 
		
		
		if (x[k] > xp[0]) {
			yking[k][0] = (W0-x[k])*yp[0][0]/(W0 - xp[0]);
			yking[k][1] = yp[0][1];
			printf("+");
		} else {
			
			i = 0;
			
			do {
				if ((x[k]-xp[i])*(x[k]-xp[i+1]) <= 0.0) {
					yking[k][0] = yp[i][0] + (yp[i+1][0]-yp[i][0])*(x[k]-xp[i])/(xp[i+1]-xp[i]);
					yking[k][1] = yp[i][1] + (yp[i+1][1]-yp[i][1])*(x[k]-xp[i])/(xp[i+1]-xp[i]);
					goto jump;
				} else {
					i++;
				}
			} while (i<kount);
		}
		
	jump:
		
		
		if (i >= kount) {
			yking[k][0] = yp[kount][0];
			yking[k][1] = yp[kount][1];
		}
		
		
		//get mass as a function of radius 
		mass[k] = -2.0*pow(yking[k][0],1.5)/yking[k][1];
		
		
		//calculate total potential energy
		if (k == 0) {
			pot = -0.5*pow(mass[k],2)/sqrt(yking[k][0]); 
		} else {
			pot -=  0.5*(pow(mass[k],2) - pow(mass[k-1],2))/(0.5*(sqrt(yking[k][0]) + sqrt(yking[k-1][0])));
		}
	}	
	
	
	
	
	/******************************
	 * DETERMINE BASIC QUANTITIES *
	 ******************************/
	
	//Total mass
	totmas = -2.0*pow(yking[M-2][0],1.5)/yking[M-2][1];
	
	//Half-mass radius
	k=0;
	
	do {
		k++;
	} while (mass[k] < 0.5*totmas);
	
	zh = yking[k][0] - (yking[k][0]-yking[k-1][0])*(mass[k]-0.5*totmas)/(mass[k]-mass[k-1]);
	*rh = sqrt(zh);
	
	//Virial radius and King radius
	*rvir = -pow(totmas,2)/(2.0*pot);
	*rking = sqrt(yking[M-2][0]);
	
	//Central velocity dispersion**2 (3-dimensional) 
	ve = sqrt(2.0*W0);
	cg1 = (-pow(ve,3)*exp(-pow(ve,2)/2.0) - 3.0*ve*exp(-pow(ve,2)/2.0) + 3.0/2.0*sqrt(2.0*PI)*erf(ve*sqrt(2.0)/2.0) - pow(ve,5)*exp(-pow(ve,2)/2.0)/5.0)/(-ve*exp(-pow(ve,2)/2.0) + sqrt(2.0*PI)*erf(ve*sqrt(2.0)/2.0)/2.0 - pow(ve,3)*exp(-pow(ve,2)/2.0)/3.0);
	
	printf("\nTheoretical values:\n");
	printf("Total mass (King units) = %g\n", totmas);
	printf("Viriral radius (King units) = %g\n", *rvir);
	printf("Half-mass radius (King units) = %g\t(Nbody units) = %g\n", *rh, *rh/ *rvir);
	printf("Edge radius (King units) = %g\t(Nbody units) = %g\n", *rking, *rking/ *rvir);
	printf("Concentration = %g\n", log10(*rking));
	printf("Core radius (King units) = %g\t(Nbody units) = %g\n", 1.0, 1.0/ *rvir);
	printf("3d velocity dispersion**2: %g (central)\t %g (global)\n", cg1, -pot/totmas);
	
	
	
	/***************************
	 * GENERATE STAR POSITIONS *
	 ***************************/
	
	printf("\nGenerating Stars:\n");	
	
	if (D>=3.0) {
		for (i=0;i<N;i++) {
		
			if ((i/1000)*1000 == i) printf("Generating orbits #%i\n", i);
		
			fmass = drand48()*mass[M-1];
		
			if (fmass < mass[0]) {
				r2 = pow(fmass/mass[0],2.0/3.0)*yking[0][0];
			} else {
				j = 0;
				do {
					j++;
					if (j>M) printf("WARNING: failing iteration\n");
				} while (mass[j] <= fmass);
			
				r2 = yking[j-1][0] + (fmass-mass[j-1])*(yking[j][0] - yking[j-1][0])/(mass[j]-mass[j-1]);
			}
		
		
			r = sqrt(r2);
			costh = 2.0*drand48()-1.0;
			sinth = sqrt(1.0-pow(costh,2));
			phi = 2.0*PI*drand48();
			xstar = r*sinth*cos(phi);
			ystar = r*sinth*sin(phi);
			zstar = r*costh;

			
			/****************
			 * CHOOSE SPEED *
			 ****************/
						
			if (r < *rking) {
				if (r < sqrt(yking[0][0])) {
					w1 = x[0];
					w = W0 - (r2/yking[0][0])*(W0 - w1);
				} else {
					j = 0;
					do {
						j++;
					} while (r > sqrt(yking[j][0]));
					wj1 = x[j-1];
					wj = x[j];
					w = wj1 + (r2-yking[j-1][0])*(wj-wj1)/(yking[j][0]-yking[j-1][0]);
				}
			} else {
				printf("radius too big\n");
			}
		
			vmax = sqrt(2.0*w);
			do {
				speed = vmax*drand48();
				fstar = pow(speed,2)*(exp(-0.5*pow(speed,2))-exp(-1.0*w));
			} while (fstar < 2.0*drand48()/exp(1.0));
		
			costh = 2.0*drand48()-1.0;
			phi = 2.0*PI*drand48();
			sinth = sqrt(1.0-pow(costh,2));
			ustar = speed*sinth*cos(phi);
			vstar = speed*sinth*sin(phi);
			wstar = speed*costh;
			mstar = star[i][0];
		
			//printf("i: %i\tr=%g\tm=%.5f\tx=%.5f y=%.5f z=%.5f\tvx=%.5f vy=%.5f vz=%.5f\n",i,sqrt(r2),mstar,xstar,ystar,zstar,ustar,vstar,wstar);
		
		
			coord[i][0] = xstar;
			coord[i][1] = ystar;
			coord[i][2] = zstar;
			coord[i][3] = ustar;
			coord[i][4] = vstar;
			coord[i][5] = wstar;
		
		}

		for (i=0;i<N;i++) {
			for (j=0;j<3;j++) {
				star[i][j+1] = 1.0*coord[i][j];
				star[i][j+4] = 1.0*coord[i][j+3];
			}
		}

	} else {
		
		fractalize(D, N, star, 1, symmetry);
		
		for (i=0;i<N;i++) {

			if ((i/1000)*1000 == i) printf("Generating orbits #%i\n", i);
			
			ri =  sqrt(pow(star[i][1],2)+pow(star[i][2],2)+pow(star[i][3],2));
			fmass = ri*mass[M-1];
			
			if (fmass < mass[0]) {
				r2 = pow(fmass/mass[0],2.0/3.0)*yking[0][0];
			} else {
				j = 0;
				do {
					j++;
					if (j>M) printf("WARNING: failing iteration\n");
				} while (mass[j] <= fmass);
				
				r2 = yking[j-1][0] + (fmass-mass[j-1])*(yking[j][0] - yking[j-1][0])/(mass[j]-mass[j-1]);
			}
			
			r = sqrt(r2);

			
			/****************
			 * CHOOSE SPEED *
			 ****************/
			
			if (r < *rking) {
				if (r < sqrt(yking[0][0])) {
					w1 = x[0];
					w = W0 - (r2/yking[0][0])*(W0 - w1);
				} else {
					j = 0;
					do {
						j++;
					} while (r > sqrt(yking[j][0]));
					wj1 = x[j-1];
					wj = x[j];
					w = wj1 + (r2-yking[j-1][0])*(wj-wj1)/(yking[j][0]-yking[j-1][0]);
				}
			} else {
				printf("radius too big\n");
			}
			
			
			vmax = sqrt(2.0*w);
			do {
				speed = vmax*drand48();
				fstar = pow(speed,2)*(exp(-0.5*pow(speed,2))-exp(-1.0*w));
			} while (fstar < 2.0*drand48()/exp(1.0));
						
			for (k=1;k<4;k++) star[i][k] *= r/ri;
			for (k=4;k<7;k++) star[i][k] *= speed;
		}
		
	}
	
		
	for (j=0;j<KMAX;j++) free (yp[j]);
	free(yp);

	for (j=0;j<N;j++) free (coord[j]);
	free(coord);
	
	return 0;
	

}

double densty(double z){
	double den = -sqrt(z)*(z+1.5)+0.75*sqrt(PI)*exp(z)*erf(sqrt(z));
	return den;
}

int odeint(double ystart0, double ystart1, double x1, double x2, double den, int *kount, double *xp, double **yp, int M, int KMAX) {
	
	double HMIN = 0.0;			//Minimum step size
	double H1 = 0.0001;			//Size of first step
	int MAXSTP = 100000;		//Maximum number of steps for integration
	double TINY = 1.0E-30;		//To avoid certain numbers get zero
	double DXSAV = 0.0001;		//Output step size for integration
	double TOL = 1.0E-12;		//Tolerance of integration
	
	
	double x;
	double h;
	int i,j;
	double y[2];
	double hdid, hnext;
	double xsav;
	double dydx[2], yscal[2];
	
	x = x1;   //King's W parameter
	if (x2-x1 >= 0.0) {
		h = sqrt(pow(H1,2));
	} else {
		h = - sqrt(pow(H1,2));
	}  //step size
	
	y[0] = ystart0;  //z^2
	y[1] = ystart1;  //2*z*dz/dW  where z is scaled radius	
	
	xsav = x-DXSAV*2.0;
	
	for (i=0;i<MAXSTP;i++) {        //limit integration to MAXSTP steps
		
		derivs(x,y,dydx,den); //find derivative
		
		for (j=0;j<2;j++) {
			yscal[j] = sqrt(pow(y[j],2))+sqrt(h*dydx[j])+TINY;  //advance y1 and y2
		}
		
		if (sqrt(pow(x-xsav,2)) > sqrt(pow(DXSAV,2))) {
			if (*kount < KMAX-1) {
				xp[*kount] = x;
				for (j=0;j<2;j++) {
					yp[*kount][j] = y[j];
				}
				*kount = *kount + 1;
				xsav = x;
			}
		}  //store x, y1 and y2 if the difference in x is smaller as the desired output step size DXSAV
		
		if (((x+h-x2)*(x+h-x1)) > 0.0) h = x2-x;
		
		rkqc(y,dydx,&x,&h,den,yscal, &hdid, &hnext, TOL);	//do a Runge-Kutta step
		
		if ((x-x2)*(x2-x1) >= 0.0) {
			ystart0 = y[0];
			ystart1 = y[1];
			
			xp[*kount] = x;
			for (j=0;j<2;j++) {
				yp[*kount][j] = y[j];
			}
			return 0;	
			*kount = *kount +1;
		}
		
		if (sqrt(pow(hnext,2)) < HMIN) {
			printf("Stepsize smaller than minimum.\n");
			return 0;
		}
		
		h = hnext;
	} 
	
	return 0;
}

int derivs(double x, double *y, double *dydx, double den){
	
	double rhox;
	
	if (x >= 0.0) {
		rhox =-sqrt(x)*(x+1.5)+0.75*sqrt(PI)*exp(x)*erf(sqrt(x));
	} else {
		rhox = 0.0;
	}
	
	dydx[0]= y[1];
	dydx[1] = 0.25*pow(y[1],2)*(6.0+9.0*y[1]*rhox/den)/y[0];	
	
	return 0;
}

int rkqc(double *y,double *dydx, double *x, double *h, double den, double *yscal, double *hdid, double *hnext, double TOL){
	
	double safety = 0.9;
	double fcor = 0.0666666667;
	double errcon = 6.0E-4;
    double pgrow = -0.20;
	double pshrnk = -0.25;
	double xsav;
	int i;
	double ysav[2],  dysav[2], ytemp[2];
	double errmax;
	double hh;
	
	xsav = *x;
	
	for (i=0;i<2;i++) {
		ysav[i] = y[i];
		dysav[i] = dydx[i];
	}
	
	do {
		hh = 0.5**h;
		rk4(xsav, ysav, dysav, hh, ytemp, den);
		
		*x = xsav + hh;
		derivs(*x,ytemp,dydx,den); //find derivative
		rk4(*x, ytemp, dydx, hh, y, den);
		
		*x = xsav + *h;
		if (*x  == xsav) {
			printf("ERROR: Stepsize not significant in RKQC.\n");
			return 0;
		}
		rk4(xsav, ysav, dysav, *h, ytemp, den);
		
		errmax = 0.0;
		for (i=0;i<2;i++) {
			ytemp[i] = y[i] - ytemp[i];
			errmax = max(errmax, sqrt(pow(ytemp[i]/yscal[i],2)));
		}
		errmax /= TOL;
		if (errmax > 1.0) *h = safety**h*(pow(errmax,pshrnk)); //if integration error is too large, decrease h
		
	} while (errmax > 1.0);
	
	
	*hdid = *h;
	if (errmax > errcon) {
		*hnext = safety**h*(pow(errmax,pgrow));//increase step size for next step
	} else {
		*hnext = 4.0**h;//if integration error is very small increase step size significantly for next step
	}
	
	for (i=0;i<2;i++) {
		y[i] += ytemp[i]*fcor;
	}
	
	return 0;
	
}

int rk4(double x, double *y, double *dydx, double h, double *yout, double den){
	double hh, h6, xh;
	double yt[2], dyt[2], dym[2];
	int i;
	
	hh=h*0.5;
	h6=h/6.0;
	xh=x+hh;
	for (i=0;i<2;i++) {
		yt[i] = y[i] + hh*dydx[i];
	}
	
	derivs(xh,yt,dyt,den); //find derivative
	for (i=0;i<2;i++) {
		yt[i] = y[i] + hh*dyt[i];
	}
	
	derivs(xh,yt,dym,den); //find derivative
	for (i=0;i<2;i++) {
		yt[i] = y[i] + h*dym[i];
		dym[i] += dyt[i];
	}
	
	derivs(x+h,yt,dyt,den); //find derivative
	for (i=0;i<2;i++) {
		yout[i] = y[i] + h6*(dydx[i]+dyt[i]+2.0*dym[i]);
	}
	
	return 0;
}

int generate_subr(int N, double S, double **star, double rtide, double rvir){
	long   i, j;
	double r1, r2, tmp, rcut, rtemp;
	double U_tot, U_tot1, K_tot;
	double UT_sub, M_sub;
	double maxim, beta;
	int    N_star;
	double M_tot;
	star1 = NULL;
	star2 = NULL;
	
	if (S<0.5) {
		printf("Setting cut-off radius of Subr model to the approximate tidal radius\n");
		rcut = rtide/rvir;		//cut-off radius for Plummer sphere = tidal radius
	} else {
		printf("Attention! Be aware that for about S>0.5 the Subr mass segregation may generate a virially hot system\nSetting cut-off radius of Subr model to twice the approximate tidal radius\n");
		rcut = 2.0*rtide/rvir;		//cut-off radius for Plummer sphere = tidal radius
	}
	M_tot = 0.;
	N_star = N;
	star1 = (struct t_star1 *)realloc(star1, N_star * sizeof(struct t_star1));
	star2 = (struct t_star2 *)realloc(star2, N_star * sizeof(struct t_star2));
	
	for (i=0; i<N_star;i++) {
		star1[i].mass = star[i][0];
		star1[i].m0 = star[i][7];
		star1[i].kw = star[i][8];
		star1[i].epoch = star[i][9];
		star1[i].spin = star[i][10];
		star1[i].rstar = star[i][11];
		star1[i].lum = star[i][12];
		star1[i].epochstar = star[i][13];
		star1[i].zstar = star[i][14];
		M_tot += star1[i].mass;
	};
		
	quick(0, N_star-1);
	
	M_sub = U_tot1 = 0.;
	for (i = 0; i < N_star; i++) {
		star2[i].U = star1[i].U_tmp = star2[i].UT_add = 0.;
		star1[i].mass /= M_tot;
		M_sub += star1[i].mass;
		star2[i].M_sub = M_sub;
		if (S > 0.) star1[i].Ui = star1[i].mass * pow(M_sub, -S);
		else star1[i].Ui = star1[i].mass;
		for (j = 0; j < i; j++) {
			tmp = star1[i].Ui * star1[j].Ui;
			star2[i].UT_add += tmp;
			star1[i].U_tmp += tmp;
			star1[j].U_tmp += tmp;
		};
		U_tot1 += star2[i].UT_add;
	};
	
	UT_sub = U_tot = K_tot = 0.;
	
	for (i = 0; i < N_star; i++) {
		if ((i/1000)*1000 == i) printf("Generating orbit #%li\n", i);
		star2[i].UT_add *= 0.5 / U_tot1;
		star2[i].UT = star1[i].U_tmp * 0.5 / U_tot1;
		UT_sub += star2[i].UT_add;
		do {
			position(i, U_tot, UT_sub, S);
			rtemp = sqrt(pow(star1[i].r[0],2)+pow(star1[i].r[1],2)+pow(star1[i].r[2],2));
		} while (rtemp>rcut);
		U_tot += star1[i].mass * star2[i].U_sub;
	};
	
	for (i = 0; i < N_star; i++) {
		maxim = 0.25 * star1[i].mass / star2[i].UT;
		find_beta(&maxim, &beta);
		if (beta > 0.)
		{
			do
			{
				r1 = maxim * drand48();
				r2 = drand48();
				r2 *= r2;
				tmp = 1. - r2;
			}
			while (r1 > r2 * exp(beta * log(tmp)));
		}
		else
		{
			do
			{
				r1 = maxim * drand48();
				r2 = sin(M_PI_2 * drand48());
				r2 *= r2;
				tmp = sqrt(1. - r2);
			}
			while (r1 > r2 * exp((2.* beta + 1.) * log(tmp)));
		};
		
		star2[i].E = r2 * (star2[i].U + star2[i].U_sub);
		
		tmp = sqrt(2. * star2[i].E);
		isorand(star2[i].v);
		star2[i].v[0] *= tmp;
		star2[i].v[1] *= tmp;
		star2[i].v[2] *= tmp;
		K_tot += star1[i].mass * star2[i].E;
	};
	
	for (i = 0; i < N_star; i++){
		//printf("%.12f  %18.14f  %18.14f  %18.14f  %16.12f  %16.12f  %16.12f\n",star1[i].mass, star1[i].r[0], star1[i].r[1], star1[i].r[2], star2[i].v[0], star2[i].v[1], star2[i].v[2]);
		star[i][0] = star1[i].mass;
		star[i][1] = star1[i].r[0];
		star[i][2] = star1[i].r[1];
		star[i][3] = star1[i].r[2];
		star[i][4] = star2[i].v[0];
		star[i][5] = star2[i].v[1];
		star[i][6] = star2[i].v[2];
		star[i][7] = star1[i].m0;
		star[i][8] = star1[i].kw;
		star[i][9] = star1[i].epoch;
		star[i][10] = star1[i].spin;
		star[i][11] = star1[i].rstar;
		star[i][12] = star1[i].lum;
		star[i][13] = star1[i].epochstar;
		star[i][14] = star1[i].zstar;
		
	};
	
	return 0;
}

void quick(int start, int stop) {
	int i, j;
	double temp, median;
	
	median = 0.5 * (star1[start].mass + star1[stop].mass);
	
	i = start;
	j = stop;
	while (i <= j)
	{
		while ((i <= stop) && (star1[i].mass > median)) i++;
		while ((j >= start) && (star1[j].mass < median)) j--;
		if (j > i)
		{
			SWAP(i,j);
			i++;
			j--;
		}
		else if (i == j)
		{
			i++;
			j--;
		};
	};
	
	if (start + 1 < j) quick(start, j);
	else if ((start < j) && (star1[start].mass < star1[j].mass)) SWAP(start, j);
	
	if (i + 1 < stop) quick(i, stop);
	else if ((i < stop) && (star1[i].mass < star1[stop].mass)) SWAP(i, stop)
}

void position(int id, double U_sub, double UT_sub, double S){
	long   i;
	double rx, ry, rz, r1;
	double rfac;
	
	rfac = pow(star2[id].M_sub, 2. * S) / (1. - S);
	
	star2[id].ntry = 0;
	
	do
	{
		star2[id].U_sub = 0.;
		star2[id].ntry++;
		
		r1 = drand48();
		r1 = exp(-2. * log(r1) / 3.);
		r1 = RSCALE * sqrt(1. / (r1 - 1.));
		
		r1 *= rfac;
		
		isorand(star1[id].r);
		star1[id].r[0] *= r1;
		star1[id].r[1] *= r1;
		star1[id].r[2] *= r1;
		star2[id].rad = r1;
		
		for (i = 0; i < id; i++)
		{
			rx = star1[i].r[0] - star1[id].r[0];
			ry = star1[i].r[1] - star1[id].r[1];
			rz = star1[i].r[2] - star1[id].r[2];
			r1 = sqrt(rx*rx + ry*ry + rz*rz);
			star1[i].U_tmp = star1[id].mass / r1;
			star2[id].U_sub += star1[i].mass / r1;
		};
		r1 = U_sub + star1[id].mass * star2[id].U_sub;
	}
	while (fabs(r1 - UT_sub) > 0.1 * (r1 + UT_sub) / sqrt(1. + id));
	
	for (i = 0; i < id; i++) star2[i].U += star1[i].U_tmp;
}

void find_beta(double *avg, double *beta){
	int i;
	double a1, a2, I1, I2, J1, J2;
	double alo, ah, Ilo, Ih, Jlo, Jh;
	double tmpavg;
	
	if (*avg > 0.75)
	{
		*beta = -0.5;
		*avg = 1.;
		return;
	};
	
	Ilo = I1 = 3. * M_PI / 16.;
	Ih  = I2 = 0.2;
	Jlo = J1 = M_PI / 4;
	Jh  = J2 = 1. / 3.;
	alo = a1 = -0.5;
	ah  = a2 = 0.;
	
	tmpavg = I2 / J2;
	i = 0;
	while (tmpavg > *avg)
	{
		if (!(i & 1))
		{
			a1 += 1.;
			I1 /= 1. + 5. / (2. * a1);
			J1 /= 1. + 3. / (2. * a1);
			Ilo = I2; Ih = I1;
			Jlo = J2; Jh = J1;
			alo = a2; ah = a1;
		}
		else
		{
			a2 += 1.;
			I2 /= 1. + 5. / (2. * a2);
			J2 /= 1. + 3. / (2. * a2);
			Ilo = I1; Ih = I2;
			Jlo = J1; Jh = J2;
			alo = a1; ah = a2;
		};
		tmpavg = Ih / Jh;
		i++;
	};
	
	*beta = alo + (ah - alo) * (*avg - Ilo / Jlo) / (tmpavg - Ilo / Jlo);
	
	if (*beta > 0.) *avg = exp(*beta * log(*beta / (1. + *beta))) / (1. + *beta);
	else *avg = 2. * exp((2.* *beta + 1.) * log((2.* *beta + 1.) / (2.* *beta + 3.))) / (2.* *beta + 3.);
}

void isorand(double *r) {
	double rad;
	
	do
	{
		r[0] = 2. * drand48() - 1.;
		r[1] = 2. * drand48() - 1.;
		r[2] = 2. * drand48() - 1.;
		rad = r[0]*r[0] + r[1]*r[1] + r[2]*r[2];
	}
	while ((rad > 1.) || (rad < 0.0001));
	
	rad = sqrt(rad);
	r[0] /= rad;
	r[1] /= rad;
	r[2] /= rad;
}

int cmpmy(double *x1, double *x2) { //smallest to largest
	if(*x1<*x2) return -1;
	return 1;
}

int cmpmy_reverse(double *x1, double *x2) { //largest to smallest
	if(*x1>*x2) return -1;
	return 1;
}

double generate_profile (int N, double **star, double Rmax, double Mtot, double *p, double *Rh, double D, int symmetry) {
	
	double Mnorm;
	double r, Rmin = 1E-05;
	int steps = 51, i, h, j; //51 steps is more than sufficient (in most cases)
	double r_array[steps], rho_array[steps], sigma_array[steps], M_array[steps], X_array[steps];
	double random[7], x,y,z, ri, vx, vy, vz;
	double r_norm, v_norm;
	
	Mnorm = M(Rmax, p);
	p[1] = Mtot/Mnorm;

	double stepsize;
	stepsize = (log10(Rmax)-log10(Rmin))/(steps-1);
	
	printf("\nIntegrating profile...\n");
	for (i=0;i<steps;i++) {
		r = pow(10.0, log10(Rmin) + stepsize*i);
		//r = Rmin + pow(1.0*i/(steps-1.0),5)*(Rmax-Rmin);
		r_array[i] = r;    //3D radius
		rho_array[i] = rho(r, p); //density(r)
		sigma_array[i] = sigma(r, p); //sigma3D(r)
		M_array[i] = M(r, p); //M(<r)
		X_array[i] = M_array[i]/Mtot;  //M(<r)/Mtot
	}
	
	printf("\n#r [pc]\t\trho_2D(r) [Msun/pc^2]\t\trho(r) [Msun/pc^3]\tsigma(r) [km/s]\t\tM(r) [Msun]\t\tX(r)\n");
	
	for (i=0;i<steps;i++){
		printf ("%f\t%.18f\t%.18f\t%.18f\t%.18f\t%.18f\n", r_array[i], rhoR(r_array[i],p), rho_array[i], sigma_array[i], M_array[i], X_array[i]);
	}
	
	
	//determine half-mass radius
	double Xi, xi;
	i = 0;
	while (X_array[i] < 0.5) {
		i++;
	};
	if (i) {
		xi = r_array[i-1];
		do {
			xi += 0.1;
			Xi = ((xi-r_array[i-1])*X_array[i]+(r_array[i]-xi)*X_array[i-1])/(r_array[i]-r_array[i-1]);
		} while ((Xi < 0.5) && (xi <= r_array[steps-1]));
		*Rh = 1.0*xi;
	} else { 
		*Rh = r_array[0];
	}
	
	
	//draw positions & velocities
	if (D>=3.0) {
	
		for (i = 0; i< N; i++){

			if ((i/1000)*1000 == i) printf("Generating orbit #%i\n", i);
			
			do {
				random[0] = drand48();
			
				j = 0;
				while (X_array[j]<random[0]) {
					j++;
				};
				if (j) {
					xi = random[0];
					Xi = ((xi-X_array[j-1])*r_array[j]+(X_array[j]-xi)*r_array[j-1])/(X_array[j]-X_array[j-1]);
					ri = 1.0*Xi;
				} else { 
					ri = r_array[0];
				}
			
				random[1] = drand48(); 
				random[2] = drand48();
			
				z = (1.0 - 2.0*random[1])*ri;
				x = sqrt(ri*ri-z*z)*cos(TWOPI*random[2]);
				y = sqrt(ri*ri-z*z)*sin(TWOPI*random[2]); 
			} while (sqrt(x*x+y*y+z*z)>Rmax); //reject particles beyond Rmax
		
			double sigma;

			j = 0;
			while (r_array[j]<ri) {
				j++;
			};
			if (j) {
				xi = ri;
				Xi = ((xi-r_array[j-1])*sigma_array[j]+(r_array[j]-xi)*sigma_array[j-1])/(r_array[j]-r_array[j-1]);
				sigma = 1.0*Xi;
			} else { 
				sigma = sigma_array[0];
			}
		
			vx = get_gauss()*sigma;
			vy = get_gauss()*sigma;
			vz = get_gauss()*sigma;
			
			star[i][1] = x;
			star[i][2] = y;
			star[i][3] = z;
			star[i][4] = vx;
			star[i][5] = vy;
			star[i][6] = vz;
		
			//printf("%f\t%f\t%f\t%f\t%f\t%f\n", x,y,z, vx,vy,vz);
		}
	} else {
		fractalize(D, N, star, 1, symmetry);
		for (i=0;i<N;i++) {

			if ((i/1000)*1000 == i) printf("Generating orbit #%i\n", i);
			
			double sigma;
			ri = sqrt(pow(star[i][1],2)+pow(star[i][2],2)+pow(star[i][3],2));

			j = 0;
			while (X_array[j]<ri) {
				j++;
			};
			if (j) {
				xi = ri;
				Xi = ((xi-X_array[j-1])*r_array[j]+(X_array[j]-xi)*r_array[j-1])/(X_array[j]-X_array[j-1]);
				r_norm = 1.0*Xi;
			} else { 
				r_norm = r_array[0];
			}
					
			j = 0;
			while (r_array[j]<r_norm) {
				j++;
			};
			if (j) {
				xi = r_norm;
				Xi = ((xi-r_array[j-1])*sigma_array[j]+(r_array[j]-xi)*sigma_array[j-1])/(r_array[j]-r_array[j-1]);
				sigma = 1.0*Xi;
			} else { 
				sigma = sigma_array[0];
			}
					
			v_norm = sigma;
			for (h=1;h<4;h++) star[i][h] *= r_norm/ri;
			for (h=4;h<7;h++) star[i][h] *= v_norm;
		}
		
	}
	
	
	return 0;
	
}

double dfridr(double (*func)(double, double*), double x, double h, double *err, double *p) {
	double CON = 1.4;
	double CON2 = (CON*CON);
	double SAFE = 2.0;
	
	int i,j;
	double errt,fac,hh,ans = 0;
	double a[10][10];
	
	double arg1, arg2;
	
	hh=h;
	a[1][1]=((*func)(x+hh, p)-(*func)(x-hh, p))/(2.0*hh);
	*err=BIGNUMBER;
	for (i=2;i<=10;i++) {
		hh /= CON;
		a[1][i]=((*func)(x+hh, p)-(*func)(x-hh, p))/(2.0*hh);
		fac=CON2;
		for (j=2;j<=i;j++) {
			a[j][i]=(a[j-1][i]*fac-a[j-1][i-1])/(fac-1.0);
			fac=CON2*fac;
			arg1 = fabs(a[j][i]-a[j-1][i]);
			arg2 = fabs(a[j][i]-a[j-1][i-1]);
			if (arg1 > arg2) errt = arg1;
			else errt = arg2;
			if (errt <= *err) {
				*err=errt;
				ans=a[j][i];
			}
		}
		if (fabs(a[i][i]-a[i-1][i-1]) >= SAFE*(*err)) break;
	}
	
	return ans;
}

double midexp(double (*func)(double, double*), double aa, double bb, int n, double *p) {
	double x,tnm,sum,del,ddel,a,b;
	static double s;
	int it,j;
	
	b=exp(-aa);
	a=0.0;
	if (n == 1) {
		x = 0.5*(a+b);
		return (s=(b-a)*func(-log(x),p)/x);
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += func(-log(x),p)/x;
			x += ddel;
			sum += func(-log(x),p)/x;
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}

double midsql(double (*func)(double, double*), double aa, double bb, int n, double *p) {
	double x,tnm,sum,del,ddel,a,b;
	static double s;
	int it,j;
	
	b=sqrt(bb-aa);
	a=0.0;
	if (n == 1) {
		x = aa+0.5*(a+b)*0.5*(a+b);
		return (s=(b-a)*2.0*x*func(x, p));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += 2.0*x*func(aa+x*x, p);
			x += ddel;
			sum += 2.0*x*func(aa+x*x, p);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}

double midsqu(double (*func)(double, double*), double aa, double bb, int n, double *p) {
	double x,tnm,sum,del,ddel,a,b;
	static double s;
	int it,j;
	
	b=sqrt(bb-aa);
	a=0.0;
	if (n == 1) {
		x = bb-0.5*(a+b)*0.5*(a+b);
		return (s=(b-a)*2.0*x*func(x, p));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += 2.0*x*func(bb-x*x, p);
			x += ddel;
			sum += 2.0*x*func(bb-x*x, p);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}

double midinf(double (*func)(double, double*), double aa, double bb, int n, double *p) {
	double x,tnm,sum,del,ddel,b,a;
	static double s;
	int it,j;
	
	b=1.0/aa;
	a=1.0/bb;
	if (n == 1) {
		x= 0.5*(a+b);
		return (s=(b-a)*func(1.0/x,p)/x/x);
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += func(1.0/x, p)/x/x;
			x += ddel;
			sum += func(1.0/x, p)/x/x;
			x += del;
		}
		return (s=(s+(b-a)*sum/tnm)/3.0);
	}
}

double midpnt(double (*func)(double, double*), double a, double b, int n, double *p) {
	double x,tnm,sum,del,ddel;
	static double s;
	int it,j;
	
	if (n == 1) {
		return (s=(b-a)*func(0.5*(a+b), p));
	} else {
		for(it=1,j=1;j<n-1;j++) it *= 3;
		tnm=it;
		del=(b-a)/(3.0*tnm);
		ddel=del+del;
		x=a+0.5*del;
		sum=0.0;
		for (j=1;j<=it;j++) {
			sum += func(x, p);
			x += ddel;
			sum += func(x, p);
			x += del;
		}
		s=(s+(b-a)*sum/tnm)/3.0;
		return s;
	}
}

double kernel (double x, double *p) {	
	double h = 1.E-06;
	double result, err;	
	result = dfridr(rhoR, x, h, &err, p)/sqrt(x*x-p[0]*p[0]);	
	return result;
}

double rhoR (double x, double *p) {	//user-defined 2d-density profile
	//	return p[1]*pow(1.0+x/p[2]*x/p[2],-0.5*p[3]); //Elson, Fall & Freeman (1987)
	return p[1]*pow(2.0,0.5*(p[3]-p[4]))*pow(x/p[2],-p[4])*pow(1.0+pow(x/p[2],p[5]),-(p[3]-p[4])/p[5]); //Nuker (Lauer et al. 1995) -> Elson, Fall & Freeman (1987) for p[4] = 0 & p[5] = 2
	
}

double rho(double r, double *p) { 
	int j = 0;
	double s,st,ost,os,stemp;
	if (r < p[2]) {
		ost = os = -1.0e30;
		for (j=1;j<=JMAX;j++) {
			if (p[4]) st=midinf(kernel,r,p[2],j,p);
			else st=midpnt(kernel,r,p[2],j,p);
			s=(4.0*st-ost)/3.0;
			if (fabs(s-os) < EPS*fabs(os)) break;
			if (s == 0.0 && os == 0.0 && j > 6) break;
			os=s;
			ost=st;
		}
		stemp = s;
		ost = os = -1.0e30;
		for (j=1;j<=JMAX;j++) {
			st=midinf(kernel,p[2],BIGNUMBER,j,p);
			s=(4.0*st-ost)/3.0;
			if (fabs(s-os) < EPS*fabs(os)) break;
			if (s == 0.0 && os == 0.0 && j > 6) break;
			os=s;
			ost=st;
		}
		s += stemp;
	} else {
		ost = os = -1.0e30;
		for (j=1;j<=JMAX;j++) {
			st=midinf(kernel,r,BIGNUMBER,j,p);
			s=(4.0*st-ost)/3.0;
			if (fabs(s-os) < EPS*fabs(os)) break;
			if (s == 0.0 && os == 0.0 && j > 6) break;
			os=s;
			ost=st;
		}
	}			
	
	return max(s * -1.0/PI, 0.0);  //no negative density!
	
}

double rho_kernel (double x, double *p) {
	double s;
	s = 4.0*PI*x*x*rho(x,p);
	return s;
}

double M(double r, double *p){	
	int j = 0;
	double s,st,ost,os;
	
	ost = os = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		st=midsql(rho_kernel,0.0,r,j,p);
		s=(4.0*st-ost)/3.0;
		if (fabs(s-os) < EPS*fabs(os)) break;
		if (s == 0.0 && os == 0.0 && j > 6) break;
		os=s;
		ost=st;
	}
	return s;
	
}

double sigma_kernel (double x, double *p) {
	double s;
	s = 1.0*rho(x,p)*M(x,p)/(x*x);
	return s;
}

double sigma(double r, double *p) { 
	int j = 0;
	double s,st,ost,os,t;
	
	ost = os = -1.0e30;
	for (j=1;j<=JMAX;j++) {
		st=midexp(sigma_kernel,r,BIGNUMBER,j,p);
		s=(4.0*st-ost)/3.0;
		if (fabs(s-os) < EPS*fabs(os)) break;
		if (s == 0.0 && os == 0.0 && j > 6) break;
		os=s;
		ost=st;
	}
	
	if (rho(r,p)) t = sqrt(G*s/rho(r,p));
	else t = 0.0;
	
	return t;
	
}

double get_gauss(void){
	double random[2],p,q;
	do {
		random[0] = 2.0*drand48()-1.0; 
		random[1] = 2.0*drand48()-1.0; 
		q = random[0]*random[0]+random[1]*random[1];
	} while (q>1.0);
	
	p = sqrt(-2.0*log(q)/q);
	return random[0]*p;
	
}

double fractalize(double D, int N, double **star, int radial, int symmetry) {
    int i, j, h, Nparent, Nparentlow;
	int Ntot = 128.0*pow(8,ceil(log(N)/log(8)));
	int Ntotorg = Ntot;
	double l = 2.0;
	double prob = pow(2.0, D-3.0);
	double scatter;
	if (radial) scatter = 0.01;
	else scatter = 0.1;
	double vx, vy, vz;
	double vscale;
	int subi;
	double morescatter = 0.1; //0.1 looks good
	
	printf("\nFractalizing initial conditions...\n\n");
	
	double **star_temp;   //temporary array for fractalized structure
	star_temp = (double **)calloc(Ntot,sizeof(double *));
	for (j=0;j<Ntot;j++){
		star_temp[j] = (double *)calloc(7,sizeof(double));
		if (star_temp[j] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}		
	
	Nparent = 0; //ur-star
	Nparentlow = Nparent;
	star_temp[Nparent][1] = 0.0;//x
	star_temp[Nparent][2] = 0.0;//y
	star_temp[Nparent][3] = 0.0;//z
	star_temp[Nparent][4] = 0.0;//vx
	star_temp[Nparent][5] = 0.0;//vy
	star_temp[Nparent][6] = 0.0;//vz
	Nparent++;
	
	i=0;
    
	while (Nparent+i*8<Ntot) {
		l /= 2.0;
		i = 0;
		for (;Nparentlow<Nparent;Nparentlow++) {
			subi = 0;
			/*if (drand48()<prob && Nparent+i<Ntot) {
				star_temp[Nparent+i][0] = 1.0;
				star_temp[Nparent+i][1] = star_temp[Nparentlow][1]+l*scatter*get_gauss();
				star_temp[Nparent+i][2] = star_temp[Nparentlow][2]+l*scatter*get_gauss();
				star_temp[Nparent+i][3] = star_temp[Nparentlow][3]+l*scatter*get_gauss();
				v = get_gauss();
				vz = (1.0 - 2.0*drand48())*v;
				vx = sqrt(v*v - vz*vz)*cos(TWOPI*drand48());
				vy = sqrt(v*v - vz*vz)*sin(TWOPI*drand48());
				star_temp[Nparent+i][4] = vx;
				star_temp[Nparent+i][5] = vy;
				star_temp[Nparent+i][6] = vz;
				i++;
				subi++;
			}*/
			if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
				star_temp[Nparent+i][0] = 1.0;
				star_temp[Nparent+i][1] = star_temp[Nparentlow][1]+l/2.0+l*scatter*get_gauss();
				star_temp[Nparent+i][2] = star_temp[Nparentlow][2]+l/2.0+l*scatter*get_gauss();
				star_temp[Nparent+i][3] = star_temp[Nparentlow][3]+l/2.0+l*scatter*get_gauss();
				vx = get_gauss();
				vy = get_gauss();
				vz = get_gauss();
				star_temp[Nparent+i][4] = vx;
				star_temp[Nparent+i][5] = vy;
				star_temp[Nparent+i][6] = vz;
				i++;
				subi++;
			}
			if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
				star_temp[Nparent+i][0] = 1.0;
				star_temp[Nparent+i][1] = star_temp[Nparentlow][1]+l/2.0+l*scatter*get_gauss();
				star_temp[Nparent+i][2] = star_temp[Nparentlow][2]+l/2.0+l*scatter*get_gauss();
				star_temp[Nparent+i][3] = star_temp[Nparentlow][3]-l/2.0+l*scatter*get_gauss();
				vx = get_gauss();
				vy = get_gauss();
				vz = get_gauss();	
				star_temp[Nparent+i][4] = vx;
				star_temp[Nparent+i][5] = vy;
				star_temp[Nparent+i][6] = vz;
				i++;
				subi++;
			}
			if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
				star_temp[Nparent+i][0] = 1.0;
				star_temp[Nparent+i][1] = star_temp[Nparentlow][1]+l/2.0+l*scatter*get_gauss();
				star_temp[Nparent+i][2] = star_temp[Nparentlow][2]-l/2.0+l*scatter*get_gauss();
				star_temp[Nparent+i][3] = star_temp[Nparentlow][3]+l/2.0+l*scatter*get_gauss();
				vx = get_gauss();
				vy = get_gauss();
				vz = get_gauss();
				star_temp[Nparent+i][4] = vx;
				star_temp[Nparent+i][5] = vy;
				star_temp[Nparent+i][6] = vz;
				i++;
				subi++;
			}
			if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
				star_temp[Nparent+i][0] = 1.0;
				star_temp[Nparent+i][1] = star_temp[Nparentlow][1]-l/2.0+l*scatter*get_gauss();
				star_temp[Nparent+i][2] = star_temp[Nparentlow][2]+l/2.0+l*scatter*get_gauss();
				star_temp[Nparent+i][3] = star_temp[Nparentlow][3]+l/2.0+l*scatter*get_gauss();
				vx = get_gauss();
				vy = get_gauss();
				vz = get_gauss();
				star_temp[Nparent+i][4] = vx;
				star_temp[Nparent+i][5] = vy;
				star_temp[Nparent+i][6] = vz;
				i++;
				subi++;
			}
			if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
				star_temp[Nparent+i][0] = 1.0;
				star_temp[Nparent+i][1] = star_temp[Nparentlow][1]+l/2.0+l*scatter*get_gauss();
				star_temp[Nparent+i][2] = star_temp[Nparentlow][2]-l/2.0+l*scatter*get_gauss();
				star_temp[Nparent+i][3] = star_temp[Nparentlow][3]-l/2.0+l*scatter*get_gauss();
				vx = get_gauss();
				vy = get_gauss();
				vz = get_gauss();
				star_temp[Nparent+i][4] = vx;
				star_temp[Nparent+i][5] = vy;
				star_temp[Nparent+i][6] = vz;
				i++;
				subi++;
			}
			if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
				star_temp[Nparent+i][0] = 1.0;
				star_temp[Nparent+i][1] = star_temp[Nparentlow][1]-l/2.0+l*scatter*get_gauss();
				star_temp[Nparent+i][2] = star_temp[Nparentlow][2]+l/2.0+l*scatter*get_gauss();
				star_temp[Nparent+i][3] = star_temp[Nparentlow][3]-l/2.0+l*scatter*get_gauss();
				vx = get_gauss();
				vy = get_gauss();
				vz = get_gauss();
				star_temp[Nparent+i][4] = vx;
				star_temp[Nparent+i][5] = vy;
				star_temp[Nparent+i][6] = vz;
				i++;
				subi++;
			}
			if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
				star_temp[Nparent+i][0] = 1.0;
				star_temp[Nparent+i][1] = star_temp[Nparentlow][1]-l/2.0+l*scatter*get_gauss();
				star_temp[Nparent+i][2] = star_temp[Nparentlow][2]-l/2.0+l*scatter*get_gauss();
				star_temp[Nparent+i][3] = star_temp[Nparentlow][3]+l/2.0+l*scatter*get_gauss();
				vx = get_gauss();
				vy = get_gauss();
				vz = get_gauss();
				star_temp[Nparent+i][4] = vx;
				star_temp[Nparent+i][5] = vy;
				star_temp[Nparent+i][6] = vz;
				i++;
				subi++;
			}
			if ((drand48()<prob && Nparent+i<Ntot) || ((Nparent == 1) && (symmetry))) {
				star_temp[Nparent+i][0] = 1.0;
				star_temp[Nparent+i][1] = star_temp[Nparentlow][1]-l/2.0+l*scatter*get_gauss();
				star_temp[Nparent+i][2] = star_temp[Nparentlow][2]-l/2.0+l*scatter*get_gauss();
				star_temp[Nparent+i][3] = star_temp[Nparentlow][3]-l/2.0+l*scatter*get_gauss();
				vx = get_gauss();
				vy = get_gauss();
				vz = get_gauss();
				star_temp[Nparent+i][4] = vx;
				star_temp[Nparent+i][5] = vy;
				star_temp[Nparent+i][6] = vz;
				i++;
				subi++;
			}
			
			//re-scaling of sub-group
			if (subi) {
				double vx, vy, vz;
				vx = 0.0; vy = 0.0; vz = 0.0;
				vscale = 0.0;
				for (j=Nparent+i-subi;j<Nparent+i;j++) {
					vx += star_temp[j][4];
					vy += star_temp[j][5];
					vz += star_temp[j][6];
				}
				vx /= 1.0*subi;
				vy /= 1.0*subi;
				vz /= 1.0*subi;
				for (j=Nparent+i-subi;j<Nparent+i;j++) {
					star_temp[j][4] -= vx;
					star_temp[j][5] -= vy;
					star_temp[j][6] -= vz;					
				}
				if (subi-1) {
					for (j=Nparent+i-subi;j<Nparent+i;j++) {
						vscale += star_temp[j][4]*star_temp[j][4]+star_temp[j][5]*star_temp[j][5]+star_temp[j][6]*star_temp[j][6];
					}
					vscale = sqrt(vscale/(subi-1));
				} else {
					vscale = 1.0;
				}
				for (j=Nparent+i-subi;j<Nparent+i;j++) {
					star_temp[j][4] = star_temp[j][4]/vscale + star_temp[Nparentlow][4]; //add bulk velocity of parent
					star_temp[j][5] = star_temp[j][5]/vscale + star_temp[Nparentlow][5];
					star_temp[j][6] = star_temp[j][6]/vscale + star_temp[Nparentlow][6];
				}				
			}
			
		}
		Nparent+=i;
	}
	Ntot = Nparent;

	double cmr[7];//centre-of-mass correction
	for (j=0; j<7; j++) cmr[j] = 0.0;
	
	for (j=0; j<Ntot; j++) {
		for (i=1;i<7;i++) 
			cmr[i] += star_temp[j][i]; 
	} 
	
	for (j=0; j<Ntot; j++) {
		for (i=1;i<7;i++)
			star_temp[j][i] -= cmr[i]/Ntot;
	}
		
	i = 0;//randomly select stars from sample with r < 1.0
	for (i=0; i<N; i++) {
		do{ 
			j = drand48()*Ntot;
		} while (!(star_temp[j][0]) || (sqrt(pow(star_temp[j][1],2)+pow(star_temp[j][2],2)+pow(star_temp[j][3],2)) > 1.0)); 
		for (h=1;h<7;h++) star[i][h] = star_temp[j][h];
		star_temp[j][0] = 0.0;
	}

	double r, r_norm, vnorm = 0.0, start[4];
	if (radial) {
		for (i=0;i<N;i++) {
			vnorm += sqrt(pow(star[i][4],2)+pow(star[i][5],2)+pow(star[i][6],2));
		}
		vnorm /= N;
		for (h=4;h<7;h++) star[i][h] /= 0.5*vnorm;
		
		for (i=0;i<N;i++) {
			r = sqrt(pow(star[i][1],2)+pow(star[i][2],2)+pow(star[i][3],2));
			r_norm = pow(r,3);
			do{ 
				for (h=1;h<4;h++) start[h] = star[i][h]*r_norm/r + pow(r_norm/r,3)*morescatter*get_gauss();
			} while (sqrt(pow(start[1],2)+pow(start[2],2)+pow(start[3],2))>1.0);
			for (h=1;h<4;h++) star[i][h] = start[h];
		}
	}
	for (j=0;j<Ntotorg;j++) free (star_temp[j]);
	free(star_temp);
	
	return 0;
}

int get_binaries(int nbin, double **star, double M, double rvir, int pairing, int *N, int adis, double amin, double amax, double Rh, int eigen, int BSE, double epoch, double Z, int remnant, int OBperiods, double msort) {		
	int i, j, k;
	double m1, m2, ecc, abin;
	double eccold, abinold, m1old, m2old;
	double pmat[3][2], rop[2], vop[2], rrel[3], vrel[3];
	double ea, mm, eadot, cosi, inc, peri, node;
	double lP, P;	
	double u1, u2;
	double q, p, x1, x2;	
	double lP1, lP2, lPmean, lPsigma;

	double zpars[20];     //metallicity parameters
	double vkick[2];	 //kick velocity for compact remnants	
	vkick[0] = 0.0;
	vkick[1] = 0.0;

	double vesc;
	if (remnant) {
		vesc = sqrt(2.0*G*M/Rh);
		if (BSE) printf("Keeping only binaries with kick velocities < escape velocity\n");
	} else {
		vesc = 1.0E10;
		if (BSE) printf("Keeping all compact remnants\n");
	}		
	
	for (i=0; i<20; i++) zpars[i] = 0;
	zcnsts_(&Z,zpars);  //get metallicity parameters
	
	if (BSE) printf("\nSetting up binary population with Z = %.4f.\n",Z);
	if (epoch) printf("\nEvolving binary population for %.1f Myr.\n",epoch);
	
	
	for (i=0; i < nbin; i++) {

		do {

			//Specify component masses
			if (BSE) {
				m1 = star[2*i][7]/M;
				m2 = star[2*i+1][7]/M;
			} else {
				m1 = star[2*i][0];
				m2 = star[2*i+1][0];
			}
		
			//Specify semi-major axis
			if (((m1*M>=msort) || (m2*M>=msort)) && (OBperiods)) {
				//derive from Sana & Evans (2011) period distribution for massive binaries
				double lPmin = 0.3, lPmax = 3.5, lPbreak = 1.0; //parameters of Sana & Evans (2011) period distribution in days (eq. 5.1)
				double Fbreak = 0.5; //fraction of binaries with periods below Pbreak
				double xperiod = drand48();
				if (xperiod <= Fbreak)
					lP = xperiod *(lPbreak-lPmin)/Fbreak + lPmin;
				else
					lP = (xperiod - Fbreak)*(lPmax-lPbreak)/(1.0-Fbreak) + lPbreak;

				P = pow(10.0,lP);//days
				P /= 365.25;//yr
				
				abin = pow((m1+m2)*M*P*P,(1.0/3.0));//AU
				abin /= 206264.806;//pc
				abin /= rvir;//Nbody units				
			} else if (adis == 0) {
				//flat semi-major axis distribution
				if (!i) printf("\nApplying flat semi-major axis distribution with amin = %g and amax = %g.\n", amin, amax);
				if (!i) amin /= rvir;
				if (!i) amax /= rvir;
				abin = amin+drand48()*(amax-amin);
			} else if (adis == 1) {
				//derive from Kroupa (1995) period distribution
				if (!i) printf("\nDeriving semi-major axis distribution from Kroupa (1995) period distribution.\n");
				double Pmin = 10, delta = 45, eta = 2.5; //parameters of Kroupa (1995) period distribution
				do {
					lP = log10(Pmin) + sqrt(delta*(exp(2.0*drand48()/eta) - 1.0));
				} while (lP > 8.43);
			
				P = pow(10.0,lP);//days
				P /= 365.25;//yr
			
				abin = pow((m1+m2)*M*P*P,(1.0/3.0));//AU
				abin /= 206264.806;//pc
				abin /= rvir;//Nbody units
			} else {
				//derive from Duquennoy & Mayor (1991) period distribution
				lPmean = 4.8;	//mean of Duquennoy & Mayor (1991) period distribution [log days]
				lPsigma = 2.3;	//full width half maximum of Duquennoy & Mayor (1991) period distribution [log days]
				if (!i) printf("\nDeriving semi-major axis distribution from Duquennoy & Mayor (1991) period distribution.\n");
				do {
					//generate two random numbers in the interval [-1,1]
					u1 = 2.0*drand48()-1.0;
					u2 = 2.0*drand48()-1.0;
					
					//combine the two random numbers
					q = u1*u1 + u2*u2;
				} while (q > 1.0);

				p = sqrt(-2.0*log(q)/q);
				x1 = u1*p;
				x2 = u2*p;
				lP1 = lPsigma*x1 + lPmean;
				lP2 = lPsigma*x2 + lPmean;

				P = pow(10,lP1);//days
				P /= 365.25;//yr
			
				abin = pow((m1+m2)*M*P*P,(1.0/3.0));//AU
				abin /= 206264.806;//pc
				abin /= rvir;//Nbody units			
			}

			if (!i && OBperiods && msort) printf("\nDeriving semi-major axis distribution for binaries with primary masses > %.3f Msun from Sana & Evans (2011) period distribution.\n",msort);
		
			//Specify eccentricity distribution
			if (!i && OBperiods && msort) printf("\nApplying thermal eccentricity distribution for low-mass systems and Sana & Evans (2011) eccentricity distribution for high-mass systems.\n");
			else if (!i) printf("\nApplying thermal eccentricity distribution.\n");	
				
			if (((m1*M>=msort) || (m2*M>=msort)) && (OBperiods)) {
				double k1, k2;
				double elim = 0.8; //maximum eccentricity for high-mass binaries
				double fcirc = 0.3; //fraction of circular orbits among high-mass binaries
				k2 = (fcirc*exp(0.5*elim)-1.0)/(exp(0.5*elim)-1.0);
				k1 = fcirc-k2;
				ecc = drand48();
				if (ecc < fcirc) ecc = 0.0;
				else {
					ecc = 2.0*log((ecc-k2)/k1);
				}
			} else {
				ecc = sqrt(drand48());   // Thermal distribution f(e)=2e 				
			}
			//ecc = 0.0;   // all circular 
		
			
			
			//Apply Kroupa (1995) eigenevolution
			eccold = ecc;
			abinold = abin;
			m1old = m1;
			m2old = m2;
			if (eigen) {
				if (!i) printf("\nApplying Kroupa (1995) eigenevolution for short-period binaries\n");
				m1*=M;//temporary re-scaling
				m2*=M;
				abin*=rvir;
				if (!(OBperiods && ((m1>=msort) || (m2>=msort)))) {
					if (m1>=m2) eigenevolution(&m1, &m2, &ecc, &abin);
					else  eigenevolution(&m2, &m1, &ecc, &abin);
				}
				m1/=M;
				m2/=M;
				abin/=rvir;
			}
		
		
			//Apply Binary Star Evolution (Hurley, Tout & Pols 2002)
			if (BSE) {
				int kw[2] = {1, 1};  //stellar type
				double mass[2];  //initial mass
				double mt[2];    //actual mass
				double r[2] = {0.0, 0.0};     //radius
				double lum[2] = {0.0, 0.0};   //luminosity
				double mc[2] = {0.0, 0.0};    //core mass
				double rc[2] = {0.0, 0.0};    //core radius
				double menv[2] = {0.0, 0.0};  //envelope mass
				double renv[2] = {0.0, 0.0};  //envelope radius
				double ospin[2] = {0.0, 0.0};		//spin
				double epoch1[2] = {0.0, 0.0};    //time spent in current evolutionary state
				double tms[2] = {0.0, 0.0};   //main-sequence lifetime
				double tphys = 0.0;       //initial age
				double tphysf = epoch;//final age
				double dtp = epoch+1; //data store value, if dtp>tphys no data will be stored
				double z = Z;      //metallicity
						
				mass[0] = m1*M;  //initial mass of primary
				mass[1] = m2*M; //initial mass of secondary
				mt[0] = mass[0];    //actual mass
				mt[1] = mass[1];    //actual mass
				vkick[0] = 0.0;
				vkick[1] = 0.0;
			
				abin *= rvir; //pc
				abin *= 206264.806; //AU
				P = sqrt(pow(abin, 3.0)/(m1+m2));//yr
				P *= 365.25;//days
			
				evolv2_(kw,mass,mt,r,lum,mc,rc,menv,renv,ospin,epoch1,tms,&tphys,&tphysf,&dtp,&z,zpars,&P,&eccold, vkick);

				P /= 365.25;//yr
			
				abin = pow((m1+m2)*P*P,(1.0/3.0));//AU
				abin /= 206264.806;//pc
				abin /= rvir;//Nbody units			
			
				m1 = mt[0]/M; //Nbody units
				m2 = mt[1]/M; //Nbody units

				star[2*i][8] = kw[0];
				star[2*i+1][8] = kw[1];
				star[2*i][9] = epoch1[0];
				star[2*i+1][9] = epoch1[1];
				star[2*i][10] = ospin[0];
				star[2*i+1][10] = ospin[1];
				star[2*i][11] = r[0];
				star[2*i+1][11] = r[1];
				star[2*i][12] = lum[0];
				star[2*i+1][12] = lum[1];
			}

			//pos & vel in binary frame
			ea = rtnewt(ecc, drand48());
			rop[0] = abin*(cos(ea) - ecc);
			rop[1] = abin*sqrt(1.0-ecc*ecc)*sin(ea);
			
			mm = sqrt((m1+m2)/pow(abin,3));
			eadot = mm/(1.0 - ecc*cos(ea));
			vop[0] = -abin*sin(ea)*eadot;
			vop[1] = abin*sqrt(1.0-ecc*ecc)*cos(ea)*eadot;
				
			//Convert to cluster frame 
			cosi = 2.0*drand48()-1.0;
			inc = acos(cosi);
			node = 2.0*PI*drand48();
			peri = 2.0*PI*drand48();
		
			pmat[0][0] = cos(peri)*cos(node) - sin(peri)*sin(node)*cosi;
			pmat[1][0] = cos(peri)*sin(node) + sin(peri)*cos(node)*cosi;
			pmat[2][0] = sin(peri)*sin(inc);
			pmat[0][1] = -sin(peri)*cos(node) - cos(peri)*sin(node)*cosi;
			pmat[1][1] = -sin(peri)*sin(node) + cos(peri)*cos(node)*cosi;
			pmat[2][1] = cos(peri)*sin(inc);
		
			for (j=0;j<3;j++) {
				rrel[j] = pmat[j][0]*rop[0] + pmat[j][1]*rop[1];
				vrel[j] = pmat[j][0]*vop[0] + pmat[j][1]*vop[1];
			}
		
			for (j=0;j<3;j++) {
				star[2*i+1][j+1] = star[2*i][j+1] + m1/(m1+m2)*rrel[j];		 //Star2 pos
				star[2*i+1][j+4] = star[2*i][j+4] + m1/(m1+m2)*vrel[j];      //Star2 vel
				star[2*i][j+1]  -= m2/(m1+m2)*rrel[j];                       //Star1 pos
				star[2*i][j+4]  -= m2/(m1+m2)*vrel[j];                       //Star1 vel
			}

			star[2*i][0] = m1;
			star[2*i+1][0] = m2;

		} while ((sqrt(vkick[0]) > vesc) && (sqrt(vkick[1]) > vesc));

		//printf("%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", m1/(m1+m2)*rrel[0], m1/(m1+m2)*rrel[1], m1/(m1+m2)*rrel[2], m1/(m1+m2)*vrel[0], m1/(m1+m2)*vrel[1], m1/(m1+m2)*vrel[2], -m2/(m1+m2)*rrel[0], -m2/(m1+m2)*rrel[1], -m2/(m1+m2)*rrel[2], -m2/(m1+m2)*vrel[0], -m2/(m1+m2)*vrel[1], -m2/(m1+m2)*vrel[2]);
		//printf("%i\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n",i,m1*M, m2*M, P, abin*rvir, ecc, m1old*M, m2old*M, abinold*rvir, eccold);
		//printf("%i\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\t%.10lf\n",i,star[2*i][4], star[2*i+1][4], star[2*i][1], star[2*i+1][1], star[2*i][1]-star[2*i+1][1], m2old*M, abinold*rvir, eccold);

		if (!star[2*i][0] || !star[2*i+1][0]) {
			nbin--; //one binary less because one component went supernova
			
			double startemp[15]; //temporarily save surviving companion
			if (!star[2*i][0]) {
				for (k=0;k<15;k++)  startemp[k] = star[2*i+1][k];
			} else {
				for (k=0;k<15;k++)  startemp[k] = star[2*i][k];
			}
			
			for (j=2*i+2;j<*N;j++) {
				for (k=0;k<15;k++) star[j-2][k] = star[j][k]; //move all remaining stars two positions up in the array
			}
			
			*N = *N-1; //reduce total number of stars by one, i.e. remove massless supernova remnant from computations

			for (k=0;k<15;k++) star[*N-1][k] = startemp[k]; //make surviving companion the last particle in array

			i--; //reduce binary counter
		} 
 
	}
	
	return 0;

}
	
void shellsort(double **array, int N, int k) {//largest up
	int i,j,l,n;
	N = N-1;
	double swap[k];
	//guess distance n
	for (n = 1; n <= N/9; n = 3*n+1);
	for (; n > 0; n /= 3) {
		for (i = n; i <= N; i++) {
			for (l=0; l<k; l++) swap[l] = array[i][l];
			for (j = i; ((j >= n) && (array[j-n][0] < swap[0])); j -= n) {
				for (l=0; l<k; l++) array[j][l] = array[j-n][l];
			}
			for (l=0; l<k; l++) array[j][l] = swap[l];
		}
	}
}

void shellsort_reverse(double **array, int N, int k) {//smallest up
	int i,j,l,n;
	N = N-1;
	double swap[k];
	//guess distance n
	for (n = 1; n <= N/9; n = 3*n+1);
	for (; n > 0; n /= 3) {
		for (i = n; i <= N; i++) {
			for (l=0; l<k; l++) swap[l] = array[i][l];
			for (j = i; ((j >= n) && (array[j-n][0] > swap[0])); j -= n) {
				for (l=0; l<k; l++) array[j][l] = array[j-n][l];
			}
			for (l=0; l<k; l++) array[j][l] = swap[l];
		}
	}
}

void shellsort_1d(double *array, int N) {//largest up
	int i,j,n;
	N = N-1;
	double swap;
	//guess distance n
	for (n = 1; n <= N/9; n = 3*n+1);
	for (; n > 0; n /= 3) {
		for (i = n; i <= N; i++) {
			swap = array[i];
			for (j = i; ((j >= n) && (array[j-n] < swap)); j -= n) {
				array[j] = array[j-n];
			}
			array[j] = swap;
		}
	}
}

void shellsort_reverse_1d(double *array, int N) {//smallest up
	int i,j,n;
	N = N-1;
	double swap;
	//guess distance n
	for (n = 1; n <= N/9; n = 3*n+1);
	for (; n > 0; n /= 3) {
		for (i = n; i <= N; i++) {
			swap = array[i];
			for (j = i; ((j >= n) && (array[j-n] > swap)); j -= n) {
				array[j] = array[j-n];
			}
			array[j] = swap;
		}
	}
}

int order(double **star, int N, double M, double msort, int pairing){
	int i,j;
	int Nhighmass;
	int columns = 15;
	double **star_temp;
	star_temp = (double **)calloc(N,sizeof(double *));
	for (j=0;j<N;j++){
		star_temp[j] = (double *)calloc(columns,sizeof(double));
		star_temp[j][0] = 0.0;
		if (star_temp[j] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}
	
	double **masses;
	masses = (double **)calloc(N,sizeof(double *));
	for (j=0;j<N;j++){
		masses[j] = (double *)calloc(2,sizeof(double));
		if (masses[j] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}
		
	//temporary scaling to astrophysical units
	for (i=0;i<N;i++) {
		//star[i][0] *= M;
		//masses[i][0] = star[i][0];
		masses[i][0] = star[i][7];
		masses[i][1] = i;
	}

	shellsort(masses, N, 2);	
	
	j = 0;
	for (i=0;i<N;i++) {
		if (masses[i][0] >= msort) {//copying to front of temporary array if massive
			star_temp[j][0] = star[(int) masses[i][1]][0];
			star_temp[j][1] = star[(int) masses[i][1]][1];
			star_temp[j][2] = star[(int) masses[i][1]][2];
			star_temp[j][3] = star[(int) masses[i][1]][3];
			star_temp[j][4] = star[(int) masses[i][1]][4];
			star_temp[j][5] = star[(int) masses[i][1]][5];
			star_temp[j][6] = star[(int) masses[i][1]][6];
			star_temp[j][7] = star[(int) masses[i][1]][7];
			star_temp[j][8] = star[(int) masses[i][1]][8];
			star_temp[j][9] = star[(int) masses[i][1]][9];
			star_temp[j][10] = star[(int) masses[i][1]][10];
			star_temp[j][11] = star[(int) masses[i][1]][11];
			star_temp[j][12] = star[(int) masses[i][1]][12];
			star_temp[j][13] = star[(int) masses[i][1]][13];
			star_temp[j][14] = star[(int) masses[i][1]][14];
			j++;
		}
	}
	if (msort) Nhighmass = j;
	
	for (i=0;i<N;i++) {
		if (masses[i][0] < msort) {//copying to random position in the back of temporary array
			do {
				j = drand48()*N;
			} while (star_temp[j][0]);
			star_temp[j][0] = star[(int) masses[i][1]][0];
			star_temp[j][1] = star[(int) masses[i][1]][1];
			star_temp[j][2] = star[(int) masses[i][1]][2];
			star_temp[j][3] = star[(int) masses[i][1]][3];
			star_temp[j][4] = star[(int) masses[i][1]][4];
			star_temp[j][5] = star[(int) masses[i][1]][5];
			star_temp[j][6] = star[(int) masses[i][1]][6];
			star_temp[j][7] = star[(int) masses[i][1]][7];
			star_temp[j][8] = star[(int) masses[i][1]][8];
			star_temp[j][9] = star[(int) masses[i][1]][9];
			star_temp[j][10] = star[(int) masses[i][1]][10];
			star_temp[j][11] = star[(int) masses[i][1]][11];
			star_temp[j][12] = star[(int) masses[i][1]][12];
			star_temp[j][13] = star[(int) masses[i][1]][13];
			star_temp[j][14] = star[(int) masses[i][1]][14];
		}
	}
	
		
	//copying back to original array and scaling back to Nbody units
	for (i=0;i<N;i++) {
		//star[i][0] = star_temp[i][0]/M;
		star[i][0] = star_temp[i][0];
		star[i][1] = star_temp[i][1];
		star[i][2] = star_temp[i][2];
		star[i][3] = star_temp[i][3];
		star[i][4] = star_temp[i][4];
		star[i][5] = star_temp[i][5];
		star[i][6] = star_temp[i][6];
		star[i][7] = star_temp[i][7];
		star[i][8] = star_temp[i][8];
		star[i][9] = star_temp[i][9];
		star[i][10] = star_temp[i][10];
		star[i][11] = star_temp[i][11];
		star[i][12] = star_temp[i][12];
		star[i][13] = star_temp[i][13];
		star[i][14] = star_temp[i][14];
	}

	for (j=0;j<N;j++) free (masses[j]);
	free(masses);

	for (j=0;j<N;j++) free (star_temp[j]);
	free(star_temp);

	if (pairing==2) segregate(star, Nhighmass, 0.0);
	
	return 0;
}

int segregate(double **star, int N, double S){
	int i,j;
	
	int columns = 15;
	double **star_temp;
	star_temp = (double **)calloc(N,sizeof(double *));
	for (j=0;j<N;j++){
		star_temp[j] = (double *)calloc(columns,sizeof(double));
		star_temp[j][0] = 0.0;
		if (star_temp[j] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}
	
	double **masses;
	masses = (double **)calloc(N,sizeof(double *));
	for (j=0;j<N;j++){
		masses[j] = (double *)calloc(2,sizeof(double));
		if (masses[j] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}
	
	for (i=0;i<N;i++) {
		masses[i][0] = star[i][0];
		masses[i][1] = i;
	}
	
	shellsort(masses, N, 2);	

	for (i=0;i<N;i++) star_temp[i][0] = 0.0;
	
	j = 0;
	int Ntemp,l;
	Ntemp = N;
	for (i=0;i<N;i++) {
		j = 1.0*(1.0-pow(drand48(),1.0-S))*Ntemp;
		l=-1;
		do {
			l++;
			if (star_temp[l][0]) {
				j++;
			}
		} while (l<j);
		star_temp[j][0] = star[(int) masses[i][1]][0];
		star_temp[j][1] = star[(int) masses[i][1]][1];
		star_temp[j][2] = star[(int) masses[i][1]][2];
		star_temp[j][3] = star[(int) masses[i][1]][3];
		star_temp[j][4] = star[(int) masses[i][1]][4];
		star_temp[j][5] = star[(int) masses[i][1]][5];
		star_temp[j][6] = star[(int) masses[i][1]][6];
		star_temp[j][7] = star[(int) masses[i][1]][7];
		star_temp[j][8] = star[(int) masses[i][1]][8];
		star_temp[j][9] = star[(int) masses[i][1]][9];
		star_temp[j][10] = star[(int) masses[i][1]][10];
		star_temp[j][11] = star[(int) masses[i][1]][11];
		star_temp[j][12] = star[(int) masses[i][1]][12];
		star_temp[j][13] = star[(int) masses[i][1]][13];
		star_temp[j][14] = star[(int) masses[i][1]][14];
		Ntemp--;
	}	
	
	//copying back to original array
	for (i=0;i<N;i++) {
		star[i][0] = star_temp[i][0];
		star[i][1] = star_temp[i][1];
		star[i][2] = star_temp[i][2];
		star[i][3] = star_temp[i][3];
		star[i][4] = star_temp[i][4];
		star[i][5] = star_temp[i][5];
		star[i][6] = star_temp[i][6];
		star[i][7] = star_temp[i][7];
		star[i][8] = star_temp[i][8];
		star[i][9] = star_temp[i][9];
		star[i][10] = star_temp[i][10];
		star[i][11] = star_temp[i][11];
		star[i][12] = star_temp[i][12];
		star[i][13] = star_temp[i][13];
		star[i][14] = star_temp[i][14];
	}
	
	for (j=0;j<N;j++) free (masses[j]);
	free(masses);
	
	for (j=0;j<N;j++) free (star_temp[j]);
	free(star_temp);
	
	return 0;
}

int energy_order(double **star, int N, int Nstars){
	int i,j;
	
	int columns = 15;
	double **star_temp;
	star_temp = (double **)calloc(N,sizeof(double *));
	for (j=0;j<N;j++){
		star_temp[j] = (double *)calloc(columns,sizeof(double));
		star_temp[j][0] = 0.0;
		if (star_temp[j] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}
	
	double **energies;
	energies = (double **)calloc(N,sizeof(double *));
	for (j=0;j<N;j++){
		energies[j] = (double *)calloc(2,sizeof(double));
		energies[j][0] = 0.0;
		if (energies[j] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}
	
	for (i=0;i<N;i++) {
		energies[i][0] = 0.5/Nstars*(star[i][4]*star[i][4]+star[i][5]*star[i][5]+star[i][6]*star[i][6])-1.0/sqrt(star[i][1]*star[i][1]+star[i][2]*star[i][2]+star[i][3]*star[i][3]);//Ekin + Epot = 0.5*1/N*v^2-1/r
		energies[i][1] = i;
	}
	shellsort_reverse(energies, N, 2);	

	for (i=0;i<N;i++) {//copying to random position in temporary array
		star_temp[i][0] = star[(int) energies[i][1]][0];
		star_temp[i][1] = star[(int) energies[i][1]][1];
		star_temp[i][2] = star[(int) energies[i][1]][2];
		star_temp[i][3] = star[(int) energies[i][1]][3];
		star_temp[i][4] = star[(int) energies[i][1]][4];
		star_temp[i][5] = star[(int) energies[i][1]][5];
		star_temp[i][6] = star[(int) energies[i][1]][6];
		star_temp[i][7] = star[(int) energies[i][1]][7];
		star_temp[i][8] = star[(int) energies[i][1]][8];
		star_temp[i][9] = star[(int) energies[i][1]][9];
		star_temp[i][10] = star[(int) energies[i][1]][10];
		star_temp[i][11] = star[(int) energies[i][1]][11];
		star_temp[i][12] = star[(int) energies[i][1]][12];
	}
	
	
	//copying back to original array
	for (i=0;i<N;i++) {
		star[i][0] = star_temp[i][0];
		star[i][1] = star_temp[i][1];
		star[i][2] = star_temp[i][2];
		star[i][3] = star_temp[i][3];
		star[i][4] = star_temp[i][4];
		star[i][5] = star_temp[i][5];
		star[i][6] = star_temp[i][6];
		star[i][7] = star_temp[i][7];
		star[i][8] = star_temp[i][8];
		star[i][9] = star_temp[i][9];
		star[i][10] = star_temp[i][10];
		star[i][11] = star_temp[i][11];
		star[i][12] = star_temp[i][12];
		star[i][13] = star_temp[i][13];
		star[i][14] = star_temp[i][14];
	}
	
	for (j=0;j<N;j++) free (energies[j]);
	free(energies);

	for (j=0;j<N;j++) free (star_temp[j]);
	free(star_temp);
	
	return 0;
}

int randomize(double **star, int N){
	int i,j;
	
	int columns = 15;
	double **star_temp;
	star_temp = (double **)calloc(N,sizeof(double *));
	for (j=0;j<N;j++){
		star_temp[j] = (double *)calloc(columns,sizeof(double));
		star_temp[j][0] = 0.0;
		if (star_temp[j] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}
	
	for (i=0;i<N;i++) {//copying randomized to temporary array
			do {
				j = drand48()*N;
			} while (star_temp[j][0]);
			star_temp[j][0] = star[i][0];
			star_temp[j][1] = star[i][1];
			star_temp[j][2] = star[i][2];
			star_temp[j][3] = star[i][3];
			star_temp[j][4] = star[i][4];
			star_temp[j][5] = star[i][5];
			star_temp[j][6] = star[i][6];
			star_temp[j][7] = star[i][7];
			star_temp[j][8] = star[i][8];
			star_temp[j][9] = star[i][9];
			star_temp[j][10] = star[i][10];
			star_temp[j][11] = star[i][11];
			star_temp[j][12] = star[i][12];
			star_temp[j][13] = star[i][13];
			star_temp[j][14] = star[i][14];
	}
	
	//copying back to original array
	for (i=0;i<N;i++) {
		star[i][0] = star_temp[i][0];
		star[i][1] = star_temp[i][1];
		star[i][2] = star_temp[i][2];
		star[i][3] = star_temp[i][3];
		star[i][4] = star_temp[i][4];
		star[i][5] = star_temp[i][5];
		star[i][6] = star_temp[i][6];
		star[i][7] = star_temp[i][7];
		star[i][8] = star_temp[i][8];
		star[i][9] = star_temp[i][9];
		star[i][10] = star_temp[i][10];
		star[i][11] = star_temp[i][11];
		star[i][12] = star_temp[i][12];
		star[i][13] = star_temp[i][13];
		star[i][14] = star_temp[i][14];
	}
	
	for (j=0;j<N;j++) free (star_temp[j]);
	free(star_temp);
	
	return 0;
}

double rtnewt (double ecc, double ma) { 
	
	double x1,x2,xacc,rtnewt,f,df,dx;
	int j,jmax;
	
	x1 = 0;
	x2 = 2*PI;
	xacc = 1E-6;
	jmax = 20;
	ma = 2*PI*ma; 
	
	rtnewt=.5*(x1+x2);
	for (j=1;j<=jmax;j++) {
		f = ma - rtnewt + ecc*sin(rtnewt);
		df = -1 + ecc*cos(rtnewt);
		dx=f/df;
		rtnewt=rtnewt-dx;
		if ((x1-rtnewt)*(rtnewt-x2)<0) 
			printf("jumped out of brackets\n");
		if(abs(dx)<xacc) return(rtnewt);
	}
	printf("RTNEWT exceeding maximum iterations\n");
	exit(-1); 
}

int eigenevolution(double *m1, double *m2, double *ecc, double *abin){
	double alpha = 28.0;
	double beta = 0.75;
	double mtot,lper,lperi,ecci,mtoti,r0,rperi,qold,qnew;
	
	*abin *= PARSEC;
	mtot = *m1+*m2;
	
	lperi = sqrt(pow(*abin,3)/mtot*4.0*PI*PI/GBIN);
	
	ecci = *ecc;
	mtoti = mtot;
	
	/* Circularisation */
	r0 = alpha*RSUN;
	rperi = *abin*(1.0-*ecc);
	alpha = -pow((r0/rperi),beta);
	
	if (*ecc > 0) {
		*ecc = exp(alpha+log(*ecc));
	}
	
	/* pre-ms mass-transfer */
	qold = *m1/ *m2;
	if (qold > 1.0) qold = 1.0/qold;
	alpha = -alpha;
	if (alpha > 1.0) {
		qnew = 1.0;
	} else {
		qnew = qold + (1.0-qold)*alpha;
	}
	
	*m1 = max(*m1,*m2);
	*m2 = qnew * *m1;
	
	mtot = *m1 + *m2;
	
	lper = lperi*pow((1.0-ecci)/(1.0-*ecc),1.5);
	lper = lper*sqrt(mtoti/mtot);
	
	*abin = pow(mtot*lper*lper*GBIN/4.0/PI/PI,1.0/3.0);
	*abin /= PARSEC;

	return 0;
}

int radial_profile(double **star, int N, double rvir, double M, int create_radial_profile, int create_cumulative_profile, int code, int *NNBMAX, double *RS0, double *Rh2D, double *Rh3D, int NNBMAX_NBODY6) {
	int i, j;
	*Rh2D = 0.0;
	*Rh3D = 0.0;
	double Mtemp;
	double **rarray;
	rarray = (double **)calloc(N,sizeof(double *));
	for (j=0;j<N;j++){
		rarray[j] = (double *)calloc(3,sizeof(double));
		if (rarray[j] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	} 
	double **rarray2D;
	rarray2D = (double **)calloc(N,sizeof(double *));
	for (j=0;j<N;j++){
		rarray2D[j] = (double *)calloc(3,sizeof(double));
		if (rarray2D[j] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	} 
	
	for (j=0; j<N; j++) {
		rarray[j][0] = rvir*sqrt(star[j][1]*star[j][1]+star[j][2]*star[j][2]+star[j][3]*star[j][3]);
		rarray[j][1] = star[j][0]*M;
		rarray[j][2] = star[j][12];
		rarray2D[j][0] = rvir*sqrt(star[j][1]*star[j][1]+star[j][2]*star[j][2]);
		rarray2D[j][1] = star[j][0]*M;
		rarray2D[j][2] = star[j][12];
	}
	shellsort_reverse(rarray, N, 3);	
	shellsort_reverse(rarray2D, N, 3);	
	
	int noofradialbins = 20;
	double rprofile[noofradialbins][9];
	double rmax = 100.0; //pc
	double rmin = 0.1;
	double stepsize;
	stepsize = (log10(rmax)-log10(rmin))/(noofradialbins-1);
		
	for (j=0;j<noofradialbins;j++) {
		rprofile[j][0] = pow(10.0, log10(rmin) + stepsize*j); //radius
		//rprofile[j][0] = 1.0*exp(log(rmax)*(j+1)/noofradialbins)-1.0; //radius
		rprofile[j][2] = 0; //number
		rprofile[j][3] = 0; //mass
		rprofile[j][4] = 0; //luminosity
		if (j == 0)
			rprofile[j][1] =  4.0/3.0*PI*pow(rprofile[j][0],3); //volume
		else
			rprofile[j][1] = 4.0/3.0*PI*(pow(rprofile[j][0],3) - pow(rprofile[j-1][0],3)); //volume
		rprofile[j][6] = 0; //number (2D)
		rprofile[j][7] = 0; //mass (2D)
		rprofile[j][8] = 0; //luminosity (2D)
		if (j == 0)
			rprofile[j][5] =  PI*pow(rprofile[j][0],2); //area
		else
			rprofile[j][5] = PI*(pow(rprofile[j][0],2) - pow(rprofile[j-1][0],2)); //area
	}
	
	j = 0; i = 0; Mtemp = 0.0;
	while ((j < noofradialbins) && (i < N)) {
		if (rarray[i][0] < rprofile[j][0]) {//3D binning in astrophysical units
			rprofile[j][2] += 1.0;
			rprofile[j][3] += rarray[i][1];
			rprofile[j][4] += rarray[i][2];
			Mtemp += rarray[i][1];
			if ((Mtemp>=0.5*M) && !(*Rh3D)) *Rh3D = rarray[i][0];
			i++;
		} else {
			j++;
		}
	}
	
	j = 0; i = 0; Mtemp = 0.0;
	while ((j < noofradialbins) && (i < N)) {
		if (rarray2D[i][0] < rprofile[j][0]) {//2D binning in astrophysical units
			rprofile[j][6] += 1.0;
			rprofile[j][7] += rarray2D[i][1];
			rprofile[j][8] += rarray2D[i][2];
			Mtemp += rarray2D[i][1];
			if ((Mtemp>=0.5*M) && !(*Rh2D)) *Rh2D = rarray2D[i][0];
			i++;
		} else {
			j++;
		}
	}
	
	if (create_radial_profile) {
		printf("\nRadial density profile:\n\n#  R [pc]   N [1/pc^3]   M [Msun/pc^3]   L [Lsun/pc^3]   N [1/pc^2]   M [Msun/pc^2]   L [Lsun/pc^2] \n");
		for (j=0;j<noofradialbins;j++) printf("%9.4f  %11.4f  %14.4f  %14.4f  %11.4f  %14.4f %15.4f\n",rprofile[j][0],rprofile[j][2]/rprofile[j][1],rprofile[j][3]/rprofile[j][1],rprofile[j][4]/rprofile[j][1],rprofile[j][6]/rprofile[j][5],rprofile[j][7]/rprofile[j][5],rprofile[j][8]/rprofile[j][5]); //print radial density profile to screen
	}
	if (create_cumulative_profile) {
		double ntemp = 0.0;
		double mtemp = 0.0;
		double ltemp = 0.0;
		double ntemp2D = 0.0;
		double mtemp2D = 0.0;
		double ltemp2D = 0.0;
		double nmax, mmax, lmax, nhalf, mhalf, lhalf;
		double nhalf2D, mhalf2D, lhalf2D;
		printf("\nCumulative profile:\n\n#  R [pc]         N     M [Msun]     L [Lsun]    N (2D)    M(2D) [Msun]    L(2D) [Lsun]\n");
		for (j=0;j<noofradialbins;j++) {
			ntemp += rprofile[j][2];
			mtemp += rprofile[j][3];
			ltemp += rprofile[j][4];
			ntemp2D += rprofile[j][6];
			mtemp2D += rprofile[j][7];
			ltemp2D += rprofile[j][8];
			printf("%9.4f  %8.0f  %11.3f  %11.3f  %8.0f  %14.4f  %14.4f\n",rprofile[j][0],ntemp,mtemp,ltemp,ntemp2D,mtemp2D,ltemp2D); //print cumulative profile to screen
		}
		nmax = ntemp;
		mmax = mtemp;
		lmax = ltemp;
		nhalf = mhalf = lhalf = 0.0;
		nhalf2D = mhalf2D = lhalf2D = 0.0;
		ntemp = mtemp = ltemp = 0.0;
		ntemp2D = mtemp2D = ltemp2D = 0.0;
		for (j=0;j<N;j++) {
			ntemp += 1.0;
			mtemp += rarray[j][1];
			ltemp += rarray[j][2];
			if (!(nhalf) && (ntemp>=0.5*nmax)) nhalf = rarray[j][0];
			if (!(mhalf) && (mtemp>=0.5*mmax)) mhalf = rarray[j][0];
			if (!(lhalf) && (ltemp>=0.5*lmax)) lhalf = rarray[j][0];
			ntemp2D += 1.0;
			mtemp2D += rarray2D[j][1];
			ltemp2D += rarray2D[j][2];
			if (!(nhalf2D) && (ntemp2D>=0.5*nmax)) nhalf2D = rarray2D[j][0];
			if (!(mhalf2D) && (mtemp2D>=0.5*mmax)) mhalf2D = rarray2D[j][0];
			if (!(lhalf2D) && (ltemp2D>=0.5*lmax)) lhalf2D = rarray2D[j][0];
		}
		printf("\nHalf-number radius = %.4f pc (%.4f pc, 2D)\nHalf-mass radius = %.4f pc (%.4f pc, 2D)\nHalf-light radius = %.4f pc (%.4f pc, 2D)\n",nhalf,nhalf2D,mhalf,mhalf2D,lhalf,lhalf2D); //print radii
		
	}
	
	
	//estimate NNBMAX and RS0 (Nbody6 only)
	if ((code == 0) || (code == 2)) {
		*NNBMAX = 2.0*sqrt(N);
		if (*NNBMAX < 30) *NNBMAX = 30;
		if (N<=*NNBMAX) *NNBMAX = 0.5*N;
		if (*NNBMAX > NNBMAX_NBODY6) *NNBMAX = NNBMAX_NBODY6;
		*RS0 = rarray[*NNBMAX][0];
		printf("\nEstimating appropriate NNBMAX = %i and RS0 = %f [pc]\n",*NNBMAX,*RS0);
	}
	
	for (j=0;j<N;j++) free (rarray[j]);
	free(rarray);
	for (j=0;j<N;j++) free (rarray2D[j]);
	free(rarray2D);

	return 0;
}

int cmd(double **star, int l, double Rgal, double *abvmag, double *vmag, double *BV, double *Teff, double *dvmag, double *dBV) {
	
	double lTeff, BC, kb;
	double bvc[8], bcc[8];
	double dbmag;
	double BCsun, abvmagsun;
	
	kb = 5.6704E-08*0.5*1.3914E9*0.5*1.3914E9/3.846E26;  //Stefan-Boltzmann constant in Lsun Rsun^-2 K^-4
	
	bvc[0] = -654597.405559323;
	bvc[1] = 1099118.61158915;
	bvc[2] = -789665.995692672;
	bvc[3] = 314714.220932623;
	bvc[4] = -75148.4728506455;
	bvc[5] = 10751.803394526;
	bvc[6] = -853.487897283685;
	bvc[7] = 28.9988730655392;
	
	bcc[0] = -4222907.80590972;
	bcc[1] = 7209333.13326442;
	bcc[2] = -5267167.04593882;
	bcc[3] = 2134724.55938336;
	bcc[4] = -518317.954642773;
	bcc[5] = 75392.2372207101;
	bcc[6] = -6082.7301194776;
	bcc[7] = 209.990478646363;
	
	BCsun = 0.11;  //sun's bolometric correction
	abvmagsun = 4.83; //sun's absolute V magnitude

	if (star[l][11] && (star[l][8]<14)) {
		*Teff = pow(star[l][12]/(4.0*PI*star[l][11]*star[l][11]*kb),0.25);
		if ((*Teff>3000.0) && (*Teff<55000.0)) {
					
			lTeff = log10(*Teff);
					
			*BV = bvc[0] + bvc[1]*lTeff + bvc[2]*pow(lTeff,2) + bvc[3]*pow(lTeff,3) + bvc[4]*pow(lTeff,4) + bvc[5]*pow(lTeff,5) + bvc[6]*pow(lTeff,6) + bvc[7]*pow(lTeff,7);
					
			BC = bcc[0] + bcc[1]*lTeff + bcc[2]*pow(lTeff,2) + bcc[3]*pow(lTeff,3) + bcc[4]*pow(lTeff,4) + bcc[5]*pow(lTeff,5) + bcc[6]*pow(lTeff,6) + bcc[7]*pow(lTeff,7);
					
			if (star[l][12]) *abvmag = -2.5*log10(star[l][12])-BC+BCsun+abvmagsun;
					
			*vmag = *abvmag + 5.0*log10(Rgal) - 5.0;
					
			
			double rand1, rand2, prand;
			
			do {
				rand1 = 2.0*drand48()-1.0;
				rand2 = 2.0*drand48()-1.0;
			} while (rand1*rand1+rand2*rand2 > 1.0);
			
			prand = sqrt(-2.0*log(rand1*rand1+rand2*rand2)/(rand1*rand1+rand2*rand2));
			*dvmag = rand1*prand*sqrt(pow(0.02,2) + pow(0.07*pow(10.0, 0.4*(*vmag-25.0)),2));
			dbmag = rand2*prand*sqrt(pow(0.02,2) + pow(0.07*pow(10.0, 0.4*(*vmag-25.0)),2));
			*dBV = *dvmag + dbmag;
			
			} else {
				*vmag = 9999.9;
				*abvmag = 9999.9;
				*BV = 9999.9;
				*dvmag = 0.0;
				*dBV = 0.0;
			}
		} else {
			*Teff = 0.0;
			*vmag = 9999.9;
			*abvmag = 9999.9;
			*BV = 9999.9;
			*dvmag = 0.0;
			*dBV = 0.0;
		}
	
	return 0;
}			

int output0(char *output, int N, int NNBMAX, double RS0, double dtadj, double dtout, double tcrit, double rvir, double mmean, int tf, int regupdate, int etaupdate, int mloss, int bin, int esc, double M, double mlow, double mup, double MMAX, double epoch, double dtplot, double Z, int nbin, double Q, double *RG, double *VG, double rtide, int gpu, double **star, int sse, int seed, double extmass, double extrad, double extdecay, double extstart){

	//Open output files
	char PARfile[50], NBODYfile[50], SSEfile[50];		
	FILE *PAR, *NBODY, *SSE12;
	sprintf(PARfile, "%s.input",output);
	PAR = fopen(PARfile,"w");
	sprintf(NBODYfile, "%s.fort.10",output);
	NBODY = fopen(NBODYfile,"w");

	int hrplot = 0;
	if (dtplot) hrplot = 1;
	if (sse) {
		sprintf(SSEfile, "%s.fort.12",output);
		SSE12 = fopen(SSEfile,"w");
		hrplot = 2;
	}		
	
	//write to .PAR file	
	fprintf(PAR,"1 5000000.0 0\n");
	fprintf(PAR,"%i 1 10 %i %i 1\n",N,seed,NNBMAX);
	fprintf(PAR,"0.02 0.02 %.8f %.8f %.8f %.8f 1.0E-03 %.8f %.8f\n",RS0,dtadj,dtout,tcrit,rvir,mmean);
	fprintf(PAR,"2 2 1 0 1 0 2 0 0 2\n");
	fprintf(PAR,"0 %i 0 %i 2 %i %i 0 %i 3\n",hrplot,tf,regupdate,etaupdate,mloss);
	fprintf(PAR,"0 %i %i 0 1 2 3 4 0 1\n",bin,esc);
	fprintf(PAR,"0 0 0 2 1 0 2 0 0 3\n");
	fprintf(PAR,"0 0 0 0 0 0 0 0 0 0\n");
	fprintf(PAR,"1.0E-5 1.0E-4 0.2 1.0 1.0E-06 0.001\n");
	fprintf(PAR,"2.350000 %.8f %.8f %i 0 %.8f %.8f %.8f\n",MMAX,mlow,nbin,Z,epoch,dtplot);
	fprintf(PAR,"%.2f 0.0 0.0 0.00000 0.125\n",Q);
	if (tf == 2) {
		fprintf(PAR,"%.8e %.8f\n",M1pointmass,sqrt(RG[0]*RG[0]+RG[1]*RG[1]+RG[2]*RG[2])/1000.0);
	} else if (tf == 3) {
		//old version:
		//fprintf(PAR,"%.6e %.6e %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",M1allen,M2allen,a2allen,b2allen,VCIRC,RCIRC,RG[0]/1000.0,RG[1]/1000.0,RG[2]/1000.0,VG[0],VG[1],VG[2]);

		//new version including bulge potential:
		fprintf(PAR,"%.6e %.6e %.6f %.6f %.6f %.6f %.6e %.6f %.6f\n",M1allen,M2allen,a2allen,b2allen,VCIRC,RCIRC, GMB, AR, GAM);
		fprintf(PAR,"%.6f %.6f %.6f %.6f %.6f %.6f\n", RG[0]/1000.0,RG[1]/1000.0,RG[2]/1000.0,VG[0],VG[1],VG[2]);
	}
	if (tf > 2)	fprintf(PAR,"%.6f %.6f %.6f %.6f\n",extmass,extrad,extdecay,extstart);

	
	
	
	//write to .fort.10 file
	int j;
	for (j=0;j<N;j++) {
		fprintf(NBODY,"%.16lf\t%.16lf %.16lf %.16lf\t%.16lf %.16lf %.16lf\n",star[j][0],star[j][1],star[j][2],star[j][3],star[j][4],star[j][5],star[j][6]);
	}
	
	//write to .fort.12 file
	if (sse) {
		for (j=0;j<N;j++) {
			fprintf(SSE12,"%.8lf\t%.0lf %.8lf %.8lf %.8lf\n",star[j][0]*M,star[j][8],star[j][7],star[j][9],star[j][10]);
			//,star[j][13],star[j][14]);
		}
	}
	
	fclose(PAR);
	fclose(NBODY);
	if (bin == 5) {
		fclose(SSE12);
		printf("\nData written to %s, %s and %s\n", PARfile, NBODYfile, SSEfile);
	} else {
		printf("\nData written to %s and %s\n", PARfile, NBODYfile);
	}
	
	return 0;
	
}

int output1(char *output, int N, double dtadj, double dtout, double tcrit, double rvir, double mmean, int tf, int regupdate, int etaupdate, int mloss, int bin, int esc, double M, double mlow, double mup, double MMAX, double epoch, double Z, int nbin, double Q, double *RG, double *VG, double rtide, int gpu, double **star){

	//Open output files
	char PARfile[50], NBODYfile[50];		
	FILE *PAR, *NBODY;
	sprintf(PARfile, "%s.PAR",output);
	PAR = fopen(PARfile,"w");
	sprintf(NBODYfile, "%s.NBODY",output);
	NBODY = fopen(NBODYfile,"w");
	
	//write to .PAR file	
	fprintf(PAR,"1 5000000.0 0\n");
	fprintf (PAR,"%i 1 10 3 8\n",N);
	fprintf(PAR,"0.02 %.8f %.8f %.8f 100.0 1000.0 1.0E-02 %.8f %.8f\n",dtadj,dtout,tcrit,rvir,mmean);
	fprintf(PAR,"2 2 1 0 1 0 2 0 0 2\n");
	fprintf(PAR,"0 0 0 %i 2 %i %i 0 %i 2\n",tf,regupdate,etaupdate,mloss);
	fprintf(PAR,"0 %i %i 0 1 2 0 0 0 2\n",bin, esc);
	fprintf(PAR,"0 0 0 0 1 0 0 2 0 3\n");
	fprintf(PAR,"0 0 0 0 0 0 0 0 0 0\n");
	fprintf(PAR,"%.8f %.8f %.8f\n",M,mlow,mup);
	fprintf(PAR,"1.0E-5 1.0E-4 0.01 1.0 1.0E-06 0.01\n");
	fprintf(PAR,"2.350000 %.8f %.8f %i %.8f %.8f 100000.0\n",MMAX,mlow, nbin, Z, epoch);	
	fprintf(PAR,"%.2f 0.0 0.0 0.00000\n",Q);
	if (tf == 1) {
		rtide = pow(1.0*M/(3.0*M1pointmass),1.0/3.0)*sqrt(RG[0]*RG[0]+RG[1]*RG[1]+RG[2]*RG[2]);
		fprintf(PAR,"%i %.8f\n",0,sqrt(1.0/(3.0*pow(rtide,3))));
		
	} else if (tf == 2) {
		fprintf(PAR,"%.8e %.8f\n",M1pointmass,sqrt(RG[0]*RG[0]+RG[1]*RG[1]+RG[2]*RG[2]));
	} else if (tf == 3) {
		fprintf(PAR,"%.6e %.6f %.6e %.6f %.6f %.6e %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",M1allen,b1allen,M2allen,a2allen,b2allen,M3allen,a3allen,RG[0]/1000.0,RG[1]/1000.0,RG[2]/1000.0,VG[0],VG[1],VG[2]);
	}
	if (gpu) fprintf(PAR,"1.0\n");
	
	//write to .NBODY file
	int j;
	for (j=0;j<N;j++) {
		fprintf(NBODY,"%.16lf\t%.16lf %.16lf %.16lf\t%.16lf %.16lf %.16lf\n",star[j][0],star[j][1],star[j][2],star[j][3],star[j][4],star[j][5],star[j][6]);
	}
	
	fclose(PAR);
	fclose(NBODY);
	printf("\nData written to %s and %s\n", PARfile, NBODYfile);
	
	return 0;
	
}

int output2(char *output, int N, int NNBMAX, double RS0, double dtadj, double dtout, double tcrit, double rvir, double mmean, int tf, int regupdate, int etaupdate, int mloss, int bin, int esc, double M, double mlow, double mup, double MMAX, double epoch, double dtplot, double Z, int nbin, double Q, double *RG, double *VG, double rtide, int gpu, double **star, int sse, int seed, double extmass, double extrad, double extdecay, double extstart){

    double GMG_temp, DISK_temp, A_temp, B_temp, VCIRC_temp, RCIRC_temp, GMB_temp, AR_temp, GAM_temp, ZDUM1_temp, ZDUM2_temp, ZDUM3_temp, ZDUM4_temp;
    
	//Open output files
	char PARfile[50], NBODYfile[50], SSEfile[50];		
	FILE *PAR, *NBODY, *SSE12;
	sprintf(PARfile, "%s.PAR",output);
	PAR = fopen(PARfile,"w");
	sprintf(NBODYfile, "%s.NBODY",output);
	NBODY = fopen(NBODYfile,"w");
	
	int hrplot = 0;
	if (dtplot) hrplot = 1;
	if (sse) {
		sprintf(SSEfile, "%s.fort.12",output);
		SSE12 = fopen(SSEfile,"w");
		hrplot = 2;
	}		
	
	
	//write to .PAR file	
	fprintf(PAR,"1 5000000.0 0\n");
	fprintf(PAR,"%i 1 10 %i %i 1\n",N,seed,NNBMAX);
	fprintf(PAR,"0.02 0.02 %.8f %.8f %.8f %.8f 1.0E-03 %.8f %.8f\n",RS0,dtadj,dtout,tcrit,rvir,mmean);
	fprintf(PAR,"2 2 1 0 1 0 2 0 0 2\n");
	fprintf(PAR,"0 %i 0 %i 2 %i %i 0 %i 3\n",hrplot,tf,regupdate,etaupdate,mloss);
	fprintf(PAR,"0 %i %i 0 1 2 3 4 0 1\n",bin, esc);
	fprintf(PAR,"0 0 0 2 1 0 2 0 0 3\n");
	fprintf(PAR,"0 0 0 0 0 0 0 0 0 0\n");
	fprintf(PAR,"1.0E-5 1.0E-4 0.2 1.0 1.0E-06 0.001\n");
	fprintf(PAR,"2.350000 %.8f %.8f %i 0 %.8f %.8f %.8f\n",MMAX,mlow,nbin,Z,epoch,dtplot);
	fprintf(PAR,"%.2f 0.0 0.0 0.00000 0.125\n",Q);

	if (tf == 2) {
		fprintf(PAR,"%.8e %.8f\n",M1pointmass,sqrt(RG[0]*RG[0]+RG[1]*RG[1]+RG[2]*RG[2])/1000.0);
	} else if (tf == 3) {
		//new version including Hernquist bulge potential + NFW:
        GMG_temp = 0.0;        //mass of central point-mass
        DISK_temp = M2_LMJ;    //mass of Myamoto disk
        A_temp = a2_LMJ;       //scale length of Myamoto disk
        B_temp = b2_LMJ;       //scale height of Myamoto disk
        VCIRC_temp = 0.0;      //desired circular velocity at RCIRC_t for logarighmic potential
        RCIRC_temp = 0.0;        //radius at which circular velocity shall be VCIRC_t
        GMB_temp = M1_LMJ;     //bulge mass
        AR_temp = b1_LMJ;      //bulge scale radius
        GAM_temp = 2.0;        //bulge profile slope (2 = Hernquist bulge)
        ZDUM1_temp = 0.0;      //smoothing length for central point mass
        ZDUM2_temp = M_NFW;    //characteristic mass of NFW halo (~M(<5.3*ZDUM3_t))
        ZDUM3_temp = R_NFW;    //scale radius of NFW profile
        ZDUM4_temp = q_halo;   //flattening of NFW or logarithmic halo potential along z-axis
		fprintf(PAR,"%.8e %.8e %.6f %.6f %.6f %.6f %.8e %.6f %.6f %.6f %.8e %.8e %.6f\n", GMG_temp, DISK_temp, A_temp, B_temp, VCIRC_temp, RCIRC_temp, GMB_temp, AR_temp, GAM_temp, ZDUM1_temp, ZDUM2_temp, ZDUM3_temp, ZDUM4_temp);
		fprintf(PAR,"%.6f %.6f %.6f %.6f %.6f %.6f\n", RG[0]/1000.0,RG[1]/1000.0,RG[2]/1000.0,VG[0],VG[1],VG[2]);
	}

	if (tf > 2) fprintf(PAR,"%.2f %.2f %.2f %.2f\n",extmass,extrad,extdecay,extstart);
	
	
	
	//write to .NBODY file
	int j;
	for (j=0;j<N;j++) {
		fprintf(NBODY,"%.16lf\t%.16lf %.16lf %.16lf\t%.16lf %.16lf %.16lf\n",star[j][0],star[j][1],star[j][2],star[j][3],star[j][4],star[j][5],star[j][6]);
	}

	//write to .fort.12 file
	if (sse) {
		for (j=0;j<N;j++) {
			fprintf(SSE12,"%.8lf\t%.0lf %.8lf %.8lf %.8lf\n",star[j][0]*M,star[j][8],star[j][7],star[j][9],star[j][10]);
					//,star[j][13],star[j][14]);
		}
	}
	
	fclose(PAR);
	fclose(NBODY);
	if (bin == 5) {
		fclose(SSE12);
		printf("\nData written to %s, %s and %s\n", PARfile, NBODYfile, SSEfile);
	} else {
		printf("\nData written to %s and %s\n", PARfile, NBODYfile);
	}
	return 0;
	
}

int output3(char *output, int N, double rvir, double rh, double mmean, double M, double epoch, double Z, double *RG, double *VG, double rtide, double **star, double Rgal){
	//Open output files
	char tablefile[20];		
	FILE *TABLE;
	sprintf(tablefile, "%s.txt",output);
	TABLE = fopen(tablefile,"w");
	
	
	//write to .txt file
	int j;
	//fprintf(TABLE,"# Star cluster generated by McLuster with the following parameters:\n");
	//fprintf(TABLE,"#\n");
	//fprintf(TABLE,"# Mass = %.3lf\t(%i stars with mean mass of %.3f)\n",M,N,mmean);
	//fprintf(TABLE,"# Half-mass radius, virial radius, tidal radius = %.3lf pc, %.3lf pc, %.3lf pc\n",rh, rvir, rtide);
	//fprintf(TABLE,"# Galactic orbit (RGx,RGy,RGz),(VGx,VGy,VGz) = (%.3lf,%.3lf,%.3lf) pc, (%.3lf,%.3lf,%.3lf) km/s\n",RG[0],RG[1],RG[2],VG[0],VG[1],VG[2]);
	//fprintf(TABLE,"# Age = %.3lf Myr\n",epoch);
	//fprintf(TABLE,"# Metallicity = %.3lf\n",Z);
	//fprintf(TABLE,"#\n");
	
#ifdef SSE	
	double abvmag, vmag, BV, Teff, dvmag, dBV;
	fprintf(TABLE,"#Mass_[Msun]\tx_[pc]\t\t\ty_[pc]\t\t\tz_[pc]\t\t\tvx_[km/s]\t\tvy_[km/s]\t\tvz_[km/s]\t\tMass(t=0)_[Msun]\tkw\tepoch1_[Myr]\tspin\t\tR_[Rsun]\tL_[Lsun]\tepoch_[Myr]\tZ\t\tM_V_[mag]\tV_[mag]\t\tB-V_[mag]\tT_eff_[K]\tdV_[mag]\td(B-V)\n");
	for (j=0;j<N;j++) {
		cmd(star, j, Rgal, &abvmag, &vmag, &BV, &Teff, &dvmag, &dBV);
		fprintf(TABLE,"%8lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%8lf\t\t%.0lf\t%8lf\t%8lf\t%8lf\t%8lf\t%8lf\t%8lf\t%8lf\t%8lf\t%8lf\t%8lf\t%8lf\t%8lf\n",star[j][0],star[j][1],star[j][2],star[j][3],star[j][4],star[j][5],star[j][6],star[j][7],star[j][8],star[j][9],star[j][10],star[j][11],star[j][12],star[j][13],star[j][14], abvmag, vmag, BV, Teff, dvmag, dBV);
	}
#else	
	fprintf(TABLE,"#Mass_[Msun]\tx_[pc]\t\t\ty_[pc]\t\t\tz_[pc]\t\t\tvx_[km/s]\t\tvy_[km/s]\t\tvz_[km/s]\n");
	for (j=0;j<N;j++) {
		fprintf(TABLE,"%8lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\t%.16lf\n",star[j][0],star[j][1],star[j][2],star[j][3],star[j][4],star[j][5],star[j][6]);
	}	
#endif
	fclose(TABLE);
	printf("\nData written to %s\n", tablefile);
	
	return 0;
	
}

int output4(char *output, int N, int NNBMAX, double RS0, double dtadj, double dtout, double tcrit, double rvir, double mmean, int tf, int regupdate, int etaupdate, int mloss, int bin, int esc, double M, double mlow, double mup, double MMAX, double epoch, double dtplot, double Z, int nbin, double Q, double *RG, double *VG, double rtide, int gpu, double **star, int sse, int seed, double extmass, double extrad, double extdecay, double extstart){
	
	//Open output files
	char PARfile[50], NBODYfile[50], SSEfile[50];		
	FILE *PAR, *NBODY, *SSE12;
	sprintf(PARfile, "%s.input",output);
	PAR = fopen(PARfile,"w");
	sprintf(NBODYfile, "%s.fort.10",output);
	NBODY = fopen(NBODYfile,"w");
	
	int hrplot = 0;
	if (dtplot) hrplot = 1;
	if (sse) {
		sprintf(SSEfile, "%s.fort.12",output);
		SSE12 = fopen(SSEfile,"w");
		hrplot = 2;
	}		
	
	//write to .PAR file	
	fprintf(PAR,"1 5000000.0 0\n");
	fprintf(PAR,"%i 1 10 %i %i 1\n",N,seed,NNBMAX);
	fprintf(PAR,"0.02 0.02 %.8f %.8f %.8f %.8f 1.0E-03 %.8f %.8f\n",RS0,dtadj,dtout,tcrit,rvir,mmean);
	fprintf(PAR,"2 2 1 0 1 0 2 0 0 2\n");
	fprintf(PAR,"-1 %i 0 %i 2 %i %i 0 %i 3\n",hrplot,tf,regupdate,etaupdate,mloss);
	fprintf(PAR,"0 %i %i 0 1 2 3 4 0 -1\n",bin,esc);
	fprintf(PAR,"0 0 0 2 1 0 2 0 0 3\n");
	fprintf(PAR,"0 1 0 1 0 0 0 0 0 0\n");
	fprintf(PAR,"1.0E-5 1.0E-4 0.2 1.0 1.0E-06 0.001\n");
	fprintf(PAR,"2.350000 %.8f %.8f %i 0 %.8f %.8f %.8f\n",MMAX,mlow,nbin,Z,epoch,dtplot);
	fprintf(PAR,"%.2f 0.0 0.0 0.00000 0.125\n",Q);
	if (tf == 2) {
		fprintf(PAR,"%.8e %.8f\n",M1pointmass,sqrt(RG[0]*RG[0]+RG[1]*RG[1]+RG[2]*RG[2])/1000.0);
	} else if (tf == 3) {
		//old version:
		//fprintf(PAR,"%.6e %.6e %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",M1allen,M2allen,a2allen,b2allen,VCIRC,RCIRC,RG[0]/1000.0,RG[1]/1000.0,RG[2]/1000.0,VG[0],VG[1],VG[2]);
		
		//new version including bulge potential:
		fprintf(PAR,"%.6e %.6e %.6f %.6f %.6f %.6f %.6e %.6f %.6f\n",M1allen,M2allen,a2allen,b2allen,VCIRC,RCIRC, GMB, AR, GAM);
		fprintf(PAR,"%.6f %.6f %.6f %.6f %.6f %.6f\n", RG[0]/1000.0,RG[1]/1000.0,RG[2]/1000.0,VG[0],VG[1],VG[2]);
	}
	if (tf > 2)	fprintf(PAR,"%.6f %.6f %.6f %.6f\n",extmass,extrad,extdecay,extstart);
	fprintf(PAR,"10000.0 2 0\n");
	fprintf(PAR,"10000.0 2 0\n");
	
	
	
	//write to .fort.10 file
	int j;
	for (j=0;j<N;j++) {
		fprintf(NBODY,"%.16lf\t%.16lf %.16lf %.16lf\t%.16lf %.16lf %.16lf\n",star[j][0],star[j][1],star[j][2],star[j][3],star[j][4],star[j][5],star[j][6]);
	}
	
	//write to .fort.12 file
	if (sse) {
		for (j=0;j<N;j++) {
			fprintf(SSE12,"%.8lf\t%.0lf %.8lf %.8lf %.8lf\n",star[j][0]*M,star[j][8],star[j][7],star[j][9],star[j][10]);
			//,star[j][13],star[j][14]);
		}
	}
	
	fclose(PAR);
	fclose(NBODY);
	if (bin == 5) {
		fclose(SSE12);
		printf("\nData written to %s, %s and %s\n", PARfile, NBODYfile, SSEfile);
	} else {
		printf("\nData written to %s and %s\n", PARfile, NBODYfile);
	}
	
	return 0;
	
}

void info(char *output, int N, double Mcl, int profile, double W0, double S, double D, double Q, double Rh, double gamma[], double a, double Rmax, double tcrit, int tf, double RG[], double VG[], int mfunc, double single_mass, double mlow, double mup, double alpha[], double mlim[], double alpha_L3, double beta_L3, double mu_L3, int weidner, int mloss, int remnant, double epoch, double Z, int prantzos, int nbin, double fbin, int pairing, double msort, int adis, double amin, double amax, int eigen, int BSE, double extmass, double extrad, double extdecay, double extstart, int code, int seed, double dtadj, double dtout, double dtplot, int gpu, int regupdate, int etaupdate, int esc, int units, int match, int symmetry, int OBperiods) {
	int i;
	char INFOfile[50];		
	FILE *INFO;
	sprintf(INFOfile, "%s.info",output);
	INFO = fopen(INFOfile,"w");

	time_t timep; 
	time(&timep);
	fprintf(INFO, "\nMcLuster model generated: %s\n\n",ctime(&timep));	
	
	fprintf(INFO, "N = %i\n",N);
	fprintf(INFO, "Mcl = %g\n",Mcl);
	fprintf(INFO, "profile = %i\n", profile);	
	fprintf(INFO, "W0 = %g\n",W0);
	fprintf(INFO, "S = %g\n",S);
	fprintf(INFO, "D = %g\n",D);
	fprintf(INFO, "Q = %g\n",Q);
	fprintf(INFO, "Rh = %g\n",Rh);
	fprintf(INFO, "gammas = %g\t%g\t%g\n",gamma[0],gamma[1],gamma[2]);
	fprintf(INFO, "a = %g\n",a);
	fprintf(INFO, "Rmax = %g\n",Rmax);
	fprintf(INFO, "tcrit = %g\n",tcrit);
	fprintf(INFO, "tf = %i\n",tf);
	fprintf(INFO, "RG = %g\t%g\t%g\n",RG[0],RG[1],RG[2]);
	fprintf(INFO, "VG = %g\t%g\t%g\n",VG[0],VG[1],VG[2]);
	fprintf(INFO, "mfunc = %i\n",mfunc);
	fprintf(INFO, "single_mass = %g\n",single_mass);
	fprintf(INFO, "mlow = %g\n",mlow);
	fprintf(INFO, "mup = %g\n",mup);
	fprintf(INFO, "alpha = "); 
	for(i=0;i<MAX_AN;i++) fprintf(INFO, "%g\t",alpha[i]); 
	fprintf(INFO, "\n");
	fprintf(INFO, "mlim = "); 
	for(i=0;i<MAX_MN;i++) fprintf(INFO, "%g\t",mlim[i]); 
	fprintf(INFO, "\n");
	fprintf(INFO, "alpha_L3 = %g\n",alpha_L3);
	fprintf(INFO, "beta_L3 = %g\n",beta_L3);
	fprintf(INFO, "mu_L3 = %g\n",mu_L3);
	fprintf(INFO, "weidner = %i\n",weidner);
	fprintf(INFO, "mloss = %i\n",mloss);
	fprintf(INFO, "remnant = %i\n",remnant);
	fprintf(INFO, "epoch = %g\n",epoch);
	fprintf(INFO, "Z = %g\n",Z);
	fprintf(INFO, "prantzos = %i\n",prantzos);
	fprintf(INFO, "nbin = %i\n",nbin);
	fprintf(INFO, "fbin = %g\n",fbin);
	fprintf(INFO, "pairing = %i\n",pairing);
	fprintf(INFO, "msort = %g\n",msort);
	fprintf(INFO, "adis = %i\n",adis);
	fprintf(INFO, "OBperiods = %i\n",OBperiods);
	fprintf(INFO, "amin = %g\n",amin);
	fprintf(INFO, "amax = %g\n",amax);
	fprintf(INFO, "eigen = %i\n",eigen);
	fprintf(INFO, "BSE = %i\n",BSE);
	fprintf(INFO, "extmass = %g\n",extmass);
	fprintf(INFO, "extrad = %g\n",extrad);
	fprintf(INFO, "extdecay = %g\n",extdecay);
	fprintf(INFO, "extstart = %g\n",extstart);
	fprintf(INFO, "code = %i\n",code);
	fprintf(INFO, "seed = %i\n",seed);
	fprintf(INFO, "dtadj = %g\n",dtadj);
	fprintf(INFO, "dtout = %g\n",dtout);
	fprintf(INFO, "dtplot = %g\n",dtplot);
	fprintf(INFO, "gpu = %i\n",gpu);
	fprintf(INFO, "regupdate = %i\n",regupdate);
	fprintf(INFO, "etaupdate = %i\n",etaupdate);
	fprintf(INFO, "esc = %i\n",esc);
	fprintf(INFO, "units = %i\n",units);
	fprintf(INFO, "match = %i\n",match);
	fprintf(INFO, "symmetry = %i\n",symmetry);

	fclose(INFO);

}

void help(double msort) {
	printf("\n Usage: mcluster -[N|M|P|W|R|r|c|g|S|D|T|Q|C|A|O|G|o|f|a|m|B|b|p|s|t|e|Z|X|V|x|u|h|?]\n");
	printf("                                                                     \n");
	printf("       -N <number> (number of stars)                                 \n");
	printf("       -M <value> (mass of cluster; specify either N or M)           \n");
	printf("       -P <0|1|2|3|-1> (density profile; 0= Plummer, 1= King (1966), \n");
	printf("                   2= Subr et al. (2007) mass-segregated,            \n");
	printf("                   3= 2-dimensional EFF/Nuker template,              \n");
	printf("                   -1= no density gradient)                          \n");
	printf("       -W <1-12> (W0 parameter for King model)                       \n");
	printf("       -R <value> (half-mass radius [pc], ignored for P = 3)         \n");
	printf("       -r <value> (scale radius of EFF/Nuker template [pc])          \n");
	printf("       -c <value> (cut-off radius of EFF/Nuker template [pc])        \n");
	printf("       -g <value> (power-law slope(s) of EFF/Nuker template; use     \n");
	printf("                   once for EFF template; use three times for Nuker  \n");
	printf("                   template (outer slope, inner slope, transition)   \n");
	printf("       -S <0.0-1.0> (degree of mass segregation; 0.0= no segregation)\n");
	printf("       -D <1.6-3.0> (fractal dimension; 3.0= no fractality)          \n");
	printf("       -T <value> (tcrit in N-body units)                            \n");
	printf("       -Q <value> (virial ratio)                                     \n");
	printf("       -C <0|1|3> (code; 0= Nbody6, 1= Nbody4, 3= table of stars)    \n");
	printf("       -A <value> (dtadj in N-body units)                            \n");
	printf("       -O <value> (deltat in N-body units)                           \n");
	printf("       -G <0|1> (GPU usage; 0= no GPU, 1= use GPU)                   \n");
	printf("       -o <name> (output name of cluster model)                      \n");
	printf("       -f <0|1|2|3|4> (IMF; 0= no IMF, 1= Kroupa (2001),             \n");
	printf("             2= user defined, 3= Kroupa (2001) with optimal sampling,\n");
	printf("             4= L3 IMF (Maschberger 2012))                           \n");
	printf("       -a <value> (IMF slope; for user defined IMF, may be used      \n"); 
	printf("                   multiple times, from low mass to high mass;       \n");
	printf("                   for L3 IMF use three times for alpha, beta and mu)\n");
	printf("       -m <value> (IMF mass limits, has to be used multiple times    \n");
	printf("                 (at least twice), from low mass to high mass [Msun])\n");
	printf("       -B <number> (number of binary systems)                        \n");
	printf("       -b <value> (binary fraction, specify either B or b)           \n");
	printf("       -p <0|1|2> (binary pairing, 0= random, 1= ordered for M>%.1f Msun,\n",msort);
	printf("                   2= random but separate pairing for M>%.1f Msun)\n",msort);
	printf("       -s <number> (seed for randomization; 0= randomize by timer)   \n");
	printf("       -t <0|1|2|3> (tidal field; 0= no tidal field, 1= near-field,  \n");
	printf("                    2= point-mass, 3= Milky-Way potential)           \n");
	printf("       -e <value> (epoch for stellar evolution [Myr])                \n");
	printf("       -Z <value> (metallicity [0.0001-0.03, 0.02 = solar])          \n");
	printf("       -X <value> (galactocentric radius vector, use 3x, [pc])       \n");
	printf("       -V <value> (cluster velocity vector, use 3x, [km/s])          \n");
	printf("       -x <value> (specify external (gas) Plummer potential, use 4x, \n");
	printf("                  1. gas mass [Msun], 2. Plummer radius [pc]         \n");
	printf("                  3. decay time for gas expulsion [Myr], 4. delay    \n");
	printf("                  time for start of gas expulsion [Myr])             \n");
	printf("       -u <0|1> (output units; 0= Nbody, 1= astrophysical)           \n");
	printf("       -h (display this help)                                        \n");
	printf("       -? (display this help)                                        \n");
	printf("                                                                     \n");
	printf(" Examples: mcluster -N 1000 -R 0.8 -P 1 -W 3.0 -f 1 -B 100 -o test1  \n");
	printf("           mcluster -f 2 -m 0.08 -a -1.35 -m 0.5 -a -2.7 -m 100.0    \n");
	printf("           mcluster -t 3 -X 8500 -X 0 -X 0 -V 0 -V 220 -V 0          \n");
	printf("           mcluster -D 1.6 -Q 0.4 -P -1                              \n");
	printf("                                                                     \n");

}	


