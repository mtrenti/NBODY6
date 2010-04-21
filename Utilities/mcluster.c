/***************************************************************************
 *   Copyright (C) 2009 by Andreas Kuepper                                 *
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
 ***************************************************************************/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include<sys/stat.h>
#include <getopt.h>

//Constants
#define G 0.0043009211 //in pc*km2*s-2*Msun
#define PI   4.0*atan(1.0)   /* PI */
#define TWOPI   8.0*atan(1.0)   /* 2PI */
#define GNBODY 1.0
#define MNBODY 1.0
#define PARSEC  3.08568E13      /* KM pro PC */

//Functions
#define max(a,b)         (a < b) ?  (b) : (a)
#define Lifetime(Mstar)  1.13E4*pow(Mstar,-3)+0.6E2*pow(Mstar,-0.75)+1.2 //Myr	[Prantzos 2007]

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

//Point-mass potential - constants:
#define M1pointmass 9.565439E+10    //solar masses

//Mass function variables
#define MAX_AN  5
#define MAX_MN  6

int generate_m1(int *N, double **star, double mlow, double mup, double *M, double *mmean, double MMAX, double Mcl);
int generate_m2(int an, double *mlim, double *alpha, double Mcl, double M_tmp, double *subcount, int *N, double *mmean, double *M, double **star, double MMAX);
double subint(double min, double max, double alpha);
double mlow(double mhigh, double alpha, double norma, double delta);
int generate_plummer(int N, double **star, double rtide, double rvir);
int generate_king(int N, double W0, double **star, double *rvir, double *rh, double *rking); 
double densty(double z);
int odeint(double ystart0, double ystart1, double x1, double x2, double den, int *kount, double *xp, double **yp, int M, int KMAX);	 
int derivs(double x, double *y, double *dydx, double den);
int rkqc(double *y,double *dydx, double *x, double *h,double den, double *yscal, double *hdid, double *hnext, double TOL);	
int rk4(double x, double *y, double *dydx, double h, double *yout, double den);
int cmpmy(float *x1, float *x2);
int cmpmy_reverse(float *x1, float *x2);
int get_binaries(int nbin, double **star, double M, double rvir, int pairing, int N, int adis, double amin, double amax);
int order(double **star, int N, double M);
double rtnewt (double ecc, double ma);
int output0(char *output, int N, int NNBMAX, double RS0, double dtadj, double dtout, double tcrit, double rvir, double mmean, int tf, int regupdate, int etaupdate, int mloss, int bin, int esc, double M, double mlow, double mup, double MMAX, double epoch, double dtplot, double Z, int nbin, double Q, double *RG, double *VG, double rtide, int gpu, double **star);
int output1(char *output, int N, double dtadj, double dtout, double tcrit, double rvir, double mmean, int tf, int regupdate, int etaupdate, int mloss, int bin, int esc, double M, double mlow, double mup, double MMAX, double epoch, double Z, int nbin, double Q, double *RG, double *VG, double rtide, int gpu, double **star);
int output2(char *output, int N, int NNBMAX, double RS0, double dtadj, double dtout, double tcrit, double rvir, double mmean, int tf, int regupdate, int etaupdate, int mloss, int bin, int esc, double M, double mlow, double mup, double MMAX, double epoch, double dtplot, double Z, int nbin, double Q, double *RG, double *VG, double rtide, int gpu, double **star);
int help();

/************************************************************
 * to do:
 * - fractalized initial conditions (Goodwin et al 2009)
 * - mass segregation
 * - Eigenevolution (Kroupa 1995)
 ************************************************************/


int main (int argv, char **argc) {

	/*******************
	 * Input variables *
	 *******************/
	
	//Basic physical parameters
	int N = 1024;		     	    //number of stars, Mcl will be set to 0 if specified!
	double Mcl = 0.0;               //total mass of the cluster, only used when N is set to 0, necessary for usage of maximum stellar mass relation of Weidner & Kroupa 2007
	int profile = 0;				//density profile; =0 Plummer sphere, =1 King profile
	double W0 = 1.0;				//King's W0 paramter [0.3-12.0]
	double Rh = 0.8;				//Half-mass radius [pc]
	double tcrit = 500.0;			//Simulation time [N-body units (Myr in Nbody6 custom)]
	int tf = 3;						//tidal field: =1 Near-field approximation, =2 point-mass galaxy, =3 Allen & Santillan (1991) MW potential (or Sverre's version of it)
	double RG[3] = {8500.0,0.0,0.0}; //Initial Galactic coordinates of the cluster [pc]
	double VG[3] = {0.0,220.0,0.0};  //Initial velocity of the cluster [km/s]
	double Q = 0.5;					//Initial virial ratio; =0.5 virial equilibrium, <0.5 collapsing, >0.5 expanding
	
	//Mass function parameters
	int mfunc = 0;					//0 = single mass stars; 1 = use Kroupa (2001) mass function; 2 = use multi power law
	double single_mass = 1.0;		//stellar mass in case of single-mass cluster 
	double alpha[MAX_AN] = {-1.35, -2.35, -2.7, 0.0, 0.0};		//alpha slopes for mfunc = 2
	double mlim[MAX_MN] = {0.08, 0.5, 4.0, 150.0, 0.0, 0.0};	//mass limits for mfunc = 2
	double mlow = 0.08;				//lower mass limit for mfunc = 1
	double mup = 100.0;				//upper mass limit for mfunc = 1
	int weidner = 0;				//Usage of Weidner & Kroupa 2007 relation for most massive star; =0 off, =1 on
	int mloss = 0;					//Stellar evolution; 0 = off, 3 = Eggleton, Tout & Hurley [KZ19]
	double epoch = 0.0;				//Star burst has been ... Myr before [e.g. 1000.0, default = 0.0]
	double Z = 0.02;				//Metallicity [0.0001-0.03, 0.02 = solar]
	double FeH = -1.41;				//Metallicity [Fe/H], only used when Z is set to 0
	int prantzos = 0;				//Usage of Prantzos 2007 relation for the life-times of stars. Set upper mass limit to Lifetime(mup) >= epoch
	
	//Binary parameters
	int nbin = 0;				    //Number of primordial binary systems
	double fbin = 0.0;				//Primordial binary fraction number of binary systems = 0.5*N*fbin, only used when nbin is set to 0 
	int pairing = 0;				//Pairing of binary component masses; 0= random pairing, 1= ordered pairing for M>5Msun
	int adis = 1;					//Semi-major axis distribution; 0= flat ranging from amin to amax, 1= based on Kroupa (1995) period distribution
	double amin = 0.0001;			//Minimum semi-major axis for adis = 0 [pc]
	double amax = 0.01;				//Maximum semi-major axis for adis = 0 [pc]
	
	//Code parameters
	int code = 2;					//Nbody version: =0 Nbody6, =1 Nbody4, =2 Nbody6 custom
	unsigned int seed = 0;			//number seed for random number generator; =0 for randomization by local time
	char *output = "test";			//name of output files
	double dtadj = 1.0;	     		//DTADJ [N-body units (Myr in Nbody6 custom)], energy-check time step
	double dtout = 1.0;		    	//DELTAT [N-body units (Myr in Nbody6 custom)], output interval, must be multiple of DTADJ
	double dtplot = 500.0;			//DTPLOT [N-body units (Myr in Nbody6 custom)], output of HRdiagnostics, should be multiple of DTOUT
	int gpu = 1;					//Use of GPU, 0= off, 1= on
	int regupdate = 1;				//Update of regularization parameters during computation; 0 = off, 0 > on
	int etaupdate = 1;				//Update of ETAI & ETAR during computation; 0 = off, 0 > on
	int esc = 1;					//Removal of escapers; 0 = no removal, 1 = regular removal at 2*R_tide
	int units = 1;				    //units of mcluster output; 0= Nbody-Units, 1= astrophysical units
	
	//McLuster internal parameters
	int check = 0;					//Make energy check at end of mcluster; =0 off, =1 on
	double Zsun = 0.02;				//Solar metallicity
	int NMAX = 500000;				//Maximum number of stars allowed in mcluster
	double upper_IMF_limit = 150.0; //Maximum stellar mass allowed in mcluster [Msun]
	int an = 0;						//counter for number of alpha slopes for mfunc = 2
	int mn = 0;						//counter for number of mass limits for mfunc = 2
	int xn = 0;						//counter for components of galactocentric radius vector
	int vn = 0;						//counter for components of cluster velocity vector 
	
	//Command line input
	int option;
	while ((option = getopt(argv, argc, "N:M:P:W:R:T:Q:C:A:O:G:o:f:a:m:B:b:S:t:X:V:h:?")) != -1) switch (option)
	{
		case 'N': N = atoi(optarg); Mcl = 0.0; break;
		case 'M': Mcl = atof(optarg); break;
		case 'P': profile = atoi(optarg); break;
		case 'W': W0 = atof(optarg); break;
		case 'R': Rh = atof(optarg); break;
		case 'T': tcrit = atof(optarg); break;
		case 'Q': Q = atof(optarg); break;
		case 'C': code = atoi(optarg); break;
		case 'A': dtadj = atof(optarg); break;
		case 'O': dtout = atof(optarg); break;
		case 'G': gpu = atoi(optarg); break;
		case 'o': output = optarg; break;
		case 'f': mfunc = atoi(optarg); break;
		case 'a' :
			if (an < MAX_AN) { alpha[an] = atof(optarg); an++; break; }
			else { printf("\nError: Number of alphas exceeded maximum limit of %d\n", MAX_AN); return 1; }
		case 'm' :
			if (mn < MAX_MN) { mlim[mn] = atof(optarg); mn++; break;
			} else { printf("\nError: Number of mass params exceded maximum limit of %d\n", MAX_MN); return 1; }
		case 'B': nbin = atoi(optarg); break;
		case 'b': fbin = atof(optarg); break;
		case 'S': seed = atoi(optarg); break;
		case 't': tf = atoi(optarg); break;
		case 'X' :
			if (xn < 3) { RG[xn] = atof(optarg); xn++; break; }
		case 'V' :
			if (vn < 3) { VG[vn] = atof(optarg); vn++; break; }
		case ':':
		case 'h':	help(); return 1;
		case '?':	help(); return 1;
	};
	
	
	/*********
	 * Start *
	 *********/
	
	printf("\n-----START----         \n"); 

	
	clock_t t1, t2;
	t1 = clock();							//start stop-watch
  
	if (seed) srand48(seed);				//initialize random number generator by seed
	else srand48((unsigned) time(NULL));	//initialize random number generator by local time
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
	double submass[MAX_AN], subcount[MAX_AN], norm[MAX_AN], N_tmp, M_tmp; //mass function parameters for mfunc = 2	
	
	if ((Mcl) && (N)) {
		printf("\nCAUTION: set either Mcl or N to 0!\nSet Mcl to 0 by default!\n");
		Mcl = 0.0;
	}	
    
	
	/***********************
	 * Generate star array *
	 ***********************/

	int columns = 7;
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
	
	if (!N && weidner && mfunc) {
		mup = upper_IMF_limit;  
		printf("\nUsing Weidner & Kroupa (2007) relation for upper stellar mass limit\n");
		
		if (Mcl < 1000.0) {
			MMAX = (log10(Mcl)*0.540563175) - 0.14120167;
			MMAX = pow(10.0,MMAX);
		} else if (Mcl < 3300.0) {
			MMAX = (log10(Mcl)*0.19186051) + 0.9058611;
			MMAX = pow(10.0,MMAX);
		} else {
			MMAX = (log10(Mcl)*0.360268003) + 0.313342031;
			MMAX = pow(10.0,MMAX);
		}
		if (MMAX > mup) MMAX = mup;
	} else {
		MMAX = upper_IMF_limit;
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

	if (mfunc == 1) {
		printf("\nMaximum stellar mass set to: %.2f\n",MMAX);
		generate_m1(&N, star, mlow, mup, &M, &mmean, MMAX, Mcl);
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
			generate_m2(an, mlim, alpha, Mcl, M_tmp, subcount, &N, &mmean, &M, star, MMAX);
		}
 	} else {
		printf("\nSetting stellar masses to %.1f solar mass\n", single_mass);
		if (!N) N = Mcl/single_mass;
		for (j=0;j<N;j++) star[j][0] = single_mass;
		mmean = single_mass;
		M = N*mmean;
		printf("\nM = %g\n", M);
	}
	for (i=0;i<N;i++) star[i][0] /= M; //scale masses to Nbody units
		
	
	
	/*************************************
	 * Generate positions and velocities *
	 *************************************/
	
	printf("\n\n-----GENERATE POSITIONS & VELOCITIES-----   \n"); 

	//evaluate approximate tidal radius assuming circular orbit
	if (tf == 3) {
		//in the case of Allen & Santillan potential, assume kappa = 1.4omega (eq. 9 in Kuepper et al. 2009)
		omega = sqrt(VG[0]*VG[0]+VG[1]*VG[1]+VG[2]*VG[2])/sqrt(RG[0]*RG[0]+RG[1]*RG[1]+RG[2]*RG[2]);
		rtide = pow(G*M/(2.0*omega*omega),1.0/3.0);
	} else { 
		//in the case of a point mass potential or near field approximation
		rtide = pow(1.0*M/(3.0*M1pointmass),1.0/3.0)*sqrt(RG[0]*RG[0]+RG[1]*RG[1]+RG[2]*RG[2]);
	}
	printf("\nApproximate tidal radius: %g (pc)\n", rtide);

	
	//generate scaled pos & vel
	if (profile == 1) {
		double rhtemp;
		printf("\nGenerating King model with parameters: N = %i\t W0 = %g\t Rh = %.3f\n",N, W0, Rh);
		generate_king(N, W0, star, &rvir, &rhtemp, &rking);
		printf ("\nrvir = %.5f\t rh = %.5f\t rking = %.5f\t rtide = %.5f (pc)\n", rvir/rhtemp*Rh, Rh, rking/rhtemp*Rh, rtide);
		rvir = rvir/rhtemp*Rh;
	} else {
		printf("\nGenerating Plummer model with parameters: N = %i\t Rh = %.3f\n",N,Rh);
		rvir = Rh/0.772764;
		rplummer = Rh/1.305;
		generate_plummer(N, star, rtide, rvir);
		printf ("\nrvir = %.5f\t rh = %.5f\t rplummer = %.5f\t rtide = %.5f (pc)\n", rvir, Rh, rplummer, rtide);
	}
		
	
	//CoM correction
	printf("\nApplying centre-of-mass correction...\n");
	for (j=0; j<7; j++) cmr[j] = 0.0;
	
	for (j=0; j<N; j++) {
		for (i=1;i<7;i++) 
			cmr[i] += star[j][0]*star[j][i]; 
	} 
		
	for (j=0; j<N; j++) {
		for (i=1;i<7;i++)
			star[j][i] -= cmr[i];
	}
	
	
	//estimate NNBMAX and RS0 (Nbody6 only)
	if ((code == 0) || (code == 2)) {
		float rarray[N];
		for (j=0; j<N; j++) {
			rarray[j] = sqrt(star[j][1]*star[j][1]+star[j][2]*star[j][2]+star[j][3]*star[j][3]);
		}
		qsort(rarray, N, sizeof(rarray[0]), (void *)cmpmy);
		NNBMAX = sqrt(N);
		RS0 = rarray[NNBMAX];
		printf("\nEstimating appropriate NNBMAX = %i and RS0 = %f\n",NNBMAX,RS0);
	}

		
	//generate binaries
	if ((!fbin) && (!nbin)) {
		printf("\nNo primordial binaries!\n");
		bin = 3; //KZ(22)
	} else {
		if (!nbin) nbin = 0.5*N*fbin;
		printf("\nCreating %i primordial binary systems, fraction: %6.2f percent.\n", nbin, 2.0*nbin/N*100);
		tscale = sqrt(rvir*rvir*rvir/(G*M));
		get_binaries(nbin, star, M, rvir, pairing, N, adis, amin, amax);
		bin = 4; //KZ(22)
	} 
	

	
	/***********
	 * Scaling * 
	 ***********/

	printf("\n\n-----SCALING-----      \n"); 
		
	//scale masses, pos & vel to astrophysical units or Nbody units
	if (units) {
		tscale = sqrt(rvir*rvir*rvir/(G*M));
		
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
	
	
	/**********
	 * Output * 
	 **********/
	
	printf("\n\n-----OUTPUT-----      \n"); 

	if (code == 0) 
		output0(output, N, NNBMAX, RS0, dtadj, dtout, tcrit, rvir, mmean, tf, regupdate, etaupdate, mloss, bin, esc, M, mlow, mup, MMAX, epoch, dtplot, Z, nbin, Q, RG, VG, rtide, gpu, star);
	else if (code == 1)
		output1(output, N, dtadj, dtout, tcrit, rvir, mmean, tf, regupdate, etaupdate, mloss, bin, esc, M, mlow, mup, MMAX, epoch, Z, nbin, Q, RG, VG, rtide, gpu, star);
	else if (code == 2)
		output2(output, N, NNBMAX, RS0, dtadj, dtout, tcrit, rvir, mmean, tf, regupdate, etaupdate, mloss, bin, esc, M, mlow, mup, MMAX, epoch, dtplot, Z, nbin, Q, RG, VG, rtide, gpu, star);
	
	
	
	/**********************
	 * Final energy check *
	 **********************/
	
	printf("\n\n-----FINISH-----  \n"); 

	if (check) {
		printf("\nMaking final energy check... (may take a while but can be aborted by pressing CTRL+c)\n");
				
		for (j=0; j<N; j++) {
			ekin += star[j][0]*((star[j][4]*star[j][4])+(star[j][5]*star[j][5])+(star[j][6]*star[j][6]));
			if (j) {
				for (i=0;i<j-1;i++) 
					epot -= star[i][0]*star[j][0]/sqrt((star[i][1]-star[j][1])*(star[i][1]-star[j][1])+(star[i][2]-star[j][2])*(star[i][2]-star[j][2])+(star[i][3]-star[j][3])*(star[i][3]-star[j][3])); 
			}
			sigma += star[j][4]*star[j][4]+star[j][5]*star[j][5]+star[j][6]*star[j][6];
		} 

		if (units) epot *= G;
		ekin *= 0.5;

		sigma = sqrt(sigma/N);
		tscale = sqrt(rvir*rvir*rvir/(G*M));
		
		printf("\nEkin = %g\t Epot = %g\t Etot = %g ", ekin, epot, ekin+epot);
		printf("\nVel.Disp. = %g\tCross.Time = %g \n", sigma, 2.0/sigma);

		if (units) printf("Vel.Disp. = %g\tCross.Time = %g (Nbody units)\n", sigma/rvir*tscale, 2.0/sigma/tscale);
		else printf("Vel.Disp. = %g\tCross.Time = %g (physical units, km/s, Myr)\n", sigma*rvir/tscale, 2.0/sigma*tscale);
	}
	
	for (j=0;j<NMAX;j++) free (star[j]);
	free(star);
	
	t2 = clock();														//stop stop-watch
	printf("\nElapsed time: %g sec\n",(double)(t2-t1)/CLOCKS_PER_SEC);	//print stopped time

	return 0;
}




/*************
 * Functions *
 *************/

int generate_m1(int *N, double **star, double mlow, double mup, double *M, double *mmean, double MMAX, double Mcl) {
	int ty, i;
	double alpha1, alpha2, c1, c2, k1, k2, xx, mth;

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
		printf("Estimated number of necessary stars: %i\n", *N);
	}
	
	
	c1 = 1.0-alpha1;  
	c2 = 1.0-alpha2; 
	*mmean = 0.0;				
	*M = 0.0;
	double mostmassive = 0.0;
	
	for (i=0; i<*N; i++) {
		do {
			xx = drand48();		
			if (xx<k1/k2)   
				star[i][0] = pow(0.5*c1*xx*k2+pow(mlow,c1),1.0/c1);
			else 
				star[i][0] = pow(c2*(xx*k2-k1)+pow(max(0.5,mlow),c2),1.0/c2);
		} while (star[i][0] > MMAX);

		if (star[i][0] > mostmassive) mostmassive = star[i][0];
		*M += star[i][0];
		if ((i==*N-1) && (*M<Mcl)) *N += 1;
	}

	printf("Total mass: %g\t(%i stars)\n",*M,*N);
	*mmean = *M/ *N;
	printf("Most massive star: %g\n",mostmassive);
	printf("Mean masses theoretical/data: %f %f\n",mth,*mmean);
	
	return 0;
}

int generate_m2(int an, double *mlim, double *alpha, double Mcl, double M_tmp, double *subcount, int *N, double *mmean, double *M, double **star, double MMAX) {
	int i, j;
	double tmp, ml, mup;
	double mostmassive = 0.0;
	*mmean = 0.0;				
	*M = 0.0;
	if (!*N) *N = 1;
	
	printf("%f\n", MMAX);
	
	for (i = 0; i < an; i++) {
			printf("# <%.2f , %.2f> .. %.2f\n", mlim[i], mlim[i+1], alpha[i]);
	}
	
	for (i = 1; i < an; i++) subcount[i] += subcount[i-1];
	
	for (i = 0; i < *N; i++) {
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

		if (star[i][0] > mostmassive) mostmassive = star[i][0];
		*M += star[i][0];
		if ((i==*N-1) && (*M<Mcl)) *N += 1;
	}

	printf("Total mass: %g\t(%i stars)\n",*M,*N);
	*mmean = *M/ *N;
	printf("Most massive star: %g\n",mostmassive);
	printf("Mean mass: %f\n",*mmean);
	
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

int generate_plummer(int N, double **star, double rtide, double rvir){
	int i,j;
	double a[9], ri, sx, sv, rcut, r2, rvirscale;
	double ke = 0.0;
	double pe = 0.0;
	
	//Scale length
	sx = 1.5*TWOPI/16.0;
	sv = sqrt(1.0/sx);
	
	printf("Setting cut-off radius of Plummer sphere to approximate tidal radius\n");
	rcut = rtide/(sx*rvir);		//cut-off radius for Plummer sphere = tidal radius in scaled length

	printf("\nGenerating Stars:\n");	

	for (i=0;i<N;i++) {
		
		if ((i/1000)*1000 == i) printf("Generating star #%i\n", i);
			
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

		if (i) {
			for (j=0;j<i-1;j++) {
				r2 = (star[i][1]-star[j][1])*(star[i][1]-star[j][1]) + (star[i][2]-star[j][2])*(star[i][2]-star[j][2]) +(star[i][3]-star[j][3])*(star[i][3]-star[j][3]) ;
				if (r2<0.000001) printf("WARNING: %i %i\tdr = %lf\n", i,j,sqrt(r2));
				pe -=  star[i][0]*star[j][0]/sqrt(r2);
			}
		}
		
		ke += star[i][0]*(pow(star[i][4],2)+pow(star[i][5],2)+pow(star[i][6],2));
	}

	rvirscale = -GNBODY*pow(MNBODY,2)/(2.0*pe);
	sx = 1.0/rvirscale;
	
	ke *= 0.5;
	sv = sqrt(4.0*ke);
		
	// Scale coordinates to N-body units
	for (i=0;i<N;i++) {
		star[i][1] *= sx;
		star[i][2] *= sx;
		star[i][3] *= sx;
		star[i][4] /= sv;
		star[i][5] /= sv;
		star[i][6] /= sv;
	}
	
	return 0;
}

int generate_king(int N, double W0, double **star, double *rvir, double *rh, double *rking){
	
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
	double vmax, vscale, speed, fstar, ustar, vstar, wstar, mstar;
	double coord[N][6];
	double pe = 0.0;
	double ke = 0.0;
	
	srand48((unsigned) time(NULL));
	
	
	
	
	
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
	
	h = pow(W0-rmin,2)/(M-1);	//step size for interpolation
	den = densty(W0);			//central density
	
	//x is King's W, y[0] is z**2, y[1] is 2*z*dz/dW, where z is scaled radius
	//so x(y[0]) is energy as function of radius2 and x(y[1]) is the derivative
	
	xstart = 0.000001;
	ystart0 = 2.0*xstart/3.0;
	ystart1 = -2.0/3.0;
	
	x1 = W0 - xstart;
	x2 = 0.0;
	
	//integrate Poisson's eqn
	printf("Integrating Poisson's equation\n",kount);

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
			pot = - 0.6*pow(mass[k],2)/sqrt(yking[k][0]);  //<<--- why is this 0.6 and not 0.5???
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
	
	for (i=0;i<N;i++) {
		
		if ((i/1000)*1000 == i) printf("Generating stars #%i\n", i);
		
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
		
		
		/****************
		 * CHOOSE SPEED *
		 ****************/
		
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
		
		//critical neighbour check
		if (i) {
			for (j=0;j<i-1;j++) {
				r2 = (coord[j][0]-xstar)*(coord[j][0]-xstar) + (coord[j][1]-ystar)*(coord[j][1]-ystar) + (coord[j][2]-zstar)*(coord[j][2]-zstar);
				if (r2<0.000001) printf("WARNING: %i %i\tdr = %lf\n", i,j,sqrt(r2));
				pe -=  mstar*star[j][0]/sqrt(r2);//1.0/sqrt(r2);
			}
		}
		ke += mstar*pow(speed,2);
		
	}

	pe = GNBODY*pe;//*pow(mstar,2);
	ke *= 0.5;
	*rvir = -GNBODY*pow(MNBODY,2)/(2.0*pe);
	vscale = sqrt(4.0*ke);
	
	
	/**********
	 * OUTPUT *
	 **********/
	
	printf("\nActual values:\n");
	printf("Edge radius (King units) = %g\t(Nbody units) = %g\n", *rking, *rking/ *rvir);
	printf("Core radius (King units) = %g\t(Nbody units) = %g\n\n", 1.0, 1.0/ *rvir);
	printf("Concentration = %g\n", log10(*rking));

	for (i=0;i<N;i++) {
		for (j=0;j<3;j++) {
			star[i][j+1] = coord[i][j]/ *rvir;
			star[i][j+4] = coord[i][j+3]/vscale;
		}
		//printf("%f  %f  %f  %f  %f  %f  %f\n", mstar, coord[i][0],coord[i][1],coord[i][2],coord[i][3],coord[i][4],coord[i][5] );
	}
	
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
	
	y[0] = ystart0;  //z2
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

int cmpmy(float *x1, float *x2) { //smallest to largest
	if(*x1<*x2) return -1;
	return 1;
}

int cmpmy_reverse(float *x1, float *x2) { //largest to smallest
	if(*x1>*x2) return -1;
	return 1;
}

int get_binaries(int nbin, double **star, double M, double rvir, int pairing, int N, int adis, double amin, double amax) {		
	int i, j;
	double m1, m2, ecc, abin;
	double pmat[3][2], rop[2], vop[2], rrel[3], vrel[3];
	double ea, mm, eadot, cosi, inc, peri, node;
	double Pmin = 10, delta = 45, eta = 2.5;
	double lP, P;	

	//Specify component pairing
	if (pairing) {
		order(star, N, M);
		printf("\nApplying ordered pairing for stars with masses > 5 Msun.\n");
	}
	
	for (i=0; i < nbin; i++) {
		//Specify component masses
		m1 = star[2*i][0];
		m2 = star[2*i+1][0];
		
		//Specify semi-major axis
		if (adis == 1) {
			//flat semi-major axis distribution
			if (!i) printf("\nApplying flat semi-major axis distribution with amin = %g and amax = %g.\n", amin, amax);
			if (!i) amin /= rvir;
			if (!i) amax /= rvir;
			abin = amin+drand48()*(amax-amin);
		} else {
			//derive from Kroupa (1995) period distribution
			if (!i) printf("\nDeriving semi-major axis distribution from Kroupa (1995) period distribution.\n");
			do {
				lP = log10(Pmin) + sqrt(delta*(exp(2.0*drand48()/eta) - 1.0));
			} while (lP > 8.43);
			
			P = pow(10,lP);//days
			P /= 365.25;//yr
			
			abin = pow((m1+m2)*M*P*P,(1.0/3.0));//AU
			abin /= 206264.806;//pc
			abin /= rvir;//Nbody units
		}			

		
		//Specify eccentricity distribution
		if (!i) printf("\nApplying thermal eccentricity distribution.\n");
		ecc = sqrt(drand48());   // Thermal distribution f(e)=2e 
		
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
	}
	
	return 0;

}

int order(double **star, int N, double M){
	int i;
	float masses[N][2];
	double mass_temp;
	
	//temporary scaling to astrophysical units
	for (i=0;i<N;i++) {
		star[i][0] *= M;
	}	
	
	for (i=0;i<N;i++) {
		masses[i][0] = star[i][0];
		masses[i][1] = i;
	}
	qsort(masses, N, sizeof(masses[0]), (void *)cmpmy_reverse);	
	
	for (i=0;i<N;i++) {
		if (masses[i][0] > 5.0) {
			mass_temp = star[i][0];
			star[i][0] = masses[i][0];
			star[(int) masses[i][1]][0] = mass_temp;
		}
	}

	//scaling back to Nbody units
	for (i=0;i<N;i++) {
		star[i][0] /= M;
	}	
	
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

int output0(char *output, int N, int NNBMAX, double RS0, double dtadj, double dtout, double tcrit, double rvir, double mmean, int tf, int regupdate, int etaupdate, int mloss, int bin, int esc, double M, double mlow, double mup, double MMAX, double epoch, double dtplot, double Z, int nbin, double Q, double *RG, double *VG, double rtide, int gpu, double **star){

	//Open output files
	char PARfile[20], NBODYfile[20];		
	FILE *PAR, *NBODY;
	sprintf(PARfile, "%s.input",output);
	PAR = fopen(PARfile,"w");
	sprintf(NBODYfile, "%s.fort.10",output);
	NBODY = fopen(NBODYfile,"w");
	
	//write to .PAR file	
	fprintf(PAR,"1 50000.0 0\n");
	fprintf(PAR,"%i 1 10 199 %i 1\n",N,NNBMAX);
	fprintf(PAR,"0.02 0.02 %.8f %.8f %.8f %.8f 1.0E-03 %.8f %.8f\n",RS0,dtadj,dtout,tcrit,rvir,mmean);
	fprintf(PAR,"2 2 1 0 1 0 2 0 0 2\n");
	fprintf(PAR,"0 0 0 %i 2 %i %i 0 %i 2\n",tf,regupdate,etaupdate,mloss);
	fprintf(PAR,"0 %i %i 0 1 2 0 0 0 2\n",bin,esc);
	fprintf(PAR,"0 0 0 0 1 0 0 2 0 3\n");
	fprintf(PAR,"1.0E-5 1.0E-4 0.01 1.0 1.0E-06 0.01\n");
	fprintf(PAR,"2.350000 %.8f %.8f %i 0 %.8f %.8f %.8f\n",mlow,MMAX,nbin,Z,epoch,dtplot);
	fprintf(PAR,"%.2f 0.0 0.0 0.00000\n",Q);
	if (tf == 2) {
		fprintf(PAR,"%.8e %.8f\n",M1pointmass,sqrt(RG[0]*RG[0]+RG[1]*RG[1]+RG[2]*RG[2]));
	} else if (tf == 3) {
		fprintf(PAR,"%.6e %.6e %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",M1allen,M2allen,a2allen,b2allen,VCIRC,RCIRC,RG[0]/1000.0,RG[1]/1000.0,RG[2]/1000.0,VG[0],VG[1],VG[2]);
	}
	fprintf(PAR,"0.0 0.0 0.0 0.0\n");
	if (gpu) fprintf(PAR,"1.0\n");
	
	
	
	//write to .NBODY file
	int j;
	for (j=0;j<N;j++) {
		fprintf(NBODY,"%.8f\t%.8f %.8f %.8f\t%.8f %.8f %.8f\n",star[j][0],star[j][1],star[j][2],star[j][3],star[j][4],star[j][5],star[j][6]);
	}
	
	fclose(PAR);
	fclose(NBODY);
	printf("\nData written to %s and %s\n", PARfile, NBODYfile);
	
	return 0;
	
}

int output1(char *output, int N, double dtadj, double dtout, double tcrit, double rvir, double mmean, int tf, int regupdate, int etaupdate, int mloss, int bin, int esc, double M, double mlow, double mup, double MMAX, double epoch, double Z, int nbin, double Q, double *RG, double *VG, double rtide, int gpu, double **star){

	//Open output files
	char PARfile[20], NBODYfile[20];		
	FILE *PAR, *NBODY;
	sprintf(PARfile, "%s.PAR",output);
	PAR = fopen(PARfile,"w");
	sprintf(NBODYfile, "%s.NBODY",output);
	NBODY = fopen(NBODYfile,"w");
	
	//write to .PAR file	
	fprintf(PAR,"1 50000.0 0\n");
	fprintf (PAR,"%i 1 10 3 8\n",N);
	fprintf(PAR,"0.02 %.8f %.8f %.8f 100.0 1000.0 1.0E-02 %.8f %.8f\n",dtadj,dtout,tcrit,rvir,mmean);
	fprintf(PAR,"2 2 1 0 1 0 2 0 0 2\n");
	fprintf(PAR,"0 0 0 %i 2 %i %i 0 %i 2\n",tf,regupdate,etaupdate,mloss);
	fprintf(PAR,"0 bin %i 0 1 2 0 0 0 2\n",bin, esc);
	fprintf(PAR,"0 0 0 0 1 0 0 2 0 3\n");
	fprintf(PAR,"0 0 0 0 0 0 0 0 0 0\n");
	fprintf(PAR,"%.8f %.8f %.8f\n",M,mlow,mup);
	fprintf(PAR,"1.0E-5 1.0E-4 0.01 1.0 1.0E-06 0.01\n");
	fprintf(PAR,"2.350000 %.8f %.8f %i %.8f %.8f 100000.0\n",mlow,MMAX, nbin, Z, epoch);	
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
		fprintf(NBODY,"%.8f\t%.8f %.8f %.8f\t%.8f %.8f %.8f\n",star[j][0],star[j][1],star[j][2],star[j][3],star[j][4],star[j][5],star[j][6]);
	}
	
	fclose(PAR);
	fclose(NBODY);
	printf("\nData written to %s and %s\n", PARfile, NBODYfile);
	
	return 0;
	
}

int output2(char *output, int N, int NNBMAX, double RS0, double dtadj, double dtout, double tcrit, double rvir, double mmean, int tf, int regupdate, int etaupdate, int mloss, int bin, int esc, double M, double mlow, double mup, double MMAX, double epoch, double dtplot, double Z, int nbin, double Q, double *RG, double *VG, double rtide, int gpu, double **star){

	//Open output files
	char PARfile[20], NBODYfile[20];		
	FILE *PAR, *NBODY;
	sprintf(PARfile, "%s.PAR",output);
	PAR = fopen(PARfile,"w");
	sprintf(NBODYfile, "%s.NBODY",output);
	NBODY = fopen(NBODYfile,"w");
	
	//write to .PAR file	
	fprintf(PAR,"1 50000.0 0\n");
	fprintf(PAR,"%i 1 10 199 %i 1\n",N,NNBMAX);
	fprintf(PAR,"0.02 0.02 %.8f %.8f %.8f %.8f 1.0E-03 %.8f %.8f\n",RS0,dtadj,dtout,tcrit,rvir,mmean);
	fprintf(PAR,"2 2 1 0 1 0 2 0 0 2\n");
	fprintf(PAR,"0 0 0 %i 2 %i %i 0 %i 2\n",tf,regupdate,etaupdate,mloss);
	fprintf(PAR,"0 %i %i 0 1 2 0 0 0 2\n",bin, esc);
	fprintf(PAR,"0 0 0 0 1 0 0 2 0 3\n");
	fprintf(PAR,"1.0E-5 1.0E-4 0.01 1.0 1.0E-06 0.01\n");
	fprintf(PAR,"2.350000 %.8f %.8f %i 0 %.8f %.8f %.8f\n",mlow,MMAX,nbin,Z,epoch,dtplot);
	fprintf(PAR,"%.2f 0.0 0.0 0.00000\n",Q);
	if (tf == 1) {
		rtide = pow(1.0*M/(3.0*M1pointmass),1.0/3.0)*sqrt(RG[0]*RG[0]+RG[1]*RG[1]+RG[2]*RG[2]);
		fprintf(PAR,"%i %.8f\n",0,sqrt(1.0/(3.0*pow(rtide,3))));
	} else if (tf == 2) {
		fprintf(PAR,"%.8e %.8f\n",M1pointmass,sqrt(RG[0]*RG[0]+RG[1]*RG[1]+RG[2]*RG[2]));
	} else if (tf == 3) {
		fprintf(PAR,"%.6e %.6f %.6e %.6f %.6f %.6e %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",M1allen,b1allen,M2allen,a2allen,b2allen,M3allen,a3allen,RG[0]/1000.0,RG[1]/1000.0,RG[2]/1000.0,VG[0],VG[1],VG[2]);
	}
	fprintf(PAR,"0.0 0.0 0.0 0.0\n");
	if (gpu) fprintf(PAR,"1.0\n");
	
	
	
	//write to .NBODY file
	int j;
	for (j=0;j<N;j++) {
		fprintf(NBODY,"%.8f\t%.8f %.8f %.8f\t%.8f %.8f %.8f\n",star[j][0],star[j][1],star[j][2],star[j][3],star[j][4],star[j][5],star[j][6]);
	}
	
	fclose(PAR);
	fclose(NBODY);
	printf("\nData written to %s and %s\n", PARfile, NBODYfile);
	
	return 0;
	
}

int help() {
	printf("\n Usage: mcluster -[N|M|P|W|R|T|Q|C|A|O|G|o|f|a|m|B|b|S|t|X|V|h|?]  \n");
	printf("                                                                     \n");
	printf("       -N <number> (number of stars)                                 \n");
	printf("       -M <value> (mass of cluster; specify either N or M)           \n");
	printf("       -P <1|2> (density profile; 1= Plummer, 2= King 1966)          \n");
	printf("       -W <1-12> (W0 parameter for King 1966 profile)                \n");
	printf("       -R <value> (half-mass radius [pc])                            \n");
	printf("       -T <value> (tcrit in N-body units)                            \n");
	printf("       -Q <value> (virial ratio)                                     \n");
	printf("       -C <1|2> (code; 1= Nbody6, 2= Nbody4)                         \n");
	printf("       -A <value> (dtadj in N-body units)                            \n");
	printf("       -O <value> (deltat in N-body units)                           \n");
	printf("       -G <0|1> (GPU usage; 0= no GPU, 1= use GPU)                   \n");
	printf("       -o <name> (name of cluster model)                             \n");
	printf("       -f <0|1|2> (IMF; 0= no IMF, 1= Kroupa (2001), 2= user defined)\n");
	printf("       -a <value> (IMF slope for user defined IMF, may be used       \n"); 
	printf("                   multiple times, from low mass to high mass)       \n");
	printf("       -m <value> (IMF mass limit for user defined IMF, may be used  \n");
	printf("                   multiple times, from low mass to high mass [Msun])\n");
	printf("       -B <number> (number of binary systems)                        \n");
	printf("       -b <value> (binary fraction, specify either B or b)           \n");
	printf("       -S <number> (seed for randomization; 0= randomize by timer)   \n");
	printf("       -t <1|2|3> (tidal field; 1= near-field, 2= point-mass,        \n");
	printf("                   3= Milky-Way potential)                           \n");
	printf("       -X <value> (component of galactocentric radius vector, use 3x)\n");
	printf("       -V <value> (component of cluster velocity vector, use 3x)     \n");
	printf("       -h (display this help)                                        \n");
	printf("       -? (display this help)                                        \n");
	printf("                                                                     \n");
	printf(" Examples: mcluster -N 1000 -R 0.8 -P 1 -f 1 -B 100 -o test1         \n");
	printf("           mcluster -f 2 -m 0.08 -a -1.35 -m 0.5 -a -2.7 -m 100.0    \n");
	printf("           mcluster -T 3 -X 8500 -X 0 -X 0 -V 0 -V 220 -V 0          \n");
	printf("                                                                     \n");

	return 0;
}

