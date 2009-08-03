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

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<string.h>
#include<sys/stat.h>

//Constants
#define G 0.0043009211 //in pc*km^2*s^-2*Msun
#define PI   4.0*atan(1.0)   /* PI */
#define TWOPI   8.0*atan(1.0)   /* 2PI */
#define GKING 1.0
#define MTOTKING 1.0

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
#define VCIRC 220			//circular velocity at RCIRC
#define RCIRC 8.500			//kpc

//Point-mass potential - constants:
#define M1pointmass 9.565439E+10    //solar masses


int generate_m(int *N, double **kubus, double mlow, double mup, double *M, double *mmean, double MMAX, double Mcl);
int generate_plummer(int N, double **kubus, double rtide, double rvir);
int generate_king(int N, double W0, double **kubus, double *rvir, double *rh, double *rking); 
double densty(double z);
int odeint(double ystart0, double ystart1, double x1, double x2, double den, int *kount, double *xp, double **yp, int M, int KMAX);	 
int derivs(double x, double *y, double *dydx, double den);
int rkqc(double *y,double *dydx, double *x, double *h,double den, double *yscal, double *hdid, double *hnext, double TOL);	
int rk4(double x, double *y, double *dydx, double h, double *yout, double den);



int main (int argc, const char * argv[]) {

	/*******************
	 * Input variables *
	 *******************/
	
	char *output = "test";			//name of output files
	unsigned int seed = 0;			//number seed for random number generator; =0 for randomization by local time

	//Physical parameters
	int N = 0;					    //number of stars, Mcl will be set to 0 if specified!
	double Mcl = 5000.0;            //total mass of the cluster, only used when N is set to 0, necessary for usage of maximum stellar mass relation of Weidner & Kroupa 2007
	int profile = 0;				//density profile; =0 Plummer sphere, =1 King profile
	double W0 = 5.0;				//King's W0 paramter [0.3-12.0]
	double Rh = 3.00;				//Half-mass radius [pc]
	double tcrit = 500.0;			//Simulation time [Tcr (Myr in Nbody6 custom)]
	int tf = 3;						//tidal field: =1 Near-field approximation, =2 point-mass galaxy, =3 Allen & Santillan (1991) MW potential (or Sverre's version of it)
	double RG[3] = {8500.0,0.0,0.0}; //Initial Galactic coordinates of the cluster [pc]
	double VG[3] = {0.0,220.0,0.0};  //Initial velocity of the cluster [km/s]

	//Code parameters
	int code = 0;					//Nbody version: =0 Nbody6, =1 Nbody4, =2 Nbody6 custom
	double dtadj = 10.0;			//DTADJ [Tcr (Myr in Nbody6 custom)], energy-check time step
	double dtout = 50.0;			//DTOUT [Tcr (Myr in Nbody6 custom)], output interval, must be multiple of DTADJ
	double dtplot = 100.0;			//DTPLOT [Myr], output of HRdiagnostics, should be multiple of DTOUT
	int gpu = 1;					//Use of GPU, 0= off, 1= on
	double RS0 = 1.0;				//Initial radius of neighbour sphere [pc], Nbody6 only
	int regupdate = 1;				//Update of regularization parameters during computation; 0 = off, 0 > on
	int etaupdate = 1;				//Update of ETAI & ETAR during computation; 0 = off, 0 > on
	int esc = 1;					//Removal of escapers; 0 = no removal, 1 = regular removal at 2*R_tide
	
	//Mass function parameters
	int mfunc = 1;					//0 = single mass stars; 1 = use Kroupa (2001) mass function
	double mlow = 0.1;				//lower mass limit
	double mup = 1.2;				//upper mass limit (will be overwritten when N is set to 0 and weidner is set to 1 with value of Weidner & Kroupa 2007)
	int weidner = 1;				//Usage of Weidner & Kroupa 2007 relation for most massive star; =0 off, =1 on
	int mloss = 3;					//Stellar evolution; 0 = off, 3 = Eggleton, Tout & Hurley [KZ19]
	double epoch = 0.0;				//Star burst has been ... Myr before [e.g. 1000.0, default = 0.0]
	double Z = 0.02;				//Metallicity [0.0001-0.03, 0.02 = solar]
	double FeH = -1.41;				//Metallicity [Fe/H], only used when Z is set to 0
	int prantzos = 1;				//Usage of Prantzos 2007 relation for the life-times of stars. Set upper mass limit to Lifetime(mup) >= epoch
	// Lifetime(Mstar) = 1.13E4*pow(Mstar,-3)+0.6E2*pow(Mstar,-0.75)+1.2; //Myr	, Prantzos 2007
	
	//Mcluster internal parameters
	int check = 1;					//Make energy check at end of mcluster; =0 off, =1 on
	double Zsun = 0.02;				//Solar metallicity
	int NMAX = 500000;				//Maximum number of stars allowed in mcluster
	double upper_IMF_limit = 150.0; //Maximum stellar mass allowed in mcluster [Msun]
	
	
	/*********
	 * Start *
	 *********/
	
	clock_t t1, t2;
	t1 = clock();							//start stop-watch
  
	if (seed) srand48(seed);				//initialize random number generator by seed
	else srand48((unsigned) time(NULL));	//initialize random number generator by local time
	int i,j;
	double M;								//Total mass [M_sun]
	double mmean;							//Mean stellar mass [M_sun]
	int NNBMAX;								//Maximum neighbour number (Nbody6 only)
	double rtide;							//Tidal radius [pc]
	double omega;							//Angular velocity of cluster around the galaxy
	double rvir;							//Virial radius [pc]
	double cmr[7];							//For CoM correction
	double rking, rplummer;					//King-, Plummer radius
	double MMAX;							//most massive star
	
	
	if ((Mcl) && (N)) {
		printf("\nCAUTION: set either Mcl or N to 0!\nSet Mcl to 0 by default...\n\n");
		Mcl = 0.0;
	}
		
	//Open output files
	char PARfile[20], NBODYfile[20];		
	FILE *PAR, *NBODY;
	sprintf(PARfile, "%s.PAR",output);
	PAR = fopen(PARfile,"w");
	sprintf(NBODYfile, "%s.NBODY",output);
	NBODY = fopen(NBODYfile,"w");
	
	
	
    /***********************
	 * Generate star array *
	 ***********************/

	int columns = 7;
	double **kubus;
	kubus = (double **)calloc(NMAX,sizeof(double *));
	for (j=0;j<NMAX;j++){
		kubus[j] = (double *)calloc(columns,sizeof(double));
		if (kubus[j] == NULL) {
			printf("\nMemory allocation failed!\n");
			return 0;
		}
	}

	
	
	/*******************************************	
	 * Evaluate Z from [Fe/H] if Z is set to 0 *
	 *******************************************/
	
	 if (!Z) {
		Z = pow(10.0, 0.977*FeH)*Zsun; //Bertelli, Bressan, Chiosi, Fagotto, Nasi, 1994, A&AS, 106, 275
		printf("Using Bertelli et al. (1994) relation to convert FeH = %.3f into Z = %.3f\n", FeH, Z);
	}
	
	
	
	/**********************************
	 * Calculate maximum stellar mass *
	 **********************************/
		 
	if (!N && weidner && mfunc) {
		mup = upper_IMF_limit;  
		printf("Using Weidner & Kroupa (2007) relation for upper stellar mass limit\n");
		
		if (Mcl < 1000.0) {
			MMAX = (log10(Mcl)*0.540563175) - 0.14120167;
			MMAX = pow(10.0,MMAX);
		} else if (Mcl < 3300.0) {
			MMAX = (log10(Mcl)*0.19186051) + 0.9058611;
			MMAX = pow(10.0,mup);
		} else {
			MMAX = (log10(Mcl)*0.360268003) + 0.313342031;
			MMAX = pow(10.0,MMAX);
		}
		if (MMAX > mup) MMAX = mup;
	} else {
		MMAX = mup;
	}
	
	if (mfunc && epoch && prantzos) {
		printf("Using Prantzos (2007) relation to reduce upper mass limit to Lifetime(mup) > epoch\n");
		while (Lifetime(MMAX) < sqrt(pow(epoch,2))) {
			MMAX -= 0.01;
		}
	}
	

	
	/*******************
	 * Generate masses *
	 *******************/
	
	if (mfunc) {
		printf("Maximum stellar mass set to: %.2f\n",MMAX);
		generate_m(&N, kubus, mlow, mup, &M, &mmean, MMAX, Mcl);
	} else {
		printf("Setting stellar masses to 1 solar mass\n");
		for (j=0;j<N;j++) kubus[j][0] = 1.0/N;
		mmean = 1.0;
		M = N*mmean;
		printf("M = %g\n", M);
	}
	
	
	
	/*************************************
	 * Generate positions and velocities *
	 *************************************/
		
	//evaluate approximate tidal radius assuming circular orbit
	if (tf == 3) {//in the case of Allen & Santillan potential, assume kappa = 1.4omega (eq. 9 in Kuepper et al. 2009)
		omega = 220.0/sqrt(RG[0]*RG[0]+RG[1]*RG[1]+RG[2]*RG[2]);
		rtide = pow(G*M/(2.0*omega*omega),1.0/3.0);
	} else { //in the case of a point mass potential or near field approximation
		rtide = pow(1.0*M/(3.0*M1pointmass),1.0/3.0)*sqrt(RG[0]*RG[0]+RG[1]*RG[1]+RG[2]*RG[2]);
	}
	printf("Approximate tidal radius: %g\n", rtide);

	
	//generate scaled pos & vel
	if (profile == 1) {
		double rhtemp;
		generate_king(N, W0, kubus, &rvir, &rhtemp, &rking);
		printf ("rvir = %.5f\t rh = %.5f\t rking = %.5f\t rtide = %.5f (pc)\n", rvir/rhtemp*Rh, Rh, rking/rhtemp*Rh, rtide);
		rvir = rvir/rhtemp*Rh;
	} else {
		printf("\nGenerating Plummer model with parameters: N = %i\t Rh = %.3f\n\n",N,Rh);
		rvir = Rh/0.772764;
		rplummer = Rh/1.305;
		generate_plummer(N, kubus, rtide, rvir);
		printf ("rvir = %.5f\t rh = %.5f\t rplummer = %.5f\t rtide = %.5f (pc)\n", rvir, Rh, rplummer, rtide);
	}
	
	
	//CoM correction
	for (j=0; j<7; j++) cmr[j] = 0.0;
	
	for (j=0; j<N; j++) {
		for (i=1;i<4;i++) 
			cmr[i] += kubus[j][0]*kubus[j][i]; 
	} 
		
	for (j=0; j<N; j++) {
		for (i=1;i<4;i++)
			kubus[j][i] -= cmr[i];
	}
	
	
	
	/**********
	 * Output *
	 **********/
		
	//write to .PAR file	
	fprintf(PAR,"1 50000.0 0\n");
	if (code == 1) {
		fprintf (PAR,"%i 1 10 3 8\n",N);
		fprintf(PAR,"0.02 %.8f %.8f %.8f 100.0 1000.0 1.0E-02 %.8f %.8f\n",dtadj,dtout,tcrit,rvir,mmean);
	} else if (code == 2) {
		NNBMAX = sqrt(N); //Guess for the neighbour number (Nbody6 only)
		fprintf(PAR,"%i 1 10 199 %i 1\n",N,NNBMAX);
		fprintf(PAR,"0.02 0.03 %.8f %.8f %.8f %.8f 1.0E-03 %.8f %.8f\n",RS0,dtadj,dtout,tcrit,rvir,mmean);
	} else {
		NNBMAX = sqrt(N); //Guess for the neighbour number (Nbody6 only)
		fprintf(PAR,"%i 1 10 199 %i 1\n",N,NNBMAX);
		fprintf(PAR,"0.02 0.03 %.8f %.8f %.8f %.8f 1.0E-03 %.8f %.8f\n",RS0,dtadj,dtout,tcrit,rvir,mmean);
	}
	fprintf(PAR,"2 2 1 0 1 0 0 0 0 2\n");
	fprintf(PAR,"0 0 0 %i 2 %i %i 0  %i 2\n",tf,regupdate,etaupdate,mloss);
	fprintf(PAR,"0 2 %i 0 1 2 -1 0 0 2\n",esc);
	fprintf(PAR,"0 0 0 0 1 0 0 2 -1 1\n");
	if (code == 1) fprintf(PAR,"0 0 0 0 0 0 0 0 0 0\n");
	if (code == 1) fprintf(PAR,"%.8f %.8f %.8f\n",M,mlow,mup);
	fprintf(PAR,"1.0E-5 1.0E-4 0.01 1.0 1.0E-06 0.01\n");
	
	if (code == 1) fprintf(PAR,"2.350000 %.8f %.8f %i %.8f 0.000000 100000.0\n",mlow,MMAX,(int) epoch,Z);
	else  fprintf(PAR,"2.350000 %.8f %.8f 0 0 %.8f %.8f %.8f\n",mlow,MMAX,Z,epoch,dtplot);
	
	fprintf(PAR,"0.5 0.0 0.0 0.00000\n");
	
	if (tf == 1) {
		if (code) {
			rtide = pow(1.0*M/(3.0*M1pointmass),1.0/3.0)*sqrt(RG[0]*RG[0]+RG[1]*RG[1]+RG[2]*RG[2]);
			fprintf(PAR,"%i %.8f\n",0,sqrt(1.0/(3.0*pow(rtide,3))));
		}
	} else if (tf == 2) {
		fprintf(PAR,"%.8e %.8f\n",M1pointmass,sqrt(RG[0]*RG[0]+RG[1]*RG[1]+RG[2]*RG[2]));
	} else if (tf == 3) {
		if (code == 2) fprintf(PAR,"%.6e %.6f %.6e %.6f %.6f %.6e %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",M1allen,b1allen,M2allen,a2allen,b2allen,M3allen,a3allen,RG[0]/1000.0,RG[1]/1000.0,RG[2]/1000.0,VG[0],VG[1],VG[2]);
		else fprintf(PAR,"%.6e %.6e %.6f %.6f %.6e %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n",M1allen,M2allen,a2allen,b2allen,VCIRC,RCIRC,RG[0]/1000.0,RG[1]/1000.0,RG[2]/1000.0,VG[0],VG[1],VG[2]);
	}
	if (code != 1) fprintf(PAR,"0.0 0.0 0.0 0.0\n");
	if (gpu) fprintf(PAR,"0.125\n");
	

	
	//write to .NBODY file
	for (j=0;j<N;j++) {
		fprintf(NBODY,"%.8f\t%.8f %.8f %.8f\t%.8f %.8f %.8f\n",kubus[j][0],kubus[j][1],kubus[j][2],kubus[j][3],kubus[j][4],kubus[j][5],kubus[j][6]);
	}

	fclose(PAR);
	fclose(NBODY);
	printf("\nData written to %s and %s\n", PARfile, NBODYfile);

	
	//energy check
	double ekin = 0.0;
	double epot = 0.0;
	
	if (check) {
		printf("\nMaking final energy check... (may take a while but can be aborted by pressing CTRL+c)\n");

		for (j=0; j<N; j++) {
			ekin += kubus[j][0]*((kubus[j][4]*kubus[j][4])+(kubus[j][5]*kubus[j][5])+(kubus[j][6]*kubus[j][6]));
			if (j) {
				for (i=0;i<j-1;i++) 
					epot -= kubus[i][0]*kubus[j][0]/sqrt((kubus[i][1]-kubus[j][1])*(kubus[i][1]-kubus[j][1])+(kubus[i][2]-kubus[j][2])*(kubus[i][2]-kubus[j][2])+(kubus[i][3]-kubus[j][3])*(kubus[i][3]-kubus[j][3])); 
			}
		} 
		ekin *= 0.5;
		
		printf("\nEkin = %g\t Epot = %g\t Etot = %g (Nbody units)\n", ekin, epot, ekin+epot);
	}
	
	
	
	for (j=0;j<NMAX;j++) free (kubus[j]);
	free(kubus);
	
	t2 = clock();														//stop stop-watch
	printf("\nElapsed time: %g sec\n",(double)(t2-t1)/CLOCKS_PER_SEC);	//print stopped time

	return 0;
}



/*************
 * Functions *
 *************/

int generate_m(int *N, double **kubus, double mlow, double mup, double *M, double *mmean, double MMAX, double Mcl) {
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
		*N = floor((Mcl-MMAX)/mth);
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
				kubus[i][0] = pow(0.5*c1*xx*k2+pow(mlow,c1),1.0/c1);
			else 
				kubus[i][0] = pow(c2*(xx*k2-k1)+pow(max(0.5,mlow),c2),1.0/c2);
		} while (kubus[i][0] > MMAX);

		if (kubus[i][0] > mostmassive) mostmassive = kubus[i][0];
		*M += kubus[i][0];
		if ((i==*N-1) && (*M<Mcl)) *N += 1;
	}

	printf("Total mass: %g\t(%i stars)\n",*M,*N);
	*mmean = *M/ *N;
	printf("Most massive star: %g\n",mostmassive);
	
	printf("Mean masses theoretical/data: %f %f\n",mth,*mmean);

	//convert masses to Nbody units
	for (i=0; i<*N; i++) {
		kubus[i][0] /= (*N**mmean);
	} 
	
	return 0;
}

int generate_plummer(int N, double **kubus, double rtide, double rvir){
	int i,j;
	double a[9], ri, sx, sv, rcut, r2, rvirscale;
	double ke = 0.0;
	double pe = 0.0;
	
	//Scale length
	sx = 1.5*TWOPI/16.0;
	sv = sqrt(1.0/sx);
	
	printf("Setting cut-off radius to approximate tidal radius\n");
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
		
			kubus[i][3] = (1.0 - 2.0*a[2])*ri;
			kubus[i][1] = sqrt(ri*ri-pow(kubus[i][3],2))*cos(TWOPI*a[3]);
			kubus[i][2] = sqrt(ri*ri-pow(kubus[i][3],2))*sin(TWOPI*a[3]); 
		} while (sqrt(pow(kubus[i][1],2)+pow(kubus[i][2],2)+pow(kubus[i][3],2))>rcut); //reject particles beyond tidal radius
		
		//velocities
		do {
			a[4] = drand48(); 
			a[5] = drand48(); 
			a[6] = pow(a[4],2)*pow(1.0 - pow(a[4],2),3.5);
		} while (0.1*a[5]>a[6]);
			
		a[8] = a[4]*sqrt(2.0)/pow(1.0 + ri*ri,0.25);
		a[6] = drand48(); 
		a[7] = drand48(); 
		
		kubus[i][6] = (1.0 - 2.0*a[6])*a[8];
		kubus[i][4] = sqrt(a[8]*a[8] - pow(kubus[i][6],2))*cos(TWOPI*a[7]);
		kubus[i][5] = sqrt(a[8]*a[8] - pow(kubus[i][6],2))*sin(TWOPI*a[7]);

		if (i) {
			for (j=0;j<i-1;j++) {
				r2 = (kubus[i][1]-kubus[j][1])*(kubus[i][1]-kubus[j][1]) + (kubus[i][2]-kubus[j][2])*(kubus[i][2]-kubus[j][2]) +(kubus[i][3]-kubus[j][3])*(kubus[i][3]-kubus[j][3]) ;
				if (r2<0.000001*rvir*rvir) printf("WARNING: %i %i\tdr = %lf\n", i,j,sqrt(r2));
				pe -=  kubus[i][0]*kubus[j][0]/sqrt(r2);
			}
		}
		
		ke += kubus[i][0]*(pow(kubus[i][4],2)+pow(kubus[i][5],2)+pow(kubus[i][6],2));
	}


	rvirscale = -GKING*pow(MTOTKING,2)/(2.0*pe);
	sx = 1.0/rvirscale;
	
	ke *= 0.5;
	sv = sqrt(4.0*ke);
		
	// Scale coordinates to N-body units
	for (i=0;i<N;i++) {
			kubus[i][1] *= sx;
			kubus[i][2] *= sx;
			kubus[i][3] *= sx;
			kubus[i][4] /= sv;
			kubus[i][5] /= sv;
			kubus[i][6] /= sv;
		}

		return 0;
}

int generate_king(int N, double W0, double **kubus, double *rvir, double *rh, double *rking){
	
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
	
	
	
	
	printf("\nGenerating King model with parameters: N = %i\t W0 = %g\n\n", N, W0);
	
	if (W0>12.0) {
		printf("W0 too large\n");
		return 0;
	} else if (W0 < 0.2) {
		printf("W0 too large\n");
		return 0;
	}
	
	
	
	
	/***************************
	 * INTERPOLATE KING VALUES *
	 ***************************/
	
	h = pow(W0-rmin,2)/(M-1);	//step size for interpolation
	den = densty(W0);			//central density
	
	//x is King's W, y[0] is z**2, y[1] is 2*z*dz/dW, where z is scaled radius
	//so x(y[0]) is energy as function of radius^2 and x(y[1]) is the derivative
	
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
	
	//Central velocity dispersion**2 (3-dimensional) *
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
		mstar = kubus[i][0];//MTOTKING/N;
		
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
				//r2 = pow(coord[j][0]-xstar,2) + pow(coord[j][1]-ystar,2) + pow(coord[j][2]-zstar,2);
				r2 = (coord[j][0]-xstar)*(coord[j][0]-xstar) + (coord[j][1]-ystar)*(coord[j][1]-ystar) + (coord[j][2]-zstar)*(coord[j][2]-zstar);
				if (r2<0.000001**rvir**rvir) printf("WARNING: %i %i\tdr = %lf\n", i,j,sqrt(r2));
				pe -=  mstar*kubus[j][0]/sqrt(r2);//1.0/sqrt(r2);
			}
		}
		ke += mstar*pow(speed,2);
		
	}
	
	pe = GKING*pe;//*pow(mstar,2);
	ke *= 0.5;
	*rvir = -GKING*pow(MTOTKING,2)/(2.0*pe);
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
			kubus[i][j+1] = coord[i][j]/ *rvir;
			kubus[i][j+4] = coord[i][j+3]/vscale;
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


