#include <iostream>
#include <cmath>
#include <cassert>

#define TMAX 8
#if 1
#include <omp.h>
#else
static inline int omp_get_num_threads(){return 1;}
static inline int omp_get_thread_num() {return 0;}
#endif

struct Jparticle{
	float x[3];
	float m;
	float v[3];
	float pad;
	Jparticle() {}
	Jparticle(double mj, double xj[3], double vj[3]){
		x[0] = xj[0];
		x[1] = xj[1];
		x[2] = xj[2];
		m    = mj;
		v[0] = vj[0];
		v[1] = vj[1];
		v[2] = vj[2];
	}
};
static Jparticle *jp_host;
// static int *nblist;
static int *nblistbuf[TMAX];
static int nbody, nbodymax;

void GPUNB_open(int nbmax){
	// std::cout << "Open GPUNB " << nbmax << std::endl;
	nbodymax = nbmax;
	jp_host = new Jparticle[nbmax];
#pragma omp parallel
	{
		int nth = omp_get_num_threads();
		assert(nth <= TMAX);
		int tid = omp_get_thread_num();
		nblistbuf[tid] = new int[nbmax];
	}
}

void GPUNB_close(){
	// std::cout << "Close GPUNB" << std::endl;
	delete [] jp_host;
#pragma omp parallel
	{
		int tid = omp_get_thread_num();
		delete [] nblistbuf[tid];
	}
}

void GPUNB_send(
		int nj,
		double mj[],
		double xj[][3],
		double vj[][3]){
	nbody = nj;
	// std::cout << "gpu send: " << nbody << " " << nbodymax << std::endl;
	assert(nbody <= nbodymax);
#pragma omp parallel for
	for(int j=0; j<nj; j++){
		jp_host[j] = Jparticle(mj[j], xj[j], vj[j]);
	}
}

void GPUNB_regf(
		int ni,
		double h2d[],
		double xid[][3],
		double vid[][3],
		double acc[][3],
		double jrk[][3],
		int lmax,
		int nbmax,
		int *listbase){
	// std::cout << " Call GPUNB_regf " << ni << std::endl;
#pragma omp parallel for
	for(int i=0; i<ni; i++){
		int tid = omp_get_thread_num();
		// std::cout << tid << " : " << i << std::endl;
		int *nblist = nblistbuf[tid];
		float ax=0, ay=0, az=0;
		float jx=0, jy=0, jz=0;
		float xi[3] = {xid[i][0], xid[i][1], xid[i][2]};
		float vi[3] = {vid[i][0], vid[i][1], vid[i][2]};
		float h2 = h2d[i];
		int nnb = 0;
		for(int j=0; j<nbody; j++){
			Jparticle &jp = jp_host[j];
			float dx = jp.x[0] - xi[0];
			float dy = jp.x[1] - xi[1];
			float dz = jp.x[2] - xi[2];
			float r2 = dx*dx + dy*dy + dz*dz;
			if(r2 < h2){
				nblist[nnb++] = j;
				continue;
			}
			float dvx = jp.v[0] - vi[0];
			float dvy = jp.v[1] - vi[1];
			float dvz = jp.v[2] - vi[2];
			float rv = dx*dvx + dy*dvy + dz*dvz;
			float rinv2 = 1.f/r2;
			float rinv1 = sqrtf(rinv2);
			rv *= -3.f * rinv2;
			float rinv3 = jp.m * rinv1 * rinv2;
			ax += rinv3 * dx;
			ay += rinv3 * dy;
			az += rinv3 * dz;
			jx += rinv3 * (dvx + rv * dx);
			jy += rinv3 * (dvy + rv * dy);
			jz += rinv3 * (dvz + rv * dz);
		}
		acc[i][0] = ax;
		acc[i][1] = ay;
		acc[i][2] = az;
		jrk[i][0] = jx;
		jrk[i][1] = jy;
		jrk[i][2] = jz;
		// store nblist in FORTRAN style
		int *nnbp = listbase + lmax * i;
		int *nblistp = nnbp + 1;
		// assert(nnb <= nbmax);
		if(nnb > nbmax){
			*nnbp = -1;
		}else{
			*nnbp = nnb;
			for(int k=0; k<nnb; k++){
				nblistp[k] = nblist[k];
			}
		}
	}
	// printf("gpu: %e %e %e %d\n", xid[0][0], acc[0][0], jrk[0][0], *listbase);
#if 1
	if(ni > 0){
		FILE *fp = fopen("Force.hst", "w");
		assert(fp);
		for(int i=0; i<ni; i++){
			int nnb =  listbase[i*lmax];
			fprintf(fp, "%d %9.2e %9.2e %9.2e %9.2e %9.2e %9.2e %d\n",
					i, acc[i][0], acc[i][1], acc[i][2], 
					   jrk[i][0], jrk[i][1], jrk[i][2], nnb);
		}
		fprintf(fp, "\n");
		fclose(fp);
		exit(1);
	}
#endif
}

extern "C" {
	void gpunb_open_( int *nbmax){
		GPUNB_open(*nbmax);
	}
	void gpunb_close_(){
		GPUNB_close();
	}
	void gpunb_send_(
			int *nj,
			double mj[],
			double xj[][3],
			double vj[][3]){
		GPUNB_send(*nj, mj, xj, vj);
	}
	void gpunb_regf_(
			int *ni,
			double h2[],
			double xi[][3],
			double vi[][3],
			double acc[][3],
			double jrk[][3],
			int *lmax,
			int *nbmax,
			int *list){ // list[][lmax]
		GPUNB_regf(*ni, h2, xi, vi, acc, jrk, *lmax, *nbmax, list);
	}
}
