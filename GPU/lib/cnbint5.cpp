#include <algorithm>
#include <cmath>
#include <cassert>
#include "v4df.h"
#include "v4sf.h"

void cnbint(
		int i,
		const double pos[][3],
		const double vel[][3],
		const double mass[],
		int nnb,
		int list[],
		double f[3],
		double fdot[3]){
	const int NBMAX = 512;
	// assert(nnb <= NBMAX);

	float xbuf[NBMAX] __attribute__ ((aligned(16)));
	float ybuf[NBMAX] __attribute__ ((aligned(16)));
	float zbuf[NBMAX] __attribute__ ((aligned(16)));
	float vxbuf[NBMAX] __attribute__ ((aligned(16)));
	float vybuf[NBMAX] __attribute__ ((aligned(16)));
	float vzbuf[NBMAX] __attribute__ ((aligned(16)));
	float mbuf[NBMAX] __attribute__ ((aligned(16)));
	assert((unsigned long)xbuf % 16 == 0);

	double xi = pos[i][0];
	double yi = pos[i][1];
	double zi = pos[i][2];
	float vxi = vel[i][0];
	float vyi = vel[i][1];
	float vzi = vel[i][2];
	v4df ax(0.0), ay(0.0), az(0.0);
	v4sf jx(0.0), jy(0.0), jz(0.0);
	for(int koff=0; koff<nnb; koff+=NBMAX){
		int nk = std::min(nnb-koff, NBMAX);
		for(int k=0; k<nk; k++){
			int j = list[k+koff];
			assert(j != i);
#if 1
			int jj = list[k+koff+4];
			__builtin_prefetch(pos[jj]);
			__builtin_prefetch(pos[jj]+2);
			__builtin_prefetch(vel[jj]);
			__builtin_prefetch(vel[jj]+2);
			__builtin_prefetch(&mass[jj]);
#endif
			double xj = pos[j][0];
			double yj = pos[j][1];
			double zj = pos[j][2];
			float vxj = vel[j][0];
			float vyj = vel[j][1];
			float vzj = vel[j][2];
			float mj = mass[j];
			xj -= xi;
			yj -= yi;
			zj -= zi;
			vxj -= vxi;
			vyj -= vyi;
			vzj -= vzi;
			xbuf[k] = xj;
			ybuf[k] = yj;
			zbuf[k] = zj;
			vxbuf[k] = vxj;
			vybuf[k] = vyj;
			vzbuf[k] = vzj;
			mbuf[k] = mj;
		}
		for(int k=nk; k%4; k++){
			xbuf[k] = 16.0f;
			ybuf[k] = 16.0f;
			zbuf[k] = 16.0f;
			vxbuf[k] = 0.0f;
			vybuf[k] = 0.0f;
			vzbuf[k] = 0.0f;
			mbuf[k] = 0.0f;
		}


		for(int k=0; k<nk; k+=4){
			v4sf dx(xbuf + k);
			v4sf dy(ybuf + k);
			v4sf dz(zbuf + k);
			v4sf dvx(vxbuf + k);
			v4sf dvy(vybuf + k);
			v4sf dvz(vzbuf + k);
			v4sf mj(mbuf + k);

			v4sf r2 = dx*dx + dy*dy + dz*dz;
			v4sf rv = dx*dvx + dy*dvy + dz*dvz;
			rv *= v4sf(-3.0f);
			// v4sf rinv1 = r2.rsqrt() & v4sf(mask);
#if 1
			v4sf rinv1 = r2.rsqrt();
#else
			v4sf rinv1 = v4sf(v4df(1.0) / v4df(r2).sqrt());
#endif
			v4sf rinv2 = rinv1 * rinv1;
			rv *= rinv2;
			v4sf rinv3 = mj * rinv1 * rinv2;
			 
			dx *= rinv3; ax += v4df(dx);
			dy *= rinv3; ay += v4df(dy);
			dz *= rinv3; az += v4df(dz);
			dvx *= rinv3; jx += dvx;
			dvy *= rinv3; jy += dvy;
			dvz *= rinv3; jz += dvz;
			dx *= rv; jx += dx;
			dy *= rv; jy += dy;
			dz *= rv; jz += dz;
		}
	}
	f[0] = ax.sum();
	f[1] = ay.sum();
	f[2] = az.sum();
	fdot[0] = (v4df(jx)).sum();
	fdot[1] = (v4df(jy)).sum();
	fdot[2] = (v4df(jz)).sum();
	assert(f[0] == f[0]);
	assert(f[1] == f[1]);
	assert(f[2] == f[2]);
	assert(fdot[0] == fdot[0]);
	assert(fdot[1] == fdot[1]);
	assert(fdot[2] == fdot[2]);
}

extern "C" {
	void cnbint_(
		int *i,
		double pos[][3],
		double vel[][3],
		double mass[],
		int *nnb,
		int *nblist,
		double f[3],
		double fdot[3]){
		cnbint(*i, pos-1, vel-1, mass-1, *nnb-1, nblist, f, fdot);
	}
}
