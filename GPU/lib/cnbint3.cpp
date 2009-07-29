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
	assert(nnb <= NBMAX);

	double xbuf[NBMAX] __attribute__ ((aligned(16)));
	double ybuf[NBMAX] __attribute__ ((aligned(16)));
	double zbuf[NBMAX] __attribute__ ((aligned(16)));
	float vxbuf[NBMAX] __attribute__ ((aligned(16)));
	float vybuf[NBMAX] __attribute__ ((aligned(16)));
	float vzbuf[NBMAX] __attribute__ ((aligned(16)));
	float mbuf[NBMAX] __attribute__ ((aligned(16)));
	for(int k=0; k<nnb; k++){
		int j = list[k];
		xbuf[k] = pos[j][0];
		ybuf[k] = pos[j][1];
		zbuf[k] = pos[j][2];
		vxbuf[k] = vel[j][0];
		vybuf[k] = vel[j][1];
		vzbuf[k] = vel[j][2];
		mbuf[k] = mass[j];
	}
	for(int k=nnb; k%4; k++){
		xbuf[k] = 16.0;
		ybuf[k] = 16.0;
		zbuf[k] = 16.0;
		vxbuf[k] = 0.0f;
		vybuf[k] = 0.0f;
		vzbuf[k] = 0.0f;
		mbuf[k] = 0.0f;
	}

	v4df xi(pos[i][0]);
	v4df yi(pos[i][1]);
	v4df zi(pos[i][2]);
	v4sf vxi(vel[i][0]);
	v4sf vyi(vel[i][1]);
	v4sf vzi(vel[i][2]);
	v4df ax(0.0), ay(0.0), az(0.0);
	v4sf jx(0.0), jy(0.0), jz(0.0);

	for(int k=0; k<nnb; k+=4){
		v4df xj(xbuf + k);
		v4df yj(ybuf + k);
		v4df zj(zbuf + k);
		v4sf vxj(vxbuf + k);
		v4sf vyj(vybuf + k);
		v4sf vzj(vzbuf + k);
		v4sf mj(mbuf + k);

		v4sf dx = v4sf::_v4sf(xj - xi);
		v4sf dy = v4sf::_v4sf(yj - yi);
		v4sf dz = v4sf::_v4sf(zj - zi);
		v4sf dvx = vxj - vxi;
		v4sf dvy = vyj - vyi;
		v4sf dvz = vzj - vzi;

		v4sf r2 = dx*dx + dy*dy + dz*dz;
		v4sf rv = dx*dvx + dy*dvy + dz*dvz;
		rv *= v4sf(-3.0f);
		// v4sf rinv1 = r2.rsqrt() & v4sf(mask);
		v4sf rinv1 = r2.rsqrt();
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
