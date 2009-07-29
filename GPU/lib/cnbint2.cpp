// #include <algorithm>
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
	typedef v4sf::_v4si v4si;
	v4df xi = pos[i][0];
	v4df yi = pos[i][1];
	v4df zi = pos[i][2];
	v4sf vxi = vel[i][0];
	v4sf vyi = vel[i][1];
	v4sf vzi = vel[i][2];
	v4df ax(0.0), ay(0.0), az(0.0);
	v4sf jx(0.0), jy(0.0), jz(0.0);
	for(int k=0; k<nnb; k+=4){
		v4si mask = {-1, -1, -1, -1};
		if(nnb-k < 4){
			if(nnb-k > 2)      mask = (v4si){-1, -1, -1,  0};
			else if(nnb-k > 1) mask = (v4si){-1, -1,  0,  0};
			else               mask = (v4si){-1,  0,  0,  0};
		}

		int j0 = list[k+0];
		int j1 = list[k+1];
		int j2 = list[k+2];
		int j3 = list[k+3];

		v4df xj (pos[j0][0], pos[j1][0], pos[j2][0], pos[j3][0]);
		v4df yj (pos[j0][1], pos[j1][1], pos[j2][1], pos[j3][1]);
		v4df zj (pos[j0][2], pos[j1][2], pos[j2][2], pos[j3][2]);
		v4sf vxj(vel[j0][0], vel[j1][0], vel[j2][0], vel[j3][0]);
		v4sf vyj(vel[j0][1], vel[j1][1], vel[j2][1], vel[j3][1]);
		v4sf vzj(vel[j0][2], vel[j1][2], vel[j2][2], vel[j3][2]);
		v4sf mj (mass[j0], mass[j1], mass[j2], mass[j3]);

		v4sf dx = v4sf::_v4sf(xj - xi);
		v4sf dy = v4sf::_v4sf(yj - yi);
		v4sf dz = v4sf::_v4sf(zj - zi);
		v4sf dvx = vxj - vxi;
		v4sf dvy = vyj - vyi;
		v4sf dvz = vzj - vzi;

		v4sf r2 = dx*dx + dy*dy + dz*dz;
		v4sf rv = dx*dvx + dy*dvy + dz*dvz;
		rv *= v4sf(-3.0f);
		v4sf rinv1 = r2.rsqrt() & v4sf(mask);
		v4sf rinv2 = rinv1 * rinv1;
		rv *= rinv2;
		// v4sf rinv3 = (mj * rinv1 * rinv2) & v4sf(mask);
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
