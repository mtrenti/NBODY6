#include <iostream>
#include <sys/time.h>
static double get_wtime(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	return tv.tv_sec + 1.e-6 * tv.tv_usec;
}

static double t0 = get_wtime();

extern "C" {
	void wtime_(int *ihour, int *imin, int *isec){
		double t = get_wtime() - t0;
		int ih = t / 3600.;
		t -= ih * 3600.;
		int im = t / 60.;
		t -= im * 60.;
		int is = t;
		*ihour = ih;
		*imin = im;
		*isec = is;
	}
}

#if 0
#include <unistd.h>
int main(){
	int ih, im, is;
	sleep(1);
	wtime_(&ih, &im, &is);
	std::cout << ih << ":" << im << ":" << is << std::endl;
	return 0;
}
#endif
