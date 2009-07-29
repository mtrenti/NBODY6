FFLAGS = -O3 -Wall
#FFLAGS = -O3 -fbounds-check
#FFLAGS = -O1
#FFLAGS += -O3 -pg
#FFLAGS = -g
FC = gfortran
CXX = g++ -O3 -Wall -march=core2 #-pg
#CXX = g++ -O3 -Wall -msse3
CUDA_PATH = /usr/local/cuda
SDK_PATH = /usr/local/NVIDIA_CUDA_SDK
NVCC = $(CUDA_PATH)/bin/nvcc --host-compilation C++ -Xcompiler "-Wall -O3 -fopenmp" -I $(SDK_PATH)/common/inc/ #-pg
RUNDIR = ./run

OBJECTS = $(SOURCE:.f=.o)

all: gpu sse

gpu: intgrt.o gpucor.o nbintp.o cnbint.o nbint.o start.o fpoly2a.o adjust.o energy2.o phicor.o cmfirr.o cmfreg.o cmfreg2.o regint2.o kspert.o swap.o scale.o wtime.o
	$(NVCC) -c ./lib/gpunb.multi.cu -o gpunb.o
	$(NVCC) -c ./lib/gpupot.gpu.cu -o gpupot.o
	cp -f *.o ../Ncode
	make gpu -C ../Ncode -f Makefile_gpu -j4
	mv -f ../Ncode/nbody6 $(RUNDIR)/nbody6.gpu

sse: intgrt.o gpucor.o nbintp.o cnbint.o nbint.o adjust.o energy2.o phicor.o cmfirr.o cmfreg.o cmfreg2.o regint2.o kspert.o swap.o start.o fpoly2a.o scale.o wtime.o
	$(CXX) -c ./lib/gpunb.sse.cpp -o gpunb.o -fopenmp
	$(CXX) -c ./lib/gpupot.sse.cpp -o gpupot.o -fopenmp
	cp -f *.o ../Ncode
	make sse -C ../Ncode -f Makefile_gpu -j4
	mv -f ../Ncode/nbody6 $(RUNDIR)/nbody6.sse

cnbint.o: ./lib/cnbint.cpp
	$(CXX) -c ./lib/cnbint.cpp

wtime.o: ./lib/wtime.cpp
	$(CXX) -c ./lib/wtime.cpp

intgrt.o: intgrt.omp.f
	$(FC) $(FFLAGS) intgrt.omp.f -c -o intgrt.o -fopenmp

start.o: start.omp.F
	$(FC) $(FFLAGS) start.omp.F -c -o start.o -fopenmp

clean:
	rm -f *.o
