FFLAGS = -O2
#FFLAGS = -g -fnonstd
FORTRAN = f77

SOURCE = \
absorb.f chfirr.f chfind.f chinit.f chlist.f chmod.f \
chpot.f chterm.f fchain.f ghost.f kcpert.f reduce.f \
reinit.f renew.f setsys.f tchain.f xcpred.f xtpert.f

OBJECTS = $(SOURCE:.f=.o)

nchain:	$(OBJECTS)

printh:
	@- \rm -f NCHAIN.TEXT
	@cath $(SOURCE) > NCHAIN.TEXT
print:
	@- \rm -f NCHAIN.TEXT
	@cat $(SOURCE) > NCHAIN.TEXT
