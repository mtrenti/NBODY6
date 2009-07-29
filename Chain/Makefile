FFLAGS = -O
#FFLAGS = -g -fnonstd
FORTRAN = f77

SOURCE = \
cfuncs.f chain.f chstab.f const.f cstab2.f cstab3.f cstab4.f \
cstab5.f derqp.f difsy1.f erel.f hpsort.f inclin.f invert.f \
ksphys.f physks.f qforce.f qpmod.f r2sort.f recoil.f redraw.f \
select.f slow.f stablc.f swcond.f switch.f transk.f transq.f \
transx.f vector.f xtf.f xtrnlu.f ycopy.f ysave.f

OBJECTS = $(SOURCE:.f=.o)

chain:	$(OBJECTS)

print:
	@- \rm -f CHAIN.TEXT
	@cat $(SOURCE) > CHAIN.TEXT

printh:
	@- \rm -f CHAIN.TEXT
	@cath $(SOURCE) > CHAIN.TEXT
