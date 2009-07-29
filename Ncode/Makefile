.KEEP_STATE:
#suggested by Keigo Nitadori 11/06/08.
#FFLAGS = -O3 -Wall -ffast-math
#FFLAGS += -O3 -fbounds-check
#FFLAGS = -O3 -pg
CFLAGS = -O3
FFLAGS = -O3
#FC = f77
#LD = f77
FC = gfortran
LD = gfortran -lstdc++ -fopenmp

VPATH = ../Chain:../Nchain

SOURCE = \
nbody6.f adjust.f assess.f bindat.f binev.f binout.f binpop.f \
block.f bodies.f brake.f brake2.f brake3.f bsetid.f chaos0.f chaos.f \
check.f checkl.f chrect.f clint.f cloud.f cloud0.f cmbody.f \
cmcorr.f cmfirr.f cmfreg.f coal.f comenv.f core.f corerd.f \
cputim.f data.f decide.f define.f deform.f degen.f delay.f \
dgcore.f dtchck.f eccmod.f ecirc.f edot.f efac2.f efac3.f \
expel.f energy.f escape.f events.f evolve.f expand.f fclose.f \
fcloud.f fcorr.f fdisk.f fhalo.f ficorr.f findj.f findm.f \
flyby.f fnuc.f fpcorr.f fpert.f fpoly1.f fpoly2.f freeze.f \
gcinit.f gcint.f giant.f giant3.f gntage.f grrad.f hcorr.f \
hiarch.f hicirc.f hidat.f higrow.f himax.f himax2.f himod.f \
hipop.f hirect.f histab.f hivel.f hmdot.f hmdot2.f hotsys.f \
hrdiag.f hrplot.f hut.f hut2.f iblock.f imf.f imfbd.f imf2.f \
impact.f induce.f inext.f input.f insert.f instar.f intgrt.f \
intide.f jacobi.f kepler.f kick.f kick2.f ksapo.f kscorr.f \
ksin2.f ksinit.f ksint.f kslist.f ksmod.f ksperi.f kspert.f \
kspoly.f kspred.f ksrect.f ksreg.f ksres.f ksres2.f ksterm.f \
kstide.f lagr.f lagr2.f levels.f magbrk.f matrix.f mdot.f merge.f \
merge2.f mix.f mloss.f mlwind.f modify.f mrenv.f mtrace.f mydump.f \
nbhist.f nbint.f nblist.f nbpot.f nbrem.f nbrest.f nbsort.f nbtide.f \
newtev.f nstab.f ntint.f offset.f orbit.f output.f peri.f permit.f \
pfac.f poti.f proto_star.f qtides.f ran2.f reflct.f regint.f \
remove.f rename.f reset.f reset2.f resolv.f rkint.f rl.f roche.f \
rpmax.f rpmax2.f rpmin.f scale.f search.f setup.f setup2.f shrink.f \
sort1.f spiral.f stability.f star.f start.f stepk.f steps.f stumpf.f \
subint.f swap.f sweep.f synch.f tail0.f tcirc.f tides.f tides2.f \
tides3.f touch.f tpert.f trdot.f trdot2.f trflow.f tstab.f tstep.f \
units.f unpert.f update.f verify.f xtrnl0.f xtrnld.f xtrnlf.f xtrnlp.f \
xtrnlt.f xtrnlv.f xvpred.f zare.f zcnsts.f zero.f zfuncs.f \
triple.f derqp3.f difsy3.f erel3.f extend.f qpmod3.f stabl3.f \
stablz.f start3.f subsys.f tperi.f trans3.f \
quad.f derqp4.f difsy4.f endreg.f erel4.f ichain.f newreg.f newsys.f \
qpmod4.f rchain.f rsort.f stabl4.f start4.f status.f trans4.f \
cfuncs.f chain.f chstab.f const.f cstab2.f cstab3.f cstab4.f cstab5.f \
derqp.f difsy1.f erel.f hpsort.f inclin.f invert.f ksphys.f physks.f \
qforce.f qpmod.f r2sort.f recoil.f redraw.f select.f slow.f stablc.f \
swcond.f switch.f transk.f transq.f transx.f vector.f xtf.f xtrnlu.f \
ycopy.f ysave.f \
absorb.f chaos2.f chdata.f chfind.f chfirr.f chinit.f chlist.f chmod.f \
chpert.f chpot.f chterm.f expel2.f fchain.f ghost.f giant2.f kcpert.f \
reduce.f reinit.f renew.f setsys.f tchain.f xcpred.f xtpert.f


OBJECTS = $(SOURCE:.f=.o)

nbody6:	$(OBJECTS)
	$(LD) $(FFLAGS) $(OBJECTS) -o nbody6

print:
	@- \rm -f NBODY6.TEXT
	@cat $(SOURCE) > NBODY6.TEXT

clean:
	rm *.o
