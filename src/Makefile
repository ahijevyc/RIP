#        Lower-level Makefile for RIP

include ../configure.rip

OBJS= \
	addorfill.o adjellip.o angle.o barb.o bes.o bogs.o \
	brnshr.o bshear.o bvfricalc.o capecalc3d.o \
	ceilingcalc.o colram.o colrln.o condheat.o \
	cpchhl.o cpchll.o \
	cpcolr.o cpmpxy.o createdtg.o cttcalc.o dbzcalc.o ddpcalc.o \
	ddx.o ddy.o derivc.o \
	drawcl.o driver.o eqthecalc.o \
	extingcalc.o fields.o fillarray.o finterp.o \
	finterp2d.o fregcalc.o frtitle.o fy.o fx.o getbasicvars.o \
	getchar.o getconvals.o getdash.o getdims.o \
	getheadinfo.o getinum.o \
	getminfo.o getpt.o getrnum.o gettimes.o getvar.o getvarinfo.o \
	hboxdraw.o hbulldraw.o \
	hchadraw.o hcondraw.o hiresmap.o hlinedraw.o \
	hmapdraw.o hodograph.o hsidsdraw.o hstrdraw.o \
	hticdraw.o htrajdraw.o hvecdraw.o igetcoind.o lennonblank.o \
	liftparcel.o locinterp.o maptform.o maptick.o \
	mconvert.o mktrjpts.o mlen.o monotonic.o netasc.o \
	numeric.o omgcalc.o pfcalc.o pltitle.o prcver.o precprob.o \
	premaptform.o profil.o pvocalc.o qgomg.o rdcolt.o \
	readspec.o relhl.o rhucalc.o \
	rip.o rotipslab.o rotpslab.o sateqthecalc.o saweli.o scondraw.o \
        seaprs.o setripheader.o smooth.o smoothcp.o smoothcz.o srflow.o \
	srflow2.o sstats.o sticdraw.o sticdraw_polar.o \
        svecdraw.o tdpcalc.o \
	tfpcalc.o thecalc.o tonpsadiabat.o \
	tserprep.o turb.o unblank.o vbardraw.o vc2dcalc.o vcondraw.o \
	velbrb.o velvctmts.o vgp.o vinterp.o \
	virtual.o virtualhyd.o viscalc.o vticdraw.o vtrajdraw.o \
	vvecdraw.o vwinddraw.o wdircalc.o wetbulbcalc.o \
	writefile.o wshear.o wspcalc.o xmapcalc.o xtodot.o binio.o v5d.o typhoon.o

# Lines from here on down should not need to be changed.  They are the
# actual rules which make uses to build a.out.

all: 	rip ripcomp ripdp_mm5 ripdp_wrf ripinterp ripshow showtraj tabdiag upscale

rip: 	$(OBJS)
	$(CF) -o $@ $(LDFLAGS) $(OBJS) $(LIBS) $(LOCAL_LIBS)

ripcomp:
	$(CF) $(FFLAGS2) -o $@ $@.f $(LOCAL_LIBS)
ripdp_mm5:	ripdp_mm5.f	fillarray.o	mconvert.o	writefile_rdp.o	\
		virtual.o
	$(CF) $(FFLAGS2) -o $@ $@.f fillarray.o mconvert.o writefile_rdp.o \
		virtual.o $(LOCAL_LIBS)
ripdp_wrf:	ripdp_wrf.f	premaptform.o	maptform.o	fillarray.o \
		mconvert.o	writefile_rdp.o	virtual.o	xtodot.o
	$(CF) $(FFLAGS3) -o $@  -L$(NETCDFLIB)  -I$(NETCDFINC) $@.f \
		premaptform.o maptform.o fillarray.o mconvert.o writefile_rdp.o \
		virtual.o xtodot.o -lnetcdf $(LOCAL_LIBS)
ripinterp:
	$(CF) $(FFLAGS2) -o $@ $@.f $(LOCAL_LIBS)
ripshow:
	$(CF) $(FFLAGS2) -o $@ $@.f $(LOCAL_LIBS)
showtraj:
	$(CF) $(FFLAGS2) -o $@ $@.f $(LOCAL_LIBS)
tabdiag:	tabdiag.f	lennonblank.f
	$(CF) $(FFLAGS2) -o $@ $@.f lennonblank.f $(LOCAL_LIBS)
upscale:
	$(CF) $(FFLAGS2) -o $@ $@.f $(LOCAL_LIBS)

addorfill.o:    comconst
adjellip.o:     comconst
binio.o:        binio.h
bogs.o:         comconst
brnshr.o:       comconst
bvfricalc.o:    comconst
capecalc3d.o:   comconst
ceilingcalc.o:  comconst
condheat.o:     comconst
cpmpxy.o:       comconst    comvctran
cttcalc.o:      comconst
dbzcalc.o:      comconst
ddpcalc.o:      comconst
ddx.o:          comconst
ddy.o:          comconst
derivc.o:       comconst
driver.o:       comconst    pointers    v5df.h      comvctran	driver.f
	$(CF) -c $(FFLAGS2) $*.f
eqthecalc.o:    comconst
extingcalc.o:   comconst
fields.o:       comconst    pointers	fields.f
	$(CF) -c $(FFLAGS2) $*.f
fregcalc.o:     comconst
fy.o:           comconst    comvctran
fx.o:           comconst    comvctran
getchar.o:      comconst
getdash.o:      comconst
getinum.o:      comconst
getpt.o:        comconst
getrnum.o:      comconst
getvar.o:      comconst
hboxdraw.o:     comconst
hbulldraw.o:    comconst
hchadraw.o:     comconst
hcondraw.o:     comconst	hcondraw.f
	$(CF) -c $(FFLAGS2) $*.f
hiresmap.o:     comconst
hlinedraw.o:    comconst
hmapdraw.o:     comconst
hsidsdraw.o:    comconst
hstrdraw.o:     comconst
hticdraw.o:     comconst
htrajdraw.o:    comconst
hvecdraw.o:     comconst
igetcoind.o:    comconst
lennonblank.o:  comconst
liftparcel.o:   comconst
locinterp.o:    comconst
netasc.o:       comconst
omgcalc.o:      comconst
pfcalc.o:       comconst
pltitle.o:      comconst
prcver.o:       comconst
precprob.o:     comconst
profil.o:       comconst
pvocalc.o:      comconst
qgomg.o:        comconst	qgomg.f
	$(CF) -c $(FFLAGS2) $*.f
rdcolt.o:       comconst
readspec.o:     comconst	readspec.f
	$(CF) -c $(FFLAGS2) $*.f
rhucalc.o:      comconst
rip.o:          comconst
sateqthecalc.o:             comconst
saweli.o:       comconst
scondraw.o:     comconst
seaprs.o:
setripheader.o:	comconst
smooth.o:       comconst
smoothcp.o:     comconst
smoothcz.o:     comconst
sstats.o:       comconst	sstats.f
	$(CF) -c $(FFLAGS2) $*.f
sticdraw.o:     comconst
sticdraw_polar.o:           comconst
svecdraw.o:     comconst
tdpcalc.o:      comconst
tfpcalc.o:      comconst
thecalc.o:      comconst
tonpsadiabat.o: comconst
tserprep.o:     comconst
turb.o:         comconst
unblank.o:      comconst
v5d.o:          binio.h     v5d.h       vis5d.h
vbardraw.o:	comconst
vc2dcalc.o:     comconst    comvctran
vcondraw.o:     comconst	vcondraw.f
	$(CF) -c $(FFLAGS2) $*.f
velbrb.o:       comconst
vinterp.o:      comconst
virtual.o:
virtualhyd.o:
viscalc.o:      comconst
vticdraw.o:     comconst
vtrajdraw.o:    comconst
vvecdraw.o:     comconst    comvctran
vwinddraw.o:    comconst    comvctran
wdircalc.o:     comconst
wetbulbcalc.o:  comconst
writefile.o:    comconst
wspcalc.o:      comconst
typhoon.o:      comconst
