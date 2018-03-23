#     Top-level Makefile for RIP
#
#
# This is the makefile for program RIP. RIP requires NCAR Graphics to run.
#
# Before you compile RIP, make sure that you either type
#
#    setenv NCARG_ROOT /usr/local/ncarg
#       (note: location of NCAR Graphics library may vary on your machine)
#
# or include the above line in your .cshrc file.
#
# To compile, edit the makefile, and uncomment the compile options
#    appropriate for your machine. Then type 'make' to compile.
#
# Other available make commands:
#
# To remove all the objects but leave the executables use the command line:
#
#       make clean
#
# To remove everything but the source files use the command line:
#
#       make clobber
#

#     Macros, these should be generic for all machines

.IGNORE:

AR 	= ar ru
RM 	= /bin/rm -f
MV 	= mv -f
LN	= ln -s
CC 	= cc
MAKE 	= make -i -f Makefile
RM_LIST	= src/*.o
RM_LIST2= rip ripcomp ripdp_mm5 ripdp_wrf ripinterp ripshow showtraj tabdiag upscale

default:
	@echo " "
	@echo "Type one oF the following:"
	@echo "	make dec 		For Run on DEC_ALPHA"
	@echo "	make linux		For Run on LINUX"
	@echo "	make sun		For Run on SUN"
	@echo "	make sun2		For Run on SUN if make sun didn't work"
	@echo "	make sun90		For Run on SUN using F90"
	@echo "	make sgi		For Run on SGI"
	@echo "	make sgi64		For Run on 64-bit SGI"
	@echo "	make ibm		For Run on IBM SP2"
	@echo "	make cray		For Run on NCAR's Cray"
	@echo "	make vpp300		For Run on Fujitsu VPP 300"
	@echo "	make vpp5000		For Run on Fujitsu VPP 5000"
	@echo "	make clean		to remove object files"
	@echo "	make clobber		to remove object files and executables"
	@echo " "
#
# Tunable parameters
#
# CF            Name of the fortran compiling system to use
# LDFLAGS       Flags to the loader
# LIBS          List of libraries
# CMD           Name of the executable
# FFLAGS        Flags to the compiler - including optimization
# FFLAGS2       Simpler compiler flags for routines that bog down the compiler
# CCFLAGS       Flags for c compiler (used for generating vis5d file)

dec:
	(cd src/ ; $(MAKE) all \
	"CF      = f90" \
        "FFLAGS  = -fast -convert big_endian -fpe2 " \
        "FFLAGS2 = -convert big_endian -fpe2 " \
	"CCFLAGS = -DLITTLE -DUNDERSCORE -c" \
        "LDFLAGS = " \
	"LIBS    = -L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lX11 -lm" )
	( $(RM) $(RM_LIST2) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrf . ; $(LN) src/ripcomp . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

#   Linux
#      Assume PGI is installed in /usr/pgi

linux:
	(cd src/ ; $(MAKE) all \
	"CF      = pgf90" \
	"FFLAGS  = -byteswapio " \
	"FFLAGS2 = -byteswapio " \
	"CCFLAGS = -DLITTLE -DUNDERSCORE -c" \
        "LDFLAGS = " \
	"LIBS    = -L$(NCARG_ROOT)/lib -lncarg -lcgm -lncarg_gks -lncarg_c -L/usr/X11R6/lib -lX11 -L$(PGI)/linux86/lib -lpgftnrtl -lpgc -L/usr/lib -lf2c" )
	( $(RM) $(RM_LIST2) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrf . ; $(LN) src/ripcomp . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

sun:
	(cd src/ ; $(MAKE) all \
	"CF      = f77" \
	"FFLAGS  = -fast -xarch=v8plus -xvector=yes -xprefetch=yes" \
	"FFLAGS2 = " \
	"CCFLAGS = -DUNDERSCORE -c" \
        "LDFLAGS = -fast -xarch=v8plus -xvector=yes -xprefetch=yes" \
	"LIBS    = -L$(NCARG_ROOT)/lib -R$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lncarg_loc -lX11 -lm -lM77 -lF77 -lsunmath" )
	( $(RM) $(RM_LIST2) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrf . ; $(LN) src/ripcomp . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

sun2:
	(cd src/ ; $(MAKE) all \
	"CF      = f77" \
	"FFLAGS  = -fast -xarch=v8plus -xvector=yes -xprefetch=yes" \
	"FFLAGS2 = " \
	"CCFLAGS = -DUNDERSCORE -c" \
        "LDFLAGS = -fast -xarch=v8plus -xvector=yes -xprefetch=yes" \
	"LIBS    = -L$(NCARG_ROOT)/lib -L/usr/openwin/lib -L/usr/dt/lib -lncarg -lncarg_gks -lncarg_c -lX11 -lm -lM77 -lF77 -lsunmath" )
	( $(RM) $(RM_LIST2) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrf . ; $(LN) src/ripcomp . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

sun90:
	(cd src/ ; $(MAKE) all \
	"CF      = f90" \
	"FFLAGS  = " \
	"FFLAGS2 = " \
	"CCFLAGS = -DUNDERSCORE -c" \
        "LDFLAGS = " \
	"LIBS    = -L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lX11 -lm -lsocket -lM77 -lF77" )
	( $(RM) $(RM_LIST2) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrf . ; $(LN) src/ripcomp . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

sgi:
	(cd src/ ; $(MAKE) all \
	"CF      = f77" \
	"FFLAGS  = -n32 -O3" \
	"FFLAGS2 = " \
	"CCFLAGS = -DUNDERSCORE -c" \
        "LDFLAGS = " \
	"LIBS    = -L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lX11" )
	( $(RM) $(RM_LIST2) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrf . ; $(LN) src/ripcomp . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

sgi64:
	(cd src/ ; $(MAKE) all \
	"CF      = f77" \
	"FFLAGS  = -64 -O3" \
	"FFLAGS2 = -64 -O1" \
	"CCFLAGS = -DUNDERSCORE -c" \
        "LDFLAGS = " \
	"LIBS    = -L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lX11" )
	( $(RM) $(RM_LIST2) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrf . ; $(LN) src/ripcomp . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

ibm:
	(cd src/ ; $(MAKE) all \
	"CF      = xlf" \
	"FFLAGS  = " \
	"FFLAGS2 = " \
	"CCFLAGS = -DLITTLE" \
        "LDFLAGS = -O3 -qarch=auto -qmaxmem=-1 -qblankpad" \
	"LIBS    = -L$(NCARG_LIB) -lncarg -lncarg_gks -lncarg_c -lX11 -lm" )
	( $(RM) $(RM_LIST2) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrf . ; $(LN) src/ripcomp . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

cray:
	(cd src/ ; $(MAKE) all \
	"CF      = ncargf77" \
	"FFLAGS  = " \
	"FFLAGS2 = " \
	"CCFLAGS = -D_CRAY" \
        "LDFLAGS = " \
	"LIBS    = " )
	( $(RM) $(RM_LIST2) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrf . ; $(LN) src/ripcomp . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

vpp300:
	(cd src/ ; $(MAKE) all \
	"CF      = frt" \
	"FFLAGS  = -sc -Kfreealloc -Karraystack3 -X9 -Sw -Wv,-Of,-te,-ilfunc,-noalias,-m3,-P255 -Oe,-Kfast" \
	"FFLAGS2 = " \
	"CCFLAGS = -DUNDERSCORE" \
        "LDFLAGS = -J" \
	"LIBS    = -L/usr/X11R6/lib -Wl,-dy /usr/X11R6/lib/libX11.so.6 -lsocket -lnsl -L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lm /nwp/nwptest/lib/libnwp.a /package/dms/librdms.a /package/dms/libgdbm.a" )
	( $(RM) $(RM_LIST2) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrf . ; $(LN) src/ripcomp . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

vpp5000:
	(cd src/ ; $(MAKE) all \
	"CF      = frt" \
	"FFLAGS  = -sc -Kfreealloc -Karraystack3 -X9 -Sw -Wv,-Of,-te,-ilfunc,-noalias,-m3,-P255 -Oe -KA32" \
	"FFLAGS2 = -KA32" \
	"CCFLAGS = -KA32 -DUNDERSCORE" \
        "LDFLAGS = -J -KA32" \
	"LIBS    = -L/usr/lib32EX -Wl,-dy -lX11 -lSM -lICE -lgen -ldl -lsocket -lw -lnsl -L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lm" )
	( $(RM) $(RM_LIST2) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrf . ; $(LN) src/ripcomp . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

clean:
	$(RM) $(RM_LIST)

clobber:	clean
	$(RM) $(RM_LIST2) ; cd src/ ; $(RM) $(RM_LIST2)
