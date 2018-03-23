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
# To compile, edit the makefile to set the NETCDFLIB and NETCDFINC variables
# for your machine's environment.  Type 'make' to see a list of supported
# architectures.  Then type 'make architecture' where architecture is one of 
# the supported architectures.
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
RM_LIST_OBJ	= *.o
RM_LIST_EXE	= rip ripcomp ripdp_mm5 ripdp_wrfnmm ripdp_wrfarw ripcut ripinterp ripshow showtraj tabdiag upscale
#

default:
	@echo " "
	@echo "Type one of the following:"
	@echo "	make dec 		to run on DEC_ALPHA"
	@echo "	make linux		to run on LINUX with PGI compiler"
	@echo "	make linuxuw		to run on LINUX (at U.Wash.) with PGI compiler"
	@echo "	make intel		to run on LINUX with INTEL compiler"
	@echo "	make mac_xlf		to run on MAC_OS_X with Xlf Compiler"
	@echo "	make mac		to run on MAC_OS_X with Absoft Compiler"
	@echo "	make sun		to run on SUN"
	@echo "	make sun2		to run on SUN if make sun didn't work"
	@echo "	make sun90		to run on SUN using F90"
	@echo "	make sgi		to run on SGI"
	@echo "	make sgi64		to run on 64-bit SGI"
	@echo "	make ibm		to run on IBM SP2"
	@echo "	make cray		to run on NCAR's Cray"
	@echo "	make vpp300		to run on Fujitsu VPP 300"
	@echo "	make vpp5000		to run on Fujitsu VPP 5000"
	@echo "	make clean		to remove object files"
	@echo "	make clobber		to remove object files and executables"
	@echo " "
#
# Tunable parameters
#
# CF            Name of the fortran compiling system to use
# LDFLAGS       Flags to the loader
# LOCAL_LIBS    List of local libraries for all executables
# LIBS          List of libraries for RIP
# CMD           Name of the executable
# FFLAGS        Flags to the compiler - including optimization
# FFLAGS2       Simpler compiler flags for routines that bog down the compiler
# FFLAGS3       Compiler flags for ripdp_wrfXXX (which depends on netcdf)
# CCFLAGS       Flags for c compiler (used for generating vis5d file)

dec:
	(cd src/ ; $(MAKE) all \
	"CF      = f90" \
        "FFLAGS  = -fast -convert big_endian -fpe2 " \
        "FFLAGS2 = -convert big_endian -fpe2 " \
        "FFLAGS3 = -fast -convert big_endian -fpe2 " \
	"CCFLAGS = -DLITTLE -DUNDERSCORE -c" \
        "LDFLAGS = " \
	"LOCAL_LIBS = " \
	"NETCDFLIB = /usr/local/netcdf/lib" \
	"NETCDFINC = /usr/local/netcdf/include" \
	"LIBS    = -L$(NCARG_ROOT)/lib -L$(NETCDFLIB) -I$(NETCDFINC) -lncarg -lncarg_gks -lncarg_c -lX11 -lm" )
	( $(RM) $(RM_LIST_EXE) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrfarw . ; $(LN) src/ripdp_wrfnmm . ; $(LN) src/ripcomp . ; $(LN) src/ripcut . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

#   Linux
#      Assume PGI is installed in /usr/pgi

linux:
	(cd src/ ; $(MAKE) all \
	"CF      = pgf90" \
	"FFLAGS  = -byteswapio " \
	"FFLAGS2 = -byteswapio " \
	"FFLAGS3 = -byteswapio " \
	"CCFLAGS = -DLITTLE -DUNDERSCORE -c" \
        "LDFLAGS = " \
	"LOCAL_LIBS = " \
	"NETCDFLIB = /usr/local/netcdf/lib" \
	"NETCDFINC = /usr/local/netcdf/include" \
	"LIBS    = -L$(NCARG_ROOT)/lib -lncarg -lcgm -lncarg_gks -lncarg_c -L/usr/X11R6/lib -lX11 -L$(PGI)/linux86/lib -lpgftnrtl -lpgc -L/usr/lib -lf2c" )
	( $(RM) $(RM_LIST_EXE) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrfarw . ; $(LN) src/ripdp_wrfnmm . ; $(LN) src/ripcomp . ; $(LN) src/ripcut . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

linuxuw:
	(cd src/ ; $(MAKE) all \
	"CF      = pgf90" \
	"FFLAGS  = -byteswapio " \
	"FFLAGS2 = -byteswapio " \
	"FFLAGS3 = -byteswapio " \
	"CCFLAGS = -DLITTLE -DUNDERSCORE -c" \
        "LDFLAGS = " \
	"LOCAL_LIBS = " \
	"NETCDFLIB = /usr/lib" \
	"NETCDFINC = /usr/include" \
	"LIBS    = -L$(NCARG_ROOT)/lib -L$(NETCDFLIB) -I$(NETCDFINC) -lncarg -lcgm -lncarg_gks -lncarg_c -L/usr/X11R6/lib -lX11 -L$(PGI)/linux86/lib -lpgftnrtl -lpgc -L/usr/lib -lf2c -lnetcdf" )
	( $(RM) $(RM_LIST_EXE) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrfarw . ; $(LN) src/ripdp_wrfnmm . ; $(LN) src/ripcomp . ; $(LN) src/ripcut . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

intel:
	(cd src/ ; $(MAKE) all \
	"CF      = ifort" \
	"FFLAGS  = -I. -convert big_endian" \
	"FFLAGS2 = -I. -convert big_endian" \
	"FFLAGS3 = -I. -convert big_endian" \
	"CCFLAGS = -I. -DLITTLE -DUNDERSCORE -c" \
        "LDFLAGS = " \
	"LOCAL_LIBS = " \
	"NETCDFLIB = /usr/local/netcdf-intel/lib" \
	"NETCDFINC = /usr/local/netcdf-intel/include" \
	"LIBS    = -L$(NCARG_ROOT)/lib -lncarg -lcgm -lncarg_gks -lncarg_c -L/usr/X11R6/lib -lX11 -L/usr/lib/gcc-lib/i386-redhat-linux/3.2.3 -L/usr/lib -lf2c" )
	( $(RM) $(RM_LIST_EXE) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrfarw . ; $(LN) src/ripdp_wrfnmm . ; $(LN) src/ripcomp . ; $(LN) src/ripcut . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

mac_xlf:
	(cd src/ ; $(MAKE) all \
	"CF      = xlf" \
	"FFLAGS  = -qextname" \
	"FFLAGS2 = -qextname" \
	"FFLAGS3 = -qextname=premaptform:maptform:fillarray:mconvert:writefile_rdp:virtual:xtodot" \
	"CCFLAGS = -DLITTLE -DUNDERSCORE" \
        "LDFLAGS = -qarch=auto -qmaxmem=-1 -qblankpad -Wl,-stack_size,10000000,-stack_addr,0xc0000000" \
	"LOCAL_LIBS = " \
	"NETCDFLIB = /usr/local/netcdf-xlf/lib" \
	"NETCDFINC = /usr/local/netcdf-xlf/include" \
	"LIBS    = -L$(NCARG_ROOT)/lib -lncarg -lcgm -lncarg_gks -lncarg_c -L/usr/X11R6/lib -lX11 -lm -L/usr/local/lib -L/opt/ibmcmp/xlf/8.1/lib/ -lxlf90 -lg2c" )
	( $(RM) $(RM_LIST_EXE) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrfarw . ; $(LN) src/ripdp_wrfnmm . ;$(LN) src/ripcomp . ; $(LN) src/ripcut . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

mac:
	(cd src/ ; $(MAKE) all \
	"CF      = pgf90" \
	"FFLAGS  = -byteswapio -g " \
	"FFLAGS2 = -byteswapio -g " \
	"FFLAGS3 = -byteswapio -g " \
	"CCFLAGS = -DLITTLE -DUNDERSCORE -c -g" \
        "LDFLAGS = " \
	"LOCAL_LIBS = -L/usr/X11R6/lib -lX11 -L/usr/local/gfortran -lgfortran" \
	"NETCDFLIB = ${NETCDF}/lib -lnetcdf" \
	"NETCDFINC = ${NETCDF}/include" \
	"LIBS    = -L$(NCARG_ROOT)/lib -lncarg -lcgm -lncarg_gks -lncarg_c -L/usr/X11R6/lib -lX11" )
	( $(RM) $(RM_LIST_EXE) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrfarw . ; $(LN) src/ripdp_wrfnmm . ;$(LN) src/ripcomp . ; $(LN) src/ripcut . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

sun:
	(cd src/ ; $(MAKE) all \
	"CF      = f77" \
	"FFLAGS  = -fast -xarch=v8plus -xvector=yes -xprefetch=yes" \
	"FFLAGS2 = " \
	"FFLAGS3 = -fast -xarch=v8plus -xvector=yes -xprefetch=yes" \
	"CCFLAGS = -DUNDERSCORE -c" \
        "LDFLAGS = -fast -xarch=v8plus -xvector=yes -xprefetch=yes" \
	"LOCAL_LIBS = " \
	"NETCDFLIB = /usr/local/netcdf/lib" \
	"NETCDFINC = /usr/local/netcdf/include" \
	"LIBS    = -L$(NCARG_ROOT)/lib -R$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lX11 -lm -lM77 -lF77 -lsunmath" )
#	"LIBS    = -L$(NCARG_ROOT)/lib -R$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lncarg_loc -lX11 -lm -lM77 -lF77 -lsunmath" )
	( $(RM) $(RM_LIST_EXE) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrfarw . ; $(LN) src/ripdp_wrfnmm . ; $(LN) src/ripcomp . ; $(LN) src/ripcut . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

sun2:
	(cd src/ ; $(MAKE) all \
	"CF      = f77" \
	"FFLAGS  = -fast -xarch=v8plus -xvector=yes -xprefetch=yes" \
	"FFLAGS2 = " \
	"FFLAGS3 = -fast -xarch=v8plus -xvector=yes -xprefetch=yes" \
	"CCFLAGS = -DUNDERSCORE -c" \
        "LDFLAGS = -fast -xarch=v8plus -xvector=yes -xprefetch=yes" \
	"LOCAL_LIBS = " \
	"NETCDFLIB = /usr/local/netcdf/lib" \
	"NETCDFINC = /usr/local/netcdf/include" \
	"LIBS    = -L$(NCARG_ROOT)/lib -L/usr/openwin/lib -L/usr/dt/lib -lncarg -lncarg_gks -lncarg_c -lX11 -lm -lM77 -lF77 -lsunmath" )
	( $(RM) $(RM_LIST_EXE) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrfarw . ; $(LN) src/ripdp_wrfnmm . ; $(LN) src/ripcomp . ; $(LN) src/ripcut . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

sun90:
	(cd src/ ; $(MAKE) all \
	"CF      = f90" \
	"FFLAGS  = " \
	"FFLAGS2 = " \
	"FFLAGS3 = " \
	"CCFLAGS = -DUNDERSCORE -c" \
        "LDFLAGS = " \
	"LOCAL_LIBS = " \
	"NETCDFLIB = /usr/local/netcdf/lib" \
	"NETCDFINC = /usr/local/netcdf/include" \
	"LIBS    = -L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lX11 -lm -lsocket -lM77 -lF77" )
	( $(RM) $(RM_LIST_EXE) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrfarw . ; $(LN) src/ripdp_wrfnmm . ; $(LN) src/ripcomp . ; $(LN) src/ripcut . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

sgi:
	(cd src/ ; $(MAKE) all \
	"CF      = f77" \
	"FFLAGS  = -n32 -O3" \
	"FFLAGS2 = " \
	"FFLAGS3 = -n32 -O3" \
	"CCFLAGS = -DUNDERSCORE -c" \
        "LDFLAGS = " \
	"LOCAL_LIBS = " \
	"NETCDFLIB = /usr/local/netcdf/lib" \
	"NETCDFINC = /usr/local/netcdf/include" \
	"LIBS    = -L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lX11" )
	( $(RM) $(RM_LIST_EXE) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrfarw . ; $(LN) src/ripdp_wrfnmm . ; $(LN) src/ripcomp . ; $(LN) src/ripcut . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

sgi64:
	(cd src/ ; $(MAKE) all \
	"CF      = f77" \
	"FFLAGS  = -64 -O3" \
	"FFLAGS2 = -64 -O1" \
	"FFLAGS3 = -64 -O3" \
	"CCFLAGS = -DUNDERSCORE -c" \
        "LDFLAGS = " \
	"LOCAL_LIBS = " \
	"NETCDFLIB = /usr/local/netcdf/lib" \
	"NETCDFINC = /usr/local/netcdf/include" \
	"LIBS    = -L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lX11" )
	( $(RM) $(RM_LIST_EXE) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrfarw . ; $(LN) src/ripdp_wrfnmm . ; $(LN) src/ripcomp . ; $(LN) src/ripcut . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

ibm:
	(cd src/ ; $(MAKE) all \
	"CF      = xlf" \
	"FFLAGS  = " \
	"FFLAGS2 = " \
	"FFLAGS3 = " \
	"CCFLAGS = -DLITTLE" \
        "LDFLAGS = -O3 -qarch=auto -qmaxmem=-1 -qblankpad" \
	"LOCAL_LIBS = " \
	"NETCDFLIB = /usr/local/netcdf/lib" \
	"NETCDFINC = /usr/local/netcdf/include" \
	"LIBS    = -L$(NCARG_LIB) -lncarg -lncarg_gks -lncarg_c -lX11 -lm" )
	( $(RM) $(RM_LIST_EXE) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrfarw . ; $(LN) src/ripdp_wrfnmm . ; $(LN) src/ripcomp . ; $(LN) src/ripcut . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

cray:
	(cd src/ ; $(MAKE) all \
	"CF      = ncargf77" \
	"FFLAGS  = " \
	"FFLAGS2 = " \
	"FFLAGS3 = " \
	"CCFLAGS = -D_CRAY" \
        "LDFLAGS = " \
	"LOCAL_LIBS = " \
	"NETCDFLIB = /usr/local/netcdf/lib" \
	"NETCDFINC = /usr/local/netcdf/include" \
	"LIBS    = " )
	( $(RM) $(RM_LIST_EXE) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrfarw . ; $(LN) src/ripdp_wrfnmm . ; $(LN) src/ripcomp . ; $(LN) src/ripcut . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

vpp300:
	(cd src/ ; $(MAKE) all \
	"CF      = frt" \
	"FFLAGS  = -sc -Kfreealloc -Karraystack3 -X9 -Sw -Wv,-Of,-te,-ilfunc,-noalias,-m3,-P255 -Oe,-Kfast" \
	"FFLAGS2 = " \
	"FFLAGS3 = -sc -Kfreealloc -Karraystack3 -X9 -Sw -Wv,-Of,-te,-ilfunc,-noalias,-m3,-P255 -Oe,-Kfast" \
	"CCFLAGS = -DUNDERSCORE" \
        "LDFLAGS = -J" \
	"LOCAL_LIBS = " \
	"NETCDFLIB = /usr/local/netcdf/lib" \
	"NETCDFINC = /usr/local/netcdf/include" \
	"LIBS    = -L/usr/X11R6/lib -Wl,-dy /usr/X11R6/lib/libX11.so.6 -lsocket -lnsl -L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lm /nwp/nwptest/lib/libnwp.a /package/dms/librdms.a /package/dms/libgdbm.a" )
	( $(RM) $(RM_LIST_EXE) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrfarw . ; $(LN) src/ripdp_wrfnmm . ; $(LN) src/ripcomp . ; $(LN) src/ripcut . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

vpp5000:
	(cd src/ ; $(MAKE) all \
	"CF      = frt" \
	"FFLAGS  = -sc -Kfreealloc -Karraystack3 -X9 -Sw -Wv,-Of,-te,-ilfunc,-noalias,-m3,-P255 -Oe -KA32" \
	"FFLAGS2 = -KA32" \
	"FFLAGS3 = -sc -Kfreealloc -Karraystack3 -X9 -Sw -Wv,-Of,-te,-ilfunc,-noalias,-m3,-P255 -Oe -KA32" \
	"CCFLAGS = -KA32 -DUNDERSCORE" \
        "LDFLAGS = -J -KA32" \
	"LOCAL_LIBS = " \
	"NETCDFLIB = /usr/local/netcdf/lib" \
	"NETCDFINC = /usr/local/netcdf/include" \
	"LIBS    = -L/usr/lib32EX -Wl,-dy -lX11 -lSM -lICE -lgen -ldl -lsocket -lw -lnsl -L$(NCARG_ROOT)/lib -lncarg -lncarg_gks -lncarg_c -lm" )
	( $(RM) $(RM_LIST_EXE) ; $(LN) src/rip . ; $(LN) src/ripdp_mm5 . ; $(LN) src/ripdp_wrfarw . ; $(LN) src/ripdp_wrfnmm . ; $(LN) src/ripcomp . ; $(LN) src/ripcut . ; $(LN) src/ripinterp . ; $(LN) src/ripshow . ; $(LN) src/showtraj . ; $(LN) src/tabdiag . ; $(LN) src/upscale . )

clean:
	$(RM) $(RM_LIST_EXE) ; cd src/ ; $(RM) $(RM_LIST_EXE)

clobber:
	$(RM) $(RM_LIST_EXE) ; cd src/ ; $(RM) $(RM_LIST_EXE) $(RM_LIST_OBJ)
