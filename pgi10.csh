#!/bin/csh

#	This is a script to have a specif session use the PGI
#	compilers on the MAC Intel (where I have the default set
#	up for g95).  To use this script, in the window/session that
#	will use PGI, source this file, then do a rehash.  Note also
#	included is the required location for the gcc lib.

# Wei:  since Michael's lib is compiled with 64-bit compiler, we need to use it too.
#       at least for ungrib.exe

#setenv PGI /usr/local/pgi_7.1-6
#setenv PGI /usr/local/pgi_10.3
setenv PGI /usr/local/pgi

#set path = ( /usr/local/pgi_7.1-6/bin-64 /usr/local/mpich2-1.0.6p1-pgi/bin \
#set path = ( /usr/local/pgi_8.0-4/bin-64 /users/duda/bin/gambel/mpich2-1.0.7/bin \
#             /usr/local/netcdf-pgi/bin $path )
set path = ( $PGI/bin \
             /usr/local/mpich2-1.2.1p1-pgi/bin \
             /usr/local/netcdf-3.6.0-p1-pgi32/bin $path )

setenv LD_LIBRARY_PATH /usr/lib:/usr/local/lib:${PGI}/lib-64
#setenv LD_LIBRARY_PATH /usr/lib:/usr/local/lib:${PGI}/lib-64:/users/duda/bin/gamble/mpich2-1.0.7/lib:/usr/local/netcdf-pgi/lib

setenv LM_LICENSE_FILE 7496@licenseb.ucar.edu:7496@licensea.ucar.edu

#setenv NETCDF /usr/local/netcdf-pgi
#setenv NETCDF /usr/local/netcdf-4.1.1-pgi32
setenv NETCDF /usr/local/netcdf-3.6.0-p1-pgi32
setenv NCARG_ROOT /usr/local/ncl

#	Need this on the LIB line in the configure.defaults:

#LIB = $(LIB_BUNDLED) $(LIB_EXTERNAL) $(LIB_LOCAL) -L/usr/lib/gcc/i686-apple-darwin8/4.0.1/x86_64 -lgcc

#	And another thing - I seem to need to use my mpich2
#	for running parallel jobs.

#/usr/local/mpich2-1.2.1p1-pgi64/bin/mpd &
#/usr/local/mpich2-1.2.1p1-pgi/bin/mpd &
#~duda/bin/gambel/mpich2-1.0.7/bin/mpd &
#~duda/bin/gambel/mpich2-1.0.7/bin/mpirun -np 4 wrf.exe &

# /stink/gill/local/bin/mpirun -np 2 wrf.exe

#	And yet another thing, MAC Intel PGI needs a new
#	WRFV3/external/io_grib1/MEL_grib1/grib_lookup.h
#	Diffs:
#+#ifdef MACOS
#+#include "/usr/lib/gcc/i686-apple-darwin8/4.0.1/include/float.h"
#+#elif
# #include <float.h>
#+#endif

