      program rip
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c   RIP is a post-processing program for the PSU/NCAR mesoscale       c
c   model, written by Mark Stoelinga at the University of             c
c   Washington and NCAR.  The name "RIP" stands for                   c
c   Read/Interpolate/Plot. The program is intended for use with       c
c   hydrostatic or nonhydrostatic mm4v8 output, mm5v0 output,         c
c   mm5v1 input or output, or mm5v3 input or output (including        c
c   pre-processor p-level data).  The data must first be              c
c   converted to the appropriate format expected by this              c
c   program, using the program RIPDP, which stands for RIP data       c
c   preparation.                                                      c
c                                                                     c
c   RIP Version 3 released Oct 2000.                                  c
c   RIP Version 2 released Dec 1997.                                  c
c   RIP was first coded informally by Mark Stoelinga in 1991.         c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   The main program does very little "work".  It simply reads
c   a data header record in order to give dimensions to
c   subroutine driver, which does most of the "work" in RIP.
c   Some of the arguments to subroutine driver are used as
c   dimensions of large arrays, which are local to
c   subroutine driver.  The purpose of this set up is
c   for dynamic memory allocation.  This capability is
c   not standard to Fortran 77, but is allowable on some f77
c   compilers (such as Cray).  It is standard to Fortran 90.
c   Therefore, RIP may be compiled with f77 compilers that allow
c   adjustable dimensioning of local (non-argument) arrays in
c   subroutines, or with any f90 compiler.
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c   Note: if you don't have a Fortran 90 compiler, or your
c   Fortran 77 compiler doesn't support adjustable dimensioning of
c   local (non-argument) arrays, then make the appropriate changes
c   in subroutine driver (see the comments at the top of subroutine
c   driver), and recompile the program.
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c   The following parameters should not need to be changed for most
c   applications:
c
c      maxtraj: the maximum number of trajectories that can be
c         calculated in trajectory calculation mode
c      maxtavl: the maximum number of available times that can be
c         read in from the .xtimes file
c
      parameter (maxtraj=1000,maxtavl=500)
c
      dimension xtimeavl(maxtavl)
c
      character argum(16)*256,fname*256,cxtimeavl(maxtavl)*9,
     &   casename*256
c
c   Namelist variables
c
      parameter (maxptimes=500)
      dimension ptimes(maxptimes),iptimes(maxptimes)
      character title*80,rootname*256,titlecolor*40,rip_root*256,
     &   ptimeunits*1
      namelist/userin/ title,rip_root,rootname,flmin,frmax,fbmin,
     &   ftmax,ptimes,iptimes,ptimeunits,tacc,mdatebf,
     &   ntextq,ntextcd,idotser,idotitle,timezone,iusdaylightrule,
     &   inearesth,iinittime,ivalidtime,
     &   itrajcalc,fcoffset,titlecolor,idescriptive,
     &   icgmsplit,maxfld,iusectrv,imakev5d
      dimension xjtraj(maxtraj),yitraj(maxtraj),zktraj(maxtraj),
     &   diag(maxtraj)
      character vctraj*1
      namelist/trajcalc/ rtim,ctim,dtfile,dttraj,
     &   vctraj,ihydrometeor,xjtraj,yitraj,zktraj
c
      dimension ptuse(maxptimes)
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32),fullsigma(128),halfsigma(128)
      character chrip(64)*64,vardesc*64,plchun*24
c
      include 'comconst'
c
c   Get command line arguments.
c
      nargum=iargc()
      do i=1,nargum
         call getarg(i,argum(i))
      enddo
c
c   Fix for machines (such as HP) that return the command name itself
c   as the first element of argum, rather than the first argument.
c
      if (argum(1)(1:4).eq.'rip_') then
         do i=1,nargum-1
            argum(i)=argum(i+1)
         enddo
         nargum=nargum-1
      endif
c
      iup=6         ! output unit# for standard rip printout
      if (nargum.eq.3.and.argum(1)(1:2).eq.'-f') then
         iup=10
         argum(1)=argum(2)
         argum(2)=argum(3)
      elseif (nargum.ne.2) then
         write(6,*) 'Usage:'
         write(6,*) '   rip [-f] model_case_name rip_case_name'
         write(6,*) 'Options:'
         write(6,*) '   -f: standard rip printout should go to a file'
         write(6,*) '       named rip_case_name.out rather than to'
         write(6,*) '       the screen.'
         write(6,*) 'Note: "rip_case_name" should be the root part of'
         write(6,*) 'the name of the ".in" file.  Do not inlcude ".in"'
         write(6,*) 'at the end of "rip_case_name" on the'
         write(6,*) 'rip command line.'
         stop
      endif
      iupt=iup
c
c   Get model case name from argument
c
      casename=argum(1)
      iendc=index(casename,' ')-1
      if (iendc.eq.-1) iendc=256
c
c   Get root name from argument
c
      rootname=argum(2)
      iendcr=index(rootname,' ')-1
      if (iendcr.eq.1) then
         write(6,*) 'rip case name all blanks.'
         stop
      elseif (iendcr.gt.236) then
         write(6,*) 'rip case name must be < or = 236 characters.'
         stop
      endif
      if (rootname(iendcr-2:iendcr).eq.'.in') then
         rootname(iendcr-2:iendcr)='   '
         iendcr=iendcr-3
      endif
c
c   Define constants.  Many are taken from Bolton (1980, MWR 108,1046-1053). 
c
      rgas=287.04  !J/K/kg
      rgasmd=.608   ! rgas_moist=rgas*(1.+rgasmd*qvp)
      cp=1004.     ! J/K/kg  Note: not using Bolton's value of 1005.7
      cpmd=.887   ! cp_moist=cp*(1.+cpmd*qvp)
      gamma=rgas/cp
      gammamd=rgasmd-cpmd  ! gamma_moist=gamma*(1.+gammamd*qvp)
      grav=9.81           ! m/s**2
      sclht=rgas*256./grav   ! 256 K is avg. trop. temp. from USSA.
      eps=0.622
      ezero=6.112  ! hPa
      xlhc0=3.1484e6   ! J/kg
      xlhctd=2370.  !
      xlhf=3.34e5
      pvc=1.e6
      rktpmps=1.94
      celkel=273.15
      eslcon1=17.67
      eslcon2=29.65
      esicon1=22.514
      esicon2=6.15e3
      thtecon1=3376. ! K
      thtecon2=2.54
      thtecon3=.81
      tlclc1=2840.
      tlclc2=3.5
      tlclc3=4.805
      tlclc4=55.
      rhoice=917.
      rhowat=1000.
      pi=4.*atan(1.)
      rpd=pi/180.
      abscoef=.145      ! cloud water absorption coefficient in m^2/g
      abscoefi=.272     ! cloud ice absorption coefficient in m^2/g
      ussalr=.0065      ! deg C per m
      rmsg=9.0e+9       ! indicates missing data or specification
c
c   Define unit numbers
c
      iuinput=7     ! input unit# for user input (namelists & plspecs)
      iuxtavl=9     ! input unit# for available xtimes file
      iustnlist=8   ! input unit# for the station list
      iudata=21     ! input unit# for the data.
      iutserstn=70  ! input unit# for the time series stations
      iutserdat=71  ! output unit# for the time series printout
      iutrajin=74   ! input unit# for trajc positns (traj plot mode)
      iutrajout=75  ! output unit# for trajc positns (traj calc mode)
      iudiagout=76  ! output unit# for trajc diagnstcs (traj calc mode)
      iuv5dout=79   ! output unit# for vis5d files (make vis5d mode)
      iuprcver=95   ! output unit# for the precip verif. results
c
c   Open the print out file
c
      if (iup.eq.10) then
         fname=rootname(1:iendcr)//'.out'
         open (unit=iup,file=fname,form='formatted',status='unknown')
      endif
c
c   Read the namelist values.
c
      fname=rootname(1:iendcr)//'.in'
      open (unit=iuinput,file=fname,form='formatted',status='old')
      fname=rootname   ! save rootname temporarily
c
      title='auto'
      rip_root='/dev/null'
      rootname=' '      ! no longer set in namelist
      mdatebf=99999999  ! no longer used
      flmin=.05
      frmax=.95
      fbmin=.10
      ftmax=.90
      do i=1,maxptimes
         ptimes(i)=9e9
         iptimes(i)=99999999
      enddo
      ptimeunits='h'
      tacc=1.
      ntextq=0
      ntextcd=0
      idotser=0
      idotitle=1
      timezone=-7.
      iusdaylightrule=1
      inearesth=0
      iinittime=1
      ivalidtime=1
      idescriptive=1
      itrajcalc=0
      imakev5d=0
      iusectrv=0
      fcoffset=0.
      titlecolor=' '
      icgmsplit=0
      maxfld=10
c
      read (iuinput,userin)
      rootname=fname
c
      if (itrajcalc.eq.1.and.imakev5d.eq.1) then
         write(6,*) 'You cannot be in trajectory calculation mode',
     &          ' (itrajcalc=1)'
         write(6,*) 'and "make vis5d data" mode (imakev5d=1) at the',
     &          ' same time.'
         write(6,*) 'You must choose one or the other.'
      endif
c
      rtim=0.0
      ctim=0.0
      dtfile=3600.0   ! seconds
      dttraj=600.0    ! seconds
      vctraj='s'
      ihydrometeor=0
      do i=1,maxtraj
         xjtraj(i)=rmsg
         yitraj(i)=rmsg
         zktraj(i)=rmsg
      enddo
c
      if (itrajcalc.eq.1) then
         read (iuinput,trajcalc)
c
c      Check if a grid is defined
c
         if (xjtraj(1).lt.0.0) call mktrjpts(xjtraj,yitraj,zktraj,
     &                                        maxtraj,rmsg)
         do i=maxtraj,1,-1
            if (xjtraj(i).ne.rmsg) then
               ntraj=i
               goto 40
            endif
         enddo
 40      continue
      else
         ntraj=1
      endif
c
c   Get rip_root from environment variable RIP_ROOT
c   if rip_root='/dev/null '
c
      if (rip_root(1:10).eq.'/dev/null ') then
         call getenv('RIP_ROOT',rip_root)
      endif
      if (rip_root(1:10).eq.'          ') then
         write(iup,*)'Either RIP_ROOT environment variable or rip_root'
         write(iup,*)'is not properly set.'
         stop
      endif
c
c   Set up lookup table for getting temperature on a pseudoadiabat.
c   (Borrow the unit number for the stationlist, just for the moment.)
c
      iendrr=index(rip_root,' ')-1
      fname=rip_root(1:iendrr)//'/psadilookup.dat'
      open (unit=iustnlist,file=fname,form='formatted',status='old')
      do i=1,14
         read(iustnlist,*)
      enddo
      read(iustnlist,*) nthte,nprs
      if (nthte.ne.150.or.nprs.ne.150) then
         write(6,*)
     &      'Number of pressure or theta_e levels in lookup table'
         write(6,*) 'file not = 150.  Check lookup table file.'
         stop
      endif
      read(iustnlist,173) (psadithte(jt),jt=1,nthte)
      read(iustnlist,173) (psadiprs(ip),ip=1,nprs)
      read(iustnlist,173) ((psaditmk(ip,jt),ip=1,nprs),jt=1,nthte)
 173  format(5e15.7)
      close(iustnlist)
c
c   Read in the available xtimes, in both character and
c   floating point arrays
c
      fname=casename(1:iendc)//'.xtimes'
      open(unit=iuxtavl,file=fname,form='formatted',status='old')
      read(iuxtavl,*)nxtavl
      if (nxtavl.gt.maxtavl) then
         write(iup,*)'There are ',nxtavl,' times in the ".xtime"'
         write(iup,*)'file, but maxtavl is only ',maxtavl,'.'
         write(iup,*)'Increase maxtavl in rip.f to be at least'
         write(iup,*)'as big as ',nxtavl,' and then recompile'
         write(iup,*)'and re-run rip.'
         stop
      endif
      do i=1,nxtavl
         read(iuxtavl,'(a9)') cxtimeavl(i)
         read(cxtimeavl(i),'(f9.5)') xtimeavl(i)
      enddo
      close (iuxtavl)
      do i=nxtavl+1,maxtavl
         cxtimeavl(i)=' '
         xtimeavl(i)=0.
      enddo
c
c   Create a file name for a file that is likely to exist,
c   just so we can read it and get the dimensions of the dataset.
c
      fname=casename(1:iendc)//'_'//cxtimeavl(1)//'_'//'ter'
      open(unit=iudata,file=fname,form='unformatted',status='unknown')
      read(iudata,err=170,end=180)
     &   vardesc,plchun,ihrip,rhrip,chrip,fullsigma,halfsigma
      miy=ihrip(4)
      mjx=ihrip(5)
      mkzh=ihrip(9)
      close (iudata)
      goto 200
c
 170  write(6,*) 'The model data header is not a format'//
     &   ' that RIP recognizes.  Stopping.'
      stop
c
 180  write(6,*) 'Unexpected EOF reached when trying to read'
      write(6,*) 'model data header.  Stopping.'
      stop
c
 200  continue
c
      rewind(iuinput)
      call driver(miy,mjx,mkzh,xtimeavl,cxtimeavl,maxtavl,nxtavl,
     &   casename,iendc,ihrip,rhrip,fullsigma,halfsigma,
     &   chrip,vardesc,plchun,title,rip_root,rootname,iendcr,
     &   ptimes,iptimes,ptuse,maxptimes,ptimeunits,tacc,ntextq,ntextcd,
     &   idotser,idotitle,timezone,iusdaylightrule,inearesth,iinittime,
     &   ivalidtime,fcoffset,
     &   titlecolor,idescriptive,icgmsplit,maxfld,itrajcalc,
     &   iusectrv,rtim,ctim,dtfile,dttraj,vctraj,ihydrometeor,
     &   xjtraj,yitraj,zktraj,diag,ntraj,imakev5d)
c
      stop
      end
