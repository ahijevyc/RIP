      program ripdp_wrf
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c   RIPDP_WRF is a data preparation program that reads in data from   c
c   the WRF modeling system and makes data files                      c
c   appropriate for input to the RIP data analysis and visulaization  c
c   program.  RIP stands for Read/Interpolate/Plot and RIPDP stands   c
c   for RIP data preparation.                                         c
c                                                                     c
c   FYI: Output unit number is 65.                                    c
c                                                                     c
c   Originally written by Mark Stoelinga, Univ. of Washington         c
c                                                                     c
c   Modified Feb 2003 for new, generalized vertical coordinate        c
c   version of RIP.  See code commented with "ccc".                   c
c                                                                     c
c   Modified March 2003 for WRF model output.                         c
c      [Wei Wang (NCAR) and Mark Stoelinga (UW)]                      c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c   The main program does very little "work".  It simply reads
c   gets some basic information from the metadata part of the
c   data file in order to give information to
c   subroutine process, which does most of the "work" in RIPDP.
c   Some of the arguments to subroutine process are used as
c   dimensions of the 2- and 3-D arrays.  The purpose of this
c   set up is for dynamic memory allocation.  This capability is
c   not standard to Fortran 77, but is allowable on some f77
c   compilers (such as Cray).  It is standard to Fortran 90.
c   Therefore, RIPDP may be compiled with f77 compilers that allow
c   adjustable dimensioning of local (non-argument) arrays in
c   subroutines, or with any f90 compiler.
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c   Note: if you don't have a Fortran 90 compiler, or your
c   Fortran 77 compiler doesn't support adjustable dimensioning of
c   local (non-argument) arrays, then make the appropriate changes
c   in subroutine process (see the comments at the top of subroutine
c   process), and recompile the program.
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c   Map transform common block
c
      common /mptf/ rpd_mptf,pi_mptf,dskmc_mptf,xlonc_mptf,rearth_mptf,
     &   ciy_mptf,cjx_mptf,yc_mptf,ihm_mptf,c1_mptf,c2_mptf,cone_mptf,
     &   conei_mptf,nproj_mptf
c
c   netcdf variables
c
      integer ncid, dimid, nf_status
      character nf_att_text*64
c
      character argum(256)*256
c
      include "netcdf.inc"
c
c   Get command line arguments.
c
      nargum=iargc()
      do i=1,nargum
         argum(i)=' '
         call getarg(i,argum(i))
      enddo
c
c   Fix for machines (such as HP) that return the command name itself
c   as the first element of argum, rather than the first argument.
c
      if (argum(1)(1:10).eq.'ripdp_wrf_') then
         do i=1,nargum-1
            argum(i)=argum(i+1)
         enddo
         nargum=nargum-1
      endif
c
      if ((argum(1)(1:3).eq.'-n '.and.nargum.lt.5).or.
     &    (argum(1)(1:3).ne.'-n '.and.nargum.lt.3)) then
         print*,'Usage:'
         print*,'  ripdp_wrf [-n namelist_file] casename [basic|all]',
     &          ' data_file_1 data_file_2 data_file_3 ...'
         print*
         print*,'Note: "namelist_file" can be any name, either with'
         print*,'an extension of your choosing or without an extension.'
         print*,'The ".in" extension is not assumed, as it is in rip.'
         print*,'"basic" tells RIPDP to process only the basic'
         print*,'variables it expects to find.  "all" tells RIP to'
         print*,'process all variables it encounters.'
         stop
      endif
c
      if (argum(1)(1:3).eq.'-n ') then
         nnl=2
         ncn=3
         ndd=4
         nsets=nargum-4
         nsetsbeg=5
      else
         nnl=0
         ncn=1
         ndd=2
         nsets=nargum-2
         nsetsbeg=3
      endif
c
c   Open first netcdf file and get model dimensions
c
      nf_status = nf_open (argum(nsetsbeg), nf_nowrite, ncid)
      call handle_err(000,nf_status)
c
      nf_att_text=' '
      nf_status = nf_get_att_text (ncid, nf_global,
     &   'TITLE', nf_att_text)
      call handle_err(001,nf_status)
      if (index(nf_att_text,'OUTPUT FROM WRF').ne.0) then
         print*,'Data is recognized as WRF model output data.'
         print*
      else
         stop
      endif
c
      nf_status = nf_inq_dimid (ncid, 'south_north_stag', dimid)
      call handle_err(002,nf_status)
      nf_status = nf_inq_dimlen (ncid, dimid, miy)
      call handle_err(003,nf_status)
c
      nf_status = nf_inq_dimid (ncid, 'west_east_stag', dimid)
      call handle_err(004,nf_status)
      nf_status = nf_inq_dimlen (ncid, dimid, mjx)
      call handle_err(005,nf_status)
c
      nf_status = nf_inq_dimid (ncid, 'bottom_top', dimid)
      call handle_err(006,nf_status)
      nf_status = nf_inq_dimlen (ncid, dimid, mkzh)
      call handle_err(007,nf_status)
c
c   Close first netcdf file
c
      nf_status = nf_close (ncid)
      call handle_err(008,nf_status)
c
c   Run a loop to open all netcdf files, for the sole purpose of
c   getting max number of model output times in any one file.
c
      nwrftimes=0
      do i=1,nsets
         ind=nsetsbeg+i-1
         nf_status = nf_open (argum(ind), nf_nowrite, ncid)
         call handle_err(009.1,nf_status)
         nf_status = nf_inq_dimid (ncid, 'Time', dimid)
         call handle_err(009.2,nf_status)
         nf_status = nf_inq_dimlen (ncid, dimid, nwrftimes_this_file)
         call handle_err(009.3,nf_status)
         nwrftimes=max(nwrftimes,nwrftimes_this_file)
         nf_status = nf_close (ncid)
         call handle_err(009.4,nf_status)
      enddo
c
      call process(miy,mjx,mkzh,argum,nnl,ncn,ndd,nsets,nsetsbeg,
     &   nwrftimes)
      stop
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine process(miy,mjx,mkzh,argum,nnl,ncn,ndd,nsets,nsetsbeg,
     &   nwrftimes)
c
c   This subroutine does most of the "work".
c
c   miy, and mjx are dot-point dimensions, in the x and y directions
c      respectively, of the domain to be analyzed.
c   mkzh is number of model levels in the domain.
c   argum carries the names of the model data files.
c   nsets is the number of files in the model dataset.
c   nsetsbeg is the element of argum that holds the first file name
c      of the model dataset.
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c   Note: if you don't have a Fortran 90 compiler, or your
c   Fortran 77 compiler doesn't support adjustable dimensioning of
c   local (non-argument) arrays, then you must do the following:
c   comment out the above subroutine declaration, and uncomment
c   the following one:
c
c      subroutine process(idumb1,idumb2,idumb3,argum,
c     &   nnl,ncn,ndd,nsets,nsetsbeg,idumb4)
c
c   Then, uncomment the following parameter statement and set the
c   parameters to the appropriate values for your model domain.
c
c      parameter (miy=118,mjx=94,mkzh=27,nwrftimes=5)
c
c   Then, recompile RIPDP to make an executable that is appropriate
c   for your domain.
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      character argum(256)*256
c
      include "netcdf.inc"
c
      dimension ter(miy,mjx),xmap(miy,mjx),xlat(miy,mjx),xlon(miy,mjx),
     &   cor(miy,mjx),xlus(miy,mjx),rtc(miy,mjx),rte(miy,mjx),
     &   tgk(miy,mjx),dmap(miy,mjx),tmk(miy,mjx,mkzh),
     &   uuu(miy,mjx,mkzh),ght(miy,mjx,mkzh),www(miy,mjx,mkzh+1),
     &   vvv(miy,mjx,mkzh),prs(miy,mjx,mkzh),sfp(miy,mjx),
     &   qvp(miy,mjx,mkzh),scr3wlev(miy,mjx,mkzh+1),
     &   scr3(miy,mjx,mkzh),scr2(miy,mjx)
c
      character varname*10,fname*256,cxtime*10,
     &   cxtimeavl(256)*10
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
c   netcdf variables
c
      real nf_uarr3(mjx,miy-1,mkzh),  nf_varr3(mjx-1,miy,mkzh),
     &     nf_tarr3(mjx-1,miy-1,mkzh),nf_warr3(mjx-1,miy-1,mkzh+1),
     &     nf_tarr2(mjx-1,miy-1)
      integer ncid, ndims, nvars, ngatts, unlimdimid, dimid,
     &   nf_status, varid, nf_att_int, iprocvarid(200)
      integer nf_ustart3(4),nf_vstart3(4),nf_tstart3(4),nf_wstart3(4),
     &   nf_tstart2(3)
      integer nf_ucount3(4),nf_vcount3(4),nf_tcount3(4),nf_wcount3(4),
     &   nf_tcount2(3),nf_tcount3s(4)
      real nf_att_real
      character nf_att_text*64, start_date*19, nf_varname*16,
     &   wrftimes(nwrftimes)*19
      real nf_znu(mkzh,nwrftimes),nf_znw(mkzh+1,nwrftimes)
      real znu(mkzh),znw(mkzh+1),znfac(mkzh)
      integer dimid_tm, dimid_we, dimid_sn, dimid_bt, dimid_sls
      integer vardimids(20), nf_att_len
c
c minfo string
c
      character minfostring*82
c
c   Namelist variables
c
      parameter (maxptimes=500)
      dimension ptimes(maxptimes),iptimes(maxptimes),ptuse(maxptimes)
      character discard(maxptimes)*16,retain(maxptimes)*16,ptimeunits*1
      namelist/userin/ ptimes,iptimes,ptimeunits,tacc,discard,retain,
     &   iexpandedout,iskpd1
c
      print*,'Welcome to your friendly RIPDP output file !'
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
      iuinput=7     ! input unit# for namelist, color table,
c                        and plspec table
      ifilecount=1
c
c   Read the namelist values.
c
      iskpd1=0
      do i=1,maxptimes
         ptimes(i)=9e9
         iptimes(i)=99999999
         discard(i)=' '
         retain(i)=' '
      enddo
      ptimeunits='h'
      tacc=1.
      iexpandedout=0
c
      if (nnl.ne.0) then
         open (unit=iuinput,file=argum(nnl),form='formatted',
     &         status='old')
         read (iuinput,userin)
      endif
      iexpanded=0
      if (iexpanded.eq.1) then
         print*,'Input data is on an expanded domain.'
         if (iexpandedout.eq.1) then
           print*,'RIPDP will process the full (expanded) domain.'
         else
           print*,'RIPDP will output the standard (unexpanded) domain.'
         endif
      endif
      tacch=tacc/3600.
c
      ibasic=0
      if (argum(ndd)(1:5).eq.'basic') ibasic=1
c
c   Determine number of names in "discard" array.
c
      do i=1,maxptimes
         if (discard(i).eq.'         ') then
            ndiscard=i-1
            goto 492
         endif
      enddo
      ndiscard=maxptimes
 492  continue
      do i=1,maxptimes
         if (retain(i).eq.'         ') then
            nretain=i-1
            goto 493
         endif
      enddo
      nretain=maxptimes
 493  continue
      if (nretain.ne.0.and.ndiscard.ne.0) then
         print*,'You cannot specify both a "discard" list and a'
         print*,'"retain" list.  Specify one or the other.'
         stop
      endif
      if (ibasic.eq.1.and.ndiscard.ne.0) then
         print*,'If you specify "basic" on the command line,'
         print*,'you should only specify a "retain" list (or no list)'
         print*,'in the namelist, but not a "discard" list.'
         stop
      endif
      if (ibasic.eq.0.and.nretain.ne.0) then
         print*,'If you specify "all" on the command line,'
         print*,'you should only specify a "discard" list (or no list)'
         print*,'in the namelist, but not a "retain" list.'
         stop
      endif
c
      iendc=index(argum(ncn),' ')-1
      if (iendc.eq.-1) iendc=256
      iendf1=iendc+12
      nxtavl=0
c
c   If using iptimes, convert the mdates in the iptimes array to
c   xtimes in the ptimes array.  Also, determine nptimes.
c
      if (ptimes(1).lt.0..or.iptimes(1).lt.0.or.(ptimes(1).eq.
     &    9e9.and.iptimes(1).eq.99999999)) then !user wants all times
         nptimes=0
         print*,'Note: RIPDP will process all times encountered.'
      elseif (ptimes(1).ne.9e9.and.iptimes(1).ne.99999999) then
         print*,'Can''t use both ptimes and iptimes.'
         stop
      elseif (iptimes(1).ne.99999999) then
         do i=1,maxptimes
            if (iptimes(i).eq.99999999) then
               nptimes=i-1
               goto 259
            else
               if (iptimes(i).lt.0) then
                  call mconvert(-iptimes(i),mhourp,1,1940)
                  ptimes(i)=-float(mhourp-mhourb)
               elseif (i.ge.2.and.iptimes(i-1).lt.0) then
                  ptimes(i)=float(iptimes(i))
               else
                  call mconvert(iptimes(i),mhourp,1,1940)
                  ptimes(i)=float(mhourp-mhourb)
               endif
            endif
         enddo
         nptimes=maxptimes
 259     continue
      else
         if (ptimeunits.eq.'h') then
            tunitfac=1.
         elseif (ptimeunits.eq.'m') then
            tunitfac=1./60.
         elseif (ptimeunits.eq.'s') then
            tunitfac=1./3600.
         endif
         do i=1,maxptimes
            if (ptimes(i).eq.9e9) then
               nptimes=i-1
               goto 261
            else
               ptimes(i)=tunitfac*ptimes(i)
            endif
         enddo
         nptimes=maxptimes
 261     continue
      endif
c
c   Process time sequences in ptimes array.
c
      ii=0
      itime=0
      ptusemax=-9e9
  400 ii=ii+1
      if (ii.gt.nptimes) goto 420
      if (ptimes(ii).ge.0.) then
         itime=itime+1
         if (itime.gt.maxptimes) then
            print*,'Number of times requested exceeds maxptimes.'
            print*,'Increase maxptimes in ripdp code, recompile,'
            print*,'and run ripdp again.'
            stop
         endif
         ptuse(itime)=ptimes(ii)
         ptusemax=max(ptuse(itime),ptusemax)
      else
         ii=ii+1
         if (ptimes(ii).gt.0.) then
            tstart=ptimes(ii-2)
            tend=-ptimes(ii-1)
            if (tend.eq.tstart) tend=tend+1.e-10
            tinc=ptimes(ii)
            tdist=tend-tstart
            isign=nint(tdist/abs(tdist))
            ntseries=int(abs(tdist)/tinc+.00001) + 1
            do i=2,ntseries
               itime=itime+1
               if (itime.gt.maxptimes) then
                  print*,
     &              'Number of times requested exceeds maxptimes.'
                  print*,
     &              'Increase maxptimes in ripdp code, recompile,'
                  print*,'and run ripdp again.'
                  stop
               endif
               ptuse(itime)=ptuse(itime-1)+isign*tinc
               ptusemax=max(ptuse(itime),ptusemax)
            enddo
         else
            print*,'Error in ptimes sequence specification.'
            stop
         endif
      endif
      goto 400
  420 nptuse=itime
      if (nptuse.eq.0) ptusemax=9e9
c
c=================================================================c
      do iwf=1,nsets         ! File loop
c=================================================================c
c
      ind=nsetsbeg+iwf-1
      nf_status = nf_open (argum(ind), nf_nowrite, ncid)
      call handle_err(009.5,nf_status)
      nf_status = nf_inq_dimid (ncid, 'Time', dimid)
      call handle_err(009.6,nf_status)
      nf_status = nf_inq_dimlen (ncid, dimid, nwrftimes_this_file)
      call handle_err(009.7,nf_status)
c
c   Get basic netcdf info
c
      nf_status = nf_inq (ncid, ndims, nvars, ngatts, unlimdimid)
      call handle_err(010,nf_status)
c
c     Initialize iprocvarid array to 0.  For each varid obtained for a
c     variable that gets processed, set that element of iprocvarid to 1.
c     This will allow ripdp to check all unprocessed variables to see if
c     a data file can be created.
c
      do i=1,200
         iprocvarid(i)=0
      enddo
c
c   Get some dimension IDs
c
      nf_status = nf_inq_dimid (ncid, 'Time', dimid_tm)
      call handle_err(011,nf_status)
      nf_status = nf_inq_dimid (ncid, 'west_east', dimid_we)
      call handle_err(012,nf_status)
      nf_status = nf_inq_dimid (ncid, 'south_north', dimid_sn)
      call handle_err(013,nf_status)
      nf_status = nf_inq_dimid (ncid, 'bottom_top', dimid_bt)
      call handle_err(014,nf_status)
      nf_status = nf_inq_dimid (ncid, 'soil_layers_stag', dimid_sls)
      call handle_err(015,nf_status)
c
c   Get array of model output dates/times
c
      nf_status = nf_inq_varid (ncid, 'Times', varid)
      call handle_err(016,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_var_text (ncid, varid, wrftimes)
      call handle_err(017,nf_status)
c
c   Get start date string and create RIP parameters for start time.
c
      nf_status = nf_get_att_text (ncid,nf_global,
     &   'START_DATE', start_date)
      call handle_err(018,nf_status)
c
      read(start_date,'(1x,6(1x,i2))')iyr,imo,idy,ihr,imn,isc
      mdateb=1000000*iyr+10000*imo+
     &   100*idy+ihr
      rhourb=imn/60.+isc/3600.
      call mconvert(mdateb,mhourb,1,1940)
c
c     Get the vertical "eta" coordinate values of both the "w" and
c     "mass" levels, and re-order in RIP's top-bottom order.
c
      nf_status = nf_inq_varid (ncid, 'ZNW', varid)
      call handle_err(019,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_var_real (ncid, varid, nf_znw)
      call handle_err(020,nf_status)
      do k=1,mkzh+1
         znw(k)=nf_znw(mkzh+1-k+1,1)
      enddo
c
      nf_status = nf_inq_varid (ncid, 'ZNU', varid)
      call handle_err(021,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_var_real (ncid, varid, nf_znu)
      call handle_err(022,nf_status)
      do k=1,mkzh
         znu(k)=nf_znu(mkzh-k+1,1)
         znfac(k)=(znu(k)-znw(k))/(znw(k+1)-znw(k))
      enddo
c
c   Set values of "count" and "start" arrays for reading netcdf data.
c
      nf_ustart3(1)=1
      nf_ustart3(2)=1
      nf_ustart3(3)=1
      nf_ustart3(4)=1     ! will be reset at each time in time loop
      nf_ucount3(1)=mjx
      nf_ucount3(2)=miy-1
      nf_ucount3(3)=mkzh
      nf_ucount3(4)=1
      nf_vstart3(1)=1
      nf_vstart3(2)=1
      nf_vstart3(3)=1
      nf_vstart3(4)=1     ! will be reset at each time in time loop
      nf_vcount3(1)=mjx-1
      nf_vcount3(2)=miy
      nf_vcount3(3)=mkzh
      nf_vcount3(4)=1
      nf_tstart3(1)=1
      nf_tstart3(2)=1
      nf_tstart3(3)=1
      nf_tstart3(4)=1     ! will be reset at each time in time loop
      nf_tcount3(1)=mjx-1
      nf_tcount3(2)=miy-1
      nf_tcount3(3)=mkzh
      nf_tcount3(4)=1
      nf_wstart3(1)=1
      nf_wstart3(2)=1
      nf_wstart3(3)=1
      nf_wstart3(4)=1     ! will be reset at each time in time loop
      nf_wcount3(1)=mjx-1
      nf_wcount3(2)=miy-1
      nf_wcount3(3)=mkzh+1
      nf_wcount3(4)=1
      nf_tstart2(1)=1
      nf_tstart2(2)=1
      nf_tstart2(3)=1     ! will be reset at each time in time loop
      nf_tcount2(1)=mjx-1
      nf_tcount2(2)=miy-1
      nf_tcount2(3)=1
c
      nf_status = nf_inq_dimid (ncid, 'soil_layers_stag', dimid)
      call handle_err(023,nf_status)
      nf_status = nf_inq_dimlen (ncid, dimid, nsoil)
      call handle_err(024,nf_status)
      nf_tcount3s(1)=mjx-1
      nf_tcount3s(2)=miy-1
      nf_tcount3s(3)=nsoil
      nf_tcount3s(4)=1
c
c   Get information from netcdf global attributes
c
      nf_status = nf_get_att_int (ncid, nf_global,
     &   'SOUTH-NORTH_GRID_DIMENSION', nf_att_int)
      call handle_err(025,nf_status)
      miycors=nf_att_int
c
      nf_status = nf_get_att_int (ncid, nf_global,
     &   'WEST-EAST_GRID_DIMENSION', nf_att_int)
      call handle_err(026,nf_status)
      mjxcors=nf_att_int
c
c      if (iexpanded.eq.1) then
c         miycors=miy
c         mjxcors=mjx
c         ioffexp=??
c         joffexp=??
c      endif
c
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'DX', nf_att_real)   ! DX is in meters
      call handle_err(027,nf_status)
      dskm=.001*nf_att_real
c
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'CEN_LAT', nf_att_real)
      call handle_err(028,nf_status)
      xlatc=nf_att_real
c
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'CEN_LON', nf_att_real)
      call handle_err(029,nf_status)
      xlonc=nf_att_real
c
      nf_status = nf_get_att_int (ncid, nf_global,
     &   'MAP_PROJ', nf_att_int)
      call handle_err(030,nf_status)
      nproj=nf_att_int
      if (nproj.gt.3.or.nproj.lt.0) then
         print*,'   Map proj. #',nproj,' is not recognized.'
         stop
      endif
c
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'TRUELAT1', nf_att_real)
      call handle_err(031,nf_status)
      truelat1=nf_att_real
c
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'TRUELAT2', nf_att_real)
      call handle_err(032,nf_status)
      truelat2=nf_att_real
c
      dskmc=dskm  ! grid spacing (in km) of coarsest domain
      dsc=dskmc*1000.
c
c   Now set up map transformation stuff
c
      if (nproj.eq.3) then   ! Mercator
         true1=0.
         true2=0.
      elseif (nproj.eq.1.and.
     &    (truelat1.gt.90..or.truelat2.gt.90.)) then  ! LC default
         true1=sign(30.,xlatc)
         true2=sign(60.,xlatc)
      elseif (nproj.eq.1) then   ! Lambert Conformal
         true1=truelat1
         true2=truelat2
      elseif (nproj.eq.2.and.truelat1.gt.90.) then   ! ST default
         true1=sign(60.,xlatc)
         true2=true1
      elseif (nproj.eq.2) then   ! Polar Sterographic
         true1=truelat1
         true2=true1
      else
         print*,'Unrecognized map projection.'
         stop
      endif
      call premaptform(dskmc,miycors,mjxcors,nproj,
     &   xlatc,xlonc,true1,true2,6)
c
c   Things for this domain
c
      ds=dskm*1000.
      refrat=dsc/ds  ! refinement ratio of this grid w.r.t. coarse grid
c
c     Note, yicorn and xjcorn, the y and x locations of the lower left
c     corner of the current domain with respect to the centered coarse
c     domain (i.e. the coarse domain centered on lat=xlatc,lon=xlonc)
c     are not set here.  They are set at each time (within the time loop),
c     to account for the future possibility of a moving nest.
c
c   Land use data set:
c
      nf_att_text=' '
      nf_status = nf_get_att_text (ncid, nf_global,
     &   'MMINLU', nf_att_text)
      call handle_err(037,nf_status)
      if (nf_att_text(1:5).eq.'USGS ') then
         ilandset=2    ! WRF SI only uses USGS 24-category land use data set
                       !   right now.
      else
         print*,'Unexpected land use data set specified.'
         stop
      endif
c
c   Other
c
      iplevdata=8   ! Anything larger than 3 means terrain-following data
c
c     Determine iice.  In RIP, iice=1 means that there are separate
c     arrays for frozen hydrometeors (as would be created by a
c     mixed-phase bulk scheme), whereas ice=0 means that frozen and
c     liquid hydrometeors are combined into a single array (e.g. the
c     cloud array contains both cloud water and ice, and the rain array
c     contains both rain and snow, as would be created by a "simple"
c     (non-mixed phase) bulk scheme.  iice will be determined by looking
c     for the presence of a snow mixing ratio array.
c     
      iice=0
      do i=varid,nvars
         nf_status = nf_inq_varname (ncid, varid, nf_varname)
         call handle_err(038,nf_status)
         if (nf_varname(1:5).eq.'QSNOW') then
            iice=1
            goto 357
         endif
      enddo
 357  continue
c
c   Write out model info to the Jim Bresch-inspired
c   ".minfo" file (only if this is the first data file)
c
      if (iwf.eq.1) then   ! start of minfo writing
c
      fname=argum(ncn)(1:iendc)//'.minfo'
      open (unit=58,file=fname,form='formatted',
     &      status='unknown')
c
      nf_att_text=' '
      nf_status = nf_get_att_text (ncid, nf_global,
     &   'TITLE', nf_att_text)
      call handle_err(039,nf_status)
      istart=index(nf_att_text,'WRF V')+4
      iend=index(nf_att_text(istart:),' ')-1+istart-1
      minfostring='Model info: '//nf_att_text(istart:iend)  ! First 18 chars.
c
      nf_status = nf_get_att_int (ncid, nf_global,
     &   'CU_PHYSICS', nf_att_int)
      call handle_err(040,nf_status)
      if (nf_att_int.eq.0) then
         minfostring(19:29)=' No Cumulus'
      elseif (nf_att_int.eq.1) then
         minfostring(19:29)=' Kain-F-Eta'
      elseif (nf_att_int.eq.2) then
         minfostring(19:29)=' Bet-Mil-Ja'
      elseif (nf_att_int.eq.99) then
         minfostring(19:29)=' Kain-Frsch'
      endif
c
      nf_status = nf_get_att_int (ncid, nf_global,
     &   'BL_PBL_PHYSICS', nf_att_int)
      call handle_err(041,nf_status)
      if (nf_att_int.eq.0) then
         minfostring(30:40)=' No PBL'
      elseif (nf_att_int.eq.1) then
         minfostring(30:40)=' MRF PBL'
      elseif (nf_att_int.eq.2) then
         minfostring(30:40)=' Mel-Yam-Ja'
      endif
c
      nf_status = nf_get_att_int (ncid, nf_global,
     &   'MP_PHYSICS', nf_att_int)
      call handle_err(042,nf_status)
      if (nf_att_int.eq.0) then
         minfostring(41:51)=' No microph'
      elseif (nf_att_int.eq.1) then
         minfostring(41:51)=' Kessler'
      elseif (nf_att_int.eq.2) then
         minfostring(41:51)=' Lin et al'
      elseif (nf_att_int.eq.3) then
         minfostring(41:51)=' NCEP simpl'
      elseif (nf_att_int.eq.4) then
         minfostring(41:51)=' NCEP mixed'
      elseif (nf_att_int.eq.5) then
         minfostring(41:51)=' Ferrier'
      elseif (nf_att_int.eq.99) then
         minfostring(41:51)=' Zhao-Carr'
      endif
c
      if (ds.ge.100000.) then
         write(minfostring(55:62),'(i3,'' km, '')') nint(.001*ds)
      elseif (ds.ge.10000.) then
         write(minfostring(55:62),'(i2,'' km, '')') nint(.001*ds)
      elseif (ds.ge.1000.) then
         write(minfostring(55:62),'(f3.1,'' km, '')') nint(.001*ds)
      else
         write(minfostring(55:62),'(i3,'' m, '')') nint(ds)
      endif
c
      write(minfostring(62:74),'(i3,'' levels, '')') mkzh
c
      nf_status = nf_get_att_real (ncid, nf_global,
     &   'DT', nf_att_real)
      call handle_err(043,nf_status)
      write(minfostring(75:81),'(i3,'' sec'')') nint(nf_att_real)
      write(58,'(a)') minfostring(1:81)
c
      endif     ! end of minfo writing
c
c   Set up stuff for RIP format data files.
c
      do i=1,32
         ihrip(i)=999999999
         rhrip(i)=9e9
         chrip(i)=' '
         chrip(i+32)=' '
      enddo
      chrip(1)=
     &'map projection (1: Lam. Conf., 2: Pol. Ster., 3: Mercator)'
      ihrip(1)=nproj
      chrip(2)=
     &'number of dot points in the y-direction (coarse domain)'
      ihrip(2)=miycors
      chrip(3)=
     &'number of dot points in the x-direction (coarse domain)'
      ihrip(3)=mjxcors
      chrip(4)=
     &'number of dot points in the y-direction (this domain)'
      ihrip(4)=miy
      chrip(5)=
     &'number of dot points in the x-direction (this domain)'
      ihrip(5)=mjx
      if (iexpanded.eq.1.and.iexpandedout.eq.0) then
         ihrip(2)=miy-2*ioffexp
         ihrip(3)=mjx-2*joffexp
         ihrip(4)=ihrip(2)
         ihrip(5)=ihrip(3)
      endif
      chrip(6)=
     &'number of dimensions of this variable (2 or 3)'
      ihrip(6)=999    ! this is set separately for each variable
      chrip(7)=
     &'grid of this variable (1: cross point dom., 0: dot point dom.)'
      ihrip(7)=999    ! this is set separately for each variable
ccc      chrip(8)=
ccc     &'vertical coordinate (0: hydrostatic sigma, 1: nonhyd. sigma)'
ccc      ihrip(8)=inhyd
ccc      chrip(9)=
ccc     &'number of half sigma levels'
      chrip(9)=
     &'number of vertical levels in the data'
      ihrip(9)=mkzh
      chrip(10)=
     &'mdateb: YYMMDDHH (truncated hour) of hour-0 for this dataset'
      ihrip(10)=mdateb
      chrip(11)=
     &'mdate: YYMMDDHH (truncated hour) of this time'
      ihrip(11)=99999999    ! this is set separately for each time
      chrip(12)=
     &'ice physics (1: sep. arrays for ice fields, 0: no sep. arrays)'
      ihrip(12)=iice
ccc      chrip(13)=
ccc     &'Program #: 1:TER. 2:DG/RG. 3:RAW. 5:INT. 6:MOD. 11:MOD.(MM5V3)'
      chrip(13)=
     &'ver. coord. type: <or=3: hgt. or prs.; >or=4: terrain-following'
      ihrip(13)=iplevdata
      chrip(14)=
     &'landuse dataset (1: old, 13-cat; 2: USGS, 24-cat; 3: SiB, 16 )'
      ihrip(14)=ilandset
c
      ijmp=32
      chrip(ijmp+1)=
     &'first true latitude (deg.)'
      rhrip(1)=truelat1
      chrip(ijmp+2)=
     &'second true latitude (deg.)'
      rhrip(2)=truelat2
      chrip(ijmp+3)=
     &'central latitude of coarse domain (deg.)'
      rhrip(3)=xlatc
      chrip(ijmp+4)=
     &'central longitude of coarse domain(deg.)'
      rhrip(4)=xlonc
      chrip(ijmp+5)=
     &'grid distance of coarse domain (km)'
      rhrip(5)=dskmc
      chrip(ijmp+6)=
     &'grid distance of this domain (km)'
      rhrip(6)=dskm
      chrip(ijmp+7)=
     &'coarse dom. y-position of lower left corner of this domain'
      rhrip(7)=9e9       ! this is set separately for each time
      chrip(ijmp+8)=
     &'coarse dom. x-position of lower left corner of this domain'
      rhrip(8)=9e9       ! this is set separately for each time
ccc      chrip(ijmp+9)=
ccc     &'pressure level (hPa) of the model top'
ccc      rhrip(9)=ptop
ccc      chrip(ijmp+10)=
ccc     &'reference sea-level pressure (Pa)'
ccc      rhrip(10)=refslp
ccc      chrip(ijmp+11)=
ccc     &'reference sea-level temperature (K)'
ccc      rhrip(11)=refslt
ccc      chrip(ijmp+12)=
ccc     &'reference lapse rate (dT/d(ln(p)), K)'
ccc      rhrip(12)=reflaps
      chrip(ijmp+13)=
     &'rhourb: diff (in h) between exact time and mdate of hour-0'
      rhrip(13)=rhourb
      chrip(ijmp+14)=
     &'rhour: diff (in h) between exact time and mdate of this data'
      rhrip(14)=9e9       ! this is set separately for each time
      chrip(ijmp+15)=
     &'xtime: exact time of this data relative to exact hour-0 (in h)'
      rhrip(15)=9e9       ! this is set separately for each time
ccc      chrip(ijmp+16)=
ccc     &'reference stratospheric constant temperature (K)'
ccc      rhrip(16)=refstratt
c
c   Initialize "max" time level. This feature causes RIPDP to
c      ignore standard data which is out of chronological order.
c
      xtimemax=-1.
c
c=================================================================c
      do iwt=1,nwrftimes_this_file     ! Time loop
c=================================================================c
c
c   Note the different time specifications:
c
c   mdateb: This refers to the truncated integer hour of the
c     beginning of the model run for model output, or of the
c     starting time for analysis data sets.  It is an 8-digit
c     integer specified as YYMMDDHH.
c
c   mhourb: This integer variable refers to the same time as
c     mdateb, but instead of the YYMMDDHH format, it is specified
c     as the number of hours since 00 UTC 1 January 1 AD, according
c     to the Gregorian calendar with full leap year specification
c     (leap year every 4 years except on century years not dvisible
c     by 400).
c
c   rhourb: This is a real number specifying the fraction of an
c     hour that the exact start time of the model forecast or
c     analysis data set exceeds the truncated integer hour of the
c     start time (mdateb/mhourb).  Typically, mesoscale models
c     are initialized precisely on the hour, so rhourb is
c     typically 0.00.  However, rhourb could, in
c     principal, be in the range 0.0 < or = rhourb < 1.0.
c
c   xtime: This is a real number referring to this particular
c     data time, and is specified as the exact number of hours
c     since the exact model start time or first analysis time,
c     i.e. since the time mhourb+rhourb.
c
c   mdate: This refers to the truncated integer hour of this
c     particular data time.  It is an 8-digit integer specified as
c     YYMMDDHH.
c
c   mhour: This integer variable refers to the same time as
c     mdate, but instead of the YYMMDDHH format, it is specified as
c     the number of hours since 00 UTC 1 January 1 AD (similar to
c     mhourb).
c
c   rhour: This is a real number specifying the fraction of an
c     hour that the exact time of this particular data time exceeds
c     the truncated integer hour of this particular data time
c     (mdate/mhour). rhour is in the range 0.0 < or = rhour < 1.0.
c
c   Set count variables that control time index in netcdf arrays
c
      nf_ustart3(4)=iwt
      nf_vstart3(4)=iwt
      nf_tstart3(4)=iwt
      nf_wstart3(4)=iwt
      nf_tstart2(3)=iwt
c
c   Create mdate, rhour, xtime.
c
      read(wrftimes(iwt),'(1x,6(1x,i2))')iyr,imo,idy,ihr,imn,isc
      mdate=1000000*iyr+10000*imo+
     &   100*idy+ihr
      rhour=imn/60.+isc/3600.
      call mconvert(mdate,mhour,1,1940)
      xtime=float(mhour-mhourb)+rhour-rhourb
c
c     Unfortunately, there is currently no information in the WRF output
c     netcdf attributes to indicate yicorn and xjcorn, the y and x
c     locations of the lower left corner of the current domain with
c     respect to the centered coarse domain (i.e. the coarse domain
c     centered on lat=xlatc,lon=xlonc).  It can only be ascertained from
c     the lat/lon data arrays.  So first, we must get latitude (XLAT)
c     and longitude (XLONG) arrays.
c
      nf_status = nf_inq_varid (ncid, 'XLAT', varid)
      call handle_err(033,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(034,nf_status)
      do j=1,mjx-1
      do i=1,miy-1
         xlat(i,j)=nf_tarr2(j,i)
      enddo
      enddo
      nf_status = nf_inq_varid (ncid, 'XLONG', varid)
      call handle_err(035,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(036,nf_status)
      do j=1,mjx-1
      do i=1,miy-1
         xlon(i,j)=nf_tarr2(j,i)
      enddo
      enddo
c
c     Use given lat/lon at cross grid point (y=1,x=1) to determine where
c     that point is in the centered coarse domain.
c
      call maptform(yicorncross,xjcorncross,xlat(1,1),xlon(1,1),-1,0.)
c
c   Set yicorn, xjcorn
c
      yicorn=yicorncross-.5/refrat
      xjcorn=xjcorncross-.5/refrat
c
c  Set RIP header values that are time-dependent
c
      ihrip(11)=mdate
      rhrip(7)=yicorn
      rhrip(8)=xjcorn
      rhrip(14)=rhour
      rhrip(15)=xtime
c
      secondspast=rhour*3600.
      print*
      print*,'****  Reading model output at'
      print*,'      forecast time=',xtime
      if (secondspast.lt..2) then
         write(6,926) mdate
      else
         write(6,927) mdate,secondspast
      endif
 926  format('        (YYMMDDHH = ',i8.8,')')
 927  format('        (YYMMDDHH = ',i8.8,' plus ',f12.5,' seconds)')
c
c   See if the time we just encountered is beyond the latest time
c   requested.  If so, then jump out of time loop here.
c
      if (xtime.gt.ptusemax+tacch) then
         print*,'   But this time is beyond the latest requested time.'
         print*,'   RIPDP is now stopping.'
         goto 1000
      endif
c
c   See if this is a time that is out of chron. order
c
      iskipit=0
      if (xtime.le.xtimemax) iskipit=1  ! this means out of ch. order
      xtimemax=max(xtimemax,xtime)
c
c   See if this is a time that the user doesn't want.
c   Note: if no ptimes or iptimes were specified (nptuse=0), then
c   it is assumed the user wants all encountered times processed
c   (unless they are out of chronological order).
c
      if (iskipit.eq.0.and.nptuse.gt.0) then
         do i=1,nptuse
            if (abs(xtime-ptuse(i)).le.tacch) goto 40
         enddo
         iskipit=2
 40      continue
      endif
c      
      if (iskipit.gt.0) then
c
         if (iskipit.eq.1) then
            print*,'   But this time is chronologically backward,',
     &        ' so we will skip it.'
         else
            print*,'   But you do not want to process this time,',
     &        ' so we will skip it.'
         endif
c
      endif
c
c   Create as much of the data file name as we know at this point.
c
      cxtime=' '
      write(cxtime,'(f10.5)')xtime
      if (cxtime(1:1).eq.' ') cxtime(1:1)='0'
      if (cxtime(2:2).eq.' ') cxtime(2:2)='0'
      if (cxtime(3:3).eq.' ') cxtime(3:3)='0'
      if (cxtime(4:4).eq.' ') cxtime(4:4)='0'
      fname=argum(ncn)(1:iendc)//'_'//cxtime//'_'
      nxtavl=nxtavl+1
      cxtimeavl(nxtavl)=cxtime
c
c   Obtain variables specifically sought by RIPDP, and write them out.
c
      print*,'Processing basic variables.'
c
c   Get U, transfer to "dot-point" grid, and write out.
c
      nf_status = nf_inq_varid (ncid, 'U', varid)
      call handle_err(044,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_ustart3,
     &   nf_ucount3, nf_uarr3)
      call handle_err(045,nf_status)
      do k=1,mkzh
         do j=1,mjx
         do i=1,miy-1
            uuu(i,j,k)=nf_uarr3(j,i,mkzh-k+1)
         enddo
         enddo
         call utodot(uuu(1,1,k),miy,mjx)
      enddo
      vardesc='Horizontal wind (y-comp.), m/s'
      plchun='m s~S~-1~N~'
      call writefile(uuu,'uuu       ',3,0,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
c   Get V, transfer to "dot-point" grid, and write out.
c
      nf_status = nf_inq_varid (ncid, 'V', varid)
      call handle_err(046,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_vstart3,
     &   nf_vcount3, nf_varr3)
      call handle_err(047,nf_status)
      do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy
            vvv(i,j,k)=nf_varr3(j,i,mkzh-k+1)
         enddo
         enddo
         call vtodot(vvv(1,1,k),miy,mjx)
      enddo
      vardesc='Horizontal wind (y-comp.), m/s'
      plchun='m s~S~-1~N~'
      call writefile(vvv,'vvv       ',3,0,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
c   Get pressure perturbation (P) and basic-state pressure (PB), combine,
c   change to hPa (from Pa), and write out.
c
      nf_status = nf_inq_varid (ncid, 'P', varid)
      call handle_err(048,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &   nf_tcount3, nf_tarr3)
      call handle_err(049,nf_status)
      do k=1,mkzh
      do j=1,mjx-1
      do i=1,miy-1
         prs(i,j,k)=nf_tarr3(j,i,mkzh-k+1)
      enddo
      enddo
      enddo
      nf_status = nf_inq_varid (ncid, 'PB', varid)
      call handle_err(050,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &   nf_tcount3, nf_tarr3)
      call handle_err(051,nf_status)
      do k=1,mkzh
      do j=1,mjx-1
      do i=1,miy-1
         prs(i,j,k)=.01*(prs(i,j,k)+nf_tarr3(j,i,mkzh-k+1))
      enddo
      enddo
      enddo
      vardesc='Pressure, hPa'
      plchun='hPa'
      call writefile(prs,'prs       ',3,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
c   If available, get water vapor mixing ratio (QVAPOR), convert it to g/kg,
c   and write out.  Otherwise, fill qvp array with 0s but don't write out.
c
      nf_status = nf_inq_varid (ncid, 'QVAPOR', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find QVAPOR.'
         call fillarray(qvp,miy*mjx*mkzh,0.)
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(052,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            qvp(i,j,k)=nf_tarr3(j,i,mkzh-k+1)*1000.  ! kg/kg -> g/kg
         enddo
         enddo
         enddo
         vardesc='Water vapor mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile(qvp,'qvp       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,
     &      iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
c
c   Get theta (T), convert it to temperature, and write out.
c
      nf_status = nf_inq_varid (ncid, 'T', varid)
      call handle_err(053,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &   nf_tcount3, nf_tarr3)
      call handle_err(054,nf_status)
      do k=1,mkzh
      do j=1,mjx-1
      do i=1,miy-1
         gammam=gamma*(1.+gammamd*.001*qvp(i,j,k))
         tmk(i,j,k)=(nf_tarr3(j,i,mkzh-k+1)+300.)*
     &      (prs(i,j,k)/1000.)**gammam
      enddo
      enddo
      enddo
      vardesc='Temperature, K'
      plchun='K'
      call writefile(tmk,'tmk       ',3,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
c     Geopotential perturbation and basic state, as well as vertical
c     velocity, are defined at vertical staggerred levels, i.e. the
c     "w" levels, similar to so-called "full sigma" levels in the MM5
c     coordinate system.  For the geopotential, since the WRF vertical
c     coordinate is mass based, it would be best to convert the
c     geopotential to an exponential height variable, then linearly
c     interpolate in WRF vert. coord. ("eta") from "w levels" to "mass
c     levels", and then convert back to geopotential height.
c
c   First read geopotential base state (PHB) into the netcdf scratch array
c   for "w-level" variables, and transfer to RIP scratch array for
c   "w-level" variables.
c
      nf_status = nf_inq_varid (ncid, 'PHB', varid)
      call handle_err(055,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_wstart3,
     &   nf_wcount3, nf_warr3)
      call handle_err(056,nf_status)
      do k=1,mkzh+1
      do j=1,mjx-1
      do i=1,miy-1
         scr3wlev(i,j,k)=nf_warr3(j,i,mkzh+1-k+1)
      enddo
      enddo
      enddo
c
c   Now read geopotential perturbation (PH) into the netcdf scratch array
c   for "w-level" variables, and add to RIP scratch array for
c   "w-level" variables.  Convert geopotential to exp(-z/H), where H
c   is a scale height.
c
      nf_status = nf_inq_varid (ncid, 'PH', varid)
      call handle_err(057,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_wstart3,
     &   nf_wcount3, nf_warr3)
      call handle_err(058,nf_status)
      do k=1,mkzh+1
      do j=1,mjx-1
      do i=1,miy-1
         scr3wlev(i,j,k)=exp(-(scr3wlev(i,j,k)+
     &      nf_warr3(j,i,mkzh+1-k+1))/(grav*sclht))
      enddo
      enddo
      enddo
c
c   Now interpolate exp(-z/H) to "mass levels", convert to height, and
c   write out.
c
      do j=1,mjx-1
      do i=1,miy-1
      do k=1,mkzh
         ght(i,j,k)=znfac(k)*scr3wlev(i,j,k+1)+
     &              (1.-znfac(k))*scr3wlev(i,j,k)
         ght(i,j,k)=-sclht*log(ght(i,j,k))
      enddo
      enddo
      enddo
      vardesc='Geopotential height, m'
      plchun='m'
      call writefile(ght,'ght       ',3,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
c   Now get vertical velocity (convert from m/s to cm/s)
c
      nf_status = nf_inq_varid (ncid, 'W', varid)
      call handle_err(059,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_wstart3,
     &   nf_wcount3, nf_warr3)
      call handle_err(060,nf_status)
      do k=1,mkzh+1
      do j=1,mjx-1
      do i=1,miy-1
         scr3wlev(i,j,k)=100.*nf_warr3(j,i,mkzh+1-k+1)
      enddo
      enddo
      enddo
c
c   Interpolate vert. vel. to "mass levels" and write out
c
      do j=1,mjx-1
      do i=1,miy-1
      do k=1,mkzh
         www(i,j,k)=znfac(k)*scr3wlev(i,j,k+1)+
     &              (1.-znfac(k))*scr3wlev(i,j,k)
      enddo
      enddo
      enddo
      vardesc='Vertical velocity, cm/s'
      plchun='cm s~S~-1~N~'
      call writefile(www,'www       ',3,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
c   Get terrain height (HGT), and write out.
c
      nf_status = nf_inq_varid (ncid, 'HGT', varid)
      call handle_err(066,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(067,nf_status)
      do j=1,mjx-1
      do i=1,miy-1
         ter(i,j)=nf_tarr2(j,i)
      enddo
      enddo
      plchun='m'
      vardesc='Terrain height AMSL, m'
      call writefile(ter,'ter       ',2,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
c   Calculate surface pressure using altimeter equation.
c
      do j=1,mjx-1
      do i=1,miy-1
         tv=virtual(tmk(i,j,mkzh),.001*qvp(i,j,mkzh))
         sfp(i,j)=prs(i,j,mkzh)*(tv/(tv+ussalr*
     *      (ght(i,j,mkzh)-ter(i,j))))**(-grav/(rgas*ussalr))
      enddo
      enddo
      vardesc='Surface pressure, hPa'
      plchun='hPa'
      call writefile(sfp,'sfp       ',2,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
c   Get mapscale factor on cross points (MAPFAC_M), and write out.
c
      nf_status = nf_inq_varid (ncid, 'MAPFAC_M', varid)
      call handle_err(068,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(069,nf_status)
      do j=1,mjx-1
      do i=1,miy-1
         xmap(i,j)=nf_tarr2(j,i)
      enddo
      enddo
      vardesc='Map factor on cross points'
      plchun='none'
      call writefile(xmap,'xmap      ',2,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
c   Transfer xmap to dot-point grid (dmap) and write out
c
      do j=1,mjx-1
      do i=1,miy-1
         dmap(i,j)=xmap(i,j)
      enddo
      enddo
      call xtodot(xmap,miy,mjx)
      vardesc='Map factor on dot points'
      plchun='none'
      call writefile(dmap,'dmap      ',2,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
c   Get coriolis parameter on cross points (F), and write out.
c
      nf_status = nf_inq_varid (ncid, 'F', varid)
      call handle_err(070,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(071,nf_status)
      do j=1,mjx-1
      do i=1,miy-1
         cor(i,j)=nf_tarr2(j,i)
      enddo
      enddo
      vardesc='Coriolis parameter, per s'
      plchun='s~S~-1~N~'
      call writefile(cor,'cor       ',2,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
c   Already have latitude (XLAT) and longitude (XLONG).  Write them out.
c
      vardesc='Latitude, degrees'
      plchun='degrees'
      call writefile(xlat,'xlat      ',2,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      vardesc='Longitude, degrees'
      plchun='degrees'
      call writefile(xlon,'xlon      ',2,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
c   Get total accumulated cumulus preciitation (RAINC), and write out.
c
      nf_status = nf_inq_varid (ncid, 'RAINC', varid)
      call handle_err(072,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(073,nf_status)
      do j=1,mjx-1
      do i=1,miy-1
         rtc(i,j)=nf_tarr2(j,i)
      enddo
      enddo
      vardesc='Cumulus precip. since h 0, mm'
      plchun='mm'
      call writefile(rtc,'rtc       ',2,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
c   Get total accumulated explicit (grid-resolved) preciitation (RAINNC),
c      and write out.
c
      nf_status = nf_inq_varid (ncid, 'RAINNC', varid)
      call handle_err(074,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(075,nf_status)
      do j=1,mjx-1
      do i=1,miy-1
         rte(i,j)=nf_tarr2(j,i)
      enddo
      enddo
      vardesc='Explicit precip. since h 0, mm'
      plchun='mm'
      call writefile(rte,'rte       ',2,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
c   Get ground temperature (and sea-surface temperature over water) (TSK),
c      and write out.
c
      nf_status = nf_inq_varid (ncid, 'TSK', varid)
      call handle_err(076,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(077,nf_status)
      do j=1,mjx-1
      do i=1,miy-1
         tgk(i,j)=nf_tarr2(j,i)
      enddo
      enddo
      vardesc='Ground/sea-surface temperature, K'
      plchun='K'
      call writefile(tgk,'tgk       ',2,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
c   Get land use (LU_INDEX), and write out.
c
      nf_status = nf_inq_varid (ncid, 'LU_INDEX', varid)
      call handle_err(078,nf_status)
      iprocvarid(varid)=1
      nf_status = nf_get_vara_real (ncid, varid, nf_tstart2,
     &   nf_tcount2, nf_tarr2)
      call handle_err(079,nf_status)
      do j=1,mjx-1
      do i=1,miy-1
         xlus(i,j)=nf_tarr2(j,i)
      enddo
      enddo
      vardesc='Land use category'
      plchun='none'
      call writefile(xlus,'xlus      ',2,1,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,
     &   iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
      print*,'Checking for hydrometeor variables.'
c
c   If available, get cloud water mixing ratio (QCLOUD), convert it to g/kg,
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'QCLOUD', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find QCLOUD.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(061,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3(i,j,k)=nf_tarr3(j,i,mkzh-k+1)*1000.  ! kg/kg -> g/kg
         enddo
         enddo
         enddo
         vardesc='Cloud water mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile(scr3,'qcw       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,
     &      iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
c
c   If available, get rain water mixing ratio (QRAIN), convert it to g/kg,
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'QRAIN', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find QRAIN.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(062,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3(i,j,k)=nf_tarr3(j,i,mkzh-k+1)*1000.  ! kg/kg -> g/kg
         enddo
         enddo
         enddo
         vardesc='Rain water mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile(scr3,'qra       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,
     &      iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
c
c   If available, get cloud ice mixing ratio (QICE), convert it to g/kg,
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'QICE', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find QICE.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(063,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3(i,j,k)=nf_tarr3(j,i,mkzh-k+1)*1000.  ! kg/kg -> g/kg
         enddo
         enddo
         enddo
         vardesc='Cloud ice mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile(scr3,'qci       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,
     &      iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
c
c   If available, get snow mixing ratio (QSNOW), convert it to g/kg,
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'QSNOW', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find QSNOW.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(064,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3(i,j,k)=nf_tarr3(j,i,mkzh-k+1)*1000.  ! kg/kg -> g/kg
         enddo
         enddo
         enddo
         vardesc='Snow mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile(scr3,'qsn       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,
     &      iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
c
c   If available, get graupel mixing ratio (QGRAUP), convert it to g/kg,
c   and write out.
c
      nf_status = nf_inq_varid (ncid, 'QGRAUP', varid)
      if (nf_status .ne. nf_noerr) then
         print*,'   Did not find QGRAUP.'
      else
         iprocvarid(varid)=1
         nf_status = nf_get_vara_real (ncid, varid, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(065,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3(i,j,k)=nf_tarr3(j,i,mkzh-k+1)*1000.  ! kg/kg -> g/kg
         enddo
         enddo
         enddo
         vardesc='Graupel mixing ratio, g/kg'
         plchun='g kg~S~-1~N~'
         call writefile(scr3,'qgr       ',3,1,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,
     &      iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      endif
c
c     Loop through other variables not specifically sought by RIPDP.
c
      if (ibasic.eq.1.and.nretain.eq.0) goto 230
c
      print*,'Processing other variables in output file.'
c
      do ivar=1,nvars
c
      nf_varname=' '
      nf_status = nf_inq_varname (ncid, ivar, nf_varname)
      call handle_err(080,nf_status)
c
c     First, jump to end of loop if variable was already processed in some
c     way shape or form.
c
      if (iprocvarid(ivar).eq.1) then
         goto 229
      endif
c
c     Next, jump to end of loop if "all" was specified on the command line
c     but the variable name is in the discard list.
c
      if (ibasic.eq.0) then
         do idis=1,ndiscard
            if (nf_varname.eq.discard(idis)) then
               goto 229
            endif
         enddo
      endif
c
c     Next, jump to end of loop if "basic" was specified on the command line
c     but the variable name is NOT in the retain list.
c
      if (ibasic.eq.1) then
         icount=0
         do iret=1,nretain
            if (nf_varname.eq.retain(iret)) icount=icount+1
         enddo
         if (icount.eq.0) goto 229
      endif
c
c     Next, check for certain combinations of variable dimensions
c     that RIP can make use of (i.e., 3D cross point array, 2D cross
c     point array, or soil-layer cross point array).
c
      nf_status = nf_inq_varndims (ncid, ivar, ndims)
      call handle_err(081,nf_status)
      nf_status = nf_inq_vardimid (ncid, ivar, vardimids)
      call handle_err(082,nf_status)
c
      if (ndims.eq.4.and.vardimids(4).eq.dimid_tm.and.
     &    vardimids(3).eq.dimid_bt.and.vardimids(2).eq.dimid_sn.and.
     &    vardimids(1).eq.dimid_we) then
         itype=1  ! 3D cross point array
      elseif (ndims.eq.3.and.vardimids(3).eq.dimid_tm.and.
     &    vardimids(2).eq.dimid_sn.and.vardimids(1).eq.dimid_we) then
         itype=2  ! 2D cross point array
      elseif (ndims.eq.4.and.vardimids(4).eq.dimid_tm.and.
     &    vardimids(3).eq.dimid_sls.and.vardimids(2).eq.dimid_sn.and.
     &    vardimids(1).eq.dimid_we) then
         itype=3  ! soil-layer cross point array
      else
         itype=0  ! none of the above
      endif
c
c     If variable does not have a combination of dimensions
c     that RIP can make use of, skip it.
c
      if (itype.eq.0) goto 229
c
      varname=nf_varname
      inname=0
      do ic=10,1,-1
         if (inname.eq.0) then
            if (varname(ic:ic).ne.' ') inname=1
         else
            if (varname(ic:ic).eq.' ') varname(ic:ic)='_'
         endif
      enddo
      nf_att_text=' '
      nf_status = nf_get_att_text (ncid, ivar,
     &   'units', nf_att_text)
      call handle_err(083,nf_status)
      plchun=nf_att_text
      nf_att_text=' '
      nf_status = nf_get_att_text (ncid, ivar,
     &   'description', nf_att_text)
      call handle_err(084,nf_status)
      nf_status = nf_inq_attlen (ncid, ivar,
     &   'description', nf_att_len)
      call handle_err(085,nf_status)
      vardesc=nf_att_text(1:nf_att_len)//', '//plchun
      icd=1
      if (itype.eq.1) then
         nf_status = nf_get_vara_real (ncid, ivar, nf_tstart3,
     &      nf_tcount3, nf_tarr3)
         call handle_err(086,nf_status)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3(i,j,k)=nf_tarr3(j,i,mkzh-k+1)
         enddo
         enddo
         enddo
         call writefile(scr3,varname,3,icd,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,
     &      iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      elseif (itype.eq.2) then
         nf_status = nf_get_vara_real (ncid, ivar, nf_tstart2,
     &      nf_tcount2, nf_tarr2)
         call handle_err(087,nf_status)
         do j=1,mjx-1
         do i=1,miy-1
            scr2(i,j)=nf_tarr2(j,i)
         enddo
         enddo
         call writefile(scr2,varname,2,icd,vardesc,plchun,
     &      fname,iendf1,ihrip,rhrip,chrip,
     &      iexpanded,
     &      iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
      elseif (itype.eq.3) then
         nf_status = nf_get_vara_real (ncid, ivar, nf_tstart3,
     &      nf_tcount3s, nf_tarr3)
         call handle_err(088,nf_status)
         isp=index(varname,' ')
         if (isp.eq.0.or.isp.eq.10) isp=9
         do k=1,nsoil
            do j=1,mjx-1
            do i=1,miy-1
               scr2(i,j)=nf_tarr3(j,i,k)
            enddo
            enddo
            write(varname(isp:isp+1),'(i2.2)') k
            vardesc=nf_att_text(1:nf_att_len)//
     &         ', layer '//varname(isp:isp+1)//', '//plchun
            call writefile(scr2,varname,2,icd,vardesc,plchun,
     &         fname,iendf1,ihrip,rhrip,chrip,
     &         iexpanded,
     &         iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
         enddo
      endif
c
 229  continue  ! point to jump to if variable is being skipped
c
      enddo     ! end of "other variables (do ivar=1,nvars)" loop
c
 230  continue  ! jump past loop
c
      print*,'Finished reading data for this time.'
c
c   Write out the available xtimes in the ".xtimes" file.  Do this after
c   every time, so that if the user kills the program in mid-run, he'll
c   still have a useful .xtimes file
c
      fname=argum(ncn)(1:iendc)//'.xtimes'
      open (unit=57,file=fname,form='formatted',status='unknown')
      write(57,*) nxtavl
      do i=1,nxtavl
         write(57,'(a10)') cxtimeavl(i)
      enddo
      close (57)
c
c=================================================================c
      enddo      ! End of time loop.
c=================================================================c
c
c=================================================================c
      enddo      ! End of file loop.
c=================================================================c
c
 1000 continue   ! destination for jumping out of both loops
c
      print*
      print*,'===================================='
      print*,' We''re outta here like Vladimir !! '
      print*,'===================================='
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine fillarray(array,ndim,val)
      dimension array(ndim)
      do i=1,ndim
         array(i)=val
      enddo
c
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine mconvert(mdate,mhour,idir,nsplityear)
c
c   mdate: an 8-digit integer specification for a date,
c      given as yymmddhh
c   mhour: an integer specificying the number of hours since
c      00 UTC 1 January 1 AD.
c
c   This routine converts an mdate to an mhour if idir=1, or vice versa
c   if idir=-1.
c
c   If idir=1, how do we know what century mdate refers to?  You
c   provide a year, called "nsplityear", denoted "aabb".  If mdate
c   is denoted "yymmddhh", then if yy >or= bb, the century is
c   assumed to be aa.  Otherwise it is assumed to be the century
c   after aa, or aa+1.
c
c   Leap year definition: every fourth year has a 29th day in February,
c      with the exception of century years not divisible by 400.
c
      dimension ndaypmo(12)
      integer yy,mm,dd,hh,aa,bb
      data ndaypmo /31,28,31,30,31,30,31,31,30,31,30,31/
c
      if (idir.eq.1) then
c
      yy=mdate/1000000
      bb=mod(nsplityear,100)
      aa=nsplityear-bb
      iyear=aa+yy
      if (yy.lt.bb) iyear=iyear+100
      iyearp=iyear-1
      idayp = iyearp*365 + iyearp/4 - iyearp/100 +iyearp/400
      mm=mod(mdate,1000000)/10000
      imonthp=mm-1
      if ((mod(iyear,4).eq.0.and.mod(iyear,100).ne.0).or.
     &    mod(iyear,400).eq.0)
     &   ndaypmo(2)=29
      do i=1,imonthp
         idayp=idayp+ndaypmo(i)
      enddo
      ndaypmo(2)=28
      dd=mod(mdate,10000)/100
      idayp=idayp+dd-1
      hh=mod(mdate,100)
      mhour=24*idayp+hh
c
      else
c
      nhour=mhour
c
c   Get an estimate of iyear that is guaranteed to be close to but
c   less than the current year
c
      iyear = max(0,nhour-48)*1.14079e-4
      ihour=24*(iyear*365+iyear/4-iyear/100+iyear/400)
 10   iyear=iyear+1
      ihourp=ihour
      ihour = 24*(iyear*365 + iyear/4 - iyear/100 +iyear/400)
      if (ihour.le.nhour) goto 10
      nhour=nhour-ihourp
      if ((mod(iyear,4).eq.0.and.mod(iyear,100).ne.0).or.
     &    mod(iyear,400).eq.0)
     &   ndaypmo(2)=29
      imo=0
      ihour=0
 20   imo=imo+1
      ihourp=ihour
      ihour=ihour+24*ndaypmo(imo)
      if (ihour.le.nhour) goto 20
      nhour=nhour-ihourp
      ndaypmo(2)=28
      iday = nhour/24 + 1
      ihour=mod(nhour,24)
      mdate=mod(iyear,100)*1000000+imo*10000+iday*100+ihour
c
      endif
c
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine writefile (var,varname,ndim,icd,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
      dimension var(miy,mjx,1+(mkzh-1)*(ndim-2))
      character varname*10,fname*256
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
c ensure data do not exceed ieee limits
c
      kk = 1+(mkzh-1)*(ndim-2)
      do k=1,kk
         do j=1,mjx-icd
         do i=1,miy-icd
           if(var(i,j,k) .gt. 3.4e+38) var(i,j,k)= 3.4e+38
           if(var(i,j,k) .lt.-3.4e+38) var(i,j,k)=-3.4e+38
           if(abs(var(i,j,k)).lt. 1.2e-38) var(i,j,k)= 0.
         enddo
         enddo
         if (icd.eq.1) then
            do i=1,miy
               var(i,mjx,k)=0.
            enddo
            do j=1,mjx-1
               var(miy,j,k)=0.
            enddo
         endif
      enddo
      iendv=index(varname,' ')-1
      if (iendv.eq.-1) iendc=10
      fname(iendf1+1:)=varname
      open(unit=65,file=fname,form='unformatted',status='unknown')
      ihrip(6)=ndim ! number of dimensions of this variable (2 or 3)
      ihrip(7)=icd  ! grid of this var. (1: cross point, 0: dot point)
ccc      write(65) vardesc,plchun,ihrip,rhrip,chrip,fullsigma,halfsigma
      write(65) vardesc,plchun,ihrip,rhrip,chrip
      if (iexpanded.eq.1.and.iexpandedout.eq.0) then
         write(65) (((var(i,j,k),i=1+ioffexp,miy-ioffexp),
     &      j=1+joffexp,mjx-joffexp),k=1,1+(mkzh-1)*(ndim-2))
      else
         write(65) var
      endif
      close (65)
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      function virtual(temp,ratmix)
c
c   This function returns virtual temperature in K, given temperature
c      in K and mixing ratio in kg/kg.
c
      parameter (eps=0.622)
c
      virtual=temp*(eps+ratmix)/(eps*(1.+ratmix))
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine handle_err(rmarker,nf_status)
      include "netcdf.inc"
      integer nf_status
      real rmarker
      if (nf_status .ne. nf_noerr) then
         write(*,*)  'NetCDF error in ripdp_wrf.  Marker = ',rmarker
         write(*,*)  '  ',nf_strerror(nf_status)
         stop 'stopped'
      endif
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine xtodot(slab,maxiy,maxjx)
c
c     This routine converts data that is on the B-grid mass grid (known
c     in MM5 lingo as "cross points") to the B-grid velocity staggered
c     grid (known in MM5 lingo as "dot points")
c
      dimension slab(maxiy,maxjx),bot(1000),rleft(1000)
c
c   Extrapolate out to top and bottom edges.
c
      do 200 j=2,maxjx-1
         bot(j)=(3.*(slab(1,j-1)+slab(1,j))-
     &      (slab(2,j-1)+slab(2,j)))/4.
         slab(maxiy,j)=(3.*(slab(maxiy-1,j-1)+slab(maxiy-1,j))-
     &      (slab(maxiy-2,j-1)+slab(maxiy-2,j)))/4.
  200 continue
c
c   Extrapolate out to left and right edges.
c
      do 300 i=2,maxiy-1
         rleft(i)=(3.*(slab(i-1,1)+slab(i,1))-
     &      (slab(i-1,2)+slab(i,2)))/4.
         slab(i,maxjx)=(3.*(slab(i-1,maxjx-1)+slab(i,maxjx-1))-
     &      (slab(i-1,maxjx-2)+slab(i,maxjx-2)))/4.
  300 continue
c
c   Extrapolate out to corners.
c
      rleft(1)=(3.*slab(1,1)-slab(2,2))/2.
      rleft(maxiy)=(3.*slab(maxiy-1,1)-slab(maxiy-2,2))/2.
      bot(maxjx)=(3.*slab(1,maxjx-1)-slab(2,maxjx-2))/2.
      slab(maxiy,maxjx)=(3.*slab(maxiy-1,maxjx-1)-
     &   slab(maxiy-2,maxjx-2))/2.
c
c   Interpolate in the interior.
c
      do 100 j=maxjx-1,2,-1
      do 100 i=maxiy-1,2,-1
         slab(i,j)=.25*(slab(i-1,j-1)+slab(i,j-1)+slab(i-1,j)+
     &      slab(i,j))
  100    continue
c
c   Put "bot" and "rleft" values into slab.
c
      do j=2,maxjx
         slab(1,j)=bot(j)
      enddo
      do i=1,maxiy
         slab(i,1)=rleft(i)
      enddo
c
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine utodot(slab,maxiy,maxjx)
c
c   This routine converts data that is on the C-grid u-velocity
c   staggered grid to the B-grid velocity staggered grid (known in
c   MM5 lingo as "dot points") 
c
      dimension slab(maxiy,maxjx),bot(2000)
c
c   Extrapolate out to top and bottom edges.
c
      do j=1,maxjx
         bot(j)=(3.*slab(1,j)-slab(2,j))/2.
         slab(maxiy,j)=(3.*slab(maxiy-1,j)-slab(maxiy-2,j))/2.
      enddo
c
c   Interpolate in the interior.
c
      do j=maxjx,1,-1
      do i=maxiy-1,2,-1
         slab(i,j)=.5*(slab(i-1,j)+slab(i,j))
      enddo
      enddo
c
c   Put "bot" values into slab.
c
      do j=1,maxjx
         slab(1,j)=bot(j)
      enddo
c
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine vtodot(slab,maxiy,maxjx)
c
c   This routine converts data that is on the C-grid v-velocity
c   staggered grid to the B-grid velocity staggered grid (known in
c   MM5 lingo as "dot points") 
c
      dimension slab(maxiy,maxjx),rleft(2000)
c
c   Extrapolate out to left and right edges.
c
      do i=1,maxiy
         rleft(i)=(3.*slab(i,1)-slab(i,2))/2.
         slab(i,maxjx)=(3.*slab(i,maxjx-1)-slab(i,maxjx-2))/2.
      enddo
c
c   Interpolate in the interior.
c
      do j=maxjx-1,2,-1
      do i=maxiy,1,-1
         slab(i,j)=.5*(slab(i,j-1)+slab(i,j))
      enddo
      enddo
c
c   Put "rleft" values into slab.
c
      do i=1,maxiy
         slab(i,1)=rleft(i)
      enddo
c
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine premaptform(dskmc,miycors,mjxcors,nproj,
     &   xlatc,xlonc,true1,true2,iup)
c
c   This routine calculates constants for routine maptform and puts
c   them in common block mptf, which should be declared in
c   (and only in) the main program and routines premaptform
c   (this routine) and maptform.
c
      common /mptf/ rpd_mptf,pi_mptf,dskmc_mptf,xlonc_mptf,rearth_mptf,
     &   ciy_mptf,cjx_mptf,yc_mptf,ihm_mptf,c1_mptf,c2_mptf,cone_mptf,
     &   conei_mptf,nproj_mptf
c
      pi_mptf=4.*atan(1.)     ! 3.1415...
      rpd_mptf=pi_mptf/180.    ! radians per degree
      rearth_mptf=6370.949  ! radius of planet, in km
      dskmc_mptf=dskmc
      xlonc_mptf=xlonc
      nproj_mptf=nproj
      ciy_mptf=.5*(1.+miycors)
      cjx_mptf=.5*(1.+mjxcors)
c
      if (nproj_mptf.eq.3) then   ! Mercator
c
      true1=0.
      true2=0.
      ihm_mptf=1
      cone_mptf=1.
      conei_mptf=1.
      c1_mptf=1.
      c2_mptf=1.
      yc_mptf=rearth_mptf*log((1.+sin(rpd_mptf*xlatc))/
     &   cos(rpd_mptf*xlatc))
c
      else   ! Lambert Comformal or Polar Stereographic
c
c   Make sure xlatc, true1, and true2 are all in same hemisphere,
c      and calculate ihm_mptf.
c
      if (xlatc.gt.0..and.true1.gt.0..and.true2.gt.0.) then
         ihm_mptf=1
      elseif (xlatc.lt.0..and.true1.lt.0..and.true2.lt.0.) then
         ihm_mptf=-1
      else
         write(iup,*)'Invalid latitude parameters for map.'
         stop
      endif
c
c   Calculate cone factor
c
      if (nproj_mptf.eq.1) then
         if (true1.ne.true2) then
            cone_mptf=log10(cos(rpd_mptf*true1)/cos(rpd_mptf*true2))/
     &           log10(tan(.25*pi_mptf-ihm_mptf*.5*rpd_mptf*true1)/
     &                 tan(.25*pi_mptf-ihm_mptf*.5*rpd_mptf*true2))
         else
            cone_mptf=cos(rpd_mptf*(90.-ihm_mptf*true1))
         endif
      elseif (nproj_mptf.eq.2) then
         cone_mptf=1.
      endif
c
c   Calculate other constants
c
      conei_mptf=1./cone_mptf
      cotrue1=ihm_mptf*90.-true1
      if (nproj_mptf.eq.1) then
         c1_mptf=rearth_mptf*sin(rpd_mptf*cotrue1)/
     &      (cone_mptf*(ihm_mptf*tan(.5*rpd_mptf*cotrue1))**cone_mptf)
         c2_mptf=tan(.5*rpd_mptf*cotrue1)*(cone_mptf/
     &      (ihm_mptf*rearth_mptf*sin(rpd_mptf*cotrue1)))**conei_mptf
         yc_mptf=-c1_mptf*(ihm_mptf*tan(.25*(ihm_mptf*pi_mptf-
     &      2.*rpd_mptf*xlatc)))**cone_mptf
      elseif (nproj_mptf.eq.2) then
         c1_mptf=1.+cos(rpd_mptf*cotrue1)
         c2_mptf=1.
         yc_mptf=-rearth_mptf*sin(.5*ihm_mptf*pi_mptf-rpd_mptf*xlatc)*
     &      c1_mptf/(1.+cos(.5*ihm_mptf*pi_mptf-rpd_mptf*xlatc))
      endif
c
      endif
c
      return
      end
c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine maptform(riy,rjx,rlat,rlon,idir,rrota)
c
c   This routine converts a coarse domain dot grid point, <riy,rjx>,
c   into a lat/lon point <rlat,rlon> if idir=1, or vice versa if
c   idir=-1. It works for Lambert Conformal (LC,1),
c   Polar Stereographic (ST,2), or Mercator (ME,3) projections,
c   with any true latitide(s).
c   It is assumed that premaptform has been called prior to this so
c   that the proper constants have been placed in the common block
c   called mptf, which should be declared in (and only in) the
c   main program and routines maptform (this routine) and premaptform.
c
      common /mptf/ rpd_mptf,pi_mptf,dskmc_mptf,xlonc_mptf,rearth_mptf,
     &   ciy_mptf,cjx_mptf,yc_mptf,ihm_mptf,c1_mptf,c2_mptf,cone_mptf,
     &   conei_mptf,nproj_mptf
c
      if (idir.eq.1) then   ! First, deal with idir=1
c
      ypoint=(riy-ciy_mptf)*dskmc_mptf+yc_mptf
      xpoint=(rjx-cjx_mptf)*dskmc_mptf
c
      if (nproj_mptf.eq.3) then
         rlat=(2.*atan(exp(ypoint/rearth_mptf))-.5*pi_mptf)/rpd_mptf
         rlon=xlonc_mptf+(xpoint/rearth_mptf)/rpd_mptf
      elseif (nproj_mptf.eq.1) then
         rlat=(.5*ihm_mptf*pi_mptf-2.*atan(c2_mptf*
     &      (sqrt(xpoint**2+ypoint**2))**conei_mptf))/rpd_mptf
         rlon=xlonc_mptf+(conei_mptf*atan2(xpoint,
     &      -ihm_mptf*ypoint))/rpd_mptf
      elseif (nproj_mptf.eq.2) then
         rlat=(.5*ihm_mptf*pi_mptf-ihm_mptf*2.*atan(sqrt(xpoint**2+
     &      ypoint**2)/(rearth_mptf*c1_mptf)))/rpd_mptf
         if(xpoint.eq.0..and.ypoint.eq.0.) then
            rlon=xlonc_mptf
         else
            rlon=xlonc_mptf+(atan2(xpoint,-ihm_mptf*ypoint))/rpd_mptf
         endif
      endif
c MGD mod
c since everything is initially in gridpoints, when we convert to lat/lon
c we should add in the rotation to the final longitude.
c going from lat/lon to gridpoints needs nothing, since lat/lon will have
c been converted, and thus already have accounted for the rotation.
      rlon=mod(rlon+900.+rrota,360.)-180.
c
      else   ! Otherwise, deal with idir=-1
c
      dlon=rlon-xlonc_mptf
      if (dlon.lt.-180.) dlon=dlon+360
      if (dlon.gt. 180.) dlon=dlon-360
      if (nproj_mptf.eq.3) then
         ypoint=rearth_mptf*log((1.+sin(rpd_mptf*rlat))/
     &      cos(rpd_mptf*rlat))
         xpoint=dlon*rpd_mptf*rearth_mptf
      elseif (nproj_mptf.eq.1) then
         ypoint=-c1_mptf*(ihm_mptf*tan(.25*(ihm_mptf*pi_mptf-
     &      2.*rpd_mptf*rlat)))**cone_mptf*cos(cone_mptf*rpd_mptf*dlon)
         xpoint=ihm_mptf*c1_mptf*(ihm_mptf*tan(.25*(ihm_mptf*pi_mptf-
     &      2.*rpd_mptf*rlat)))**cone_mptf*sin(cone_mptf*rpd_mptf*dlon)
      elseif (nproj_mptf.eq.2) then
         ypoint=-rearth_mptf*sin(.5*ihm_mptf*pi_mptf-rpd_mptf*rlat)*
     &      c1_mptf/(1.+cos(.5*ihm_mptf*pi_mptf-rpd_mptf*rlat))*
     &      cos(rpd_mptf*dlon)
         xpoint=ihm_mptf*rearth_mptf*sin(.5*ihm_mptf*pi_mptf-rpd_mptf*
     &      rlat)*c1_mptf/(1.+cos(.5*ihm_mptf*pi_mptf-rpd_mptf*rlat))*
     &      sin(rpd_mptf*dlon)
      endif
      riy=(ypoint-yc_mptf)/dskmc_mptf+ciy_mptf
      rjx=xpoint/dskmc_mptf+cjx_mptf
c
      endif
c
      return
      end
