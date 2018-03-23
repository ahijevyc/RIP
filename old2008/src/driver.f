c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine driver(miy,mjx,mkzh,mabpl,morpl,xtimeavl,
     &     cxtimeavl,ncxc,maxtavl,nxtavl,casename,iendc,
     &     title,rip_root,rootname,iendcr,ptimes,iptimes,
     &     ptuse,maxptimes,ptimeunits,tacc,ntextq,ntextcd,
     &     idotser,idotitle,timezone,iusdaylightrule,
     &     inearesth,iinittime,ivalidtime,fcoffset,
     &     titlecolor,idescriptive,icgmsplit,maxfld,
     &     itrajcalc,iusectrv,rtim,ctim,dtfile,dttraj,
     &     vctraj,ihydrometeor,xjtraj,yitraj,zktraj,diag,
     &     ntraj,imakev5d)
c
c   This subroutine was formerly the main program of RIP.  It does
c   most of the "work".
c
c   miy, and mjx are dot-point dimensions, in the x and y directions
c      respectively, of the domain to be analyzed.
c   mkzh is number of vertical levels in the domain.
c   mabpl and morpl are the maximum number of abscissa and ordinate
c      gridpoints, respectively, in a plot
c   xtimeavl is an array containing the xtimes that are available
c      for processing in this dataset.
c   cxtimeavl is the same as xtimeavl, but in character form.  It is
c      used for constructing file names.
c   ncxc is the number of characters in cxtimeavl (either 9 or 10,
c      depending on whether an older or newer version of RIPDP was run).
c   maxtavl is the dimension of xtimeavl and cxtimeavl.
c   nxtavl is the number of actual values in xtimeavl and cxtimeavl.
c   casename is a character string containing the case name for the
c      dataset.
c   iendc is number of meaningful characters in the variable
c      'casename.'
c   title,rip_root,rootname,ptimes,iptimes,ptimeunits,tacc,ntextq,ntextcd,
c      idotser,idotitle,timezone,iusdaylightrule,inearesth,
c      iinittime,ivalidtime,fcoffset,titlecolor,idescriptive,
c       and icgmsplit are namelist variables
c      that used to be read in in subroutine driver, but are now
c      read in the main program and passed to driver.
c   ptuse is a variable that will hold all the actual times requested
c      (i.e. with the time series in ptimes/iptimes expanded out)
c   maxptimes is the dimension of ptimes, iptimes, and ptuse
c   iendcr is number of meaningful characters in the variable
c      'rootname.'
c   maxfld is also a namelist variable, used to dimension the work
c      array, and is typically set to a value around 8 - 10.
c   itrajcalc is a flag that determines whether RIP is being run
c      in trajectory calculation mode (0=no, 1=yes)
c   rtim,ctim,dtfile,dttraj,vctraj,ihydrometeor,xjtraj,yitraj,
c      and zktraj are namelist variables for trajectory calculation.
c   diag is an array to hold diagnostic quantities, in trajectory
c      calculation mode
c   ntraj is the number of trajectories specified in xjtraj,yitraj,
c      zktraj, and diag
c   iusectrv is a flag to use contrived fields.
c   imakev5d is a flag that determines whether RIP is being run
c      in "make vis5d data" mode
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
c   Note: if you don't have a Fortran 90 compiler, or your
c   Fortran 77 compiler doesn't support adjustable dimensioning of
c   local (non-argument) arrays, then you must do the following:
c   comment out the above subroutine declaration, and uncomment
c   the following one:
c
c      subroutine driver(miydumb,mjxdumb,mkzhdumb,
c     &     mabpldumb,morpldumb,xtimeavl,
c     &     cxtimeavl,ncxc,maxtavl,nxtavl,casename,iendc,
c     &     title,rip_root,rootname,iendcr,ptimes,iptimes,
c     &     ptuse,maxptimes,ptimeunits,tacc,ntextq,ntextcd,
c     &     idotser,idotitle,timezone,iusdaylightrule,
c     &     inearesth,iinittime,ivalidtime,fcoffset,
c     &     titlecolor,idescriptive,icgmsplit,maxflddumb,
c     &     itrajcalcdumb,iusectrv,rtim,ctim,dtfile,dttraj,
c     &     vctraj,ihydrometeor,xjtraj,yitraj,zktraj,diag,
c     &     ntraj,imakev5d)
c
c   Then, uncomment the following include statement.
c
c       include 'parinc'
c
c   Also check rotuines qgomg.f, and saweli.f, for
c   instructions on similar adjustments that need to be made for
c   those routines.
c
c   Then, make a file called 'parinc' that has the following line,
c   but without the comment character:
c
c       parameter(miy=91,mjx=101,mkzh=27,maxfld=10,itrajcalc=1)
c       parameter (mabpl=700, morpl=350)
c
c   All of these parameters are described above.  Set them for your
c   particular domain and RIP requirements.  Note that maxfld and
c   itrajcalc must have the same values here as they do in the
c   &userin namelist.
c
c   Finally, run make to create a new executable.
c
c   Once these changes are made, you only need to edit 'parinc'
c   and run make to create an executable for a different size domain.
c
c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c
      dimension xtimeavl(maxtavl)
      character cxtimeavl(maxtavl)*10,casename*256
      dimension ptimes(maxptimes),iptimes(maxptimes),ptuse(maxptimes)
      character title*80,rootname*256,titlecolor*40,rip_root*256,
     &   ptimeunits*1
      dimension xjtraj(ntraj),yitraj(ntraj),zktraj(ntraj),diag(ntraj)
      character vctraj*1
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
c   Other parameters (which shouldn't need to be changed
c      for most applications):
c
c   maxpl, maxlev, maxfr are the maximum number of plots, plot levels,
c      and frames respectively (needed for dimensioning).
c   maxbuf is the maximum number of gflash buffers that can be saved
c   maxcosq is the max number of colors that can be defined for a
c      color sequence
c   maxtserv, maxtsers, and maxtsert are the maximum number of
c      time series variables, stations, and times, respectively.
c   niam,ncs are dimensions for work arrays for NCAR Graphics AREAS
c      routines
c
      parameter (maxpl=1500, maxlev=50, maxfr=300)
      parameter (maxbuf=220, maxcosq=30)
      parameter (maxtserv=10,maxtsers=30,maxtsert=19)
      parameter (niam=2000000,ncs=200000)
c
      dimension uuu(miy,mjx,mkzh), vvv(miy,mjx,mkzh),
     &   tmk(miy,mjx,mkzh), qvp(miy,mjx,mkzh), prs(miy,mjx,mkzh),
     &   ght(miy,mjx,mkzh), www(miy,mjx,mkzh),
     &   dmap(miy,mjx), xmap(miy,mjx), ter(miy,mjx),
     &   cor(miy,mjx), sfp(miy,mjx),
     &   unorth(miy,mjx), vnorth(miy,mjx),
     &   vc3d(miy,mjx,mkzh), wk(miy,mjx,mkzh*maxfld), prssou(mkzh),
     &   pslab1(mabpl,morpl), pslab2(mabpl,morpl),
     &   vcground(mabpl)
c
      dimension nltf(maxfr), numppf(maxfr), intim(maxfr),
     &   rtime(100,maxfr), ixwin(2,maxpl), iywin(2,maxpl),
     &   rcrag(2,maxpl), rcrbg(2,maxpl), indwk(3,maxpl),
     &   icdwk(maxpl), icong(maxpl), iconl(maxpl), icozr(maxpl),
     &   ilwll(maxpl), ilwng(maxpl), ilwnl(maxpl), ilwzr(maxpl),
     &   idall(maxpl), idang(maxpl), idanl(maxpl), idazr(maxpl),
     &   ilcnl(maxpl), ilczr(maxpl), ilcbr(maxpl), ipwlb(maxpl),
     &   iorlb(maxpl), igdir(maxpl), ihvbr(maxpl),
     &   ipwhl(maxpl), ipwbr(maxpl), ifclb(maxpl), ifcnl(maxpl),
     &   ifczr(maxpl), ifchl(maxpl), ilclo(maxpl), ifclo(maxpl),
     &   rwdbr(maxpl), ismcp(maxpl), ismth(maxpl), rcint(maxpl),
     &   rcbeg(maxpl), idash(maxpl), icolr(maxpl), icoll(maxpl),
     &   ilcll(maxpl), ilchl(maxpl), rtslb(maxpl), rtshl(maxpl),
     &   imjsk(maxpl), icomg(maxpl), icosq(maxcosq,maxpl),
     &   rcosq(maxcosq,maxpl), incon(maxpl), rcend(maxpl),
     &   rslcg(2,maxpl), ilinw(maxpl), rlevl(maxlev,maxpl),
     &   rlavl(maxlev,maxpl), rlevs(maxlev,maxpl), inlvs(maxpl),
     &   ioulw(maxpl), iouco(maxpl), imfco(6,maxpl), rsepa(32,maxpl),
     &   idimn(maxpl), rvwin(2,maxpl), ixavg(maxpl),
     &   iintv(maxpl), rvcmx(maxpl), ivvnx(maxpl), rvvms(maxpl),
     &   incsq(maxpl), raxlg(maxpl), raxld(maxpl), raxlv(maxpl),
     &   raxtg(maxpl), raxtd(maxpl), raxtv(maxpl), raddf(maxpl),
     &   rstrm(2,maxpl), rstmv(2,maxpl), rrfst(4,maxpl),
     &   rtjsp(3,50,maxpl),itjns(maxpl),iqgsm(maxpl),
     &   itjid(30,maxpl),itjni(maxpl),rtjar(2,maxpl),rtjst(maxpl),
     &   rtjen(maxpl),rtjti(maxpl),rdiff(maxpl),rrota(maxpl),
     &   pavprof(1000)
c
      dimension ixwingf(2,maxbuf), iywingf(2,maxbuf),
     &   yicorngf(maxbuf), xjcorngf(maxbuf), icolrgf(maxbuf),
     &   ilinwgf(maxbuf), idashgf(maxbuf), rtslbgf(maxbuf),
     &   rcintgf(maxbuf), ioulwgf(maxbuf), ioucogf(maxbuf),
     &   imfcogf(6,maxbuf), iwhatgf(maxbuf), raxlggf(maxbuf),
     &   raxtggf(maxbuf), rcragvc(2), rcrbgvc(2), rslcgprv(2),
     &   fred(0:255), fgreen(0:255), fblue(0:255),
     &   iam(niam), xcs(ncs), ycs(ncs),
     &   tserdat(maxtsert,maxtserv,maxtsers), tseryi(maxtsers),
     &   tserxj(maxtsers), mdatetser(maxtsert),
     &   rhourtser(maxtsert),rrotagf(maxbuf)
c
      logical lnobr(maxpl),lnozr(maxpl),lnolb(maxpl),
     &   lnmsg(maxpl),lmult(maxpl),lnttl(maxpl),lverf(maxpl),
     &   lchfl(maxpl),larng(maxpl),lhide(maxpl),lnmin(maxpl),
     &   lredo(maxpl),lgrad(maxpl),llapl(maxpl),lhadv(maxpl),
     &   lhodo(maxpl),lmand(maxpl),lnogd(maxpl),ldfrl(maxpl),
     &   lnsmm(maxpl),lnvlb(maxpl),lsndg(maxpl),lbogs(maxpl),
     &   lplrs(maxpl)
      logical lhodogf(maxbuf),lsndggf(maxbuf),lmandgf(maxbuf),lnogdvc,
     &   ldiffsuccess
      character cfeld(3,maxpl)*10,cptyp(maxpl)*2,ccmth(maxpl)*4,
     &   cvcor(maxpl)*1,cdum*80,alphabet*52,cfulb(maxpl)*5,
     &   csloc(2,maxpl)*20,ccrsa(2,maxpl)*20,ccrsb(2,maxpl)*20,
     &   csout(maxpl)*58,ctjfl(maxpl)*256,cdiff(maxpl)*256,
     &   csave(maxpl)*10,ctitl(maxpl)*82,engplttl(maxpl)*36,
     &   conam(0:255)*40,cv5nm(maxpl)*8,
     &   titlestr*82,cgmname*80,vcncheck*1,cnohl(maxpl)*1,
     &   vc3dtyp*1,vc2dtyp*1,unwk(maxpl)*24,
     &   cmllmgf(maxbuf)*5,coutygf(maxbuf)*32,coudsgf(maxbuf)*5,
     &   cmllm(maxpl)*5,couty(maxpl)*32,couds(maxpl)*5,
     &   icaoid*4,locdesc*44,csids(40,maxpl)*20,minfo*256,
     &   tserlocpoint(maxtsers)*58,
     &   tserlocdesc(maxtsers)*44,tservname(maxtserv)*82
c
      dimension nsids(maxpl)
c
      dimension uusv(1+(miy-1)*itrajcalc,1+(mjx-1)*itrajcalc,
     &               1+(mkzh-1)*itrajcalc),
     &   vvsv(1+(miy-1)*itrajcalc,1+(mjx-1)*itrajcalc,
     &        1+(mkzh-1)*itrajcalc),
     &   ehsv(1+(miy-1)*itrajcalc,1+(mjx-1)*itrajcalc,   
     &         1+(mkzh-1)*itrajcalc),                    ! eh: "expon. height"
     &   ehdotsv(1+(miy-1)*itrajcalc,1+(mjx-1)*itrajcalc,   
     &         1+(mkzh-1)*itrajcalc)                     ! d(eh)/dt
c
      integer fcasthr, fchour(maxtsert)
c
c   Array for call to GKS routine GSASF
c
      dimension iasf(13)
      data iasf / 1,1,1,1,1,1,1,1,1,1,1,1,1 /
c
      include 'comconst'
      include 'comvctran'
      include 'pointers'
c
c   NCAR Graphics common blocks
c
      dimension mconcp(1000),icoindcp(1000)
      common /cpack/ mconcp,icoindcp,nconarea,icpfchl,icpfclo,
     &   icpfclb,icpfcnl,icpfczr
      common /emap/ llcolor,lllinw,llndot,mpfillco(6),llmask,ioutype
      common /str03/  inita , initb , arowl , iterp , iterc , igflg
     +             ,  imsg , uvmsg , icyc , displ , dispc , cstop
      common / stpar /
     +                iud1       ,ivd1       ,ipd1       ,
     +                ixd1       ,ixdm       ,iyd1       ,iydn       ,
     +                ixm1       ,iym1       ,ixm2       ,iym2       ,
     +                iwkd       ,iwku       ,iset       ,ierr       ,
     +                ixin       ,iyin       ,imsk       ,icpm       ,
     +                nlvl       ,ipai       ,ictv       ,wdlv       ,
     +                uvmn       ,uvmx       ,pmin       ,pmax       ,
     +                ithn       ,iplr       ,isst       ,
     +                iclr(64)           ,tvlu(64)
c
c   Map transform common block
c
      common /mptf/ rpd_mptf,pi_mptf,dskmc_mptf,xlonc_mptf,rearth_mptf,
     &   ciy_mptf,cjx_mptf,yc_mptf,ihm_mptf,c1_mptf,c2_mptf,cone_mptf,
     &   conei_mptf,nproj_mptf
c
c   Stationlist common block
c
      common /sl/ igot_sl,ns_sl,icao_sl,iwmo_sl,slat_sl,slon_sl,loc_sl
      dimension iwmo_sl(15000),slat_sl(15000),slon_sl(15000)
      character icao_sl(15000)*4,loc_sl(15000)*44
c
c some architectures (our SGI-64) don't initialize these properly, so do
c it here.
c
      COMMON /VEC1/   ASH        ,EXT        ,ICTRFG     ,ILAB       ,
     +                IOFFD      ,IOFFM      ,ISX        ,ISY        ,
     +                RMN        ,RMX        ,SIDE       ,SIZE       ,
     +                XLT        ,YBT        ,ZMN        ,ZMX
C
      COMMON /VEC2/   BIG        ,INCX       ,INCY
c
c   Vis5d variables
c
      parameter (maxv5delements=2000000)
      include 'v5df.h'
      external v5dcreate,v5dclose,v5dwrite
      integer nr, nc, nl(MAXVARS)
      integer numtimes
      integer numvars
      character varname(MAXVARS)*10
      integer datestamp(MAXTIMES)
      integer timestamp(MAXTIMES)
      integer compress
      integer projection
      real proj_args(100)
      integer vertical
      real vert_args(MAXLEVELS)
      pointer(i_v5darray,v5darray(maxv5delements))
c MGD vars
      real fbmino,ftmaxo,flmino,frmaxo,vdif
c
      write(iup,*)'Welcome to your friendly RIP output file !'
c     call flush(iup)
      igotit_sl=0
c
c velvct common variables
c
      incx = 1
      incy = 1
      ext = .25
      ictrfg = 1
      ilab = 0
      ioffd = 0
      ioffm = 0
      rmn = 160.
      rmx = 6400.
      side = 0.90
      size = 256.
      xlt = 0.05
      ybt = 0.05
c
c   Convert tacc from seconds to hours
c
      tacch=tacc/3600.
c
c   Constants for gks
c
      ier=6             ! error output
      iwkidcgm=1        ! workstation id for cgm file output
      iucgm=3           ! fortran unit number for cgm file output
      iwtypecgm=1       ! workstation type (1=cgm file output)
      iwkidgf=9         ! workstation id for gflash utility
      iugf=4            ! fortran unit number for gflash utility
      iwtypegf=3        ! workstation type (3=gflash utility)
      alphabet=
     &   'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'
c
c   Initialize the gflash buffer counter
c
      ngfbuf=0
c
c   For gflash buffers, in addition to other specific parameters,
c   the variable iwhatgf indicates the general type of background
c   saved in a particular buffer:
c       1: map background
c       2: skewt background
c       3: horizontal tick marks
c       4: polar skewt background
c
c   Plot settup: The plots (except for axis labels and title) will
c      be confined to a "usable rectangle" with left and right edges
c      and top and bottom at flmin, frmax, ftmax, and fbmin respect-
c      ively (in the fractional coordinate system).  The size of the
c      possibly non-square plot will be adjusted to fill as much of the
c      "usable rectangle" as possible.  The plot title will appear in
c      the space above the "usable rectangle", and the axis labels will
c      appear in the left and bottom margins.
c
c   Read the color table
c
      call rdcolt (rip_root,nco,conam,fred,fgreen,fblue)
c
c   Get information from header record that was read in main program
c
      nxt=1
      write(iup,*) 'Get header info.'
      call getheadinfo(casename,iendc,xtimeavl,cxtimeavl,
     &   nxt,maxtavl,ncxc,nproj,miycors,mjxcors,mdateb,mhourb,iice,
     &   iplevdata,true1,true2,xlatc,
     &   xlonc,dskmc,dskm,yicorn,xjcorn,rhourb,dsc,ds,refrat,iup)
c
c   Set values of rip file header variables
c
      call setripheader(ihrip,rhrip,chrip,miy,mjx,mkzh)
c
      if (itrajcalc.eq.0) then
c
c      If using iptimes, convert the mdates in the iptimes array to
c         xtimes in the ptimes array.  Also, determine nptimes.
c
         if (ptimes(1).lt.0..or.iptimes(1).lt.0.or.(ptimes(1).eq.
     &       9e9.and.iptimes(1).eq.99999999)) then !user wants all times
            nptimes=0
            write(iup,*)'Note: RIP will plot all available times.'
         elseif (ptimes(1).ne.9e9.and.iptimes(1).ne.99999999) then
            write(iup,*)'Can''t use both ptimes and iptimes.'
            stop
         endif
         if (iptimes(1).ne.99999999) then
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
 259        continue
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
 261        continue
         endif
c
c      Process time sequences in ptimes array.
c
         ii=0
         itime=0
  100    ii=ii+1
         if (ii.gt.nptimes) goto 120
         if (ptimes(ii).ge.0.) then
            itime=itime+1
            if (itime.gt.maxptimes) then
               write(iup,*)
     &            'Number of times requested exceeds maxptimes.'
               write(iup,*)'Increase maxptimes in rip code, recompile,'
               write(iup,*)'and run rip again.'
               stop
            endif
            ptuse(itime)=ptimes(ii)
         else
            ii=ii+1
            if (ptimes(ii).gt.0.) then
               tstart=ptimes(ii-2)
               tend=-ptimes(ii-1)
               tinc=ptimes(ii)
               tdist=tend-tstart
               isign=nint(tdist/abs(tdist))
               ntseries=int(abs(tdist)/tinc+.00001) + 1
               do i=2,ntseries
                  itime=itime+1
                  if (itime.gt.maxptimes) then
                     write(iup,*)
     &                 'Number of times requested exceeds maxptimes.'
                     write(iup,*)
     &                 'Increase maxptimes in rip code, recompile,'
                     write(iup,*)'and run rip again.'
                     stop
                  endif
                  ptuse(itime)=ptuse(itime-1)+isign*tinc
               enddo
            else
               write(iup,*)'Error in ptimes sequence specification.'
               stop
            endif
         endif
         goto 100
  120    nptuse=itime
c
      elseif (itrajcalc.eq.1) then
c
c      First set ptuse.  ptimes or iptimes that were set in the namelist
c         are ignored (i.e. they don't matter)
c
         nptuse=nint(abs(rtim-ctim)*3600./dtfile)+1
         do i=1,nptuse
            ptuse(i)=rtim+(i-1.)*(ctim-rtim)/(nptuse-1.)
         enddo
c
c      Test to make sure all specified times are available.
c
         iunav=0
         do i=1,nptuse
            do j=1,nxtavl
               if (abs(xtimeavl(j)-ptuse(i)).le.tacch) goto 46
            enddo
            iunav=iunav+1
            write(iup,*)'Time ',ptuse(i),' is not available.'
 46         continue
         enddo
         if (iunav.gt.0) then
           write(iup,*)'Change your user input file and run rip again.'
           stop
         endif
c
      endif
      call premaptform(dskmc,miycors,mjxcors,nproj,
     &   xlatc,xlonc,true1,true2,iup)
c
c   Zero the vertical coordinate array
c
      call fillarray(vc3d,miy*mjx*mkzh,0.)
c
c   Read the plot specifications.
c
      call readspec      (nfr,numppf,nltf,intim,rtime,cfeld,icomg,
     &   cptyp,incon,icolr,icosq,rcosq,icoll,ilcll,ilchl,rtslb,rtshl,
     &   icong,iconl,icozr,ilwll,ilwng,ilwnl,ilwzr,idall,
     &   idang,idanl,idazr,ilcnl,ilczr,
     &   ilcbr,ipwlb,iorlb,ipwhl,ipwbr,ifclb,ifcnl,ifczr,ifchl,
     &   ilclo,ifclo,ccmth,rwdbr,ihvbr,ccrsa,ccrsb,csloc,rsepa,
     &   imjsk,nptuse,ptuse,ixwin,iywin,lhide,csave,lredo,lnogd,
     &   ilinw,lnobr,lnozr,cnohl,lnolb,lnmsg,lnttl,lverf,cvcor,lnmin,
     &   rlevs,inlvs,rlevl,rlavl,rvwin,ixavg,idash,ismth,ismcp,iintv,
     &   rvcmx,ivvnx,rvvms,rcint,rcbeg,rcend,lmult,larng,lchfl,lhodo,
     &   lmand,lsndg,cdiff,rdiff,ldfrl,lbogs,
     &   raxlg,raxld,raxlv,raxtg,raxtd,raxtv,rstrm,rrfst,raddf,
     &   cfulb,conam,incsq,cmllm,couty,iqgsm,iopenverf,
     &   ctjfl,ctitl,cv5nm,rtjsp,itjns,itjid,itjni,
     &   rtjar,rtjst,rtjen,rtjti,lgrad,llapl,lhadv,igdir,
     &   couds,ioulw,iouco,imfco,fred,fgreen,fblue,nco,
     &   csids,nsids,lnsmm,lnvlb,itrajcalc,imakev5d,
     &   maxcosq,maxfr,maxlev,maxpl,miy,mjx,mkzh,rrota,lplrs)
c
c   Open the precip verification file if necessary
c
      if (iopenverf.gt.0)
     &   open (unit=iuprcver,file=rootname(1:iendcr)//'.prcver',
     &   form='formatted',status='unknown')
c
c   Set up time series stuff if needed.
c
      if (idotser.eq.1) call tserprep(tseryi,tserxj,
     &   tserlocpoint,tserlocdesc,ntsers,ntsert,maxtsers,
     &   rip_root,rrota(ipl),miy,mjx)
c
c   Get "minfo" information
c
      call getminfo(casename,iendc,igotminfo,minfo,iendmi,iup)
c
c   Set trivial values for some plotting quantities if
c   calculating trajectories or making vis5d files.
c
      if (itrajcalc.eq.1.or.imakev5d.eq.1) then
c
c   For all plot spec statements, force cvcor to be the same as vctraj.
c   For all frames, force nltf (number of levels in this frame) to be
c   1.  Whatever vcor values were set in the plspec table are
c   overwritten (i.e. they don't matter), and levl values are also
c   ignored.
c
      if (imakev5d.eq.1) nv5dlevels=nltf(1)
      do ifr=1,nfr
         ipl=ifr
         if (itrajcalc.eq.1) then
            cvcor(ipl)=vctraj
         elseif (imakev5d.eq.1) then
            cvcor(ipl)='z'
         endif
         nltf(ifr)=1
         if (numppf(ifr).ne.1) then
            write(iup,*)'For trajectory calculation or producing',
     &         ' Vis5d data there should only be'
            write(iup,*)'ONE plot specification statement per frame.'
            write(iup,*)'For frame #',ifr,' there are ',numppf(ifr),
     &         ' PSSs.'
            stop
         endif
      enddo
c
      endif
c
      if (itrajcalc.eq.1) then
c
c   Just in case your dttraj doesn't divide evenly into your dtfile:
c
      ntrajtime=nint(dtfile/dttraj)
      dttraj=dtfile/ntrajtime
c
      ntrajdir=1
      if (rtim.gt.ctim) ntrajdir=-1
c
c   Open files
c
      open (unit=iutrajout,file=rootname(1:iendcr)//'.traj',
     &   form='unformatted',status='unknown')
      open (unit=iudiagout,file=rootname(1:iendcr)//'.diag',
     &   form='unformatted',status='unknown')
c
      endif
c
c   Initialize vis5d file
c
      if (imakev5d.eq.1) then
c
c   Use "cgmname" to hold vis5d file name
c
      cgmname=rootname(1:iendcr)//'.v5d'
c
      nr=iywin(2,1)-iywin(1,1)  ! nr and nc are number of cross points,
      nc=ixwin(2,1)-ixwin(1,1)  ! whereas xwin and ywin are dot points
      numtimes=nptuse
      numvars=nfr
c
c   Set number of levels for each variable
c
      do iv=1,numvars
         nl(iv)=nv5dlevels
      enddo
      do iv=numvars+1,MAXVARS
         nl(iv)=IMISSING
      enddo
c
      do iv=1,numvars
         if (varname(iv).eq.'uuu       ') then
            varname(iv)='U         '
         elseif (varname(iv).eq.'vvv       ') then
            varname(iv)='V         '
         elseif (varname(iv).eq.'www       ') then
            varname(iv)='W         '
         elseif (cv5nm(iv).eq.'samasvar') then
            varname(iv)=cfeld(1,iv)
         else
            varname(iv)=cv5nm(iv)
         endif
      enddo
      do iv=numvars+1,MAXVARS
         varname(iv)=' '
      enddo
c
c   Set timestamp and datestamp arrays
c
      do it=1,numtimes
c
c      First determine if this time is available.
c
         iavail=0
         do i=1,nxtavl
            if (abs(xtimeavl(i)-ptuse(it)).le.tacch) then
               iavail=1
               nxt=i
               goto 145
            endif
         enddo
 145     continue
         if (iavail.ne.1) then
            write(iup,*)'In "make vis5d" mode: requested time ',
     &         ptuse(it),' is not available. Stopping.'
            stop
         endif
         xtime=xtimeavl(nxt)
         hrspastmdateb=rhourb+xtime
         mhourtrunc=mhourb+int(hrspastmdateb)
         rhourtrunc=hrspastmdateb-float(mhourtrunc-mhourb)
         if (rhourtrunc.lt.0.0.or.rhourtrunc.ge.1.0) then
            write(iup,*)'In "make vis5d" mode: problem with'
            write(iup,*)'setting rhourtrunc.'
            stop
         endif
         call mconvert(mdatetrunc,mhourtrunc,-1,1940)
         iyy=mdatetrunc/1000000
         mdateyearstart=iyy*1000000+10100
         call mconvert(mdateyearstart,mhouryearstart,1,1940)
         ijd=1+(mhourtrunc-mhouryearstart)/24
         ihh=mod(mhourtrunc,24)
         imm=int(rhourtrunc*60.)
         iss=(rhourtrunc*3600.-imm*60.)
         if (iss.lt.0.0.or.iss.ge.60.) then
            write(iup,*)'In "make vis5d" mode: problem with'
            write(iup,*)'setting iss.   iss=',iss
            stop
         endif
         timestamp(it)=ihh*10000+imm*100+iss
         datestamp(it)=iyy*1000+ijd
      enddo
      do it=numtimes+1,MAXTIMES
         timestamp(it)=IMISSING
         datestamp(it)=IMISSING
      enddo
c
      compress=1
c
c   Set up map projection
c
      if (nproj.eq.1) then ! Lambert Conformal
         projection=2
         if (true1*true2.le.0.0) then
            write(iup,*)'Invalid true lats. for LC map proj.'
            stop
         endif
         proj_args(1) = max(true1,true2)
         proj_args(2) = min(true1,true2)
         if (true1.gt.0.) then
            call maptform(yipole,xjpole,90.,0.,-1,rrota(ipl))
         else
            call maptform(yipole,xjpole,-90.,0.,-1,rrota(ipl))
         endif
         yipole=1.+(yipole-yicorn)*refrat
         xjpole=1.+(xjpole-xjcorn)*refrat
         proj_args(3) = iywin(2,1)-yipole+0.5
         proj_args(4) = xjpole-.5-(ixwin(1,1)-1)
         proj_args(5) = -xlonc ! west long. is positive in Vis5d
         proj_args(6) = dskm
         do ipa=7,100
            proj_args(ipa) = MISSING
         enddo
      elseif (nproj.eq.2) then ! Polar Stereographic
         write(iup,*)'Sorry, can''t do Polar Stereographic'
         write(iup,*)'projection for Vis5d data conversion.'
         stop
      elseif (nproj.eq.3) then ! Mercator
         write(iup,*)'Sorry, can''t do Mercator projection'
         write(iup,*)'for Vis5d data conversion.'
         stop
      endif
c
c   Set level values
c      
      vertical=2
      do il=1,nv5dlevels
         vert_args(il)=rlevl(il,1)
      enddo
      do il=nv5dlevels+1,MAXLEVELS
         vert_args(il)=MISSING
      enddo
c
      iv5derr=v5dcreate( cgmname, numtimes, numvars, nr, nc, nl,
     &   varname, timestamp, datestamp, compress, projection,
     &   proj_args, vertical, vert_args )
      if (iv5derr.eq.0) then
         write(iup,*)'v5dcreate failed for some unknown reason.'
         stop
      endif
c
      iv5dcount=0
c
      endif
c
c   LOOP THROUGH TIME LEVELS.
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
c---------------------------------------------------------------------c
      do 1000 ipltime=1,nptuse     ! TIME LOOP
c---------------------------------------------------------------------c
c
c   First determine if this time is available.
c
      iavail=0
      do i=1,nxtavl
         if (abs(xtimeavl(i)-ptuse(ipltime)).le.tacch) then
            iavail=1
            nxt=i
            goto 40
         endif
      enddo
 40   continue
      write(iup,*)
      if (iavail.eq.1) then
         write(iup,*)'Requested time ',ptuse(ipltime),' is available.'
c        call flush(iup)
      else
         write(iup,*)'Requested time ',ptuse(ipltime),
     &      ' is not available.'
c        call flush(iup)
         goto 1000
      endif
      ifirstplot=1
      xtime=xtimeavl(nxt)
      hrspastmdateb=rhourb+xtime
      mhour=mhourb+int(hrspastmdateb)
      fcasthr=nint(hrspastmdateb)
      rhour=hrspastmdateb-float(mhour-mhourb)
      call mconvert(mdate,mhour,-1,1940)
c
c   Call getheadinfo again at each time, for the sole purpose of
c   obtaining possibly new values of yicorn and xjcorn, if the domain
c   has moved.
c
      call getheadinfo(casename,iendc,xtimeavl,cxtimeavl,
     &   nxt,maxtavl,ncxc,nproj,miycors,mjxcors,mdateb,mhourb,iice,
     &   iplevdata,true1,true2,xlatc,
     &   xlonc,dskmc,dskm,yicorn,xjcorn,rhourb,dsc,ds,refrat,iup)
c
c   Read the basic fields.
c
      call getbasicvars(casename,iendc,cxtimeavl,xtimeavl,nxt,
     &   ncxc,maxtavl,uuu,vvv,tmk,qvp,www,prs,ght,sfp,ter,
     &   dmap,xmap,cor,miy,mjx,mkzh)
c
c   This is where one would insert a "contrived fields" routine,
c   if one were so inclined.
c
c
c   Update the time-dependent elements of the rip file header arrays
c
      ihrip(11)=mdate
      rhrip(14)=rhour
      rhrip(15)=xtime
      rhrip(7)=yicorn
      rhrip(8)=xjcorn
c
c   Create unorth,vnorth (requires yicorn, xjcorn)
c
      do j=1,mjx
         rj1=xjcorn+(j-1.)/refrat
      do i=1,miy
         ri1=yicorn+(i-1.)/refrat
         call maptform(ri1,rj1,rlat,rlon,1,rrota(ipl))
         rlat2=min(89.9999,rlat+.1)
         call maptform(ri2,rj2,rlat2,rlon,-1,rrota(ipl))
         unn=rj2-rj1
         vnn=ri2-ri1
         dnorth=sqrt(unn**2+vnn**2)
         unorth(i,j)=unn/dnorth
         vnorth(i,j)=vnn/dnorth
      enddo
      enddo
c
c   Create pavprof, average vertical profile of pressure, for use with
c   the vvms keyword in cross-section and sounding vector plots.
c
      nstep=max(1,nint(sqrt((miy*mjx)/100.)))
      ntot=((miy-2)/nstep+1)*((mjx-2)/nstep+1)
      do k=1,mkzh
         pavprof(k)=0.
         do j=1,mjx-1,nstep
         do i=1,miy-1,nstep
            pavprof(k)=pavprof(k)+prs(i,j,k)
         enddo
         enddo
         pavprof(k)=pavprof(k)/float(ntot)
      enddo
c
c   Increment the time series time counter
c
      if (idotser.eq.1) then
         fchour(ntsert+1)=fcasthr
         ntsert=ntsert+1
         mdatetser(ntsert)=mdate
         rhourtser(ntsert)=rhour
         ntserv=0
      endif
c
c -------------------  P L O T    S E C T I O N  ---------------------
c
c   We now have all the necesary information to begin plotting
c   for this time level, so we will begin the plotting loop.
c   First make the cgm file name.
c
      if ((ipltime.eq.1.or.icgmsplit.eq.1).and.itrajcalc.eq.0.and.
     &     imakev5d.eq.0.and.nfr.gt.0) then
c
      if (icgmsplit.eq.1) then
         cgmname=
     &      rootname(1:iendcr)//'.cgm'//alphabet(ipltime:ipltime)
      else
         cgmname=rootname(1:iendcr)//'.cgm'
      endif
c
c   Open gks, name the cgm output file for this time, and open the
c      metafile workstation to unit 3 (instead of 2, as would be done
c      by a call to opngks)
c
      call gopks (ier,isz)
      call gesc (-1391,1,cgmname,1,1,cdum)
      call gopwk (iwkidcgm,iucgm,iwtypecgm)
      call gacwk (iwkidcgm)
c
c   Open gflash workstation to unit 4
c
      call gopwk (iwkidgf,iugf,iwtypegf)
c
c   Make sure all GKS aspect source flags are set to "individual", and
c   area fill is solid.
c
      call gsasf(iasf)
      call gsfais(1)
c
c   Assign the colors.
c
      do ico=0,nco
         call gscr (1,ico,fred(ico),fgreen(ico),fblue(ico))
      enddo
c      icomax=nco  ! Setting this here causes rip to run out of
c                  ! available colors (max=255)
c
c   Assign special value for CONPACK
c
      call cpsetr('SPV',rmsg)
      ivcs=0            ! needed for vertical coord. trans. in x-secs.
c
      endif
c
c   Reset the alternate vertical coordinate identifier, and identifiers
c      of cross-section end points for alternate vertical coordinate
c
      vc3dtyp='?'
      vc2dtyp='?'
      rcragvc(1)=rmsg
      rcragvc(2)=rmsg
      rcrbgvc(1)=rmsg
      rcrbgvc(2)=rmsg
      lnogdvc=.false.
c
c   Initialize previous sounding location.
c
      rslcgprv(1)=99999.
      rslcgprv(2)=99999.
c
c---------------------------------------------------------------------c
      do 950 ifr=1,nfr     ! FRAME LOOP
c---------------------------------------------------------------------c
c
c      if (ifr.gt.1) then   ! for testing purposes
c         fbmin=fbmin-.01
c         frmax=frmax-.01
c         ngfbuf=0
c      endif
      iplstrt=1
      do i=1,ifr-1
         iplstrt=iplstrt+numppf(i)
      enddo
      iplend=iplstrt+numppf(ifr)-1
      do itt=1,intim(ifr)
         if (abs(xtime-rtime(itt,ifr)).le.tacch) then
            goto 616
         endif
      enddo
      goto 950
  616 continue
c
c   Zero the work array.
c
      incwk=1
      maxslab=mkzh*maxfld
      call fillarray(wk,miy*mjx*maxslab,0.)
c
c---------------------------------------------------------------------c
      do 900 ilev=1,nltf(ifr)     ! LEVEL LOOP
c---------------------------------------------------------------------c
c
      if (itrajcalc.eq.1) then
         write(iup,*)'Working on trajectory diagnostic quantity ',ifr
      elseif (imakev5d.eq.1) then
         write(iup,*)'Working on Vis5D variable ',ifr
      else
         write(iup,*)
         write(iup,*)'Frame number ',ifr,',  Level number ',ilev
         write(iup,*)'------------------------------------------'
      endif
c     call flush(iup)
c
c   Initialize the "ceiling" of available space for text at the top of
c   the screen, and the "floor" of available space for messages and
c   labels at the bottom of the screen.  toptextclg will decrease as
c   new text is added at the top.  bottextfloor will increase as new
c   text/label bars are added at the bottom.
c
      toptextclg=.992
      bottextfloor=.008
      prfmax = -99999.
      prfmin =  99999.
      icomax=nco  ! this is the correct place to reset icomax=nco
c
      if (itrajcalc.eq.0.and.imakev5d.eq.0) then
c
c      Set text quality, font, function-code character,
c      and color for frame title.
c
         call pcseti ('QU',ntextq)
         call pcseti ('CD',ntextcd)
         call pcsetc ('FC','~')
         do i=0,nco
            if (titlecolor.eq.conam(i)) then
               call gsplci(i)
               call gstxci(i)
               goto 620
            endif
         enddo
  620    continue
c
c   Adjust mdateb, rhourb, and xtime based on user-specified
c      forecast offset.
c
         xtimefc=xtime-fcoffset
         rhourbfc=rhourb+fcoffset
         if (rhourbfc.ge.0) then
            mhadd=int(rhourbfc)
         else
            isafe=int(-rhourbfc)+2
            mhadd=int(rhourbfc+float(isafe))-isafe
         endif
         rhourbfc=rhourbfc-float(mhadd)
         call mconvert(mdatebfc,mhourb+mhadd,-1,1940)
c
c      Make title at top of frame.
c
         call frtitle(title,casename,iendc,rootname,
     &      mdatebfc,rhourbfc,xtimefc,timezone,iusdaylightrule,
     &      inearesth,idotitle,iinittime,ivalidtime,toptextclg)
c
c      Print "minfo" information at bottom of frame
c
         if (igotminfo.eq.1) then
            idominfo=1
            do ipl=iplstrt,iplend
               if ((cptyp(ipl)(1:1).eq.'s').or.lnmin(ipl)) then
                  idominfo=0
                  goto 179
               endif
            enddo
 179        continue
            if (idominfo.eq.1) then
               chsize=.009
               ypos=bottextfloor+.5*chsize
               call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
               call plchhq (.1,ypos,minfo(1:iendmi),chsize,0.,-1.)
               bottextfloor=bottextfloor+1.9*chsize
            endif
         endif
c
         call gsplci(1)
         call gstxci(1)
      endif
c
c---------------------------------------------------------------------c
      do 700 ipl=iplstrt,iplend     ! PLOT OVERLAY LOOP
c---------------------------------------------------------------------c
c
c      icomax=nco  ! resetting this here causes the problem with multiple
c                  ! color-filled overlays in one frame
c
c   Several things that need to be done only if ilev=1.
c
      if (ilev.eq.1) then
c
c   Interpret cross section location info.
c
      call locinterp(ccrsa(1,ipl),rcrag(2,ipl),rcrag(1,ipl),
     &   rlat,rlon,iwmo,icaoid,locdesc,rip_root,rrota(ipl))
      call locinterp(ccrsb(1,ipl),rcrbg(2,ipl),rcrbg(1,ipl),
     &   rlat2,rlon2,iwmo,icaoid,locdesc,rip_root,rrota(ipl))
c
c   Interpret sounding location info.
c
      call locinterp(csloc(1,ipl),rslcg(2,ipl),rslcg(1,ipl),
     &   rlat,rlon,iwmo,icaoid,locdesc,rip_root,rrota(ipl))
      if ((cptyp(ipl)(1:1).eq.'s'.and.cptyp(ipl)(2:2).ne.'b').or.
     &    (cptyp(ipl)(1:1).eq.'p') ) then
         csout(ipl)=' '
         write(csout(ipl),213)'x,y=      ,         lat,lon=',
     &      rlat,',',rlon
         if (rlat .lt. 0.) then
           ins = 1
         else
           ins = 0
         endif
         itdone=0
         if (icaoid.ne.'XXXX') then
            write(csout(ipl)(43:),214)'  stn=',icaoid
            itdone=1
         endif
         if (iwmo.ne.99999) then
            if (itdone.eq.0) then
               write(csout(ipl)(43:),215)'   stn=',iwmo
            else
               write(csout(ipl)(53:),216)',',iwmo
            endif
         endif
         sxgn=1.+(rslcg(2,ipl)-xjcorn)*refrat
         sygn=1.+(rslcg(1,ipl)-yicorn)*refrat
         if (sxgn.le..5.or.sxgn.ge.mjx-.5.or.
     &       sygn.le..5.or.sygn.ge.miy-.5) then
            write(iup,*)'The sounding requested at the following'
            write(iup,*)'location is outside the cross-point domain:'
            write(iup,*) csout(ipl)
            write(iup,*)'Please correct input file and re-execute RIP.'
            stop
         endif
  213    format(a28,f6.2,a1,f7.2)
  214    format(a6,a4)
  215    format(a7,i5)
  216    format(a1,i5)
      endif
c
c   Interpret storm speed.
c
      if (rstrm(1,ipl).ne.rmsg) then
         if (rstrm(2,ipl).eq.rmsg) then
            dy=rcrbg(1,ipl)-rcrag(1,ipl)
            dx=rcrbg(2,ipl)-rcrag(2,ipl)
            hypot=sqrt(dy*dy+dx*dx)
            rstmv(2,ipl)=rstrm(1,ipl)*dx/hypot
            rstmv(1,ipl)=rstrm(1,ipl)*dy/hypot
         else
            rstmv(2,ipl)=rstrm(2,ipl)
            rstmv(1,ipl)=rstrm(1,ipl)
         endif
      else
         rstmv(1,ipl)=0.
         rstmv(2,ipl)=0.
      endif
c
c   Calculate some cross section parameters.
c
      ydist=refrat*(rcrbg(1,ipl)-rcrag(1,ipl))
      xdist=refrat*(rcrbg(2,ipl)-rcrag(2,ipl))
      xseclen=sqrt(ydist**2+xdist**2)
      distmax=max(abs(xdist),abs(ydist))
      nscrs=2*nint(distmax)+1
c      nscrs=nint(distmax)+1
      nscross=nscrs  ! for vctran common block
      mkzhcross=mkzh  ! for vctran common block
c
c   Put the appropriate 2- or 3-d arrays into
c      the work 3-d array. This part is similar
c      to routine FIELDS in program SIGMA.
c
      if (cptyp(ipl)(2:2).ne.'b'.and.cptyp(ipl)(2:2).ne.'t') then
         call fields(cfeld,wk,indwk,icdwk,rlevl,rlavl,unwk,
     &      uuu,vvv,tmk,qvp,prs,ght,www,sfp,dmap,xmap,ter,
     &      cor,unorth,vnorth,rstmv,rrfst,pslab1,pslab2,
     &      incwk,ipl,iplstrt,idimn,rcrag,ismcp,
     &      rcrbg,cptyp,mdate,lverf,ydist,xdist,xseclen,
     &      nscrs,raddf,csave,lredo,ihrip,rhrip,rsepa,
     &      chrip,vardesc,lgrad,llapl,lhadv,
     &      igdir,iqgsm,plchun,casename,iendc,engplttl,ctjfl,rtjst,
     &      rtjen,cdiff,rdiff,ldfrl,xtime,tacch,
     &      xtimeavl,cxtimeavl,ncxc,maxtavl,nxtavl,nxt,ldiffsuccess,
     &      maxslab,maxlev,maxpl,miy,mjx,mkzh,mabpl,morpl)
      endif
c
      endif  ! end of stuff to do only if ilev=1
c
c   Check to make sure "level" is not out of bounds if vcor=s.
c
      if (cvcor(ipl).eq.'s') then
         rlevl(ilev,ipl)=min(float(mkzh),max(1.,rlevl(ilev,ipl)))
         isw=0
         if (rlavl(ilev,ipl).lt.0.) then
            isw=1
            rlavl(ilev,ipl)=-rlavl(ilev,ipl)
         endif
         rlavl(ilev,ipl)=min(float(mkzh),max(1.,rlavl(ilev,ipl)))
         if (isw.eq.1) rlavl(ilev,ipl)=-rlavl(ilev,ipl)
      endif
c
c   We now have the desired array(s) in the wk array so we are
c      ready so start the plotting. First, set text quality and font.
c
      if (itrajcalc.eq.0.and.imakev5d.eq.0) then
         call pcseti ('QU',ntextq)
         call pcseti ('CD',ntextcd)
      endif
c
c   Set flag that informs vinterp that a differenced field is being
c   plotted, so it knows not to do its fancy extrapolation below
c   ground.
c
      idiffflag=0
      if (cdiff(ipl)(1:5).ne.'none '.or.rdiff(ipl).ne.rmsg)
     &   idiffflag=1
c
c   Read trajectory header records if plotting trajectories.
c
      ifltrack=0
      if (cptyp(ipl)(2:2).eq.'t') then
         if (index(ctjfl(ipl),'fltrack').eq.0) then
            open (unit=iutrajin,file=ctjfl(ipl),form='unformatted',
     &            status='old')
            read (iutrajin)
            read (iutrajin) rtim,ctim,dttraj,ntrajplt
            ifltrack=0
         else
c
c         In flight track file, the following format issues are important:
c         The first 5 lines should look like this (excluding comment column):
cFlt NNNN - YYYY-MM-DD
c Time           tans-lat        tans-lon        tans-alt
c UTC            deg     deg     meters 
cHH:MM:SS          LL.LL         -LLL.LL          AAAA.A
c
c         Subsequent lines repeat the format of the 4th line above.
c         The important things for proper reading of data is that YYYY-MM-DD
c         (the UTC date of the first data point) should be in the first line
c         in the exact columns shown; and HH:MM:SS (UTC times of data points)
c         should be in the fourth and subsequent lines, in the exact columns
c         shown.  The lat, lon, and altitude data needn't be in those exact
c         columns, but should be space or tab delimited from each other and
c         from the time (they are read from the string in columns 9-31, using
c         list-directed I/O).  Data points should be every 1 second.
c         Missing values of lat, lon, or altitude should be -999.
c
            open (unit=iutrajin,file=ctjfl(ipl),form='formatted',
     &            status='old')
            read(iutrajin,'(12x,3(1x,i2))')nyear,nmonth,nday
            read(iutrajin,*)
            read(iutrajin,*)
            do i=1,1000
               read(iutrajin,'(3(i2,1x))')nh,nm,ns
               if (mod(ns,10).eq.0) then
                  nhour=nh
                  nmin=nm
                  nsec=ns
                  goto 280
               endif
            enddo
 280        iflt1=i
            do i=1,1000000
               read(iutrajin,'(3(i2,1x))',end=287)nh,nm,ns
               if (mod(ns,10).eq.0) iflt2=i
            enddo
 287        continue
            rewind (iutrajin)
            do i=1,3+iflt1-1
               read(iutrajin,*)
            enddo
            mdatetrajb=nyear*1000000+nmonth*10000+nday*100+nhour
            rhourtrajb=nmin/60.+nsec/3600.
            call mconvert(mdatetrajb,mhourtrajb,1,1940)
            rtim=float(mhourtrajb-mhourb)+rhourtrajb-rhourb
            dttraj=10. ! data avialable every second, but we're coarsening
c                      ! to every 10 seconds (about 1 km distance).
            ctim=rtim+float(iflt2-iflt1)/3600.
            ntrajplt=1
            ifltrack=1
         endif
         ntrajtimeplt=nint(abs(rtim-ctim)/dttraj*3600) + 1
      endif
c
c   Make title for individual plot underneath general title.
c
      if (cptyp(ipl)(2:2).ne.'b'.and.itrajcalc.eq.0.and.
     &    imakev5d.eq.0) then
         call pltitle(ctitl,cfeld,cptyp,rlevl,rlavl,cvcor,
     &      rstmv,lgrad,llapl,lhadv,rcrag,rcrbg,rslcg,ismth,
     &      ismcp,toptextclg,ilev,ixavg,icomg,csout,idimn,raddf,
     &      ifltrack,rtjst,rtjen,rtjti,titlestr,
     &      cdiff,rdiff,ldfrl,xtime,itjns,rtim,ctim,lnttl,lnsmm,lnvlb,
     &      ldiffsuccess,idescriptive,engplttl,maxlev,maxpl,mkzh,ipl)
      elseif (itrajcalc.eq.0.and.imakev5d.eq.0) then
         if (cfeld(1,ipl)(1:4).eq.'map ') then
            write(iup,*)'map background'
         elseif (cfeld(1,ipl)(1:4).eq.'tic ') then
            write(iup,*)'tic marks'
         elseif (cfeld(1,ipl)(1:4).eq.'box ') then
            write(iup,*)'box'
         elseif (cfeld(1,ipl)(1:5).eq.'line ') then
            write(iup,*)'line'
         elseif (cfeld(1,ipl)(1:7).eq.'bullet ') then
            write(iup,*)'bullet'
         elseif (cfeld(1,ipl)(1:5).eq.'vbar ') then
            write(iup,*)'vbar'
         elseif (cfeld(1,ipl)(1:5).eq.'sids ') then
            write(iup,*)'station ids'
         endif
      endif
c
      if (raddf(ipl).ne.0.0) goto 700
c
c   Put alternate vertical coord. in vc3d.
c
      vcncheck=cvcor(ipl)
      if (vcncheck.eq.'l'.or.vcncheck.eq.'x') vcncheck='p'
      if (vcncheck.ne.vc3dtyp) then
         if (vcncheck.eq.'s') then
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               vc3d(i,j,k)=float(k)
            enddo
            enddo
            enddo
         elseif (vcncheck.eq.'z') then
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               vc3d(i,j,k)=exp(-ght(i,j,k)/sclht)
            enddo
            enddo
            enddo
         elseif (vcncheck.eq.'p') then
            call addorfill(prs,vc3d,miy,mjx,mkzh,3,1,1.,0.)
         elseif (vcncheck.eq.'t') then
            call thecalc(prs,tmk,qvp,vc3d,miy,mjx,mkzh)
            call monotonic(vc3d,prs,1,.01,0,cor,miy,mjx,mkzh)
         elseif (vcncheck.eq.'e') then
            call eqthecalc(qvp,tmk,prs,vc3d,miy,mjx,mkzh)
            call monotonic(vc3d,prs,1,.01,0,cor,miy,mjx,mkzh)
         elseif (vcncheck.eq.'q') then
            ifree=incwk
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
            call thecalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
            call pvocalc(xmap,uuu,vvv,cor,scr3a,prs,
     &         vc3d,miy,mjx,mkzh)
            call monotonic(vc3d,prs,1,.001,1,cor,miy,mjx,mkzh)
         else
            write(iup,*)'Don''t understand vertical coordinate ',
     &         vcncheck
            stop
         endif
         vc3dtyp=vcncheck
      endif
c
c   Set up default values of vertical coordinate limits
c   for vertical cross sections (vtickinc is also used
c   for horizontal trajectory plots, and isense is used
c   for horizontal trajectory plots).
c
      if (cvcor(ipl).eq.'s') then ! these are in k-level
         vtickinc=1.
         defvv1=float(mkzh)
         defvv2=1.
         isense=-1  ! (k-index decreases w/ height)
      elseif (cvcor(ipl).eq.'p'.or.cvcor(ipl).eq.'l'.or.
     &        cvcor(ipl).eq.'x') then ! these are in hPa
         vtickinc=10.
         defvv1=1050.
         defvv2=100.
         isense=-1
      elseif (cvcor(ipl).eq.'z') then ! these are in km
         vtickinc=.1
         defvv1=0.
         defvv2=15.
         isense=1
      elseif (cvcor(ipl).eq.'t') then ! these are in K
         vtickinc=1.
         defvv1=260.
         defvv2=400.
         isense=1   ! theta increases w/ height
      elseif (cvcor(ipl).eq.'e') then ! these are in K
         vtickinc=1.
         defvv1=260.
         defvv2=400.
         isense=1
      elseif (cvcor(ipl).eq.'q') then ! these are in PVU
         vtickinc=.1
         defvv1=-.5
         defvv2=5.5
         isense=1
      endif
c
c   Next, calculate set limits for vertical coordinate
c
      if (rvwin(1,ipl).eq.rmsg) then
         vv1=defvv1
      else
         vv1=rvwin(1,ipl)
      endif
      if (rvwin(2,ipl).eq.rmsg) then
         vv2=defvv2
      else
         vv2=rvwin(2,ipl)
      endif
c
c   Adjust vv1,vv2 so that they are divisible by vtickinc,
c   then make set limits
c
      vv1=nint(vv1/vtickinc)*vtickinc
      vv2=nint(vv2/vtickinc)*vtickinc
      if (cvcor(ipl).eq.'l') then
         set1=alog(vv1)
         set2=alog(vv2)
      elseif (cvcor(ipl).eq.'x') then
         set1=(vv1)**gamma
         set2=(vv2)**gamma
      else
         set1=vv1
         set2=vv2
      endif
c
c   Do the trajectory calculation stuff if asked for.
c
      if (itrajcalc.eq.1) then
c
      if (ipltime.eq.1.and.ifr.eq.1) then
c
c      Create "sv" variables at initial time.  Regardless of what vertical
c      coordinate was used to specify the trajectory initial points, all
c      trajectory calculations will be done in terms of exponential
c      height [exp(-z/H)] (z in meters), and the vertical velocity used
c      will be d[exp(-z/H)]/dt.
c
         yicornsv=yicorn  !   Account for the possibility of a
         xjcornsv=xjcorn  !   moving domain.
c
         do k=1,mkzh
         do j=1,mjx
         do i=1,miy
            uusv(i,j,k)=uuu(i,j,k)
            vvsv(i,j,k)=vvv(i,j,k)
         enddo
         enddo
         enddo
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            ehsv(i,j,k)=exp(-ght(i,j,k)/sclht)
c           Note, www is in cm/s, but we want v.v. in m/s
            ehdotsv(i,j,k)=-.01*www(i,j,k)*ehsv(i,j,k)/sclht
         enddo
         enddo
         enddo
c
c      If ihydrometeor=1, include average hydrometeor vertical velocity,
c      i.e., particle vertical velocity = vertical air velocity plus
c      average particle terminal fall speed.
c
         if (ihydrometeor.eq.1) then
            ifree=incwk
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3e,wk,maxslab)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,
     &         ncxc,'qra       ',miy,mjx,mkzh,maxtavl,
     &         3,0,scr3c,istat)
            if (istat.eq.-1) call fillarray(scr3c,miy*mjx*mkzh,0.)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,
     &         ncxc,'qgr       ',miy,mjx,mkzh,maxtavl,
     &         3,0,scr3d,istat)
            if (istat.eq.-1) call fillarray(scr3d,miy*mjx*mkzh,0.)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,
     &         ncxc,'qsn       ',miy,mjx,mkzh,maxtavl,
     &         3,0,scr3e,istat)
            if (istat.eq.-1) call fillarray(scr3e,miy*mjx*mkzh,0.)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               vrain=0.0
               vgrap=0.0
               vsnow=0.0
               dens=prs(i,j,k)*100./
     &              (rgas*virtual(tmk(i,j,k),qvp(i,j,k))) ! in kg/m3
               if (scr3c(i,j,k).gt.0.0) 
     &              vrain=20.82*(dens*scr3c(i,j,k)*.001)**.2
               if (scr3d(i,j,k).gt.0.0)
     &              vgrap=3.964*(dens*scr3d(i,j,k)*.001)**.0925
               if (scr3e(i,j,k).gt.0.0)
     &              vsnow=1.987*(dens*scr3e(i,j,k)*.001)**.1025
               vdom=(scr3c(i,j,k)+scr3d(i,j,k)+scr3e(i,j,k))
               if (vdom.gt.0.0) then
                  vratr=scr3c(i,j,k)/vdom
                  vratg=scr3d(i,j,k)/vdom
                  vrats=scr3e(i,j,k)/vdom
                  vfall=vgrap*vratg + vrain*vratr + vsnow*vrats
               else
                  vfall=0.0
               endif
               ehdotsv(i,j,k) = ehdotsv(i,j,k) - vfall*ehsv(i,j,k)/sclht
            enddo
            enddo
            enddo
         endif
c
c      Convert initial yitraj and xjtraj values to coarse domain
c      coordinates, and convert initial zktraj values to
c      exponential height.
c
         do itr=1,ntraj
            if (xjtraj(itr).ne.rmsg) then
               yitraj(itr)=yicorn+(yitraj(itr)-1.)/refrat
               xjtraj(itr)=xjcorn+(xjtraj(itr)-1.)/refrat
            endif
            if (cvcor(1).eq.'z') then
               zktraj(itr)=exp(-1000.*zktraj(itr)/sclht) !zktraj was given in km
            else
               zktraj(itr)=finterp(vc3d,ehsv,1,miy,mjx,mkzh,yitraj(itr),
     &            xjtraj(itr),zktraj(itr),refrat,yicorn,xjcorn,rmsg,iup)
            endif
         enddo
c
c      Write out header information (standard rip header and trajectory
c         header) to trajectory position file.
c
         vardesc='trajectories'
         plchun=' '
         ihrip(6)=1  ! number of dimensions is not relevant
c                      to trajectories
         ihrip(7)=0  ! trajectory positions are always w.r.t
c                      dot-point grid
         ihrip(11)=0   ! mdate is not relevant for trajectory files
         rhrip(14)=0.  ! rhour is not relevant for trajectory files
         rhrip(15)=0.  ! xtime is not relevant for trajectory files
         write (iutrajout)vardesc,plchun,ihrip,rhrip,chrip
         write (iutrajout)rtim,ctim,dttraj,ntraj
c
c      Write out initial trajectory positions to file
c
         write(iutrajout)(yitraj(itr),itr=1,ntraj),
     &      (xjtraj(itr),itr=1,ntraj),(zktraj(itr),itr=1,ntraj)
c
c      Write out header information (standard rip header and trajectory
c         header) to trajectory diagnostics file.  Note: in trajcalc
c         mode, the variable nfr, which is the number of
c         frames specified in the plspec table, is also the number
c         if diagnostic variables to be calculated along trajectories.
c         Also, the variable dtfile is the relevant timestep for
c         diagnostics (instead of dttraj) because diagnostics are
c         calculated only at each model output time, not at all the
c         in-between times that trajectories are calculated at.
c
         vardesc='diagnostics'
         write(iudiagout)vardesc,plchun,ihrip,rhrip,chrip
         write(iudiagout)rtim,ctim,dtfile,ntraj,nfr
c
      endif
c
      if (ipltime.eq.1) then
c
c      Calculate and write out diagnostic quantity at first file time.
c
         do itr=1,ntraj
            if (xjtraj(itr).ne.rmsg) then
               if (idimn(ipl).eq.3) then
                  diag(itr)=finterp(ehsv,wk(1,1,indwk(1,ipl)),
     &               icdwk(ipl),miy,mjx,mkzh,yitraj(itr),xjtraj(itr),
     &               zktraj(itr),refrat,yicorn,xjcorn,rmsg,iup)
               else
                  diag(itr)=finterp2d(wk(1,1,indwk(1,ipl)),
     &               icdwk(ipl),miy,mjx,yitraj(itr),xjtraj(itr),
     &               refrat,yicorn,xjcorn,rmsg,iup)
               endif
            else
               diag(itr)=rmsg
            endif
         enddo
         write(iudiagout)(diag(itr),itr=1,ntraj)
c
      endif
c
      if (ipltime.gt.1.and.ifr.eq.1) then
c
c      Do time loop to integrate trajectory positions
c
c      Put current exponential height ("eh") into scr3a, and current
c      vertical velocity ("ehdot") into scr3b.
c
         ifree=incwk
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3a(i,j,k)=exp(-ght(i,j,k)/sclht)
c           Note, www is in cm/s, but we want v.v. in m/s
            scr3b(i,j,k)=-.01*www(i,j,k)*ehsv(i,j,k)/sclht
         enddo
         enddo
         enddo
c
c      If ihydrometeor=1, include average hydrometeor vertical velocity,
c      i.e., particle vertical velocity = vertical air velocity plus
c      average particle terminal fall speed.
c
         if (ihydrometeor.eq.1) then
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3e,wk,maxslab)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,
     &         ncxc,'qra       ',miy,mjx,mkzh,maxtavl,
     &         3,0,scr3c,istat)
            if (istat.eq.-1) call fillarray(scr3c,miy*mjx*mkzh,0.)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,
     &         ncxc,'qgr       ',miy,mjx,mkzh,maxtavl,
     &         3,0,scr3d,istat)
            if (istat.eq.-1) call fillarray(scr3d,miy*mjx*mkzh,0.)
            call getvar(casename,iendc,cxtimeavl,xtimeavl,nxt,
     &         ncxc,'qsn       ',miy,mjx,mkzh,maxtavl,
     &         3,0,scr3e,istat)
            if (istat.eq.-1) call fillarray(scr3e,miy*mjx*mkzh,0.)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               vrain=0.0
               vgrap=0.0
               vsnow=0.0
               dens=prs(i,j,k)*100./
     &              (rgas*virtual(tmk(i,j,k),qvp(i,j,k))) ! in kg/m3
               if (scr3c(i,j,k).gt.0.0) 
     &              vrain=20.82*(dens*scr3c(i,j,k)*.001)**.2
               if (scr3d(i,j,k).gt.0.0)
     &              vgrap=3.964*(dens*scr3d(i,j,k)*.001)**.0925
               if (scr3e(i,j,k).gt.0.0)
     &              vsnow=1.987*(dens*scr3e(i,j,k)*.001)**.1025
               vdom=(scr3c(i,j,k)+scr3d(i,j,k)+scr3e(i,j,k))
               if (vdom.gt.0.0) then
                  vratr=scr3c(i,j,k)/vdom
                  vratg=scr3d(i,j,k)/vdom
                  vrats=scr3e(i,j,k)/vdom
                  vfall=vgrap*vratg + vrain*vratr + vsnow*vrats
               else
                  vfall=0.0
               endif
               scr3b(i,j,k) = scr3b(i,j,k) - vfall*ehsv(i,j,k)/sclht
            enddo
            enddo
            enddo
         endif
c
c      Integrate the trajectories
c
         dsci=1./dsc
         do itm=1,ntrajtime
            fac1=(itm-1)*dttraj/dtfile
            fac2=itm*dttraj/dtfile
         do itr=1,ntraj
            if (xjtraj(itr).ne.rmsg) then
c
c            First iteration
c
               uutr2=finterp(scr3a,uuu,0,miy,mjx,mkzh,
     &            yitraj(itr),xjtraj(itr),zktraj(itr),
     &            refrat,yicorn,xjcorn,rmsg,iup)
               uutr1=finterp(ehsv,uusv,0,miy,mjx,mkzh,
     &            yitraj(itr),xjtraj(itr),zktraj(itr),
     &            refrat,yicornsv,xjcornsv,rmsg,iup)
               vvtr2=finterp(scr3a,vvv,0,miy,mjx,mkzh,
     &            yitraj(itr),xjtraj(itr),zktraj(itr),
     &            refrat,yicorn,xjcorn,rmsg,iup)
               vvtr1=finterp(ehsv,vvsv,0,miy,mjx,mkzh,
     &            yitraj(itr),xjtraj(itr),zktraj(itr),
     &            refrat,yicornsv,xjcornsv,rmsg,iup)
               ehdottr2=finterp(scr3a,scr3b,1,miy,mjx,mkzh,
     &            yitraj(itr),xjtraj(itr),zktraj(itr),
     &            refrat,yicorn,xjcorn,rmsg,iup)
               ehdottr1=finterp(ehsv,ehdotsv,1,miy,mjx,mkzh,
     &            yitraj(itr),xjtraj(itr),zktraj(itr),
     &            refrat,yicornsv,xjcornsv,rmsg,iup)
               if (uutr2.eq.rmsg.or.uutr1.eq.rmsg.or.
     &             vvtr2.eq.rmsg.or.vvtr1.eq.rmsg.or.
     &             ehdottr2.eq.rmsg.or.ehdottr1.eq.rmsg) then
                  xjtraj(itr)=rmsg
                  yitraj(itr)=rmsg
                  zktraj(itr)=rmsg
                  goto 169
               endif
               uutr=fac1*uutr2+(1.-fac1)*uutr1
               vvtr=fac1*vvtr2+(1.-fac1)*vvtr1
               ehdottr=fac1*ehdottr2+(1.-fac1)*ehdottr1
               xjnew=xjtraj(itr)+uutr*dsci*ntrajdir*dttraj
               yinew=yitraj(itr)+vvtr*dsci*ntrajdir*dttraj
               zknew=zktraj(itr)+ehdottr*ntrajdir*dttraj
c
c            Make sure zknew is within vertical domain.
c
               zknewlim=finterp(scr3a,scr3a,1,miy,mjx,mkzh,
     &            yinew,xjnew,zknew,refrat,yicorn,xjcorn,rmsg,iup)
               zknew=zknewlim
c
c            Second iteration
c
               uutr2=finterp(scr3a,uuu,0,miy,mjx,mkzh,
     &            yinew,xjnew,zknew,
     &            refrat,yicorn,xjcorn,rmsg,iup)
               uutr1=finterp(ehsv,uusv,0,miy,mjx,mkzh,
     &            yinew,xjnew,zknew,
     &            refrat,yicornsv,xjcornsv,rmsg,iup)
               vvtr2=finterp(scr3a,vvv,0,miy,mjx,mkzh,
     &            yinew,xjnew,zknew,
     &            refrat,yicorn,xjcorn,rmsg,iup)
               vvtr1=finterp(ehsv,vvsv,0,miy,mjx,mkzh,
     &            yinew,xjnew,zknew,
     &            refrat,yicornsv,xjcornsv,rmsg,iup)
               ehdottr2=finterp(scr3a,scr3b,1,miy,mjx,mkzh,
     &            yinew,xjnew,zknew,
     &            refrat,yicorn,xjcorn,rmsg,iup)
               ehdottr1=finterp(ehsv,ehdotsv,1,miy,mjx,mkzh,
     &            yinew,xjnew,zknew,
     &            refrat,yicornsv,xjcornsv,rmsg,iup)
               if (uutr2.eq.rmsg.or.uutr1.eq.rmsg.or.
     &             vvtr2.eq.rmsg.or.vvtr1.eq.rmsg.or.
     &             ehdottr2.eq.rmsg.or.ehdottr1.eq.rmsg) then
                  xjtraj(itr)=rmsg
                  yitraj(itr)=rmsg
                  zktraj(itr)=rmsg
                  goto 169
               endif
               uutr=.5*((fac2*uutr2+(1.-fac2)*uutr1)+uutr)
               vvtr=.5*((fac2*vvtr2+(1.-fac2)*vvtr1)+vvtr)
               ehdottr=.5*((fac2*ehdottr2+(1.-fac2)*ehdottr1)+ehdottr)
               xjtraj(itr)=xjtraj(itr)+uutr*dsci*ntrajdir*dttraj
               yitraj(itr)=yitraj(itr)+vvtr*dsci*ntrajdir*dttraj
               zktraj(itr)=zktraj(itr)+ehdottr*ntrajdir*dttraj
c
c            Make sure zktraj is within vertical domain.
c
               zktrajlim=finterp(scr3a,scr3a,1,miy,mjx,mkzh,
     &            yitraj(itr),xjtraj(itr),zktraj(itr),
     &            refrat,yicorn,xjcorn,rmsg,iup)
               zktraj(itr)=zktrajlim
 169           continue
            endif
         enddo
            write(iutrajout)(yitraj(itr),itr=1,ntraj),
     &         (xjtraj(itr),itr=1,ntraj),(zktraj(itr),itr=1,ntraj)
         enddo
c
c      Put current u,v,eh,ehdot into sv arrays
c
         do k=1,mkzh
         do j=1,mjx
         do i=1,miy
            uusv(i,j,k)=uuu(i,j,k)
            vvsv(i,j,k)=vvv(i,j,k)
            ehsv(i,j,k)=scr3a(i,j,k)
            ehdotsv(i,j,k)=scr3b(i,j,k)
         enddo
         enddo
         enddo
         yicornsv=yicorn
         xjcornsv=xjcorn
c   
      endif
c
      if (ipltime.gt.1) then
c
c      Calculate and write out diagnostic quantity
c
         do itr=1,ntraj
            if (xjtraj(itr).ne.rmsg) then
               if (idimn(ipl).eq.3) then
                  diag(itr)=finterp(ehsv,wk(1,1,indwk(1,ipl)),
     &               icdwk(ipl),miy,mjx,mkzh,yitraj(itr),xjtraj(itr),
     &               zktraj(itr),refrat,yicorn,xjcorn,rmsg,iup)
               else
                  diag(itr)=finterp2d(wk(1,1,indwk(1,ipl)),
     &               icdwk(ipl),miy,mjx,yitraj(itr),xjtraj(itr),
     &               refrat,yicorn,xjcorn,rmsg,iup)
               endif
            else
               diag(itr)=rmsg
            endif
         enddo
         write(iudiagout)(diag(itr),itr=1,ntraj)
c
      endif
c
      goto 700
c
      endif  ! end trajectory calculation stuff
c
c   Add variables to the Vis5d file
c
      if (imakev5d.eq.1) then
c
      ifree=incwk
      nwholelevsneeded=(nr*nc*nv5dlevels)/(miy*mjx)+1
      ifreenew=ifree+nwholelevsneeded
      if (ifreenew.gt.maxslab) then
         write(iup,*)'In creating Vis5d data array, the program is'
         write(iup,*)'trying to allocate more work space than is'
         write(iup,*)'available.  Try increasing maxfld and run'
         write(iup,*)'rip again.'
         stop
      endif
      i_v5darray=loc(wk(1,1,ifree))
      kend=nv5dlevels
      do kz=1,kend
         if (idimn(ipl).eq.3) then
            if (icdwk(ipl).eq.0) then
               niy=miy
               njx=mjx
            else
               niy=miy-1
               njx=mjx-1
            endif
            call vinterp('z',vert_args(kz),1,1,icdwk(ipl),vc3d,
     &         tmk,qvp,prs,ght,ter,sfp,lhide(ipl),
     &         idiffflag,cfeld(1,ipl),wk(1,1,indwk(1,ipl)),
     &         pslab1,mabpl,morpl,njx,niy,miy,mjx,mkzh)
         else
            do i=1,miy-icdwk(ipl)
            do j=1,mjx-icdwk(ipl)
               pslab1(j,i)=wk(i,j,indwk(1,ipl))
            enddo
            enddo
         endif
         if (icdwk(ipl).eq.0) then
            do i=1,miy-1
            do j=1,mjx-1
               if (pslab1(j,i).ne.rmsg.and.pslab1(j+1,i).ne.rmsg.and.
     &          pslab1(j,i+1).ne.rmsg.and.pslab1(j+1,i+1).ne.rmsg) then
                  pslab1(j,i)=.25*(pslab1(j,i)+pslab1(j+1,i)+
     &                             pslab1(j,i+1)+pslab1(j+1,i+1))
               else
                  pslab1(j,i)=rmsg
               endif
            enddo
            enddo
         endif
         do jrip=ixwin(1,1),ixwin(2,1)-1
            jv5d=jrip-ixwin(1,1)+1
         do irip=iywin(1,1),iywin(2,1)-1
            iv5d=iywin(2,1)-irip
            iv5delem=(kz-1)*(nr*nc)+(jv5d-1)*nr+iv5d
            v5darray(iv5delem)=pslab1(jrip,irip)
            if (v5darray(iv5delem).eq.rmsg) then
               v5darray(iv5delem)=MISSING
            endif
            if (varname(ifr).eq.'W         ') v5darray(iv5delem)=
     &         v5darray(iv5delem)*.01  ! Vis5d wants W in m/s
         enddo
         enddo
      enddo
      write(iup,*)'Writing to Vis5D data file: ',cfeld(1,ipl)
      iv5derr=v5dwrite(ipltime,ifr,v5darray)
      if (iv5derr.eq.0) then
         write(iup,*)'v5dwrite failed for some unknown reason.'
         write(iup,*)'Variable = ',cfeld(1,ipl)
         stop
      endif
      iv5dcount=iv5dcount+1
c
      goto 700
c
      endif  ! end of writing to vis5d file
c
c   Decide what to plot.
c
c PLOT DATA
c MGD begin mod
c Square up the window so that everything will rotate nicely by a quarter
c turn - this is not necessary if the data array is square, but we can't 
c assume that...
      if(rrota(ipl) .eq. 90. .or. rrota(ipl) .eq. -90.) then
        fbmino = fbmin
        ftmaxo = ftmax
        flmino = flmin
        frmaxo = frmax
        if(ixwin(2,ipl)-ixwin(1,ipl) .ge. 
     &     iywin(2,ipl)-iywin(1,ipl)) then
          vdif = ftmax - fbmin
          flmin = .5-vdif/2.
          frmax = .5+vdif/2.
        else
          vdif = frmax - flmin
          fbmin = (fbmino+ftmaxo)/2.-vdif/2.
          ftmax = (fbmino+ftmaxo)/2.+vdif/2.
        endif
      endif
c MGD end mod
      if (cptyp(ipl)(1:1).eq.'h') then           ! horizontal plot
         if (cfeld(1,ipl)(1:4).eq.'map ') then           ! map
            call hmapdraw(ngfbuf,ixwin,ixwingf,iywin,iywingf,
     &         yicorngf,xjcorngf,icolr,icolrgf,ilinw,ilinwgf,
     &         idash,idashgf,rtslb,rtslbgf,rcint,rcintgf,cmllm,
     &         cmllmgf,couty,coutygf,couds,coudsgf,ioulw,ioulwgf,
     &         iouco,ioucogf,imfco,imfcogf,iwhatgf,iam,xcs,ycs,
     &         rip_root,niam,ncs,maxbuf,maxpl,ipl,rrota,rrotagf)
         elseif (cfeld(1,ipl)(1:4).eq.'tic ') then       ! tic marks
            call hticdraw(ngfbuf,ilinw,ixwin,iywin,raxlg,raxtg,icolr,
     &         rtslb,ilinwgf,ixwingf,iywingf,raxlggf,raxtggf,icolrgf,
     &         rtslbgf,iwhatgf,maxbuf,maxpl,ipl,rrota,rrotagf)
         elseif (cfeld(1,ipl)(1:4).eq.'box ') then       ! box
            call hboxdraw(ilinw,ixwin,iywin,icolr,rcrag,rcrbg,
     &         maxpl,ipl,rrota)
         elseif (cfeld(1,ipl)(1:5).eq.'line ') then       ! line
            call hlinedraw(ilinw,ixwin,iywin,icolr,rcrag,rcrbg,
     &         maxpl,ipl,rrota)
         elseif (cfeld(1,ipl)(1:7).eq.'bullet ') then       ! bullet
            call hbulldraw(rtslb,ctitl,ixwin,iywin,icolr,rcrag,
     &         maxpl,ipl,rrota)
         elseif (cfeld(1,ipl)(1:5).eq.'sids ') then       ! station ids
            do nn = 1, nsids(ipl)
              csloc(1,ipl) = csids(nn,ipl)
              csloc(2,ipl) = 'missing             '
              call locinterp(csloc(1,ipl),gridx,gridy,
     &          rlat,rlon,iwmo,icaoid,locdesc,rip_root,rrota(ipl))
              rcrag(1,ipl) = gridy
              rcrag(2,ipl) = gridx
              call hsidsdraw(rtslb,ilinw,ixwin,iywin,icolr,rcrag,
     &              icaoid,maxpl,ipl,rrota)
            enddo
         elseif (cptyp(ipl)(2:2).eq.'h') then    ! characters
            call hchadraw(ixwin,iywin,rlevl,wk(1,1,indwk(1,ipl)),icdwk,
     &         ilev,iintv,icosq,rcosq,lchfl,icolr,incsq,pslab1,maxcosq,
     &         mabpl,morpl,maxlev,maxpl,miy,mjx,mkzh,ipl,rrota)
         elseif (cptyp(ipl)(2:2).eq.'c') then    ! contours
c            Save plot title for time series information.
            if (idotser.eq.1) then
               ntserv=ntserv+1
               if (ntserv .gt. maxtserv) then
                 write(iup,*) 'Number of time series variables exceeds'
                 write(iup,*) 'maxtserv. Increase maxtserv in driver.f'
                 write(iup,*) 'and recompile.'
                 stop
               endif
               if (ipltime.eq.1) tservname(ntserv)=titlestr
            endif
            call hcondraw(xtime,ilinw,vc3d,tmk,qvp,
     &         prs,ght,ter,sfp,
     &         ixwin,iywin,ismth,rcint,rcbeg,rcend,lmult,larng,
     &         idash,rlevl,rlavl,cnohl,lnolb,lnobr,lnozr,incon,
     &         bottextfloor,cfeld,cvcor,wk(1,1,indwk(1,ipl)),
     &         icdwk,unwk,ilev,icolr,icoll,ilcll,ilchl,rtslb,
     &         rtshl,imjsk,icomg,lnmsg,icong,iconl,icozr,idimn,
     &         lhide,lgrad,lhadv,ilwll,ilwng,ilwnl,ilwzr,idall,
     &         idang,idanl,idazr,ilcnl,ilczr,
     &         ilcbr,ipwlb,iorlb,ipwhl,ipwbr,ifclb,ifcnl,
     &         ifczr,ifchl,ilclo,ifclo,ccmth,rwdbr,ihvbr,
     &         idotser,tseryi,tserxj,tserdat,ntsers,ntsert,ntserv,
     &         icosq,rcosq,incsq,fred,fgreen,fblue,nco,icomax,pslab1,
     &         iam,xcs,ycs,niam,ncs,idiffflag,
     &         maxtserv,maxtsers,maxtsert,maxcosq,mabpl,morpl,
     &         maxlev,maxpl,miy,mjx,mkzh,ipl,rrota)
         elseif (cptyp(ipl)(2:2).eq.'v') then  ! vectors of <f1,f2>
            call hvecdraw(ilinw,vc3d,tmk,qvp,
     &         prs,ght,ter,sfp,
     &         icolr,ixwin,iywin,ismth,rvcmx,cfulb,unwk,lhide,icomg,
     &         iintv,rlevl,rlavl,cfeld,cvcor,idimn,idiffflag,
     &         wk(1,1,indwk(1,ipl)),wk(1,1,indwk(2,ipl)),icdwk,ilev,
     &         lnmsg,bottextfloor,pslab1,pslab2,mabpl,morpl,maxlev,
     &         maxpl,miy,mjx,mkzh,ipl,rrota)
         elseif (cptyp(ipl)(2:2).eq.'s') then    ! streamlines
            ifree=incwk
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
            call hstrdraw(ilinw,iintv,scr3a,vc3d,tmk,qvp,
     &         prs,ght,ter,sfp,
     &         ixwin,iywin,ismth,icolr,lhide,
     &         rlevl,rlavl,cfeld,cvcor,idimn,wk(1,1,indwk(1,ipl)),
     &         wk(1,1,indwk(2,ipl)),icdwk,ilev,idiffflag,
     &         pslab1,pslab2,mabpl,morpl,maxlev,maxpl,
     &         miy,mjx,mkzh,ipl,rrota)
         elseif (cptyp(ipl)(2:2).eq.'t') then    ! trajectories
            ifree=incwk
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               scr3a(i,j,k)=exp(-ght(i,j,k)/sclht)
            enddo
            enddo
            enddo
            call htrajdraw(ilinw,vc3d,scr3a,ixwin,iywin,
     &         ifltrack,idash,cnohl,lnolb,cvcor,icolr,icomg,
     &         ilcll,ilchl,rtslb,rtshl,
     &         icong,ilwng,idang,lnmsg,rvwin,rtjsp,itjns,itjid,itjni,
     &         rtjar,cfeld,rtjst,rtjen,rtjti,rstmv,rtim,ctim,dttraj,
     &         ntrajplt,ntrajtimeplt,vtickinc,isense,xtime,maxpl,
     &         miy,mjx,mkzh,ipl,rrota)
            close (iutrajin)
         endif
      elseif (cptyp(ipl)(1:1).eq.'v') then  ! vertical cross-section plots
c
c      First, calulate new vertical coordinate for x-sec, if needed
c
         if (rcrag(1,ipl).ne.rcragvc(1).or.
     &       rcrag(2,ipl).ne.rcragvc(2).or.
     &       rcrbg(1,ipl).ne.rcrbgvc(1).or.
     &       rcrbg(2,ipl).ne.rcrbgvc(2).or.
     &       cvcor(ipl).ne.vc2dtyp.or.
     &       lnogd(ipl).neqv.lnogdvc) then
            call vc2dcalc(vc3d,prs,cvcor(ipl),rcrag(1,ipl),
     &         rcrbg(1,ipl),nscrs,sfp,vcground,mabpl,miy,mjx,mkzh)
            rcragvc(1)=rcrag(1,ipl)
            rcragvc(2)=rcrag(2,ipl)
            rcrbgvc(1)=rcrbg(1,ipl)
            rcrbgvc(2)=rcrbg(2,ipl)
            vc2dtyp=cvcor(ipl)
            lnogdvc=lnogd(ipl)
         endif
c
c      OK, make vertical plot
c
         if (cfeld(1,ipl)(1:4).eq.'tic ') then           ! tic marks
            call vticdraw(ilinw,icolr,xseclen,raxtd,raxld,cvcor,
     &         rtslb,nscrs,raxlv,raxtv,lnogd,vcground,
     &         rlat,rlon,rlat2,rlon2,
     &         vv1,vv2,set1,set2,vtickinc,mabpl,maxpl,ipl)
         elseif (cfeld(1,ipl)(1:5).eq.'vbar ') then      ! vertical bar
            call vbardraw(ilinw,idash,icolr,set1,set2,xdist,ydist,
     &         rcrag,rslcg,xseclen,maxpl,ipl)
         elseif (cptyp(ipl)(2:2).eq.'c') then    ! contours
            call vcondraw(ilinw,rcrag,rcrbg,ismth,rcint,rcbeg,rcend,
     &         lmult,larng,idash,ixavg,cnohl,lnolb,lnobr,lnozr,icomg,
     &         incon,bottextfloor,wk(1,1,indwk(1,ipl)),
     &         icdwk,unwk,icolr,icoll,ilcll,
     &         ilchl,rtslb,rtshl,imjsk,lnmsg,cfeld,
     &         icong,iconl,icozr,ilwll,ilwng,ilwnl,ilwzr,idall,
     &         idang,idanl,idazr,ilcnl,ilczr,ihvbr,
     &         ilcbr,ipwlb,iorlb,ipwhl,ipwbr,ifclb,ifcnl,
     &         ifczr,ifchl,ilclo,ifclo,ccmth,idimn,rwdbr,
     &         nscrs,set1,set2,xdist,ydist,xseclen,icosq,
     &         rcosq,incsq,fred,fgreen,fblue,nco,icomax,
     &         pslab1,iam,xcs,ycs,niam,ncs,
     &         maxcosq,mabpl,morpl,maxpl,
     &         miy,mjx,mkzh,ipl)
         elseif (cptyp(ipl)(2:2).eq.'v') then
c         Vectors of paral. comp. of <f1,f2,f3>
c
c         Check to make sure the chosen form of vertical velocity
c         is consistent with the chosen vertical coordinate
c         for the cross section.
c
            if (cvcor(ipl).eq.'p') then
               if (cfeld(3,ipl)(1:3).ne.'omg'.and.
     &             cfeld(3,ipl)(1:5).ne.'qgomg'.and.
     &             cfeld(3,ipl)(1:5).ne.'qmomg'.and.
     &             cfeld(3,ipl)(1:5).ne.'seomf'.and.
     &             cfeld(3,ipl)(1:5).ne.'seomb'.and.
     &             cfeld(3,ipl)(1:5).ne.'smomf'.and.
     &             cfeld(3,ipl)(1:5).ne.'smomb') then
                  write(iup,*)'The only allowable vertical velocities',
     &               ' for vectors in a'
                  write(iup,*)'cross section with pressure as the',
     &               ' vertical coordinate are:'
                  write(iup,*)'   omg, qgomg, qmomg, seomf, seomb,'
                  write(iup,*)'   smomf, and smomb.'
                  stop
               endif
            elseif (cvcor(ipl).eq.'z') then
               if (cfeld(3,ipl)(1:3).ne.'www') then
                  write(iup,*)'The only allowable vertical velocities',
     &               ' for vectors in a'
                  write(iup,*)'cross section with height as the',
     &               ' vertical coordinate are:'
                  write(iup,*)'   www'
                  stop
               endif
            else
               write(iup,*)'Circulation vectors can only be drawn'
               write(iup,*)'in a cross section with height or'
               write(iup,*)'pressure as the vertical coordinate,'
               write(iup,*)'because those are the only two'
               write(iup,*)'vertical coordinates that have a'
               write(iup,*)'compatible form of vertical velocity'
               write(iup,*)'(www and omg) available in RIP.'
               stop
            endif
            call vvecdraw(ilinw,wk(1,1,indwk(3,ipl)),rcrag,rcrbg,
     &         ismth,icolr,rvcmx,rvvms,pavprof,ivvnx,ixavg,
     &         wk(1,1,indwk(1,ipl)),wk(1,1,indwk(2,ipl)),icdwk,icomg,
     &         lnmsg,nscrs,set1,set2,xdist,ydist,xseclen,cvcor,
     &         cfeld,vv1,vv2,bottextfloor,pslab1,pslab2,mabpl,morpl,
     &         maxpl,miy,mjx,mkzh,ipl)
         elseif (cptyp(ipl)(2:2).eq.'w') then
c           Vectors showing horiz. wind, defined by <f1,f2>
            call vwinddraw(ilinw,rcrag,rcrbg,ismth,unwk,icomg,
     &         icolr,rvcmx,cfulb,rvvms,pavprof,ivvnx,ixavg,
     &         wk(1,1,indwk(1,ipl)),wk(1,1,indwk(2,ipl)),icdwk,
     &         lnmsg,nscrs,set1,set2,xdist,ydist,xseclen,cvcor,
     &         vv1,vv2,bottextfloor,unorth,vnorth,pslab1,pslab2,
     &         mabpl,morpl,maxpl,miy,mjx,mkzh,ipl,rrota)
         elseif (cptyp(ipl)(2:2).eq.'t') then  ! trajectories
            ifree=incwk
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               scr3a(i,j,k)=exp(-ght(i,j,k)/sclht)
            enddo
            enddo
            enddo
            call vtrajdraw(ilinw,vc3d,scr3a,idash,ifltrack,
     &         cnohl,lnolb,cvcor,icolr,ilcll,ilchl,rtslb,rtshl,
     &         rtjsp,itjns,itjid,itjni,rtjar,cfeld,rtjst,
     &         rtjen,rtjti,rstmv,set1,set2,xdist,ydist,rcrag,xseclen,
     &         rtim,ctim,dttraj,ntrajplt,ntrajtimeplt,xtime,maxpl,
     &         miy,mjx,mkzh,ipl,rrota(ipl))
            close (iutrajin)
         endif
      elseif (cptyp(ipl)(1:1).eq.'s') then       ! skewt
         flminsou=0.  ! You can change these, but the positioning
         frmaxsou=1.  ! of the stats and hodograph information will
         fbminsou=0.  ! get screwed up. 
         ftmaxsou=.9
         if (cfeld(1,ipl)(1:4).eq.'tic ') then   ! skewt grid
            igray=igetcoind('light.gray',conam,nco)
            if (lplrs(ipl)) then
               call sticdraw_polar(icolr,igray,ilinw,ngfbuf,icolrgf,
     &            ilinwgf,lhodo,lmand,lsndg,lhodogf,lmandgf,lsndggf,
     &            iwhatgf,maxbuf,flminsou,frmaxsou,fbminsou,ftmaxsou,
     &            maxpl,ipl)
            else
               call sticdraw(icolr,igray,ilinw,ngfbuf,icolrgf,ilinwgf,
     &            lhodo,lmand,lsndg,lhodogf,lmandgf,lsndggf,iwhatgf,
     &            maxbuf,flminsou,frmaxsou,fbminsou,ftmaxsou,maxpl,ipl)
            endif
         elseif (cptyp(ipl)(2:2).eq.'c') then    ! contours
            if (lplrs(ipl)) then
               call scondraw_polar(ilinw,prs,rslcg,idash,icolr,icdwk,
     &            wk(1,1,indwk(1,ipl)),ipl,rslcgprv,prssou,
     &            flminsou,frmaxsou,fbminsou,ftmaxsou,maxpl,miy,mjx,
     &            mkzh)
            else
               call scondraw(ilinw,prs,rslcg,idash,icolr,icdwk,
     &            wk(1,1,indwk(1,ipl)),ipl,rslcgprv,prssou,
     &            flminsou,frmaxsou,fbminsou,ftmaxsou,maxpl,miy,mjx,
     &            mkzh)
            endif
            rslcgprv(1)=rslcg(1,ipl)
            rslcgprv(2)=rslcg(2,ipl)
         elseif (cptyp(ipl)(2:2).eq.'v') then     ! wind vectors
            if (lplrs(ipl)) then
               call svecdraw_polar(ilinw,prs,rslcg,icolr,icdwk,
     &            wk(1,1,indwk(1,ipl)),wk(1,1,indwk(2,ipl)),
     &            ipl,rslcgprv,unorth,vnorth,rvvms,pavprof,prssou,
     &            cfulb,lhodo,ins,flminsou,frmaxsou,fbminsou,ftmaxsou,
     &            icomax,maxpl,miy,mjx,mkzh)
            else
               call svecdraw(ilinw,prs,rslcg,icolr,icdwk,
     &            wk(1,1,indwk(1,ipl)),wk(1,1,indwk(2,ipl)),
     &            ipl,rslcgprv,unorth,vnorth,rvvms,pavprof,prssou,
     &            cfulb,lhodo,ins,flminsou,frmaxsou,fbminsou,ftmaxsou,
     &            icomax,maxpl,miy,mjx,mkzh)
            endif
            if (lbogs(ipl)) then
               call bogs(uuu,vvv,tmk,qvp,prs,ght,unorth,vnorth,ter,
     &                 mdate,rlat,rlon,rslcg,ipl,maxpl,miy,mjx,mkzh)
            endif
	    if (lsndg(ipl)) call sstats(ilinw,prs,rslcg,icolr,
     &         icdwk,uuu,vvv,tmk,qvp,
     &         ipl,rslcgprv,unorth,vnorth,prssou,
     &         ins,flminsou,frmaxsou,fbminsou,ftmaxsou,
     &         maxpl,miy,mjx,mkzh)
            rslcgprv(1)=rslcg(1,ipl)
            rslcgprv(2)=rslcg(2,ipl)
         endif
      elseif (cptyp(ipl)(1:2).eq.'pc') then   ! vertical profile contour
	 if (rcbeg(ipl) .gt. -rmsg) prfmin = rcbeg(ipl)
	 if (rcend(ipl) .lt.  rmsg) prfmax = rcend(ipl)
         call profil (wk(1,1,indwk(1,ipl)),maxpl,ipl,prfmax,prfmin,
     &      vc3d,cvcor(ipl),rslcg,set1,set2,unwk,csout,icomg,ilinw,
     &      miy,mjx,mkzh)
      endif
c MGD begin mod
c reset the window paramters if we have changed them
      if(rrota(ipl) .eq. 90. .or. rrota(ipl) .eq. -90.) then
        ftmax = ftmaxo
        fbmin = fbmino
        flmin = flmino
        frmax = frmaxo
      endif
c MGD end mod
c
  700 continue      ! End of plot loop
c
      if (itrajcalc.eq.0.and.imakev5d.eq.0) call frame
c
c   If horizontal plots, check if all plots were 2D.  If so, then jump
c   out of level loop after first iteration, just in case rip thinks it
c   needs to plot more than one level (because levl might have more than
c   one value, carried over from a previous frame).
c
      igetout=1
      do ipl=iplstrt,iplend
         if (cptyp(ipl)(1:1).eq.'h'.and.idimn(ipl).eq.3) igetout=0
      enddo
      if (igetout.eq.1) goto 950
c
  900 continue      ! End of level loop
  950 continue      ! End of frame loop
c
      if (icgmsplit.eq.1.and.
     &    itrajcalc.eq.0.and.imakev5d.eq.0.and.nfr.gt.0) then
c
c      Close the GFLASH workstation
c
         call gclwk (iwkidgf)
c
c      Close the metafile workstation and GKS
c
         call gdawk (iwkidcgm)
         call gclwk (iwkidcgm)
         call gclks
      endif
c
 1000 continue       ! End of time loop.
c
c   Write out the time series info.
c
      if (idotser.eq.1) then
         open (unit=iutserdat,file=rootname(1:iendcr)//'.tser',
     &      form='formatted',status='unknown')
         write(iutserdat,1100)ntsers,ntserv,ntsert
         write(iutserdat,*)
         do is=1,ntsers
            if (tserlocdesc(is)(1:10).eq.'          ') then
               write(iutserdat,'(a)') tserlocpoint(is)
            else
               write(iutserdat,1102) tserlocpoint(is),tserlocdesc(is)
            endif
            do iv=1,ntserv
               write(iutserdat,1104)tservname(iv)
               do it=1,ntsert
                  write(iutserdat,1106)mdatetser(it),fchour(it),
     &               rhourtser(it),tserdat(it,iv,is)
               enddo
            enddo
         enddo
      endif
c
 1100 format('There are ',i5,' locations, ',i3,' variables, and ',
     &        i3,' times.')
 1102 format(a,' (',a,')')
 1104 format(3x,a82)
 1106 format('      mdate = ',i8.8,','i2,', rhour = ',f7.4,
     &   '   value = ',e15.8)
c
      if (icgmsplit.ne.1.and.itrajcalc.eq.0.and.imakev5d.eq.0) then
c
c      Close the GFLASH workstation
c
         call gclwk (iwkidgf)
c
c      Close the metafile workstation and GKS
c
         call gdawk (iwkidcgm)
         call gclwk (iwkidcgm)
         call gclks
      endif
c
c   Close vis5d file
c
      if (imakev5d.eq.1) then
         if (iv5dcount.ne.numtimes*numvars) then
            write(iup,*)'Vis5d: number of times or number of'
            write(iup,*)'variables processed is not consistent'
            write(iup,*)'with what was expected.'
            write(iup,*)'numtimes,numvars,iv5dcount=',
     &         numtimes,numvars,iv5dcount
         endif
         iv5derr=v5dclose()
      endif
c
      write(iup,*)
      write(iup,*)'===================================='
      write(iup,*)' We''re outta here like Vladimir !! '
      write(iup,*)'===================================='
      return
      end
