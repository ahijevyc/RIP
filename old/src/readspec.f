c    
c Modification history -- Dave Ahijevych (MMM)
c
c
c 20020211:
c	added reals rfdmn and rfdmx
c
c 20020118:
c	added integer rmply
c		multiply next field by next field to this power
c		0.0	= do nothing
c		1.0	= multiply by next field
c		-1.0	= divide by next field
c		0.5	= multiply by square root of next field
c
c*********************************************************************c
c                                                                     c
      subroutine readspec(nfr,numppf,nltf,intim,rtime,cfeld,icomg,
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
     &   raxlg,raxld,raxlv,raxtg,raxtd,raxtv,rstrm,rrfst,raddf,rmply,
     &   rfdmn,rfdmx,
     &   cfulb,conam,incsq,sigh,cmllm,couty,iqgsm,iopenverf,
     &   ctjfl,ctitl,cv5nm,rtjsp,itjns,itjid,itjni,
     &   rtjar,rtjst,rtjen,rtjti,lgrad,llapl,lhadv,igdir,
     &   couds,ioulw,iouco,imfco,fred,fgreen,fblue,nco,
     &   csids,nsids,lnsmm,lnvlb,itrajcalc,imakev5d,
     &   maxcosq,maxfr,maxlev,maxpl,miy,mjx,mkzh,rrota,lplrs)
c
c   This subroutine reads information from the plot specification file
c   and transfers that information to the plot specification arrays.
c
      dimension numppf(maxfr),nltf(maxfr),intim(maxfr),incsq(maxpl),
     &   rtime(100,maxfr),ixwin(2,maxpl),iywin(2,maxpl),
     &   ilinw(maxpl),incon(maxpl),ptuse(100),
     &   icolr(maxpl),icoll(maxpl),ilcll(maxpl),ilchl(maxpl),
     &   rtslb(maxpl),rtshl(maxpl),imjsk(maxpl),icomg(maxpl),
     &   icosq(maxcosq,maxpl),rcosq(maxcosq,maxpl),rsepa(32,maxpl),
     &   rlevs(maxlev,maxpl),inlvs(maxpl),rlevl(maxlev,maxpl),
     &   rlavl(maxlev,maxpl),rvwin(2,maxpl),ixavg(maxpl),idash(maxpl),
     &   ismth(maxpl),ismcp(maxpl),iintv(maxpl),
     &   rvcmx(maxpl),ivvnx(maxpl),igdir(maxpl),
     &   icong(maxpl),iconl(maxpl),icozr(maxpl),ilwll(maxpl),
     &   ilwng(maxpl),ilwnl(maxpl),ilwzr(maxpl),idall(maxpl),
     &   idang(maxpl),idanl(maxpl),idazr(maxpl),
     &   ilcnl(maxpl),ilczr(maxpl),ilcbr(maxpl),ipwlb(maxpl),
     &   iorlb(maxpl),rdiff(maxpl),
     &   ipwhl(maxpl),ipwbr(maxpl),ifclb(maxpl),ifcnl(maxpl),
     &   ifczr(maxpl),ifchl(maxpl),ihvbr(maxpl),
     &   ilclo(maxpl),ifclo(maxpl),rwdbr(maxpl),raddf(maxpl),
     &   rmply(maxpl),rfdmn(maxpl),rfdmx(maxpl),
     &   rvvms(maxpl),rcint(maxpl),rstrm(2,maxpl),rrfst(4,maxpl),
     &   rcbeg(maxpl),rcend(maxpl),iqgsm(maxpl),
     &   raxlg(maxpl),
     &   raxld(maxpl),raxlv(maxpl),raxtg(maxpl),
     &   raxtd(maxpl),raxtv(maxpl),
     &   rtjsp(3,50,maxpl),itjns(maxpl),
     &   itjid(30,maxpl),itjni(maxpl),rtjar(2,maxpl),rtjst(maxpl),
     &   rtjen(maxpl),rtjti(maxpl),
     &   sigh(mkzh),ioulw(maxpl),iouco(maxpl),imfco(6,maxpl),
     &   fred(0:255),fgreen(0:255),fblue(0:255),rrota(maxpl)
      logical lnobr(maxpl),lnozr(maxpl),lnolb(maxpl),
     &   lchfl(maxpl),lmult(maxpl),larng(maxpl),lnmsg(maxpl),
     &   lnttl(maxpl),lverf(maxpl),lhide(maxpl),lnogd(maxpl),
     &   lredo(maxpl),lgrad(maxpl),llapl(maxpl),lhadv(maxpl),
     &   lhodo(maxpl),lmand(maxpl),lnsmm(maxpl),ldfrl(maxpl),
     &   lnvlb(maxpl),lsndg(maxpl),lnmin(maxpl),lbogs(maxpl),
     &   lplrs(maxpl)
      character cfeld(3,maxpl)*10, cptyp(maxpl)*2, cvcor(maxpl)*1,
     &   conam(0:255)*40,cfulb(maxpl)*5,ccmth(maxpl)*4,csave(maxpl)*10,
     &   ctjfl(maxpl)*256,ctitl(maxpl)*82,cv5nm(maxpl)*8,
     &   csloc(2,maxpl)*20,ccrsa(2,maxpl)*20,ccrsb(2,maxpl)*20,
     &   cmllm(maxpl)*5,couty(maxpl)*32,couds(maxpl)*5,
     &   csids(40,maxpl)*20,cdiff(maxpl)*256,cnohl(maxpl)*1
      dimension nsids(maxpl)
c
      dimension inlvl(5000)
      character string*500,c4*4,cospec*40,
     &   stringt*500,cvcorprev*1
c
      logical numeric
      external numeric
c
      include 'comconst'
c
c   First, find the colors defined as black and white (other than
c   the default foreground and background colors).  These will be used
c   later in the routine, in the section that sets the values of
c   icomg.
c
      igotblack=0
      igotwhite=0
      do i=2,nco
         if (fred(i).eq.0.0.and.fgreen(i).eq.0.0.and.
     &       fblue(i).eq.0.0) then
            igotblack=1
            icoindblack=i
         endif
         if (fred(i).eq.1.0.and.fgreen(i).eq.1.0.and.
     &       fblue(i).eq.1.0) then
            igotwhite=1
            icoindwhite=i
         endif
      enddo
      if (igotblack.eq.0.or.igotwhite.eq.0) then
         write(iup,*)'in readspec, couldn''t find black or white',
     &      ' color.'
         write(iup,*)'igotblack,igotwhite=',igotblack,igotwhite
      endif
c
c   Find correct location in input file.
c
 5    read(iuinput,'(a240)',end=500) stringt
      if (index(stringt,'=============').eq.0) goto 5
c
c   Skip next 2 lines of the plspec table.
c
      read(iuinput,*)
      read(iuinput,*)
      linec = 3  ! line counter
c
c   Begin read loop.
c
      ifr=0      ! frame number
      ipl=0      ! plot number
      ipltf=0    ! number of plots in this frame
      itimedef=0 ! time range not defined for 1st frame (use all times)
      cvcorprev='?'
      lenstr=len(string)
      lenstrt=len(stringt)
c
   15 string=' '
      istartrd=1
      inewrvwin=0
   20 stringt=' '
      read(iuinput,'(a240)',end=500) stringt
      linec = linec + 1
      if (stringt(1:32).eq.'                                ') goto 500
      if (index(stringt,'#').ne.0) goto 20
      do i=1,lenstrt   ! convert tabs to blanks
         if (ichar(stringt(i:i)).eq.9) stringt(i:i)=' '
      enddo
      call unblank(stringt)
      string(istartrd:lenstr)=stringt(1:lenstr-istartrd+1)
      indcont=index(stringt,'>')
      if (indcont.ne.0) then    ! continuation character found
         istartrd=istartrd+indcont-1
         goto 20
      endif
      c4=string(1:4)
c
      if (c4.eq.'time') then  ! time specifiers
         ipos=1
         itimedef=1
         ipos=ipos+5
         intim(ifr+1)=0
   22    intim(ifr+1)=intim(ifr+1)+1  !  number of time specifiers
         call getrnum(string,ipos,mkzh,rtime(intim(ifr+1),ifr+1))
         if (string(ipos-1:ipos-1).eq.',') goto 22
         goto 15
      endif
c
      if (c4.eq.'====') then
c
c      Check if empty frame.
c
         if (ipltf.eq.0) then
            write(iup,*)'Empty frame - check specification file.'
            goto 15
         endif
c
c      Do frame stuff.
c
         ifr=ifr+1
         numppf(ifr)=ipltf
         ipltf=0
         if (itimedef.eq.0) then
            intim(ifr)=nptuse
            do 23 itt=1,nptuse
               rtime(itt,ifr)=ptuse(itt)
   23       continue
         endif
         itimedef=0
         goto 15
c
      endif
c
      ipl=ipl+1
      if (ipl.gt.maxpl) then
         write(iup,*)'ipl exceeds maxpl.'
	 write(iup,*) 'ipl = ',ipl,' maxpl = ',maxpl
         stop
      endif
      ipltf=ipltf+1
c
c   Set defaults. Most of the defaults are reset with every plot, so
c      they MUST be respecified with each plot (i.e. they don't
c      carry over from the previous plot). Horizontal window,
c      vert.coor., level specifiers, cross-section endpoints, and
c      cross-section pressure window carry over
c
      if (ipl.eq.1) then
c MGD
         rrota(ipl)=0.
         ixwin(1,ipl)=1       ! left window limit (in x direction)
         ixwin(2,ipl)=mjx     ! right window limit (in x direction)
         iywin(1,ipl)=1       ! bottom window limit (in y direction)
         iywin(2,ipl)=miy     ! top window limit (in y direction)
         csloc(1,ipl)='3                   '
     &                        ! loc of sounding
         csloc(2,ipl)='3                   '
     &                        ! loc of sounding
         ccrsa(1,ipl)='3                   '
     &                        ! loc of cross-sec end pt. A (left)
         ccrsa(2,ipl)='3                   '
     &                        ! loc of cross-sec end pt. A (left)
         ccrsb(1,ipl)='10                  '
     &                        ! loc of cross-sec end pt. B (right)
         ccrsb(2,ipl)='10                  '
     &                        ! loc of cross-sec end pt. B (right)
         cvcor(ipl)='s'       ! v. coor. (s,p,l,x,z,t,e, or q)
         inlvs(ipl)=1         ! number of level specifiers
         rlevs(1,ipl)=float(mkzh) ! level specifiers
         rvwin(1,ipl)=rmsg    ! lower window limit in vertical dir.
         rvwin(2,ipl)=rmsg    ! upper window limit in vertical dir.
      else
         ixwin(1,ipl)=ixwin(1,ipl-1)
         ixwin(2,ipl)=ixwin(2,ipl-1)
         iywin(1,ipl)=iywin(1,ipl-1)
         iywin(2,ipl)=iywin(2,ipl-1)
         csloc(1,ipl)=csloc(1,ipl-1)
         csloc(2,ipl)=csloc(2,ipl-1)
         ccrsa(1,ipl)=ccrsa(1,ipl-1)
         ccrsa(2,ipl)=ccrsa(2,ipl-1)
         ccrsb(1,ipl)=ccrsb(1,ipl-1)
         ccrsb(2,ipl)=ccrsb(2,ipl-1)
         cvcor(ipl)=cvcor(ipl-1)
         inlvs(ipl)=inlvs(ipl-1)
         do 25 ilev=1,inlvs(ipl)
            rlevs(ilev,ipl)=rlevs(ilev,ipl-1)
   25    continue
         rvwin(1,ipl)=rvwin(1,ipl-1)
         rvwin(2,ipl)=rvwin(2,ipl-1)
      endif
      cfeld(1,ipl)='dm1       '      ! plot field 1
      cfeld(2,ipl)='dm2       '      ! plot field 2
      cfeld(3,ipl)='dm3       '      ! plot field 3
      cptyp(ipl)='xx'         ! plot type ('xx' indicates no type requested)
      ilinw(ipl)=1            ! line width, normal (thinnest) = 1
      lnobr(ipl)=.false.      ! no label bar for filled contour plots
      rwdbr(ipl)=-1.          ! width of strip in which to confine label bar
c                             !   (incl. labels) (-1: let rip figure it out)
      ihvbr(ipl)=-1           ! label bar orient. flag: -1 = let rip choose;
c                             !   0=hor; 1=ver
      lnozr(ipl)=.false.      ! zero contour suppressed
      cnohl(ipl)=' '          ! highs and/or lows suppressed
      lnolb(ipl)=.false.      ! contour labels suppressed
      lnmsg(ipl)=.false.      ! no message for contour, vector plots
      lnmin(ipl)=.false.      ! no model info message at bottom of frame
      lhide(ipl)=.false.      ! don't show stuff in below-ground areas
      lredo(ipl)=.false.      ! force calculation of the field
      lnogd(ipl)=.false.      ! show pseudo-ground sfc in x-secs (p-lev data)
      lnttl(ipl)=.false.      ! no plot title
      lverf(ipl)=.false.      ! verify the field against observations
      lgrad(ipl)=.false.      ! calc. and plot hor. grad. of the field
      llapl(ipl)=.false.      ! calc. and plot hor. lapl. of the field
      lhadv(ipl)=.false.      ! calc. and plot hor. advec. of the field
      cmllm(ipl)='none '      ! l/l mask-out: 'none' 'land' or 'water'
      couty(ipl)='PS'         ! map outline set (for 'OU' in EZMAP), or
c                             ! if not set to acceptable value for 'OU',
c                             ! use hires map info from this file in rip_root
      couds(ipl)='dot  '      ! map outlines: dot or solid
      csave(ipl)='dontsave  ' ! name of file for saving field
      ioulw(ipl)=999999       ! outl wth or dot spc (999999 -> 1 or 12)
      iouco(ipl)=999999       ! map outl color (999999 -> same as colr)
      do i=1,6
         imfco(i,ipl)=999999  ! map fill colors (999999 -> don't fill)
      enddo
      do i=1,40
         csids(i,ipl)=' '     ! station IDs to be shown on plot
      enddo
      nsids(ipl)=0
      ixavg(ipl)=0            ! num. of gr. pts. to average x-secs
      idash(ipl)=70           ! dash pattern
      ismth(ipl)=0            ! number of smoothing passes
      ismcp(ipl)=0            ! number of constant-prs smoothing passes
      iqgsm(ipl)=0            ! # smth passes in qgomg.f (p-lev data)
      iintv(ipl)=1            ! grid interv. betw. hor. vecs. or chars.
      rvcmx(ipl)=0.           ! mag. of vector of length intv
      ivvnx(ipl)=20           ! approx. # of vert. vecs. in hor. direc.
      rvvms(ipl)=0.0          ! min. delta-sigma betw. vert. vecs.
      cfulb(ipl)='5mps '      ! magnitude of full barb
      rcint(ipl)=rmsg         ! contour interval
      rcbeg(ipl)=-rmsg        ! starting contour
      rcend(ipl)=rmsg         ! ending contour
      lmult(ipl)=.false.      ! multiplicative contouring flag
      larng(ipl)=.false.      ! auto colr range from min to max contour
      lhodo(ipl)=.false.      ! plot a hodograph on a skewt
      lbogs(ipl)=.false.      ! print out bogus sndg in little-r format
      lmand(ipl)=.false.      ! plot lines at mandatory levels on sndg bkg
      lsndg(ipl)=.false.      ! plot sounding analysis parameters on a skewt
      lchfl(ipl)=.false.      ! fill boxes instead of printing values
      ccmth(ipl)='cont'       ! 'cont', 'fill', or 'both'
      incon(ipl)=20           ! nice number of contours
      icolr(ipl)=1            ! main foreground color
      lnsmm(ipl)=.false.      ! no smoothing message on the plot title
      lnvlb(ipl)=.false.      ! no vertical level information on the plot title
      do 28 i=1,maxcosq
         icosq(i,ipl)=1       ! colors for color sequence
         rcosq(i,ipl)=0.      ! values for color sequence
   28 continue
      incsq(ipl)=0            ! number of ramp colors
      cdiff(ipl)='none'       ! data set to subtract (model run differencing)
      rdiff(ipl)=rmsg         ! data time to subtract (time differencing)
      ldfrl(ipl)=.false.      ! rdiff is rel. to curr. time instead of absolute
      icozr(ipl)=999999       ! color for zero contour
      ilczr(ipl)=999999       ! label color for zero contour
      icoll(ipl)=999999       ! color for labeled contours
      ilcll(ipl)=999999       ! label color for labeled contours
      icong(ipl)=999999       ! color for negative unlabeled contours
      iconl(ipl)=999999       ! color for negative labeled contours
      ilcnl(ipl)=999999       ! label color for neg labeled contours
      ilchl(ipl)=999999       ! label color for hi/lo
      ilclo(ipl)=999999       ! label color for lows
      ilwzr(ipl)=999999       ! line width for zero contour
      ilwll(ipl)=999999       ! line width for labeled contours
      ilwng(ipl)=999999       ! line width for neg unlabeled contours
      ilwnl(ipl)=999999       ! line width for neg labeled contours
      idazr(ipl)=999999       ! dash pattern for zero contour
      idall(ipl)=999999       ! dash pattern for labeled contours
      idang(ipl)=999999       ! dash pattern for neg unlabeled contours
      idanl(ipl)=999999       ! dash pattern for neg labeled contours
      iorlb(ipl)=1            ! orientation of contour labels
      ipwlb(ipl)=1            ! perim width of contour and hi/lo labels
      ipwhl(ipl)=999999       ! perimeter width for hi/lo label boxes
      ifclb(ipl)=999999       ! fill color for labels
      ifcnl(ipl)=999999       ! fill color for neg contour label boxes
      ifczr(ipl)=999999       ! fill color for zero contour label boxes
      ifchl(ipl)=999999       ! fill color for hi/lo label boxes
      ifclo(ipl)=999999       ! fill color for low label boxes
      rtslb(ipl)=.01          ! text sz of contour and high/low labels
      rtshl(ipl)=rmsg         ! text sz for hi/los
      ilcbr(ipl)=1            ! color of label bar labs and perim lines
      ipwbr(ipl)=1            ! width of label bar perimeter lines
      imjsk(ipl)=3            ! # of unlab. conts. bet. lab. conts.
      raxlg(ipl)=rmsg         ! intv betw labld ticks (in gr pts)
      raxld(ipl)=rmsg         ! intv betw labld ticks (in km)
      raxlv(ipl)=rmsg         ! intv betw labld ticks (in km,hPa,K,sig)
      raxtg(ipl)=rmsg         ! intv betw small ticks (in gr pts)
      raxtd(ipl)=rmsg         ! intv betw small ticks (in km)
      raxtv(ipl)=rmsg         ! intv betw small ticks (in km,hPa,K,sig)
      igdir(ipl)=362          ! direc. of grad., if "grad" or "lapl" are set
      rstrm(1,ipl)=rmsg       ! y-vel. of storm in m/s (for rel. wind)
      rstrm(2,ipl)=rmsg       ! x-vel. of storm in m/s (for rel. wind)
      rrfst(1,ipl)=0.0        ! ref. state SLP (Pa) for pert. field plots
      rrfst(2,ipl)=0.0        ! ref. state SLT (K) for pert. field plots
      rrfst(3,ipl)=0.0        ! ref. state [dT/d(ln p)] for pert plots
      rrfst(4,ipl)=0.0        ! ref. state strat. T (K) for pert plots
      raddf(ipl)=0.0          ! mult fctr for adding this field to next
      rmply(ipl)=0.0          ! power fctr for multiplying this field by next
      rfdmn(ipl)=rmsg         ! lower limit of field values
      rfdmx(ipl)=rmsg         ! upper limit of field values
      rsepa(1,ipl)=200.       ! imaxlit parameter for Saw.-El.
      rsepa(2,ipl)=1.         ! errmin parameter for Saw.-El.
      rsepa(3,ipl)=1.6        ! alphor parameter for Saw.-El.
      rsepa(4,ipl)=0.         ! ixaverage parameter for Saw.-El.
      rsepa(5,ipl)=0.         ! smfac parameter for Saw.-El.
      rsepa(6,ipl)=15.        ! imaxbig parameter for Saw.-El.
      rsepa(7,ipl)=80.        ! rhithresh parameter for Saw.-El.
      rsepa(8,ipl)=2.0        ! smfstb parameter for Saw.-El.
      rsepa(9,ipl)=1.0        ! switch for SG (1) or QG (2) version of LHS
      do i=10,32
         rsepa(i,ipl)=rmsg
      enddo
      ctjfl(ipl)='junk'       ! name trajec. position info. file
      ctitl(ipl)='auto'       ! title for each plot ('auto': rip makes it)
      cv5nm(ipl)='samasvar'   ! var. name for Vis5D; samasvar-same as feld
      do i=1,50
      do j=1,3
         rtjsp(j,i,ipl)=rmsg  ! storm pos for trj, j=1,2,3-> xtime,y,x
      enddo
      enddo
      itjns(ipl)=0            ! # of storm pos vals (0-> grnd-rel trjs)
      itjid(1,ipl)=1          ! ID numbers of trajectories to plot
      itjni(ipl)=0            ! number of values in itjid
      rtjar(1,ipl)=.003       ! arrow width at bottom of vert. window
      rtjar(2,ipl)=.035       ! arrow width at top of vert. window
      rtjst(ipl)=rmsg         ! trajectory start time (xtime in hours)
      rtjen(ipl)=rmsg         ! trajectory end time (xtime in hours)
      rtjti(ipl)=1.0          ! traj time intv (hrs., for arrw hds or swarms)

                              ! "plrs" stands for "polar skew-t".  With
                              ! this logical flag set in the plot 
                              ! specification for a sounding, the skew-t 
                              ! (and background) will be shifted 30 K, so
                              ! that polar temperatures will fit on the
                              ! background.
      lplrs(ipl) = .false.    ! by default, skew-t does not have polar offset

c
c   Assign plot characteristics.
c
      ipos=1
 30   c4=string(ipos:ipos+3)
      iopenverf=0
c
      if (c4.eq.'feld') then
         ipos=ipos+5
         call getchar(string,ipos,cfeld(1,ipl),0)
         if (string(ipos-1:ipos-1).eq.',') then
            call getchar(string,ipos,cfeld(2,ipl),0)
         endif
         if (string(ipos-1:ipos-1).eq.',') then
            call getchar(string,ipos,cfeld(3,ipl),0)
         endif
c MGD begin mod
      elseif (c4.eq.'rota') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,rrota(ipl))
c MGD end mod
c KWM begin mod
      elseif (c4.eq.'plrs') then
         ipos=ipos+5
         lplrs(ipl) = .TRUE.
c KWM end mod
      elseif (c4.eq.'ptyp') then
         ipos=ipos+5
         cptyp(ipl)=string(ipos:ipos+1)
         ipos=ipos+3
      elseif (c4.eq.'xwin') then
         ipos=ipos+5
         call getinum(string,ipos,ixwin(1,ipl))
         call getinum(string,ipos,ixwin(2,ipl))
      elseif (c4.eq.'ywin') then
         ipos=ipos+5
         call getinum(string,ipos,iywin(1,ipl))
         call getinum(string,ipos,iywin(2,ipl))
      elseif (c4.eq.'crsa') then
         ipos=ipos+5
         nterm=0
   33    nterm=nterm+1
         if (nterm.gt.2) then
            write(iup,*)'Too many terms for crsa.'
            write(iup,*)'ipl=',ipl
            stop
         endif
         call getchar(string,ipos,ccrsa(nterm,ipl),0)
         if (string(ipos-1:ipos-1).eq.',') goto 33
         if (nterm.eq.1) ccrsa(2,ipl)='missing             '
      elseif (c4.eq.'crsb') then
         ipos=ipos+5
         nterm=0
   35    nterm=nterm+1
         if (nterm.gt.2) then
            write(iup,*)'Too many terms for crsb.'
            write(iup,*)'ipl=',ipl
            stop
         endif
         call getchar(string,ipos,ccrsb(nterm,ipl),0)
         if (string(ipos-1:ipos-1).eq.',') goto 35
         if (nterm.eq.1) ccrsb(2,ipl)='missing             '
      elseif (c4.eq.'sloc') then
         ipos=ipos+5
         nterm=0
   37    nterm=nterm+1
         if (nterm.gt.2) then
            write(iup,*)'Too many terms for sloc.'
            write(iup,*)'ipl=',ipl
            stop
         endif
         call getchar(string,ipos,csloc(nterm,ipl),0)
         if (string(ipos-1:ipos-1).eq.',') goto 37
         if (nterm.eq.1) csloc(2,ipl)='missing             '
      elseif (c4.eq.'sids') then
         ipos=ipos+5
         nterm=0
  137    nterm=nterm+1
         if (nterm.gt.20) then
            write(iup,*)'Too many terms for idns.'
            write(iup,*)'ipl=',ipl
            stop
         endif
         call getchar(string,ipos,csids(nterm,ipl),0)
         nsids(ipl) = nterm
         if (string(ipos-1:ipos-1).eq.',') goto 137
      elseif (c4.eq.'linw') then
         ipos=ipos+5
         call getinum(string,ipos,ilinw(ipl))
      elseif (c4.eq.'nobr') then
         ipos=ipos+5
         lnobr(ipl)=.true.
      elseif (c4.eq.'wdbr') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,rwdbr(ipl))
      elseif (c4.eq.'hvbr') then
         ipos=ipos+5
         call getinum(string,ipos,ihvbr(ipl))
      elseif (c4.eq.'nozr') then
         ipos=ipos+5
         lnozr(ipl)=.true.
      elseif (c4.eq.'nohl') then
         ipos=ipos+5
         if (string(ipos-1:ipos-1).eq.';'.or.
     &       string(ipos-1:ipos-1).eq.' ') then
            cnohl(ipl)='B'  ! for "Both H and L labels suppressed"
         elseif (string(ipos-1:ipos).eq.'=;'.or.
     &           string(ipos-1:ipos).eq.'= ') then
            cnohl(ipl)='B'  ! for "Both H and L labels suppressed"
            ipos=ipos+1
         else
            call getchar(string,ipos,cnohl(ipl),0)
         endif
      elseif (c4.eq.'nolb') then
         ipos=ipos+5
         lnolb(ipl)=.true.
      elseif (c4.eq.'nmsg') then
         ipos=ipos+5
         lnmsg(ipl)=.true.
      elseif (c4.eq.'nmin') then
         ipos=ipos+5
         lnmin(ipl)=.true.
      elseif (c4.eq.'nsmm') then
         ipos=ipos+5
         lnsmm(ipl)=.true.
      elseif (c4.eq.'nvlb') then
         ipos=ipos+5
         lnvlb(ipl)=.true.
      elseif (c4.eq.'hide') then
         ipos=ipos+5
         lhide(ipl)=.true.
      elseif (c4.eq.'redo') then
         ipos=ipos+5
         lredo(ipl)=.true.
      elseif (c4.eq.'nogd') then
         ipos=ipos+5
         lnogd(ipl)=.true.
      elseif (c4.eq.'nttl') then
         ipos=ipos+5
         lnttl(ipl)=.true.
      elseif (c4.eq.'verf') then
         ipos=ipos+5
         lverf(ipl)=.true.
         iopenverf=1
      elseif (c4.eq.'vcor') then
         ipos=ipos+5
         cvcor(ipl)=string(ipos:ipos)
         ipos=ipos+2
      elseif (c4.eq.'levs') then
         ipos=ipos+5
         inlvs(ipl)=0
   40    inlvs(ipl)=inlvs(ipl)+1
         call getrnum(string,ipos,mkzh,rlevs(inlvs(ipl),ipl))
         if (string(ipos-1:ipos-1).eq.',') goto 40
      elseif (c4.eq.'sepa') then
         ipos=ipos+5
         nsepa=0
   44    nsepa=nsepa+1
         if (nsepa.gt.32) then
            write(iup,*)'Max allowed values for sepa is 32.'
            write(iup,*)'Currently, only 8 actually get used.'
            write(iup,*)'Stopping.'
            stop
         endif
         if (string(ipos:ipos).eq.','.or.string(ipos:ipos).eq.';'.or.
     &       string(ipos:ipos).eq.' ') then
            ipos=ipos+1
         else
            call getrnum(string,ipos,mkzh,rsepa(nsepa,ipl))
         endif
         if (string(ipos-1:ipos-1).eq.',') goto 44
      elseif (c4.eq.'strm') then
         ipos=ipos+5
         nstrm=0
 41      nstrm=nstrm+1
         if (nstrm.gt.2) then
            write(iup,*)'No more than 2 values should be supplied'//
     &      ' for storm velocity.'
            stop
         endif
         call getrnum(string,ipos,mkzh,rstrm(nstrm,ipl))
         if (string(ipos-1:ipos-1).eq.',') goto 41
         if (nstrm.eq.2) then
c
c         Values were given in the order u_strm,v_strm, but should be
c         stored as v_strm,u_strm.
c
            rstrmsv=rstrm(1,ipl)
            rstrm(1,ipl)=rstrm(2,ipl)
            rstrm(2,ipl)=rstrmsv
         endif
      elseif (c4.eq.'rfst') then
         ipos=ipos+5
         nrfst=0
 42      nrfst=nrfst+1
         call getrnum(string,ipos,mkzh,rrfst(nrfst,ipl))
         if (string(ipos-1:ipos-1).eq.',') goto 42
         if (nrfst.ne.4) then
            write(iup,*) 'Exactly 4 values should be supplied',
     &                   ' for reference state:'
            write(iup,*) 'slp, slt, lapse, and strat. T.  You ',
     &                   'supplied ',nrfst,'.'
            stop
         endif
      elseif (c4.eq.'vwin') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,rvwin(1,ipl))
         if (string(ipos-1:ipos-1).eq.',') then
            call getrnum(string,ipos,mkzh,rvwin(2,ipl))
         endif
         inewrvwin=1
      elseif (c4.eq.'oulw') then
         ipos=ipos+5
         call getinum(string,ipos,ioulw(ipl))
      elseif (c4.eq.'ouco') then
         ipos=ipos+5
         call getchar(string,ipos,cospec,0)
         iouco(ipl)=igetcoind(cospec,conam,nco)
      elseif (c4.eq.'mfco') then
         ipos=ipos+5
         nmfco=0
 43      nmfco=nmfco+1
         if (nmfco.gt.6) then
            write(iup,*)'You are specifying more than 6 map fill',
     &         ' colors.'
            stop
         endif
         call getchar(string,ipos,cospec,0)
         imfco(nmfco,ipl)=igetcoind(cospec,conam,nco)
         if (string(ipos-1:ipos-1).eq.',') goto 43
      elseif (c4.eq.'xavg') then
         ipos=ipos+5
         call getinum(string,ipos,ixavg(ipl))
      elseif (c4.eq.'dash') then
         ipos=ipos+5
         call getinum(string,ipos,idash(ipl))
      elseif (c4.eq.'smth') then
         ipos=ipos+5
         call getinum(string,ipos,ismth(ipl))
      elseif (c4.eq.'smcp') then
         ipos=ipos+5
         call getinum(string,ipos,ismcp(ipl))
      elseif (c4.eq.'qgsm') then
         ipos=ipos+5
         call getinum(string,ipos,iqgsm(ipl))
      elseif (c4.eq.'intv') then
         ipos=ipos+5
         call getinum(string,ipos,iintv(ipl))
      elseif (c4.eq.'vcmx') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,rvcmx(ipl))
      elseif (c4.eq.'vvnx') then
         ipos=ipos+5
         call getinum(string,ipos,ivvnx(ipl))
      elseif (c4.eq.'vvms') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,rvvms(ipl))
      elseif (c4.eq.'cint') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,rcint(ipl))
      elseif (c4.eq.'cbeg') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,rcbeg(ipl))
      elseif (c4.eq.'cend') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,rcend(ipl))
      elseif (c4.eq.'mult') then
         ipos=ipos+5
         lmult(ipl)=.true.
      elseif (c4.eq.'arng') then
         ipos=ipos+5
         larng(ipl)=.true.
      elseif (c4.eq.'hodo') then
         ipos=ipos+5
         lhodo(ipl)=.true.
      elseif (c4.eq.'bogs') then
         ipos=ipos+5
         lbogs(ipl)=.true.
      elseif (c4.eq.'mand') then
         ipos=ipos+5
         lmand(ipl)=.true.
      elseif (c4.eq.'sndg') then
         ipos=ipos+5
         lsndg(ipl)=.true.
      elseif (c4.eq.'chfl') then
         ipos=ipos+5
         lchfl(ipl)=.true.
      elseif (c4.eq.'grad') then
         ipos=ipos+5
         lgrad(ipl)=.true.
      elseif (c4.eq.'lapl') then
         ipos=ipos+5
         llapl(ipl)=.true.
      elseif (c4.eq.'hadv') then
         ipos=ipos+5
         lhadv(ipl)=.true.
      elseif (c4.eq.'ncon') then
         ipos=ipos+5
         call getinum(string,ipos,incon(ipl))
      elseif (c4.eq.'mllm') then
         ipos=ipos+5
         call getchar(string,ipos,cmllm(ipl),0)
      elseif (c4.eq.'outy') then
         ipos=ipos+5
         call getchar(string,ipos,couty(ipl),0)
      elseif (c4.eq.'ouds') then
         ipos=ipos+5
         call getchar(string,ipos,couds(ipl),0)
      elseif (c4.eq.'save') then
         ipos=ipos+5
         if (string(ipos-1:ipos-1).eq.';'.or.
     &       string(ipos-1:ipos-1).eq.' ') then
            csave(ipl)='sameasfeld'
         else
            call getchar(string,ipos,csave(ipl),0)
         endif
      elseif (c4.eq.'fulb') then
         ipos=ipos+5
         call getchar(string,ipos,cfulb(ipl),0)
      elseif (c4.eq.'colr') then
         ipos=ipos+5
         call getchar(string,ipos,cospec,0)
         icolr(ipl)=igetcoind(cospec,conam,nco)
      elseif (c4.eq.'cosq') then
         ipos=ipos+5
         incsq(ipl)=0
   45    incsq(ipl)=incsq(ipl)+1
         call getrnum(string,ipos,mkzh,rcosq(incsq(ipl),ipl))
         call getchar(string,ipos,cospec,0)
         icosq(incsq(ipl),ipl)=igetcoind(cospec,conam,nco)
         if (string(ipos-1:ipos-1).eq.',') goto 45
      elseif (c4.eq.'diff') then
         ipos=ipos+5
 48      do i=ipos,lenstr
            if (string(i:i).eq.';'.or.string(i:i).eq.','.or.
     &            string(i:i).eq.' ') then
               goto 220
            endif
         enddo
  220    ilast=i-1
         if (numeric(string(ipos:ilast))) then
            call getrnum(string,ipos,mkzh,rdiff(ipl))
         elseif (ilast-ipos.eq.2.and.string(ipos:ilast).eq.'rel') then
            ldfrl(ipl)=.true.
            ipos=ilast+2
         else
            call getchar(string,ipos,cdiff(ipl),0)
         endif
         if (string(ipos-1:ipos-1).eq.',') goto 48
      elseif (c4.eq.'coll') then
         ipos=ipos+5
         call getchar(string,ipos,cospec,0)
         icoll(ipl)=igetcoind(cospec,conam,nco)
      elseif (c4.eq.'lcll') then
         ipos=ipos+5
         call getchar(string,ipos,cospec,0)
         ilcll(ipl)=igetcoind(cospec,conam,nco)
      elseif (c4.eq.'lchl') then
         ipos=ipos+5
         call getchar(string,ipos,cospec,0)
         ilchl(ipl)=igetcoind(cospec,conam,nco)
      elseif (c4.eq.'lclo') then
         ipos=ipos+5
         call getchar(string,ipos,cospec,0)
         ilclo(ipl)=igetcoind(cospec,conam,nco)
      elseif (c4.eq.'tslb') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,rtslb(ipl))
      elseif (c4.eq.'tshl') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,rtshl(ipl))
      elseif (c4.eq.'cong') then
         ipos=ipos+5
         call getchar(string,ipos,cospec,0)
         icong(ipl)=igetcoind(cospec,conam,nco)
      elseif (c4.eq.'conl') then
         ipos=ipos+5
         call getchar(string,ipos,cospec,0)
         iconl(ipl)=igetcoind(cospec,conam,nco)
      elseif (c4.eq.'cozr') then
         ipos=ipos+5
         call getchar(string,ipos,cospec,0)
         icozr(ipl)=igetcoind(cospec,conam,nco)
      elseif (c4.eq.'lwll') then
         ipos=ipos+5
         call getinum(string,ipos,ilwll(ipl))
      elseif (c4.eq.'lwng') then
         ipos=ipos+5
         call getinum(string,ipos,ilwng(ipl))
      elseif (c4.eq.'lwnl') then
         ipos=ipos+5
         call getinum(string,ipos,ilwnl(ipl))
      elseif (c4.eq.'lwzr') then
         ipos=ipos+5
         call getinum(string,ipos,ilwzr(ipl))
      elseif (c4.eq.'dall') then
         ipos=ipos+5
         call getinum(string,ipos,idall(ipl))
      elseif (c4.eq.'dang') then
         ipos=ipos+5
         call getinum(string,ipos,idang(ipl))
      elseif (c4.eq.'danl') then
         ipos=ipos+5
         call getinum(string,ipos,idanl(ipl))
      elseif (c4.eq.'dazr') then
         ipos=ipos+5
         call getinum(string,ipos,idazr(ipl))
      elseif (c4.eq.'lcnl') then
         ipos=ipos+5
         call getchar(string,ipos,cospec,0)
         ilcnl(ipl)=igetcoind(cospec,conam,nco)
      elseif (c4.eq.'lczr') then
         ipos=ipos+5
         call getchar(string,ipos,cospec,0)
         ilczr(ipl)=igetcoind(cospec,conam,nco)
      elseif (c4.eq.'lcbr') then
         ipos=ipos+5
         call getchar(string,ipos,cospec,0)
         ilcbr(ipl)=igetcoind(cospec,conam,nco)
      elseif (c4.eq.'orlb') then
         ipos=ipos+5
         call getinum(string,ipos,iorlb(ipl))
      elseif (c4.eq.'pwlb') then
         ipos=ipos+5
         call getinum(string,ipos,ipwlb(ipl))
      elseif (c4.eq.'pwhl') then
         ipos=ipos+5
         call getinum(string,ipos,ipwhl(ipl))
      elseif (c4.eq.'pwbr') then
         ipos=ipos+5
         call getinum(string,ipos,ipwbr(ipl))
      elseif (c4.eq.'fclb') then
         ipos=ipos+5
         call getchar(string,ipos,cospec,0)
         ifclb(ipl)=igetcoind(cospec,conam,nco)
      elseif (c4.eq.'fcnl') then
         ipos=ipos+5
         call getchar(string,ipos,cospec,0)
         ifcnl(ipl)=igetcoind(cospec,conam,nco)
      elseif (c4.eq.'fczr') then
         ipos=ipos+5
         call getchar(string,ipos,cospec,0)
         ifczr(ipl)=igetcoind(cospec,conam,nco)
      elseif (c4.eq.'fchl') then
         ipos=ipos+5
         call getchar(string,ipos,cospec,0)
         ifchl(ipl)=igetcoind(cospec,conam,nco)
      elseif (c4.eq.'fclo') then
         ipos=ipos+5
         call getchar(string,ipos,cospec,0)
         ifclo(ipl)=igetcoind(cospec,conam,nco)
      elseif (c4.eq.'mjsk') then
         ipos=ipos+5
         call getinum(string,ipos,imjsk(ipl))
      elseif (c4.eq.'cmth') then
         ipos=ipos+5
         call getchar(string,ipos,ccmth(ipl),0)
      elseif (c4.eq.'axlg') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,raxlg(ipl))
      elseif (c4.eq.'axld') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,raxld(ipl))
      elseif (c4.eq.'axlv') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,raxlv(ipl))
      elseif (c4.eq.'axtg') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,raxtg(ipl))
      elseif (c4.eq.'axtd') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,raxtd(ipl))
      elseif (c4.eq.'axtv') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,raxtv(ipl))
      elseif (c4.eq.'addf') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,raddf(ipl))
      elseif (c4.eq.'mply') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,rmply(ipl))
      elseif (c4.eq.'fdmn') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,rfdmn(ipl))
      elseif (c4.eq.'fdmx') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,rfdmx(ipl))
      elseif (c4.eq.'gdir') then
         ipos=ipos+5
         call getinum(string,ipos,igdir(ipl))
      elseif (c4.eq.'tjfl') then
         ipos=ipos+5
         call getchar(string,ipos,ctjfl(ipl),0)
      elseif (c4.eq.'titl') then
         ipos=ipos+5
         call getchar(string,ipos,ctitl(ipl),1)
      elseif (c4.eq.'v5nm') then
         ipos=ipos+5
         call getchar(string,ipos,cv5nm(ipl),0)
      elseif (c4.eq.'tjsp') then
         ipos=ipos+5
         itjns(ipl)=0
 49      itjns(ipl)=itjns(ipl)+1
         call getrnum(string,ipos,mkzh,rtjsp(1,itjns(ipl),ipl)) !time
         call getrnum(string,ipos,mkzh,rtjsp(3,itjns(ipl),ipl)) !x
         call getrnum(string,ipos,mkzh,rtjsp(2,itjns(ipl),ipl)) !y
         if (string(ipos-1:ipos-1).eq.',') goto 49
      elseif (c4.eq.'tjid') then
         ipos=ipos+5
         itjni(ipl)=0
 50      itjni(ipl)=itjni(ipl)+1
         call getinum(string,ipos,itjid(itjni(ipl),ipl))
         if (string(ipos-1:ipos-1).eq.',') goto 50
      elseif (c4.eq.'tjar') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,rtjar(2,ipl))
         if (string(ipos-1:ipos-1).eq.',') then
            rtjar(1,ipl)=rtjar(2,ipl)
            call getrnum(string,ipos,mkzh,rtjar(2,ipl))
         endif
      elseif (c4.eq.'tjst') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,rtjst(ipl))
      elseif (c4.eq.'tjen') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,rtjen(ipl))
      elseif (c4.eq.'tjti') then
         ipos=ipos+5
         call getrnum(string,ipos,mkzh,rtjti(ipl))
      else
         write(iup,*)'Don''t recognize the plot specifier.'
         write(iup,*)'   For ipl,ipos= ',ipl,ipos
         write(iup,*)'   at line ',linec
         write(iup,'(a10,a4,a1)') '    c4 = "',c4,'"'
         write(iup,*)'   in the following string:'
         write(iup,'(a)') string
         stop
      endif
c
c   The RIP convention is that the final plot specifier in
c   a plot specification statement not be followed by a
c   semicolon.  However, users often do it anyway, so as a
c   magnanimous gesture, rip will attempt to forgive this
c   error and recognize the end of the PSS, even though
c   there is a semicolon.
c
      if (string(ipos-1:ipos-1).eq.';'.and.
     &    string(ipos:ipos).ne.' ') goto 30
c
c   Warn user if no ptyp was specified, and stop.
c
      if (cptyp(ipl).eq.'xx') then
         if (itrajcalc.eq.0.and.imakev5d.eq.0) then
           write(iup,*)'Warning: you must always specifiy a plot type'
           write(iup,*)'(i.e. ptyp=something) [except when doing a'
           write(iup,*)'trajectory calculation run (itrajcalc=1)'
           write(iup,*)'or a Vis5D data creation run (imakev5d=1)].'
           write(iup,*)'You didn''t for plot spec. statment number ',ipl
           write(iup,*)'RIP will assume you want a hor. contour plot.'
         endif
         cptyp(ipl)='hc'
      endif
c
c   Set contour parameters to default values if they weren't
c   set explicitly.  Also set these for trajectories, since
c   a few of them get used for trajs.
c
      if (cptyp(ipl)(2:2).eq.'c'.or.cptyp(ipl)(2:2).eq.'t') then
c
         if (icozr(ipl).eq.999999) icozr(ipl)=icolr(ipl)
         if (ilczr(ipl).eq.999999) ilczr(ipl)=icozr(ipl)
         if (icoll(ipl).ne.999999.and.icong(ipl).eq.999999) then
            icong(ipl)=icolr(ipl)
            if (iconl(ipl).eq.999999) iconl(ipl)=icoll(ipl)
         elseif (icoll(ipl).eq.999999.and.icong(ipl).ne.999999) then
            icoll(ipl)=icolr(ipl)
            if (iconl(ipl).eq.999999) iconl(ipl)=icong(ipl)
         elseif (icoll(ipl).eq.999999.and.icong(ipl).eq.999999) then
            icoll(ipl)=icolr(ipl)
            icong(ipl)=icolr(ipl)
            if (iconl(ipl).eq.999999) iconl(ipl)=icolr(ipl)
         endif
         if (ilcll(ipl).eq.999999) ilcll(ipl)=icoll(ipl)
         if (ilcnl(ipl).eq.999999) ilcnl(ipl)=iconl(ipl)
         if (ilchl(ipl).eq.999999) ilchl(ipl)=icolr(ipl)
         if (ilclo(ipl).eq.999999) ilclo(ipl)=ilchl(ipl)
c
         if (ilwzr(ipl).eq.999999) ilwzr(ipl)=ilinw(ipl)
         if (ilwll(ipl).ne.999999.and.ilwng(ipl).eq.999999) then
            ilwng(ipl)=ilinw(ipl)
            if (ilwnl(ipl).eq.999999) ilwnl(ipl)=ilwll(ipl)
         elseif (ilwll(ipl).eq.999999.and.ilwng(ipl).ne.999999) then
            ilwll(ipl)=ilinw(ipl)
            if (ilwnl(ipl).eq.999999) ilwnl(ipl)=ilwng(ipl)
         elseif (ilwll(ipl).eq.999999.and.ilwng(ipl).eq.999999) then
            ilwll(ipl)=ilinw(ipl)
            ilwng(ipl)=ilinw(ipl)
            if (ilwnl(ipl).eq.999999) ilwnl(ipl)=ilinw(ipl)
         endif
c
         if (idazr(ipl).eq.999999) idazr(ipl)=idash(ipl)
         if (idall(ipl).ne.999999.and.idang(ipl).eq.999999) then
            idang(ipl)=idash(ipl)
            if (idanl(ipl).eq.999999) idanl(ipl)=idall(ipl)
         elseif (idall(ipl).eq.999999.and.idang(ipl).ne.999999) then
            idall(ipl)=idash(ipl)
            if (idanl(ipl).eq.999999) idanl(ipl)=idang(ipl)
         elseif (idall(ipl).eq.999999.and.idang(ipl).eq.999999) then
            idall(ipl)=idash(ipl)
            idang(ipl)=idash(ipl)
            if (idanl(ipl).eq.999999) idanl(ipl)=idash(ipl)
         endif
c
         if (ipwhl(ipl).eq.999999) ipwhl(ipl)=ipwlb(ipl)
c
         if (ifcnl(ipl).eq.999999) ifcnl(ipl)=ifclb(ipl)
         if (ifczr(ipl).eq.999999) ifczr(ipl)=ifclb(ipl)
         if (ifchl(ipl).eq.999999) ifchl(ipl)=ifclb(ipl)
         if (ifclo(ipl).eq.999999) ifclo(ipl)=ifchl(ipl)
c
         if (rtshl(ipl).eq.rmsg  ) rtshl(ipl)=rtslb(ipl)
c
      endif
c
c   Set map parameters to default values if they weren't
c   set explicitly
c
      if (cfeld(1,ipl)(1:3).eq.'map') then
         if (iouco(ipl).eq.999999) iouco(ipl)=icolr(ipl)
         if (couds(ipl).eq.'dot  ') then
            if (ioulw(ipl).eq.999999) ioulw(ipl)=12
         else
            if (ioulw(ipl).eq.999999) ioulw(ipl)=1
         endif
         if (imfco(1,ipl).ne.999999) then
            do i=2,6
               if (imfco(i,ipl).eq.999999) imfco(i,ipl)=imfco(i-1,ipl)
            enddo
         endif
         if (rcint(ipl).eq.rmsg) rcint(ipl)=10.
         if (idash(ipl).eq.70) idash(ipl)=32
      endif
c
c   Set a message color different from icolr if the value of
c   icolr ends up being close to the background color.
c   First, check to make sure that the first element of the array
c   that holds the color indices has a value of zero, which is the
c   background color index.
c
      diff=abs(fred(icolr(ipl))-fred(0))+
     &     abs(fgreen(icolr(ipl))-fgreen(0))+
     &     abs(fblue(icolr(ipl))-fblue(0))
      if (diff.lt..3) then
         bgtot=fred(0)+fgreen(0)+fblue(0)
         if (bgtot.lt.1.5) then
            icomg(ipl)=icoindwhite
         else
            icomg(ipl)=icoindblack
         endif
      else
         icomg(ipl)=icolr(ipl)
      endif
c
c   Set reference state parameters for perturbation fields
c   to default values if they weren't set explicitly
c
      if (rrfst(1,ipl).eq.0.) then
         rrfst(1,ipl)=refslp
         rrfst(2,ipl)=refslt
         rrfst(3,ipl)=reflaps
         rrfst(4,ipl)=refstratt
      endif
c
c   Check ixwin and iywin, and set to extreme values if the
c   values given are outside the extremes.
c
      if (ixwin(1,ipl).lt.1.or.ixwin(1,ipl).gt.mjx) ixwin(1,ipl)=1
      if (ixwin(2,ipl).lt.1.or.ixwin(2,ipl).gt.mjx) ixwin(2,ipl)=mjx
      if (iywin(1,ipl).lt.1.or.iywin(1,ipl).gt.miy) iywin(1,ipl)=1
      if (iywin(2,ipl).lt.1.or.iywin(2,ipl).gt.miy) iywin(2,ipl)=miy
c
c   Force rvwin to reset to missing values if cvcor has changed from
c   ipl-1 value.
c
      if (cvcor(ipl).ne.cvcorprev.and.inewrvwin.eq.0) then
         rvwin(1,ipl)=rmsg
         rvwin(2,ipl)=rmsg
      endif
      cvcorprev=cvcor(ipl)
c
c   Interpret the level info
c
      if (cptyp(ipl)(1:1).eq.'h') then
c
c      Deal with level sequences, averaging, and differencing
c
         ii=0
         ilev=0
  100    ii=ii+1
         if (ii.gt.inlvs(ipl)) goto 120
         if (rlevs(ii,ipl).ge.0.) then
            ilev=ilev+1
            rlevl(ilev,ipl)=rlevs(ii,ipl)
            rlavl(ilev,ipl)=rlevl(ilev,ipl)
         else
            ii=ii+1
            if (rlevs(ii,ipl).gt.0.) then
               rstart=rlevs(ii-2,ipl)
               rend=-rlevs(ii-1,ipl)
               rinc=rlevs(ii,ipl)
               rdist=rend-rstart
               isign=nint(rdist/abs(rdist))
               nlseries=int(abs(rdist)/rinc+.00001) + 1
               do 110 i=2,nlseries
                  ilev=ilev+1
                  rlevl(ilev,ipl)=rlevl(ilev-1,ipl)+isign*rinc
                  rlavl(ilev,ipl)=rlevl(ilev,ipl)
  110          continue
            elseif (rlevs(ii,ipl).eq.0.) then
               rlavl(ilev,ipl)=-rlevs(ii-1,ipl)
            elseif (rlevs(ii,ipl).eq.-1.) then
               rlavl(ilev,ipl)=rlevs(ii-1,ipl)
            else
               write(iup,*)'Error in level-series or',
     &            ' level-averaging spec.'
               stop
            endif
         endif
         goto 100
  120    inlvl(ipl)=ilev
c
c      Convert sigma value specifications to nearest k indices.
c
         if (cvcor(ipl).eq.'s') then
            do 130 i=1,inlvl(ipl)
               if (rlevl(i,ipl).lt.1.0) then
                  if (rlevl(i,ipl).ge..5*(sigh(mkzh)+
     &                sigh(mkzh-1))) then
                     rlevl(i,ipl)=mkzh
                  elseif (rlevl(i,ipl).le..5*(sigh(1)+sigh(2))) then
                     rlevl(i,ipl)=1.0
                  else
                     do 122 k=2,mkzh-1
                        if (rlevl(i,ipl).ge..5*(sigh(k-1)+sigh(k)).and.
     &                     rlevl(i,ipl).le..5*(sigh(k+1)+sigh(k))) then
                           rlevl(i,ipl)=k
                           goto 124
                        endif
  122                continue
  124                continue
                  endif
               else
                  rlevl(i,ipl)=min(rlevl(i,ipl),float(mkzh))
               endif
               isign=nint(rlavl(i,ipl)/abs(rlavl(i,ipl)))
               rlavl(i,ipl)=isign*rlavl(i,ipl)
               if (rlavl(i,ipl).lt.1.0) then
                  if (rlavl(i,ipl).ge..5*(sigh(mkzh)+
     &                sigh(mkzh-1))) then
                     rlavl(i,ipl)=mkzh
                  elseif (rlavl(i,ipl).le..5*(sigh(1)+sigh(2))) then
                     rlavl(i,ipl)=1.0
                  else
                     do 126 k=2,mkzh-1
                        if (rlavl(i,ipl).ge..5*(sigh(k-1)+sigh(k)).and.
     &                     rlavl(i,ipl).le..5*(sigh(k+1)+sigh(k))) then
                           rlavl(i,ipl)=k
                           goto 128
                        endif
  126                continue
  128                continue
                  endif
               else
                  rlavl(i,ipl)=min(rlavl(i,ipl),float(mkzh))
               endif
               rlavl(i,ipl)=isign*rlavl(i,ipl)
  130       continue
         endif
      else
         rlevl(1,ipl)=1.
         rlavl(1,ipl)=1.
         inlvl(ipl)=1
      endif
c
      goto 15
c
  500 nfr=ifr
c
      if (nfr.eq.0) then
c
c      Have to make at least ONE plot...
c
         iuinputsv=iuinput
         iuinput=83
         open(unit=iuinput,file='tmp.in',status='unknown',
     &      form='formatted')
         write(iuinput,'(a)') '=========================='//
     &      '================================================='
         write(iuinput,'(a)') '----------------------    '//
     &      'Plot Specification Table    ---------------------'
         write(iuinput,'(a)') '=========================='//
     &      '================================================='
         write(iuinput,'(a)')'feld=ght'
         write(iuinput,'(a)')'==========='
         close (iuinput)
         open(unit=iuinput,file='tmp.in',status='old',
     &      form='formatted')
         read(iuinput,*)
         read(iuinput,*)
         read(iuinput,*)
	 linec = linec + 3
         goto 15
      endif
c
c   Adjust level info for frames.
c
      iplef=0
      do 600 ifr=1,nfr
         iplbf=iplef+1
         iplef=iplbf+numppf(ifr)-1
         mlevf=1
         do 550 i=iplbf,iplef
            mlevf=max(mlevf,inlvl(i))
  550    continue
         do 567 i=iplbf,iplef
            do 565 j=inlvl(i)+1,mlevf
               rlevl(j,i) = rlevl(j-1,i)
               rlavl(j,i) = rlavl(j-1,i)
  565       continue
  567    continue
         nltf(ifr)=mlevf
  600 continue
c
      if (iuinput.eq.83) iuinput=iuinputsv
c
      return
      end