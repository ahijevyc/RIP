c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hcondraw(xtime,ilinw,sigf,vc3d,tmk,qvp,
     &         prs,ght,ter,pstx,prs_tsf,sigh,iprog,
     &         ixwin,iywin,ismth,rcint,rcbeg,rcend,lmult,larng,
     &         idash,rlevl,rlavl,cnohl,lnolb,lnobr,lnozr,incon,
     &         bottextfloor,cfeld,cvcor,work,icdwk,unwk,ilev,
     &         icolr,icoll,ilcll,ilchl,rtslb,rtshl,imjsk,icomg,
     &         lnmsg,icong,iconl,icozr,idimn,lhide,lgrad,lhadv,
     &         ilwll,ilwng,ilwnl,ilwzr,idall,
     &         idang,idanl,idazr,ilcnl,ilczr,
     &         ilcbr,ipwlb,iorlb,ipwhl,ipwbr,ifclb,ifcnl,
     &         ifczr,ifchl,ilclo,ifclo,ccmth,rwdbr,ihvbr,
     &         idotser,tseryi,tserxj,tserdat,ntsers,ntsert,ntserv,
     &         icosq,rcosq,incsq,fred,fgreen,fblue,nco,icomax,pslab1,
     &         iam,xcs,ycs,niam,ncs,idiffflag,
     &         maxtserv,maxtsers,maxtsert,maxcosq,mabpl,morpl,
     &         maxlev,maxpl,miy,mjx,mkzh,ipl,rrota)
c
      parameter (maxcon=80)
c
      dimension sigf(mkzh+1),vc3d(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   qvp(miy,mjx,mkzh),ght(miy,mjx,mkzh),prs(miy,mjx,mkzh),
     &   ter(miy,mjx),pstx(miy,mjx),sigh(mkzh),ixwin(2,maxpl),
     &   iywin(2,maxpl),ismth(maxpl),prs_tsf(miy,mjx),
     &   rcint(maxpl),rcbeg(maxpl),rcend(maxpl),
     &   idash(maxpl),rlevl(maxlev,maxpl),incsq(maxpl),
     &   rlavl(maxlev,maxpl),ilinw(maxpl),incon(maxpl),icomg(maxpl),
     &   icolr(maxpl),icoll(maxpl),ilcll(maxpl),ilchl(maxpl),
     &   rtslb(maxpl),rtshl(maxpl),imjsk(maxpl),idimn(maxpl),
     &   icosq(maxcosq,maxpl),rcosq(maxcosq,maxpl),
     &   icong(maxpl),iconl(maxpl),icozr(maxpl),ilwll(maxpl),
     &   ilwng(maxpl),ilwnl(maxpl),ilwzr(maxpl),idall(maxpl),
     &   idang(maxpl),idanl(maxpl),idazr(maxpl),
     &   ilcnl(maxpl),ilczr(maxpl),
     &   ilcbr(maxpl),ipwlb(maxpl),iorlb(maxpl),
     &   ipwhl(maxpl),ipwbr(maxpl),ifclb(maxpl),ifcnl(maxpl),
     &   ifczr(maxpl),ifchl(maxpl),
     &   ilclo(maxpl),ifclo(maxpl),rwdbr(maxpl),ihvbr(maxpl),
     &   work(miy,mjx,mkzh),icdwk(maxpl),
     &   tserdat(maxtsert,maxtserv,maxtsers),tseryi(maxtsers),
     &   tserxj(maxtsers),fred(0:255),fgreen(0:255),fblue(0:255),
     &   pslab1(mabpl,morpl),rrota(maxpl)
      dimension iam(niam),xcs(ncs),ycs(ncs),temp(mabpl,morpl),
     &   temp2(morpl,mabpl)
      logical lnolb(maxpl),lnobr(maxpl),lnozr(maxpl),
     &   lmult(maxpl),larng(maxpl),lnmsg(maxpl),lhide(maxpl),
     &   lgrad(maxpl),lhadv(maxpl),lhidel
      character cfeld(3,maxpl)*10,cvcor(maxpl)*1,ccmth(maxpl)*4,
     &   unwk(maxpl)*24,cnohl(maxpl)*1
c
      dimension redcosq(500),greencosq(500),bluecosq(500)
      dimension rwrk(5000),iwrk(5000),rect(4),
     &   valcon(maxcon),majcon(maxcon),lbbarcon(maxcon),
     &   iaia(500),igia(500),lfin(maxcon+1)
      character messg*126,llbs(maxcon+2)*24
      real ch,cv,uh,uv
      REAL psmin, psmax
      external drawcl,cpcolr
c
      parameter (nsp=6)
      dimension spx(nsp),spy(nsp),spt(nsp)
c      data (spt(i),i=1,nsp) / 24.,25.,27.,30.,33.,36. /
c      data (spx(i),i=1,nsp) / 85.,87.,90.,95.,98.,99. /
c      data (spy(i),i=1,nsp) / 56.,57.,58.,69.,70.,73. /
      data (spt(i),i=1,nsp) / 0.,1.,3.,6.,9.,12. /
      data (spx(i),i=1,nsp) / 37.,37.,39.,41.,43.,46. /
      data (spy(i),i=1,nsp) / 42.,42.,44.,47.,51.,56. /
c
      include 'comconst'
c
      dimension mconcp(1000),icoindcp(1000)
      common /cpack/ mconcp,icoindcp,nconarea,icpfchl,icpfclo,
     &   icpfclb,icpfcnl,icpfczr
c
c   Make proper set call, as well as number of values to
c      plot in each direction
c
      xintervs=ixwin(2,ipl)-ixwin(1,ipl)
      yintervs=iywin(2,ipl)-iywin(1,ipl)
      faspect=(ftmax-fbmin)/(frmax-flmin)
      aspect=yintervs/xintervs
      if (aspect.lt.faspect) then
         fl=flmin
         fr=frmax
         fextra=.5*((ftmax-fbmin)-aspect*(frmax-flmin))
         fb=fbmin+fextra
         ft=ftmax-fextra
      else
         fb=fbmin
         ft=ftmax
         fextra=.5*((frmax-flmin)-1./aspect*(ftmax-fbmin))
         fl=flmin+fextra
         fr=frmax-fextra
      endif
      if (icdwk(ipl).eq.0) then
         ul=1.
         ur=1.+xintervs
         ub=1.
         ut=1.+yintervs
         niy=nint(ut)
         njx=nint(ur)
      else
         ul=.5
         ur=.5+xintervs
         ub=.5
         ut=.5+yintervs
         niy=int(ut)
         njx=int(ur)
      endif
      if(rrota(ipl) .ne. 90. .and. rrota(ipl) .ne. -90.) then
        call set(fl,fr,fb,ft,ul,ur,ub,ut,1)
c MGD begin mod
c need to reverse dimensions of set() window if rotating by (-)90
      else
        cv = (ft-fb)/2.
        ch = (fr-fl)/2.
        call set((fr+fl)/2.-cv,(fr+fl)/2.+cv,
     &           (ft+fb)/2.-ch,(ft+fb)/2.+ch,ub,ut,ul,ur,1)
      endif
c MGD end mod
c
c   Put appropriate data into horizontal slab.
c
      if ((cvcor(ipl).eq.'s'.and.rlevl(ilev,ipl).eq.
     &   rlavl(ilev,ipl)).or.idimn(ipl).eq.2) then
         rmaxv=-rmsg
         rminv=rmsg
         do 90 j=1,njx
            jj=j+ixwin(1,ipl)-1
            do 80 i=1,niy
               ii=i+iywin(1,ipl)-1
               pslab1(j,i)=work(ii,jj,nint(rlevl(ilev,ipl)))
               rmaxv=max(rmaxv,pslab1(j,i))
               rminv=min(rminv,pslab1(j,i))
   80       continue
   90    continue
      elseif (cvcor(ipl).eq.'s'.and.rlevl(ilev,ipl).ne.
     &      rlavl(ilev,ipl).and.rlavl(ilev,ipl).ge.0) then
         call fillarray(pslab1,mabpl*morpl,0.)
         lev1=min(nint(rlevl(ilev,ipl)),nint(rlavl(ilev,ipl)))
         lev2=max(nint(rlevl(ilev,ipl)),nint(rlavl(ilev,ipl)))
         sigtot=sigf(lev2+1)-sigf(lev1)
         do 125 k=lev1,lev2
         do 120 j=1,njx
            jj=j+ixwin(1,ipl)-1
            do 115 i=1,niy
               ii=i+iywin(1,ipl)-1
               pslab1(j,i)=pslab1(j,i)+work(ii,jj,k)*
     &          (sigf(k+1)-sigf(k))/sigtot
  115       continue
  120    continue
  125    continue
      elseif (cvcor(ipl).eq.'s'.and.rlevl(ilev,ipl).ne.
     &      rlavl(ilev,ipl).and.rlavl(ilev,ipl).lt.0) then
         lev1=nint(rlevl(ilev,ipl))
         lev2=nint(-rlavl(ilev,ipl))
         do j=1,njx
            jj=j+ixwin(1,ipl)-1
            do i=1,niy
               ii=i+iywin(1,ipl)-1
               pslab1(j,i)=work(ii,jj,lev1)-work(ii,jj,lev2)
            enddo
         enddo
      else
         if (lgrad(ipl).or.lhadv(ipl)) then
            lhidel=.true.
         else
            lhidel=lhide(ipl)
         endif
         call vinterp(cvcor(ipl),rlevl(ilev,ipl),ixwin(1,ipl),
     &      iywin(1,ipl),icdwk(ipl),vc3d,tmk,qvp,
     &      prs,ght,ter,pstx,sigh,sigf,prs_tsf,
     &      lhidel,idiffflag,cfeld(1,ipl),work,
     &      pslab1,iprog,mabpl,morpl,njx,niy,miy,mjx,mkzh)
      endif
      PSMAX  = -500.
      PSMIN  =  500.
      pssum  =    0.
      pssum2 =    0.
      DO 150 J=1,NJX
      DO 150 I=1,NIY
         PSMAX  = MAX(PSMAX,PSLAB1(J,I))
         PSMIN  = MIN(PSMIN,PSLAB1(J,I))
         pssum  = pssum  + PSLAB1(J,I)
         pssum2 = pssum2 + PSLAB1(J,I)**2
  150 CONTINUE
      WRITE(IUP,*)'njx,niy=',njx,niy
      WRITE(IUP,*)'MAX, MIN TEMP.=',PSMAX,PSMIN
      WRITE(iup,*)'sum,sum2,mean,rmse=',pssum,pssum2,
     & pssum/njx/niy, SQRT(pssum2/njx/niy)
c
c   Smooth field if called for.
c
      call smooth(pslab1,ismth(ipl),mabpl,njx,niy)
c
c   Calculate time series data if called for.
c
      do 300 istn=1,idotser*ntsers
         posy=tseryi(istn)-iywin(1,ipl)+1.-.5*icdwk(ipl)
         posx=tserxj(istn)-ixwin(1,ipl)+1.-.5*icdwk(ipl)
         if (posx.le.1..or.posx.ge.float(njx).or.
     &       posy.le.1..or.posy.ge.float(niy)) goto 290
         jl=int(posx)
         jr=jl+1
         ib=int(posy)
         it=ib+1
         ratlr=posx-jl
         ratbt=posy-ib
         if (pslab1(jl,it).eq.rmsg.or.
     &       pslab1(jr,it).eq.rmsg.or.
     &       pslab1(jl,ib).eq.rmsg.or.
     &       pslab1(jr,ib).eq.rmsg) goto 290
         tserdat(ntsert,ntserv,istn)=
     +      (1.-ratlr)*(   ratbt)*pslab1(jl,it)+
     +      (   ratlr)*(   ratbt)*pslab1(jr,it)+
     +      (1.-ratlr)*(1.-ratbt)*pslab1(jl,ib)+
     +      (   ratlr)*(1.-ratbt)*pslab1(jr,ib)
         goto 300
  290    tserdat(ntsert,ntserv,istn)=rmsg
  300 continue
c      if (idotser.eq.1) write(iup,*)'   rmaxv,rminv=',rmaxv,rminv
cc
cc   Do low center thing
cc
c      if (lnobr(ipl)) then
c         do isp=1,nsp
c            if (abs(spt(isp)-xtime).le..15) then
c               ispg=isp
c               goto 543
c            endif
c         enddo
c         stop 333
c 543     continue
c         dist=800./dskm
c         pouter=0.
c         do ideg=0,350,10
c            xring=spx(ispg)+dist*cos(rpd*ideg)
c            yring=spy(ispg)+dist*sin(rpd*ideg)
c            ixring=nint(min(max(xring,1.),mjx-1.))
c            iyring=nint(min(max(yring,1.),miy-1.))
c            pouter=pouter+pslab1(ixring,iyring)
c         enddo
c         pouter=pouter/36.
cc         write(iup,*) 'x,y=',nint(spx(ispg)),nint(spy(ispg))
c         pcentral=pslab1(nint(spx(ispg)),nint(spy(ispg)))
c         write(iup,*) 'time,pcentral,pouter,pdiff=',
c     &      spt(ispg),pcentral,pouter,pcentral-pouter
c      endif
      if (rwdbr(ipl).ge.100.) then  ! one-time code - ignore this stuff
         gradient=rcosq(1,ipl)  !  in hPa/km
         direction=rcosq(2,ipl)  ! degrees
         valinterior=rcosq(3,ipl)
         xinterior=41.
         yinterior=47.
c mod this too...
         do j=1,mjx
         do i=1,miy
            gval=valinterior+gradient*dskm*(
     &         (j-xinterior)*cos(rpd*direction)+
     &         (i-yinterior)*sin(rpd*direction))
            if (rwdbr(ipl).eq.100.) then
               pslab1(j,i)=gval
            else
               pslab1(j,i)=pslab1(j,i)+gval
            endif
         enddo
         enddo
      endif
c
c   Determine number of contours and contour values.
c
      call getconvals(rcbeg(ipl),rcend(ipl),rcint(ipl),
     &   incon(ipl),lmult(ipl),imjsk(ipl),pslab1,mabpl,njx,niy,rmsg,
     &   maxcon,valcon,majcon,cintuse,numcon,iup)
c
      call cpseti('SET',0)   ! we'll use our own set call
      call cpseti('CLS',0)   ! we'll use our own cont. values
c
c   Some settings for hi/lo and contour labels
c
      call cpseti('NSD',-4)
      call cpseti('NOF',7)
      call cpsetr ('PC6',.6)
c
c   Filling:
c
      if ((ccmth(ipl).eq.'fill'.or.ccmth(ipl).eq.'both').and.
     &    numcon.gt.0) then

c MGD begin mod
c for ccw rotation, x becomes y and inverted y becomes x
c we just copy in this manner to a new array with reversed dimensions
      if(rrota(ipl) .eq. 90.) then
        do i=1,njx
        do j=1,niy
          temp2(niy-j+1,i) = pslab1(i,j)
        enddo
        enddo
c for cw rotation, inverted x becomes y, and y becomes x
      elseif(rrota(ipl) .eq. -90.) then
        do i=1,njx
        do j=1,niy
          temp2(j,njx-i+1) = pslab1(i,j)
        enddo
        enddo
c for 180 rotation, invert both x and y
      elseif(rrota(ipl) .eq. -180. .or. rrota(ipl) .eq. 180.) then
        do i=1,njx
        do j=1,niy
          temp(njx-i+1,niy-j+1) = pslab1(i,j)
        enddo
        enddo
      endif

c 180 degree rotations yield same dimension of matrix, so go ahead 
c and copy data back to the old one
      if(rrota(ipl) .ne. 90. .and. rrota(ipl) .ne. -90. .and. rrota(ipl)
     &   .ne. 0.) then
      do i=1,njx
      do j=1,niy
        pslab1(i,j) = temp(i,j)
      enddo
      enddo
      endif
c MGD end mod

c
c   Assign the red, green, and blue fractions to the color sequence
c
      do i=1,incsq(ipl)
         do j=0,nco
            if (icosq(i,ipl).eq.j) then
               redcosq(i)=fred(j)
               greencosq(i)=fgreen(j)
               bluecosq(i)=fblue(j)
               goto 320
            endif
         enddo
         redcosq(i)=1.
         greencosq(i)=1.
         bluecosq(i)=1.
 320     continue
      enddo
c
c   Colors are assigned to ranges of values bounded by adjacent
c   contour values, as well as the semi-infinite ranges below
c   the lowest and above the highest contours.  Therefore, numcon+1
c   colors need to be defined. In order to choose a color from the
c   color sequence, each range must be assigned a specific value.
c   That value will be set equal to the average of the two
c   surrounding contour values.  The value for the range
c   <minus infinity to valcon(1)> will be assigned a value of the
c   lowest contour minus one half the difference between the lowest
c   and second lowest contours.  A similar definition applies to
c   the value for the range <valcon(numcon) to plus infinity>.
c   Based on the specific value assigned to each range, a color will
c   be defined based on a linear interpolation from the color sequence.
c
      nconarea=numcon+1
      mm=incsq(ipl)
      if (numcon.gt.1) then
         valmin=valcon(1)-.5*(valcon(2)-valcon(1))
         valmax=valcon(numcon)+.5*(valcon(numcon)-valcon(numcon-1))
      else
c
c      If there is only one contour value generated, then
c      set valmin to the average of all grid values that are less than
c      the single contour value, and valmax to the average of all
c      grid values that are greater than the single contour value, 
c
         numgt=0
         avggt=0.
         numlt=0
         avglt=0.
c MGD mod
         if(rrota(ipl).ne.90. .or. rrota(ipl).ne.-90.) then
         do j=1,njx
         do i=1,niy
            if (pslab1(j,i).ne.rmsg.and.
     &          pslab1(j,i).gt.valcon(1)) then
               numgt=numgt+1
               avggt=avggt+pslab1(j,i)
            elseif (pslab1(j,i).ne.rmsg.and.
     &              pslab1(j,i).lt.valcon(1)) then
               numlt=numlt+1
               avglt=avglt+pslab1(j,i)
            endif
         enddo
         enddo
c MGD begin mod
c avggt and avglt should be done for our array of reversed dimensions
c if we are working with either a 90 or -90 degree rotation. Same as
c normal case, but using transposed array and dimensions.
         else
         do j=1,niy
         do i=1,njx
            if (temp2(j,i).ne.rmsg.and.
     &          temp2(j,i).gt.valcon(1)) then
               numgt=numgt+1
               avggt=avggt+temp2(j,i)
            elseif (temp2(j,i).ne.rmsg.and.
     &              temp2(j,i).lt.valcon(1)) then
               numlt=numlt+1
               avglt=avglt+temp2(j,i)
            endif
         enddo
         enddo
         endif
         if (numgt.gt.0) then
            valmax=avggt/numgt
	 else
            valmax=valcon(1)
         endif
         if (numlt.gt.0) then
            valmin=avglt/numlt
	 else
            valmin=valcon(1)
         endif
      endif
c
      do i=1,numcon+1
         red=-1.
c
c      Note, in the following code, if icosc<0, that means the user
c      chose "transparent" as the color, which means the value of
c      "red" should remain -1.
c
         if (i.eq.1) then
            val=valmin
         elseif (i.eq.numcon+1) then
            val=valmax
         else
            val=.5*(valcon(i-1)+valcon(i))
         endif
         if (larng(ipl)) val=100.*(val-valmin)/(valmax-valmin)
         if (val.le.rcosq(1,ipl)) then
            if (icosq(1,ipl).gt.0) then
               red=redcosq(1)
               green=greencosq(1)
               blue=bluecosq(1)
            endif
         elseif (val.ge.rcosq(mm,ipl)) then
            if (icosq(mm,ipl).gt.0) then
               red=redcosq(mm)
               green=greencosq(mm)
               blue=bluecosq(mm)
            endif
         else
            do isq=1,incsq(ipl)-1
               if (val.ge.rcosq(isq,ipl).and.
     &             val.le.rcosq(isq+1,ipl)) then
                  if (icosq(isq,ipl).gt.0.and.
     &                icosq(isq+1,ipl).gt.0) then
                     fac=(val-rcosq(isq,ipl))/
     &                  (rcosq(isq+1,ipl)-rcosq(isq,ipl))
                     red=fac*redcosq(isq+1)+(1.-fac)*redcosq(isq)
                     green=fac*greencosq(isq+1)+(1.-fac)*greencosq(isq)
                     blue=fac*bluecosq(isq+1)+(1.-fac)*bluecosq(isq)
                  endif
                  goto 330
               endif
            enddo
 330        continue
         endif
c
c      Assign color index to icoindcp
c
         if (red.lt.0.) then
            icoindcp(i)=-1 ! transparent
         else
c
c         First check if color already exists
c
            do ico=2,icomax
               if (red.eq.fred(ico).and.green.eq.fgreen(ico).and.
     &             blue.eq.fblue(ico)) then
                  icoindcp(i)=ico
                  goto 357
               endif
            enddo
c
c         If not, set the new color and assign its index to icoincdcp
c
            icomax=icomax+1
            if (icomax.eq.256) then
               write(iup,*) 'hcondraw is trying to define color number'
               write(iup,*) '256, but 255 is the limit.  Try increasing'
               write(iup,*) 'the contour interval, or reduce the number'
               write(iup,*) 'of colors defined in the color table.'
               stop
            endif
            call gscr(1,icomax,red,green,blue)
            fred(icomax)=red
            fgreen(icomax)=green
            fblue(icomax)=blue
            icoindcp(i)=icomax
         endif
 357     continue
      enddo
c
c  Increase number of vertical strips for efficiency:
c
      call cpseti('NVS',4)
c
      call cpseti('NCL',numcon)
      do i=1,numcon
         call cpseti('PAI',i)
         call cpsetr('CLV',valcon(i))
         call cpseti('CLU',1)
         call cpseti('AIB',i)
         call cpseti('AIA',i+1)
      enddo
c
      if(rrota(ipl) .ne. 90. .and. rrota(ipl) .ne. -90.) then
      call cprect(pslab1,mabpl,njx,niy,rwrk,5000,iwrk,5000)
      call arinam(iam,niam)   ! initialize the area map
      call cpclam(pslab1,rwrk,iwrk,iam)  ! put cont. lines in area map
      call arscam(iam,xcs,ycs,ncs,iaia,igia,500,cpcolr)
c MGD begin mod
c if (-)90 degree rotation, then we need to use our array with reversed
c dimensions
      else
      call cprect(temp2,morpl,niy,njx,rwrk,5000,iwrk,5000)
      call arinam(iam,niam)   ! initialize the area map
      call cpclam(temp2,rwrk,iwrk,iam)  ! put cont. lines in area map
      call arscam(iam,xcs,ycs,ncs,iaia,igia,500,cpcolr)
      endif
c MGD end mod
c
c   Make label bar
c
      if (.not.lnobr(ipl)) then
c
         call lbseti('CBL',ilcbr(ipl))
         call lbseti('CLB',ilcbr(ipl))
         call lbsetr('WBL',float(ipwbr(ipl)))
         call lbsetr('WLB',float(ipwbr(ipl)))
c
c      Decide whether to orient the bar horizontally or vertically
c
c      Determine anticipated available space at bottom and on right.
c      This needn't be exact, it's just for roughly determining whether a
c      horizontal or vertical label bar will more efficiently use
c      the available space.  For bottom, assume .02 is needed for
c      tick labels (typically present) and perhaps 3 additional fields
c      will be plotted.  For right, assume .045 is needed for lat/lon
c      labels (typically present).
c
         avlbottom=fb-bottextfloor-.02-3.*.0152
         avlright=1.-fr-.045
c
         if (rwdbr(ipl).eq.-1.) then
            reqbottom=.036
            reqright=.07
         else
            reqbottom=rwdbr(ipl)
            reqright=rwdbr(ipl)
         endif
c
         ihov=ihvbr(ipl)
         if (ihov.eq.-1) then
            if (avlbottom-reqbottom.gt.avlright-reqright) then
               ihov=0
            else
               ihov=1
            endif
         endif
         if (ihov.eq.0) then
            xleb=fl
            xreb=fr
            ybeb=bottextfloor
            yteb=ybeb+reqbottom
            wsfb=1.
            hsfb=.45
         else
            xleb=1.-reqright
            xreb=1.
            ybeb=fb
            yteb=ft
            wsfb=.231
            hsfb=1.
         endif
         iftp=0
         lbab=1
c
c      Get label bar labels and colors
c
         lstep=numcon/20+1  ! don't overcrowd label. skip some values
c
c      Use getconvals to decide which contours (color transitions) should
c      be labeled on the label bar.
c
         call getconvals(rcbeg(ipl),rcend(ipl),rcint(ipl),
     &      incon(ipl),lmult(ipl),lstep-1,pslab1,mabpl,njx,niy,rmsg,
     &      maxcon,valcon,lbbarcon,cintuse,numcon,iup)
c
         nlbs=1
         nbox=1
         if (icoindcp(1).ge.0) then
            lfin(1)=icoindcp(1)
         else
            lfin(1)=0
         endif
         llbs(1)=' '
         do i = 1,numcon
            nlbs=nlbs+1
            nbox=nbox+1
            if ((lbbarcon(i).eq.-2.or.lbbarcon(i).eq.0.or.
     &           lbbarcon(i).eq.2).and.nlbs.le.numcon+2-lstep) then
               call cpseti ('PAI',i)          ! set internal array index
               call cpsetr ('ZDV',valcon(i))  ! give value to a converter
               call cpgetc ('ZDV',llbs(nlbs)) ! get a string back
            else
               llbs(nlbs)=' '
            endif
c            icii=min(i+1+lstep/2,nconarea)
c            if (icoindcp(icii).ge.0) then
c               lfin(nbox)=icoindcp(icii)
            if (icoindcp(nbox).ge.0) then
               lfin(nbox)=icoindcp(nbox)
            else
               lfin(nbox)=0
            endif
         enddo
         nlbs=nlbs+1
         llbs(nlbs)=unwk(ipl)
c
         call lblbar(ihov,xleb,xreb,ybeb,yteb,nbox,wsfb,hsfb,lfin,
     &      iftp,llbs,nlbs,lbab)
         if (ihov.eq.0) bottextfloor=bottextfloor+1.2*reqbottom
c
      endif
c
      endif
c
c   Contours:
c
      if ((ccmth(ipl).eq.'cont'.or.ccmth(ipl).eq.'both').and.
     &    numcon.gt.0) then
c
      call cpsetc('ILT',' ') ! no CONPACK informational label
      call cpseti('LBC',-1)  ! use current fill color for lab. boxes
c
c   Set up high/low stuff
c
      if (cnohl(ipl).ne.'B'.and.cnohl(ipl).ne.'b') then
c         call cpsetc('HIT','H~PRL~0~PRU~ $ZDV$')
c         call cpsetc('LOT','L~PRL~0~PRU~ $ZDV$')
         if (cnohl(ipl).eq.'H'.or.cnohl(ipl).eq.'h') then
            call cpsetc('HIT',' ')
         else
            call cpsetc('HIT','H~V-1QH-50~ $ZDV$')
         endif
         if (cnohl(ipl).eq.'L'.or.cnohl(ipl).eq.'l') then
            call cpsetc('LOT',' ')
         else
            call cpsetc('LOT','L~V-1QH-50~ $ZDV$')
         endif
         call cpseti('HIC',ilchl(ipl))
         call cpseti('LOC',ilclo(ipl))
         rtshlvp=rtshl(ipl)/(fr-fl)
c         rtshlvp=rtshl(ipl)
         call cpsetr('HLS',rtshlvp)
         if (ipwhl(ipl).ne.0) then
            call cpsetr('HLL',float(ipwhl(ipl)))
         else
            call cpsetr('HLL',1.)
         endif
         if (ifchl(ipl).ne.999999) then
            icpfchl=ifchl(ipl)
            if (ifclo(ipl).ne.999999) then
               icpfclo=ifclo(ipl)
            else
               icpfclo=ifchl(ipl)
            endif
         endif
         if (ipwhl(ipl).eq.0.and.ifchl(ipl).eq.999999) then
            ihlb=0
         elseif (ipwhl(ipl).ne.0.and.ifchl(ipl).eq.999999) then
            ihlb=1
         elseif (ipwhl(ipl).eq.0.and.ifchl(ipl).ne.999999) then
            ihlb=2
         else
            ihlb=3
         endif
         call cpseti('HLB',ihlb)
      else
         call cpsetc('HLT',' ')
      endif
c
c   Set up contour line stuff
c
c MGD begin mod
c this is all pretty much the same as for filling method, but it 
c works out better if we do the rotations separately for each
c if you want this explained again, see the filling section comments
      if(rrota(ipl) .eq. 90.) then
        do i=1,njx
        do j=1,niy
          temp2(niy-j+1,i) = pslab1(i,j)
c         temp2(niy-j,i) = pslab1(i,j)
        enddo
        enddo
      elseif(rrota(ipl) .eq. -90.) then
        do i=1,njx
        do j=1,niy
          temp2(j,njx-i+1) = pslab1(i,j)
c         temp2(j,njx-i) = pslab1(i,j)
        enddo
        enddo
      elseif(rrota(ipl) .eq. -180. .or. rrota(ipl) .eq. 180.) then
        do i=1,njx
        do j=1,niy
          temp(njx-i+1,niy-j+1) = pslab1(i,j)
c         temp(njx-i,niy-j) = pslab1(i,j)
        enddo
        enddo
      endif

      if(rrota(ipl) .ne. 90. .and. rrota(ipl) .ne. -90. .and. rrota(ipl)
     &   .ne. 0.) then
      do i=1,njx
      do j=1,niy
        pslab1(i,j) = temp(i,j)
      enddo
      enddo
      endif

      if(rrota(ipl) .ne. 90. .and. rrota(ipl) .ne. -90.) then
        call cprect(pslab1,mabpl,njx,niy,rwrk,5000,iwrk,5000)
      else
        call cprect(temp2,morpl,niy,njx,rwrk,5000,iwrk,5000)
      endif
c MGD end mod
      if (iorlb(ipl).eq.1) then
         call cpseti('LLP',3)   ! penalty scheme for line labels
         call cpseti('LLO',0)   ! all labels at angle=0 (hor. orient.)
      elseif (iorlb(ipl).eq.2) then
         call cpseti('LLP',3)   ! penalty scheme for line labels
         call cpseti('LLO',1)   ! all labels oriented along contour
      elseif (iorlb(ipl).eq.3) then
         call cpseti('LLP',1)   ! labels drawn with DASHPAT,
c                                 like in old CONREC
      else
         write(iup,*)'For ipl=',ipl,', invalid value for orlb.'
         write(iup,*)'Should be 1, 2, or 3.'
         stop
      endif
      rtslbvp=rtslb(ipl)/(fr-fl)
c      rtslbvp=rtslb(ipl)
      call cpsetr('LLS',rtslbvp)
      call cpseti('NCL',numcon)
      if (.not.lnolb(ipl)) then
         if (ipwlb(ipl).ne.0) then
            call cpsetr('LLL',float(ipwlb(ipl)))
         else
            call cpsetr('LLL',1.)
         endif
         if (ifclb(ipl).ne.999999) then
            icpfclb=ifclb(ipl)
            if (ifcnl(ipl).ne.999999) then
               icpfcnl=ifcnl(ipl)
            else
               icpfcnl=ifclb(ipl)
            endif
            if (ifczr(ipl).ne.999999) then
               icpfczr=ifczr(ipl)
            else
               icpfczr=ifclb(ipl)
            endif
         endif
         if (ipwlb(ipl).eq.0.and.ifclb(ipl).eq.999999) then
            illb=0
         elseif (ipwlb(ipl).ne.0.and.ifclb(ipl).eq.999999) then
            illb=1
         elseif (ipwlb(ipl).eq.0.and.ifclb(ipl).ne.999999) then
            illb=2
         else
            illb=3
         endif
         call cpseti('LLB',illb)
      endif
      do i=1,numcon
         call cpseti('PAI',i)
         call cpsetr('CLV',valcon(i))
         mconcp(i)=majcon(i)
         if (majcon(i).eq.1) then
            call cpseti('CLU',1)
            call cpseti('CLC',icolr(ipl))
            call cpsetr('CLL',float(ilinw(ipl)))
            call getdash(idash(ipl),ndot)
            call cpseti('CLD',ndot)
         elseif (majcon(i).eq.2) then
            if (.not.lnolb(ipl)) then
               call cpseti('CLU',3)
               call cpseti('LLC',ilcll(ipl))
            else
               call cpseti('CLU',1)
            endif
            call cpseti('CLC',icoll(ipl))
            call cpsetr('CLL',float(ilwll(ipl)))
            call getdash(idall(ipl),ndot)
            call cpseti('CLD',ndot)
         elseif (majcon(i).eq.-1) then
            call cpseti('CLU',1)
            call cpseti('CLC',icong(ipl))
            call cpsetr('CLL',float(ilwng(ipl)))
            call getdash(idang(ipl),ndot)
            call cpseti('CLD',ndot)
         elseif (majcon(i).eq.-2) then
            if (.not.lnolb(ipl)) then
               call cpseti('CLU',3)
               call cpseti('LLC',ilcnl(ipl))
            else
               call cpseti('CLU',1)
            endif
            call cpseti('CLC',iconl(ipl))
            call cpsetr('CLL',float(ilwnl(ipl)))
            call getdash(idanl(ipl),ndot)
            call cpseti('CLD',ndot)
         elseif (majcon(i).eq.0) then
            if (.not.lnozr(ipl)) then
               if (.not.lnolb(ipl)) then
                  call cpseti('CLU',3)
                  call cpseti('LLC',ilczr(ipl))
               else
                  call cpseti('CLU',1)
               endif
               call cpseti('CLC',icozr(ipl))
               call cpsetr('CLL',float(ilwzr(ipl)))
               call getdash(idazr(ipl),ndot)
               call cpseti('CLD',ndot)
            else
               call cpseti('CLU',0)
            endif
         endif
         call cpgeti('CLU',iclu)
      enddo

c
c MGD begin mod
c same deal as with filling - just use our own array with reversed order
c of dimensions for (-)90 degree rotated grids
      if(rrota(ipl) .ne. 90. .and. rrota(ipl) .ne. -90.) then
        call arinam(iam,niam)   ! initialize the area map
        call cplbam(pslab1,rwrk,iwrk,iam)  ! put label boxes in area map
        call cpcldm(pslab1,rwrk,iwrk,iam,drawcl)   ! draw contour lines
        call cplbdr(pslab1,rwrk,iwrk)   ! make labels
      else
        call arinam(iam,niam)   ! initialize the area map
        call cplbam(temp2,rwrk,iwrk,iam)  ! put label boxes in area map
        call cpcldm(temp2,rwrk,iwrk,iam,drawcl)   ! draw contour lines
        call cplbdr(temp2,rwrk,iwrk)   ! make labels
      endif
c MGD end mod
c
c   Write message at bottom
c
      if (.not.lnmsg(ipl)) then
c
      messg ='                                           LOW=123'//
     &       '456789012  HIGH=123456789012  INTERVAL= 1234567890'//
     &       '12  SCALE=123456789012'
c
c     messg ='CONTOURS:  UNITS=123456789012345678901234  LOW=123'//
c    &       '456789012  HIGH=123456789012  INTERVAL= 1234567890'//
c    &       '12  SCALE=123456789012'
C             12345678901234567890123456789012345678901234567890
C                      1         2         3         4         5
      do iii=24,1,-1
         if (unwk(ipl)(iii:iii).ne.' ') then
            ilch=iii
            goto 412
         endif
      enddo
      ilch=1
 412  continue
      ibwk=1+(24-ilch)
      messg(ibwk:ibwk+16)='CONTOURS:  UNITS='
      write(messg(ibwk+17:ibwk+16+ilch),'(a)') unwk(ipl)(1:ilch)
      write(messg(48:59),'(g12.5)') valcon(1)
      write(messg(67:78),'(g12.5)') valcon(numcon)
      if (.not.lmult(ipl)) then
         write(messg(91:102),'(g12.5)') cintuse
      else
         write(messg(90:102),'(a1,g12.5)') 'X',cintuse
      endif
      call cpgetr('SFU',ash)
      write(messg(111:122),'(g12.5)') ash
      nchar=122
      if (ash .eq. 1.) nchar = 102
      call gqclip (ierr,iclp,rect)
      call gsclip (0)
      call gstxci(icomg(ipl))
      call gsplci(icomg(ipl))
      chsize=.008
      ypos=bottextfloor+.5*chsize
      call pcgeti ('QU',ntextqq)
      call pcseti ('QU',0)
      call plchhq (cfux(.5),cfuy(ypos),messg(ibwk:nchar),chsize,0.,0.)
      call pcseti ('QU',ntextqq)
      call gsclip (iclp)
      bottextfloor=bottextfloor+1.9*chsize
c
      endif
c
      endif
c
      return
      end
