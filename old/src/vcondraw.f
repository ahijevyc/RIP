c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine vcondraw(ilinw,rcrag,rcrbg,ismth,rcint,rcbeg,rcend,
     &         lmult,larng,idash,ixavg,cnohl,lnolb,lnobr,lnozr,icomg,
     &         incon,bottextfloor,work,icdwk,unwk,icolr,icoll,ilcll,
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
c
      parameter (maxcon=80)
c
      dimension rcrag(2,maxpl),rcrbg(2,maxpl),ismth(maxpl),
     &   ilinw(maxpl),rcint(maxpl),rcbeg(maxpl),rcend(maxpl),
     &   incon(maxpl),idimn(maxpl),icomg(maxpl),
     &   idash(maxpl),ixavg(maxpl),icolr(maxpl),
     &   icoll(maxpl),ilcll(maxpl),ilchl(maxpl),
     &   rtslb(maxpl),rtshl(maxpl),imjsk(maxpl),
     &   icong(maxpl),iconl(maxpl),icozr(maxpl),ilwll(maxpl),
     &   ilwng(maxpl),ilwnl(maxpl),ilwzr(maxpl),idall(maxpl),
     &   idang(maxpl),idanl(maxpl),idazr(maxpl),
     &   ilcnl(maxpl),ilczr(maxpl),
     &   ilcbr(maxpl),ipwlb(maxpl),iorlb(maxpl),
     &   ipwhl(maxpl),ipwbr(maxpl),ifclb(maxpl),ifcnl(maxpl),
     &   ifczr(maxpl),ifchl(maxpl),
     &   ilclo(maxpl),ifclo(maxpl),rwdbr(maxpl),ihvbr(maxpl),
     &   incsq(maxpl),icosq(maxcosq,maxpl),rcosq(maxcosq,maxpl),
     &   fred(0:255),fgreen(0:255),fblue(0:255),
     &   work(miy,mjx,mkzh),icdwk(maxpl),
     &   pslab1(mabpl,morpl)
      dimension iam(niam),xcs(ncs),ycs(ncs)
      logical lmult(maxpl),larng(maxpl),lnolb(maxpl),
     &   lnobr(maxpl),lnozr(maxpl),lnmsg(maxpl)
      character unwk(maxpl)*24,ccmth(maxpl)*4,cfeld(3,maxpl)*10,
     &   cnohl(maxpl)*1
c
      dimension redcosq(500),greencosq(500),
     &   bluecosq(500)
      dimension rwrk(5000),iwrk(5000),rect(4),
     &   valcon(maxcon),majcon(maxcon),lbbarcon(maxcon),
     &   iaia(500),igia(500),lfin(maxcon+1)
      character messg*126,llbs(maxcon+2)*24
      external drawcl,cpcolr
c
      include 'comconst'
c
      dimension mconcp(1000),icoindcp(1000)
      common /cpack/ mconcp,icoindcp,nconarea,icpfchl,icpfclo,
     &   icpfclb,icpfcnl,icpfczr
c
      fb=fbmin
      ft=ftmax
      fl=flmin
      fr=frmax
c
c   Make set call for conrec.
c
      call set(fl,fr,fb,ft,1.,float(nscrs),set1,set2,1)
c
c   If field is Sawyer-Eliassen streamfunction, data is already
c   arranged into cross-section array.
c
      if (cfeld(1,ipl)(1:2).eq.'se'.or.cfeld(1,ipl)(1:2).eq.'sm') then
         do k=1,mkzh
         do ls=1,nscrs
            n1d=(k-1)*nscrs+ls
            kk=1+(n1d-1)/(miy*mjx)
            n1dleft=n1d-(kk-1)*miy*mjx
            jj=1+(n1dleft-1)/miy
            ii=n1dleft-(jj-1)*miy
            pslab1(ls,k)=work(ii,jj,kk)
         enddo
         enddo
         goto 431
      endif
c
c   Set up x-sec averaging parameters
c
      cosangle=xdist/xseclen
      sinangle=ydist/xseclen
      rnavg=1./(2.*ixavg(ipl)+1.)
c
      do 50 k=1,mkzh
      do 50 ls=1,nscrs
         pslab1(ls,k)=0.
   50 continue
c
c   Interpolate gridded data to x-section.
c
      caxgn=1.+(rcrag(2,ipl)-xjcorn)*refrat
      caygn=1.+(rcrag(1,ipl)-yicorn)*refrat
      cbxgn=1.+(rcrbg(2,ipl)-xjcorn)*refrat
      cbygn=1.+(rcrbg(1,ipl)-yicorn)*refrat
      do 103 islab=-ixavg(ipl),ixavg(ipl)
         xj1t=caxgn+islab*sinangle
         xj2t=cbxgn+islab*sinangle
         yi1t=caygn-islab*cosangle
         yi2t=cbygn-islab*cosangle
         if (xj1t.le.1.5.or.xj1t.ge.mjx-.5.or.
     &       xj2t.le.1.5.or.xj2t.ge.mjx-.5.or.
     &       yi1t.le.1.5.or.yi1t.ge.miy-.5.or.
     &       yi2t.le.1.5.or.yi2t.ge.miy-.5) then
            write(iup,*)'xj1t,xj2t,yi1t,yi2t=',xj1t,xj2t,yi1t,yi2t
            write(iup,*)'Cross sec. endpoints must be between 1.5'
            write(iup,*)'and (miy-.5) or (mjx-.5).'
            stop
         endif
c
         if (idimn(ipl).eq.2) then
            kend=1  ! only do the 2-D slab
         else
            kend=mkzh
         endif
c
      do 100 k=1,kend
         kp=k-kend+mkzh
      do 100 ls=1,nscrs
         posx=xj1t+(ls-1.)/(nscrs-1.)*(xj2t-xj1t) -
     &        icdwk(ipl)*.5
         posy=yi1t+(ls-1.)/(nscrs-1.)*(yi2t-yi1t) -
     &        icdwk(ipl)*.5
         jl=int(posx)
         jr=jl+1
         ib=int(posy)
         it=ib+1
         ratlr=posx-jl
         ratbt=posy-ib
         if (work(it,jl,k).eq.rmsg.or.
     &       work(it,jr,k).eq.rmsg.or.
     &       work(ib,jl,k).eq.rmsg.or.
     &       work(ib,jr,k).eq.rmsg.or.
     &       pslab1(ls,kp).eq.rmsg) then
            pslab1(ls,kp)=rmsg
         else
            wk1=work(it,jl,k)
            wk2=work(it,jr,k)
            wk3=work(ib,jl,k)
            wk4=work(ib,jr,k)
            pslab1(ls,kp)=pslab1(ls,kp)+ (
     &                   (1.-ratlr)*(   ratbt)*wk1+
     +                   (   ratlr)*(   ratbt)*wk2+
     +                   (1.-ratlr)*(1.-ratbt)*wk3+
     +                   (   ratlr)*(1.-ratbt)*wk4 )*rnavg
         endif
  100 continue
c
  103 continue
c
 431  continue
c
c   Smooth field if called for.
c
      call smooth(pslab1,ismth(ipl),mabpl,nscrs,mkzh)
c
c   Determine number of contours and contour values.
c
      call getconvals(rcbeg(ipl),rcend(ipl),rcint(ipl),
     &   incon(ipl),lmult(ipl),imjsk(ipl),pslab1,mabpl,nscrs,mkzh,rmsg,
     &   maxcon,valcon,majcon,cintuse,numcon,iup)
c
      call cpseti('SET',0)   ! we'll use our own set call
      call cpseti('CLS',0)   ! we'll use our own cont. values
      call cpseti('MAP',1)   ! enable coord. transformation via CPMPXY
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
      valmin=valcon(1)-.5*(valcon(2)-valcon(1))
      valmax=valcon(numcon)+.5*(valcon(numcon)-valcon(numcon-1))
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
      call cprect(pslab1,mabpl,nscrs,mkzh,rwrk,5000,iwrk,5000)
      call arinam(iam,niam)   ! initialize the area map
      call cpclam(pslab1,rwrk,iwrk,iam)  ! put cont. lines in area map
      call arscam(iam,xcs,ycs,ncs,iaia,igia,500,cpcolr)
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
c      the available space.  For bottom, assume .05 is needed for
c      tick and axis labels (typically present) and perhaps 3 additional fields
c      will be plotted.  For right, assume no space is needed.
c
         avlbottom=fb-bottextfloor-.05-3.*.0152
         avlright=1.-fr
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
     &      incon(ipl),lmult(ipl),lstep-1,pslab1,mabpl,nscrs,mkzh,rmsg,
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
      call cprect(pslab1,mabpl,nscrs,mkzh,rwrk,5000,iwrk,5000)
      if (iorlb(ipl).eq.1) then
         call cpseti('LLP',3)   ! penalty scheme for line labels
         call cpseti('LLO',0)   ! all labels at angle=0 (hor orientatn)
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
      call arinam(iam,niam)   ! initialize the area map
      call cplbam(pslab1,rwrk,iwrk,iam)  ! put label boxes in area map
      call cpcldm(pslab1,rwrk,iwrk,iam,drawcl)   ! draw contour lines
      call cplbdr(pslab1,rwrk,iwrk)   ! make labels
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
      call cpseti('MAP',0)   ! disable coord. transformation via CPMPXY
c
      return
      end