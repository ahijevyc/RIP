c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hvecdraw(ilinw,sigf,vc3d,tmk,qvp,
     &   prs,ght,ter,pstx,prs_tsf,sigh,iprog,
     &   icolr,ixwin,iywin,ismth,rvcmx,cfulb,unwk,lhide,icomg,
     &   iintv,rlevl,rlavl,cfeld,cvcor,idimn,idiffflag,
     &   work1,work2,icdwk,ilev,
     &   lnmsg,bottextfloor,pslab1,pslab2,mabpl,morpl,maxlev,
     &   maxpl,miy,mjx,mkzh,ipl,rrota)
c
      dimension sigf(mkzh+1),vc3d(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   qvp(miy,mjx,mkzh),ght(miy,mjx,mkzh),prs(miy,mjx,mkzh),
     &   ter(miy,mjx),pstx(miy,mjx),sigh(mkzh),prs_tsf(miy,mjx),
     &   ixwin(2,maxpl),idimn(maxpl),icomg(maxpl),
     &   iywin(2,maxpl),ismth(maxpl),ilinw(maxpl),rvcmx(maxpl),
     &   iintv(maxpl),rlevl(maxlev,maxpl),rlavl(maxlev,maxpl),
     &   work1(miy,mjx,mkzh),work2(miy,mjx,mkzh),icdwk(maxpl),
     &   icolr(maxpl),pslab1(mabpl,morpl),pslab2(mabpl,morpl),
     &   rrota(maxpl),tempu(mabpl,morpl),tempv(mabpl,morpl),
     &   tempu2(morpl,mabpl),tempv2(morpl,mabpl)
      logical lnmsg(maxpl),lhide(maxpl)
      character cfeld(3,maxpl)*10,cvcor(maxpl)*1,cfulb(maxpl)*5,
     &   unwk(maxpl)*24
c
      dimension vecskip(2)
      character string*48
      real cv,ch
c
      include 'comconst'
c
      vecskip(1)=rmsg
      vecskip(2)=rmsg
c
c   Set line width
c
      lwidth=ilinw(ipl)*1000
      call setusv('LW',lwidth)
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
c MGD begin mod
c if we are dealing with a (-)90 degree rotation, then the corners
c of the set() window should reflect a rotated rectangle (accomodate
c the reversal of dimensions)
      if(rrota(ipl) .ne. 90. .and. rrota(ipl) .ne. -90.) then
        call set(fl,fr,fb,ft,ul,ur,ub,ut,1)
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
         do 90 j=1,njx
            jj=j+ixwin(1,ipl)-1
            do 90 i=1,niy
               ii=i+iywin(1,ipl)-1
               pslab1(j,i)=work1(ii,jj,nint(rlevl(ilev,ipl)))
               pslab2(j,i)=work2(ii,jj,nint(rlevl(ilev,ipl)))
   90    continue
      elseif (cvcor(ipl).eq.'s'.and.rlevl(ilev,ipl).ne.
     &      rlavl(ilev,ipl).and.rlavl(ilev,ipl).ge.0) then
         call fillarray(pslab1,mabpl*morpl,0.)
         call fillarray(pslab2,mabpl*morpl,0.)
         lev1=min(nint(rlevl(ilev,ipl)),nint(rlavl(ilev,ipl)))
         lev2=max(nint(rlevl(ilev,ipl)),nint(rlavl(ilev,ipl)))
         sigtot=sigf(lev2+1)-sigf(lev1)
         do 125 k=lev1,lev2
         do 120 j=1,njx
            jj=j+ixwin(1,ipl)-1
            do 120 i=1,niy
               ii=i+iywin(1,ipl)-1
               pslab1(j,i)=pslab1(j,i)+work1(ii,jj,k)*
     &          (sigf(k+1)-sigf(k))/sigtot
               pslab2(j,i)=pslab2(j,i)+work2(ii,jj,k)*
     &          (sigf(k+1)-sigf(k))/sigtot
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
               pslab1(j,i)=work1(ii,jj,lev1)-work1(ii,jj,lev2)
               pslab2(j,i)=work2(ii,jj,lev1)-work2(ii,jj,lev2)
            enddo
         enddo
      else
         call vinterp(cvcor(ipl),rlevl(ilev,ipl),ixwin(1,ipl),
     &      iywin(1,ipl),icdwk(ipl),vc3d,tmk,qvp,
     &      prs,ght,ter,pstx,sigh,sigf,prs_tsf,
     &      lhide(ipl),idiffflag,cfeld(1,ipl),work1,
     &      pslab1,iprog,mabpl,morpl,njx,niy,miy,mjx,mkzh)
         call vinterp(cvcor(ipl),rlevl(ilev,ipl),ixwin(1,ipl),
     &      iywin(1,ipl),icdwk(ipl),vc3d,tmk,qvp,
     &      prs,ght,ter,pstx,sigh,sigf,prs_tsf,
     &      lhide(ipl),idiffflag,cfeld(2,ipl),work2,
     &      pslab2,iprog,mabpl,morpl,njx,niy,miy,mjx,mkzh)
      endif

c MGD begin mod
c ahh, rotating arrays of vectors... 
c for cw rotation, each of the u and v arrays is transposed, then
c inverted vertically. but, for such a rotation, the v component
c must be negated to reflect the fact that positive x axis has become
c negative y axis 
      if(rrota(ipl) .eq. -90.) then
        do i=1,njx
        do j=1,niy
          tempu2(j,njx-i+1) = pslab2(i,j)
          tempv2(j,njx-i+1) = -1*pslab1(i,j)
        enddo
        enddo
c similar situation, but with things reversed a bit
      elseif(rrota(ipl) .eq. 90.) then
        do i=1,njx
        do j=1,niy
          tempu2(niy-j+1,i) = -1*pslab2(i,j)
          tempv2(niy-j+1,i) = pslab1(i,j)
c         if(sqrt(tempu2(niy-j+1,i)**2 + tempv2(niy-j+1,i)**2)
c    &        .gt. 200) 
c    &     write(6,*) 'BIGGIE: ',niy-j+1,i,tempu2(niy-j+1,i),
c    &     tempv2(niy-j+1,i) 
        enddo
        enddo
c for 180 degree rotations, need to invert both x and y of both 
c arrays, and also negate the values of the components.
      elseif(rrota(ipl) .eq. -180. .or. rrota(ipl) .eq. 180.) then
        do i=1,njx
        do j=1,niy
          tempu(njx-i+1,niy-j+1) = -1*pslab1(i,j)
          tempv(njx-i+1,niy-j+1) = -1*pslab2(i,j)
        enddo
        enddo
      endif

c dimensions for 180 degree rotated fields will be the same, so keep
c the data in old array
      if(rrota(ipl) .eq. 180. .or. rrota(ipl) .eq. -180.) then
      do i=1,njx
      do j=1,niy
        pslab1(i,j) = tempu(i,j)
        pslab2(i,j) = tempv(i,j)
      enddo
      enddo
      endif
c MGD end mod

c      print*,'Enter C,D:'
c      read*,ccc,ddd
c      do j=1,njx
c      do i=1,niy
c         xx=float(j-njx/2)
c         yy=float(i-niy/2)
c         pslab1(j,i)=(ddd-ccc)*xx
c         pslab2(j,i)=(-ddd-ccc)*yy
c      enddo
c      enddo
c
c   Smooth data if necessary
c
c MGD begin mod
c niy and njx arguements used to be miy and mjx, respectively, but
c this caused some problems when smoothing rotated vector fields,
c especially when intv was >= 2
      if(rrota(ipl) .ne. 90. .and. rrota(ipl) .ne. -90.) then
        call smooth(pslab1,ismth(ipl),mabpl,njx,niy)
        call smooth(pslab2,ismth(ipl),mabpl,njx,niy)
      else
        call smooth(tempu2,ismth(ipl),morpl,niy,njx)
        call smooth(tempv2,ismth(ipl),morpl,niy,njx)
      endif
c MGD end mod
c
c   Put in special values where we don't want vectors
c
      vmagmax=0.
      intp=iintv(ipl)
      if (intp.gt.0) then  ! regular staggering
         intph=0
      else   ! "diamond pattern" staggering
         intp=-intp
         intph=intp/2
      endif
      if(rrota(ipl) .ne. 90. .and. rrota(ipl) .ne. -90.) then
      do j=1,mjx
      do i=1,miy
         if ( ((mod(i-1,intp).eq.0.and.mod(j-1,intp).eq.0).or.
     &         (mod(i-1+intph,intp).eq.0.and.
     &          mod(j-1+intph,intp).eq.0)) .and.
     &       pslab1(j,i).ne.rmsg.and.pslab2(j,i).ne.rmsg) then
            vmag=sqrt(pslab1(j,i)**2+pslab2(j,i)**2)
            vmagmax=max(vmagmax,vmag)
         else
            pslab1(j,i)=rmsg
            pslab2(j,i)=rmsg
         endif
      enddo
      enddo
c MGD begin mod
c if we have rotated (-)90, do same as we normally would, but with
c our own array of reversed dimensions.
      else
      do j=1,miy
      do i=1,mjx
         if ( ((mod(i-1,intp).eq.0.and.mod(j-1,intp).eq.0).or.
     &         (mod(i-1+intph,intp).eq.0.and.
     &          mod(j-1+intph,intp).eq.0)) .and.
     &       tempu2(j,i).ne.rmsg.and.tempv2(j,i).ne.rmsg) then
            vmag=sqrt(tempu2(j,i)**2+tempv2(j,i)**2)
            vmagmax=max(vmagmax,vmag)
         else
            tempu2(j,i)=rmsg
            tempv2(j,i)=rmsg
         endif
      enddo
      enddo
      endif
c MGD end mod
c
c   If vmagmax=0, then why the hell are we doing this?
c
      if (vmagmax.eq.0.) goto 200
c
      call getusv('XF',ixpau)
      ixpau=2**ixpau
      gskip=float(abs(iintv(ipl)))
      if (iintv(ipl).lt.0) gskip=gskip*.7071
      rpaubetv=(fr-fl)/xintervs*ixpau*gskip
      if (rvcmx(ipl).gt.0) then
         rpaupervmag=rpaubetv/rvcmx(ipl)
      elseif (rvcmx(ipl).eq.0) then
         rpaupervmag=rpaubetv/vmagmax
      endif
      if (rvcmx(ipl).ge.0) then
         rpaumax=rpaupervmag*vmagmax
         npaumax=int(rpaumax)+1
         vmagpmax=vmagmax*npaumax/rpaumax
      endif

c
c   Call velvct
c
      if (rvcmx(ipl).ge.0) then   ! vectors
         imxvpl=0
c MGD use whichever array is appropriate for the rotation at hand
      if(rrota(ipl) .ne. 90. .and. rrota(ipl) .ne. -90.) then
         call velvctmts(pslab1,mabpl,pslab2,mabpl,mjx,miy,
     &         0.,vmagpmax,1,npaumax,4,vecskip,imxvpl,icolr(ipl))
      else 
c if (-)90 rotation, then need to use our arrays with interchanged 
c dimensions...
         call velvctmts(tempu2,morpl,tempv2,morpl,miy,mjx,
     &         0.,vmagpmax,1,npaumax,4,vecskip,imxvpl,icolr(ipl))
      endif
c
         if (.not.lnmsg(ipl)) then
            call setusv('LW',1000)
            call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
            call gsplci(icomg(ipl))
            call gstxci(icomg(ipl))
            chsize=.008
            ypos=bottextfloor+.5*chsize
            xleft=.3
            string=' '
            write(string,'(a16,f5.1,1x,a24)')
     &         'MAXIMUM VECTOR: ',vmagmax,unwk(ipl)
            call pcseti('TE',1)
            call pcgeti ('QU',ntextqq)
            call pcseti ('QU',0)
            call plchhq(xleft,ypos,string,chsize,0.,-1.)
            call pcseti ('QU',ntextqq)
            call pcgetr('DR',distr)
            call pcseti('TE',0)
            vecmaxfcor=vmagmax*rpaupervmag/ixpau
            xstart=xleft+distr+.04
            xend=xstart+vecmaxfcor
            call setusv('LW',lwidth)
            call line(xstart,ypos,xend,ypos)
            dxarrow=.24*vecmaxfcor*cos(.45)
            dyarrow=.24*vecmaxfcor*sin(.45)
            call line(xend,ypos,xend-dxarrow,ypos+dyarrow)
            call line(xend,ypos,xend-dxarrow,ypos-dyarrow)
            call gsplci(1)
            call gstxci(1)
            bottextfloor=bottextfloor+1.9*chsize
         endif
c
      else                  ! wind barbs
c
c      Convert to appropriate units so that a full barb represents
c         the desired magnitude.   The velbrb routine always assumes
c         that a full barb = 10 units.  Hence, if cfulb=10mps, the
c         conversion factor is 1.  If cfulb=5mps, the conversion factor
c         is 2.  If cfulb=10kts, the conversion factor is 1.94.
c
         string=' '
         if (index(cfulb(ipl),'10mps').ne.0) then
            barbfac=1.
            write(string,'(a41)')
     &         'BARB VECTORS:  FULL BARB = 10 m s~S~-1~N~'
            nch=41
         elseif (index(cfulb(ipl),'5mps').ne.0) then
            barbfac=2.
            write(string,'(a40)')
     &         'BARB VECTORS:  FULL BARB = 5 m s~S~-1~N~'
            nch=40
         else
            barbfac=rktpmps
            write(string,'(a33)')
     &         'BARB VECTORS:  FULL BARB = 10 kts'
            nch=33
         endif
c
c MGD begin mod
c basically same as usual, but set up two parallel sections - uses
c the one that utilizes the appropriate array (depending on rotation)
      if(rrota(ipl) .ne. 90. .and. rrota(ipl) .ne. -90.) then
         do 140 j=1,mjx
         do 140 i=1,miy
            if (pslab1(j,i).ne.vecskip(1))
     &         pslab1(j,i)=barbfac*pslab1(j,i)
            if (pslab2(j,i).ne.vecskip(2))
     &         pslab2(j,i)=barbfac*pslab2(j,i)
  140    continue
         call gsplci(icolr(ipl))
         call gstxci(icolr(ipl))
         call velbrb(pslab1,mabpl,pslab2,mabpl,mjx,miy,
     &       .4*gskip,4,vecskip,iywin(1,ipl),ixwin(1,ipl),0,rrota(ipl))
         else
         do 240 j=1,mjx
         do 240 i=1,miy
            if (tempu2(j,i).ne.vecskip(1))
     &         tempu2(j,i)=barbfac*tempu2(j,i)
            if (tempv2(j,i).ne.vecskip(2))
     &         tempv2(j,i)=barbfac*tempv2(j,i)
  240    continue
         call gsplci(icolr(ipl))
         call gstxci(icolr(ipl))
         call velbrb(tempu2,morpl,tempv2,morpl,miy,mjx,
     &       .4*gskip,4,vecskip,ixwin(1,ipl),iywin(1,ipl),0,rrota(ipl))
         endif
c MGD end mod
         if (.not.lnmsg(ipl)) then
            call setusv('LW',1000)
            call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
            chsize=.008
            ypos=bottextfloor+.5*chsize
            call pcgeti ('QU',ntextqq)
            call pcseti ('QU',0)
            call gsplci(icomg(ipl))
            call gstxci(icomg(ipl))
            call plchhq(.5,ypos,string(1:nch),chsize,0.,0.)
            call pcseti ('QU',ntextqq)
            bottextfloor=bottextfloor+1.9*chsize
         endif
         call gsplci(1)
         call gstxci(1)
      endif
c
  200 continue
      call setusv('LW',1000)
c
      return
      end
