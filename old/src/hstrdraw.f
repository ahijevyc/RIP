c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hstrdraw(ilinw,iintv,sigf,workspc,vc3d,tmk,qvp,
     &   prs,ght,ter,pstx,prs_tsf,sigh,iprog,
     &   ixwin,iywin,ismth,icolr,lhide,
     &   rlevl,rlavl,cfeld,cvcor,idimn,work1,work2,icdwk,ilev,
     &   idiffflag,pslab1,pslab2,mabpl,morpl,maxlev,maxpl,
     &   miy,mjx,mkzh,ipl,rrota)
c
      dimension sigf(mkzh+1),workspc(miy,mjx,mkzh),vc3d(miy,mjx,mkzh),
     &   tmk(miy,mjx,mkzh),qvp(miy,mjx,mkzh),ght(miy,mjx,mkzh),
     &   ter(miy,mjx),pstx(miy,mjx),sigh(mkzh),prs(miy,mjx,mkzh),
     &   icolr(maxpl),idimn(maxpl),prs_tsf(miy,mjx),
     &   ixwin(2,maxpl),iywin(2,maxpl),ismth(maxpl),ilinw(maxpl),
     &   iintv(maxpl),rrota(maxpl),
     &   rlevl(maxlev,maxpl),rlavl(maxlev,maxpl),
     &   work1(miy,mjx,mkzh),work2(miy,mjx,mkzh),icdwk(maxpl),
     &   pslab1(mabpl,morpl),pslab2(mabpl,morpl)
      character cfeld(3,maxpl)*10,cvcor(maxpl)*1
      logical lhide(maxpl)
c
      include 'comconst'
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
      imsg=1            ! needed for strmln common block
      uvmsg=rmsg        ! needed for strmln common block
c
c   Set color and common block variable that controls line width,
c      and common block variable that controls density.
c
      write(6,*) 'HSTRDRAW'
      call gsplci(icolr(ipl))
      wdlv=float(ilinw(ipl))
      inita=iintv(ipl)
      initb=inita
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
      call set(fl,fr,fb,ft,ul,ur,ub,ut,1)
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
      call smooth(pslab1,ismth(ipl),mabpl,njx,niy)
      call smooth(pslab2,ismth(ipl),mabpl,njx,niy)
c
c   Call strmln
c
      call strmln(pslab1,pslab2,workspc,mabpl,njx,niy,1,istrerr)
      call setusv('LW',1000)
      call gsplci(1)
c
      return
      end
