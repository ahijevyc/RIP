c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine contrivefg(cor,dmap,xmap,ter,pstd,pstx,prs,ght,
     &   tmk,qvp,uuu,vvv,www,sigh,rootname,iendcr,miy,mjx,mkzh)
c
c   This routine generates artificial model output fields for testing,
c   etc.
c
      dimension cor(miy,mjx),dmap(miy,mjx),xmap(miy,mjx),ter(miy,mjx),
     &   pstd(miy,mjx),pstx(miy,mjx),prs(miy,mjx,mkzh),
     &   qvp(miy,mjx,mkzh),
     &   ght(miy,mjx,mkzh),tmk(miy,mjx,mkzh),uuu(miy,mjx,mkzh),
     &   vvv(miy,mjx,mkzh),www(miy,mjx,mkzh),sigh(mkzh)
      character rootname*256
c
      include 'comconst'
c
      open (unit=49,file=rootname(1:iendcr)//'.ctin',
     &   form='formatted',status='old')
      read(49,*)graddir
      read(49,*)div
      read(49,*)f1
      read(49,*)f2
c
c   First set up terrain, xmap, dmap, cor
c
      do j=1,mjx-1
      do i=1,miy-1
         if (i.ge.33.and.i.le.58) then
            ter(i,j)=800.*(sin(pi*(i-33.)/25.))**2
         else
            ter(i,j)=0.
         endif
         if (i.ge.32.and.i.le.45.and.j.ge.34.and.j.le.59) then
            ter(i,j)=max(ter(i,j),800.*(sin(pi*(j-34.)/25.))**2)
         endif
         if (i.ge.20.and.i.le.32.and.j.ge.34.and.j.le.59) then
            ter(i,j)=800.*(sin(pi*(j-34.)/25.))**2*
     &                    (sin(pi*(i-19.5)/25.))**2
         endif
         xmap(i,j)=1.
         dmap(i,j)=1.
         cor(i,j)=1.e-4
      enddo
      enddo
      call xtodot(dmap,miy,mjx)
c
c   Next, define pstx, pstd based on reference state.
c
      cc1=rgas/grav*(-.5)*reflaps
      cc2=rgas/grav*(reflaps*log(.01*refslp)-refslt)
      cc3=rgas/grav*(refslt-.5*reflaps*log(.01*refslp))*log(.01*refslp)
c
      do j = 1, mjx-1
      do i = 1, miy-1
c
c      The following assumes that all terrain is below the reference
c         tropop.
c
         xxx=(-cc2-sqrt(cc2*cc2-4.*cc1*(cc3-ter(i,j))))/(2.*cc1)
         pstx(i,j)=exp(xxx)-ptop
         pstd(i,j)=pstx(i,j)
      enddo
      enddo
c
      call xtodot(pstd,miy,mjx)
c
c   Next define 3D fields.
c
      xjmid=.5*(1.+mjx)
      yimid=.5*(1.+miy)
      gang=(90.-graddir)*rpd
      width=500. ! km
      hwidth=width/2.
      widthtran=200. ! km
      hwplust=hwidth+widthtran
      thsealevmid=290. ! K
      dthdxsealev=2. ! K/100km
      dthdz=4. ! K/km
      do k=1,mkzh
      do j = 1, mjx-1
      do i = 1, miy-1
         alnpref=log(sigh(k)*pstx(i,j)+ptop)
         ght(i,j,k)=cc1*alnpref*alnpref+cc2*alnpref+cc3
         dthdx=dthdxsealev*(12000.-ght(i,j,k))/12000.
         dthmax=.01*dthdx/widthtran*(hwplust*(hwplust-hwidth)+
     &       .5*(hwidth*hwidth-hwplust*hwplust)+hwidth*widthtran)
         distnzone=((j+.5-xjmid)*cos(gang)+(i+.5-yimid)*sin(gang))*dskm
         adz=abs(distnzone)
         if (adz.le.hwidth) then
            dtheta=adz*.01*dthdx
         elseif(adz.le.hwidth+widthtran) then
            dtheta=.01*dthdx/widthtran*(hwplust*(adz-hwidth)+
     &         .5*(hwidth*hwidth-adz*adz)+hwidth*widthtran)
         else
            dtheta=dthmax
         endif
         thetamid=thsealevmid+ght(i,j,k)*.001*dthdz
         if (distnzone.ge.0.) then
            theta=thetamid+dtheta
            www(i,j,k)=dtheta
         else
            theta=thetamid-dtheta
            www(i,j,k)=-dtheta
         endif
         prs(i,j,k)=sigh(k)*pstx(i,j)+ptop
         tmk(i,j,k)=theta*(prs(i,j,k)/1000.)**gamma
         es = ezero * exp( eslcon1*(tmk(i,j,k)-celkel)/
     &      (tmk(i,j,k)-eslcon2) )
         qvp(i,j,k)=eps*es/(prs(i,j,k)-es)
         gammam=gamma*(1.+gammamd*qvp(i,j,k))
         tmk(i,j,k)=theta*(prs(i,j,k)/1000.)**gammam
         uuu(i,j,k)=(.5*(div+f1)*(j-xjmid)+.5*f2*(i-yimid))*ds
         vvv(i,j,k)=(.5*f2*(j-xjmid)+.5*(div-f1)*(i-yimid))*ds
      enddo
      enddo
      enddo
      do j = 1, mjx
      do i = 1, miy
      do k=1,mkzh
         uuu(i,j,k)=(.5*(div+f1)*(j-xjmid)+.5*f2*(i-yimid))*ds
         vvv(i,j,k)=(.5*f2*(j-xjmid)+.5*(div-f1)*(i-yimid))*ds
      enddo
      enddo
      enddo
c
      return
      end
