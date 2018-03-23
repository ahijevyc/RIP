c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine vervcalc(pstd,dmap,sigf,xmap,pstx,sigh,
     &   uuu,vvv,qvp,tmk,www,prs,verv,ind,miy,mjx,mkzh)
c
      dimension pstd(miy,mjx),dmap(miy,mjx),sigf(mkzh+1),
     &   xmap(miy,mjx),pstx(miy,mjx),sigh(mkzh),uuu(miy,mjx,mkzh),
     &   vvv(miy,mjx,mkzh),qvp(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   www(miy,mjx,mkzh),prs(miy,mjx,mkzh),
     &   verv(miy,mjx,mkzh)
c
      dimension hmfd(100),sgd(100)
c
      include 'comconst'
c
c   Ind: (1) - w (dz/dt) in cm/s
c        (2) - omega (dp/dt) in dPa/s
c        (3) - sigmadot (d sigma/dt) in /s
c
      ds8=8.*ds
      ds2=2.*ds
c
      do 1000 j=1,mjx-1
c
      jp1=min(j+1,mjx-1)
      jm1=max(j-1,1)
      xjfact=2./(jp1-jm1)
c
      do 1000 i=1,miy-1
c
      ip1=min(i+1,miy-1)
      im1=max(i-1,1)
      yifact=2./(ip1-im1)
c
      if (inhyd.eq.0) then
c
c   Calculate pstar tendency and horizontal mass flux div.
c
      pstten=0.
      do 200 k=1,mkzh
       hmfd(k)=
     &  (pstd(i+1,j+1)*( uuu(i+1,j+1,k)+vvv(i+1,j+1,k))/dmap(i+1,j+1)+
     &   pstd(i  ,j+1)*( uuu(i  ,j+1,k)-vvv(i  ,j+1,k))/dmap(i  ,j+1)+
     &   pstd(i+1,j  )*(-uuu(i+1,j  ,k)+vvv(i+1,j  ,k))/dmap(i+1,j  )+
     &   pstd(i  ,j  )*(-uuu(i  ,j  ,k)-vvv(i  ,j  ,k))/dmap(i  ,j  ))/
     &   ds2
         pstten=pstten-hmfd(k)*(sigf(k+1)-sigf(k))*xmap(i,j)*xmap(i,j)
  200 continue
c
c   Calculate sigmadot
c
      sgd(1)=0.
      do 400 k=2,mkzh+1
         sgd(k)=sgd(k-1)-(pstten+xmap(i,j)*xmap(i,j)*hmfd(k-1))*
     &      (sigf(k)-sigf(k-1))/pstx(i,j)
  400 continue
c
c   Check sigmadot at ground is less than ~1 cm/s
c
      sgdcms=sgd(mkzh+1)*1.e6
      if (abs(sgdcms).gt.1.) then
         write(iup,*)'Surface sigmadot too big!  i,j,sgdcms= ',
     &      i,j,sgdcms
         stop
      endif
c
c   Interpolate sigmadot to half levels.
c
      do 420 k=1,mkzh
         verv(i,j,k)=.5*(sgd(k)+sgd(k+1))
  420 continue
c
      if (ind.eq.1.or.ind.eq.2) then
c
c      Calculate omega.
c
         do 700 k=1,mkzh
            verv(i,j,k)=(pstx(i,j)*verv(i,j,k)+sigh(k)*(pstten+(
     &       (uuu(i,j,k)+uuu(i+1,j,k)+uuu(i,j+1,k)+uuu(i+1,j+1,k))*
     &       (pstx(i,jp1)-pstx(i,jm1))*xjfact+
     &       (vvv(i,j,k)+vvv(i+1,j,k)+vvv(i,j+1,k)+vvv(i+1,j+1,k))*
     &       (pstx(ip1,j)-pstx(im1,j))*yifact)/ds8*xmap(i,j)))*
     &       1000.  !dPa/sec
  700    continue
      endif
c
      if (ind.eq.1) then
c
c      Calculate approximate w (ignoring slope and vertical velocity
c         of pressure surface). No pert. needed here (hydrostatic).
c
         do 710 k=1,mkzh
            verv(i,j,k)=-rgas*virtual(tmk(i,j,k),qvp(i,j,k))/
     &         (grav*(sigh(k)*pstx(i,j)+ptop))*verv(i,j,k)*
     &         .1  ! verv was in dPa/sec, want it in cm/sec
  710    continue
      endif
c
      elseif (inhyd.eq.1) then
c
c   w already exists.
c
      do 800 k=1,mkzh
         verv(i,j,k)=www(i,j,k)  ! Note: www is in cm/s
  800 continue
c
      if (ind.eq.2) then
c
c      Calculate approximate omega.
c
         do 810 k=1,mkzh
            verv(i,j,k)=-grav*prs(i,j,k)/
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k)))*verv(i,j,k)*
     &         10. ! verv was in cm/sec, want it in dPa/sec
  810    continue
      endif
c
      if (ind.eq.3) then
c
c      Calculate sigmadot.
c
         do 820 k=1,mkzh
            refprspas=100.*(sigh(k)*pstx(i,j)+ptop) ! No pert. needed.
            reftmk=refslt+reflaps*log(refprspas/1.e5)
            reftmk=max(reftmk,refstratt)
            refrho=refprspas/(rgas*reftmk)
            verv(i,j,k)=-refrho*grav*.01*verv(i,j,k)/(100.*pstx(i,j))-
     &         sigh(k)/pstx(i,j)*
     &         ((uuu(i,j,k)+uuu(i+1,j,k)+uuu(i,j+1,k)+uuu(i+1,j+1,k))*
     &          (pstx(i,jp1)-pstx(i,jm1))*xjfact+
     &          (vvv(i,j,k)+vvv(i+1,j,k)+vvv(i,j+1,k)+vvv(i+1,j+1,k))*
     &          (pstx(ip1,j)-pstx(im1,j))*yifact  )/ds8*xmap(i,j)
c       The ".01" is because verv was in cm/sec, so it should
c       first be converted to m/sec
  820    continue
      endif
c
      endif
c
 1000 continue
      return
      end
