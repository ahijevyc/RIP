c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine pvocalc(sigh,xmap,uuu,vvv,cor,theta,prs,pv,
     &   miy,mjx,mkzh)
c
      dimension sigh(mkzh),xmap(miy,mjx),uuu(miy,mjx,mkzh),
     &   vvv(miy,mjx,mkzh),cor(miy,mjx),theta(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),pv(miy,mjx,mkzh)
c
      include 'comconst'
c
      icn=0
      do 200 k=1,mkzh
c
      kp1=min(k+1,mkzh)
      km1=max(k-1,1)
      dsig=sigh(kp1)-sigh(km1)
         do 100 j=1,mjx-1
         jp1=min(j+1,mjx-1)
         jm1=max(j-1,1)
         do 100 i=1,miy-1
         ip1=min(i+1,miy-1)
         im1=max(i-1,1)
            dss=ds/xmap(i,j)
            dudy=.5*(uuu(i+1,j,k)+uuu(i+1,j+1,k)-
     &               uuu(i,j,k)-uuu(i,j+1,k))/dss
            dvdx=.5*(vvv(i+1,j+1,k)+vvv(i,j+1,k)-
     &               vvv(i+1,j,k)-vvv(i,j,k))/dss
            avort=dvdx-dudy+cor(i,j)
            dudsi=.25*(uuu(i,j,kp1)+uuu(i+1,j,kp1)+
     &                 uuu(i,j+1,kp1)+uuu(i+1,j+1,kp1)-
     &                 uuu(i,j,km1)-uuu(i+1,j,km1)-
     &                 uuu(i,j+1,km1)-uuu(i+1,j+1,km1))/dsig
            dvdsi=.25*(vvv(i,j,kp1)+vvv(i+1,j,kp1)+
     &                 vvv(i,j+1,kp1)+vvv(i+1,j+1,kp1)-
     &                 vvv(i,j,km1)-vvv(i+1,j,km1)-
     &                 vvv(i,j+1,km1)-vvv(i+1,j+1,km1))/dsig
            dthdsi=(theta(i,j,kp1)-theta(i,j,km1))/dsig
            dsidp=.01*dsig/(prs(i,j,kp1)-prs(i,j,km1))
            dy=dss*(ip1-im1)
            dx=dss*(jp1-jm1)
            dthdx=(theta(i,jp1,k)-theta(i,jm1,k))/dx
            dthdy=(theta(ip1,j,k)-theta(im1,j,k))/dy
            pvort=-grav*dsidp*(dthdsi*avort-dvdsi*dthdx+dudsi*dthdy)
c
c         Convert from si to familiar units:
c
            pv(i,j,k)=pvort*1.e6
c
  100    continue
  200 continue
      return
      end
