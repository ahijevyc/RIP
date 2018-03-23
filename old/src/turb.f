c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine turb(sigh,pstx,prs,uvel,vvel,ght,xmap,dmap,tt,qv,
     &   scr2a,scr2b,arr,miy,mjx,mkzh)
c
c compute the turbulence index, TI2, from Ellrod and Knapp, 
c W&F, March 1992.
c
      dimension sigh(mkzh),pstx(miy,mjx),prs(miy,mjx,mkzh),
     &   tt(miy,mjx,mkzh),ght(miy,mjx,mkzh),uvel(miy,mjx,mkzh),
     &   vvel(miy,mjx,mkzh),arr(miy,mjx,mkzh), qv(miy,mjx,mkzh),
     &   xmap(miy,mjx),dmap(miy,mjx),scr2a(miy,mjx),scr2b(miy,mjx)
      real d1(miy,mjx,mkzh), d2(miy,mjx,mkzh), d3(miy,mjx,mkzh),
     &    d4(miy,mjx,mkzh), d5(miy,mjx,mkzh)
c
      dimension wsp(mkzh)
c
      include 'comconst'
c
c     write(6,*) 'beginning turb, miy = ',miy,' mjx = ',mjx
c     write(6,*) 'mkzh = ',mkzh
      do 201 j=1,mjx-1
      do 200 i=1,miy-1
c
      do 150 k=1,mkzh
c        p = prs(i,j,k)
c
         wsp(k)=.25*sqrt((uvel(i,j,k)+uvel(i+1,j,k)+
     &                    uvel(i,j+1,k)+uvel(i+1,j+1,k))**2+
     &                   (vvel(i,j,k)+vvel(i+1,j,k)+
     &                    vvel(i,j+1,k)+vvel(i+1,j+1,k))**2)
c     if (i.eq.38.and.j.eq.38) write(6,*) 'k = ',k,' spd = ',
c    &   wsp(k)
  150 continue
c
      do 160 k=1,mkzh
         kp1=min(k+1,mkzh)
         km1=max(k-1,1)
         if (inhyd.eq.1) then
            refprspas=100.*(sigh(k)*pstx(i,j)+ptop) ! No pp needed.
            reftmk=refslt+reflaps*log(refprspas/1.e5)
            reftmk=max(reftmk,refstratt)
            refrho=refprspas/(rgas*reftmk)
            dsdz=-refrho*grav/(pstx(i,j)*100.)
         else
            dsdz=-(sigh(k)+ptop/pstx(i,j))*grav/(rgas*
     &         virtual(tt(i,j,k),qv(i,j,k)))
c              No pp needed here (hyd.)
         endif
c           shr=(wsp(kp1)-wsp(km1))/(sigh(kp1)-sigh(km1))*dsdz
            d5(i,j,k)=sqrt((uvel(i,j,km1)-uvel(i,j,kp1))**2 + 
     &  (vvel(i,j,km1)-vvel(i,j,kp1))**2)/(sigh(km1)-sigh(kp1))*dsdz
c     if (i.eq.38.and.j.eq.38) write(6,*) 'k = ',k,' shear = ',
c    &   d5(i,j,k),' old = ',shr
  160 continue
c
  200 continue
  201 continue
c     write(6,*) 'call first derivcz'
      call derivcz(uvel,0,ght,xmap,dmap,pstx,qv,tt,sigh,
     &   scr2a,scr2b,d1,0,'x',miy,mjx,mkzh)
c     write(6,*) 'call 2 derivcz'
      call derivcz(vvel,0,ght,xmap,dmap,pstx,qv,tt,sigh,
     &   scr2a,scr2b,d2,0,'y',miy,mjx,mkzh)
c     write(6,*) 'call 3 derivcz'
      call derivcz(vvel,0,ght,xmap,dmap,pstx,qv,tt,sigh,
     &   scr2a,scr2b,d3,0,'x',miy,mjx,mkzh)
c     write(6,*) 'call 4 derivcz'
      call derivcz(uvel,0,ght,xmap,dmap,pstx,qv,tt,sigh,
     &   scr2a,scr2b,d4,0,'y',miy,mjx,mkzh)
      do 302 k = 1, mkzh
      do 301 j = 1, mjx-1
      do 300 i = 1, miy-1
        arr(i,j,k) = d5(i,j,k) * (sqrt((d1(i,j,k)-d2(i,j,k))**2 + 
     &     (d3(i,j,k)+d4(i,j,k))**2) - (d1(i,j,k)+d2(i,j,k))) * 1.e7
c     if (i.eq.38.and.j.eq.38) write(6,*) 'k = ',k,' index = ',
c    &   arr(i,j,k)
  300 continue
  301 continue
  302 continue
      return
      end
