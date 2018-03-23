c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine pvocalc(xmap,uuu,vvv,cor,theta,prs,pv,
     &   www,qv,miy,mjx,mkzh,xtime,ght,ter)
      parameter (pvtr=1.5,thtc=310.,mm=20,k1=29,ibox=1,
     &  k3=17,k2=24,kw=21,rd=287.,p0=1.e3,zpi=3.14159,lv1=2.5e6,
     &  tht1=350.,tht2=370.,i1=90,i2=150,j1=170,j2=290,
C     &  u_domspeed=14.,klevmax=13,iv=2)  !5 July
C     &  u_domspeed=10.,klevmax=13,iv=2)   !SI run
C     &  u_domspeed=8.5,klevmax=13,iv=2)   !Phys run
     &  u_domspeed=0.,klevmax=13,iv=2)   !Phys run
      parameter (dtmkmax=10.,dtmkmin=-10.)
      parameter (dqvmax=0.005,dqvmin=-0.005)
c
      dimension xmap(miy,mjx),uuu(miy,mjx,mkzh),
     &   vvv(miy,mjx,mkzh),cor(miy,mjx),theta(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),pv(miy,mjx,mkzh),www(miy,mjx,mkzh),
     &   pvd(miy,mjx,mkzh),tmk3(miy,mjx,mkzh),
     &   qv(miy,mjx,mkzh),vor(miy,mjx,mkzh),prd(miy,mjx,mkzh)
      dimension up(miy,mjx,klevmax),vp(miy,mjx,klevmax),
     &   wp(miy,mjx,klevmax),vorp(miy,mjx,klevmax),
     &   thp(miy,mjx,klevmax)
      dimension pb(miy),pvar(miy,10),div(miy,mjx,mkzh),
     &   prtrh(mjx),vorh(mjx),itop(miy)
      dimension pvx(mkzh),pvy(mkzh),pvav(mkzh),prav(mkzh),
     &   pvth(miy,mjx),pvrad(mjx/2),divrad(miy/2),icrad(mjx/2),
     &   wrad(miy/2),pvav1(miy),vav(miy),wx(mkzh),wy(mkzh),
     &   wav(mkzh),qrx(mkzh),qry(mkzh),qrav(mkzh),
     &   qrav1(36),wav1(36),icav1(36),prlev(klevmax)
      real pstr(mkzh),ptil(mkzh),nstr(mkzh),ntil(mkzh)
      integer pstrd(mkzh,100),ptild(mkzh,100),nstrd(mkzh,100),
     &   ntild(mkzh,100),nptend(mkzh),nntend(mkzh),
     &   nptends(mkzh),nptendt(mkzh),nntends(mkzh),nntendt(mkzh)
c
c---below is where the new kbmin3 array is dimensioned (ST) 
c
      integer ipvor(mkzh),invor(mkzh),kbmin3(miy,mjx,mkzh)
      integer ioutc(miy,mjx,mkzh),ioutd(miy,mjx,mkzh)
      dimension dv1x(miy,mjx),dv2x(miy,mjx),dv1y(miy,mjx),dv2y(miy,mjx)
      real disv1(mkzh,200),disv2(mkzh,200)
      real t1(mkzh),p1(mkzh),q1(mkzh)
C      real ubav1(mkzh),vbav1(mkzh)
      integer kbmin(miy,mjx),kmsemax(miy,mjx)
      real hfx(miy,mjx),qfx(miy,mjx),pblh(miy,mjx),zmsemax(miy,mjx),
     &     zbmin(miy,mjx),tbmin(miy,mjx),tpmin(miy,mjx),tlcl(miy,mjx),
     &     ter(miy,mjx),ght(miy,mjx,mkzh)
      character getfromcore(20)*10

      data prlev /350, 400, 450, 500, 550, 600, 650, 700, 750, 800,
     &  850, 900, 950/
      data (getfromcore(i),i=1,2) /'ght       ','ter       '/
C
      integer icav(miy)
c
      include 'comconst'
c

      open(unit=50,file='bamex3_0530_pvo_all_kbmin.out',status='new')
      do kboss =20,mkzh

      read(16)kbmin
      read(16)tbmin
      read(16)zbmin
      read(16)tlcl
      read(16)tpmin
      read(17)zmsemax
      read(17)kmsemax
      read(17)ter
      close(16)
      close(17)
      close(18)
      read(18)hfx
      print*,hfx(miy-1,mjx-1)
      read(18)qfx
      print*,qfx(miy-1,mjx-1)
      read(18)pblh
      print*,pblh(miy-1,mjx-1)
      close(18)
      icn=0
c
c---initialize new kbmin3 array to be used in future budget calcs (ST)
c
      do i=1,miy
        do j=1,mjx
          do k=1,mkzh
            kbmin3(i,j,k)=0
          end do
        end do
      end do
c
      do k=1,mkzh
       pstr(k)=0.
       nstr(k)=0.
       ptil(k)=0.
       ntil(k)=0.
       ipvor(k)=0.
       invor(k)=0.
       nptend(k)=0
       nntend(k)=0
       nptends(k)=0
       nntends(k)=0
       nptendt(k)=0
       nntendt(k)=0
C       ubav1(k)=0.
C       vbav1(k)=0.
       do ibin=1,100
        pstrd(k,ibin)=0
        nstrd(k,ibin)=0
        ptild(k,ibin)=0
        ntild(k,ibin)=0
       enddo
      do i=1,miy
      do j=1,mjx
       ioutc(i,j,k)=0
       ioutd(i,j,k)=0
       ! set kmsemax = kboss
       kmsemax(i,j) = kboss
       ! tie zmsemax to kboss
       zmsemax(i,j) = ght(i,j,kboss)
      enddo
      enddo
      enddo
      do 200 k=1,mkzh
       do i=1,200
        disv1(k,i)=0.
        disv2(k,i)=0.
       enddo
c
      kp1=min(k+1,mkzh)
      km1=max(k-1,1)
      pav=0.
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
            dudx=.5*(uuu(i,j+1,k)+uuu(i+1,j+1,k)-
     &               uuu(i,j,k)-uuu(i+1,j,k))/dss
            dvdy=.5*(vvv(i+1,j,k)+vvv(i+1,j+1,k)-
     &               vvv(i,j,k)-vvv(i,j+1,k))/dss
            avort=dvdx-dudy+cor(i,j)
            dp=prs(i,j,kp1)-prs(i,j,km1)
            dudp=.25*(uuu(i,j,kp1)+uuu(i+1,j,kp1)+
     &                uuu(i,j+1,kp1)+uuu(i+1,j+1,kp1)-
     &                uuu(i,j,km1)-uuu(i+1,j,km1)-
     &                uuu(i,j+1,km1)-uuu(i+1,j+1,km1))/dp
            dvdp=.25*(vvv(i,j,kp1)+vvv(i+1,j,kp1)+
     &                vvv(i,j+1,kp1)+vvv(i+1,j+1,kp1)-
     &                vvv(i,j,km1)-vvv(i+1,j,km1)-
     &                vvv(i,j+1,km1)-vvv(i+1,j+1,km1))/dp
            dthdp=(theta(i,j,kp1)-theta(i,j,km1))/dp
            dy=dss*(ip1-im1)
            dx=dss*(jp1-jm1)
            dthdx=(theta(i,jp1,k)-theta(i,jm1,k))/dx
            dthdy=(theta(ip1,j,k)-theta(im1,j,k))/dy
            pvort=-grav*(dthdp*avort-dvdp*dthdx+dudp*dthdy)
            div(i,j,k)=1.e4*(dudx+dvdy)
            vor(i,j,k)=1.e4*(dvdx-dudy+cor(i,j))
            pav=pav+ prs(i,j,k)/float(miy)/float(mjx)
c cp comes in as 1004.000
            tmk3(i,j,k)=theta(i,j,k)*((prs(i,j,k)/p0)**(rd/cp))
c
c         Convert to PVU: 1.e-2 to go from "per hPa" to "per Pa" (SI),
c         and 1.e6 to go from SI to PVU
c
            pv(i,j,k)=pvort*1.e4
c
  100    continue
C         write(43,*)k,pav
  200 continue
c
      do 270 j=1,mjx
      do 270 i=1,miy
      do 270 k=1,18
       pv(i,j,k)=0.
 270  continue
      tdiffmax=0.
      do 271 j=3,mjx-2
      do 271 i=3,miy-2
       ibeg=max(i-ibox,3)
       iend=min(i+ibox,miy-2)
       jbeg=max(j-ibox,3)
       jend=min(j+ibox,mjx-2)
       pblhav=0.
       hfxav=0.
       qfxav=0.
       rhoav=0.
       iibl=0
       isumt=0
       if ((i.eq.3).and.(j.eq.3)) then
        print*,kmsemax(i,j),zmsemax(i,j),pblh(i,j),ter(i,j)
       endif
       do k=1,mkzh
        t1(k)=tmk3(i,j,k)
        q1(k)=qv(i,j,k)
        p1(k)=prs(i,j,k)
       enddo
       kk=kmsemax(i,j)
       dt=0.
       dqv=0.
C this call to tlift is when bmin (tl) amd bmax (tm) are retrieved.
       if (dtmk .ne. 0) print*,'i,j,dtmk=',i,j,dtmk
       call tlift(t1,dtmk,zmsemax(i,j),p1,q1,dqv,tl,tm,i,j,kk,mkzh,kb)
c
c---setting levels of kbmin for future budget calculations (ST) 
c
       kbmin3(i,j,kboss)=kb
c
       kk=kmsemax(i,j)
       do 274 ii=ibeg,iend
       do 274 jj=jbeg,jend
        kkmax=min(kmsemax(i,j)+1,mkzh)
        kkmin=max(kmsemax(i,j)-1,1)
        dss=ds/xmap(ii,jj)
        dens=1.e2*prs(ii,jj,kk)/(rd*tmk3(ii,jj,kk))
        omega=-0.01*www(ii,jj,kk)*grav*dens  !hydrostatic approx: OK for averaging
        ucav=0.25*(uuu(ii+1,jj,kk)+uuu(ii,jj+1,kk)+uuu(ii+1,jj+1,kk)+
     &    uuu(ii,jj,kk))
        vcav=0.25*(vvv(ii+1,jj,kk)+vvv(ii,jj+1,kk)+vvv(ii+1,jj+1,kk)+
     &    vvv(ii,jj,kk))
C
        tadvx=-ucav*(tmk3(ii,jj+1,kk)-tmk3(ii,jj-1,kk))/(2.*dss)
        tadvy=-vcav*(tmk3(ii+1,jj,kk)-tmk3(ii-1,jj,kk))/(2.*dss)
        dpdx=1.e2*(prs(ii,jj+1,kk)-prs(ii,jj-1,kk))/(2.*dss)
        dpdy=1.e2*(prs(ii+1,jj,kk)-prs(ii-1,jj,kk))/(2.*dss)
        dtdp=0.01*(tmk3(ii,jj,kkmax)-tmk3(ii,jj,kkmin))/
     &     (prs(ii,jj,kkmax)-prs(ii,jj,kkmin))
        tadvcorx=ucav*dtdp*dpdx
        tadvcory=vcav*dtdp*dpdy
        tadvp=-omega*dtdp
        tadvad=-0.01*www(ii,jj,kk)*grav/cp
C
        qadvx=-ucav*(qv(ii,jj+1,kk)-qv(ii,jj-1,kk))/(2.*dss)
        qadvy=-vcav*(qv(ii+1,jj,kk)-qv(ii-1,jj,kk))/(2.*dss)
        dqdp=0.01*(qv(ii,jj,kkmax)-qv(ii,jj,kkmin))/
     &     (prs(ii,jj,kkmax)-prs(ii,jj,kkmin))
        qadvp=-omega*dqdp
        qadvcorx=ucav*dqdp*dpdx
        qadvcory=vcav*dqdp*dpdy
C
        pblhav=pblhav + pblh(ii,jj)-ter(ii,jj)
        hfxav=hfxav + hfx(ii,jj)
        qfxav=qfxav + qfx(ii,jj)
        rhoav=rhoav + dens
        if (zmsemax(ii,jj).le.pblh(ii,jj)+ter(ii,jj)) then
         iibl=iibl + 1
        endif
        pv(i,j,1)=pv(i,j,1) + tadvx+tadvy+tadvcorx+tadvcory
        pv(i,j,2)=pv(i,j,2) + tadvp+tadvad
        pv(i,j,3)=pv(i,j,3) + (qadvx+qadvy+qadvcorx+qadvcory)*lv1/cp
        pv(i,j,4)=pv(i,j,4) + qadvp*lv1/cp
C Now compute advection of T at level of max neg buoyancy
        kkmaxb=min(kb+1,mkzh)
        kkminb=max(kb-1,1)
        dens=1.e2*prs(ii,jj,kb)/(rd*tmk3(ii,jj,kb))
        omega=-0.01*www(ii,jj,kb)*grav*dens  !hydrostatic approx: OK for averaging
        ucav=0.25*(uuu(ii+1,jj,kb)+uuu(ii,jj+1,kb)+uuu(ii+1,jj+1,kb)+
     &    uuu(ii,jj,kb))
        vcav=0.25*(vvv(ii+1,jj,kb)+vvv(ii,jj+1,kb)+vvv(ii+1,jj+1,kb)+
     &    vvv(ii,jj,kb))
        tadvx=-ucav*(tmk3(ii,jj+1,kb)-tmk3(ii,jj-1,kb))/(2.*dss)
        tadvy=-vcav*(tmk3(ii+1,jj,kb)-tmk3(ii-1,jj,kb))/(2.*dss)
        dpdx=1.e2*(prs(ii,jj+1,kb)-prs(ii,jj-1,kb))/(2.*dss)
        dpdy=1.e2*(prs(ii+1,jj,kb)-prs(ii-1,jj,kb))/(2.*dss)
        dtdp=0.01*(tmk3(ii,jj,kkmaxb)-tmk3(ii,jj,kkminb))/
     &     (prs(ii,jj,kkmaxb)-prs(ii,jj,kkminb))
        tadvcorx=ucav*dtdp*dpdx
        tadvcory=vcav*dtdp*dpdy
        tadvp=-omega*dtdp
        tadvad=-0.01*www(ii,jj,kb)*grav/cp
C
        pv(i,j,13)=pv(i,j,13) + tadvx+tadvy+tadvcorx+tadvcory
        pv(i,j,14)=pv(i,j,14) + tadvp+tadvad
        isumt=isumt + 1
 274   continue
       pv(i,j,1)=3.6e3*pv(i,j,1)/float(isumt)
       pv(i,j,2)=3.6e3*pv(i,j,2)/float(isumt)
       pv(i,j,3)=3.6e3*pv(i,j,3)/float(isumt)
       pv(i,j,4)=3.6e3*pv(i,j,4)/float(isumt)
       pv(i,j,13)=3.6e3*pv(i,j,13)/float(isumt)
       pv(i,j,14)=3.6e3*pv(i,j,14)/float(isumt)
       pv(i,j,5)=pv(i,j,1)+pv(i,j,2)+pv(i,j,3)+pv(i,j,4)
       rhoav=rhoav/float(isumt)
c
c---Fix to make sure that flux forcing is only calculated through
c---pbl depth since parcels are now lifted from multiple levels (ST)
c
       if (ght(i,j,kboss).le.pblh(i,j)) then
       pv(i,j,6)=3.6e3*hfxav*float(iibl)/
     &      (cp*float(isumt)*rhoav*pblhav)
       pv(i,j,7)=3.6e3*lv1*qfxav*float(iibl)/
     &      (cp*float(isumt)*rhoav*pblhav)
       pv(i,j,7)=3.6e3*qfxav*float(iibl)/
     &      (cp*float(isumt)*rhoav*pblhav)
       else
       pv(i,j,6)=0.
       pv(i,j,7)=0.
       end if
c
c---end of fix
c
       pv(i,j,8)=10.*float(iibl)/float(isumt)
       pv(i,j,9)=pv(i,j,5)+pv(i,j,6)+pv(i,j,7)
c       pv(i,j,7)=zbmin(i,j)
       pv(i,j,8)=tm
       pv(i,j,12)=tl
c       pv(i,j,17)=kmsemax(i,j)
       pv(i,j,18)=zmsemax(i,j)-ter(i,j)
C
       kk=kmsemax(i,j)
c       dtmk=pv(i,j,1)+pv(i,j,2)+pv(i,j,6)
       dtmk=pv(i,j,1)+pv(i,j,2)
c       dtmk=pv(i,j,1)
c       dtmk=pv(i,j,1)
       if (dtmk.gt.dtmkmax) dtmk=dtmkmax
       if (dtmk.lt.dtmkmin) dtmk=dtmkmin
       dqv=0.
       call tlift(t1,dtmk,zmsemax(i,j),p1,q1,dqv,tla,tma,
     &     i,j,kk,mkzh,tba)
       if ( abs(tla-tl).gt.tdiffmax) then
        idmax=i
        jdmax=j
        tdiffmax=abs(tla-tl)
        tlmax=tl
        tlamax=tla
       endif
       pv(i,j,10)=tla-tl  !difference in lifted parcel temp from changing temp
C
       dtmk=0.
c       dqv=(pv(i,j,3)+pv(i,j,4)+pv(i,j,7))*cp/lv1
       dqv=(pv(i,j,3)+pv(i,j,4))*cp/lv1
c       dqv=pv(i,j,3)*cp/lv1
c       dqv=pv(i,j,4)
       if (dqv.gt.dqvmax) dqv=dqvmax
       if (dqv.lt.dqvmin) dqv=dqvmin
       call tlift(t1,dtmk,zmsemax(i,j),p1,q1,dqv,tlb,tma,
     &    i,j,kk,mkzh,tbb)
       pv(i,j,11)=tlb-tl
C
 931   format(2i5,3f10.2)
c       pv(i,j,5)=pv(i,j,1)+pv(i,j,2)+pv(i,j,6)
       pv(i,j,5)=pv(i,j,1)+pv(i,j,2)
c       pv(i,j,9)=pv(i,j,3)+pv(i,j,4)+pv(i,j,7)
       pv(i,j,9)=pv(i,j,3)+pv(i,j,4)
       pv(i,j,15)=-(pv(i,j,13)+pv(i,j,14))
       pv(i,j,16)=pv(i,j,10)+pv(i,j,11)+pv(i,j,15)
c       pv(i,j,16)=pv(i,j,10)+pv(i,j,11)-pv(i,j,13)
c       pv(i,j,16)=pv(i,j,10)-pv(i,j,13)
       if ((i.eq.300).and.(j.eq.300)) then
        print*,(pv(i,j,nnk),nnk=1,15)
       endif
 271  continue
C
      print*,'got here.'
      write(6,939)idmax,jdmax,tdiffmax,tl1max,tlmax
 939  format(2i5,3f11.4)
      do 252 j=2,mjx-1
      do 252 i=2,miy-1
       do k=1,mkzh
        prd(i,j,k)=0.25*(prs(i-1,j,k)+prs(i,j-1,k)+prs(i-1,j-1,k)+
     &       prs(i,j,k))
       enddo
C Vertical interpolation to pressure coordinates
       klev=1
       do k=2,mkzh
C        print*,i,j,k,prs(i,j,k),klev,prlev(klev)
        if (prs(i,j,k).lt.prlev(klev)) then
         if (k.eq.mkzh) then
C          print*,'below ground, w.',i,j,k,prs(i,j,k),klev,prlev(klev)
          do kk=klev,klevmax
           wp(i,j,kk)=0.
           thp(i,j,kk)=thp(i,j,kk-1)
           vorp(i,j,kk)=vorp(i,j,klev)
           ioutc(i,j,kk)=1
          enddo
          goto 35
         else
          goto 34
         endif
        endif
        if ((prs(i,j,k).ge.prlev(klev)).and.
     &      (prs(i,j,k-1).le.prlev(klev))) then
         wp(i,j,klev)=www(i,j,k-1) + (prlev(klev)-prs(i,j,k-1))*
     &         (www(i,j,k)-www(i,j,k-1))/(prs(i,j,k)-prs(i,j,k-1))
         vorp(i,j,klev)=vor(i,j,k-1) + (prlev(klev)-prs(i,j,k-1))*
     &         (vor(i,j,k)-vor(i,j,k-1))/(prs(i,j,k)-prs(i,j,k-1))
         thp(i,j,klev)=theta(i,j,k-1) + (prlev(klev)-prs(i,j,k-1))*
     &         (theta(i,j,k)-theta(i,j,k-1))/(prs(i,j,k)-prs(i,j,k-1))
         klev=klev + 1
         if (klev.gt.klevmax) goto 35
         if (k.eq.mkzh) then
          do kk=klev,klevmax
           wp(i,j,kk)=0.
           thp(i,j,kk)=thp(i,j,kk-1)
           vorp(i,j,kk)=vorp(i,j,klev)
           ioutc(i,j,kk)=1
          enddo
          goto 35
         endif
        endif
 34     continue
       enddo
 35    klev=1
       do k=2,mkzh
C        print*,i,j,k,prd(i,j,k),klev,prlev(klev)
        if (prd(i,j,k).lt.prlev(klev)) then
         if (k.eq.mkzh) then
C          print*,'below ground, u,v.',i,j,k,prd(i,j,k)
          do kk=klev,klevmax
           up(i,j,kk)=up(i,j,klev-1)
           vp(i,j,kk)=vp(i,j,klev-1)
           ioutd(i,j,kk)=1
          enddo
          goto 38
         else
          goto 37
         endif
        endif
        if ((prd(i,j,k).ge.prlev(klev)).and.
     &      (prd(i,j,k-1).le.prlev(klev))) then
         up(i,j,klev)=uuu(i,j,k-1) + (prlev(klev)-prd(i,j,k-1))*
     &         (uuu(i,j,k)-uuu(i,j,k-1))/(prd(i,j,k)-prd(i,j,k-1))
         vp(i,j,klev)=vvv(i,j,k-1) + (prlev(klev)-prd(i,j,k-1))*
     &         (vvv(i,j,k)-vvv(i,j,k-1))/(prd(i,j,k)-prd(i,j,k-1))
         klev=klev + 1
	 if (klev.gt.klevmax) goto 38
	 if (k.eq.mkzh) then
          do kk=klev,klevmax
           up(i,j,kk)=up(i,j,klev-1)
           vp(i,j,kk)=vp(i,j,klev-1)
           ioutd(i,j,kk)=1
          enddo
          goto 38
         endif
        endif
 37     continue
       enddo
 38    continue
 252  continue
c
C End Vertical Interpolation


C
C       Write to file='pvocalc.out'
C       In do-loops i=imin,imax and j=jmin,jmax
C       define output subdomain with imin, imax, jmin, jmax.
C   
        do j=459,513
         do i=360,414
c        do j=418,472
c         do i=449,503
c        do j=380,430
c          do i=405,430
c        do j=440,490
c          do i=297,322
c        do j=295,320
c          do i=435,485


c            I wasn't sure how many fields or which pv levels to output.
c            These must be defined below. - Ahijevych
c            Also make sure the number of fields matches the number
c             before F12.3 in format. 
c
             write(50,990) kboss,kbmin3(i,j,kboss),j,i,pv(i,j,18),
     &          pv(i,j,10),pv(i,j,11),pv(i,j,13),pv(i,j,14),pv(i,j,12)
          enddo
        enddo
      enddo ! kboss loop
 990  format (4(i3,1x),6(f10.4))
      close(unit=50)
      write(19) kbmin3
      return
      end
