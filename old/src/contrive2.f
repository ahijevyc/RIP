c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine contrive2(cor,dmap,xmap,ter,pstd,pstx,prs,ght,
     &   tmk,qvp,uuu,vvv,www,sigh,rootname,iendcr,miy,mjx,mkzh)
c
c   This routine generates artificial model output fields for testing,
c   etc.
c
      parameter(mkzslab=250,ztopslab=22.,vagfac=4.)
c
      dimension cor(miy,mjx),dmap(miy,mjx),xmap(miy,mjx),ter(miy,mjx),
     &   pstd(miy,mjx),pstx(miy,mjx),prs(miy,mjx,mkzh),
     &   qvp(miy,mjx,mkzh),
     &   ght(miy,mjx,mkzh),tmk(miy,mjx,mkzh),uuu(miy,mjx,mkzh),
     &   vvv(miy,mjx,mkzh),www(miy,mjx,mkzh),sigh(mkzh)
      dimension refhp(miy,mjx),refhu(miy,mjx),refhv(miy,mjx),
     &   refhth(miy,mjx)
      dimension slabth(miy,mkzslab)
      character rootname*256
      real m1,m2,m3,m4,m5,m6,m7,m8
c
      include 'comconst'
c
      open (unit=49,file=rootname(1:iendcr)//'.ctin',
     &   form='formatted',status='old')
      read(49,*)zjet
      read(49,*)slthetacold
      read(49,*)dthetadzt
      read(49,*)dthetadzs
      read(49,*)dthetadyf
      read(49,*)dthetadx
      read(49,*)frontslope
      read(49,*)frontwidth
      read(49,*)sttpshgt
      read(49,*)sttpsslope
      read(49,*)strattranslope
      read(49,*)iycenter
      read(49,*)nsmooth
      read(49,*)def
      read(49,*)refh
      read(49,*)refhpmid
c
c   First set up terrain, xmap, dmap, cor
c
      do j=1,mjx-1
      do i=1,miy-1
c         if (i.ge.33.and.i.le.58) then
c            ter(i,j)=800.*(sin(pi*(i-33.)/25.))**2
c         else
            ter(i,j)=0.
c         endif
c         if (i.ge.32.and.i.le.45.and.j.ge.34.and.j.le.59) then
c            ter(i,j)=max(ter(i,j),800.*(sin(pi*(j-34.)/25.))**2)
c         endif
c         if (i.ge.20.and.i.le.32.and.j.ge.34.and.j.le.59) then
c            ter(i,j)=800.*(sin(pi*(j-34.)/25.))**2*
c     &                    (sin(pi*(i-19.5)/25.))**2
c         endif
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
c   Next define ght based on reference state.
c
      do k=1,mkzh
      do j = 1, mjx-1
      do i = 1, miy-1
         alnpref=log(sigh(k)*pstx(i,j)+ptop)
         ght(i,j,k)=cc1*alnpref*alnpref+cc2*alnpref+cc3
      enddo
      enddo
      enddo
c
c   Next, calculate theta.  First, determine line equations.  Assume
c   y-axis points toward colder air, and x-axis points in direction
c   of main tropopause-level jet.
c
c   1: southern frontal boundary
c
      m1=frontslope
      b1=0.
      dthetadzf=dthetadzt-dthetadyf/frontslope
      if (dthetadzf.ge.dthetadzs) then
         write(iup,*)'dthetadzf,dthetadzs=',dthetadzf,dthetadzs
         write(iup,*)'Increase strat. stability.'
         stop
      endif
c
c   2: northern frontal boundary
c
      m2=frontslope
      b2=-frontslope*frontwidth
c
c   3: stratosphere/front transition
c
      m3=dthetadyf/(dthetadzs-dthetadzf)
      b3=0.
c
c   4: polar tropopause (flat)
c
      m4=0.
      b4=m3*b2/(m3-m2)
c
c   6: sloped subtropical tropopause
c
      m6=sttpsslope
      b6=0.
c
c   7: flat subtropical tropopause
c
      m7=0.
      b7=sttpshgt-zjet
c
c   5: northern stratospheric transition
c
      m5=strattranslope
      dthetadysts=(dthetadzt-dthetadzs)/
     &   (1./m6-1./m5)
      dthetadzsts=dthetadzt-dthetadysts/m6
      b5=0.
c
c   8: southern stratospheric transition
c
      m8=strattranslope
      b8=b7-m8*(b7-b6)/m6
c
      zgrel=-zjet
      yya=(b4-b2)/m2
      yyb=0.
      yyc=(b7-b6)/m6
      yyd=(zgrel-b2)/m2
      yye=(zgrel-b1)/m1
c
c   Initialize the field
c
      dy=dskm
      dx=dskm
      rimidx=.5*miy
      rjmidx=.5*mjx
      imidx=nint(rimidx)
      jmidx=nint(rjmidx)
      rimidd=.5*(miy+1)
      rjmidd=.5*(mjx+1)
      imidd=nint(rimidd)
      jmidd=nint(rjmidd)
      zbottom=m2*min(yyc-yyd-20.,0.)
      if (zbottom.le.-1000.) then
         write(iup,*)'Warning: zbottom = ',zbottom,' km.'
         write(iup,*)'This is very low, and may cause inaccuracies.'
      endif
      zbottomrel=zbottom-zjet
      yyd=(zbottomrel-b2)/m2
      yye=(zbottomrel-b1)/m1
      dz=.005*dy
      mkzslabuse=1+nint(ztopslab/dz)
      if (mkzslabuse.gt.mkzslab) then
         write(iup,*)'Increase mkzslab to at least ',mkzslabuse
         stop
      endif
      ztopslabuse=(mkzslabuse-1)*dz
c
c   Calculate theta
c
      bottomthetacold=slthetacold+zbottom*dthetadzt
      do i=1,miy-1
         y=(i-iycenter)*dy
         if (y.ge.yyd) then
            bottomtheta=bottomthetacold
         elseif (y.le.yyd.and.y.ge.yye) then
            bottomtheta=bottomthetacold+dthetadyf*(y-yyd)
         else
            bottomtheta=bottomthetacold+dthetadyf*(yye-yyd)
         endif
         z1=m1*y+b1
         z2=m2*y+b2
         z4=m4*y+b4
         z3=m3*y+b3
         z5=m5*y+b5
         z6=m6*y+b6
         z7=m7*y+b7
         z8=m8*y+b8
         do k=1,mkzslabuse
            z=zgrel+(k-1)*dz
            if (y.ge.yya) then
               if (z.le.z4) then
                  slabth(i,k)=bottomtheta+
     &               (z-zbottomrel)*dthetadzt
               elseif (z.le.z5) then
                  slabth(i,k)=bottomtheta+
     &               (z4-zbottomrel)*dthetadzt+
     &               (z-z4)*dthetadzs
               elseif (z.le.z8) then
                  slabth(i,k)=bottomtheta+
     &               (z4-zbottomrel)*dthetadzt+
     &               (z5-z4)*dthetadzs+
     &               (z-z5)*dthetadzsts
               else
                  slabth(i,k)=bottomtheta+
     &               (z4-zbottomrel)*dthetadzt+
     &               (z5-z4)*dthetadzs+
     &               (z8-z5)*dthetadzsts+
     &               (z-z8)*dthetadzs
               endif
            elseif (y.ge.yyb) then
               if (z.le.z2) then
                  slabth(i,k)=bottomtheta+
     &               (z-zbottomrel)*dthetadzt
               elseif (z.le.z3) then
                  slabth(i,k)=bottomtheta+
     &               (z2-zbottomrel)*dthetadzt+
     &               (z-z2)*dthetadzf
               elseif (z.le.z5) then
                  slabth(i,k)=bottomtheta+
     &               (z2-zbottomrel)*dthetadzt+
     &               (z3-z2)*dthetadzf+
     &               (z-z3)*dthetadzs
               elseif (z.le.z8) then
                  slabth(i,k)=bottomtheta+
     &               (z2-zbottomrel)*dthetadzt+
     &               (z3-z2)*dthetadzf+
     &               (z5-z3)*dthetadzs+
     &               (z-z5)*dthetadzsts
               else
                  slabth(i,k)=bottomtheta+
     &               (z2-zbottomrel)*dthetadzt+
     &               (z3-z2)*dthetadzf+
     &               (z5-z3)*dthetadzs+
     &               (z8-z5)*dthetadzsts+
     &               (z-z8)*dthetadzs
               endif
            elseif (y.ge.yyc) then
               if (z.le.z2) then
                  slabth(i,k)=bottomtheta+
     &               (z-zbottomrel)*dthetadzt
               elseif (z.le.z1) then
                  slabth(i,k)=bottomtheta+
     &               (z2-zbottomrel)*dthetadzt+
     &               (z-z2)*dthetadzf
               elseif (z.le.z6) then
                  slabth(i,k)=bottomtheta+
     &               (z2-zbottomrel)*dthetadzt+
     &               (z1-z2)*dthetadzf+
     &               (z-z1)*dthetadzt
               elseif (z.le.z8) then
                  slabth(i,k)=bottomtheta+
     &               (z2-zbottomrel)*dthetadzt+
     &               (z1-z2)*dthetadzf+
     &               (z6-z1)*dthetadzt+
     &               (z-z6)*dthetadzsts
               else
                  slabth(i,k)=bottomtheta+
     &               (z2-zbottomrel)*dthetadzt+
     &               (z1-z2)*dthetadzf+
     &               (z6-z1)*dthetadzt+
     &               (z8-z6)*dthetadzsts+
     &               (z-z8)*dthetadzs
               endif
            elseif (y.ge.yyd) then
               if (z.le.z2) then
                  slabth(i,k)=bottomtheta+
     &               (z-zbottomrel)*dthetadzt
               elseif (z.le.z1) then
                  slabth(i,k)=bottomtheta+
     &               (z2-zbottomrel)*dthetadzt+
     &               (z-z2)*dthetadzf
               elseif (z.le.z7) then
                  slabth(i,k)=bottomtheta+
     &               (z2-zbottomrel)*dthetadzt+
     &               (z1-z2)*dthetadzf+
     &               (z-z1)*dthetadzt
               else
                  slabth(i,k)=bottomtheta+
     &               (z2-zbottomrel)*dthetadzt+
     &               (z1-z2)*dthetadzf+
     &               (z7-z1)*dthetadzt+
     &               (z-z7)*dthetadzs
               endif
            elseif (y.ge.yye) then
               if (z.le.z1) then
                  slabth(i,k)=bottomtheta+
     &               (z-zbottomrel)*dthetadzf
               elseif (z.le.z7) then
                  slabth(i,k)=bottomtheta+
     &               (z1-zbottomrel)*dthetadzf+
     &               (z-z1)*dthetadzt
               else
                  slabth(i,k)=bottomtheta+
     &               (z1-zbottomrel)*dthetadzf+
     &               (z7-z1)*dthetadzt+
     &               (z-z7)*dthetadzs
               endif
            else
               if (z.le.z7) then
                  slabth(i,k)=bottomtheta+
     &               (z-zbottomrel)*dthetadzt
               else
                  slabth(i,k)=bottomtheta+
     &               (z7-zbottomrel)*dthetadzt+
     &               (z-z7)*dthetadzs
               endif
            endif
         enddo
      enddo
c
      call smooth(slabth,nsmooth,miy,miy-1,mkzslabuse,rmsg)
c
c   Interpolate theta to model levels, and include x-gradient of theta
c
      do k=1,mkzh
      do j=1,mjx-1
         dtheta=(j-jmidx)*dx*dthetadx
      do i=1,miy-1
         zmod=ght(i,j,k)/1000.
         if (zmod.le.0.) then
            tmk(i,j,k)=slabth(i,1)+dtheta
         elseif (zmod.ge.ztopslabuse) then
            tmk(i,j,k)=slabth(i,mkzslabuse)+dtheta
         else
            do kslab=1,mkzslabuse-1
               z=(kslab-1)*dz
               z1=z+dz
               if (zmod.ge.z.and.zmod.le.z1) then
                  tmk(i,j,k)=((zmod-z)*slabth(i,kslab+1)+
     &                        (z1-zmod)*slabth(i,kslab))/dz+dtheta
                  goto 43
               endif
            enddo
            write(iup,*)'couldn''t interpolate contrived theta.'
            write(iup,*)'i,j,k=',i,j,k
            stop
 43         continue
         endif
      enddo
      enddo
      enddo
c
c   Calculate theta at reference height
c
      refhold=refh
      refhpmidold=refhpmid
      refh=dz*nint(refh/dz)
      refhpmid=refhpmidold*
     &   exp(-grav*1000.*(refh-refhold)/(rgas*celkel))
      write(iup,*)'refh,refhpmid adjusted to ',refh,refhpmid
      refhm=1000.*refh
      kslabrefh=1+nint(refh/dz)
      do j=1,mjx-1
         dtheta=(j-jmidx)*dx*dthetadx
         do i=1,miy-1
            refhth(i,j)=slabth(i,kslabrefh)+dtheta
         enddo
      enddo
c
c   Next, make wind at reference height level, which is just
c   a geostrophic deformation pattern.
c
      do j = 1, mjx
         x=(j-rjmidd)*ds
      do i = 1, miy
         y=(i-rimidd)*ds
         refhu(i,j)=.5*def*x
         refhv(i,j)=-.5*def*y
      enddo
      enddo
c
c   Next, make pressure at reference height
c
      do j = 1, mjx-1
      do i = 1, miy-1
         refhp(i,j)=refhpmid ! first estimate of refhp
      enddo
      enddo
      do iter=1,5 ! iterative refinement of refhp
         do j=2,mjx-1
            avgcor=.5*(cor(1,j-1)+cor(1,j))
            avgrefhv=.5*(refhv(1,j)+refhv(2,j))
            avgrefhp=.5*(refhp(1,j-1)+refhp(1,j))
            avgrefht=.5*(refhth(1,j-1)+refhth(1,j))*
     &         (1000./avgrefhp)**gamma
            refhp(1,j)=refhp(1,j-1)*
     &         exp(avgcor*avgrefhv/(rgas*avgrefht)*ds)
         enddo
         do j=1,mjx-1
         do i=2,miy-1
            avgcor=.5*(cor(i-1,j)+cor(i,j))
            avgrefhu=.5*(refhu(i,j)+refhu(i,j+1))
            avgrefhp=.5*(refhp(i-1,j)+refhp(i,j))
            avgrefht=.5*(refhth(i-1,j)+refhth(i,j))*
     &         (1000./avgrefhp)**gamma
            refhp(i,j)=refhp(i-1,j)*
     &         exp(-avgcor*avgrefhu/(rgas*avgrefht)*ds)
         enddo
         enddo
         refhperrormid=refhp(imidx,jmidx)-refhpmid
         do j=1,mjx-1
         do i=1,miy-1
            refhp(i,j)=refhp(i,j)-refhperrormid
         enddo
         enddo
      enddo
c
c   Calculate pressure at all levels
c
      const=grav*gamma*1000.**gamma/rgas
      do j = 1, mjx-1
      do i = 1, miy-1
c
c   First get two sigma levels surrounding reference level
c
      do k=1,mkzh-1
         if (refhm.le.ght(i,j,k).and.refhm.ge.ght(i,j,k+1)) then
            kabove=k
            kbelow=k+1
            goto 344
         endif
      enddo
      write(iup,*)'couldn''t find two sigma levels surrounding'
      write(iup,*)'reference level, i,j=',i,j
      stop
 344  continue
c
c   Get p at levels above
c
      prs(i,j,kabove)=(refhp(i,j)**gamma-2.*const*
     &   (ght(i,j,kabove)-refhm)/
     &   (tmk(i,j,kabove)+refhth(i,j)))**(1./gamma)
      do k=kabove-1,1,-1
         prs(i,j,k)=(prs(i,j,k+1)**gamma-2.*const*
     &      (ght(i,j,k)-ght(i,j,k+1))/
     &      (tmk(i,j,k)+tmk(i,j,k+1)))**(1./gamma)
      enddo
c
c   Get p at levels below
c      
      prs(i,j,kbelow)=(refhp(i,j)**gamma-2.*const*
     &   (ght(i,j,kbelow)-refhm)/
     &   (tmk(i,j,kbelow)+refhth(i,j)))**(1./gamma)
      do k=kbelow+1,mkzh
         prs(i,j,k)=(prs(i,j,k-1)**gamma-2.*const*
     &      (ght(i,j,k)-ght(i,j,k-1))/
     &      (tmk(i,j,k)+tmk(i,j,k-1)))**(1./gamma)
      enddo
c
      enddo
      enddo
c
c   Convert theta to temperature, set www and qvp to 0.
c
      do k=1,mkzh
      do j = 1, mjx-1
      do i = 1, miy-1
         tmk(i,j,k)=tmk(i,j,k)*(prs(i,j,k)/1000.)**gamma
         www(i,j,k)=0.
         qvp(i,j,k)=0.
      enddo
      enddo
      enddo
c
c   Calculate geostrophic winds aloft.  Use refhu and refhv for scratch
c   space, since they're not needed anymore.
c
      call derivcz(prs,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &   refhu,refhv,uuu,1,'y',miy,mjx,mkzh)
      call derivcz(prs,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &   refhu,refhv,vvv,1,'x',miy,mjx,mkzh)
      fbar=cor(miy/2,mjx/2)
      do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            tv=virtual(tmk(i,j,k),qvp(i,j,k))
            uuu(i,j,k)=-rgas*tv*uuu(i,j,k)/(prs(i,j,k)*fbar)
            vag=vagfac*(-(prs(i,j,k)-600.)/500.)
            vvv(i,j,k)= rgas*tv*vvv(i,j,k)/(prs(i,j,k)*fbar)+vag
         enddo
         enddo
         call xtodot(uuu(1,1,k),miy,mjx)
         call xtodot(vvv(1,1,k),miy,mjx)
      enddo
c
      return
      end
