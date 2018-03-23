c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine derivcp(arr,icdin,prs,xmap,dmap,pstx,
     &   qvp,tmk,sigh,scr2a,scr2b,deriv,icdout,dir,miy,mjx,mkzh)
c
c   Calculates a horizontal derivative (in either the x or y
c      direction) holding p constant.  Both input and output can
c      be on either dot or cross grid.
c
      dimension arr(miy,mjx,mkzh),prs(miy,mjx,mkzh),xmap(miy,mjx),
     &   dmap(miy,mjx),pstx(miy,mjx),qvp(miy,mjx,mkzh),
     &   tmk(miy,mjx,mkzh),sigh(mkzh),scr2a(miy,mjx),scr2b(miy,mjx),
     &   deriv(miy,mjx,mkzh)
      character*1 dir
c
      include 'comconst'
c
c   Big k-loop
c
      do 200 k=1,mkzh
c
c      First, calculate d(sigma)/dp (per hPa), and put into scr2a
c
         do j=1,mjx-1
         do i=1,miy-1
            if (inhyd.eq.1) then
               refprshpa=sigh(k)*pstx(i,j)+ptop !No pert needed
               reftmk=refslt+reflaps*log(refprshpa/1000.)
c
c            Note, the expression below assumes dp/dz=-rho*g (hydrostatic
c            relationship), which is technically inconsistent with the fact
c            that inhyd=1.  We assume the error is negligible for the purpose
c            of transforming to pressure coordinates.
c
               scr2a(i,j)=refprshpa*virtual(tmk(i,j,k),qvp(i,j,k))/
     &            (pstx(i,j)*prs(i,j,k)*reftmk)
            else
               scr2a(i,j)=1./pstx(i,j)
            endif
         enddo
         enddo
c
c      Put it on dot points if output grid is dot point.
c
         if (icdout.eq.0) call xtodot(scr2a,miy,mjx)
c
c      Put the vertical derivative of arr into scr2b.
c
         kp1=min(mkzh,k+1)
         km1=max(1,k-1)
         dsig=sigh(kp1)-sigh(km1)
         do j=1,mjx-icdin
         do i=1,miy-icdin
            scr2b(i,j)=(arr(i,j,kp1)-arr(i,j,km1))/dsig
         enddo
         enddo
c
c      Switch grids if necessary.
c
         if (icdin.eq.1.and.icdout.eq.0) then
            call xtodot(scr2b,miy,mjx)
         elseif (icdin.eq.0.and.icdout.eq.1) then
            do j=1,mjx-1
            do i=1,miy-1
               scr2b(i,j)=.25*
     &            (scr2b(i,j)+scr2b(i+1,j)+scr2b(i,j+1)+scr2b(i+1,j+1))
            enddo
            enddo
         endif
c
c      Combine scr2a and scr2b
c
         do j=1,mjx-icdout
         do i=1,miy-icdout
            scr2a(i,j)=scr2a(i,j)*scr2b(i,j)
         enddo
         enddo
c
c      Get horizontal derivatives of arr and prs.
c
         if (dir.eq.'x') then
            call ddx(arr(1,1,k),icdin,deriv(1,1,k),xmap,dmap,icdout,
     &         miy,mjx)
            call ddx(prs(1,1,k),    1,        scr2b,xmap,dmap,icdout,
     &         miy,mjx)
         else
            call ddy(arr(1,1,k),icdin,deriv(1,1,k),xmap,dmap,icdout,
     &         miy,mjx)
            call ddy(prs(1,1,k),    1,        scr2b,xmap,dmap,icdout,
     &         miy,mjx)
         endif
c
c   Combine terms.
c
         do j=1,mjx-icdout
         do i=1,miy-icdout
            deriv(i,j,k)=deriv(i,j,k)-scr2a(i,j)*scr2b(i,j)
         enddo
         enddo
c
  200 continue
c
      return
      end
