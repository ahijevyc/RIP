c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine dpcalc(sigh,sigf,pstx,prs,dp,miy,mjx,mkzh)
c
c   Calculates the pressure thickness of sigma layers
c
      dimension sigh(mkzh),sigf(mkzh+1),pstx(miy,mjx),
     &   prs(miy,mjx,mkzh),dp(miy,mjx,mkzh),pfull(200)
c
      include 'comconst'
c
      do j=1,mjx-1
      do i=1,miy-1
         pfull(1)=prs(i,j,1)-sigh(1)*pstx(i,j)
         pfull(mkzh+1)=prs(i,j,mkzh)+(1.-sigh(mkzh))*pstx(i,j)
         do k=2,mkzh
            pfull(k)=((sigh(k)-sigf(k))*prs(i,j,k-1)+
     &         (sigf(k)-sigh(k-1))*prs(i,j,k))/(sigh(k)-sigh(k-1))
         enddo
         do k=1,mkzh
            dp(i,j,k)=pfull(k+1)-pfull(k)
         enddo
      enddo
      enddo
c
      return
      end
