c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine pfcalc(sigh,sigf,pstx,prs,pf,miy,mjx,mkzh)
c
c   Calculates the pressure at full sigma levels. It is
c   assumed that pressure at the model bottom [i.e., the
c   last full sigma level] is pstx + ptop + p-prime(mkzh),
c   or pr(mkzh) + (1-sigh(mkzh))*pstx, and that the
c   pressure at the model top [i.e., the 1st full sigma
c   level] is ptop + p_prime(1), or pr(1)-sigh(1)*pstx.
c   The array only contains mkzh levels, so the pressure at
c   the model top is excluded.  The kth value of pf is the
c   pressure at the (k+1)th full sigma level, which is the
c   lower bounding level for the kth sigma layer.  Thus, pf(mkzh)
c   would be surface pressure.
c
      dimension sigh(mkzh),sigf(mkzh+1),pstx(miy,mjx),
     &   prs(miy,mjx,mkzh),pf(miy,mjx,mkzh)
c
      include 'comconst'
c
      do j=1,mjx-1
      do i=1,miy-1
         pf(i,j,mkzh)=prs(i,j,mkzh)+(1.-sigh(mkzh))*pstx(i,j)
         do k=2,mkzh
            pf(i,j,k-1)=((sigh(k)-sigf(k))*prs(i,j,k-1)+
     &         (sigf(k)-sigh(k-1))*prs(i,j,k))/(sigh(k)-sigh(k-1))
         enddo
      enddo
      enddo
c
      return
      end
