c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine slpbmcalc(qvp,tmk,pstx,prs,sigh,sigf,ter,slp,miy,mjx,
     &   mkzh)
c
c   This routine calculates SLP using a modified form of the
c   Benjamin and Miller (1990) method.  Instead of using 700-hPa
c   temperature, it uses the temperature at a sigma level (the same
c   sigma level for all points).  The sigma level that is chosen is
c   that which would be closest to 850 hPa for an ocean point with
c   surface pressure of 1000 hPa.
c
      dimension qvp(miy,mjx,mkzh),tmk(miy,mjx,mkzh),pstx(miy,mjx),
     &   prs(miy,mjx,mkzh),sigh(mkzh),sigf(mkzh+1),ter(miy,mjx),
     &   slp(miy,mjx)
c
      include 'comconst'
c
      expon=rgas*ussalr/grav
      exponi=1./expon
c
c   Get desired sigma level to use temperature from.
c
      sigc=(850.-ptop)/(1000.-ptop)
      do k=1,mkzh
         if (sigc.ge.sigf(k).and.sigc.le.sigf(k+1)) kupper=k
      enddo
c
c   Create sea-level pressure.
c
      do 190 j=1,mjx-1
      do 190 i=1,miy-1
         psfc=prs(i,j,mkzh)+(1.-sigh(mkzh))*pstx(i,j)
         pupper=prs(i,j,kupper)
         t0=tmk(i,j,kupper)*(psfc/pupper)**expon
         t0v=virtual(t0,qvp(i,j,mkzh))
         slp(i,j)=psfc*(1.+ussalr/t0v*ter(i,j))**exponi
 190  continue
      return
      end
