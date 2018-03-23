c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine slpcalc(qvp,tmk,pstx,prs,sigh,ter,slp,miy,mjx,mkzh)
c
      dimension qvp(miy,mjx,mkzh),tmk(miy,mjx,mkzh),pstx(miy,mjx),
     &   prs(miy,mjx,mkzh),sigh(mkzh),ter(miy,mjx),slp(miy,mjx)
c
      include 'comconst'
c
      expon=rgas*ussalr/grav
      exponi=1./expon
c
c   Create sea-level pressure.
c   Given: at lowest half sigma layer (LHSL): Tv, p
c          at surface: z, p
c          at sea level: z (=0)
c
c   First, get z at LHSL using altimeter equation between LHSL and
c   surface. Then, get p at sea level using alt. eqn. between LHSL
c   and sea level. Assume p. pert. at surface equals p. pert. at LHSL.
c
      do 190 j=1,mjx-1
      do 190 i=1,miy-1
         tvlhsl=virtual(tmk(i,j,mkzh),qvp(i,j,mkzh))
         prslhsl=prs(i,j,mkzh)
         psurf=prs(i,j,mkzh)+(1.-sigh(mkzh))*pstx(i,j)
         ghtlhsl=ter(i,j)+tvlhsl/ussalr*
     &      ((psurf/prslhsl)**expon - 1.)
         slp(i,j)=prslhsl*(1.+ussalr/tvlhsl*ghtlhsl)**exponi
 190  continue
      return
      end
