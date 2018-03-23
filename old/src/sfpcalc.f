c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine sfpcalc(pstx,sigh,prs,sfp,miy,mjx,mkzh)
c
      dimension pstx(miy,mjx),prs(miy,mjx,mkzh),sfp(miy,mjx),
     &   sigh(mkzh)
c
      include 'comconst'
c
c   Create surface pressure. Assume pres. pert. at surface equals
c   pres. pert. at LHSL.
c
      do 190 j=1,mjx-1
      do 190 i=1,miy-1
c         sfp(i,j)=100.*sin(2.*pi*float(j)/40.)*
c     &      (.5+.5*sin(2.*pi*float(i)/340.))
         sfp(i,j)=prs(i,j,mkzh)+(1.-sigh(mkzh))*pstx(i,j)
 190  continue
      return
      end
