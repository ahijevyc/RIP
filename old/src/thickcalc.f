c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine thickcalc(qvp,tmk,ght,prs,sigh,sigf,ter,pstx,
     &   thick,p1,p2,miy,mjx,mkzh)
c
      dimension qvp(miy,mjx,mkzh),tmk(miy,mjx,mkzh),ght(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),sigf(mkzh+1),ter(miy,mjx),pstx(miy,mjx),
     &   thick(miy,mjx),sigh(mkzh)
c
      include 'comconst'
c
      expon=rgas*ussalr/grav
      exponi=1./expon
c
c   Get desired sigma level on which to base subterranean
c   temperature profile.
c
      sigc=(850.-ptop)/(1000.-ptop)
      do k=1,mkzh
         if (sigc.ge.sigf(k).and.sigc.le.sigf(k+1)) kupper=k
      enddo
c
c   Calculate dp
c
      do ip=1,2
         pr=p1*10.               ! hm to m
         if (ip.eq.2) pr=p2*10.  ! kPa to hPa
c
      do j=1,mjx-1
      do i=1,miy-1
         if (pr.lt.prs(i,j,1)) then
            ghtpr=rmsg
            goto 35
         endif
         do k=1,mkzh-1
            if (pr.ge.prs(i,j,k).and.pr.le.prs(i,j,k+1)) then
               expz=exp(-ght(i,j,k)/sclht)
               expzp1=exp(-ght(i,j,k+1)/sclht)
               expzpr=((pr-prs(i,j,k))*expzp1+(prs(i,j,k+1)-pr)*expz)/
     &                (prs(i,j,k+1)-prs(i,j,k))
               ghtpr=-log(expzpr)*sclht
               goto 35
            endif
         enddo
c
c      If you've gotten to here, then the level is below the lowest
c      half sigma level.
c
         psfc=prs(i,j,mkzh)+(1.-sigh(mkzh))*pstx(i,j)
         if (pr.le.psfc) then
            expz=exp(-ght(i,j,mkzh)/sclht)
            expzp1=exp(-ter(i,j)/sclht)
            expzpr=((pr-prs(i,j,mkzh))*expzp1+(psfc-pr)*expz)/
     &             (psfc-prs(i,j,mkzh))
            ghtpr=-log(expzpr)*sclht
         else
c
c         If you've gotten to here, then the level is below ground.
c
            tsbg=tmk(i,j,kupper)*(psfc/prs(i,j,kupper))**expon
            tsbgv=virtual(tsbg,qvp(i,j,mkzh))
            ghtpr=ter(i,j)+tsbgv/ussalr*(1.-(pr/psfc)**expon)
         endif
 35      continue
         if (ip.eq.1) then
            thick(i,j)=ghtpr
         else
            thick(i,j)=.1*abs(thick(i,j)-ghtpr)  ! m to dam
         endif
      enddo
      enddo
c
      enddo
c
      return
      end
