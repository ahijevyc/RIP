c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine dpbhcalc(qvp,tmk,ght,prs,sigf,sigh,pstx,
     &   ter,dpbh,h1,h2,miy,mjx,mkzh)
c
c   This routine calculates the pressure differencfe between two
c   height levels.
c
      dimension qvp(miy,mjx,mkzh),tmk(miy,mjx,mkzh),ght(miy,mjx,mkzh),
     &   prs(miy,mjx,mkzh),sigf(mkzh+1),sigh(mkzh),ter(miy,mjx),
     &   pstx(miy,mjx),dpbh(miy,mjx)
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
      do ih=1,2
         hh=h1*100.               ! hm to m
         if (ih.eq.2) hh=h2*100.  ! hm to m
c
      do j=1,mjx-1
      do i=1,miy-1
         if (hh.gt.ght(i,j,1)) then
            prshh=rmsg
            goto 35
         endif
         do k=1,mkzh-1
            if (hh.le.ght(i,j,k).and.hh.ge.ght(i,j,k+1)) then
               exphh=exp(-hh/sclht)
               expz=exp(-ght(i,j,k)/sclht)
               expzp1=exp(-ght(i,j,k+1)/sclht)
               prshh=((exphh-expz)*prs(i,j,k+1)+
     &                (expzp1-exphh)*prs(i,j,k))/(expzp1-expz)
               goto 35
            endif
         enddo
c
c      If you've gotten to here, then the level is below the lowest
c      half sigma level.
c
         sfp=prs(i,j,mkzh)+(1.-sigh(mkzh))*pstx(i,j)
         if (hh.ge.ter(i,j)) then
c
c         Level is between lowest half sigma level and ground.
c
            exphh=exp(-hh/sclht)
            expz=exp(-ght(i,j,mkzh)/sclht)
            expzp1=exp(-ter(i,j)/sclht)
            prshh=((exphh-expz)*sfp+
     &          (expzp1-exphh)*prs(i,j,mkzh))/(expzp1-expz)
         else
c
c         Level is below ground.
c
            pupper=prs(i,j,kupper)
            t0=tmk(i,j,kupper)*(sfp/pupper)**expon
            t0v=virtual(t0,qvp(i,j,mkzh))
            prshh=sfp*(1.+ussalr/t0v*(ter(i,j)-hh))**exponi
         endif
 35      continue
         if (ih.eq.1) then
            dpbh(i,j)=prshh
         else
            dpbh(i,j)=dpbh(i,j)-prshh
         endif
      enddo
      enddo
c
      enddo
c
      return
      end
