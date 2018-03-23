c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine cttcalc(pstx,sigf,sigh,tmk,qcw,qci,ctt,
     &   miy,mjx,mkzh)
c
c   This routine calculates cloud top brightness temperatures in Cels.
c   It assumes:
c      1) zenith angle is 0
c      2) brightness temperature is roughly the temperature at
c         unity optical depth into the cloud
c      3) cloud absorption coefficient is constant
c
      dimension pstx(miy,mjx),sigf(mkzh+1),sigh(mkzh),
     &   tmk(miy,mjx,mkzh),qcw(miy,mjx,mkzh),ctt(miy,mjx),
     &   qci(miy,mjx,mkzh)
c
      include 'comconst'
c
c   Create cloud-top temperature.
c
      do 190 j=1,mjx-1
      do 190 i=1,miy-1
         opdepthd=0.
         k=0
c
c      Integrate downward from model top, calculating path at full
c      sigma levels.
c
   20    opdepthu=opdepthd
         k=k+1
         dp=100.*(pstx(i,j)*(sigf(k+1)-sigf(k))) !ignoring pres. pert.
         if (iice.eq.0) then
            if (tmk(i,j,k).lt.celkel) then
c             Note: abscoefi is m**2/g, qcw is g/kg,
c                   so no convrsion needed
               opdepthd=opdepthu+abscoefi*qcw(i,j,k)*dp/grav
            else
               opdepthd=opdepthu+abscoef*qcw(i,j,k)*dp/grav
            endif
         else
            opdepthd=opdepthd+(abscoef*qcw(i,j,k)+
     &                        abscoefi*qci(i,j,k))*dp/grav
         endif
         if (opdepthd.lt.1..and.k.lt.mkzh) then
            goto 20
         elseif (opdepthd.lt.1..and.k.eq.mkzh) then
            sigctt=sigh(mkzh)
         else
            fac=(1.-opdepthu)/(opdepthd-opdepthu)
            sigctt=sigf(k)+fac*(sigf(k+1)-sigf(k))
            sigctt=min(sigh(mkzh),max(sigh(1),sigctt))
         endif
         do 30 k=2,mkzh
            if (sigctt.ge.sigh(k-1).and.sigctt.le.sigh(k)) then
               fac=(sigctt-sigh(k-1))/(sigh(k)-sigh(k-1))
               ctt(i,j)=tmk(i,j,k-1)+
     &            fac*(tmk(i,j,k)-tmk(i,j,k-1))-celkel
               goto 40
            endif
   30    continue
   40    continue
 190  continue
      return
      end
