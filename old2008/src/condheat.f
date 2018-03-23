c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine condheat(tmk,qvp,omg,ifreeze,prs,miy,mjx,mkzh,dthdt)
c
c      *** uses temperature and vertical velocity to compute
c      ***    d/dt (theta) due to condensation from lifting,
c      ***    in kelvins per day
c
      dimension prs(miy,mjx,mkzh),tmk(miy,mjx,mkzh),qvp(miy,mjx,mkzh),
     &   omg(miy,mjx,mkzh),dthdt(miy,mjx,mkzh)
c
      include 'comconst'
c
      rvap=rgas/eps
c
      do k=1,mkzh
      do j=1,mjx-1
      do i=1,miy-1
         q = qvp(i,j,k)
         t = tmk(i,j,k)
         p = prs(i,j,k)
         e = q*p/(eps+q)
c
c   es should be wrt ice if temperature is below freezing,
c      and latent heat should include that of fusion
c
         xlhc=xlhc0-xlhctd*t
         if (t.lt.celkel.and.ifreeze.eq.1) then
            es=ezero*exp(esicon1-esicon2/t)
            xlatent=xlhc+xlhf
         else
            es=ezero*exp(eslcon1*(t-celkel)/(t-eslcon2))
            xlatent=xlhc
         endif
         rhu=100.*(e*(p-es))/(es*(p-e))
         if (omg(i,j,k).lt.0..and.rhu.gt.99.) then
            ws=eps*es/(p-es)
c
c         The following pseudoadiabatic lapse rate can be found in several
c         text books (Bluestein, Cotton and Anthes, Rogers and Yau, etc.).
c         It is probably not exactly consistent with Bolton's
c         relationships for theta_e and e_s, which are used more consistently
c         throughout the rip code.
c
            rmalr=(grav/cp)*(1.+xlatent*ws/(rgas*tmk(i,j,k)))/
     &         (1.+xlatent*xlatent*eps*ws/
     &            (rgas*cp*tmk(i,j,k)*tmk(i,j,k)))
            omghPaps=omg(i,j,k)*.001   ! dPa/s to hPa/s
            gammam=gamma*(1.+gammamd*q)
            cpm=cp*(1.+cpmd*q)
            rgasm=rgas*(1.+rgasmd*q)
            theta=tmk(i,j,k)*(1000./p)**gammam
            dthdt(i,j,k)=(rgasm*theta*omghPaps)/(grav*p)*
     &         (rmalr-grav/cpm)*3600.  ! K/s to K/hr
         else
            dthdt(i,j,k)=0.
         endif
      enddo
      enddo
      enddo
      return
      end
