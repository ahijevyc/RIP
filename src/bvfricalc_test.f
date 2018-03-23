c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine bvfricalc(prs,ght,tmk,velscr,qvp,qls,
     &   uvel,vvel,work,cosa,sina,cphase,field,isat,arr,miy,mjx,mkzh)
c
c   Note: qvp is vapor mix. rat., qls is sum of all liquid
c   and solid water mixing ratios (both in kg/kg)
c
      dimension prs(miy,mjx,mkzh), ght(miy,mjx,mkzh),
     &   tmk(miy,mjx,mkzh),qvp(miy,mjx,mkzh),qls(miy,mjx,mkzh),
     &   velscr(miy,mjx,mkzh),uvel(miy,mjx,mkzh),
     &   vvel(miy,mjx,mkzh),work(miy,mjx),arr(miy,mjx,mkzh),
     &   alapse(miy,mjx,mkzh),thte(miy,mjx,mkzh),ci(miy,mjx,mkzh)
      character field*10
c
      dimension ws(100),wt(100),th(100),rh(100)
c
      include 'comconst'
c
      if (field(1:2).ne.'bv') then
         if (field(7:9).ne.'   ') then
c            read(field(7:9),'(i3)') nkmsmth
c            numpas=300+nint(nkmsmth/dskm)
            read(field(7:9),'(i3)') numpas
            write(iup,*)'Smoothing vel. with numpas = ',numpas
         else
            numpas=0
         endif
         do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
              velscr(i,j,k)=.25*(cosa*(uvel(i,j,k)+uvel(i+1,j,k)+
     &                          uvel(i,j+1,k)+uvel(i+1,j+1,k))+
     &                    sina*(vvel(i,j,k)+vvel(i+1,j,k)+
     &                          vvel(i,j+1,k)+vvel(i+1,j+1,k)))
            enddo
            enddo
c
c         Smooth this level
c
            call smooth(velscr(1,1,k),work,numpas,miy,miy-1,mjx-1)
         enddo
      endif
c
c---calculate needed arrays for bvfsq array based on eq.(6) 
c---of Kirshbaum and Durran (2004)
c
c---theta-e, specific heat of ice, moist adiabatic lapse rate
c
      ethc=2825.42
      do 1000 k = 1, mkzh
      do 1000 j = 1, mjx-1
      do 1000 i = 1, miy-1
        kp1=min(k+1,mkzh)
        km1=max(k-1,1)
        q=max(qvp(i,j,k),1.e-15)
        t=tmk(i,j,k)
        p=prs(i,j,k)
        e=q*p/(eps+q)
        tlcl=tlclc1/(log(t**tlclc2/e)-tlclc3)+tlclc4
        thte(i,j,k)=t*(1000./p)**(gamma*(1.+gammamd*q))*
     &              exp((2825.42*q)/t)
        ci(i,j,k)=2106.+7.3*(tmk(i,j,k)-273.15)
        atmkp1=tonpsadiabat(thte(i,j,k),prs(i,j,kp1))
        atmkm1=tonpsadiabat(thte(i,j,k),prs(i,j,km1))
        alapse(i,j,k)=(atmkp1-atmkm1)/(ght(i,j,kp1)-ght(i,j,km1))
 1000 continue
      do 200 j=1,mjx-1
      do 200 i=1,miy-1
c
c      
      do 150 k=1,mkzh
         t = tmk(i,j,k)
         gammam=gamma*(1.+gammamd*qvp(i,j,k))
         th(k)=t*(1000./prs(i,j,k))**gammam
c
c   es should be wrt ice if state=ice and temp. is below freezing.
c
         if (field(6:6).eq.'i'.and.t.lt.celkel) then
            es = ezero * exp( esicon1-esicon2/t )
         else
            es = ezero * exp( eslcon1*(t-celkel)/(t-eslcon2) )
         endif
         w = qvp(i,j,k)
         wsat = eps*es/(prs(i,j,k)-es)
         rh(k)=100.*w/wsat
         ws(k) = wsat
         wt(k) = wsat + qls(i,j,k)
  150 continue
c
      do 160 k=1,mkzh
         kp1=min(k+1,mkzh)
         km1=max(k-1,1)
         kk=max(min(k,mkzh-1),2)
         if ((field(6:6).eq.'l'.or.field(6:6).eq.'i').and.
     &        (rh(k).gt.90..or.isat.eq.1)) then
c
c----Calculation of equation (6) from Kirshbaum and Durran (2004)
c
         coeff=1./(1.+wt(k))
         deltaz=ght(i,j,kp1)-ght(i,j,km1)
         term1=alapse(i,j,k)*((cp/thte(i,j,k)*(thte(i,j,kp1)-
     &         thte(i,j,km1))/deltaz)+2./deltaz*((qls(i,j,k)*ci(i,j,k)/
     &         thte(i,j,k))*(thte(i,j,kp1)-thte(i,j,km1))+
     &         (ci(i,j,k)*log(thte(i,j,k))*(qls(i,j,kp1)-qls(i,j,km1))
     &         +(qls(i,j,k)*log(thte(i,j,k))*(ci(i,j,kp1)-ci(i,j,km1)
     &          )))))
         term2=(ci(i,j,k)*alapse(i,j,k)*log(tmk(i,j,k))+grav)*
     &         (qls(i,j,kp1)-qls(i,j,km1))/deltaz
         bvfsq=coeff*(term1-term2)
c
c         The following is from Durran and Klemp (1982)
c
c            xlhc=xlhc0-xlhctd*tmk(i,j,k)
c            bvfsq=grav*((1.+xlhc*ws(k)/(rgas*tmk(i,j,k)))/
c     &         (1.+eps*xlhc**2*ws(k)/(cp*rgas*tmk(i,j,k)**2))*
c     &         ((th(kp1)-th(km1))/(ght(i,j,kp1)-ght(i,j,km1))/th(k)+
c     &         xlhc/(cp*tmk(i,j,k))*(ws(kp1)-ws(km1))/
c     &         (ght(i,j,kp1)-ght(i,j,km1)))-(wt(kp1)-wt(km1))/
c     &         (ght(i,j,kp1)-ght(i,j,km1)))
         else
            bvfsq=grav*(th(kp1)-th(km1))/(ght(i,j,kp1)-ght(i,j,km1))/
     &         th(k)
         endif
         if (field(1:5).eq.'bvfsq') then
            arr(i,j,k)=bvfsq  ! per second squared
         elseif (field(1:5).eq.'richn') then
            dudz=(velscr(i,j,kp1)-velscr(i,j,km1))/
     &         (ght(i,j,kp1)-ght(i,j,km1))
            arr(i,j,k)=bvfsq/(dudz*dudz)  ! dim'less
c
c         Use a logarithimc function for large Ri
c
            if (arr(i,j,k).gt.1.0) arr(i,j,k)=1.0+log(arr(i,j,k))
         elseif (field(1:4).eq.'spsq') then
            arr(i,j,k)=0.0
            umc=velscr(i,j,k)-cphase
            if (abs(umc).lt..001) umc=.001
            umci=1./umc
            if (field(5:5).eq.'a'.or.field(5:5).eq.'t')
     &         arr(i,j,k)=arr(i,j,k)+bvfsq*umci*umci
            if (field(5:5).eq.'b'.or.field(5:5).eq.'t') then
               d2udz2=((velscr(i,j,kk+1)-velscr(i,j,kk))/
     &              (ght(i,j,kk+1)-ght(i,j,kk))-
     &              (velscr(i,j,kk)-velscr(i,j,kk-1))/
     &              (ght(i,j,kk)-ght(i,j,kk-1)))/
     &             (.5*(ght(i,j,kk+1)-ght(i,j,kk-1)))
               arr(i,j,k)=arr(i,j,k)-d2udz2*umci
            endif
            arr(i,j,k)=1.e6*arr(i,j,k)  ! per km squared
         endif
  160 continue
c
  200 continue
      return
      end
