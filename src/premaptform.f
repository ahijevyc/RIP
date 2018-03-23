c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine premaptform(dskmc,miycors,mjxcors,nproj,
     &   xlatc,xlonc,true1,true2,iup)
c
c     This routine calculates constants for routine maptform and puts
c     them in common block mptf, which should be declared in (and only
c     in) the main program and routines premaptform (this routine),
c     maptform, and mapfaccalc.
c
      common /mptf/ rpd_mptf,pi_mptf,dskmc_mptf,xlonc_mptf,rearth_mptf,
     &   ciy_mptf,cjx_mptf,yc_mptf,ihm_mptf,c1_mptf,c2_mptf,cone_mptf,
     &   conei_mptf,nproj_mptf
c
      pi_mptf=4.*atan(1.)     ! 3.1415...
      rpd_mptf=pi_mptf/180.    ! radians per degree
      rearth_mptf=6370.949  ! radius of planet, in km
      dskmc_mptf=dskmc
      xlonc_mptf=xlonc
      nproj_mptf=nproj
      ciy_mptf=.5*(1.+miycors)
      cjx_mptf=.5*(1.+mjxcors)
c
      if (nproj_mptf.eq.1) then   ! Lambert Comformal
c
c   Make sure xlatc, true1, and true2 are all in same hemisphere,
c      and calculate ihm_mptf.
c
      if (xlatc.gt.0..and.true1.gt.0..and.true2.gt.0.) then
         ihm_mptf=1
      elseif (xlatc.lt.0..and.true1.lt.0..and.true2.lt.0.) then
         ihm_mptf=-1
      else
         write(iup,*)'Invalid latitude parameters for LC map.'
         stop
      endif
c
c   Calculate cone factor
c
      if (true1.ne.true2) then
         cone_mptf=log10(cos(rpd_mptf*true1)/cos(rpd_mptf*true2))/
     &        log10(tan(.25*pi_mptf-ihm_mptf*.5*rpd_mptf*true1)/
     &              tan(.25*pi_mptf-ihm_mptf*.5*rpd_mptf*true2))
      else
         cone_mptf=cos(rpd_mptf*(90.-ihm_mptf*true1))
      endif
c
c   Calculate other constants
c
      conei_mptf=1./cone_mptf
      cotrue1=ihm_mptf*90.-true1
      c1_mptf=rearth_mptf*sin(rpd_mptf*cotrue1)/
     &   (cone_mptf*(ihm_mptf*tan(.5*rpd_mptf*cotrue1))**cone_mptf)
      c2_mptf=tan(.5*rpd_mptf*cotrue1)*(cone_mptf/
     &   (ihm_mptf*rearth_mptf*sin(rpd_mptf*cotrue1)))**conei_mptf
      yc_mptf=-c1_mptf*(ihm_mptf*tan(.25*(ihm_mptf*pi_mptf-
     &   2.*rpd_mptf*xlatc)))**cone_mptf
c
      elseif (nproj_mptf.eq.2) then   ! Polar Stereographic
c
c   Make sure xlatc and true1 are both in same hemisphere,
c      and calculate ihm_mptf.
c
      if (xlatc.gt.0..and.true1.gt.0.) then
         ihm_mptf=1
      elseif (xlatc.lt.0..and.true1.lt.0.) then
         ihm_mptf=-1
      else
         write(iup,*)'Invalid latitude parameters for PS map.'
         stop
      endif
c
c   Calculate cone factor
c
      cone_mptf=1.
c
c   Calculate other constants
c
      conei_mptf=1./cone_mptf
      cotrue1=ihm_mptf*90.-true1
      c1_mptf=1.+cos(rpd_mptf*cotrue1)
      c2_mptf=1. !not used
      yc_mptf=-rearth_mptf*sin(.5*ihm_mptf*pi_mptf-rpd_mptf*xlatc)*
     &   c1_mptf/(1.+cos(.5*ihm_mptf*pi_mptf-rpd_mptf*xlatc))
c
      elseif (nproj_mptf.eq.0.or.nproj_mptf.eq.3) then   ! Mercator
c
c       (Note: nproj=0 means "idealized" (i.e. no map), but we'll treat
c        it as Mercator so nothing in RIP goes haywire.)
c
      ihm_mptf=1
      cone_mptf=1.
      conei_mptf=1.
      c1_mptf=cos(rpd_mptf*true1)
      c2_mptf=1.  !not used
      yc_mptf=rearth_mptf*log((1.+sin(rpd_mptf*xlatc))/
     &   cos(rpd_mptf*xlatc))
c
      endif
c
      return
      end
