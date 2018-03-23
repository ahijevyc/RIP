c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine maptform(riy,rjx,rlat,rlon,idir,rrota)
c
c   This routine converts a coarse domain dot grid point, <riy,rjx>,
c   into a lat/lon point <rlat,rlon> if idir=1, or vice versa if
c   idir=-1. It works for Lambert Conformal (LC,1),
c   Polar Stereographic (ST,2), or Mercator (ME,3) projections,
c   with any true latitide(s).
c   It is assumed that premaptform has been called prior to this so
c   that the proper constants have been placed in the common block
c   called mptf, which should be declared in (and only in) the
c   main program and routines maptform (this routine) and premaptform.
c
      common /mptf/ rpd_mptf,pi_mptf,dskmc_mptf,xlonc_mptf,rearth_mptf,
     &   ciy_mptf,cjx_mptf,yc_mptf,ihm_mptf,c1_mptf,c2_mptf,cone_mptf,
     &   conei_mptf,nproj_mptf
c
      if (idir.eq.1) then   ! First, deal with idir=1
c
      ypoint=(riy-ciy_mptf)*dskmc_mptf+yc_mptf
      xpoint=(rjx-cjx_mptf)*dskmc_mptf
c
      if (nproj_mptf.eq.3) then
         rlat=(2.*atan(exp(ypoint/rearth_mptf))-.5*pi_mptf)/rpd_mptf
         rlon=xlonc_mptf+(xpoint/rearth_mptf)/rpd_mptf
      elseif (nproj_mptf.eq.1) then
         rlat=(.5*ihm_mptf*pi_mptf-2.*atan(c2_mptf*
     &      (sqrt(xpoint**2+ypoint**2))**conei_mptf))/rpd_mptf
         rlon=xlonc_mptf+(conei_mptf*atan2(xpoint,
     &      -ihm_mptf*ypoint))/rpd_mptf
      elseif (nproj_mptf.eq.2) then
         rlat=(.5*ihm_mptf*pi_mptf-ihm_mptf*2.*atan(sqrt(xpoint**2+
     &      ypoint**2)/(rearth_mptf*c1_mptf)))/rpd_mptf
         if(xpoint.eq.0..and.ypoint.eq.0.) then
            rlon=xlonc_mptf
         else
            rlon=xlonc_mptf+(atan2(xpoint,-ihm_mptf*ypoint))/rpd_mptf
         endif
      endif
c MGD mod
c since everything is initially in gridpoints, when we convert to lat/lon
c we should add in the rotation to the final longitude.
c going from lat/lon to gridpoints needs nothing, since lat/lon will have
c been converted, and thus already have accounted for the rotation.
      rlon=mod(rlon+900.+rrota,360.)-180.
c
      else   ! Otherwise, deal with idir=-1
c
      dlon=rlon-xlonc_mptf
      if (dlon.lt.-180.) dlon=dlon+360
      if (dlon.gt. 180.) dlon=dlon-360
      if (nproj_mptf.eq.3) then
         ypoint=rearth_mptf*log((1.+sin(rpd_mptf*rlat))/
     &      cos(rpd_mptf*rlat))
         xpoint=dlon*rpd_mptf*rearth_mptf
      elseif (nproj_mptf.eq.1) then
         ypoint=-c1_mptf*(ihm_mptf*tan(.25*(ihm_mptf*pi_mptf-
     &      2.*rpd_mptf*rlat)))**cone_mptf*cos(cone_mptf*rpd_mptf*dlon)
         xpoint=ihm_mptf*c1_mptf*(ihm_mptf*tan(.25*(ihm_mptf*pi_mptf-
     &      2.*rpd_mptf*rlat)))**cone_mptf*sin(cone_mptf*rpd_mptf*dlon)
      elseif (nproj_mptf.eq.2) then
         ypoint=-rearth_mptf*sin(.5*ihm_mptf*pi_mptf-rpd_mptf*rlat)*
     &      c1_mptf/(1.+cos(.5*ihm_mptf*pi_mptf-rpd_mptf*rlat))*
     &      cos(rpd_mptf*dlon)
         xpoint=ihm_mptf*rearth_mptf*sin(.5*ihm_mptf*pi_mptf-rpd_mptf*
     &      rlat)*c1_mptf/(1.+cos(.5*ihm_mptf*pi_mptf-rpd_mptf*rlat))*
     &      sin(rpd_mptf*dlon)
      endif
      riy=(ypoint-yc_mptf)/dskmc_mptf+ciy_mptf
      rjx=xpoint/dskmc_mptf+cjx_mptf
c
      endif
c
      return
      end