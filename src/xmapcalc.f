c                                                                     c
c*********************************************************************c
c                                                                     c
      function xmapcalc(riy,rjx,isw)
c
c     This function calculates the map scale factor at a coarse domain
c     dot grid point, <riy,rjx> (if isw=1), or at a latitude in degrees
c     (given in riy, if isw=2).  It works for Lambert Conformal (LC,1),
c     Polar Stereographic (ST,2), or Mercator (ME,3) projections, with
c     any true latitide(s).  It is assumed that premaptform has been
c     called prior to this so that the proper constants have been placed
c     in the common block called mptf, which should be declared in (and
c     only in) the main program and routines mapfaccalc, (this routine),
c     maptform, and premaptform.
c
      common /mptf/ rpd_mptf,pi_mptf,dskmc_mptf,xlonc_mptf,rearth_mptf,
     &   ciy_mptf,cjx_mptf,yc_mptf,ihm_mptf,c1_mptf,c2_mptf,cone_mptf,
     &   conei_mptf,nproj_mptf
c
      if (isw.eq.1) then   ! Get latitude
         call maptform(riy,rjx,rlat,rlon,1)
      else                 ! Assume riy holds latitude
         rlat=riy
      endif
      rlatr=rpd_mptf*rlat
c
      if (nproj_mptf.eq.1) then
         angle=.25*(pi_mptf-2.*ihm_mptf*rlatr)
         xmapcalc=ihm_mptf*cone_mptf*c1_mptf/rearth_mptf*
     &      tan(angle)**cone_mptf/sin(2.*angle)
      elseif (nproj_mptf.eq.2) then
         xmapcalc=c1_mptf/(1.+cos(.5*ihm_mptf*pi_mptf-rlatr))
      elseif (nproj_mptf.eq.0.or.nproj_mptf.eq.3) then ! Mercator
c
c       (Note: nproj=0 means "idealized" (i.e. no map), but we'll treat
c        it as Mercator so nothing in RIP goes haywire.)
c
         xmapcalc=1./cos(rlatr)
      endif
c
      return
      end
