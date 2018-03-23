c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine vinterp(cvcord,rlevel,ixjs,iyis,icdwork,vc3d,tmk,qvp,
     &   prs,ght,ter,pstx,sigh,sigf,prs_tsf,lhide,idiffflag,field,
     &   work,pslab,iprog,mabpl,morpl,njx,niy,miy,mjx,mkzh)
c
c   Assumption: Everything varies linearly with sigma ( and therefore
c      with pressure) expept height (Z). Instead, exp(-Z/H) varies
c      linearly with sigma, where H is a constant scale height.
c      Therefore, when interpolating Z to a pressure or isentropic
c      surface, first convert it to exp(-Z/H), then interpolate,
c      then convert back. Also, when interpolating everything else
c      to Z surfaces, instead of using plain old Z as the vert. coord.,
c      use exp(-Z/H).
c
      dimension vc3d(miy,mjx,mkzh),tmk(miy,mjx,mkzh),qvp(miy,mjx,mkzh),
     &   work(miy,mjx,mkzh),pslab(mabpl,morpl),
     &   ght(miy,mjx,mkzh),ter(miy,mjx),pstx(miy,mjx),sigh(mkzh),
     &   sigf(mkzh+1),prs(miy,mjx,mkzh),prs_tsf(miy,mjx)
      logical lhide, dosfcprs
      character cvcord*1,cvc*1,field*10
c
      include 'comconst'
c
      expon=rgas*ussalr/grav
      exponi=1./expon
      iwarn = 0
c
      cvc=cvcord
      if (cvcord.eq.'l'.or.cvcord.eq.'x') cvc='p'
c
c   Trivial cases: field and vert. coord. are the same.
c
      itriv=0
      if (cvc.eq.'z'.and.field(1:4).eq.'ght ') then
         itriv=1
      elseif (cvc.eq.'p'.and.field(1:4).eq.'prs ') then
         itriv=2
      endif
c
c   Special cases: special extrapolation based on US standard atm.
c   will be performed.
c
      icase=0
      if (cvc.eq.'p'.or.cvc.eq.'z') then
         if (field(1:4).eq.'ght '.or.field(1:4).eq.'prs ') then
            icase=1
         elseif (field(1:4).eq.'tmk ') then
            icase=2
         elseif (field(1:4).eq.'tmc ') then
            icase=3
         elseif (field(1:4).eq.'the ') then
            icase=4
         elseif (field(1:4).eq.'eth ') then
            icase=5
         endif
      endif
c
c   Set flag for using surface pressure as vlev.  Should be true if
c   data is pressure-level data, requested vert. coord. is pressure,
c   and requested level is 1001.
c
      dosfcprs=.false.
      if (iprog.le.3.and.cvc.eq.'p'.and.nint(rlevel).eq.1001)
     &   dosfcprs=.true.
c
c   Get desired sigma level on which to base subterranean
c   temperature profile.
c
      if (icase.gt.0) then
         sigc=(850.-ptop)/(1000.-ptop)
         do k=1,mkzh
            if (sigc.ge.sigf(k).and.sigc.le.sigf(k+1)) kupper=k
         enddo
      endif
c
      if (cvc.eq.'z') then
         vlev=exp(-rlevel*1000./sclht)
      else
         vlev=rlevel
      endif
c
      do j=1,njx
      jj=j+ixjs-1
      jjph=min(jj,mjx-1)
      jjmh=max(jj-1,1)
      do i=1,niy
      ii=i+iyis-1
      iiph=min(ii,miy-1)
      iimh=max(ii-1,1)
c
c      if (lhide) then  ! just for testing purposes
c         pslab(j,i)=vc3d(ii,jj,1)
c         goto 333
c      endif
c
c   Set vlev to true surface pressure if flag is true.
c
      if (dosfcprs) vlev=prs_tsf(ii,jj)
c
c   Trival cases
c
      if (itriv.eq.1) then
         pslab(j,i)=-sclht*log(vlev)
         goto 333
      elseif (itriv.eq.2) then
         pslab(j,i)=vlev
         goto 333
      endif
c
      ifound=0
      do k=1,mkzh-1
         if (icdwork.eq.0) then
            vcp1=.25*(vc3d(iimh,jjmh,k+1)+vc3d(iimh,jjph,k+1) +
     &                vc3d(iiph,jjmh,k+1)+vc3d(iiph,jjph,k+1) )
            vcp0=.25*(vc3d(iimh,jjmh,k)+vc3d(iimh,jjph,k) +
     &                vc3d(iiph,jjmh,k)+vc3d(iiph,jjph,k) )
         else
            vcp1=vc3d(ii,jj,k+1)
            vcp0=vc3d(ii,jj,k)
         endif
         if ( (vlev.ge.vcp0.and.vlev.le.vcp1) .or.
     &        (vlev.le.vcp0.and.vlev.ge.vcp1)     ) then
            if (work(ii,jj,k).eq.rmsg.or.work(ii,jj,k+1).eq.rmsg)then
               pslab(j,i)=rmsg
               goto 105
            endif
            if (field(1:4).eq.'ght ') then
               valp0=exp(-work(ii,jj,k)/sclht)
               valp1=exp(-work(ii,jj,k+1)/sclht)
            else
               valp0=work(ii,jj,k)
               valp1=work(ii,jj,k+1)
            endif
            pslab(j,i)=(vlev-vcp0)*(valp1-valp0)/(vcp1-vcp0)+
     &                 valp0
            if (field(1:4).eq.'ght ') then
               pslab(j,i)=-sclht*log(pslab(j,i))
            endif
  105       ifound=1
            goto 115
         endif
      enddo
 115  continue
c
      if (ifound.eq.1) goto 333
c
c   Grid point is outside (above or below) model domain
c   Set to missing value if "hide" was specified, so grid point
c   will be excluded from plot.
c
      if (lhide) then
         pslab(j,i)=rmsg
         goto 333
      endif
c
      if (icdwork.eq.0) then
         vclhsl=.25*(vc3d(iimh,jjmh,mkzh)+
     &      vc3d(iimh,jjph,mkzh)+vc3d(iiph,jjmh,mkzh)+
     &      vc3d(iiph,jjph,mkzh) )
         vctophsl=.25*(vc3d(iimh,jjmh,1)+
     &      vc3d(iimh,jjph,1)+vc3d(iiph,jjmh,1)+
     &      vc3d(iiph,jjph,1) )
      else
         vclhsl=vc3d(ii,jj,mkzh)
         vctophsl=vc3d(ii,jj,1)
         ttlhsl=tmk(ii,jj,mkzh)
         qvlhsl=qvp(ii,jj,mkzh)
      endif
      diff=vctophsl-vclhsl
      isign=nint(diff/abs(diff))
c
c   If "hide" not specified, set value for grid points above model top
c   to value at highest half sigma level.
c
      if (isign*vlev.ge.isign*vctophsl) then
         pslab(j,i)=work(ii,jj,1)
         iwarn=iwarn+1
         goto 333
      endif
c
c   Only remaining possibility is that the specified level is below
c   lowest half sigma level.  If lowest half-sigma-level value is missing,
c   set interpolated value to missing.
c
      if (work(ii,jj,mkzh).eq.rmsg) then
         pslab(j,i)=rmsg
         goto 333
      endif
c
c   In most other cases, set the value to that
c   at the lowest half sigma level.
c
      pslab(j,i)=work(ii,jj,mkzh)
c
c   For the special cases of pressure on height levels or height on
c   pressure levels, or temperature-related variables on pressure or
c   height levels, perform a special extrapolation based on
c   US Standard Atmosphere.
c
c   Note, this special extrapolation will not be performed if this is
c   a differenced field, because the differenced information will
c   be lost.  In this case, the lowest half sigma level value should
c   be retained.
c
      if (idiffflag.eq.1) goto 333
c
      if (icase.gt.0) then
c
      plhsl=prs(ii,jj,mkzh)
      zlhsl=ght(ii,jj,mkzh)
      ezlhsl=exp(-zlhsl/sclht)
      tlhsl=tmk(ii,jj,mkzh)
      psurf=plhsl+(1.-sigh(mkzh))*pstx(ii,jj)
      zsurf=ter(ii,jj)
      ezsurf=exp(-zsurf/sclht)
c
      if ((cvc.eq.'z'.and.vlev.gt.ezsurf).or.
     &    (cvc.eq.'p'.and.vlev.lt.psurf)) then
c
c      We are below the lowest half sigma level but above the ground.
c
         if (cvc.eq.'p') then
            plev=vlev
            ezlev=((plev-plhsl)*ezsurf+(psurf-plev)*ezlhsl)/
     &          (psurf-plhsl)
            zlev=-sclht*log(ezlev)
            if (field(1:4).eq.'ght ') then
               pslab(j,i)=zlev
               goto 333
            endif
         elseif (cvc.eq.'z') then
            ezlev=vlev
            zlev=-sclht*log(ezlev)
            plev=((ezlev-ezlhsl)*psurf+(ezsurf-ezlev)*plhsl)/
     &          (ezsurf-ezlhsl)
            if (field(1:4).eq.'prs ') then
               pslab(j,i)=plev
               goto 333
            endif
         endif
c
      else
c
c      We are below ground.
c
         tsurfbg=tmk(ii,jj,kupper)*(psurf/prs(ii,jj,kupper))**expon
         tsurfbgv=virtual(tsurfbg,qvp(ii,jj,mkzh))
         if (cvc.eq.'p') then
            plev=vlev
            zlev=zsurf+tsurfbgv/ussalr*(1.-(vlev/psurf)**expon)
            if (field(1:4).eq.'ght ') then
               pslab(j,i)=zlev
               goto 333
            endif
         elseif (cvc.eq.'z') then
            zlev=-sclht*log(vlev)
            plev=psurf*(1.+ussalr/tsurfbgv*(zsurf-zlev))**exponi
            if (field(1:4).eq.'prs ') then
               pslab(j,i)=plev
               goto 333
            endif
         endif
c
      endif
c
      if (icase.ge.2) then
c
c      One of the temperature fields was requested.  Extrapolate from
c      T at the lowest half sigma level using US standard lapse rate.
c      Assume qv at the requested level is the same as at the lowest
c      half sigma level.
c
         tlev=tlhsl+(zlhsl-zlev)*ussalr
         qvlev=qvp(ii,jj,mkzh)
         if (icase.eq.2) then
            pslab(j,i)=tlev
         elseif (icase.eq.3) then
            pslab(j,i)=tlev-celkel
         elseif (icase.eq.4) then
            gammam=gamma*(1.+gammamd*qvlev)
            pslab(j,i)=tlev*(1000./plev)**gammam
         elseif (icase.eq.5) then
            e=qvlev*plev/(eps+qvlev)
            tlcl=tlclc1/(log(tlev**tlclc2/e)-tlclc3)+tlclc4
            pslab(j,i)=tlev*(1000./plev)**(gamma*(1.+gammamd*qvlev))*
     &         exp((thtecon1/tlcl-thtecon2)*
     &         qvlev*(1.+thtecon3*qvlev))
         endif
c
      endif
c
      endif  ! end of special extrapolation
c
 333  continue
c
      enddo
      enddo
c
      if (iwarn .gt. 0) then
	write(iup,*) 'Above model top at ',iwarn,' points --',
     &             ' using model top data'
      endif
      return
      end
