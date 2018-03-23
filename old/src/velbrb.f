c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine velbrb(u,lu,v,lv,m,n,rlength,ispv,spv,ixwin,iywin,
     &  ihv,rrota)
c
c   Everything is as in velvct except:
c      There is no hi and lo specification.
c      NSET is hardwired to 1
c      RLENGTH is the length, in grid spaces, of a one-full-barb
c      wind arrow.  If the grid is irregular, a grid length is taken
c      as an x-direction grid length in the center of the grid.
c
      dimension u(lu,n),v(lv,n),spv(2)
      include 'comconst'
c
c   Determine vector length in fractional coordinates.
c
      if ( xlatc .ge. 0. ) then   
        ins = 0
      else
	ins = 1
      endif
      xhalf=.5*m
      yhalf=.5*n
      gridfrac=cufx(fx(xhalf+.5,yhalf))-cufx(fx(xhalf-.5,yhalf))
      do 1000 iy=1,n
      do 1000 ix=1,m
c
      if (ispv.eq.1.and.u(ix,iy).eq.spv(1)) goto 999
      if (ispv.eq.2.and.v(ix,iy).eq.spv(2)) goto 999
      if (ispv.eq.3.and.(u(ix,iy).eq.spv(1).or.
     &    v(ix,iy).eq.spv(2))) goto 999
      if (ispv.eq.4.and.u(ix,iy).eq.spv(1).and.
     &    v(ix,iy).eq.spv(2)) goto 999
c
      rix=fx(float(ix),float(iy))
      riy=fy(float(ix),float(iy))
      if (riy .eq. float(iy) .and. ihv .eq. 0 ) then
	xpp = xjcorn + (rix+(ixwin-1) -1.)/refrat
	ypp = yicorn + (riy+(iywin-1) -1.)/refrat
cwrite(6,*) 'rix = ',rix,' riy = ',riy,' ix = ',ix,' iy = ',iy
cwrite(6,*) 'xpp = ',xpp,' ypp = ',ypp,' xjcorn = ',xjcorn,
c    & ' yicorn = ',yicorn
	call maptform(ypp,xpp,rlat,rlon,1,rrota)
cwrite(6,*) 'rlat = ',rlat
	if (rlat .lt. 0.) then
	  ins = 1
	else
	  ins = 0
	endif
      endif
      call barb(rix,riy,u(ix,iy),v(ix,iy),rlength*gridfrac,ins)
c
 999  continue
c
 1000 continue
      return
      end
