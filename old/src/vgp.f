c--------------------------------------------------------------------
      subroutine vgp (u, v, ght, ter, cape, g, miy, mjx, mkzh)
c Routine to compute vorticity-generation potential
c Compute this after doing cape (it reads the cape in from a scratch
c array to save time).
      dimension u(miy,mjx,mkzh), v(miy,mjx,mkzh), ght(miy,mjx,mkzh)
     & ,ter(miy,mjx), cape(miy,mjx), g(miy,mjx)
      parameter (pi=3.14159265, dtr=pi/180., dpr=180./pi)
c     rewind(58)
c     read(58,*) iii,jjj
      do j = 1, mjx-1
	do i = 1, miy-1
c find the index nearest 3 km AGL
	  k3 = 0
	  do 6 k = mkzh, 2, -1
	    if (((ght(i,j,k) - ter(i,j)) .gt. 3000.)) then
	      k3 = k
	      go to 8
	    endif
    6     continue
    8     continue
c calculate the 0-3 km hodograph length
	  totshr = 0.
	  depth=ght(i,j,k3)-ght(i,j,mkzh)
	  do k = mkzh-1, k3, -1
	    su = abs(u(i,j,k+1) - u(i,j,k))
	    sv = abs(v(i,j,k+1) - v(i,j,k))
	    totshr = sqrt(su*su + sv*sv) + totshr
	  enddo
	  totshr = totshr/depth
	  g(i,j) = sqrt(cape(i,j)) * totshr
c         if (i.eq.iii .and. j.eq.jjj) then
c           write(6,*) 'cape = ',cape(i,j),' totshr = ',totshr
c           write(6,*) 'k3 = ',k3,' ght(i,j,k3) = ',ght(i,j,k3),
c    & ' elev = ',ter(i,j)
c           write(6,*) 'depth = ',depth,' vgp = ',g(i,j)
c         endif
        enddo
      enddo
      end
