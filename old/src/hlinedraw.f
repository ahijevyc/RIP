c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hlinedraw(ilinw,ixwin,iywin,icolr,rcrag,rcrbg,
     &         maxpl,ipl,rrota)
c
      dimension ixwin(2,maxpl),iywin(2,maxpl),rcrag(2,maxpl),
     &   rcrbg(2,maxpl),icolr(maxpl),ilinw(maxpl),rrota(maxpl)
      real cv,ch,vm,hm
c
      include 'comconst'
c
c   Set line width, color
c
      lwidth=ilinw(ipl)*1000
      call setusv('LW',lwidth)
      call gsplci(icolr(ipl))
      call gstxci(icolr(ipl))
c
c   Make proper set call, as well as number of values to
c      plot in each direction
c
c MGD begin mod
c if window is rectangular (assume it is) then set needs to be called
c with different parameters to accomodate rotated rectangle
      if (rrota(ipl) .eq. 90. .or. rrota(ipl) .eq. -90.) then
        yintervs=ixwin(2,ipl)-ixwin(1,ipl)
        xintervs=iywin(2,ipl)-iywin(1,ipl)
        cv = (ftmax-fbmin)/2.
        ch = (frmax-flmin)/2.
        vm = (ftmax+fbmin)/2.
        hm = (frmax+flmin)/2.
        fbmin=vm-ch
        ftmax=vm+ch
        flmin=hm-cv
        frmax=hm+cv
c MGD end mod
      else
        xintervs=ixwin(2,ipl)-ixwin(1,ipl)
        yintervs=iywin(2,ipl)-iywin(1,ipl)
      endif
      faspect=(ftmax-fbmin)/(frmax-flmin)
      aspect=yintervs/xintervs
      if (aspect.lt.faspect) then
         fl=flmin
         fr=frmax
         fextra=.5*((ftmax-fbmin)-aspect*(frmax-flmin))
         fb=fbmin+fextra
         ft=ftmax-fextra
      else
         fb=fbmin
         ft=ftmax
         fextra=.5*((frmax-flmin)-1./aspect*(ftmax-fbmin))
         fl=flmin+fextra
         fr=frmax-fextra
      endif
c MGD begin mod
c switch x and y dimensions
      if (rrota(ipl) .eq. 90. .or. rrota(ipl) .eq. -90.) then
        ul=float(iywin(1,ipl))
        ur=float(iywin(2,ipl))
        ub=float(ixwin(1,ipl))
        ut=float(ixwin(2,ipl))
c MGD end mod
      else
        ul=float(ixwin(1,ipl))
        ur=float(ixwin(2,ipl))
        ub=float(iywin(1,ipl))
        ut=float(iywin(2,ipl))
      endif
      call set(fl,fr,fb,ft,ul,ur,ub,ut,1)
c
c   Draw the line
c
      caxgn=1.+(rcrag(2,ipl)-xjcorn)*refrat
      caygn=1.+(rcrag(1,ipl)-yicorn)*refrat
      cbxgn=1.+(rcrbg(2,ipl)-xjcorn)*refrat
      cbygn=1.+(rcrbg(1,ipl)-yicorn)*refrat
      write(6,*) 'line rota:', rrota(ipl)
c MGD begin mod
c pretty straight forward transposition/inversion of line 
c end-point coordinated
      if(rrota(ipl) .eq. 90.) then
        call line(iywin(2,ipl)+iywin(1,ipl)-caygn+1,caxgn,
     &    iywin(2,ipl)+iywin(1,ipl)-cbygn+1,cbxgn)
      elseif(rrota(ipl) .eq. -90.) then
        call line(caygn,ixwin(2,ipl)+ixwin(1,ipl)-caxgn+1,cbygn,
     &    ixwin(2,ipl)+ixwin(1,ipl)-cbxgn+1)
      elseif(rrota(ipl) .eq. 180. .or. rrota(ipl) .eq. -180.) then
        call line(ixwin(2,ipl)+ixwin(1,ipl)-caxgn+1,iywin(2,ipl)+
     &    iywin(1,ipl)-caygn+1,ixwin(2,ipl)+ixwin(1,ipl)-cbxgn+1,
     &    iywin(2,ipl)+iywin(1,ipl)-cbygn+1)
c MGD end mod
      else
        call line(caxgn,caygn,cbxgn,cbygn)
      endif
c
      call setusv('LW',1000)
      call gsplci(1)
      call gstxci(1)
c
      return
      end
