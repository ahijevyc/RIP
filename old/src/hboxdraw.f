c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hboxdraw(ilinw,ixwin,iywin,icolr,rcrag,rcrbg,
     &         maxpl,ipl,rrota)
c
      dimension ixwin(2,maxpl),iywin(2,maxpl),rcrag(2,maxpl),
     &   rcrbg(2,maxpl),icolr(maxpl),ilinw(maxpl),rrota(maxpl)
c
      include 'comconst'
c
c   Set line width, color
c
      write(6,*) 'HBOXDRAW'
      lwidth=ilinw(ipl)*1000
      call setusv('LW',lwidth)
      call gsplci(icolr(ipl))
      call gstxci(icolr(ipl))
c
c   Make proper set call, as well as number of values to
c      plot in each direction
c
      xintervs=ixwin(2,ipl)-ixwin(1,ipl)
      yintervs=iywin(2,ipl)-iywin(1,ipl)
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
      ul=float(ixwin(1,ipl))
      ur=float(ixwin(2,ipl))
      ub=float(iywin(1,ipl))
      ut=float(iywin(2,ipl))
      call set(fl,fr,fb,ft,ul,ur,ub,ut,1)
c
c   Draw the box
c
      caxgn=1.+(rcrag(2,ipl)-xjcorn)*refrat
      caygn=1.+(rcrag(1,ipl)-yicorn)*refrat
      cbxgn=1.+(rcrbg(2,ipl)-xjcorn)*refrat
      cbygn=1.+(rcrbg(1,ipl)-yicorn)*refrat
      call line(caxgn,caygn,caxgn,cbygn)
      call line(caxgn,cbygn,cbxgn,cbygn)
      call line(cbxgn,cbygn,cbxgn,caygn)
      call line(cbxgn,caygn,caxgn,caygn)
c
      call setusv('LW',1000)
      call gsplci(1)
      call gstxci(1)
c
      return
      end
