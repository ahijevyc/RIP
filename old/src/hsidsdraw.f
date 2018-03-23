c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hsidsdraw(rtslb,ixwin,iywin,icolr,rcrag,
     &   icaoid, maxpl,ipl,rrota)
c
      dimension ixwin(2,maxpl),iywin(2,maxpl),rcrag(2,maxpl),
     &   icolr(maxpl),rtslb(maxpl),rrota(maxpl)
      character icaoid*4
      real cv,ch,vm,hm,flmino,frmaxo,fbmino,ftmaxo
c
c routine to plot station IDs on the map background (ptyp=hb)
c
      include 'comconst'
c
c   Make proper set call, as well as number of values to
c      plot in each direction
c
c MGD begin mod
c modify min and max so that when aspect ratio is calculated, we 
c get useful value for a rotated rectangular window
      if (rrota(ipl) .eq. 90. .or. rrota(ipl) .eq. -90.) then
        yintervs=ixwin(2,ipl)-ixwin(1,ipl)
        xintervs=iywin(2,ipl)-iywin(1,ipl)
        cv = (ftmax-fbmin)/2.
        ch = (frmax-flmin)/2.
        vm = (ftmax+fbmin)/2.
        hm = (frmax+flmin)/2.
        fbmino=fbmin
        ftmaxo=ftmax
        flmino=flmin
        frmaxo=frmax
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
c different corner coordinates for call to set()
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
c   Draw the text
c
      setfrac=(fr-fl)/(ur-ul)*50.
      chsiz=max(.004,setfrac*rtslb(ipl))
      call gstxci(icolr(ipl))
      call gsplci(icolr(ipl))
      caxgn=1.+(rcrag(2,ipl)-xjcorn)*refrat
      caygn=1.+(rcrag(1,ipl)-yicorn)*refrat

c MGD begin mod
c transpose/invert x and y coordinates of text as necessary...
      if(rrota(ipl).eq.90.) then
        call plchhq(iywin(2,ipl)+iywin(1,ipl)-caygn+1,caxgn,
     &    icaoid,chsiz,0.,0.)
      elseif(rrota(ipl).eq.-90.) then
        call plchhq(caygn,ixwin(2,ipl)+ixwin(1,ipl)-caxgn+1,
     &    icaoid,chsiz,0.,0.)
      elseif(rrota(ipl).eq.-180..or.rrota(ipl).eq.180.) then
        call plchhq(ixwin(2,ipl)+ixwin(1,ipl)-caxgn+1,
     &  iywin(2,ipl)+iywin(1,ipl)-caygn+1,icaoid,chsiz,0.,0.)
c MGD end mod
      else
        call plchhq(caxgn,caygn,icaoid,chsiz,0.,0.)
      endif
c
      call gsplci(1)
      call gstxci(1)
      call gspmci(1)
      call gsfaci(1)
c
c MGD begin mod
c reset the window parameters if we have changed them
      if(rrota(ipl) .eq. 90. .or. rrota(ipl) .eq. -90.) then
        fbmino=fbmin
        ftmaxo=ftmax
        flmino=flmin
        frmaxo=frmax
      endif
c MGD end mod
      return
      end
