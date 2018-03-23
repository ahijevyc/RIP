c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hbulldraw(rtslb,ctitl,ixwin,iywin,icolr,rcrag,
     &   maxpl,ipl,rrota)
c
      dimension ixwin(2,maxpl),iywin(2,maxpl),rcrag(2,maxpl),
     &   icolr(maxpl),rtslb(maxpl),rrota(maxpl)
      character ctitl(maxpl)*82
c
      include 'comconst'
c
c   Make proper set call, as well as number of values to
c      plot in each direction
c
c MGD begin mod
c reverse dimensions for call to set for (-)90 degree rotations
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
c   Draw the bullet
c
      caxgn=1.+(rcrag(2,ipl)-xjcorn)*refrat
      caygn=1.+(rcrag(1,ipl)-yicorn)*refrat
      if (ctitl(ipl)(1:8).eq.'auto    ') then
         size=xintervs/(fr-fl)*rtslb(ipl)
c MGD begin mod
c transpose/invert x and y coordinates of the bullet
         if(rrota(ipl) .eq. 90.) then
           call ngdots(iywin(2,ipl)-caygn+1,caxgn,1,size,icolr(ipl))
         elseif(rrota(ipl) .eq. -90.) then
           call ngdots(caygn,ixwin(2,ipl)-caxgn+1,1,size,icolr(ipl))
         elseif(rrota(ipl) .eq. 180. .or. rrota(ipl) .eq. -180.) then
           call ngdots(ixwin(2,ipl)-caxgn+1,iywin(2,ipl)-caygn+1,
     &                 1,size,icolr(ipl))
c MGD end mod
         else
           call ngdots(caxgn,caygn,1,size,icolr(ipl))
         endif
      else
         if (ctitl(ipl)(82:82).ne.' ') then
            iendct=82
         else
            do ich=82,1,-1
               if (ctitl(ipl)(ich:ich).ne.' ') then
                  iendct=ich
                  goto 39
               endif
            enddo
            iendct=0
 39         continue
         endif
         if (iendct.lt.1) goto 30
         call gstxci(icolr(ipl))
         call gsplci(icolr(ipl))
c MGD begin mod
c same deal, but for the text
         if(rrota(ipl) .eq. 90.) then
           call plchhq(iywin(2,ipl)+iywin(1,ipl)-caygn+1,caxgn,
     &       ctitl(ipl)(1:iendct),rtslb(ipl),0.,0.)
         elseif(rrota(ipl) .eq. -90.) then
           call plchhq(caygn,ixwin(2,ipl)+ixwin(1,ipl)-caxgn+1,
     &       ctitl(ipl)(1:iendct),rtslb(ipl),0.,0.)
         elseif(rrota(ipl) .eq. 180. .or. rrota(ipl) .eq. -180.) then
           call plchhq(ixwin(2,ipl)+ixwin(1,ipl)-caxgn+1,
     &       iywin(2,ipl)+iywin(1,ipl)-caygn+1,ctitl(ipl)(1:iendct),
     &       rtslb(ipl),0.,0.)
c MGD end mod
         else
           call plchhq(caxgn,caygn,ctitl(ipl)(1:iendct),
     &                 rtslb(ipl),0.,0.)
         endif
      endif
 30   continue
c
      call gsplci(1)
      call gstxci(1)
      call gspmci(1)
      call gsfaci(1)
c
      return
      end
