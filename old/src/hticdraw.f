c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hticdraw(ngfbuf,ilinw,ixwin,iywin,raxlg,raxtg,icolr,
     &         rtslb,ilinwgf,ixwingf,iywingf,raxlggf,raxtggf,icolrgf,
     &         rtslbgf,iwhatgf,maxbuf,maxpl,ipl,rrota,rrotagf)
c
      dimension ixwin(2,maxpl),iywin(2,maxpl),raxtg(maxpl),
     &   icolr(maxpl),raxlg(maxpl),ilinw(maxpl),rrotagf(maxbuf),
     &   rtslb(maxpl),ixwingf(2,maxbuf),rtslbgf(maxbuf),
     &   iywingf(2,maxbuf),raxtggf(maxbuf),raxlggf(maxbuf),
     &   ilinwgf(maxbuf),icolrgf(maxbuf),iwhatgf(maxbuf),rrota(maxpl)
c
      character axlab*3
      real temp,cv,ch,vm,hm
      real fbmino,ftmaxo,flmino,frmaxo
c
      include 'comconst'
c
c   Check if this tic background is identical to any of those saved
c      in the GFLASH buffers, and if so, flush that buffer to
c      this frame.
c
      do 10 ib=1,ngfbuf
         if (iwhatgf(ib).ne.3) goto 10
         if ( ixwin(1,ipl).eq.ixwingf(1,ib).and.ixwin(2,ipl).eq.
     &        ixwingf(2,ib).and.iywin(1,ipl).eq.iywingf(1,ib).and.
     &        iywin(2,ipl).eq.iywingf(2,ib).and.raxtg(ipl).eq.
     &        raxtggf(ib).and.raxlg(ipl).eq.raxlggf(ib).and.
     &        ilinw(ipl).eq.ilinwgf(ib).and.icolr(ipl).eq.
     &        icolrgf(ib).and.rtslb(ipl).eq.
     &        rtslbgf(ib)) then
c MGD begin mod
c If rotations are not the same for two sets of ticks, then the ticks 
c are not really the same, so we must redraw...
c        if((ipl .gt. 1 .and. rrota(ipl) .eq. rrota(ib)) .or.
c    &      ipl .eq. 1) then
         if(rrota(ipl) .eq. rrotagf(ib)) then
            call gflas3(ib)
            goto 120
         endif
c MGD end mod
         endif
   10 continue
c
c   If identical tic background was not found in GFLASH buffer,
c      draw the tic background to a new GFLASH buffer, flush this
c      new buffer to the current frame, and save the parameters for
c      this new tic background.
c
      ngfbuf=ngfbuf+1
      if (ngfbuf.gt.maxbuf) then
         write(iup,*)'In hticdraw: too many gflash buffers.'
         stop
      endif
      call gflas1(ngfbuf)
c
c   Set line width, color
c
      lwidth=ilinw(ipl)*1000
      call setusv('LW',lwidth)
      call gsplci(icolr(ipl))
      call gstxci(icolr(ipl))
c
c   Make set call and draw perimeter.
c
c MGD begin mod
c account for non-square grids being rotated - gotta fit in new window
c with same center as before
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
      else
c MGD end mod
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
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call line(fl,fb,fl,ft)
      call line(fr,fb,fr,ft)
      call line(fl,ft,fr,ft)
      call line(fl,fb,fr,fb)
c
c   Draw tick marks and label axes.
c
c MGD begin mod
c basically, if grid is rotated by (-)90, then coordinates of gridpoints
c also get transposed and possibly inverted...
      if (rrota(ipl) .eq. 90. .or. rrota(ipl) .eq. -90.) then
      xdist=(.6*rtslb(ipl)+.01)
      ydist=.01
      tickl=.005
      gint=(fr-fl)/xintervs
      if (raxlg(ipl).eq.rmsg) then
         iaxlg=10
      else
         iaxlg=nint(raxlg(ipl))
      endif
      if (raxtg(ipl).eq.rmsg) then
         iaxtg=1
      else
         iaxtg=nint(raxtg(ipl))
      endif
c
      if (iaxlg.ne.0) then
         do j=iywin(1,ipl),iywin(2,ipl)
            xtick=fl+(j-iywin(1,ipl))*gint
            if(mod(j,iaxlg).eq.0) then
               if(rrota(ipl) .eq. -90.) then
                 write(axlab,910) j
c if 90, then direction of increase for x labels will reverse
               elseif(rrota(ipl) .eq. 90.) then
                 write(axlab,910) iywin(2,ipl)+iywin(1,ipl)-j
               endif
               iaxs=3-int(alog10(float(j)+.5))
               call setusv('LW',1000)
c              call plchhq(xtick,fb-ydist,axlab(iaxs:),
               call plchhq(xtick,fb-ydist,axlab,
     &            rtslb(ipl),0.,0.)
               call setusv('LW',lwidth)
               call line(xtick,fb,xtick,fb+2.*tickl)
               call line(xtick,ft,xtick,ft-2.*tickl)
            endif
         enddo
         do i=ixwin(1,ipl),ixwin(2,ipl)
            ytick=fb+(i-ixwin(1,ipl))*gint
            if(mod(i,iaxlg).eq.0) then
c if -90 rotation, then direction of increase for y labels will reverse
               if(rrota(ipl) .eq. -90.) then
                 write(axlab,910) ixwin(2,ipl)+ixwin(1,ipl)-i
               elseif(rrota(ipl) .eq. 90.) then
                 write(axlab,910) i
               endif
               call setusv('LW',1000)
               call plchhq(fl-xdist,ytick,axlab,rtslb(ipl),0.,1.)
               call setusv('LW',lwidth)
               call line(fl,ytick,fl+2.*tickl,ytick)
               call line(fr,ytick,fr-2.*tickl,ytick)
            endif
         enddo
      endif
      if (iaxtg.ne.0) then
         do j=iywin(1,ipl),iywin(2,ipl)
            xtick=fl+(j-iywin(1,ipl))*gint
            if(mod(j,iaxtg).eq.0) then
               call line(xtick,fb,xtick,fb+tickl)
               call line(xtick,ft,xtick,ft-tickl)
            endif
         enddo
         do i=ixwin(1,ipl),ixwin(2,ipl)
            ytick=fb+(i-ixwin(1,ipl))*gint
            if(mod(i,iaxtg).eq.0) then
               call line(fl,ytick,fl+tickl,ytick)
               call line(fr,ytick,fr-tickl,ytick)
            endif
         enddo
      endif
c MGD end mod
      else
      ydist=(.6*rtslb(ipl)+.01)
      xdist=.01
      tickl=.005
      gint=(fr-fl)/xintervs
      if (raxlg(ipl).eq.rmsg) then
         iaxlg=10
      else
         iaxlg=nint(raxlg(ipl))
      endif
      if (raxtg(ipl).eq.rmsg) then
         iaxtg=1
      else
         iaxtg=nint(raxtg(ipl))
      endif
c
      if (iaxlg.ne.0) then
         do j=ixwin(1,ipl),ixwin(2,ipl)
            xtick=fl+(j-ixwin(1,ipl))*gint
            if(mod(j,iaxlg).eq.0) then
c MGD mod
c if 180 degree rotation, then both labels need to be inverted
               if(rrota(ipl) .eq. -180. .or. 
     &                rrota(ipl) .eq. 180.) then
                 write(axlab,910) ixwin(2,ipl)+ixwin(1,ipl)-j
               else
                 write(axlab,910) j
               endif
               iaxs=3-int(alog10(float(j)+.5))
               call setusv('LW',1000)
c              call plchhq(xtick,fb-ydist,axlab(iaxs:),
               call plchhq(xtick,fb-ydist,axlab,
     &            rtslb(ipl),0.,0.)
               call setusv('LW',lwidth)
               call line(xtick,fb,xtick,fb+2.*tickl)
               call line(xtick,ft,xtick,ft-2.*tickl)
            endif
         enddo
         do i=iywin(1,ipl),iywin(2,ipl)
            ytick=fb+(i-iywin(1,ipl))*gint
            if(mod(i,iaxlg).eq.0) then
c if 180 degree rotation, then both labels need to be inverted
               if(rrota(ipl) .eq. -180. .or. 
     &                rrota(ipl) .eq. 180.) then
                 write(axlab,910) iywin(2,ipl)+iywin(1,ipl)-i
               else
                 write(axlab,910) i
               endif
               call setusv('LW',1000)
               call plchhq(fl-xdist,ytick,axlab,rtslb(ipl),0.,1.)
               call setusv('LW',lwidth)
               call line(fl,ytick,fl+2.*tickl,ytick)
               call line(fr,ytick,fr-2.*tickl,ytick)
            endif
         enddo
      endif
      if (iaxtg.ne.0) then
         do j=ixwin(1,ipl),ixwin(2,ipl)
            xtick=fl+(j-ixwin(1,ipl))*gint
            if(mod(j,iaxtg).eq.0) then
               call line(xtick,fb,xtick,fb+tickl)
               call line(xtick,ft,xtick,ft-tickl)
            endif
         enddo
         do i=iywin(1,ipl),iywin(2,ipl)
            ytick=fb+(i-iywin(1,ipl))*gint
            if(mod(i,iaxtg).eq.0) then
               call line(fl,ytick,fl+tickl,ytick)
               call line(fr,ytick,fr-tickl,ytick)
            endif
         enddo
      endif
      endif
  910 format(i3)
c
      call gflas2
      call gflas3(ngfbuf)
      call setusv('LW',1000)
      call gsplci(1)
      call gstxci(1)
      iywingf(1,ngfbuf)=iywin(1,ipl)
      iywingf(2,ngfbuf)=iywin(2,ipl)
      ixwingf(1,ngfbuf)=ixwin(1,ipl)
      ixwingf(2,ngfbuf)=ixwin(2,ipl)
      raxtggf(ngfbuf)=raxtg(ipl)
      raxlggf(ngfbuf)=raxlg(ipl)
      ilinwgf(ngfbuf)=ilinw(ipl)
      icolrgf(ngfbuf)=icolr(ipl)
      rtslbgf(ngfbuf)=rtslb(ipl)
c MGD mod for buffering
      rrotagf(ib)=rrota(ipl)
      iwhatgf(ngfbuf)=3
c      write(iup,*)'   Made new tic background number ',ngfbuf,
c     &   ' for plot ',ipl
c
  120 continue
c MGD begin mod
c reset the values of the parameters we changed for -90 and 90 rotations
      if (rrota(ipl) .eq. 90. .or. rrota(ipl) .eq. -90.) then
        fbmin=fbmino
        ftmax=ftmaxo
        flmin=flmino
        frmax=frmaxo
      endif
c MGD end mod
      return
      end
