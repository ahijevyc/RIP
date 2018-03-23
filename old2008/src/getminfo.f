c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine getminfo(casename,iendc,igotminfo,minfo,iendmi,iup)
c
      character casename*(*),fname*256
      character minfo*256
c
c   Read in the minfo string from the .minfo file.
c
      igotminfo=0
      fname=casename(1:iendc)//'.minfo'
      open(unit=25,file=fname,err=1021,form='formatted',status='old')
      read(25,'(a)',end=1021,err=1021) minfo
      do im=132,1,-1
         if (minfo(im:im).ne.' ') then
            iendmi=im
            goto 109
         endif
      enddo
      iendmi=0
 109  continue
      if (iendmi.gt.0) then
         write(iup,*) minfo(1:iendmi)
c        call flush(iup)
         igotminfo=1
      endif
      goto 1022
 1021 write(iup,*) 'Model info file not available, so model info'
      write(iup,*) 'won''t be plotted at bottom of frame. Proceeding.'
 1022 continue
      close (25)
c
      return
      end
