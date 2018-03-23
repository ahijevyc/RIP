c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine frtitle(title,casename,iendc,rootname,
     &   mdateb,rhourb,xtime,timezone,iusdaylightrule,inearesth,
     &   idotitle,iinittime,ivalidtime,toptextclg,iprog)
c
      character title*80,str*80,casename*256,rootname*256,
     &   blank*80
      character*22 dtgi,dtgfz,dtgfl
c
      blank=' '
c
c   Write title at top of plot.
c
      call set(0.,1.,0.,1.,0.,1.,0.,1.,1)
      call createdtg(mdateb,rhourb,xtime,timezone,
     &   iusdaylightrule,inearesth,dtgi,dtgfz,dtgfl)
c
      if (title(1:10).eq.'none      '.or.idotitle.eq.0) then
         str=' '
      elseif (title(1:10).eq.'auto      ') then
         ibegc=1
         ibegcr=1
         do i=256,2,-1
            if (ibegc.eq.1.and.casename(i-1:i-1).eq.'/') ibegc=i
            if (ibegcr.eq.1.and.rootname(i-1:i-1).eq.'/') ibegcr=i
         enddo
         ilenc=min(iendc-ibegc+1,15)
         ilencr=31-ilenc
         str='Dataset: '//casename(ibegc:ibegc+ilenc-1)//
     &       '  RIP: '//rootname(ibegcr:ibegcr+ilencr-1)
      else
         str=title
      endif
c
      ifirstpart=0
      if (str.ne.blank) then
         chsize=.012
         ypos=toptextclg-.5*chsize
         call plchhq(.005,ypos,str(1:47),chsize,0.,-1.)
         toptextclg=toptextclg-1.9*chsize
         ifirstpart=1
      endif
c
      if (iinittime.eq.1) then
         str=' '
c         if (iprog.eq.6.or.iprog.eq.11) then
c            write(str,'(a6,a22)') 'Init: ',dtgi
c            iendst=28
c         elseif (iprog.eq.2.or.iprog.eq.3.or.iprog.eq.5) then
c            write(str,'(a7,a22)') 'Valid: ',dtgfz
c            iendst=29
c         endif
         write(str,'(a6,a22)') 'Init: ',dtgi
         iendst=28
         if (ifirstpart.eq.0) then
            chsize=.012
            ypos=toptextclg-.5*chsize
         endif
         if (inearesth.eq.1) iendst=iendst-2
         call plchhq(.995,ypos,str(1:iendst),chsize,0.,1.)
         if (ifirstpart.eq.0) toptextclg=toptextclg-1.9*chsize
      endif
c
c  Now, the valid time line, if asked for
c
c      if (ivalidtime.eq.1.and.(iprog.eq.6.or.iprog.eq.11)) then
      if (ivalidtime.eq.1) then
         str=' '
         ypos=toptextclg-.006
         if (iprog.eq.6.or.iprog.eq.11) then
            write(str,'(a6)') 'Fcst: '
            iendst=6
         else
            write(str,'(a4)') 'T + '
            iendst=4
         endif
         if (inearesth.eq.0) then
            if (xtime.ge.0..or.str(1:4).eq.'Fcst') then
               write(str(iendst+1:),'(f6.2,a2)') xtime,' h'
            else
               str(3:3)='-'
               write(str(iendst+1:),'(f6.2,a2)') -xtime,' h'
            endif
         else
            ixtime=nint(xtime)
            if (ixtime.ge.0.or.str(1:4).eq.'Fcst') then
               write(str(iendst+1:),'(i3,a2)') ixtime,' h'
            else
               str(3:3)='-'
               write(str(iendst+1:),'(i3,a2)') -ixtime,' h'
            endif
         endif
         call plchhq(.005,ypos,str(1:13),.012,0.,-1.)
         if (inearesth.eq.0) then
            write(str,'(a7,a22,a2,a22,a1)')
     &         'Valid: ',dtgfz,' (',dtgfl,')'
            call plchhq(.995,ypos,str(1:54),.012,0.,1.)
         else
            write(str,'(a7,a20,a2,a20,a1)')
     &         'Valid: ',dtgfz(1:20),' (',dtgfl(1:20),')'
            call plchhq(.995,ypos,str(1:50),.012,0.,1.)
         endif
         toptextclg=toptextclg-.021
      endif
c
      return
      end
