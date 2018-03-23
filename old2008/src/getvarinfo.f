c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine getvarinfo(casename,iendc,cxtimeavl,xtimeavl,nxt,
     &   ncxc,varname,maxtavl,ndim,icd,vardesc,plchun,
     &   iabort,istat,iup)
c
      dimension xtimeavl(maxtavl)
      character fname*256,varname*10,cxtimeavl(maxtavl)*10,
     &   casename*(*)
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
c     If you are attempting to re-code rip to read data directly from
c     native model output file(s), this would be the appropriate place
c     to first try to obtain info on the requested variable from that file(s).
c     However, if that fails, you should still have the routine attempt
c     to get info on the variable from a rip data file, in case it was
c     previously saved to such a file by rip.
c
      fname=casename(1:iendc)//'_'//
     &      cxtimeavl(nxt)(1:ncxc)//'_'//varname
      open(unit=25,err=30,file=fname,form='unformatted',
     &   status='old')
      read (25)
     &   vardesc,plchun,ihrip,rhrip,chrip
      ndim=ihrip(6)
      icd=ihrip(7)
      close (25)
      istat=1
      return
c
 30   istat=-1
      if (iabort.eq.1) then
         write(iup,*)'Fatal: couldn''t find file named'
         write(iup,*)'   ',fname
         stop
      endif
      return
      end
