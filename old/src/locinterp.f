c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine locinterp(string,gridx,gridy,rlat,rlon,iwmo,tlc,
     &   rip_root,rrota)
c
      character string(2)*20,tlc*4,fnm*256,threelc*4,rip_root*256
c
      logical numeric
      external numeric
c
      include 'comconst'
c
      if (string(1).eq.'missing             '.and.
     &    string(2).eq.'missing             ') then ! Everyth. mssng.
         write(iup,*)'In locinterp, location specification is'
         write(iup,*)'completely missing.  Stopping.'
         stop
      elseif (string(2).eq.'missing             ') then ! WMO or 3LC
         if (numeric(string(1)(1:5))) then ! WMO
            read(string(1)(1:7),'(i5)') iwmo
            tlc='XXXX'
         else ! 3LC
            read(string(1)(1:4),'(a4)') tlc
	    if(string(1)(4:4).eq. ' ') tlc = 'K'//string(1)(1:3) 
            iwmo=99999
         endif
         iendci=index(rip_root,' ')-1
         fnm=rip_root(1:iendci)//'/stationlist'
         open (unit=iustnlist,file=fnm,form='formatted',status='old')
         read(iustnlist,*)
         read(iustnlist,*)
  200    read(iustnlist,'(23x,a4,5x,f6.2,f8.2,7x,i5)',end=205)
     &      threelc,rlat,rlon,nwmo
	 if (threelc(1:1) .eq. ' ') threelc(1:1) = 'K'
         if (iwmo.eq.nwmo.or.tlc.eq.threelc) then
            call maptform(gridy,gridx,rlat,rlon,-1,rrota)
            if (nwmo.ne.00000) iwmo=nwmo
            if (threelc.ne.'K---') tlc=threelc
            goto 210
         endif
         goto 200
  205    write(iup,*)'Couldn''t find station requested for'
         write(iup,*)'   iwmo=',iwmo,'  and tlc=#',tlc,'#'
         stop
  210    continue
         close (iustnlist)
      else ! Either x/y or lat/lon specification.
         iwmo=99999
         tlc='XXXX'
         inlat1=index(string(1),'lat')
         inlon1=index(string(1),'lon')
         inlat2=index(string(2),'lat')
         inlon2=index(string(2),'lon')
         if (inlat1.ne.0.or.inlon1.ne.0) then ! Lat/lon specification.
            if (inlat1.ne.0) then ! Lat is first.
               read(string(1)(1:inlat1-1),fmt=*) rlat
               read(string(2)(1:inlon2-1),fmt=*) rlon
            else ! Lon is first.
               read(string(1)(1:inlon1-1),fmt=*) rlon
               read(string(2)(1:inlat2-1),fmt=*) rlat
            endif
            call maptform(gridy,gridx,rlat,rlon,-1,rrota)
         else ! x/y specification
            read(string(1)(1:12),fmt=*) gridxn
            read(string(2)(1:12),fmt=*) gridyn
            gridx=(gridxn-1.)/refrat+xjcorn
            gridy=(gridyn-1.)/refrat+yicorn
            call maptform(gridy,gridx,rlat,rlon,1,rrota)
         endif
      endif
      return
      end
