c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hmapdraw(ngfbuf,ixwin,ixwingf,iywin,iywingf,
     &         yicorngf,xjcorngf,icolr,icolrgf,ilinw,ilinwgf,
     &         idash,idashgf,rtslb,rtslbgf,rcint,rcintgf,cmllm,
     &         cmllmgf,couty,coutygf,couds,coudsgf,ioulw,ioulwgf,
     &         iouco,ioucogf,imfco,imfcogf,iwhatgf,iam,xcs,ycs,
     &         rip_root,niam,ncs,maxbuf,maxpl,ipl,rrota,rrotagf)
c
      dimension ixwin(2,maxpl),ixwingf(2,maxbuf),iywin(2,maxpl),
     &   iywingf(2,maxbuf),yicorngf(maxbuf),xjcorngf(maxbuf),
     &   icolr(maxpl),icolrgf(maxbuf),ilinw(maxpl),
     &   ilinwgf(maxbuf),idash(maxpl),idashgf(maxbuf),
     &   rtslb(maxpl),rtslbgf(maxbuf),rcint(maxpl),
     &   rcintgf(maxbuf),ioulw(maxpl),ioulwgf(maxbuf),
     &   iouco(maxpl),ioucogf(maxbuf),imfco(6,maxpl),
     &   imfcogf(6,maxbuf),iwhatgf(maxbuf), rrotagf(maxbuf)
      dimension iam(niam),xcs(ncs),ycs(ncs),rrota(maxpl)
      dimension icol(5),icsf(5)
c
      character chproj*2,cmllm(maxpl)*5,cmllmgf(maxbuf)*5,
     &   couty(maxpl)*32,coutygf(maxbuf)*32,couds(maxpl)*5,
     &   coudsgf(maxbuf)*5,rip_root*256
c
      dimension iai(2),iag(2)
      real llclat,llclon,urclat,urclon
c
      parameter (nama=300000, ncra=50000)
      dimension iama(nama),xcra(ncra),ycra(ncra)
      common /mapspc/iama,xcra,ycra
c
      common /emap/ llcolor,lllinw,llndot,mpfillco(6),llmask,ioutype
c
      external colram
      external colrln
c
      include 'comconst'
c
c   Check if this map is identical to any of those saved in the
c      GFLASH buffers, and if so, flush that buffer to this frame.
c
      do 10 ib=1,ngfbuf
         if (iwhatgf(ib).ne.1) goto 10
         if ( ixwin(1,ipl).eq.ixwingf(1,ib).and.    ! x-window
     &        ixwin(2,ipl).eq.ixwingf(2,ib).and.    ! x-window
     &        iywin(1,ipl).eq.iywingf(1,ib).and.    ! y-window
     &        iywin(2,ipl).eq.iywingf(2,ib).and.    ! y-window
     &        yicorn.eq.yicorngf(ib).and.           ! y-corner point in crs dom
     &        xjcorn.eq.xjcorngf(ib).and.           ! x-corner point in crs dom
     &        icolr(ipl).eq.icolrgf(ib).and.  ! color of lat/lon labels, lines
     &        ilinw(ipl).eq.ilinwgf(ib).and.  ! line width of lat/lon lines
     &        idash(ipl).eq.idashgf(ib).and.  ! dash pattern of lat/lon lines
     &        rtslb(ipl).eq.rtslbgf(ib).and.  ! text size of lat/lon labels
     &        rcint(ipl).eq.rcintgf(ib).and.  ! interval of lat/lon lines
     &        cmllm(ipl).eq.cmllmgf(ib).and.  ! masking of lat/lon lines
     &        couty(ipl).eq.coutygf(ib).and.  ! map outline data set type
     &        couds(ipl).eq.coudsgf(ib).and.  ! dot or solid map outlines
     &        ioulw(ipl).eq.ioulwgf(ib).and.  ! line width of solid map outlines
     &        iouco(ipl).eq.ioucogf(ib).and.  ! color of map outlines
     &        imfco(1,ipl).eq.imfcogf(1,ib).and.  !  map fill color 1,2,3,etc.
     &        imfco(2,ipl).eq.imfcogf(2,ib).and.
     &        imfco(3,ipl).eq.imfcogf(3,ib).and.
     &        imfco(4,ipl).eq.imfcogf(4,ib).and.
     &        imfco(5,ipl).eq.imfcogf(5,ib).and.
     &        imfco(6,ipl).eq.imfcogf(6,ib).and.
     &        rrota(ipl).eq.rrotagf(ib)) then     ! rotation angle for pol st
            call gflas3(ib)
            goto 20
         endif
   10 continue
c
c   If identical map was not found in GFLASH buffer, draw a new one.
c
      ngfbuf=ngfbuf+1
      if (ngfbuf.gt.maxbuf) then
         write(iup,*)'In hmapdraw: too many gflash buffers.'
         stop
      endif
      ib=ngfbuf
      call gflas1(ib)
c
c   This should save unrotated window lat and lon for shifting
c   the window after projection is rotated about origin
c   (projection is rotated, window is not, then window is shifted to 
c   new position to achieve a rotated map outline)
c
      yllc=yicorn+(iywin(1,ipl)-1.)/refrat
      xllc=xjcorn+(ixwin(1,ipl)-1.)/refrat
      call maptform(yllc,xllc,llclat,llclon,1,0.)
      call maptform(yllc,xllc,rlatllc,rlonllc,1,rrota(ipl))
      yurc=yicorn+(iywin(2,ipl)-1.)/refrat
      xurc=xjcorn+(ixwin(2,ipl)-1.)/refrat
      call maptform(yurc,xurc,urclat,urclon,1,0.)
      call maptform(yurc,xurc,rlaturc,rlonurc,1,rrota(ipl))
c
c   Change xlonc, but since LC projections dont work, this is irrelevant
c
      if(nproj .eq. 1) then
        xlonc=xlonc-rrota(ipl)
        plon=xlonc
      else
        plon=xlonc
      endif
c
      if (nproj.eq.1) then
         chproj='LC'
         plat=true1
         rota=true2
c        plat=90.
c        rota=plat
      elseif (nproj.eq.2) then
         chproj='ST'
         plat=sign(90.,xlatc)
c
c   This is the goodie that rotates the outline by rrota
c
         rota=-1*rrota(ipl)
      elseif (nproj.eq.3) then
         chproj='ME'
         plat=0.
         rota=0.
      endif
c
c   Start ezmap.
c
      call mappos (flmin,frmax,fbmin,ftmax)
      call maproj (chproj,plat,plon,rota)
      call mapset ('CO',llclat,llclon,urclat,urclon)
      call mapint
c
c   Various initial stuff
c
      call mapsti ('MV',2)  ! minimum vector, in plotting units (out of 32768)
      call mapsti ('VS',5)  ! number of vertical strip divisions for area map
      call mapsti ('LA',0)  ! label meridians and poles?
      call mapsti ('PE',0)  ! draw perimeter?
      call gsplci(iouco(ipl))
      call gspmci(iouco(ipl))
      call gstxci(iouco(ipl))
      call getusv ('LW',lwsv)
      if (couds(ipl).eq.'dot  ') then
         if (couty(ipl)(1:3).eq.'NO '.or.couty(ipl)(1:3).eq.'CO '.or.
     &       couty(ipl)(1:3).eq.'US '.or.couty(ipl)(1:3).eq.'PS '.or.
     &       couty(ipl)(1:3).eq.'PO '.or.couty(ipl)(1:7).eq.'Earth..'
     &      ) then
            call mapsti ('DO',1)
            call mapsti ('DD',ioulw(ipl))
         else
            write(iup,*) '   Warning: dotted outlines only work with'
            write(iup,*) '   outline choices NO, CO, US, PS, PO,'
            write(iup,*) '   and Earth..n, not with RANGS/GSHSS (RG)'
            write(iup,*) '   or custom map outline files.'
         endif
      elseif (couds(ipl).eq.'solid') then
         call mapsti ('DO',0)
         call setusv ('LW',ioulw(ipl)*1000)
      else
         write(iup,*)'In hmapdraw, dot/solid specifier is bad.'
         write(iup,*)'ipl,couds=',ipl,'  #',couds(ipl),'#'
         stop
      endif
c
c   Define area filling and masking parameters for common blocks
c
      if (couty(ipl)(1:3).eq.'NO '.or.couty(ipl)(1:3).eq.'CO '.or.
     &    couty(ipl)(1:3).eq.'US '.or.couty(ipl)(1:3).eq.'PS '.or.
     &    couty(ipl)(1:3).eq.'PO '.or.couty(ipl)(1:7).eq.'Earth..'.or.
     &    couty(ipl)(1:2).eq.'RG') then
         do i=1,6
            mpfillco(i)=imfco(i,ipl)
         enddo
      else
         if (imfco(1,ipl).ne.999999) then
            write(iup,*) '   Warning: area filling only works with'
            write(iup,*) '   outline choices NO, CO, US, PS, PO,'
            write(iup,*) '   Earth..n, and RANGS/GSHSS (RG), not'
            write(iup,*) '   with custom map outline files.'
         endif
         imfco(1,ipl)=999999
      endif
      if (couty(ipl)(1:3).eq.'NO '.or.couty(ipl)(1:3).eq.'CO '.or.
     &    couty(ipl)(1:3).eq.'US '.or.couty(ipl)(1:3).eq.'PS '.or.
     &    couty(ipl)(1:3).eq.'PO '.or.couty(ipl)(1:7).eq.'Earth..') then
         if (cmllm(ipl).eq.'none ') then
            llmask=0
         elseif (cmllm(ipl).eq.'land ') then
            llmask=1
         elseif (cmllm(ipl).eq.'water') then
            llmask=-1
         else
            write(iup,*)'In hmapdraw, mask specifier is bad.'
            write(iup,*)'ipl,cmllm=',ipl,cmllm(ipl)
            stop
         endif
      else
         if (cmllm(ipl).ne.'none ') then
            write(iup,*) '   Warning: masking of lat/lon lines by land'
            write(iup,*) '   or water only works with outline choices'
            write(iup,*) '   NO, CO, US, PS, PO, and Earth..n,'
            write(iup,*) '   not with RANGS/GSHSS (RG) or custom'
            write(iup,*) '   map outline files.'
         endif
         llmask=0
      endif
c
c   Do things that are specific to outline type, and plot map
c
      if (couty(ipl)(1:3).eq.'NO '.or.couty(ipl)(1:3).eq.'CO '.or.
     &    couty(ipl)(1:3).eq.'US '.or.couty(ipl)(1:3).eq.'PS '.or.
     &    couty(ipl)(1:3).eq.'PO ') then  ! traditional outline data sets
         call mapstc ('OU',couty(ipl)(1:2))
         ioutype=1
         if (mpfillco(1).ne.999999.or.llmask.ne.0) then
            call arinam (iam,niam) ! Initialize the area map.
            call mapbla (iam) ! Add edges to the area map.
            call arpram (iam,0,0,0) ! Pre-process the area map.
            if (mpfillco(1).ne.999999) then   ! color the map
               call arscam (iam,xcs,ycs,ncs,iai,iag,2,colram)
            endif
         endif
         call plotit(0,0,0)  ! Flush plotit's buffers
         call maplot
      elseif (couty(ipl)(1:7).eq.'Earth..') then  ! Earth.. outline data sets
         if  (couty(ipl)(9:9).eq.'L') then
            read(couty(ipl)(10:10),'(i1)') iearthlev
         else
            iearthlev=4
         endif
         ioutype=2
         if (mpfillco(1).ne.999999.or.llmask.ne.0) then
            call arinam (iam,niam) ! Initialize the area map.
            call mplnam (couty(ipl)(1:8),iearthlev,iam) ! Add edges to area map.
            call arpram (iam,0,0,0) ! Pre-process the area map.
            if (mpfillco(1).ne.999999) then   ! color the map
               call arscam (iam,xcs,ycs,ncs,iai,iag,2,colram)
            endif
         endif
         call plotit(0,0,0)  ! Flush plotit's buffers
         call mplndr(couty(ipl)(1:8),iearthlev)
      elseif (couty(ipl)(1:2).eq.'RG') then  ! RANGS/GSHHS outline data sets
         write(iup,*) 'To use the RG data set you must have'
         write(iup,*) 'NCAR Graphics Version 4.3.  If you have it,'
         write(iup,*) 'comment out this print statement and the'
         write(iup,*) 'subsequent stop statement in routine hmapdraw.,'
         write(iup,*) 'and uncomment the following section of code.'
         write(iup,*) 'Then recompile and rerun rip.'
         stop
c         if (couty(ipl)(3:3).ne.' ') then
c            read(couty(ipl)(3:3),'(i1)') irgres  ! 0 (finest) to 4 (coarsest)
c         else
c            call mdrgdl(irgres) ! Get suggested resolution of RANGS/GSHHS data.
c         endif
c         call mpseti('II',64064)
c         if (mpfillco(1).ne.999999) then
c            do i=1,5
c               icol(i)=iouco(ipl)
c               icsf(i)=mpfillco(i)
c            enddo
c            call mdrgsc(icol,icsf)
c            call mdrgsf (irgres,xcra,ncra,iama,nama)
c            call mdrgol (irgres,xcra,ncra)
c         else
c            do i=1,5
c               icol(i)=iouco(ipl)
c               icsf(i)=iouco(ipl)
c            enddo
c            call mdrgsc(icol,icsf)
c            call mdrgol (irgres,xcra,ncra)
c         endif
      else  ! Custom outline data sets
         call hiresmap(rip_root,couty(ipl),iutserstn)
      endif
c
      call setusv ('LW',lwsv)
c
c   Lat/lon labels
c
      if (rtslb(ipl).gt.0..and.rcint(ipl).gt.0.) then
         call gsplci(icolr(ipl))
         call gspmci(icolr(ipl))
         call gstxci(icolr(ipl))
         if(nproj .eq. 1 .or. rrota(ipl) .eq. 0.) then
         call maptick(refrat,yicorn,xjcorn,iywin(1,ipl),ixwin(1,ipl),
     &      nint(rcint(ipl)),2,rtslb(ipl),iup,rrota(ipl))
         endif
      endif
c
c   Lat/lon lines:
c
      if (rcint(ipl).gt.0.0) then
         call mapsti ('GR',nint(rcint(ipl)))
         call getdash(idash(ipl),llndot)
         call getdash(10,llndotch)
         llcolor=icolr(ipl)
         lllinw=ilinw(ipl)
         if (llmask.ne.0) then
            call mapgrm (iam,xcs,ycs,ncs,iai,iag,2,colrln)
         else
            call mapgrd
         endif
      endif
c
      call gflas2
      call gflas3(ib)
c
      iywingf(1,ib)=iywin(1,ipl)
      iywingf(2,ib)=iywin(2,ipl)
      ixwingf(1,ib)=ixwin(1,ipl)
      ixwingf(2,ib)=ixwin(2,ipl)
      yicorngf(ib)=yicorn
      xjcorngf(ib)=xjcorn
      icolrgf(ib)=icolr(ipl)
      ilinwgf(ib)=ilinw(ipl)
      idashgf(ib)=idash(ipl)
      rtslbgf(ib)=rtslb(ipl)
      rcintgf(ib)=rcint(ipl)
      cmllmgf(ib)=cmllm(ipl)
      coutygf(ib)=couty(ipl)
      coudsgf(ib)=couds(ipl)
      ioulwgf(ib)=ioulw(ipl)
      ioucogf(ib)=iouco(ipl)
      imfcogf(1,ib)=imfco(1,ipl)
      imfcogf(2,ib)=imfco(2,ipl)
      imfcogf(3,ib)=imfco(3,ipl)
      imfcogf(4,ib)=imfco(4,ipl)
      imfcogf(5,ib)=imfco(5,ipl)
      imfcogf(6,ib)=imfco(6,ipl)
      rrotagf(ib)=rrota(ipl)
      iwhatgf(ib)=1
c
   20 continue
c
c   Set color back to white
c
      call gsplci (1)
      call gspmci (1)
      call gstxci (1)
      call gsfaci (1)
      return
      end
