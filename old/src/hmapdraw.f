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
c
      character chproj*2,cmllm(maxpl)*5,cmllmgf(maxbuf)*5,
     &   couty(maxpl)*32,coutygf(maxbuf)*32,couds(maxpl)*5,
     &   coudsgf(maxbuf)*5,rip_root*256
c
      dimension iai(2),iag(2)
      character fnm*256
      real llclat,llclon,urclat,urclon
c
      common /emap/ llcolor,lllinw,llndot,mpfillco(6),llmask
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
         if ( ixwin(1,ipl).eq.ixwingf(1,ib).and.ixwin(2,ipl).eq.
     &        ixwingf(2,ib).and.iywin(1,ipl).eq.iywingf(1,ib).and.
     &        iywin(2,ipl).eq.iywingf(2,ib).and.
     &        yicorn.eq.yicorngf(ib).and.
     &        xjcorn.eq.xjcorngf(ib).and.
     &        icolr(ipl).eq.icolrgf(ib).and.
     &        ilinw(ipl).eq.ilinwgf(ib).and.
     &        idash(ipl).eq.idashgf(ib).and.
     &        rtslb(ipl).eq.rtslbgf(ib).and.
     &        rcint(ipl).eq.rcintgf(ib).and.
     &        cmllm(ipl).eq.cmllmgf(ib).and.
     &        couty(ipl).eq.coutygf(ib).and.
     &        couds(ipl).eq.coudsgf(ib).and.
     &        ioulw(ipl).eq.ioulwgf(ib).and.
     &        iouco(ipl).eq.ioucogf(ib).and.
     &        imfco(1,ipl).eq.imfcogf(1,ib).and.
     &        imfco(2,ipl).eq.imfcogf(2,ib).and.
     &        imfco(3,ipl).eq.imfcogf(3,ib).and.
     &        imfco(4,ipl).eq.imfcogf(4,ib).and.
     &        imfco(5,ipl).eq.imfcogf(5,ib).and.
     &        imfco(6,ipl).eq.imfcogf(6,ib)      ) then
c MGD begin mod
c don't keep old map outline if rotations are different
c maps aren't identical if their rotations are different...
         if(rrota(ipl) .eq. rrotagf(ib)) then
            call gflas3(ib)
            goto 20
         endif
c MGD end mod
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
c MGD begin mod 
c this should save unrotated window lat and lon for shifting
c the window after projection is rotated about origin
c (projection is rotated, window is not, then window is shifted to 
c  new position to achieve a rotated map outline)
      yllc=yicorn+(iywin(1,ipl)-1.)/refrat
      xllc=xjcorn+(ixwin(1,ipl)-1.)/refrat
      call maptform(yllc,xllc,llclat,llclon,1,0.)
      call maptform(yllc,xllc,rlatllc,rlonllc,1,rrota(ipl))
      yurc=yicorn+(iywin(2,ipl)-1.)/refrat
      xurc=xjcorn+(ixwin(2,ipl)-1.)/refrat
      call maptform(yurc,xurc,urclat,urclon,1,0.)
      call maptform(yurc,xurc,rlaturc,rlonurc,1,rrota(ipl))
c
c change xlonc, but since LC projections dont work, this is irrelevant
      if(nproj .eq. 1) then
        xlonc=xlonc-rrota(ipl)
        plon=xlonc
      else
        plon=xlonc
      endif
c MGD end mod
      if (nproj.eq.1) then
         chproj='LC'
         plat=true1
         rota=true2
c        plat=90.
c        rota=plat
      elseif (nproj.eq.2) then
         chproj='ST'
         plat=sign(90.,xlatc)
c MGD this is the goodie that rotates the outline by rrota
         rota=-1*rrota(ipl)
      elseif (nproj.eq.3) then
         chproj='ME'
         plat=0.
         rota=0.
      endif
c
c   Do various "mapst" calls
c
      ihires=0
      if (couty(ipl)(1:3).ne.'NO '.and.couty(ipl)(1:3).ne.'CO '.and.
     &    couty(ipl)(1:3).ne.'US '.and.couty(ipl)(1:3).ne.'PS '.and.
     &    couty(ipl)(1:3).ne.'PO '.and.couty(ipl)(1:3).ne.'M1 '.and.
     &    couty(ipl)(1:3).ne.'M3 ') then
         ihires=1
      endif
      if (ihires.eq.0) then
         call mapstc ('OU',couty(ipl)(1:2))
      else
         call mapstc ('OU','PS')
      endif
      call mapsti ('MV',2)
      call mapsti ('G1',1)
      call mapsti ('G2',2)
      call mapsti ('VS',5)
      call mapsti ('LA',0)
      call mapsti ('PE',0)
c
c   Start ezmap.
c
      call mappos (flmin,frmax,fbmin,ftmax)
      call maproj (chproj,plat,plon,rota)
      call mapset ('CO',llclat,llclon,urclat,urclon)
      call mapint
c
      if ((imfco(1,ipl).ne.999999.or.cmllm(ipl).ne.'none ').and.
     &    (couty(ipl)(1:3).eq.'NO '.or.couty(ipl)(1:3).eq.'CO '.or.
     &     couty(ipl)(1:3).eq.'US '.or.couty(ipl)(1:3).eq.'PS '.or.
     &     couty(ipl)(1:3).eq.'PO ')) then
         call arinam (iam,niam) ! Initialize the area map.
         call mapbla (iam) ! Add edges to the area map.
         call arpram (iam,0,0,0) ! Pre-process the area map.
      elseif (imfco(1,ipl).ne.999999.or.cmllm(ipl).ne.'none ') then
         write(iup,*) '   Warning: area filling and lat/lon'
         write(iup,*) '   masking only work with standard NCAR'
         write(iup,*) '   Graphics outline choices (i.e., with'
         write(iup,*) '   outy=NO, CO, US, PS, or PO).'
      endif
c
      if (imfco(1,ipl).ne.999999.and.
     &    (couty(ipl)(1:3).eq.'NO '.or.couty(ipl)(1:3).eq.'CO '.or.
     &     couty(ipl)(1:3).eq.'US '.or.couty(ipl)(1:3).eq.'PS '.or.
     &     couty(ipl)(1:3).eq.'PO ')) then
c
c      Color the map.
c
         do i=1,6
            mpfillco(i)=imfco(i,ipl)
         enddo
         call arscam (iam,xcs,ycs,ncs,iai,iag,2,colram)
c
      endif
c
      call plotit(0,0,0)  ! Flush plotit's buffers
c
c   Draw the map outlines.
c
      call gsplci(iouco(ipl))
      call gspmci(iouco(ipl))
      call gstxci(iouco(ipl))
      call getusv ('LW',lwsv)
      if (couds(ipl).eq.'dot  ') then
         call mapsti ('DO',1)
         call mapsti ('DD',ioulw(ipl))
      elseif (couds(ipl).eq.'solid') then
         call mapsti ('DO',0)
         call setusv ('LW',ioulw(ipl)*1000)
      else
         write(iup,*)'In hmapdraw, dot/solid specifier is bad.'
         write(iup,*)'ipl,couds=',ipl,'  #',couds(ipl),'#'
         stop
      endif
      if (ihires.eq.1) then
c
c      Use HIRES map routine "borrowed" from GRAPH
c
        iendci=index(rip_root,' ')-1
        fnm=rip_root(1:iendci)//'/'//couty(ipl)
        iendfnm=index(fnm,' ')-1
        open (unit=83,file=fnm,form='formatted',status='old',err=33)
	call hiresmap(83)
        close (83)
        goto 34
 33     continue
        write(iup,*)'Could not find the hires map data file'
        write(iup,*)'specified by rip_root and outy.'
        write(iup,*)'Looked for file "',fnm(1:iendfnm),'"'
 34     continue
      else
        call maplot
c       CALL MPLNDR ('Earth..1',3)   ! ncarg 4.1.0  higher resolution map
      endif
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
         if (llmask.ne.0.and.
     &       (couty(ipl)(1:3).eq.'NO '.or.couty(ipl)(1:3).eq.'CO '.or.
     &        couty(ipl)(1:3).eq.'US '.or.couty(ipl)(1:3).eq.'PS '.or.
     &        couty(ipl)(1:3).eq.'PO ')) then
            call mapgrm (iam,xcs,ycs,ncs,iai,iag,2,colrln)
         else
            call mapgrd
         endif
      endif
c
c     call supcon (40.583,-105.083, xjim, yjim)
c     call pwritx(xjim, yjim, 'FCL',3,12, 0, 1)
c     call supcon (40.01,-105.25, xjim, yjim)
c     call pwritx(xjim, yjim, 'BOU',3,12, 0, 1)
c     call supcon (32.84,-106.14, xjim, yjim)
c     call pwritx(xjim, yjim, 'HMN',3,12, 0, 1)
c     call supcon (35.04,-106.62, xjim, yjim)
c     call pwritx(xjim, yjim, 'ABQ',3,12, 0, 1)

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
c MGD mod 
c for buffering...
      rrotagf(ib)=rrota(ipl)
      iwhatgf(ib)=1
c      write(iup,*)'   Made new map background number ',ib,
c     &   ' for plot ',ipl
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
