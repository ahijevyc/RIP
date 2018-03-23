      subroutine netasc(prs1,prs2,sigf,sigh,casename,iendc,
     &   ctjfl,rtjst,rtjen,asc,miy,mjx,mkzh)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                      c
c   This subroutine reads in a gridswarm trajectory position file and  c
c   calculates the neta ascent for a specified time period at all      c
c   domain grid points that are within the gridswarm grid.             c

c The program assumes:
c
c 1) The entire trajectory file is one 3D gridswarm, i.e. trajectory
c    initial points lie in a rectangular grid.
c 2) The gridswarm points lie on model cross points.  This means the
c    initial trajectory x and y locations should all be an
c    integer + 0.5.
c 3) The grid has a horizontal grid spacing that is some integer
c    multiple of the model grid space.  Interpolation will be used
c    to get net ascent at the model points in between the trajectory
c    grid points.
c 4) The gridswarm vertical levels must be half sigma values of
c    the model grid, though the order does not matter, and not
c    all levels need to be included.  Interpolation will be used
c    to get net ascent at the model levels in between the trajectory
c    levels.
c 5) The trajectories are arranged such that the x-index of the
c    gridswarm varies most rapidly, then the y, and then the z.
c
      parameter (maxtraj=1000,maxtrajtime=220)
      dimension prs1(miy,mjx,mkzh),prs2(miy,mjx,mkzh),
     &   asc(miy,mjx,mkzh),sigf(mkzh+1),sigh(mkzh)
      dimension stortr(maxtrajtime,maxtraj,3),asctr(maxtraj),
     &   sight(300),kgot(300)
      character fname*256,casename*256,chtime*10,ctjfl*256
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32),fullsigma(128),halfsigma(128)
      character chrip(64)*64,vardesc*64,plchun*24
c
      include 'comconst'
c
c   Read in data from trajectory file
c
      open (unit=iutrajin,file=ctjfl,form='unformatted',status='old')
      read (iutrajin)
      read (iutrajin) rtim,ctim,dttraj,ntraj
      ntrajtime=nint(abs(rtim-ctim)/dttraj*3600) + 1
      if (rtim.lt.ctim) then
         itm1=1
         itm2=ntrajtime
         itmi=1
      else
         itm1=ntrajtime
         itm2=1
         itmi=-1
      endif
      if (rtjst.eq.rmsg) then
         starttim=min(rtim,ctim)
      else
         starttim=max(rtjst,min(rtim,ctim))
      endif
      if (rtjen.eq.rmsg) then
         endtim=max(rtim,ctim)
      else
         endtim=min(rtjen,max(rtim,ctim))
      endif
      do itm=itm1,itm2,itmi
         read(iutrajin) (stortr(itm,itr,1),itr=1,ntraj),
     &       (stortr(itm,itr,2),itr=1,ntraj),
     &       (stortr(itm,itr,3),itr=1,ntraj)
      enddo
      close (iutrajin)
      itma1=1+nint((starttim-min(rtim,ctim))/dttraj*3600)
      itma2=itma1+nint((endtim-starttim)/dttraj*3600)
c
c   Read in pressure data at earliest time in net ascent window.
c
c      write(iup,*)'starttim=',starttim
c      write(iup,*)'1000.+starttim=',1000.+starttim
      write(chtime,'(f10.5)') 1000.+starttim
      fname=casename(1:iendc)//'_'//chtime(2:10)//'_'//'prs'
c      write(iup,*)'fname=',fname,'   iudata=',iudata
      open(unit=iudata,file=fname,form='unformatted',status='old')
      read(iudata,err=170,end=180)
     &   vardesc,plchun,ihrip,rhrip,chrip,fullsigma,halfsigma
      il=ihrip(4)
      jl=ihrip(5)
      kl=ihrip(9)
      if (il.ne.miy.or.jl.ne.mjx.or.kl.ne.mkzh) then
         write(iup,*)'Parameters miy,mjx,mkzh must be set'
         write(iup,*)'to domain dimensions of ',il,jl,kl
      endif
      read(iudata,err=190,end=195)
     &   (((prs1(i,j,k),i=1,miy),j=1,mjx),k=1,mkzh)
      close (iudata)
c
c   Read in pressure data at latest time in net ascent window.
c
      write(chtime,'(f10.5)') 1000.+endtim
      fname=casename(1:iendc)//'_'//chtime(2:10)//'_'//'prs'
      open(unit=iudata,file=fname,form='unformatted',status='unknown')
      read(iudata,err=170,end=180)
     &   vardesc,plchun,ihrip,rhrip,chrip,fullsigma,halfsigma
      il=ihrip(4)
      jl=ihrip(5)
      kl=ihrip(9)
      if (il.ne.miy.or.jl.ne.mjx.or.kl.ne.mkzh) then
         write(iup,*)'Parameters miy,mjx,mkzh must be set'
         write(iup,*)'to domain dimensions of ',il,jl,kl
      endif
      read(iudata,err=190,end=195)
     &   (((prs2(i,j,k),i=1,miy),j=1,mjx),k=1,mkzh)
      close (iudata)
c
c   Get information from rip header record
c
      sigf(mkzh+1)=fullsigma(mkzh+1)
      do k=1,mkzh
         sigh(k)=halfsigma(k)
         sigf(k)=fullsigma(k)
      enddo
      miycors=ihrip(2)
      mjxcors=ihrip(3)
      inhyd=ihrip(8)
      dskmc=rhrip(5)
      dskm=rhrip(6)
      yicorn=rhrip(7)
      xjcorn=rhrip(8)
      ptop=rhrip(9)
      refslp=rhrip(10)
      refslt=rhrip(11)
      reflaps=rhrip(12)
      refstratt=rhrip(16)
      dsc=dskmc*1000.
      ds=dskm*1000.
      refrat=dsc/ds
c
c   Get pressure at trjectory locations.
c   Also, convert x and y values (which are in coarse dom. grid points)
c   to grid points in the current domain.
c
c      write(iup,*)'rtim,ctim,starttim,endtim=',
c     &   rtim,ctim,starttim,endtim
c      write(iup,*)'itm1,itm2,itmi,itma1,itma2=',
c     &   itm1,itm2,itmi,itma1,itma2
      do itr=1,ntraj
         if (stortr(itm1,itr,1).eq.rmsg.or.
     &       stortr(itm1,itr,2).eq.rmsg.or.
     &       stortr(itm1,itr,3).eq.rmsg) then
            write(iup,*)'problem with position of a',
     &         ' trajectory at rtim.'
            write(iup,*)'itr,st1,st2,st3=',itr,stortr(itm1,itr,1),
     &         stortr(itm1,itr,2),stortr(itm1,itr,3)
            stop
         elseif (itm1.ne.itma1.and.itm1.ne.itma2) then
            stortr(itm1,itr,1)=
     &         1.+refrat*(stortr(itm1,itr,1)-yicorn)
            stortr(itm1,itr,2)=
     &         1.+refrat*(stortr(itm1,itr,2)-xjcorn)
         endif
         if (stortr(itma1,itr,1).ne.rmsg) then
            pr1=finterp(prs1,1,sigh,miy,mjx,mkzh,
     &         stortr(itma1,itr,1),stortr(itma1,itr,2),
     &         stortr(itma1,itr,3),refrat,yicorn,xjcorn,rmsg,iup)
            if (pr1.eq.rmsg) then
               stortr(itma1,itr,1)=rmsg
               stortr(itma1,itr,2)=rmsg
               stortr(itma1,itr,3)=rmsg
            else
               stortr(itma1,itr,1)=
     &            1.+refrat*(stortr(itma1,itr,1)-yicorn)
               stortr(itma1,itr,2)=
     &            1.+refrat*(stortr(itma1,itr,2)-xjcorn)
            endif
         endif
         if (stortr(itma2,itr,1).ne.rmsg) then
            pr2=finterp(prs2,1,sigh,miy,mjx,mkzh,
     &         stortr(itma2,itr,1),stortr(itma2,itr,2),
     &         stortr(itma2,itr,3),refrat,yicorn,xjcorn,rmsg,iup)
            if (pr2.eq.rmsg) then
               stortr(itma2,itr,1)=rmsg
               stortr(itma2,itr,2)=rmsg
               stortr(itma2,itr,3)=rmsg
            else
               stortr(itma2,itr,1)=
     &            1.+refrat*(stortr(itma2,itr,1)-yicorn)
               stortr(itma2,itr,2)=
     &            1.+refrat*(stortr(itma2,itr,2)-xjcorn)
            endif
         endif
         if (stortr(itma1,itr,1).ne.rmsg.and.
     &       stortr(itma2,itr,1).ne.rmsg) then
            asctr(itr)=pr2-pr1
         else
            asctr(itr)=rmsg
         endif
      enddo
c
c   Figure out grid space and dimensions of gridswarm.
c
c      do itr=1,150
c         write(iup,*)'itr,xval,yval=',stortr(itm1,itr,2)
c      enddo
      iytest=abs(nint(100.*(stortr(itm1,1,1)-nint(stortr(itm1,1,1)))))
      jxtest=abs(nint(100.*(stortr(itm1,1,2)-nint(stortr(itm1,1,2)))))
      if (iytest.ne.50.or.jxtest.ne.50) then
         write(iup,*)'In netasc: trj. swarm corner not on cross point.'
         write(iup,*)'x,y vals. = ',stortr(itm1,1,2),stortr(itm1,1,1)
         write(iup,*)'jxtest,iytest=',jxtest,iytest
         stop
      endif
      iorientest=nint(100.*(stortr(itm1,2,1)-stortr(itm1,1,1)))
      if (iorientest.ne.0) then
         write(iup,*)'In netasc: trj. swarm not properly oriented.'
         write(iup,*)'First two y vals. = ',stortr(itm1,1,1),
     &      stortr(itm1,2,1)
         stop
      endif
      ids=nint(stortr(itm1,2,2)-stortr(itm1,1,2))
      iycorn=nint(stortr(itm1,1,1)-.5)
      jxcorn=nint(stortr(itm1,1,2)-.5)
      igotnx=0
      do itr=2,ntraj
         if (igotnx.eq.0.and.
     &       stortr(itm1,itr,2).lt.stortr(itm1,itr-1,2)) then
            igotnx=1
            nx=itr-1
         endif
         if (stortr(itm1,itr,1).lt.stortr(itm1,itr-1,1)-.01) then
            ny=(itr-1)/nx !prob
            goto 33
         endif
      enddo
 33   nkzh=ntraj/(nx*ny)
c      write(iup,*)'Grid swarm: ids,nx,ny,nkzh,iycorn,jxcorn=',
c     &   ids,nx,ny,nkzh,iycorn,jxcorn
      do k=1,nkzh
         sight(k)=stortr(itm1,1+(k-1)*nx*ny,3)
c         write(iup,*)'   k,sight=',k,sight(k)
      enddo
c
c   Fill the model array with the values from the gridswarm.
c      
      do k=1,mkzh
         kgot(k)=0
         do j=1,mjx-1
         do i=1,miy-1
            asc(i,j,k)=rmsg
         enddo
         enddo
      enddo
c
      itr=0
      ibeg=iycorn
      iend=ibeg+(ny-1)*ids
      jbeg=jxcorn
      jend=jbeg+(nx-1)*ids
      do kt=1,nkzh
         do k=1,mkzh
            if (abs(sight(kt)-sigh(k)).lt..00001) then
               kgot(k)=1
               goto 37
            endif
         enddo
         write(iup,*)'Couldn''t find matching value in model sigma',
     &      ' array.'
         write(iup,*)'kt,sight=',kt,sight(kt)
         do k=1,mkzh
            write(iup,*)'k,sigh=',k,sigh(k)
         enddo
         stop
 37      continue
         do i=ibeg,iend,ids
         do j=jbeg,jend,ids
            itr=itr+1
            asc(i,j,k)=asctr(itr)
         enddo
         enddo
      enddo
c
c   Fill in the missing values in between
c
      do k=1,mkzh
         if (kgot(k).eq.1) then
            do j=jbeg,jend-ids,ids
            do i=ibeg,iend-ids,ids
               if (asc(i,j,k).ne.rmsg.and.
     &             asc(i+ids,j,k).ne.rmsg.and.
     &             asc(i,j+ids,k).ne.rmsg.and.
     &             asc(i+ids,j+ids,k).ne.rmsg) then
                  do jj=j,j+ids
                  do ii=i,i+ids
                     rj=(jj-j)/float(ids)
                     ri=(ii-i)/float(ids)
                     asc(ii,jj,k)=( ri  )*( rj  )*asc(i+ids,j+ids,k)+
     &                            (1.-ri)*( rj  )*asc(i    ,j+ids,k)+
     &                            ( ri  )*(1.-rj)*asc(i+ids,j    ,k)+
     &                            (1.-ri)*(1.-rj)*asc(i    ,j    ,k)
                  enddo
                  enddo
               endif
            enddo
            enddo
         endif
      enddo
c
      do j=jbeg,jend
      do i=ibeg,iend
         do k=2,mkzh-1
            if (asc(i,j,k).eq.rmsg) then
               kabove=-1
               kbelow=-1
               do kk=k-1,1,-1
                  if (asc(i,j,kk).ne.rmsg) then
                     kabove=kk
                     goto 73
                  endif
               enddo
 73            continue
               do kk=k+1,mkzh
                  if (asc(i,j,kk).ne.rmsg) then
                     kbelow=kk
                     goto 75
                  endif
               enddo
 75            continue
               if (kabove.ne.-1.and.kbelow.ne.-1) then
                  denom=sigh(kbelow)-sigh(kabove)
                  r1=(sigh(k)-sigh(kabove))/denom
                  r2=(sigh(kbelow)-sigh(k))/denom
                  asc(i,j,k)=r1*asc(i,j,kbelow)+r2*asc(i,j,kabove)
               endif
            endif
         enddo
      enddo
      enddo
c
      return
c
 170  write(iup,*)'In netasc: The model data header is not a format'//
     &   ' that RIP recognizes.  Stopping.'
      stop
c
 180  write(iup,*)'In netasc: Unexpected EOF reached when trying'
      write(iup,*)'to read model data header.  Stopping.'
      stop
c
 190  write(iup,*)'In netasc: Error in reading the model data array.'//
     &   ' Stopping.'
      stop
c
 195  write(iup,*)'In netasc: Unexpected EOF encountered while'
      write(iup,*)'reading the model data array. Stopping.'
      stop
c
      end
