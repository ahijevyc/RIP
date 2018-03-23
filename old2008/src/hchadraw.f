c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine hchadraw(ixwin,iywin,rlevl,work,icdwk,ilev,
     &         iintv,icosq,rcosq,lchfl,icolr,incsq,pslab1,maxcosq,
     &         mabpl,morpl,maxlev,maxpl,miy,mjx,mkzh,ipl,rrota)
c
      dimension ixwin(2,maxpl),iywin(2,maxpl),
     &   icosq(maxcosq,maxpl),rlevl(maxlev,maxpl),iintv(maxpl),
     &   work(miy,mjx,mkzh),icdwk(maxpl),
     &   rcosq(maxcosq,maxpl),incsq(maxpl),icolr(maxpl),
     &   pslab1(mabpl,morpl),rrota(maxpl),temp(mabpl,morpl),
     &   temp2(morpl,mabpl)
      logical lchfl(maxpl)
      real ch,cv
c
      dimension xra(4),yra(4),dst(10),ind(10)
      character pchar*1, pchar2*2
c
      include 'comconst'
c
c   Make proper set call, as well as number of values to
c      plot in each direction
c
      xintervs=ixwin(2,ipl)-ixwin(1,ipl)
      yintervs=iywin(2,ipl)-iywin(1,ipl)
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
      if (icdwk(ipl).eq.0) then
         ul=1.
         ur=1.+xintervs
         ub=1.
         ut=1.+yintervs
         niy=nint(ut)
         njx=nint(ur)
      else
         ul=.5
         ur=.5+xintervs
         ub=.5
         ut=.5+yintervs
         niy=int(ut)
         njx=int(ur)
      endif
c MGD begin mod
c just switch around corners when calling set() if rotating by (-)90
c because we have to assume grid is a rectangle
      if(rrota(ipl) .ne. 90. .and. rrota(ipl) .ne. -90.) then
        call set(fl,fr,fb,ft,ul,ur,ub,ut,1)
      else
c MGD end mod
        cv = (ft-fb)/2.
        ch = (fr-fl)/2.
        call set((fr+fl)/2.-cv,(fr+fl)/2.+cv,
     &           (ft+fb)/2.-ch,(ft+fb)/2.+ch,ub,ut,ul,ur,1)
      endif
c
c   Put appropriate data into horizontal slab.
c
      do 90 j=1,njx
         jj=j+ixwin(1,ipl)-1
         do 90 i=1,niy
            ii=i+iywin(1,ipl)-1
            pslab1(j,i)=work(ii,jj,nint(rlevl(ilev,ipl)))
   90 continue

c MGD begin mod
c transpose/invert the data grids... 
c for ccw rotation, x becomes y and inverted y becomes x
      if(rrota(ipl) .eq. 90.) then
        do i=1,njx
        do j=1,niy
          temp2(niy-j+1,i) = pslab1(i,j)
        enddo
        enddo
c for cw rotation, inverted x becomes y and y becomes x
      elseif(rrota(ipl) .eq. -90.) then
        do i=1,njx
        do j=1,niy
          temp2(j,njx-i+1) = pslab1(i,j)
        enddo
        enddo
c for 180 rotation, invert both x and y
      elseif(rrota(ipl) .eq. -180. .or. rrota(ipl) .eq. 180.) then
        do i=1,njx
        do j=1,niy
          temp(njx-i+1,niy-j+1) = pslab1(i,j)
        enddo
        enddo
      endif

c if 180 degree rotation, dimensions of array are the same, so use
c the original one
      if(rrota(ipl) .ne. 90. .and. rrota(ipl) .ne. -90. .and. rrota(ipl)
     &   .ne. 0.) then
      do i=1,njx
      do j=1,niy
        pslab1(i,j) = temp(i,j)
      enddo
      enddo
      endif
c MGD end mod

c
c   Plot the characters
c
      setfrac=(fr-fl)/(ur-ul)*50.
      if (iintv(ipl).gt.0) then
         chsiz=max(.004,.012*setfrac*iintv(ipl))
      else
         chsiz=max(.004,.012*setfrac*iintv(ipl)*.5)
      endif
      chsiz2=max(.004,.7*chsiz)
c
      if (.not.lchfl(ipl)) then
c
      intp=iintv(ipl)
      if (intp.gt.0) then  ! regular staggering
         intph=0
      else   ! "diamond pattern" staggering
         intp=-intp
         intph=intp/2
      endif
c MGD mod
      if(rrota(ipl) .ne. 90. .and. rrota(ipl) .ne. -90.) then
      do j=1,njx
      do i=1,niy
         if ( ((mod(i-1,intp).eq.0.and.mod(j-1,intp).eq.0).or.
     &         (mod(i-1+intph,intp).eq.0.and.
     &          mod(j-1+intph,intp).eq.0)) .and.
     &       pslab1(j,i).ne.rmsg) then
            ips=nint(min(pslab1(j,i),100.4))
            if (ips.gt.99.or.ips.lt.0) goto 100
            rj=j
            ri=i
            do 95 ir = 1,incsq(ipl)
               if (nint(rcosq(ir,ipl)).eq.ips) then
                  call gstxci(icosq(ir,ipl))
                  call gsplci(icosq(ir,ipl))
                  goto 96
               endif
   95       continue
            call gstxci(icolr(ipl))
            call gsplci(icolr(ipl))
   96       continue        
            if (ips.lt.10) then
               write(pchar,'(i1)') ips
               call plchhq(rj,ri,pchar,chsiz,0.,0.)
            else
               write(pchar2,'(i2)') ips
               call plchhq(rj,ri,pchar2,chsiz2,0.,0.)
            endif
         endif
  100    continue
      enddo
      enddo
c MGD begin mod
c basically the same as before, except using array with reversed
c dimensions... (temp2 instead of pslab1)
      else
      do j=1,niy
      do i=1,njx
         if ( ((mod(i-1,intp).eq.0.and.mod(j-1,intp).eq.0).or.
     &         (mod(i-1+intph,intp).eq.0.and.
     &          mod(j-1+intph,intp).eq.0)) .and.
     &       temp2(j,i).ne.rmsg) then
            ips=nint(min(temp2(j,i),100.4))
            if (ips.gt.99.or.ips.lt.0) goto 200
            rj=j
            ri=i
            do 195 ir = 1,incsq(ipl)
               if (nint(rcosq(ir,ipl)).eq.ips) then
                  call gstxci(icosq(ir,ipl))
                  call gsplci(icosq(ir,ipl))
                  goto 196
               endif
  195       continue
            call gstxci(icolr(ipl))
            call gsplci(icolr(ipl))
  196       continue        
            if (ips.lt.10) then
               write(pchar,'(i1)') ips
               call plchhq(rj,ri,pchar,chsiz,0.,0.)
            else
               write(pchar2,'(i2)') ips
               call plchhq(rj,ri,pchar2,chsiz2,0.,0.)
            endif
         endif
  200    continue
      enddo
      enddo
      endif
c MGD end mod
c
      else
c
      if(rrota(ipl) .ne. 90. .and. rrota(ipl) .ne. -90.) then
      do 120 j=1,njx
      do 120 i=1,niy
         ips=nint(min(pslab1(j,i),100.4))
         if (ips.gt.99.or.ips.lt.0) goto 120
         rj=j
         ri=i
         do 105 ir = 1,incsq(ipl)
            if (nint(rcosq(ir,ipl)).eq.ips) then
               ici=icosq(ir,ipl)
               if (ici.eq.-1) go to 120 ! Transparent
               goto 106
            endif
  105    continue
         ici=1
  106    continue        
         xra(1)=rj-.5
         xra(2)=rj+.5
         xra(3)=rj+.5
         xra(4)=rj-.5
         yra(1)=ri-.5
         yra(2)=ri-.5
         yra(3)=ri+.5
         yra(4)=ri+.5
         call sfsgfa(xra,yra,4,dst,10,ind,10,ici)
  120 continue
c MGD begin mod
c same story - (-)90 rotations use differently dimensioned array
c temp2(j,i) instead of pslab(i,j) 
      else
      do 220 j=1,niy
      do 220 i=1,njx
         ips=nint(min(temp2(j,i),100.4))
         if (ips.gt.99.or.ips.lt.0) goto 220
         rj=j
         ri=i
         do 205 ir = 1,incsq(ipl)
            if (nint(rcosq(ir,ipl)).eq.ips) then
               ici=icosq(ir,ipl)
               if (ici.eq.-1) go to 220 ! Transparent
               goto 206
            endif
  205    continue
         ici=1
  206    continue        
         xra(1)=rj-.5
         xra(2)=rj+.5
         xra(3)=rj+.5
         xra(4)=rj-.5
         yra(1)=ri-.5
         yra(2)=ri-.5
         yra(3)=ri+.5
         yra(4)=ri+.5
         call sfsgfa(xra,yra,4,dst,10,ind,10,ici)
  220 continue
      endif
c MGD end mod
c
      endif
c
      call gstxci(1)
      call gsplci(1)
      return
      end
