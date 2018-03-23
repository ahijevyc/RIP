c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine readbasic(uuu,vvv,tmk,qvp,www,prs,ght,prs_tsf,ter,
     &   dmap,xmap,cor,xlus,pstx,pstd,sigh,sigf,ptop,inhyd,iprog,
     &   fname,iendf1,iudata,miy,mjx,mkzh)
c
      dimension uuu(miy,mjx,mkzh),vvv(miy,mjx,mkzh),tmk(miy,mjx,mkzh),
     &   qvp(miy,mjx,mkzh),www(miy,mjx,mkzh),prs(miy,mjx,mkzh),
     &   ght(miy,mjx,mkzh),prs_tsf(miy,mjx),ter(miy,mjx),
     &   dmap(miy,mjx),xmap(miy,mjx),cor(miy,mjx),xlus(miy,mjx),
     &   pstx(miy,mjx),pstd(miy,mjx),sigf(mkzh+1),sigh(mkzh)
      character fname*256
c
c   Read files that should always be there, regardless of iprog
c
      call readdat(iudata,fname,iendf1,'ter       ',
     &   miy,mjx,mkzh,2,1,ter,istat)
      call readdat(iudata,fname,iendf1,'dmap      ',
     &   miy,mjx,mkzh,2,1,dmap,istat)
      call readdat(iudata,fname,iendf1,'xmap      ',
     &   miy,mjx,mkzh,2,1,xmap,istat)
      call readdat(iudata,fname,iendf1,'cor       ',
     &   miy,mjx,mkzh,2,1,cor,istat)
      call readdat(iudata,fname,iendf1,'xlus      ',
     &   miy,mjx,mkzh,2,1,xlus,istat)
      call readdat(iudata,fname,iendf1,'pstx      ',
     &   miy,mjx,mkzh,2,1,pstx,istat)
      call readdat(iudata,fname,iendf1,'pstd      ',
     &   miy,mjx,mkzh,2,1,pstd,istat)
c
      if (iprog.eq.1) then
c
c   Fill other arrays with reasonable info if this is TERRAIN data
c
      do j=1,mjx
      do i=1,miy
         uuu(i,j,1)=0.
         vvv(i,j,1)=0.
         tmk(i,j,1)=273.
         qvp(i,j,1)=0.
         www(i,j,1)=0.
         prs(i,j,1)=sigh(1)*pstx(i,j)+ptop
         ght(i,j,1)=5000.
      enddo
      enddo
c
      else
c
c   Read files that will only be there for iprog > 1
c   (i.e. anything other than TERRAIN output)
c
      call readdat(iudata,fname,iendf1,'uuu       ',
     &   miy,mjx,mkzh,3,1,uuu,istat)
      call readdat(iudata,fname,iendf1,'vvv       ',
     &   miy,mjx,mkzh,3,1,vvv,istat)
      call readdat(iudata,fname,iendf1,'tmk       ',
     &   miy,mjx,mkzh,3,1,tmk,istat)
      call readdat(iudata,fname,iendf1,'qvp       ',
     &   miy,mjx,mkzh,3,0,qvp,istat)
      if (istat.eq.-1) then
         call fillarray(qvp,miy*mjx*mkzh,0.)
      else ! convert from g/kg to kg/kg
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            qvp(i,j,k)=.001*qvp(i,j,k)
         enddo
         enddo
         enddo
      endif
      if (iprog.ge.5) then
         call readdat(iudata,fname,iendf1,'prs       ',
     &      miy,mjx,mkzh,3,1,prs,istat)
         do j=1,mjx-1
         do i=1,miy-1
            prs_tsf(i,j)=prs(i,j,mkzh)+(1.-sigh(mkzh))*pstx(i,j)
         enddo
         enddo
      elseif (iprog.eq.2.or.iprog.eq.3) then
         do k=1,mkzh
            pressure=sigh(k)*pstx(1,1)+ptop
            do j=1,mjx-1
            do i=1,miy-1
               prs(i,j,k)=pressure
            enddo
            enddo
         enddo
         call readdat(iudata,fname,iendf1,'prs_tsf   ',
     &      miy,mjx,mkzh,2,1,prs_tsf,istat)
      endif
      call readdat(iudata,fname,iendf1,'ght       ',
     &   miy,mjx,mkzh,3,1,ght,istat)
      if (inhyd.eq.1) then
         call readdat(iudata,fname,iendf1,'www       ',
     &      miy,mjx,mkzh,3,1,www,istat)
      else
         call vervcalc(pstd,dmap,sigf,xmap,pstx,sigh,
     &      uuu,vvv,qvp,tmk,www,prs,www,1,miy,mjx,mkzh)
      endif
c
      endif
c
      return
      end
