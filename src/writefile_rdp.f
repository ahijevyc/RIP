c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine writefile_rdp (var,varname,ndim,icd,vardesc,plchun,
     &   fname,iendf1,ihrip,rhrip,chrip,
     &   iexpanded,iexpandedout,ioffexp,joffexp,miy,mjx,mkzh)
c
      dimension var(miy,mjx,1+(mkzh-1)*(ndim-2))
      character varname*10,fname*256
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32)
      character chrip(64)*64,vardesc*64,plchun*24
c
c ensure data do not exceed ieee limits
c
      kk = 1+(mkzh-1)*(ndim-2)
      do k=1,kk
         do j=1,mjx-icd
         do i=1,miy-icd
           if(var(i,j,k) .gt. 3.4e+38) var(i,j,k)= 3.4e+38
           if(var(i,j,k) .lt.-3.4e+38) var(i,j,k)=-3.4e+38
           if(abs(var(i,j,k)).lt. 1.2e-38) var(i,j,k)= 0.
         enddo
         enddo
         if (icd.eq.1) then
            do i=1,miy
               var(i,mjx,k)=0.
            enddo
            do j=1,mjx-1
               var(miy,j,k)=0.
            enddo
         endif
      enddo
      iendv=index(varname,' ')-1
      if (iendv.eq.-1) iendc=10
      fname(iendf1+1:)=varname
      open(unit=65,file=fname,form='unformatted',status='unknown')
      ihrip(6)=ndim ! number of dimensions of this variable (2 or 3)
      ihrip(7)=icd  ! grid of this var. (1: cross point, 0: dot point)
      write(65) vardesc,plchun,ihrip,rhrip,chrip
      if (iexpanded.eq.1.and.iexpandedout.eq.0) then
         write(65) (((var(i,j,k),i=1+ioffexp,miy-ioffexp),
     &      j=1+joffexp,mjx-joffexp),k=1,1+(mkzh-1)*(ndim-2))
      else
         write(65) var
      endif
      close (65)
      return
      end
