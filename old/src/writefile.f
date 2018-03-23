c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine writefile (var,varname,numpas,ndim,icd,vardesc,
     &   plchun,fname,iendf1,ihrip,rhrip,chrip,fullsigma,halfsigma,
     &   miy,mjx,mkzh)
c
      dimension var(miy,mjx,1+(mkzh-1)*(ndim-2))
      character varname*10,fname*256,temp*4
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32),fullsigma(128),halfsigma(128)
      character chrip(64)*64,vardesc*64,plchun*24
c
      include 'comconst'
c
      iendv=index(varname,' ')-1
      if (iendv.eq.-1) iendv=10
      if (numpas.gt.0) then
         write(temp,'(i4)') 1000+numpas
         fname(iendf1+1:)=varname(1:iendv)//temp(2:4)
      else
         fname(iendf1+1:)=varname(1:iendv)
      endif
      open(unit=iudata,file=fname,form='unformatted',status='unknown')
      ihrip(6)=ndim ! number of dimensions of this variable (2 or 3)
      ihrip(7)=icd  ! grid of this var. (1:cross point, 0:dot point)
      write(iudata) vardesc,plchun,ihrip,rhrip,chrip,
     &   fullsigma,halfsigma
      write(iudata) var
      close (iudata)
      return
      end
