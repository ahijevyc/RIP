      program ripinterp
c
c   This program reads in model output (in rip-format files)
c   from a coarse domain and from a fine domain, and creates a new file
c   which has the data from the coarse domain file interpolated
c   (bi-linearly) to the fine domain.  The header and data dimensions
c   of the new file will be that of the fine domain, and the case name
c   used in the file name will be the same as that of the fine domain
c   file that was read in.
c
      character argum(16)*256
c
c   RIP header variables
c
      dimension ihrip1(32),rhrip1(32)
      character chrip1(64)*64,vardesc1*64,plchun1*24
      dimension ihrip2(32),rhrip2(32)
      character chrip2(64)*64,vardesc2*64,plchun2*24
c
      dimension arr1(200,200,40),arr2(200,200,40)
c
c   Get command line arguments.
c
      nargum=iargc()
      if (nargum.ne.3) then
         print*,
     &      'Usage: ripinterp coarse_file fine_file new_fine_file'
         stop
      endif
      do i=1,nargum
         call getarg(i,argum(i))
      enddo
c
c   Fix for machines (such as HP) that return the command name itself
c   as the first element of argum, rather than the first argument.
c
      if (argum(1)(1:10).eq.'ripinterp_') then
         do i=1,nargum-1
            argum(i)=argum(i+1)
         enddo
         nargum=nargum-1
      endif
c
c   Get information from file name and from header record for domain 1
c   (the coarser domain)
c
      ixtpos=index(argum(1),'_')+1
      read(argum(1)(ixtpos:ixtpos+8),'(f9.5)')xtime1
      open (unit=11,file=argum(1),form='unformatted',status='old')
      read(11,err=170,end=180)
     &   vardesc1,plchun1,ihrip1,rhrip1,chrip1
      goto 190
 170  print*,'Error in reading file 1.'
      stop
 180  print*,'Unexpected EOF encountered in reading file 1.'
      stop
 190  continue
      il1=ihrip1(4)
      jl1=ihrip1(5)
      kl1=ihrip1(9)
      nproj1=ihrip1(1)
      miycors1=ihrip1(2)
      mjxcors1=ihrip1(3)
      ndim1=ihrip1(6)
      icd1=ihrip1(7)
      mdateb1=ihrip1(10)
      mdate1=ihrip1(11)
      iice1=ihrip1(12)
      truelat11=rhrip1(1)
      truelat21=rhrip1(2)
      xlatc1=rhrip1(3)
      xlonc1=rhrip1(4)
      dskmc1=rhrip1(5)
      dskm1=rhrip1(6)
      yicorn1=rhrip1(7)
      xjcorn1=rhrip1(8)
      rhourb1=rhrip1(13)
c
c   Get information from file name and from header record for domain 2
c   (the finer domain)
c
      ixtpos=index(argum(2),'_')+1
      read(argum(2)(ixtpos:ixtpos+8),'(f9.5)')xtime2
      open (unit=12,file=argum(2),form='unformatted',status='old')
      read(12,err=270,end=280)
     &   vardesc2,plchun2,ihrip2,rhrip2,chrip2
      goto 290
 270  print*,'Error in reading file 2.'
      stop
 280  print*,'Unexpected EOF encountered in reading file 2.'
      stop
 290  continue
      il2=ihrip2(4)
      jl2=ihrip2(5)
      kl2=ihrip2(9)
      nproj2=ihrip2(1)
      miycors2=ihrip2(2)
      mjxcors2=ihrip2(3)
      ndim2=ihrip2(6)
      icd2=ihrip2(7)
      mdateb2=ihrip2(10)
      mdate2=ihrip2(11)
      iice2=ihrip2(12)
      truelat12=rhrip2(1)
      truelat22=rhrip2(2)
      xlatc2=rhrip2(3)
      xlonc2=rhrip2(4)
      dskmc2=rhrip2(5)
      dskm2=rhrip2(6)
      yicorn2=rhrip2(7)
      xjcorn2=rhrip2(8)
      rhourb2=rhrip2(13)
c
c   Check compatibility of domains
c
      ier=0
      if (kl2.ne.kl1) then
         print*,'Error: kl1,kl2=',kl1,kl2
         ier=ier+1
      endif
      if (nproj2.ne.nproj1) then
         print*,'nproj1,nproj2=',nproj1,nproj2
         ier=ier+1
      endif
      if (miycors2.ne.miycors1) then
         print*,'miycors1,miycors2=',miycors1,miycors2
         ier=ier+1
      endif
      if (mjxcors2.ne.mjxcors1) then
         print*,'mjxcors1,mjxcors2=',mjxcors1,mjxcors2
         ier=ier+1
      endif
      if (ndim2.ne.ndim1) then
         print*,'ndim1,ndim2=',ndim1,ndim2
         ier=ier+1
      endif
      if (icd2.ne.icd1) then
         print*,'icd1,icd2=',icd1,icd2
         ier=ier+1
      endif
      icd=icd1
      if (mdate2.ne.mdate1) then
         write(6,'(''mdate1,mdate2='',2(i8.8,1x))')mdate1,mdate2
         ier=ier+1
      endif
      if (truelat12.ne.truelat11) then
         print*,'truelat11,truelat12=',truelat11,truelat12
         ier=ier+1
      endif
      if (truelat22.ne.truelat21) then
         print*,'truelat21,truelat22=',truelat21,truelat22
         ier=ier+1
      endif
      if (xlatc2.ne.xlatc1) then
         print*,'xlatc1,xlatc2=',xlatc1,xlatc2
         ier=ier+1
      endif
      if (xlonc2.ne.xlonc1) then
         print*,'xlonc1,xlonc2=',xlonc1,xlonc2
         ier=ier+1
      endif
      if (dskmc2.ne.dskmc1) then
         print*,'dskmc1,dskmc2=',dskmc1,dskmc2
         ier=ier+1
      endif
      if (ier.gt.0) stop
      irefrat=nint(1000.*dskm1/dskm2)
      if (irefrat.ne.1000.and.irefrat.ne.3000.and.
     &    irefrat.ne.5000.and.irefrat.ne.7000.and.
     &    irefrat.ne.9000.and.irefrat.ne.11000) then
         print*,'Refinement ratio must be an odd number from 1 to 11.'
         stop
      endif
      irefrat=irefrat/1000
c
c  Determine position of lower left corner of domain 2 in domain 1.
c
      irefrat1=nint(dskmc1/dskm1)
      yicorn=(yicorn2-yicorn1)*irefrat1+1
      xjcorn=(xjcorn2-xjcorn1)*irefrat1+1
      ychk=abs(mod(irefrat*yicorn,1.))
      xchk=abs(mod(irefrat*xjcorn,1.))
      if (.not.((ychk.lt..01.or.ychk.gt..99).and.
     &          (xchk.lt..01.or.xchk.gt..99))) then
         print*,'Domains are not grid-point coincident.'
         stop
      endif
      print*,'dskmc1,dskm1,irefrat1,yicorn1,xjcorn1,yicorn2,xjcorn2='
      print*, dskmc1,dskm1,irefrat1,yicorn1,xjcorn1,yicorn2,xjcorn2
      print*,'Refinement ratio, yicorn, xjcorn = ',
     &         irefrat,yicorn,xjcorn
c
c   Read the data
c
      klread=1+(ndim1-2)*(kl1-1)
      read(11) (((arr1(i,j,k),i=1,il1),j=1,jl1),k=1,klread)
      read(12) (((arr2(i,j,k),i=1,il2),j=1,jl2),k=1,klread)
c
c   Start interpolation loop
c
      do k=1,klread
      do j2=1,jl2-icd
         rj1=(j2-1.+.5*icd)/irefrat+xjcorn-.5*icd
         j1m=int(rj1)
         j1p=j1m+1
         facj=rj1-j1m
      do i2=1,il2-icd
         ri1=(i2-1.+.5*icd)/irefrat+yicorn-.5*icd
         i1m=int(ri1)
         i1p=i1m+1
         faci=ri1-i1m
         arr2(i2,j2,k)=(1.-facj)*(1.-faci)*arr1(i1m,j1m,k)+
     &                 (   facj)*(1.-faci)*arr1(i1m,j1p,k)+
     &                 (1.-facj)*(   faci)*arr1(i1p,j1m,k)+
     &                 (   facj)*(   faci)*arr1(i1p,j1p,k)
      enddo
      enddo
      enddo
c
c   Write the new file
c
      open (unit=21,file=argum(3),form='unformatted',status='unknown')
      write(21)
     &   vardesc1,plchun1,ihrip2,rhrip2,chrip2
      write(21) (((arr2(i,j,k),i=1,il2),j=1,jl2),k=1,klread)
c
      stop
      end
