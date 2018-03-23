      subroutine getheadinfo(ihrip,rhrip,fullsigma,halfsigma,chrip,
     &   vardesc,plchun,sigf,sigh,nproj,miycors,mjxcors,inhyd,
     &   mdateb,mhourb,iice,iprog,ilandset,iwater,truelat1,truelat2,
     &   xlatc,xlonc,dskmc,dskm,yicorn,xjcorn,ptop,
     &   refslp,refslt,reflaps,refstratt,rhourb,dsc,ds,refrat,mkzh)
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32),fullsigma(128),halfsigma(128)
      character chrip(64)*64,vardesc*64,plchun*24
c
      dimension sigf(mkzh+1),sigh(mkzh)

      sigf(mkzh+1)=fullsigma(mkzh+1)
      do k=1,mkzh
         sigh(k)=halfsigma(k)
         sigf(k)=fullsigma(k)
      enddo
      nproj=ihrip(1)
      miycors=ihrip(2)
      mjxcors=ihrip(3)
      inhyd=ihrip(8)
      mdateb=ihrip(10)
      call mconvert(mdateb,mhourb,1,1940)
      iice=ihrip(12)
      iprog=ihrip(13)
      if (iprog.lt.1.or.iprog.gt.50) iprog=6
      ilandset=ihrip(14)
      if (ilandset.lt.1.or.ilandset.gt.10) ilandset=1
      if (ilandset.eq.1) then
         iwater=7
      elseif (ilandset.eq.2) then
         iwater=16
      elseif (ilandset.eq.3) then
         iwater=15
      endif
      truelat1=rhrip(1)
      truelat2=rhrip(2)
      xlatc=rhrip(3)
      xlonc=rhrip(4)
      dskmc=rhrip(5)
      dskm=rhrip(6)
      yicorn=rhrip(7)
      xjcorn=rhrip(8)
      ptop=rhrip(9)
      refslp=rhrip(10)
      refslt=rhrip(11)
      reflaps=rhrip(12)
      refstratt=rhrip(16)
      if (refstratt.gt.500.or.refstratt.lt.0.) refstratt=0.1
      rhourb=rhrip(13)
      dsc=dskmc*1000.
      ds=dskm*1000.
      refrat=dsc/ds
c
      return
      end
