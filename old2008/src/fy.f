c                                                                     c
c*********************************************************************c
c                                                                     c
      function fy(x,y)
c
      include 'comconst'
      include 'comvctran'
c
      if (ivcs.eq.1) then
         reml=x-int(x)
         remk=y-int(y)
         l=max(1,int(x))
         k=max(1,int(y))
         lpl=min(int(x)+1,nscross)
         kpl=min(int(y)+1,mkzhcross)
         fy=(   reml)*(   remk)*vc2d(lpl,kpl)+
     &      (1.-reml)*(   remk)*vc2d(l,kpl)+
     &      (   reml)*(1.-remk)*vc2d(lpl,k)+
     &      (1.-reml)*(1.-remk)*vc2d(l,k)
      else
         fy=y
      endif
      return
      end
