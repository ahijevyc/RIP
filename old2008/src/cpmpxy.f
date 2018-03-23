c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine cpmpxy(imap,xinp,yinp,xotp,yotp)
c
c   This is a user-supplied subroutine that transforms each grid
c   point's vertical level index to the chosen vertical coordinate.
c
      include 'comconst'
      include 'comvctran'
c
      if (imap.ne.0) then
         reml=xinp-int(xinp)
         remk=yinp-int(yinp)
         l=max(1,int(xinp))
         k=max(1,int(yinp))
         lpl=min(int(xinp)+1,nscross)
         kpl=min(int(yinp)+1,mkzhcross)
         yotp=(   reml)*(   remk)*vc2d(lpl,kpl)+
     &      (1.-reml)*(   remk)*vc2d(l,kpl)+
     &      (   reml)*(1.-remk)*vc2d(lpl,k)+
     &      (1.-reml)*(1.-remk)*vc2d(l,k)
         xotp=xinp
      else
         yotp=yinp
         xotp=xinp
      endif
      return
      end
