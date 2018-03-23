c                                                                     c
c*********************************************************************c
c                                                                     c
      real function finterp2d(arr,icd,miy,mjx,yi,xj,
     &   refrat,yicorn,xjcorn,rmsg)
c
c   Note, this function assumes that yi and xj are with respect to the
c   coarse domain dot-point grid.  The source array can be either
c   dot or cross (icd must be set accordingly).
c
      dimension arr(miy,mjx)
c
      xjthisdom=1.+(xj-xjcorn)*refrat
      yithisdom=1.+(yi-yicorn)*refrat
      if ((icd.eq.0.and.
     &      (xjthisdom.lt.1..or.xjthisdom.gt.float(mjx).or.
     &       yithisdom.lt.1..or.yithisdom.gt.float(miy)))   .or.
     &    (icd.eq.1.and.
     &      (xjthisdom.lt.1.5.or.xjthisdom.gt.float(mjx)-.5.or.
     &       yithisdom.lt.1.5.or.yithisdom.gt.float(miy)-.5))) then
         finterp2d=rmsg
         return
      endif
      jarrg=int(xjthisdom-icd*.5)
      jarrgp=jarrg+1
      iarrg=int(yithisdom-icd*.5)
      iarrgp=iarrg+1
      ratx=xjthisdom-icd*.5-jarrg
      raty=yithisdom-icd*.5-iarrg
      if (arr(iarrgp,jarrg).ne.rmsg.and.
     +    arr(iarrgp,jarrgp).ne.rmsg.and.
     +    arr(iarrg,jarrg).ne.rmsg.and.
     +    arr(iarrg,jarrgp).ne.rmsg) then
         finterp2d=(1.-ratx)*(   raty)*arr(iarrgp,jarrg)+
     +           (   ratx)*(   raty)*arr(iarrgp,jarrgp)+
     +           (1.-ratx)*(1.-raty)*arr(iarrg,jarrg)+
     +           (   ratx)*(1.-raty)*arr(iarrg,jarrgp)
      else
         finterp2d=rmsg
      endif
      return
      end
