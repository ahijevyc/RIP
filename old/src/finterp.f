c                                                                     c
c*********************************************************************c
c                                                                     c
      real function finterp(arr,icd,sigh,miy,mjx,mkzh,yi,xj,zk,
     &   refrat,yicorn,xjcorn,rmsg,iup)
c
c   Note, this function assumes that yi and xj are with respect to the
c   coarse domain dot-point grid.  The source array can be either
c   dot or cross (icd must be set accordingly).
c
      dimension arr(miy,mjx,mkzh),sigh(mkzh)
c
      xjthisdom=1.+(xj-xjcorn)*refrat
      yithisdom=1.+(yi-yicorn)*refrat
      if ((icd.eq.0.and.
     &      (xjthisdom.lt.1..or.xjthisdom.gt.float(mjx).or.
     &       yithisdom.lt.1..or.yithisdom.gt.float(miy)))   .or.
     &    (icd.eq.1.and.
     &      (xjthisdom.lt.1.5.or.xjthisdom.gt.float(mjx)-.5.or.
     &       yithisdom.lt.1.5.or.yithisdom.gt.float(miy)-.5))) then
         finterp=rmsg
         return
      endif
      jarrg=int(xjthisdom-icd*.5)
      jarrgp=jarrg+1
      iarrg=int(yithisdom-icd*.5)
      iarrgp=iarrg+1
      ratx=xjthisdom-icd*.5-jarrg
      raty=yithisdom-icd*.5-iarrg
      do kch=1,mkzh-1
         k=kch
         kp=k+1
         if (zk.ge.sigh(k).and.zk.le.sigh(kp)) then
            ratz=(zk-sigh(k))/(sigh(kp)-sigh(k))
            goto 100
         endif
      enddo
      if (zk.ge.sigh(mkzh).and.zk.le.1.050) then
         k=mkzh-1
         kp=mkzh
         ratz=(zk-sigh(k))/(sigh(kp)-sigh(k)) ! linearly extrapolate
         goto 100
      endif
      write(iup,*)'In finterp, point is outside of valid sigma range.'
      write(iup,*)'zk=',zk
      stop
 100  continue
      finterp=0.
      totwgt=0.
      if (arr(iarrgp,jarrg,kp).ne.rmsg) then
         wgt=(1.-ratx)*(   raty)*(   ratz)
         finterp=finterp+wgt*arr(iarrgp,jarrg,kp)
         totwgt=totwgt+wgt
      endif
      if (arr(iarrgp,jarrgp,kp).ne.rmsg) then
         wgt=(   ratx)*(   raty)*(   ratz)
         finterp=finterp+wgt*arr(iarrgp,jarrgp,kp)
         totwgt=totwgt+wgt
      endif
      if (arr(iarrg,jarrg,kp).ne.rmsg) then
         wgt=(1.-ratx)*(1.-raty)*(   ratz)
         finterp=finterp+wgt*arr(iarrg,jarrg,kp)
         totwgt=totwgt+wgt
      endif
      if (arr(iarrg,jarrgp,kp).ne.rmsg) then
         wgt=(   ratx)*(1.-raty)*(   ratz)
         finterp=finterp+wgt*arr(iarrg,jarrgp,kp)
         totwgt=totwgt+wgt
      endif
      if (arr(iarrgp,jarrg,k).ne.rmsg) then
         wgt=(1.-ratx)*(   raty)*(1.-ratz)
         finterp=finterp+wgt*arr(iarrgp,jarrg,k)
         totwgt=totwgt+wgt
      endif
      if (arr(iarrgp,jarrgp,k).ne.rmsg) then
         wgt=(   ratx)*(   raty)*(1.-ratz)
         finterp=finterp+wgt*arr(iarrgp,jarrgp,k)
         totwgt=totwgt+wgt
      endif
      if (arr(iarrg,jarrg,k).ne.rmsg) then
         wgt=(1.-ratx)*(1.-raty)*(1.-ratz)
         finterp=finterp+wgt*arr(iarrg,jarrg,k)
         totwgt=totwgt+wgt
      endif
      if (arr(iarrg,jarrgp,k).ne.rmsg) then
         wgt=(   ratx)*(1.-raty)*(1.-ratz)
         finterp=finterp+wgt*arr(iarrg,jarrgp,k)
         totwgt=totwgt+wgt
      endif
      if (totwgt.ne.0.) then
         finterp=finterp/totwgt
      else
         finterp=rmsg
      endif
      return
      end
