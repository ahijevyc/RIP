c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine vc2dcalc(vc3d,cvcor,sigh,rcrag,rcrbg,nscrs,
     &   prs_tsf,pstx,idoground_tsf,vcground,mabpl,miy,mjx,mkzh)
c
      dimension vc3d(miy,mjx,mkzh),rcrag(2),rcrbg(2),vcground(mabpl),
     &   sigh(mkzh),prs_tsf(miy,mjx),pstx(miy,mjx)
      character cvcor*1
c
      include 'comconst'
      dimension vc2d(500,250)
      common /vctran/ nscross,mkzhcross,ivcs,vc2d
c
c   Interpolate vc3d to x-section.
c
      caxgn=1.+(rcrag(2)-xjcorn)*refrat
      caygn=1.+(rcrag(1)-yicorn)*refrat
      cbxgn=1.+(rcrbg(2)-xjcorn)*refrat
      cbygn=1.+(rcrbg(1)-yicorn)*refrat
      xj1t=caxgn
      xj2t=cbxgn
      yi1t=caygn
      yi2t=cbygn
      if (xj1t.le.1.5.or.xj1t.ge.mjx-.5.or.
     &    xj2t.le.1.5.or.xj2t.ge.mjx-.5.or.
     &    yi1t.le.1.5.or.yi1t.ge.miy-.5.or.
     &    yi2t.le.1.5.or.yi2t.ge.miy-.5) then
         write(iup,*)'In vc2dclac: Cross sec. endpoints must be'
         write(iup,*)'between 1.5 and (miy-.5) or (mjx-.5).'
         stop
      endif
c
      do ls=1,nscrs
         if (cvcor.ne.'s'.or.idoground_tsf.eq.1) then
            posx=xj1t+(ls-1.)/(nscrs-1.)*(xj2t-xj1t)-.5
            posy=yi1t+(ls-1.)/(nscrs-1.)*(yi2t-yi1t)-.5
            jl=int(posx)
            jr=jl+1
            ib=int(posy)
            it=ib+1
            ratlr=posx-jl
            ratbt=posy-ib
         endif
         do k=1,mkzh
            if (cvcor.ne.'s') then
               if (vc3d(it,jl,k).eq.rmsg.or.
     &             vc3d(it,jr,k).eq.rmsg.or.
     &             vc3d(ib,jl,k).eq.rmsg.or.
     &             vc3d(ib,jr,k).eq.rmsg) then
                  write(iup,*)'vc3d has an rmsg value.'
                  stop
               else
                  wk1=vc3d(it,jl,k)
                  wk2=vc3d(it,jr,k)
                  wk3=vc3d(ib,jl,k)
                  wk4=vc3d(ib,jr,k)
               endif
               vc2d(ls,k)=
     &            (1.-ratlr)*(   ratbt)*wk1+
     +            (   ratlr)*(   ratbt)*wk2+
     +            (1.-ratlr)*(1.-ratbt)*wk3+
     +            (   ratlr)*(1.-ratbt)*wk4 
            else
               vc2d(ls,k)=sigh(k)
            endif
         enddo
         if (idoground_tsf.eq.1) then
            if (prs_tsf(it,jl).eq.rmsg.or.
     &          prs_tsf(it,jr).eq.rmsg.or.
     &          prs_tsf(ib,jl).eq.rmsg.or.
     &          prs_tsf(ib,jr).eq.rmsg) then
               write(iup,*)'prs_tsf has an rmsg value.'
               stop
            else
               pt1=prs_tsf(it,jl)
               pt2=prs_tsf(it,jr)
               pt3=prs_tsf(ib,jl)
               pt4=prs_tsf(ib,jr)
            endif
c
c         Interpolate true surface pressure to x-sec, put in vcground
c
            vcground(ls)=
     &         (1.-ratlr)*(   ratbt)*pt1+
     +         (   ratlr)*(   ratbt)*pt2+
     +         (1.-ratlr)*(1.-ratbt)*pt3+
     +         (   ratlr)*(1.-ratbt)*pt4 
c
c         Convert to sigma.  Since idoground_tsf=1 only for pressure-level
c         data sets, sigma should be hydrostatic form.  Also, for p-level
c         data, pstx is constant, so just use lower left corner value.
c
            vcground(ls)=(vcground(ls)-ptop)/pstx(1,1)
         endif
      enddo
c
c   Calculate the vertical corrdinate at the ground surface,
c   by linear interpolation/extrapolation
c
      if (idoground_tsf.eq.0) then
         do ls=1,nscrs
            vcground(ls)=vc2d(ls,mkzh)+(vc2d(ls,mkzh)-
     &         vc2d(ls,mkzh-1))*(1.-sigh(mkzh))/
     &         (sigh(mkzh)-sigh(mkzh-1))
         enddo
      else
         do ls=1,nscrs
            if (vcground(ls).le.sigh(1)) then ! above model top - extrap.
               vcground(ls)=vc2d(ls,1)-
     &            (vc2d(ls,2)-vc2d(ls,1))*
     &            (sigh(1)-vcground(ls))/(sigh(2)-sigh(1))
            elseif (vcground(ls).ge.sigh(mkzh)) then ! below bottom - extrap.
               vcground(ls)=vc2d(ls,mkzh)+
     &            (vc2d(ls,mkzh)-vc2d(ls,mkzh-1))*
     &            (vcground(ls)-sigh(mkzh))/(sigh(mkzh)-sigh(mkzh-1))
            else  ! somewhere within model sigma levels
               do k=mkzh,2,-1
                  if (vcground(ls).le.sigh(k).and.
     &                vcground(ls).ge.sigh(k-1)) then
                     vcground(ls)=
     &                  ((vcground(ls)-sigh(k-1))*vc2d(ls,k)+
     &                   (sigh(k)-vcground(ls))*vc2d(ls,k-1))/
     &                  (sigh(k)-sigh(k-1))
                     goto 107
                  endif
               enddo
 107           continue
            endif
         enddo
      endif
c
c   Special things to be done for particular
c      vertical coordinates
c
      do ls=1,nscrs
         if (cvcor.eq.'z') then ! convert to height (km)
            vcground(ls)=-.001*sclht*alog(vcground(ls))

         elseif (cvcor.eq.'l') then ! convert to log pressure
            vcground(ls)=alog(vcground(ls))
         elseif (cvcor.eq.'x') then ! convert to exner function
            vcground(ls)=(vcground(ls))**gamma
         endif
         do k=1,mkzh
            if (cvcor.eq.'z') then ! convert to height (km)
               vc2d(ls,k)=-.001*sclht*alog(vc2d(ls,k))

            elseif (cvcor.eq.'l') then ! convert to log pressure
               vc2d(ls,k)=alog(vc2d(ls,k))
            elseif (cvcor.eq.'x') then ! convert to exner function
               vc2d(ls,k)=(vc2d(ls,k))**gamma
            endif
         enddo
      enddo
c
      return
      end
