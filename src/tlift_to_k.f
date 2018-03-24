      SUBROUTINE tlift_to_k (t, dt, z0, p, q, dq, tl, tm, i, j, kk, 
     &mkzh,
     &kmaxb, k_to_lift_to)
      DIMENSION t(mkzh), q(mkzh), p(mkzh), z(mkzh), tb(mkzh)

c Copied from tlift.f_dave 20130306.

c Instead of assigning the minimum buoyancy to tl, 
c assign the buoyancy at the prescribed level k_to_lift_to.

c Unlike tlift, variables kminb, tm, and kmaxb are not defined in tlift_to_k.


! Intent is optional, but it helps prevent bugs.
! intent(in) means the variable can enter, but can't be changed.
! intent(out) means the variable is set inside the procedure and sent 
!             back to the main program with any initial values ignored.
! intent(inout) is the default. It means the variable comes in with
!               a value and leaves with a value.
      real,    intent(in) :: t, dt, z0, p, q, dq
      real,    intent(out):: tl, tm ! slightly different than tlift (see above)
      integer, intent(in) :: i, j, kk, mkzh, k_to_lift_to ! input
      real c1

! k_to_lift_to has to be a valid level. 
! In pvocalc it was initialized to 0 if kbmin3 was 
! not already created in the form of file fort.19.
! So make sure k_to_lift_to is not 0.
      if (k_to_lift_to.eq.0) then 
        write(*,*) 'tlift_to_k: k_to_lift_to is undefined.'
        write(*,*)'Did you write kbmin3 to fort.19?'
        stop
      endif
c      open (unit=50,file="tliftnew_debug.out",status="unknown")
      rd=287.
c      pt=400.
      pt=300.
!highest level to look for neg buoyancy
      cp=1005.
      lv=2.5e6
      grav=9.81
      p00=1.e3
      dp=10.
      bsat=243.5
      asat=17.67
      t00=273.15
! 
!---Average parcel through specificed pressure depth pavd
! 
      pavd=50.
      pavd2=25.
      lvl=kk
   10 lvl=lvl+1
      dp=p(lvl)-p(kk)
      IF (lvl.gt.mkzh.or.dp.gt.pavd2) THEN
        lvlbot=lvl-1
        pbot=p(lvlbot)
        tbot=t(lvlbot)
        qbot=q(lvlbot)
      ELSE
        GO TO 10
      END IF
      ptop=pbot-pavd
      lvl=lvlbot
   20 lvl=lvl-1
      IF (p(lvl).lt.ptop) THEN
        delta=(ptop-p(lvl+1))/(p(lvl)-p(lvl+1))
        ttop=t(lvl+1)+(t(lvl)-t(lvl+1))*delta
        qtop=q(lvl+1)+(q(lvl)-q(lvl+1))*delta
        pave=(pbot+ptop)/2.
        tave=(tbot+ttop)/2.
        qave=(qbot+qtop)/2.
      ELSE
        GO TO 20
      END IF
! 
!----end parcel average
! 
      t0=tave+dt
      q0=max(qave+dq,1.e-4)
      p0=pave
! 
!      t0=t(kk) + dt
!      q0=max(q(kk)+dq,1.e-4)
!      p0=p(kk)
      z(kk)=z0
      es00=6.11
      tminb=10.
      tmaxb=-100.
! Compute geopotential height
      DO k=kk-1,1,-1
        tv1=t(k)*(1.+0.61*q(k))
        tv2=t(k+1)*(1.+0.61*q(k+1))
        z(k)=z(k+1)-rd*0.5*(tv1+tv2)*(p(k)-p(k+1))/(grav*0.5*(p(k)+p(k+
     &   1)))
      END DO
! Find z_lcl relative to parcel with max MSE
! Check to see if parcel already saturated
      tlc=t0-t00
      esl=es00*exp(asat*tlc/(bsat+tlc))
      qsl=0.622*esl/(p0-esl)
      ! if saturated 
      ! new temperature t0 and mixing ratio q0 blends effects of dt and dq.
      IF (q0.gt.qsl) THEN
        ! get saturation temperature for q0
        esl = p0*q0/(0.622+q0)! get water vapor pressure esl from pressure and mixing ratio 
        c1 = alog(esl/es00) ! just an abbreviation to use below
        tlc = bsat*c1/(asat - c1) 
        ! remember to add 273.15 (t00) to tlc.
        ! blend saturation temperature for q0 and new t0
        t0 = ((tlc+t00) + t0)/2. 

        ! get blended saturation mixing ratio using temperature above.
        esl=es00*exp(asat*(t0-t00)/(bsat+t0-t00))
        q0=0.622*esl/(p0-esl)
      END IF
! 
      the0=t0*(p0/p00)**(-rd/cp)
      zmse0=cp*t0+grav*z0+lv*q0
      tlk=t00+(bsat/asat)*alog(q0*p0/(0.622*es00))
      iter=0
   40 tlc=tlk-t00
!     pl treated like lifted condensation level in IF-block below
      pl=p00*((the0/tlk)**(-cp/rd))
      esl=es00*exp(asat*tlc/(bsat+tlc))
      qsl=0.622*esl/(pl-esl)
      tlk2=tlk+((bsat+tlc)/asat)*(q0-qsl)/q0
      IF (abs(tlk2-tlk).gt.0.01) THEN
!iterate for LCL
        iter=iter+1
        tlk=tlk2
        IF (iter.lt.10) GO TO 40
      END IF
      zl=z0-cp*(tlk2-t0)/grav
!height of lcl
      DO k1=kk,1,-1
        IF (p(k1).le.pt) GO TO 60
      END DO
! 
   60 tp=the0*(p(kk)/p00)**(rd/cp)
      t2a=tp
      if (j.eq.395.and.i.eq.400) write(6,*) 'k1 = ',k1
      DO k2=kk-1,k1,-1
        IF (p(k2).gt.pl) THEN
          tp=the0*(p(k2)/p00)**(rd/cp)
          tb(k2)=tp*(1.+0.61*q0)-t(k2)*(1.+0.61*q(k2))
          t2a=tp
          zmsep=cp*tp+grav*z(k2)+lv*q0
c          IF (tb(k2).lt.tminb) THEN
c            tminb=tb(k2)
c            kminb=k2
c          END IF
c          IF (tb(k2).gt.tmaxb) THEN
c            tmaxb=tb(k2)
c            kmaxb=k2
c          END IF

c Instead of defining tminb by the minimum buoyancy, 
c return the buoyancy at the prescribed level k_to_lift_to.
c kminb, tmaxb, and kmaxb are not used.
c If you change this again remember there are 2 places.
          IF (k2.eq.k_to_lift_to) tmaxb=tb(k2)
c
          if (j.eq.395.and.i.eq.400) then
          write (6,990) kk,iter,k2,tmaxb
          end if
c
        ELSE
!saturated
! if here, parcel is saturated at level of max negative buoyancy
          iter=0
   70     t2c=t2a-273.15
          es2a=6.11*exp(asat*t2c/(bsat+t2c))
          q2a=0.622*es2a/(p(k2)-es2a)
          zmsep=cp*t2a+grav*z(k2)+lv*q2a
          t2=t2a+0.25*(zmse0-zmsep)/cp
          q2=q2a+0.25*(zmse0-zmsep)/lv
          iter=iter+1
          IF (iter.gt.10) THEN
!        print*,'too many iterations.',i,j,t2,t2a
            GO TO 90
          END IF
          IF (abs(t2a-t2).gt.0.01) THEN
            q2a=q2
            t2a=t2
            GO TO 70
            tp=t2
          END IF
   90     CONTINUE
          tb(k2)=t2*(1.+0.61*q2)-t(k2)*(1.+0.61*q(k2))


c          IF (tb(k2).lt.tminb) THEN
c            tminb=tb(k2)
c            kminb=k2
c          END IF
c          IF (tb(k2).gt.tmaxb) THEN
c            tmaxb=tb(k2)
c            kmaxb=k2
c          END IF

c Instead of defining tminb by the minimum buoyancy, 
c return the buoyancy at the prescribed level k_to_lift_to.
c kminb, tmaxb, and kmaxb are not used.
c If you change this again remember there are 2 places.
          IF (k2.eq.k_to_lift_to) tmaxb=tb(k2)
c
          if (j.eq.395.and.i.eq.400) then
          write (6,990) kk,iter,k2,tmaxb
          end if
c
        END IF
      END DO
      tl=tminb
      tm=tmaxb
 990  format (3(i3,1x),f10.4)
c      close (unit=50)
      RETURN
      END
  