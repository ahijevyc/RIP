c--------------------------------------------------------------------
      subroutine wshear (iw,
     +                      u, v, ght, prsf, ter, wk, miy, mjx, mkzh,
     +                      vlev_bot, vlev_top,
     +                      cvcor, cvcorunits)
c	
c	Everything here was originally copied from brnshr.f of the RIP_20011208
c	package.
c	Modified by Dave Ahijevych (MMM/NCAR).
c
c
c	Description and date of code changes:
c
c
c	20011213
c
c
c
c	Changed cvcor from a 4-character string to a
c	single character.
c
c	cvcor		: CHARACTER
c
c		Possible values:
c
c		'a'	= meters AGL
c		'z'	= geopotential height (m)
c		'p'	= pressure (mb or hPa)
c		's'	= index of sigma surface
c
c
c
c	Added cvcorunits to argument list.
c	Now this variable is determined in the 
c	calling procedure.
c
c	This variable was originally called cunits,
c	but changed to cvcorunits to match its 
c	name in fields.f.
c
c	Added sigma coordinates option.
c
c
c	20011211
c
c	Added argument cvcor to end of dummy argument list.
c	This is a 4-character string that describes the vertical
c	coordinate system of the vlev_bot and vlev_top 
c	arguments.
c
c	cvcor		: CHARACTER*4
c
c		Possible values:
c
c		'zagl'	= meters AGL
c		'z   '	= geopotential height (m)
c		'p   '	= pressure (mb or hPa)
c
c
c	Also modified the CALL wshear( )
c	line in fields.f to accommodate
c	the additional argument.
c
c
c	Added REAL vlev.  It holds the value of the 
c	current vertical level in the loop.
c
c	vlev		: REAL
c			In the vertical coordinate loop,
c			holds the current 
c			height, pressure, or whatever.
c
c
c
c	Added CHARACTER string cunits
c
c	cunits		: CHARACTER*16
c			Units of vertical coordinate system
c			Currently only used for nice-looking
c			debugging messages.
c
c
c	Added REAL dmax
c
c	dmax		: REAL
c			If the vertical separation between 
c			vcor_top or vcor_bot and the
c			first sigma layer beneath it exceeds
c			dmax, issue a warning.
c                       Only used for debugging purposes.
c
c
c	pfsf		: REAL
c			3-dimensional array of 
c			pressure on sigma surfaces
c			(mb or hPa)

c
c
c	20011208
c
c	Changed subroutine name from brnshr to wshear
c
c	Added two arguments to the end of the argument list.
c	These are the bottom and top level to use in the 
c	vertical wind shear calculation.
c	Before this, the bottom and top level were hard coded
c	to 500. and 6000. m AGL.
c
c	vlev_bot	: REAL
c			bottom level in wind shear calculation
c
c	vlev_top	: REAL
c			top level in wind shear calculation
c
c	Added variable declaration for REALs vlev_bot and vlev_top.
c
c	Added test to ensure the bottom level is below the 
c	top level.
c
c	Renamed k500 k_bot.
c	Renamed k6000 k_top.
c	
c
c
c
c	Added argument iw to the front of the dummy argument list
c	and its explicit integer type declaration.
c	This variable works the same as iw in bshear.f.
c	It determines whether the u-component or the
c	v-component of the wind will be assigned to 
c	2-dimensional work field wk.
c	wk is the output field.
c
c	iw		: INTEGER
c		1		= return u-component
c		other value	= return v-component
c
c


      INTEGER    iw
      REAL       vlev_bot, vlev_top
      REAL       vlev
      REAL       dmax
      CHARACTER  cvcor, cvcorunits*16



c
c   This program originally calculated the
c   Bulk Richardson Number Shear, as found in Stensrud et al. (1997)
c
      dimension u(miy,mjx,mkzh), v(miy,mjx,mkzh), ght(miy,mjx,mkzh),
     &   prsf(miy, mjx,mkzh), ter(miy,mjx), wk(miy,mjx)
c

c	This was something in the original code related to
c	the storm motion vector, I think. -- Dave
c      dimension ucross(500),vcross(500)
c
      include 'comconst'
c



c
c	Make sure the vertical coordinate is one we recognize.
c
c	Also assign cvcorunits and dmax, depending on the 
c	vertical coordinate system.
c
c	If the vertical coordinate system is pressure,
c	negate vlev_bot and vlev_top.
c	This makes the inequality tests work for both
c	pressure and height coordinate systems.
c

      dmax   = 0.
      IF (cvcor.EQ.'a') THEN
         dmax   = 300.
c         WRITE(iup,*)"In subroutine wshear..."
c         WRITE(iup,*)"Using height above ground level (m)"//
c     +               " as vertical coordinate"
      ELSE IF (cvcor.EQ.'z') THEN
         dmax   = 300.
c         WRITE(iup,*)"In subroutine wshear..."
c         WRITE(iup,*)"Using geopotential height (m)"//
c     +               " as vertical coordinate"
      ELSE IF (cvcor.EQ.'p') THEN
         dmax   = 50.
         IF (vlev_top.GT.0.) vlev_top = -1. * vlev_top
         IF (vlev_bot.GT.0.) vlev_bot = -1. * vlev_bot
      ELSE IF (cvcor.EQ.'s') THEN
         IF (vlev_bot.GT.mkzh) THEN
            WRITE(iup,*) 'In subroutine wshear'
            WRITE(iup,*) 'Your request exceeds the number '//
     +       'of available sigma layers (max ', mkzh, ')'
            WRITE(iup,*)'Stopping'
            STOP
         ENDIF
         IF (vlev_top.GT.0.) vlev_top = -1. * vlev_top
         IF (vlev_bot.GT.0.) vlev_bot = -1. * vlev_bot
      ELSE
         WRITE(iup,*)"In subroutine wshear..."
         WRITE(iup,*)"Vertical coordinate not recognized: "//cvcor
         WRITE(iup,*)"Stopping."
         STOP
      ENDIF 


c
c	Make sure the 1st vertical level in the 
c	argument list is lower than the 2nd argument
c	level.
c	This subroutine assumes the lower vertical 
c	level is listed first.
c
c	Note: If the vertical coordinate system is
c	pressure, i.e. cvcor = 'p',
c	vlev_bot and vlev_top must first be negated.
c	Same with sigma coordinates, because I 
c	believe the sigma level index increases downward.
c
c

      IF (vlev_bot .GE. vlev_top) THEN
         WRITE(iup,*) "In subroutine wshear..."
         WRITE(iup,*) "Requested bottom level ", vlev_bot
         WRITE(iup,*) "the same as or higher than the top level ",
     +              vlev_bot, "!"
         WRITE(iup,*) "Stopping"
         STOP
      ENDIF


      do j = 1, mjx-1
      do i = 1, miy-1
c         sdh = 0.
c         su = 0.
c         sv = 0.
c
c	This originally was hard-coded to:
c      Find the first level below 500 m AGL and the first level
c      below 6000 m AGL
c
c	But now it uses the arguments vlev_bot and vlev_top instead.
c
c
         k_bot = 0
         k_top = 0
         do k = 1,mkzh


c
c	As the vertical coordinate system dictates,
c	assign the current vertical coordinate value
c	to vlev.
c
 
            IF (cvcor.EQ.'a') THEN
c		Use height above ground level
               vlev = ght(i,j,k)-ter(i,j)
            ELSE IF (cvcor.EQ.'z') THEN
c		Use ordinary geopotential height
               vlev = ght(i,j,k)
            ELSE IF (cvcor.EQ.'p') THEN

c		Use pressure on sigma surfaces.
c
c		Since pressure increases with sigma
c		levels (the opposite of what happens
c		with height), negate vlev, vlev_top, 
c		and vlev_bot so the inequality
c		tests that find the closest sigma 
c		layer beneath vlev_top and vlev_bot
c		work correctly for both pressure
c		coordinates and height coordinates.

               vlev = -1. * prsf(i,j,k)
c               WRITE(iup,*) 'prsf=',prsf(i,j,k)
               IF (vlev_top.GT.0.) vlev_top = -1. * vlev_top
               IF (vlev_bot.GT.0.) vlev_bot = -1. * vlev_bot
            ELSE IF (cvcor.EQ.'s') THEN
               vlev = -1 * k
               IF (vlev_top.GT.0.) vlev_top = -1. * vlev_top
               IF (vlev_bot.GT.0.) vlev_bot = -1. * vlev_bot
            ELSE
               WRITE(iup,*)"In subroutine wshear..."
               WRITE(iup,'(1X,"at (i,j) (",I4,","I4,")")') i,j
               WRITE(iup,*) "Couldn't get current vertical coordinate"
               WRITE(iup,*) "Strange cvcor: "//cvcor
               WRITE(iup,*) "Stopping."
               STOP
            ENDIF 



c		This originally used 6000. instead of vlev_top.
c		This originally had vlev.le.vlev_top
            if ((vlev.le.vlev_top).and.(k_top.eq.0)) then
               k_top = k

c		Warn the user if the closest sigma layer below
c		vlev_top is more than dmax cvcorunits away.

               IF ((vlev - vlev_top).GT. dmax) THEN
                  WRITE(iup,*) "In subroutine wshear"
                  WRITE(iup,'(1X,"at (i,j) (",I4,","I4,")")') i,j
                  WRITE(iup,'(1X,A,F7.0,A," is ",F6.1,A," away.")')
     +             " the closest sigma level below ",
     +             vlev_top, cvcorunits, vlev - vlev_top, cvcorunits
               ENDIF


            endif
c		This originally used 500. instead of vlev_bot.
c		This originally had vlev.lt.vlev_bot
            if ((vlev.le.vlev_bot).and.(k_bot.eq.0)) then
               k_bot = k

c		Warn the user if the closest sigma layer below
c		vlev_bot is more than dmax cvcorunits away.

               IF ((vlev - vlev_bot).GT. dmax) THEN
                  WRITE(iup,*) "In subroutine wshear"
                  WRITE(iup,'(1X,"at (i,j) (",I4,","I4,")")') i,j
                  WRITE(iup,'(1X,A,F7.0,A," is ",F6.1,A," away.")')
     +             " the closest sigma level below ", 
     +             vlev_bot, cvcorunits, vlev - vlev_bot, cvcorunits
               ENDIF


            endif
         enddo
         if (k_top.eq.0) then
            write(iup,*)'In wshear, couldn''t get the first'
            write(iup,*)' sigma layer below ',vlev_top, cvcorunits
            stop
         endif
         if (k_bot.eq.0) then
            write(iup,*)'In wshear, couldn''t get the first'
            write(iup,*)' sigma layer below ',vlev_bot, cvcorunits
            stop
         endif


 
c
c	This was something in the original code related to
c	the storm motion vector, I think. -- Dave
c      dimension ucross(500),vcross(500)
c
c         ucross(k)=.25*(u(i,j,k)+u(i+1,j,k)+
c     &               u(i,j+1,k)+u(i+1,j+1,k))
c         vcross(k)=.25*(v(i,j,k)+v(i+1,j,k)+
c     &               v(i,j+1,k)+v(i+1,j+1,k))
c
c
c      Calculate a 0-6 km AGL mean wind
c
c         do k = mkzh, k_top, -1
c            dh = prsf(i,j,k) - prsf(i,j,k-1)
c            sdh = sdh + dh
c            su = su + dh*ucross(k)
c            sv = sv + dh*vcross(k)
c         enddo
c         ua = su / sdh
c         va = sv / sdh
c
c      Calculate a 0-500 m AGL mean wind
c
c         sdh = 0.
c         su = 0.
c         sv = 0.
c         do k = mkzh, k_bot, -1
c            dh = prsf(i,j,k) - prsf(i,j,k-1)
c            sdh = sdh + dh
c            su = su + dh*ucross(k)
c            sv = sv + dh*vcross(k)
c         enddo
c         u0 = su / sdh
c         v0 = sv / sdh

c	This originally was in the code...
c   Bulk Richardson Number Shear:
c
c         wk(i,j) = 0.5 * ( (ua - u0)**2 + (va - v0)**2 )



c	Added for testing purposes...
c         IF (i.EQ.6 .AND. j.EQ.51) THEN
c            WRITE(iup,*) "u,v top = ",u(i,j,k_top),v(i,j,k_top)
c            WRITE(iup,*) "u,v bot = ",u(i,j,k_bot),v(i,j,k_bot)
c            WRITE(iup,*) "u,v shr = ",u(i,j,k_top) - u(i,j,k_bot),
c     +                                v(i,j,k_top) - v(i,j,k_bot)
c         ENDIF




c
c	Depending on the value of iw, assign
c	either the u-component or the v-component
c	of the wind shear to wk.
c	
         
         IF (iw .EQ. 1) THEN
             wk(i,j) = u(i,j,k_top) - u(i,j,k_bot)
         ELSE
             wk(i,j) = v(i,j,k_top) - v(i,j,k_bot)
         ENDIF

c
      enddo
      enddo
c
      return
      end
