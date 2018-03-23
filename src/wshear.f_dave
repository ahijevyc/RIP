c--------------------------------------------------------------------
      subroutine wshear (iw,
     +                      u, v, ght, prsf, ter, wk, miy, mjx, mkzh,
     +                      vlev_bot, vlev_top,
     +                      cvcor, cvcorunits)



c  Calculates the wind at vlev_top minus vlev_bot.
c  Uses the closest model vertical level.
c  The requested levels may not match a model level exactly if 
c  the requested levels vlev_bot and vlev_top are in terms of  
c  height or pressure, but the program will use the closest model levels
c  at or beneath each one. 
c  If the difference between the requested vertical level
c  and the closest model level is greater than dmax, then a warning
c  is produced.  No interpolation to the requested level is performed.
c 
c
c	
c	Everything here was originally copied from brnshr.f of the RIP_20011208
c	package.
c	Modified by Dave Ahijevych (MMM/NCAR).
c
c
c	Description and date of code changes:
c
c  20100426
c  Return wind shear magnitude if iw=3.  
c
c  20100424
c
c  Replaced references to sigma level with model level vertical (k) index. 
c  This is the way the RIP User's guide refers to it.
c  Added debug options
c
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
c		's'	= model vertical level (k) index 
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
c	Added model vertical level (k) coordinates option.
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
c			first model vertical level beneath it exceeds
c			dmax, issue a warning.
c                       Only used for debugging purposes.
c
c
c	pfsf		: REAL
c			3-dimensional array of 
c			pressure on model surfaces
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
      INTEGER    debug


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

      debug = 0
      IF (debug.GE.1) THEN
         WRITE(iup,*)"In subroutine wshear..."
         WRITE(iup,*)"iw =", iw, " miy=", miy, " mjx=", mjx
         WRITE(iup,*)"mkzh =", mkzh, " vlev_bot=", vlev_bot
         WRITE(iup,*)"vlev_top=", vlev_top, " cvcor=", cvcor
         WRITE(iup,*)"cvcorunits=", cvcorunits
      ENDIF 

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
         IF (debug.GE.1) WRITE(iup,*)"Using height above ground "//
     +               "level (m) as vertical coordinate"
      ELSE IF (cvcor.EQ.'z') THEN
         dmax   = 300.
         IF (debug.GE.1) WRITE(iup,*)"Using geopotential height (m)"//
     +               " as vertical coordinate"
      ELSE IF (cvcor.EQ.'p') THEN
         IF (debug.GE.1) WRITE(iup,*)"Using pressure (hPa)"//
     +               " as vertical coordinate"
         dmax   = 50.
         IF (vlev_top.GT.0.) vlev_top = -1. * vlev_top
         IF (vlev_bot.GT.0.) vlev_bot = -1. * vlev_bot
      ELSE IF (cvcor.EQ.'s') THEN
         WRITE(iup,*)"Using k-index"//
     +               " as vertical coordinate"
         IF (vlev_bot.GT.mkzh) THEN
            WRITE(iup,*) 'In subroutine wshear'
            WRITE(iup,*) 'Your request exceeds the number '//
     +       'of vertical levels (max ', mkzh, ')'
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
c	Same with k-index, because I 
c	believe the k-index index increases downward.
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
               IF(debug.ge.2) WRITE(iup,*) "vlev=",vlev
            ELSE IF (cvcor.EQ.'p') THEN

c		Use pressure on model surfaces.
c
c		Since pressure increases with k-index
c		(the opposite of what happens
c		with height), negate vlev, vlev_top, 
c		and vlev_bot so the inequality
c		tests that find the closest model 
c		layer beneath vlev_top and vlev_bot
c		work correctly for both pressure
c		coordinates and height coordinates.

               IF (debug.ge.2) WRITE(iup,*) 'negating prsf'
               vlev = -1. * prsf(i,j,k)
               IF (debug.ge.2) WRITE(iup,*) 'prsf=',prsf(i,j,k)
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

c		Warn the user if the closest model level below
c		vlev_top is more than dmax cvcorunits away.

               IF ((vlev - vlev_top).GT. dmax) THEN
                  WRITE(iup,*) "In subroutine wshear"
                  WRITE(iup,'(1X,"at (i,j) (",I4,","I4,")")') i,j
                  WRITE(iup,'(1X,A,F7.0,A," is ",F6.1,A," away.")')
     +             " the closest model level below ",
     +             vlev_top, cvcorunits, vlev - vlev_top, cvcorunits
               ENDIF


            endif
c		This originally used 500. instead of vlev_bot.
c		This originally had vlev.lt.vlev_bot
            if ((vlev.le.vlev_bot).and.(k_bot.eq.0)) then
               k_bot = k

c		Warn the user if the closest model level below
c		vlev_bot is more than dmax cvcorunits away.

               IF ((vlev - vlev_bot).GT. dmax) THEN
                  WRITE(iup,*) "In subroutine wshear"
                  WRITE(iup,'(1X,"at (i,j) (",I4,","I4,")")') i,j
                  WRITE(iup,'(1X,A,F7.0,A," is ",F6.1,A," away.")')
     +             " the closest model level below ", 
     +             vlev_bot, cvcorunits, vlev - vlev_bot, cvcorunits
               ENDIF


            endif
         enddo
         if (k_top.eq.0) then
            write(iup,*)'In wshear, couldn''t find the first'
            write(iup,*)' model level below ',vlev_top, cvcorunits
            stop
         endif
         if (k_bot.eq.0) then
            write(iup,*)'In wshear, couldn''t find the first'
            write(iup,*)' model level below ',vlev_bot, cvcorunits
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
         IF (debug.GE.1 .AND. i.EQ.6 .AND. j.EQ.51) THEN
            WRITE(iup,*) "u,v top = ",u(i,j,k_top),v(i,j,k_top)
            WRITE(iup,*) "u,v bot = ",u(i,j,k_bot),v(i,j,k_bot)
            WRITE(iup,*) "u,v shr = ",u(i,j,k_top) - u(i,j,k_bot),
     +                                v(i,j,k_top) - v(i,j,k_bot)
         ENDIF




c
c	Depending on the value of iw, assign
c	either the u-component or the v-component
c	of the wind shear to wk.
c	
         
         IF (iw .EQ. 1) THEN
             wk(i,j) = u(i,j,k_top) - u(i,j,k_bot)
         ELSEIF (iw .EQ. 2) THEN
             wk(i,j) = v(i,j,k_top) - v(i,j,k_bot)
         ELSEIF (iw .EQ. 3) THEN
             wk(i,j) = SQRT((u(i,j,k_top) - u(i,j,k_bot))**2. +
     &                      (v(i,j,k_top) - v(i,j,k_bot))**2.)
         ELSE
             WRITE(iup,*) 'Bad iw value: ', iw
             STOP
         ENDIF

c
      enddo
      enddo
c
      return
      end
