c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine getrnum(string,ipos,mkzh,rval)
      character string*(*)
c
      include 'comconst'
c
      lenstr=len(string)
      do 10 i=ipos,lenstr
         if (string(i:i).eq.';'.or.string(i:i).eq.','.or.
     &         string(i:i).eq.' '.or.string(i:i+1).eq.'fb') then
            goto 20
         endif
   10 continue
   20 ilast=i-1
      read(string(ipos:ilast),fmt=*,err=50) rval
      ibrel=0
      if (string(i:i+1).eq.'fb') then
         ibrel=1
         ilast=ilast+2
      endif
c
c   If reading a level specifier intended to be taken as a vertical level
c   index in bottom-to-top order (specified by appended "fb"), then
c   convert to vertical level index in standard top-to-bottom order.
c
      if (ibrel.eq.1) then
         rval=float(mkzh)+1.-rval
      elseif (ibrel.eq.-1) then
         rval=-(float(mkzh)+1.-rval)
      endif
      ipos=ilast+2
      return
   50 write(iup,*)'   Can''t read list-directed format from "',
     &   string(ipos:ilast),'"'
      stop
      end
