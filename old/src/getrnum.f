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
     &         string(i:i).eq.' ') then
            goto 20
         endif
   10 continue
   20 ilast=i-1
      ibrel=0
      if (string(ipos:ipos).eq.'b') then
         ibrel=1
         ipos=ipos+1
      elseif (string(ipos:ipos).eq.'-'.and.
     &        string(ipos+1:ipos+1).eq.'b') then
         ibrel=-1
         ipos=ipos+2
      endif
      read(string(ipos:ilast),fmt=*,err=50) rval
c
c   If reading a level specifier intended to be taken as a sigma level
c   index in bottom-to-top order (specified by the preceding "b"), then
c   convert to sigma level index in standard top-to-bottom order.
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
