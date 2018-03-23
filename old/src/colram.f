c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine colram (xcs,ycs,ncs,iai,iag,nai)
      dimension xcs(*),ycs(*),iai(*),iag(*)
      common /emap/ llcolor,lllinw,llndot,mpfillco(6),llmask
c
c   This is a user-supplied subroutine that is used by NCAR Graphics
c   for making color-filled maps.
c
        itm=1
        do 101 i=1,nai
          if (iai(i).lt.0) itm=0
  101   continue
        if (itm.ne.0) then
          itm=0
          do 102 i=1,nai
            if (iag(i).eq.1) itm=iai(i)
  102     continue
          if (itm.gt.0) then
c
c  set area fill color index.
c
            call gsfaci(mpfillco(mapaci(itm)))
c
            call gfa (ncs-1,xcs,ycs)
          endif
        endif
        return
      end
