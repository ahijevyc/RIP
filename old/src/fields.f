c
c
c	Everything here was originally copied from fields.f of the 
c	RIP_20011208 package.
c	Modified by Dave Ahijevych (MMM/NCAR).
c
c	Description and date of code changes:
c
c
c
c
c	20020219
c	When using keyword mply, a negative value implies division.
c	This opens the door to divide-by-zero error.
c	Therefore, ensure that you do not divide by zero.
c
c	I should program this to STOP with an error if this occurs, but
c	I cheat by setting undefined quotients equal to zero, just to get
c	through the program.  There is no "missing value" or "bad value" or
c	"undefined value."  Otherwise I would set the quotient equal to one
c	of them.
c
c	20020211
c
c	Added rfdmn and rfdmx (field min and max)
c	for thresholding field.
c	Thresholding is done after all adding and multiplying, but
c	before saving.
c
c
c	20020112
c
c	Added ismth to allow saving of smoothed fields.
c	This early smoothing is only done if cfeld starts with 'beta'.
c
c	Normally, the graphing subroutine does the smoothing
c	after the save process (in this subroutine).
c
c	A scratch field is filled with smoothed data and saved here.
c	The working field is unaltered.  This is because we don't
c	want to smooth once here, and then again in the graphing
c	routine.
c
c
c	20020105
c
c	Added rmply (analogous to raddf) for multiplying by 
c	another field raised to a particular power.
c	Just as a negative raddf value means subtraction,
c	a negative rmply value means division.
c
c
c	20011213
c
c	Added variations to the feld name for vertical 
c	wind shear.
c
c	ushxMMMNNN and vshxMMMNNN
c
c	The 4th character indicates the vertical 
c	coordinate system used in the 
c	6-digit level specification MMMNNN.
c	The character in the 'x' position   
c	follows the naming convention of
c	vcor, with an additional option 
c	'a' representing height 
c	above ground level (AGL). 
c
c	For the u-component of the vertical shear, 
c	the first 3 characters of feld are 'ush'.
c	Similarly, for the v-component of shear, 
c	the first 3 characters of feld are 'vsh'.
c
c
c	Examples:
c	
c	usha005030 = u-component of wind shear
c			from 0.5 to 3 km AGL 
c	ushz010060 = u-component of wind shear
c			from 1 to 6 km MSL 
c	vshp085070 = v-component of wind shear
c			from 850 to 700mb
c	vshs015005 = v-component of wind shear
c			from sigma level index 15 to 5
c		(like pressure, sigma increases downward)
c
c
c	As seen in the examples, 
c	heights (AGL and MSL) are specified
c	in hectometers, pressures in kPa,  
c	and sigma levels in index values.
c
c
c
c
c
c	Added to argument list for subroutine wshear.
c
c	cvcorunits	: CHARACTER*16
c			string (16 char max) holding
c			units of the
c			vertical coordinate system used 
c			in the wind shear
c			layer specification.
c
c	Added variable cvcorunits to the character
c	declaration section.
c
c	Changed the vertical coordinate argument to
c	subroutine wshear from a 4-character
c	to 1-character variable.
c
c
c	Added 's' (sigma surfaces) option.
c
c	Converted the requested sigma level indices
c	to actual sigma values via variable sigh.
c	The sigma values may be used in the plot label.
c
c	Right now, only sigma level indices may be used
c	in ushsMMMNNN and usfsMMMNNN.
c	In the future, I might add 'from the bottom'
c	notation 'b' and allow the user to specify the 
c	second from the bottom sigma level as 'b02'.
c
c 
c	20011211
c
c	I allowed different vertical
c	coordinates for the vertical wind shear calculation.
c	Therefore I added another argument to the end of the
c	argument list for subroutine wshear.
c	This is a character indicating the
c	vertical coordinate system used to request
c	the wind shear layer.
c
c	Current choices:
c		'a' = height (m AGL)
c		'z' = geopotential height (m MSL)
c		'p' = pressure (mb or hPa)
c		's' = sigma-surface (index)
c
c	Future choices:
c		't' = potential temp surface (K)
c
c
c	Changed z1 and z2 to vlev_bot vlev_top so they
c	are more generic names (z usually implies
c	geopotential height).
c
c
c
c	20011208
c	Added to large ELSE-IF block to accommodate
c	requests for wind shear stuff.
c	Allowed the bottom and top levels to be requested the
c	same way as for thickness calculations.
c
c	One difference, however, is that the levels may be 
c	requested in different units.
c
c
c      
c
c*********************************************************************c
c                                                                     c
      subroutine fields(cfeld,wk,indwk,icdwk,rlevl,rlavl,unwk,
     &   uuu,vvv,tmk,qvp,prs,ght,www,pstx,prs_tsf,dmap,xmap,ter,
     &   pstd,cor,xlus,sigh,unorth,vnorth,rstmv,rrfst,pslab1,
     &   fname,iendf1,sigf,incwk,ipl,iplstrt,idimn,rcrag,ismcp,ismth,
     &   rcrbg,cptyp,mdate,lverf,ydist,xdist,xseclen,
     &   nscrs,raddf,rmply,rfdmn,rfdmx,csave,lredo,ihrip,rhrip,rsepa,
     &   fullsigma,halfsigma,chrip,vardesc,lgrad,llapl,lhadv,
     &   igdir,iqgsm,plchun,casename,iendc,engplttl,ctjfl,rtjst,
     &   rtjen,iprog,cdiff,rdiff,ldfrl,xtime,tacch,
     &   xtimeavl,cxtimeavl,maxtavl,nxtavl,nxt,ldiffsuccess,
     &   maxslab,maxlev,maxpl,miy,mjx,mkzh,mabpl,morpl)
c
c   Compute fields to plot. When adding a new field, the first call to
c   getpt must be for the output field (pl2 or pl3), then call getpt
c   for the scratch arrays.
c
      dimension uuu(miy,mjx,mkzh), vvv(miy,mjx,mkzh),
     &   tmk(miy,mjx,mkzh), qvp(miy,mjx,mkzh), prs(miy,mjx,mkzh),
     &   ght(miy,mjx,mkzh), www(miy,mjx,mkzh), pstx(miy,mjx),
     &   dmap(miy,mjx), xmap(miy,mjx), ter(miy,mjx), pstd(miy,mjx),
     &   cor(miy,mjx), xlus(miy,mjx), prs_tsf(miy,mjx)
c
      dimension wk(miy,mjx,maxslab),indwk(3,maxpl),icdwk(maxpl),
     &   sigh(mkzh),pslab1(mabpl,morpl),
     &   unorth(miy,mjx),vnorth(miy,mjx),rstmv(2,maxpl),
     &   sigf(mkzh+1),rcrag(2,maxpl),rcrbg(2,maxpl),
     &   idimn(maxpl),igdir(maxpl),rsepa(32,maxpl),
     &   rlevl(maxlev,maxpl),rlavl(maxlev,maxpl),rrfst(4,maxpl),
     &   raddf(maxpl),rmply(maxpl),rfdmn(maxpl),rfdmx(maxpl),
     &   ismcp(maxpl),
     &   ismth(maxpl),iqgsm(maxpl),
     &   rtjst(maxpl),rtjen(maxpl),rdiff(maxpl),
     &   xtimeavl(maxtavl),xtimeavl_sv(maxtavl)
c
      character cfeld(3,maxpl)*10,cptyp(maxpl)*2,
     &   csave(maxpl)*10,ctjfl(maxpl)*256,engplttl(maxpl)*36,
     &   cdiff(maxpl)*256,unwk(maxpl)*24,fname*256,varname*10,
     &   casename*256,cxtimeavl(maxtavl)*9,fname_sv*256,
     &   casename_sv*256,cxtimeavl_sv(maxtavl)*9,
     &   cvcorunits*16
      logical lverf(maxpl),lredo(maxpl),lgrad(maxpl),
     &   llapl(maxpl),lhadv(maxpl),ldfrl(maxpl),ldiffsuccess
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32),fullsigma(128),halfsigma(128)
      character chrip(64)*64,vardesc*64,plchun*24
c
      character getfromcore(50)*10
c
      logical numeric
      external numeric
c
      include 'comconst'
c
      include 'pointers'
c
c   Following are fields that are already in core memory (i.e.,
c   in arrays by the same name) and therefore do not need to be
c   calculated or read from a file.
c
      data (getfromcore(i),i=1,15) /'uuu       ','vvv       ',
     &    'tmk       ','prs       ','qvp       ','ght       ',
     &    'ter       ','dmap      ','xmap      ','xlus      ',
     &    'cor       ','pstx      ','pstd      ','www       ',
     &    'prs_tsf   '/
c
      do 2000 ifld=1,3
c
      if (ifld.eq.2.and.
     &    (cptyp(ipl)(2:2).eq.'c'.or.cptyp(ipl)(2:2).eq.'h')) goto 2000
      if (ifld.eq.3.and.cptyp(ipl).ne.'vv') goto 2000
c
c   Don't read or calculate field if already in work array
c
      do iplchk=iplstrt,ipl-1
         if ((cfeld(ifld,ipl).eq.cfeld(ifld,iplchk)).and.
     &       (ismcp(ipl).eq.ismcp(iplchk)).and.
     &       (ismth(ipl).eq.ismth(iplchk)).and.
     &       (lgrad(ipl).eqv.lgrad(iplchk)).and.
     &       (llapl(ipl).eqv.llapl(iplchk)).and.
     &       (lhadv(ipl).eqv.lhadv(iplchk)).and.
     &       igdir(ipl).eq.igdir(iplchk).and.
     &       iqgsm(ipl).eq.iqgsm(iplchk).and.
     &       rdiff(ipl).eq.rdiff(iplchk).and.
     &       cdiff(ipl).eq.cdiff(iplchk)) then
            if (iplchk.gt.iplstrt.and.raddf(iplchk-1).ne.0..and.
     &          raddf(iplchk).eq.0.) goto 291
            if (iplchk.gt.iplstrt.and.rmply(iplchk-1).ne.0..and.
     &          rmply(iplchk).eq.0.) goto 291
            if (rsepa(1,ipl).eq.rsepa(1,iplchk).and.
     &          rsepa(2,ipl).eq.rsepa(2,iplchk).and.
     &          rsepa(3,ipl).eq.rsepa(3,iplchk).and.
     &          rsepa(4,ipl).eq.rsepa(4,iplchk).and.
     &          rsepa(5,ipl).eq.rsepa(5,iplchk).and.
     &          rsepa(6,ipl).eq.rsepa(6,iplchk).and.
     &          rsepa(7,ipl).eq.rsepa(7,iplchk).and.
     &          rsepa(8,ipl).eq.rsepa(8,iplchk).and.
     &          rsepa(9,ipl).eq.rsepa(9,iplchk)) then
               indwk(ifld,ipl)=indwk(ifld,iplchk)
               icdwk(ipl)=icdwk(iplchk)
               unwk(ipl)=unwk(iplchk)
               engplttl(ipl)=engplttl(iplchk)
               idimn(ipl)=3
               if (idimn(iplchk).eq.2) idimn(ipl)=2
               goto 1980
            endif
         endif
      enddo
 291  continue
c
c   Some preliminary stuff
c
      idiffpass=0
      idimn(ipl)=3
      engplttl(ipl)=' '
 35   ifree=incwk
      idiffpass=idiffpass+1
c
c   If the field is one that is not in core memory, it may be in
c   a file, either because it was processed by ripdp, or because
c   it was calculated by RIP and saved.  Therefore, for fields
c   not in core memory (i.e., not in the "getfromcore" list),
c   check to see if a file exists and if so, read it from the
c   file (unless the user has set "redo" because he/she specifically
c   wants the field to be re-calculated).
c
      if (lredo(ipl)) goto 9
      do i=1,50
         if (getfromcore(i).eq.'          ') goto 7
         if (cfeld(ifld,ipl).eq.getfromcore(i)) goto 9
      enddo
 7    continue
      fname(iendf1+1:)=cfeld(ifld,ipl)
      open (unit=iudata,err=9,file=fname,form='unformatted',
     &   status='old')
      write(iup,*)'    Reading field ',cfeld(ifld,ipl),' from a file.'
      read (iudata)
     &   vardesc,plchun,ihrip,rhrip,chrip,fullsigma,halfsigma
      idimn(ipl)=ihrip(6)
      indwk(ifld,ipl)=incwk
      icdwk(ipl)=ihrip(7)
      unwk(ipl)=plchun
      do i=35,1,-1
         if (vardesc(i:i+1).eq.', ') then
            iendv=i-1
            goto 71
         endif
      enddo
      if (vardesc(36:36).eq.',') then
         iendv=35
      else
         iendv=36
      endif
 71   continue
      engplttl(ipl)=vardesc(1:iendv)
      if (idimn(ipl).eq.2) then
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         read (iudata) pl2
      elseif (idimn(ipl).eq.3) then
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         read (iudata) pl3
      else
         write(iup,*)'Dimension of data is not 2 or 3.  Stopping.'
         stop
      endif
      close (iudata)
      goto 1970
 9    continue
c
c   Load or calculate new field.
c
      if (cfeld(ifld,ipl)(1:4).eq.'uuu ') then! u-velocity, m/s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call addorfill(uuu,pl3,miy,mjx,mkzh,3,0,1.,0.)
         if (rstmv(2,ipl).ne.0.0) then
            do k=1,mkzh
            do j=1,mjx
            do i=1,miy
               pl3(i,j,k)=uuu(i,j,k)-rstmv(2,ipl)
            enddo
            enddo
            enddo
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Horizontal wind (x-comp.)'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'vvv ') then! v-velocity, m/s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call addorfill(vvv,pl3,miy,mjx,mkzh,3,0,1.,0.)
         if (rstmv(1,ipl).ne.0.0) then
            do k=1,mkzh
            do j=1,mjx
            do i=1,miy
               pl3(i,j,k)=vvv(i,j,k)-rstmv(1,ipl)
            enddo
            enddo
            enddo
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Horizontal wind (y-comp.)'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tmc ') then! temperature, deg C
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do k=1,mkzh
         do j=1,mjx
         do i=1,miy
            pl3(i,j,k)=tmk(i,j,k)-celkel
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Temperature'
         unwk(ipl)='~S~o~N~C'
      elseif (cfeld(ifld,ipl)(1:7).eq.'tptslr ') then
c        Temperature pert. from standard lapse rate, K (or deg. C)
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do k=1,mkzh
         do j=1,mjx
         do i=1,miy
            trefer=celkel+max(-55.,15.-.0065*ght(i,j,k))
            pl3(i,j,k)=tmk(i,j,k)-trefer
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Temperature Perturb. (from USSALR)'
         unwk(ipl)='K (or ~S~o~N~C)'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tmk ') then! temperature, K
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call addorfill(tmk,pl3,miy,mjx,mkzh,3,1,1.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Temperature'
         unwk(ipl)='K'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tmf ') then! temperature, deg F
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do k=1,mkzh
         do j=1,mjx
         do i=1,miy
            pl3(i,j,k)=(tmk(i,j,k)-celkel)*1.8 + 32.0
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Temperature'
         unwk(ipl)='~S~o~N~F'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tpt ') then! tmp. prtrb., K (or deg C)
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         alnslp=log(.01*rrfst(1,ipl))
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            alnp=log(prs(i,j,k))
            reftmk=rrfst(2,ipl)+rrfst(3,ipl)*(alnp-alnslp)
            reftmk=max(reftmk,rrfst(4,ipl))
            pl3(i,j,k)=tmk(i,j,k)-reftmk
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Temperature pert. (from MM5 ref. state)'
         unwk(ipl)='K (or ~S~o~N~C)'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tdp ') then! dewpoint, deg C
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call tdpcalc(qvp,prs,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Dewpoint temperature'
         unwk(ipl)='~S~o~N~C'
      elseif (cfeld(ifld,ipl)(1:4).eq.'twb ') then! wet-blb tmp., deg C
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call wetbulbcalc(prs,tmk,qvp,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Wetbulb temperature'
         unwk(ipl)='~S~o~N~C'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tvk ') then! virtual temp., K
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do k = 1, mkzh
         do j = 1, mjx-1
         do i = 1, miy-1
            pl3(i,j,k)=virtual(tmk(i,j,k),qvp(i,j,k))
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Virtual temperature'
         unwk(ipl)='K'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tvc ') then! virtual temp., deg. C
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do k = 1, mkzh
         do j = 1, mjx-1
         do i = 1, miy-1
            pl3(i,j,k)=virtual(tmk(i,j,k),qvp(i,j,k))-celkel
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Virtual temperature'
         unwk(ipl)='~S~o~N~C'
      elseif (cfeld(ifld,ipl)(1:4).eq.'qvp ') then! vapor mix rat, g/kg
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call addorfill(qvp,pl3,miy,mjx,mkzh,3,1,1000.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Water vapor mixing ratio'
         unwk(ipl)='g kg~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:3).eq.'qra') then! rain mix rat, g/kg
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call readdat(iudata,fname,iendf1,'qra       ',
     &      miy,mjx,mkzh,3,1,pl3,istat)
         if (iice.eq.0.and.cfeld(ifld,ipl)(4:4).ne.'b') then
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (tmk(i,j,k).lt.celkel) pl3(i,j,k)=0.
            enddo
            enddo
            enddo
         elseif (iice.eq.0.and.cfeld(ifld,ipl)(4:4).eq.'b') then
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
            call wetbulbcalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
            call precprob(tmk,scr3a,ght,ter,scr3b,miy,mjx,mkzh)
            do j=1,mjx-1
            do i=1,miy-1
               probsnow=scr3b(i,j,3)
            do k=1,mkzh
               if (probsnow.ge.50..or.(probsnow.lt.50..and.
     &             k.ne.mkzh.and.tmk(i,j,k).lt.celkel)) 
     &            pl3(i,j,k)=0.
            enddo
            enddo
            enddo
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Rain water mixing ratio'
         unwk(ipl)='g kg~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:3).eq.'qsn') then! snow mix rat, g/kg
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         if (iice.eq.1) then
            call readdat(iudata,fname,iendf1,'qsn       ',
     &         miy,mjx,mkzh,3,1,pl3,istat)
         elseif (iice.eq.0.and.cfeld(ifld,ipl)(4:4).ne.'b') then
            call readdat(iudata,fname,iendf1,'qra       ',
     &         miy,mjx,mkzh,3,1,pl3,istat)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (tmk(i,j,k).ge.celkel) pl3(i,j,k)=0.
            enddo
            enddo
            enddo
         elseif (iice.eq.0.and.cfeld(ifld,ipl)(4:4).eq.'b') then
            call readdat(iudata,fname,iendf1,'qra       ',
     &         miy,mjx,mkzh,3,1,pl3,istat)
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
            call wetbulbcalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
            call precprob(tmk,scr3a,ght,ter,scr3b,miy,mjx,mkzh)
            do j=1,mjx-1
            do i=1,miy-1
               probsnow=scr3b(i,j,3)
            do k=1,mkzh
               if (.not.(probsnow.ge.50..or.(probsnow.lt.50..and.
     &             k.ne.mkzh.and.tmk(i,j,k).lt.celkel))) 
     &            pl3(i,j,k)=0.
            enddo
            enddo
            enddo
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Snow mixing ratio'
         unwk(ipl)='g kg~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'qpr ') then! qsn+qra+qgra, g/kg
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call readdat(iudata,fname,iendf1,'qra       ',
     &      miy,mjx,mkzh,3,1,pl3,istat)
         if (iice.eq.1) then
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
            call readdat(iudata,fname,iendf1,'qsn       ',
     &         miy,mjx,mkzh,3,1,scr3a,istat)
            call addorfill(scr3a,pl3,miy,mjx,mkzh,3,1,1.,1.)
            call readdat(iudata,fname,iendf1,'qgr       ',
     &         miy,mjx,mkzh,3,0,scr3a,istat)
            if (istat.ge.0)
     &         call addorfill(scr3a,pl3,miy,mjx,mkzh,3,1,1.,1.)
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Total precipitation mixing ratio'
         unwk(ipl)='g kg~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'snf ') then! % snow=(qsn+qgra)/qpr
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call readdat(iudata,fname,iendf1,'qra       ',
     &      miy,mjx,mkzh,3,1,pl3,istat)
         if (iice.eq.0) then
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (pl3(i,j,k).le..0001e-3) then
                  pl3(i,j,k)=rmsg
               elseif (tmk(i,j,k).ge.celkel) then
                  pl3(i,j,k)=0.
               else
                  pl3(i,j,k)=100.
               endif
            enddo
            enddo
            enddo
         elseif (iice.eq.1) then
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
            call readdat(iudata,fname,iendf1,'qsn       ',
     &         miy,mjx,mkzh,3,1,scr3a,istat)
            call readdat(iudata,fname,iendf1,'qgr       ',
     &         miy,mjx,mkzh,3,0,scr3b,istat)
            if (istat.ge.0)
     &         call addorfill(scr3b,scr3a,miy,mjx,mkzh,3,1,1.,1.)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               qtot=pl3(i,j,k)+scr3a(i,j,k)
               if (qtot.le..0001e-3) then
                  pl3(i,j,k)=rmsg
               else
                  pl3(i,j,k)=100.*scr3a(i,j,k)/qtot
               endif
            enddo
            enddo
            enddo
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Frozen precipitation fraction'
         unwk(ipl)='%'
      elseif (cfeld(ifld,ipl)(1:4).eq.'dbz ') then! reflectivity, dBZ
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call readdat(iudata,fname,iendf1,'qra       ',
     &      miy,mjx,mkzh,3,1,pl3,istat)
         if (iice.eq.1) then
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
            call readdat(iudata,fname,iendf1,'qsn       ',
     &         miy,mjx,mkzh,3,1,scr3a,istat)
            call readdat(iudata,fname,iendf1,'qgr       ',
     &         miy,mjx,mkzh,3,0,scr3b,istat)
         endif
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            dens=prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k)))      ! density
            if (iice.eq.0) then
               if (tmk(i,j,k).ge.celkel) then
                  pl3(i,j,k)=((.001*pl3(i,j,k)*dens)**1.75)*
     $               3.630803e-9 * 1.e18                  ! Z for rain
               else
                  pl3(i,j,k)=((.001*pl3(i,j,k)*dens)**1.75)*
     $               2.18500e-10 * 1.e18                  ! Z for snow
               endif
            elseif (iice.eq.1) then
               pl3(i,j,k)=((.001*pl3(i,j,k)*dens)**1.75)*
     $               3.630803e-9 * 1.e18                  ! Z for rain
               pl3(i,j,k)= pl3(i,j,k)+
     $            ((.001*scr3a(i,j,k)*dens)**1.75)*
     $               2.18500e-10 * 1.e18                  ! Z for snow
               if (istat.ge.0) pl3(i,j,k)= pl3(i,j,k)+
     $            ((.001*scr3b(i,j,k)*dens)**1.75)*
     $               1.033267e-9 * 1.e18                  ! Z for graup
            endif
            if (pl3(i,j,k).ne.0.)
     $         pl3(i,j,k)=10*log10(pl3(i,j,k)) ! dBZ
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Reflectivity'
         unwk(ipl)='dBZ'
      elseif (cfeld(ifld,ipl)(1:5).eq.'pcpw '.or.  ! precip'bl water, mm
     &        cfeld(ifld,ipl)(1:7).eq.'intcld '.or. ! integ. cld. hyds., mm
     &        cfeld(ifld,ipl)(1:7).eq.'intpcp ') then! integ. pcp. hyds., mm
         if (cfeld(ifld,ipl)(1:5).eq.'pcpw ') then
            iv=1
            ic=1
            ip=1
            engplttl(ipl)='Precipitable water'
         elseif (cfeld(ifld,ipl)(1:7).eq.'intcld ') then
            iv=0
            ic=1
            ip=0
            engplttl(ipl)='Column-integ. cloud hydrometeors'
         elseif (cfeld(ifld,ipl)(1:7).eq.'intpcp ') then
            iv=0
            ic=0
            ip=1
            engplttl(ipl)='Column-integ. precip. hydrometeors'
         endif
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call fillarray(scr3a,miy*mjx*mkzh,0.)
         call fillarray(scr3b,miy*mjx*mkzh,0.)
         call addorfill(qvp,scr3a,miy,mjx,mkzh,3,1,1.,1.)
         if (iv.eq.1) call addorfill(qvp,scr3b,miy,mjx,mkzh,
     &                               3,1,1.,1.)
         call readdat(iudata,fname,iendf1,'qra       ',
     &      miy,mjx,mkzh,3,1,scr3c,istat)
         call addorfill(scr3c,scr3a,miy,mjx,mkzh,3,1,.001,1.)
         if (ip.eq.1) call addorfill(scr3c,scr3b,miy,mjx,mkzh,
     &                               3,1,.001,1.)
         call readdat(iudata,fname,iendf1,'qcw       ',
     &      miy,mjx,mkzh,3,1,scr3c,istat)
         call addorfill(scr3c,scr3a,miy,mjx,mkzh,3,1,.001,1.)
         if (ic.eq.1) call addorfill(scr3c,scr3b,miy,mjx,mkzh,
     &                               3,1,.001,1.)
         if (iice.eq.1) then
            call readdat(iudata,fname,iendf1,'qci       ',
     &         miy,mjx,mkzh,3,1,scr3c,istat)
            call addorfill(scr3c,scr3a,miy,mjx,mkzh,3,1,.001,1.)
            if (ic.eq.1) call addorfill(scr3c,scr3b,miy,mjx,mkzh,
     &                                  3,1,.001,1.)
            call readdat(iudata,fname,iendf1,'qsn       ',
     &         miy,mjx,mkzh,3,1,scr3c,istat)
            call addorfill(scr3c,scr3a,miy,mjx,mkzh,3,1,.001,1.)
            if (ip.eq.1) call addorfill(scr3c,scr3b,miy,mjx,mkzh,
     &                                  3,1,.001,1.)
            call readdat(iudata,fname,iendf1,'qgr       ',
     &         miy,mjx,mkzh,3,0,scr3c,istat)
            if (istat.ge.0) then
               call addorfill(scr3c,scr3a,miy,mjx,mkzh,3,1,.001,1.)
               if (ip.eq.1) call addorfill(scr3c,scr3b,miy,mjx,mkzh,
     &                                     3,1,.001,1.)
            endif
         endif
         call pfcalc(sigh,sigf,pstx,prs,scr3c,miy,mjx,mkzh)
         do j=1,mjx-1
         do i=1,miy-1
            pl2(i,j)=0.
            do k=1,mkzh
               if (k.eq.1) then
                  pabove=prs(i,j,1)-sigh(1)*pstx(i,j)
               else
                  pabove=scr3c(i,j,k-1)
               endif
               pl2(i,j)=pl2(i,j)+scr3b(i,j,k)/(1.+scr3a(i,j,k))*
     &            100.*(scr3c(i,j,k)-pabove) !mix rats are kg/kg, prs's are hPa
            enddo
            pl2(i,j)=pl2(i,j)/grav  ! dens of water cancels with conv to mm
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='mm'
c      elseif (cfeld(ifld,ipl)(1:9).eq.'qnz6log10') then! nuc.rate,cm^3/s,log10
c         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
c         call readdat(iudata,fname,iendf1,'qnz6      ',
c     &      miy,mjx,mkzh,3,1,pl3,istat)
c         do k=1,mkzh
c         do j=1,mjx-1
c         do i=1,miy-1
c            pl3(i,j,k)=alog10(pl3(i,j,k))
c         enddo
c         enddo
c         enddo
c         indwk(ifld,ipl)=incwk
c         icdwk(ipl)=1
c         engplttl(ipl)='Nuc. rate'
c         unwk(ipl)='cm~S~3~N~ s~S~-1~N~, log'
      elseif (cfeld(ifld,ipl)(1:3).eq.'qcw') then! cld wat mix r., g/kg
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call readdat(iudata,fname,iendf1,'qcw       ',
     &      miy,mjx,mkzh,3,1,pl3,istat)
         if (iice.eq.0.and.cfeld(ifld,ipl)(4:4).ne.'b') then
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (tmk(i,j,k).lt.celkel) pl3(i,j,k)=0.
            enddo
            enddo
            enddo
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Cloud water mixing ratio'
         unwk(ipl)='g kg~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:3).eq.'qci') then! cld ice mix r., g/kg
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         if (iice.eq.1) then
            call readdat(iudata,fname,iendf1,'qci       ',
     &         miy,mjx,mkzh,3,1,pl3,istat)
         elseif (iice.eq.0) then
            call readdat(iudata,fname,iendf1,'qcw       ',
     &         miy,mjx,mkzh,3,1,pl3,istat)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (tmk(i,j,k).ge.celkel) pl3(i,j,k)=0.
            enddo
            enddo
            enddo
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Cloud ice mixing ratio'
         unwk(ipl)='g kg~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'qcl ') then! qcw+qci, g/kg
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call readdat(iudata,fname,iendf1,'qcw       ',
     &      miy,mjx,mkzh,3,1,pl3,istat)
         if (iice.eq.1) then
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
            call readdat(iudata,fname,iendf1,'qci       ',
     &         miy,mjx,mkzh,3,1,scr3a,istat)
            call addorfill(scr3a,pl3,miy,mjx,mkzh,3,1,1.,1.)
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Total cloud mixing ratio'
         unwk(ipl)='g kg~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'ght ') then! geop. height, m
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call addorfill(ght,pl3,miy,mjx,mkzh,3,1,1.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Geopotential height'
         unwk(ipl)='m'
      elseif (cfeld(ifld,ipl)(1:4).eq.'sig ') then! sigma
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=sigh(k)
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Sigma'
         unwk(ipl)='none'
      elseif (cfeld(ifld,ipl)(1:4).eq.'prs ') then! pressure, hPa
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call addorfill(prs,pl3,miy,mjx,mkzh,3,1,1.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Pressure'
         unwk(ipl)='hPa'
      elseif (cfeld(ifld,ipl)(1:4).eq.'rho ') then! density, kg/m**3
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=prs(i,j,k)*100./
     &         (rgas*virtual(tmk(i,j,k),qvp(i,j,k)))
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Density'
         unwk(ipl)='kg m~S~-3~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'vacc') then! vert. acc., m/s**2
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
c
         if (inhyd.eq.0) then
            call fillarray(pl3,miy*mjx*mkzh,0.)
         else
            call readdat(iudata,fname,iendf1,'qra       ',
     &         miy,mjx,mkzh,3,1,pl3,istat)
            call readdat(iudata,fname,iendf1,'qcw       ',
     &         miy,mjx,mkzh,3,1,scr3a,istat)
            call addorfill(scr3a,pl3,miy,mjx,mkzh,
     &         3,1,1.,1.)
            if (iice.eq.1) then
               call readdat(iudata,fname,iendf1,'qsn       ',
     &            miy,mjx,mkzh,3,1,scr3a,istat)
               call addorfill(scr3a,pl3,miy,mjx,mkzh,
     &            3,1,1.,1.)
               call readdat(iudata,fname,iendf1,'qci       ',
     &            miy,mjx,mkzh,3,1,scr3a,istat)
               call addorfill(scr3a,pl3,miy,mjx,mkzh,
     &            3,1,1.,1.)
               call readdat(iudata,fname,iendf1,'qgr       ',
     &            miy,mjx,mkzh,3,0,scr3a,istat)
               if (istat.ge.0) call addorfill(scr3a,pl3,miy,mjx,mkzh,
     &            3,1,1.,1.)
            endif
            do k=1,mkzh
               kp1=min(k+1,mkzh)
               km1=max(k-1,1)
            do j=1,mjx-1
            do i=1,miy-1
               hydrometmr=.001*pl3(i,j,k) ! g/kg to kg/kg
               refprs=sigh(k)*pstx(i,j)+ptop
               reftmk=refslt+reflaps*alog(refprs/(.01*refslp))
               reftmk=max(reftmk,refstratt)
               refrho=refprs/rgas/reftmk
               rho=prs(i,j,k)/(rgas*virtual(tmk(i,j,k),qvp(i,j,k)))
               pl3(i,j,k)=grav*refrho/rho*(1./pstx(i,j)*
     &            ((prs(i,j,kp1)-prs(i,j,km1))/
     &            (sigh(kp1)-sigh(km1))-pstx(i,j))+
     &            (virtual(tmk(i,j,k),qvp(i,j,k))-reftmk)/tmk(i,j,k)-
     &            reftmk*(prs(i,j,k)-refprs)/(tmk(i,j,k)*refprs))-
     &            grav*hydrometmr
            enddo
            enddo
            enddo
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Vertical acceleration'
         unwk(ipl)='m s~S~-2~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'ppt ') then! press pertbn, hPa
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         cc1=rgas/grav*(-.5)*rrfst(3,ipl)
         cc2=rgas/grav*(rrfst(3,ipl)*log(.01*rrfst(1,ipl))-
     &      rrfst(2,ipl))
         cc3=rgas/grav*(rrfst(2,ipl)-.5*rrfst(3,ipl)*
     &      log(.01*rrfst(1,ipl)))*log(.01*rrfst(1,ipl))
         do j=1,mjx-1
         do i=1,miy-1
         do k=1,mkzh
            pl3(i,j,k)=prs(i,j,k)-exp(
     &         (-cc2-sqrt(cc2*cc2-4.*cc1*(cc3-ght(i,j,k))))/(2.*cc1) )
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Pressure pert. (from MM5 ref. state)'
         unwk(ipl)='hPa'
      elseif (cfeld(ifld,ipl)(1:4).eq.'the ') then! pot. temperature, K
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call thecalc(prs,tmk,qvp,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Potential temperature'
         unwk(ipl)='K'
      elseif (cfeld(ifld,ipl)(1:9).eq.'the_sfan ') then! sfc anal. of theta, K
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call readdat(iudata,fname,iendf1,'tmk_sfan  ',
     &      miy,mjx,mkzh,2,1,pl2,istat)
         call readdat(iudata,fname,iendf1,'qvp_sfan  ',
     &      miy,mjx,mkzh,2,1,scr2a,istat)
         do j=1,mjx-1
         do i=1,miy-1
            gammam=gamma*(1.+gammamd*.001*scr2a(i,j))
            pl2(i,j)=pl2(i,j)*(1000./prs_tsf(i,j))**gammam
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Surface potential temperature'
         unwk(ipl)='K'
      elseif (cfeld(ifld,ipl)(1:4).eq.'eth ') then! eqv. pot. temp., K
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call eqthecalc(qvp,tmk,prs,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Equivalent potential temperature'
         unwk(ipl)='K'
      elseif (cfeld(ifld,ipl)(1:4).eq.'rcum'.or.
     &        cfeld(ifld,ipl)(1:4).eq.'rexp'.or.
     &        cfeld(ifld,ipl)(1:4).eq.'rtot') then  ! rainfall, mm
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         if (cfeld(ifld,ipl)(1:4).eq.'rcum') then
            icum=1
            iexp=0
            engplttl(ipl)='Cumulus precip.'
            iendeng=15
         elseif (cfeld(ifld,ipl)(1:4).eq.'rexp') then
            icum=0
            iexp=1
            engplttl(ipl)='Explicit precip.'
            iendeng=16
         elseif (cfeld(ifld,ipl)(1:4).eq.'rtot') then
            icum=1
            iexp=1
            engplttl(ipl)='Total precip.'
            iendeng=13
         endif
         iendcf=4
         if (cfeld(ifld,ipl)(iendcf+1:iendcf+2).eq.'sh') then
c
c         Given time is interpretted as "since hour x"
c
            irel=0
            ispos=iendcf+3
            iepos=index(cfeld(ifld,ipl),' ')-1
            if (iepos.eq.-1) iepos=10
         else
c
c         Given time is interpretted as "in past x hours"
c
            irel=1
            ispos=iendcf+1
            iepos=iendcf+index(cfeld(ifld,ipl)(iendcf+1:),'h')-1
         endif
         if (.not.numeric(cfeld(ifld,ipl)(ispos:iepos))) then
            write(iup,*)'Processing rainfall specifier in fields.f.'
            write(iup,*)'Non-numeric characters found where time (in h)'
            write(iup,*)'is supposed to be specified.'
            write(iup,*)'feld was given as ',cfeld(ifld,ipl)
            stop
         endif
         read(cfeld(ifld,ipl)(ispos:iepos),fmt=*)raintime
         if (irel.eq.0) then
            xtime_get=raintime
            engplttl(ipl)(iendeng+1:)=' since h '//
     &         cfeld(ifld,ipl)(ispos:iepos)
         else
            xtime_get=max(0.0,xtime-raintime)
            engplttl(ipl)(iendeng+1:)=' in past '//
     &         cfeld(ifld,ipl)(ispos:iepos)//' h'
         endif
c
c      Get current rainfall
c
         do j=1,mjx-1
         do i=1,miy-1
            pl2(i,j)=0.
         enddo
         enddo
         if (icum.eq.1) then
            call readdat(iudata,fname,iendf1,'rtc       ',
     &         miy,mjx,mkzh,2,1,scr2a,istat)
            call addorfill(scr2a,pl2,miy,mjx,mkzh,2,1,1.,1.)
         endif
         if (iexp.eq.1) then
            call readdat(iudata,fname,iendf1,'rte       ',
     &         miy,mjx,mkzh,2,1,scr2a,istat)
            call addorfill(scr2a,pl2,miy,mjx,mkzh,2,1,1.,1.)
         endif
c
c      Subtract previous rainfall if xtime_get is greater than ~.1 seconds
c
         if (xtime_get.ge.3e-5) then
c
c         Check if requested time is available
c
            iavail=0
            do i=1,nxtavl
               if (abs(xtimeavl(i)-xtime_get).le.tacch) then
                  iavail=1
                  nxt_get=i
                  goto 57
               endif
            enddo
 57         continue
            if (iavail.eq.0) then
               write(iup,*)'   Requested previous time ',
     &            xtime_get,' for'
               write(iup,*)'   acccumulated rainfall is not available.'
               write(iup,*)'   RIP will plot requested rainfall type'
               write(iup,*)'   since h 0.'
               goto 59
            endif
c
c         Build file name for previous time
c
            fname=casename(1:iendc)//'_'//cxtimeavl(nxt_get)//'_'
c
            if (icum.eq.1) then
               call readdat(iudata,fname,iendf1,'rtc       ',
     &            miy,mjx,mkzh,2,1,scr2a,istat)
               call addorfill(scr2a,pl2,miy,mjx,mkzh,2,1,-1.,1.)
            endif
            if (iexp.eq.1) then
               call readdat(iudata,fname,iendf1,'rte       ',
     &            miy,mjx,mkzh,2,1,scr2a,istat)
               call addorfill(scr2a,pl2,miy,mjx,mkzh,2,1,-1.,1.)
            endif
            fname=casename(1:iendc)//'_'//cxtimeavl(nxt)//'_'
         endif
 59      continue
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='mm'
      elseif (cfeld(ifld,ipl)(1:4).eq.'ter ') then! terrain, m
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         if (iprog.eq.2.or.iprog.eq.3) then
            call readdat(iudata,fname,iendf1,'ter_tsf   ',
     &         miy,mjx,mkzh,2,1,pl2,istat)
         else
            call addorfill(ter,pl2,miy,mjx,mkzh,2,1,1.,0.)
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Terrain height AMSL'
         unwk(ipl)='m'
      elseif (cfeld(ifld,ipl).eq.'ter_pseudo') then! height of pseudo-terrain
c                                                  ! for p-level data, m
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call addorfill(ter,pl2,miy,mjx,mkzh,2,1,1.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Pseudo-terrain height AMSL'
         unwk(ipl)='m'
      elseif (cfeld(ifld,ipl)(1:5).eq.'xmap ') then! map fac (x-points)
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call addorfill(xmap,pl2,miy,mjx,mkzh,2,1,1.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Map factor on cross points'
         unwk(ipl)='none'
      elseif (cfeld(ifld,ipl)(1:5).eq.'dmap ') then! map fac (dot points)
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call addorfill(dmap,pl2,miy,mjx,mkzh,2,0,1.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Map factor on dot points'
         unwk(ipl)='none'
      elseif (cfeld(ifld,ipl)(1:4).eq.'cor ') then! Coriolis, per s
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call addorfill(cor,pl2,miy,mjx,mkzh,2,1,1.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Coriolis parameter'
         unwk(ipl)='s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:5).eq.'xlus ') then! land use
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call addorfill(xlus,pl2,miy,mjx,mkzh,2,1,1.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Land use category'
         unwk(ipl)='none'
      elseif (cfeld(ifld,ipl)(1:4).eq.'sfp '.or.
     &        cfeld(ifld,ipl)(1:8).eq.'prs_tsf ') then ! prs at true sfc, hPa
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         do j=1,mjx-1
         do i=1,miy-1
            pl2(i,j)=prs_tsf(i,j)
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         if (iprog.le.3) then
            engplttl(ipl)='Pressure at true surface'
         else
            engplttl(ipl)='Surface pressure'
         endif
         unwk(ipl)='hPa'
      elseif (cfeld(ifld,ipl)(1:5).eq.'pstx ') then! p-star (x-points), hPa
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call addorfill(pstx,pl2,miy,mjx,mkzh,2,1,1.,0.)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='P-star on cross points'
         unwk(ipl)='hPa'
      elseif (cfeld(ifld,ipl)(1:5).eq.'pstd ') then! p-star (dot-points), hPa
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call addorfill(pstx,pl2,miy,mjx,mkzh,2,0,1.,0.)
         call xtodot(pl2,miy,mjx)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='P-star on dot points'
         unwk(ipl)='hPa'
      elseif (cfeld(ifld,ipl)(1:4).eq.'frgm'.or.
     &        cfeld(ifld,ipl)(1:4).eq.'frem') then
c
c      Frontogenesis (Miller form), K per 100 km per day
c
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3e,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3f,wk,maxslab)
c
c      Put theta or theta-e in scr3a
c
         if (cfeld(ifld,ipl)(3:3).eq.'g') then
            call thecalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
         elseif (cfeld(ifld,ipl)(3:3).eq.'e') then
            call eqthecalc(qvp,tmk,prs,scr3a,miy,mjx,mkzh)
         endif
c
         call derivcz(scr3a,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
         call derivcz(scr3a,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,scr3c,1,'y',miy,mjx,mkzh)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3d(i,j,k)=sqrt(scr3b(i,j,k)*scr3b(i,j,k)+
     &         scr3c(i,j,k)*scr3c(i,j,k))
         enddo
         enddo
         enddo
         if (cfeld(ifld,ipl)(5:7).eq.'had') then
            call derivcz(scr3d,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3e,1,'x',miy,mjx,mkzh)
            call derivcz(scr3d,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3f,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               ucross=.25*(uuu(i,j,k)+uuu(i+1,j,k)+uuu(i,j+1,k)+
     &                     uuu(i+1,j+1,k))
               vcross=.25*(vvv(i,j,k)+vvv(i+1,j,k)+vvv(i,j+1,k)+
     &                     vvv(i+1,j+1,k))
               pl3(i,j,k)=-(ucross*scr3e(i,j,k)+vcross*scr3f(i,j,k))
            enddo
            enddo
            enddo
            engplttl(ipl)='Miller frontogenesis (hor. adv.)'
         elseif (cfeld(ifld,ipl)(5:7).eq.'vad') then
            do k=1,mkzh
               kp1=min(mkzh,k+1)
               km1=max(1,k-1)
            do j=1,mjx-1
            do i=1,miy-1
               pl3(i,j,k)=
     &            -.01*www(i,j,k)*(scr3d(i,j,kp1)-scr3d(i,j,km1))/
     &            (ght(i,j,kp1)-ght(i,j,km1))
            enddo
            enddo
            enddo
            engplttl(ipl)='Miller frontogenesis (ver. adv.)'
         elseif (cfeld(ifld,ipl)(5:7).eq.'div') then
            call derivcz(uuu,0,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3e,1,'x',miy,mjx,mkzh)
            call derivcz(vvv,0,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3f,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (scr3d(i,j,k).gt.0.) then
                  pl3(i,j,k)=-.5/scr3d(i,j,k)*(scr3b(i,j,k)*
     &               scr3b(i,j,k)+scr3c(i,j,k)*scr3c(i,j,k))*
     &               (scr3e(i,j,k)+scr3f(i,j,k))
               else
                  pl3(i,j,k)=0.
               endif
            enddo
            enddo
            enddo
            engplttl(ipl)='Miller frontogenesis (div. term)'
         elseif (cfeld(ifld,ipl)(5:7).eq.'def') then
            call derivcz(uuu,0,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3e,1,'x',miy,mjx,mkzh)
            call derivcz(vvv,0,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3f,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (scr3d(i,j,k).gt.0.) then
                  pl3(i,j,k)=-.5/scr3d(i,j,k)*(scr3b(i,j,k)*
     &               scr3b(i,j,k)-scr3c(i,j,k)*scr3c(i,j,k))*
     &               (scr3e(i,j,k)-scr3f(i,j,k))
               else
                  pl3(i,j,k)=0.
               endif
            enddo
            enddo
            enddo
            call derivcz(uuu,0,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3e,1,'y',miy,mjx,mkzh)
            call derivcz(vvv,0,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3f,1,'x',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (scr3d(i,j,k).gt.0.) then
                  pl3(i,j,k)=pl3(i,j,k)-1./scr3d(i,j,k)*scr3b(i,j,k)*
     &               scr3c(i,j,k)*(scr3e(i,j,k)+scr3f(i,j,k))
               else
                  pl3(i,j,k)=0.
               endif
            enddo
            enddo
            enddo
            engplttl(ipl)='Miller frontogenesis (def. term)'
         elseif (cfeld(ifld,ipl)(5:7).eq.'dia') then
            call vervcalc(pstd,dmap,sigf,xmap,pstx,sigh,
     &         uuu,vvv,qvp,tmk,www,prs,scr3e,2,miy,mjx,mkzh)
            call condheat(tmk,qvp,scr3e,0,prs,miy,mjx,mkzh,pl3)
            call derivcz(pl3,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3e,1,'x',miy,mjx,mkzh)
            call derivcz(pl3,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3f,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               if (scr3d(i,j,k).gt.0.) then
                  pl3(i,j,k)=1./scr3d(i,j,k)*2.778e-4*(scr3b(i,j,k)*
     &               scr3e(i,j,k)+scr3c(i,j,k)*scr3f(i,j,k))
               else
                  pl3(i,j,k)=0.
               endif
            enddo
            enddo
            enddo
            engplttl(ipl)='Miller frontogenesis (dia. term)'
         elseif (cfeld(ifld,ipl)(5:7).eq.'til') then
            call derivcz(www,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3e,1,'x',miy,mjx,mkzh)
            call derivcz(www,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3f,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
               kp1=min(mkzh,k+1)
               km1=max(1,k-1)
            do j=1,mjx-1
            do i=1,miy-1
               if (scr3d(i,j,k).gt.0.) then
                  pl3(i,j,k)=-1./scr3d(i,j,k)*.01*(scr3a(i,j,kp1)-
     &               scr3a(i,j,km1))/(ght(i,j,kp1)-ght(i,j,km1))*
     &               (scr3b(i,j,k)*scr3e(i,j,k)+
     &                scr3c(i,j,k)*scr3f(i,j,k))
               else
                  pl3(i,j,k)=0.
               endif
            enddo
            enddo
            enddo
            engplttl(ipl)='Miller frontogenesis (til. term)'
         endif
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=pl3(i,j,k)*3.6e8 ! K/m/s to K/100km/hour
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='K (100 km h)~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:5).eq.'frg2d'.or.
     &        cfeld(ifld,ipl)(1:5).eq.'fre2d') then
c
c      Frontogenesis (uni-directional form), K per 100 km per day
c
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3e,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3f,wk,maxslab)
c
c      Put theta or theta-e in scr3a
c
         if (cfeld(ifld,ipl)(3:3).eq.'g') then
            call thecalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
         elseif (cfeld(ifld,ipl)(3:3).eq.'e') then
            call eqthecalc(qvp,tmk,prs,scr3a,miy,mjx,mkzh)
         endif
         cosa=xdist/xseclen
         sina=ydist/xseclen
c
c      For some terms, put thermal gradient into scr3a instead
c
         if (cfeld(ifld,ipl)(6:8).ne.'dia'.and.
     &       cfeld(ifld,ipl)(6:8).ne.'dii'.and.
     &       cfeld(ifld,ipl)(6:8).ne.'til') then
c
         call derivcz(scr3a,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
         call derivcz(scr3a,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,scr3c,1,'y',miy,mjx,mkzh)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            if (cfeld(ifld,ipl)(6:8).ne.'shr') then
c              thermal grad. along cross-sec:
               scr3a(i,j,k)=cosa*scr3b(i,j,k)+sina*scr3c(i,j,k)
            else
c              thermal grad. into cross-sec:
               scr3a(i,j,k)=sina*scr3b(i,j,k)-cosa*scr3c(i,j,k)
            endif
         enddo
         enddo
         enddo
c
         endif
c
c      Smooth scr3a (on constant pressure), since it gets used in
c         almost every term.
c
         call smoothcp(scr3a,1,scr3b,prs,pslab1,iqgsm(ipl),
     &      miy,mjx,mkzh,mabpl,morpl)
c
         if (cfeld(ifld,ipl)(6:8).eq.'aad') then
            call derivcz(scr3a,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
            call derivcz(scr3a,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3c,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               scr3d(i,j,k)=.25*(uuu(i,j,k)+uuu(i+1,j,k)+uuu(i,j+1,k)+
     &                     uuu(i+1,j+1,k))-rstmv(2,ipl)
               scr3e(i,j,k)=.25*(vvv(i,j,k)+vvv(i+1,j,k)+vvv(i,j+1,k)+
     &                     vvv(i+1,j+1,k))-rstmv(1,ipl)
            enddo
            enddo
            enddo
            call smoothcp(scr3d,1,scr3a,prs,pslab1,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            call smoothcp(scr3e,1,scr3a,prs,pslab1,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               pl3(i,j,k)=-(cosa*scr3d(i,j,k)+sina*scr3e(i,j,k))*
     &                     (cosa*scr3b(i,j,k)+sina*scr3c(i,j,k))
            enddo
            enddo
            enddo
            engplttl(ipl)='2D frontogenesis (adv. along c.s.)'
         elseif (cfeld(ifld,ipl)(6:8).eq.'iad') then
            call derivcz(scr3a,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
            call derivcz(scr3a,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3c,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               scr3d(i,j,k)=.25*(uuu(i,j,k)+uuu(i+1,j,k)+uuu(i,j+1,k)+
     &                     uuu(i+1,j+1,k))-rstmv(2,ipl)
               scr3e(i,j,k)=.25*(vvv(i,j,k)+vvv(i+1,j,k)+vvv(i,j+1,k)+
     &                     vvv(i+1,j+1,k))-rstmv(1,ipl)
            enddo
            enddo
            enddo
            call smoothcp(scr3d,1,scr3a,prs,pslab1,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            call smoothcp(scr3e,1,scr3a,prs,pslab1,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               pl3(i,j,k)=-(sina*scr3d(i,j,k)-cosa*scr3e(i,j,k))*
     &                     (sina*scr3b(i,j,k)-cosa*scr3c(i,j,k))
            enddo
            enddo
            enddo
            engplttl(ipl)='2D frontogenesis (adv. into c.s.)'
         elseif (cfeld(ifld,ipl)(6:8).eq.'vad') then
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               scr3b(i,j,k)=www(i,j,k)
            enddo
            enddo
            enddo
            call smoothcp(scr3b,1,scr3c,prs,pslab1,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            do k=1,mkzh
               kp1=min(mkzh,k+1)
               km1=max(1,k-1)
            do j=1,mjx-1
            do i=1,miy-1
               pl3(i,j,k)=
     &            -.01*scr3b(i,j,k)*(scr3a(i,j,kp1)-scr3a(i,j,km1))/
     &            (ght(i,j,kp1)-ght(i,j,km1))
            enddo
            enddo
            enddo
            engplttl(ipl)='2D frontogenesis (ver. adv.)'
         elseif (cfeld(ifld,ipl)(6:8).eq.'cnf') then
            do k=1,mkzh
            do j=1,mjx
            do i=1,miy
               scr3b(i,j,k)=cosa*uuu(i,j,k)+sina*vvv(i,j,k) !wind along x-sec.
            enddo
            enddo
            enddo
            call smoothcp(scr3b,1,scr3c,prs,pslab1,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            call derivcz(scr3b,0,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3c,1,'x',miy,mjx,mkzh)
            call derivcz(scr3b,0,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3d,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               duadxa=cosa*scr3c(i,j,k)+sina*scr3d(i,j,k)
               pl3(i,j,k)=-scr3a(i,j,k)*duadxa
            enddo
            enddo
            enddo
            engplttl(ipl)='2D frontogenesis (cnf. term)'
         elseif (cfeld(ifld,ipl)(6:8).eq.'shr') then
            do k=1,mkzh
            do j=1,mjx
            do i=1,miy
               scr3b(i,j,k)=sina*uuu(i,j,k)-cosa*vvv(i,j,k) !wind into x-sec.
            enddo
            enddo
            enddo
            call smoothcp(scr3b,0,scr3c,prs,pslab1,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            call derivcz(scr3b,0,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3c,1,'x',miy,mjx,mkzh)
            call derivcz(scr3b,0,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3d,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               duidxi=sina*scr3c(i,j,k)-cosa*scr3d(i,j,k)
               pl3(i,j,k)=-scr3a(i,j,k)*duidxi
            enddo
            enddo
            enddo
            engplttl(ipl)='2D frontogenesis (shr. term)'
         elseif (cfeld(ifld,ipl)(6:8).eq.'til') then
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               scr3b(i,j,k)=www(i,j,k)
            enddo
            enddo
            enddo
            call smoothcp(scr3b,1,scr3c,prs,pslab1,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            call derivcz(scr3b,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3c,1,'x',miy,mjx,mkzh)
            call derivcz(scr3b,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3d,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
               kp1=min(mkzh,k+1)
               km1=max(1,k-1)
            do j=1,mjx-1
            do i=1,miy-1
               pl3(i,j,k)=-.01*(scr3a(i,j,kp1)-
     &            scr3a(i,j,km1))/(ght(i,j,kp1)-ght(i,j,km1))*
     &            (cosa*scr3c(i,j,k)+sina*scr3d(i,j,k))
            enddo
            enddo
            enddo
            engplttl(ipl)='2D frontogenesis (til. term)'
         elseif (cfeld(ifld,ipl)(6:8).eq.'dia') then
            call vervcalc(pstd,dmap,sigf,xmap,pstx,sigh,
     &         uuu,vvv,qvp,tmk,www,prs,scr3a,2,miy,mjx,mkzh)
            call condheat(tmk,qvp,scr3a,0,prs,miy,mjx,mkzh,scr3b)
            call smoothcp(scr3b,1,scr3c,prs,pslab1,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            call derivcz(scr3b,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3c,1,'x',miy,mjx,mkzh)
            call derivcz(scr3b,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3d,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               pl3(i,j,k)=2.778e-4*(cosa*scr3c(i,j,k)+
     &            sina*scr3d(i,j,k)) ! K/m/h to K/m/s
            enddo
            enddo
            enddo
            engplttl(ipl)='2D frontogenesis (dia. term)'
         elseif (cfeld(ifld,ipl)(6:8).eq.'dii') then
            call vervcalc(pstd,dmap,sigf,xmap,pstx,sigh,
     &         uuu,vvv,qvp,tmk,www,prs,scr3a,2,miy,mjx,mkzh)
            call condheat(tmk,qvp,scr3a,1,prs,miy,mjx,mkzh,scr3b)
            call smoothcp(scr3b,1,scr3c,prs,pslab1,iqgsm(ipl),
     &         miy,mjx,mkzh,mabpl,morpl)
            call derivcz(scr3b,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3c,1,'x',miy,mjx,mkzh)
            call derivcz(scr3b,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3d,1,'y',miy,mjx,mkzh)
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               pl3(i,j,k)=2.778e-4*(cosa*scr3c(i,j,k)+
     &            sina*scr3d(i,j,k)) ! K/m/h to K/m/s
            enddo
            enddo
            enddo
            engplttl(ipl)='2D frontogenesis (dia.(ice) term)'
         endif
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=pl3(i,j,k)*3.6e8 ! K/m/s to K/100km/hour
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='K (100 km h)~S~-1~N~'
         if (cfeld(ifld,ipl)(3:3).eq.'e') then
            engplttl(ipl)(4:16)='theta-e frgn.'
         endif
      elseif (cfeld(ifld,ipl)(1:5).eq.'ugeo '.or.
     &        cfeld(ifld,ipl)(1:6).eq.'uageo ') then ! u_geos. or
c                                                      u_ageos., m/s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call derivcz(prs,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,pl3,1,'y',miy,mjx,mkzh)
         fbar=cor(miy/2,mjx/2)
         do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               tv=virtual(tmk(i,j,k),qvp(i,j,k))
               pl3(i,j,k)=-rgas*tv*pl3(i,j,k)/(prs(i,j,k)*fbar)
            enddo
            enddo
            call xtodot(pl3(1,1,k),miy,mjx)
            if (cfeld(ifld,ipl)(1:6).eq.'uageo ') then
               do j=1,mjx
               do i=1,miy
                  pl3(i,j,k)=uuu(i,j,k)-pl3(i,j,k)
               enddo
               enddo
            else
               if (rstmv(2,ipl).ne.0.0) then
                  do j=1,mjx
                  do i=1,miy
                     pl3(i,j,k)=pl3(i,j,k)-rstmv(2,ipl)
                  enddo
                  enddo
               endif
            endif
         enddo
         if (cfeld(ifld,ipl)(1:6).eq.'uageo ') then
            engplttl(ipl)='Ageostrophic wind (x-comp.)'
         else
            engplttl(ipl)='Geostrophic wind (x-comp.)'
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:5).eq.'vgeo '.or.
     &        cfeld(ifld,ipl)(1:6).eq.'vageo ') then ! v_geos. or
c                                                      v_ageos., m/s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call derivcz(prs,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,pl3,1,'x',miy,mjx,mkzh)
         fbar=cor(miy/2,mjx/2)
         do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               tv=virtual(tmk(i,j,k),qvp(i,j,k))
               pl3(i,j,k)=rgas*tv*pl3(i,j,k)/(prs(i,j,k)*fbar)
            enddo
            enddo
            call xtodot(pl3(1,1,k),miy,mjx)
            if (cfeld(ifld,ipl)(1:6).eq.'vageo ') then
               do j=1,mjx
               do i=1,miy
                  pl3(i,j,k)=vvv(i,j,k)-pl3(i,j,k)
               enddo
               enddo
            else
               if (rstmv(1,ipl).ne.0.0) then
                  do j=1,mjx
                  do i=1,miy
                     pl3(i,j,k)=pl3(i,j,k)-rstmv(1,ipl)
                  enddo
                  enddo
               endif
            endif
         enddo
         if (cfeld(ifld,ipl)(1:6).eq.'uageo ') then
            engplttl(ipl)='Ageostrophic wind (y-comp.)'
         else
            engplttl(ipl)='Geostrophic wind (y-comp.)'
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:6).eq.'xptgeo') then ! geo.wind prlel. to x-sec
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call derivcz(prs,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,scr3a,1,'y',miy,mjx,mkzh)
         call derivcz(prs,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
         icen=nint(.5*(rcrag(2,ipl)+rcrbg(2,ipl)))
         jcen=nint(.5*(rcrag(1,ipl)+rcrbg(1,ipl)))
         fbar=cor(icen,jcen)
         cosa=xdist/xseclen
         sina=ydist/xseclen
         do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               tv=virtual(tmk(i,j,k),qvp(i,j,k))
               ugeo=-rgas*tv*scr3a(i,j,k)/(prs(i,j,k)*fbar)-rstmv(2,ipl)
               vgeo= rgas*tv*scr3b(i,j,k)/(prs(i,j,k)*fbar)-rstmv(1,ipl)
               pl3(i,j,k)=cosa*ugeo+sina*vgeo
            enddo
            enddo
            call xtodot(pl3(1,1,k),miy,mjx)
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Geostrophic wind along cross section'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:7).eq.'xptageo') then !ageo. wind prll. to xsec
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call derivcz(prs,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,scr3a,1,'y',miy,mjx,mkzh)
         call derivcz(prs,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
         icen=nint(.5*(rcrag(2,ipl)+rcrbg(2,ipl)))
         jcen=nint(.5*(rcrag(1,ipl)+rcrbg(1,ipl)))
         fbar=cor(icen,jcen)
         cosa=xdist/xseclen
         sina=ydist/xseclen
         do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               tv=virtual(tmk(i,j,k),qvp(i,j,k))
               ugeo=-rgas*tv*scr3a(i,j,k)/(prs(i,j,k)*fbar)
               vgeo= rgas*tv*scr3b(i,j,k)/(prs(i,j,k)*fbar)
               pl3(i,j,k)=cosa*(uuu(i,j,k)-ugeo)+sina*(vvv(i,j,k)-vgeo)
            enddo
            enddo
            call xtodot(pl3(1,1,k),miy,mjx)
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Ageostrophic wind along cross section'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:6).eq.'xntgeo') then ! geo.wind norm. to x-sec
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call derivcz(prs,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,scr3a,1,'y',miy,mjx,mkzh)
         call derivcz(prs,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
         icen=nint(.5*(rcrag(2,ipl)+rcrbg(2,ipl)))
         jcen=nint(.5*(rcrag(1,ipl)+rcrbg(1,ipl)))
         fbar=cor(icen,jcen)
         cosa=xdist/xseclen
         sina=ydist/xseclen
         do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               tv=virtual(tmk(i,j,k),qvp(i,j,k))
               ugeo=-rgas*tv*scr3a(i,j,k)/(prs(i,j,k)*fbar)-rstmv(2,ipl)
               vgeo= rgas*tv*scr3b(i,j,k)/(prs(i,j,k)*fbar)-rstmv(1,ipl)
               pl3(i,j,k)=-sina*ugeo+cosa*vgeo
            enddo
            enddo
            call xtodot(pl3(1,1,k),miy,mjx)
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Geostrophic wind into cross section'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:7).eq.'xntageo') then !ageo. wind norm. to xsec
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call derivcz(prs,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,scr3a,1,'y',miy,mjx,mkzh)
         call derivcz(prs,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
         icen=nint(.5*(rcrag(2,ipl)+rcrbg(2,ipl)))
         jcen=nint(.5*(rcrag(1,ipl)+rcrbg(1,ipl)))
         fbar=cor(icen,jcen)
         cosa=xdist/xseclen
         sina=ydist/xseclen
         do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               tv=virtual(tmk(i,j,k),qvp(i,j,k))
               ugeo=-rgas*tv*scr3a(i,j,k)/(prs(i,j,k)*fbar)
               vgeo= rgas*tv*scr3b(i,j,k)/(prs(i,j,k)*fbar)
               pl3(i,j,k)=-sina*(uuu(i,j,k)-ugeo)+cosa*(vvv(i,j,k)-vgeo)
            enddo
            enddo
            call xtodot(pl3(1,1,k),miy,mjx)
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Ageostrophic wind into cross section'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:3).eq.'slp') then! sea-lev press, hPa
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         if (cfeld(ifld,ipl)(4:7).eq.'old ') then! simple SLP, hPa
            call slpcalc(qvp,tmk,pstx,prs,sigh,ter,pl2,
     &         miy,mjx,mkzh)
         elseif (cfeld(ifld,ipl)(1:6).eq.'slpbm '.or.
     &           cfeld(ifld,ipl)(1:4).eq.'slp ') then! B&M SLP, hPa
c
c         Inspired by Benjamin and Miller (199?).  See subroutine
c         for details.
c
            call slpbmcalc(qvp,tmk,pstx,prs,sigh,sigf,ter,pl2,
     &         miy,mjx,mkzh)
         elseif (cfeld(ifld,ipl)(1:6).eq.'slpgr ') then! GRAPH SLP, hPa
c
c         This is the version in Graph
c
            call sfpcalc(pstx,sigh,prs,scr2a,miy,mjx,mkzh)
            call seaprs(tmk,prs,ter,scr2a,miy,mjx,mkzh,pl2,iup)
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Sea-level pressure'
         unwk(ipl)='hPa'
      elseif (cfeld(ifld,ipl)(1:4).eq.'dpbh') then! pressure diff.
c                                                   betw. htl levels, hPa
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         read(cfeld(ifld,ipl)(5:10),'(2f3.0)') h1,h2
         call dpbhcalc(qvp,tmk,ght,prs,sigf,sigh,pstx,ter,pl2,h1,h2,
     &      miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Press. diff. bet. '//cfeld(ifld,ipl)(5:6)//
     &      '.'//cfeld(ifld,ipl)(7:7)//' and '//cfeld(ifld,ipl)(8:9)//
     &      '.'//cfeld(ifld,ipl)(10:10)//' km'
         unwk(ipl)='hPa'
      elseif (cfeld(ifld,ipl)(1:4).eq.'thck') then! thickness
c                                                   betw. prs levels, dam
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         read(cfeld(ifld,ipl)(5:10),'(2f3.0)') p1,p2
         call thickcalc(qvp,tmk,ght,prs,sigh,sigf,ter,pstx,
     &      pl2,p1,p2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)=cfeld(ifld,ipl)(5:7)//'0'//' to '//
     &      cfeld(ifld,ipl)(8:10)//'0 hPa thickness'
         unwk(ipl)='dam'
      elseif (cfeld(ifld,ipl)(1:4).eq.'ctt ') then! cld-top temp, deg C
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call readdat(iudata,fname,iendf1,'qcw       ',
     &      miy,mjx,mkzh,3,1,scr3a,istat)
         if (iice.eq.1) then
            call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
            call readdat(iudata,fname,iendf1,'qci       ',
     &         miy,mjx,mkzh,3,1,scr3b,istat)
         endif
         call cttcalc(pstx,sigf,sigh,tmk,scr3a,scr3b,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Cloud-top temperature'
         unwk(ipl)='~S~o~N~C'
      elseif (cfeld(ifld,ipl)(1:5).eq.'cap3 '.or.   ! 3D CAPE, J/kg
     &        cfeld(ifld,ipl)(1:5).eq.'cin3 ') then ! 3D conv. inhib., J/kg
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call pfcalc(sigh,sigf,pstx,prs,scr3c,miy,mjx,mkzh)
         call capecalc3d(prs,tmk,qvp,ght,ter,scr3c,scr3a,scr3b,
     &      miy,mjx,mkzh,1)
         if (cfeld(ifld,ipl)(1:4).eq.'cap3') then
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               pl3(i,j,k)=scr3a(i,j,k)
            enddo
            enddo
            enddo
            engplttl(ipl)='CAPE'
         elseif (cfeld(ifld,ipl)(1:4).eq.'cin3') then
            do k=1,mkzh
            do j=1,mjx-1
            do i=1,miy-1
               pl3(i,j,k)=scr3b(i,j,k)
            enddo
            enddo
            enddo
            engplttl(ipl)='Convective inhibition'
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='J kg~S~-1~N~'
c
c      In the following: if "save" was requested, it will be done at the end
c      of this routine, but we might as well also save the other field
c      that was calculated (either cap3 or cin3), because they take so
c      long to calculate. 
c
         if (csave(ipl).ne.'dontsave  '.and..not.lgrad(ipl).and..not.
     &       llapl(ipl).and..not.lhadv(ipl).and.ismcp(ipl).eq.0
     &       .and.ismth(ipl).eq.0) then
            if (cfeld(ifld,ipl)(1:5).eq.'cap3 ') then
               vardesc='Conv. Inh., J/kg'
               plchun=unwk(ipl)
               ndimen=3
               call writefile(scr3b,'cin3      ',
     &            0,ndimen,icdwk(ipl),vardesc,plchun,fname,iendf1,
     &            ihrip,rhrip,chrip,fullsigma,halfsigma,miy,mjx,mkzh)
            elseif (cfeld(ifld,ipl)(1:5).eq.'cin3 ') then
               vardesc='CAPE, J/kg'
               plchun=unwk(ipl)
               ndimen=3
               call writefile(scr3a,'cap3      ',
     &            0,ndimen,icdwk(ipl),vardesc,plchun,fname,iendf1,
     &            ihrip,rhrip,chrip,fullsigma,halfsigma,miy,mjx,mkzh)
            endif
         endif
      elseif (cfeld(ifld,ipl)(1:5).eq.'mcap '.or.   ! Max CAPE, J/kg
     &        cfeld(ifld,ipl)(1:5).eq.'mcin ') then ! assoc. CIN, J/kg
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call pfcalc(sigh,sigf,pstx,prs,scr3c,miy,mjx,mkzh)
         call capecalc3d(prs,tmk,qvp,ght,ter,scr3c,scr3a,scr3b,
     &      miy,mjx,mkzh,0)
         if (cfeld(ifld,ipl)(1:4).eq.'mcap') then
            do j=1,mjx-1
            do i=1,miy-1
               pl2(i,j)=scr3a(i,j,mkzh)
            enddo
            enddo
            engplttl(ipl)='CAPE (for parcel with max theta-e)'
         elseif (cfeld(ifld,ipl)(1:4).eq.'mcin') then
            do j=1,mjx-1
            do i=1,miy-1
               pl2(i,j)=scr3b(i,j,mkzh)
            enddo
            enddo
            engplttl(ipl)='CIN (for parcel with max theta-e)'
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='J kg~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'lcl '.or.
     &        cfeld(ifld,ipl)(1:4).eq.'lfc ') then! LCL or LFC, meters AGL
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
c
c      Put LCL (LFC) in the mkzh-1 (mkzh-2) slab of scr3b
c
         call pfcalc(sigh,sigf,pstx,prs,scr3c,miy,mjx,mkzh)
         call capecalc3d(prs,tmk,qvp,ght,ter,scr3c,scr3a,scr3b,
     &      miy,mjx,mkzh,0)
c
c      Now put it in pl2
c
         if (cfeld(ifld,ipl)(1:4).eq.'lfc ') then
            kget=mkzh-2
            engplttl(ipl)='LFC (for parcel with max theta-e)'
         else
            kget=mkzh-1
            engplttl(ipl)='LCL (for parcel with max theta-e)'
         endif
         do j = 1, mjx-1
         do i = 1, miy-1
            pl2(i,j)=scr3b(i,j,kget)
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='m (AGL)'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:4).eq.'tmcl') then!tmp. of lifted parcel, K
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         read(cfeld(ifld,ipl)(5:6),'(i2)') kpar !nth level from bottom
         kpar=mkzh-kpar+1   ! sigma level
         call liftparcel(prs,tmk,qvp,ght,pl3,kpar,0,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Lifted T from sig=.000'
         write(engplttl(ipl)(19:22),'(f4.3)') sigh(kpar)
         unwk(ipl)='~S~o~N~C'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tvcl') then! T_v of lifted parc., K
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         read(cfeld(ifld,ipl)(5:6),'(i2)') kpar !nth level from bottom
         kpar=mkzh-kpar+1   ! sigma level
         call liftparcel(prs,tmk,qvp,ght,pl3,kpar,1,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Lifted Tv from sig=.000'
         write(engplttl(ipl)(20:23),'(f4.3)') sigh(kpar)
         unwk(ipl)='~S~o~N~C'
      elseif (cfeld(ifld,ipl)(1:5).eq.'ethmx') then
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2c,wk,maxslab)
         call eqthecalc(qvp,tmk,prs,scr3a,miy,mjx,mkzh)
         do j=1,mjx-1
         do i=1,miy-1
            eqthetamax=0.
            do k=1,mkzh
               if (ght(i,j,k)-ter(i,j).lt.3000..and.
     &             scr3a(i,j,k).gt.eqthetamax) then
                  kmax=k
                  eqthetamax=scr3a(i,j,k)
               endif
            enddo
            pmax=prs(i,j,kmax)
            if (cfeld(ifld,ipl)(6:6).eq.'v') then
               pl2(i,j)=eqthetamax
            elseif (cfeld(ifld,ipl)(6:6).eq.'p') then
               pl2(i,j)=pmax
            endif
         enddo
         enddo
         if (cfeld(ifld,ipl)(6:6).eq.'v') then! max eth, K
            engplttl(ipl)='Max theta-e below 3000 m AGL'
            unwk(ipl)='K'
         elseif (cfeld(ifld,ipl)(6:6).eq.'p') then! p-lev. of max eth, hPa
            engplttl(ipl)='Prs. at max theta-e below 3000 m AGL'
            unwk(ipl)='hPa'
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)=' '
         unwk(ipl)='K'
      elseif (cfeld(ifld,ipl)(1:4).eq.'freg') then! flight regul.
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
         if (cfeld(ifld,ipl)(6:6).eq.'b') then! Bocchieri distinction
            call wetbulbcalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
            call precprob(tmk,scr3a,ght,ter,scr3b,miy,mjx,mkzh)
            call addorfill(scr3b(1,1,3),scr2c,miy,mjx,mkzh,2,1,1.,0.)
         endif
         call readdat(iudata,fname,iendf1,'qcw       ',
     &      miy,mjx,mkzh,3,1,scr3a,istat)
         call readdat(iudata,fname,iendf1,'qra       ',
     &      miy,mjx,mkzh,3,1,scr3b,istat)
         if (iice.eq.1) then
            call readdat(iudata,fname,iendf1,'qci       ',
     &         miy,mjx,mkzh,3,1,scr3c,istat)
            call readdat(iudata,fname,iendf1,'qsn       ',
     &         miy,mjx,mkzh,3,1,scr3d,istat)
         endif
         if (cfeld(ifld,ipl)(5:5).eq.'c'.or.
     &       cfeld(ifld,ipl)(5:5).eq.'b') then! get ceiling
            clgfac=1.
            call ceilingcalc(qvp,scr3a,scr3b,scr3c,scr3d,tmk,sigh,
     &         sigf,pstx,prs,scr2c,cfeld(ifld,ipl)(6:6),clgfac,scr2a,
     &         miy,mjx,mkzh)
         else
            call fillarray(scr2a,miy*mjx,rmsg)
         endif
         if (cfeld(ifld,ipl)(5:5).eq.'v'.or.
     &       cfeld(ifld,ipl)(5:5).eq.'b') then! get visibility
            call viscalc(qvp,scr3a,scr3b,scr3c,scr3d,tmk,prs,
     &         scr2c,cfeld(ifld,ipl)(6:6),scr2b,
     &         miy,mjx,mkzh)
         else
            call fillarray(scr2b,miy*mjx,rmsg)
         endif
         call fregcalc(scr2a,scr2b,cfeld(ifld,ipl)(5:5),pl2,miy,mjx)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Flight regulation category'
         unwk(ipl)='none'
      elseif (cfeld(ifld,ipl)(1:5).eq.'frego') then! obs flight regul.
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call readdat(iudata,fname,iendf1,'clgo      ',
     &      miy,mjx,mkzh,2,1,scr2a,istat)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call readdat(iudata,fname,iendf1,'viso      ',
     &      miy,mjx,mkzh,2,1,scr2b,istat)
         do j=1,mjx-1
         do i=1,miy-1
            scr2a(i,j)=.3048*scr2a(i,j)  ! ft. to m
            scr2b(i,j)=1.609*scr2b(i,j)  ! miles to km
         enddo
         enddo
         call fregcalc(scr2a,scr2b,cfeld(ifld,ipl)(6:6),
     &      pl2,miy,mjx)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Obs. flight regulation category'
         unwk(ipl)='none'
      elseif (cfeld(ifld,ipl)(1:5).eq.'fregf') then! FT flight regul.
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call readdat(iudata,fname,iendf1,'clgf      ',
     &      miy,mjx,mkzh,2,1,scr2a,istat)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call readdat(iudata,fname,iendf1,'visf      ',
     &      miy,mjx,mkzh,2,1,scr2b,istat)
         do j=1,mjx-1
         do i=1,miy-1
            scr2a(i,j)=.3048*scr2a(i,j)  ! ft. to m
            scr2b(i,j)=1.609*scr2b(i,j)  ! miles to km
         enddo
         enddo
         call fregcalc(scr2a,scr2b,cfeld(ifld,ipl)(6:6),
     &      pl2,miy,mjx)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Fcst. flight regulation category'
         unwk(ipl)='none'
      elseif (cfeld(ifld,ipl)(1:3).eq.'vis'.and.
     &        cfeld(ifld,ipl)(4:4).ne.'o'.and.
     &        cfeld(ifld,ipl)(4:4).ne.'f') then! hor. vis., km
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
         if (cfeld(ifld,ipl)(4:4).eq.'b') then! Bocchieri distinction
            call wetbulbcalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
            call precprob(tmk,scr3a,ght,ter,scr3b,miy,mjx,mkzh)
            call addorfill(scr3b(1,1,3),scr2a,miy,mjx,mkzh,2,1,1.,0.)
         endif
         call readdat(iudata,fname,iendf1,'qcw       ',
     &      miy,mjx,mkzh,3,1,scr3a,istat)
         call readdat(iudata,fname,iendf1,'qra       ',
     &      miy,mjx,mkzh,3,1,scr3b,istat)
         if (iice.eq.1) then
            call readdat(iudata,fname,iendf1,'qci       ',
     &         miy,mjx,mkzh,3,1,scr3c,istat)
            call readdat(iudata,fname,iendf1,'qsn       ',
     &         miy,mjx,mkzh,3,1,scr3d,istat)
         endif
         call viscalc(qvp,scr3a,scr3b,scr3c,scr3d,tmk,prs,
     &      scr2a,cfeld(ifld,ipl)(4:4),pl2,miy,mjx,mkzh)
         do j=1,mjx-1
         do i=1,miy-1
            if (ter(i,j).le.0.0.and.nint(xlus(i,j)).eq.iwater) then
               pl2(i,j)=rmsg
            endif
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Visibility'
         unwk(ipl)='km'
      elseif (cfeld(ifld,ipl)(1:5).eq.'viso ') then! anal obs vis., km
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call readdat(iudata,fname,iendf1,'viso      ',
     &      miy,mjx,mkzh,2,1,pl2,istat)
         do j=1,mjx-1
         do i=1,miy-1
            if (ter(i,j).le.0.0.and.nint(xlus(i,j)).eq.iwater) then
               pl2(i,j)=rmsg
            else
               pl2(i,j)=1.609*pl2(i,j)  ! miles to km
            endif
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Obs. visibility'
         unwk(ipl)='km'
      elseif (cfeld(ifld,ipl)(1:5).eq.'visf ') then! anal FT vis., km
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call readdat(iudata,fname,iendf1,'visf      ',
     &      miy,mjx,mkzh,2,1,pl2,istat)
         do j=1,mjx-1
         do i=1,miy-1
            if (ter(i,j).le.0.0.and.nint(xlus(i,j)).eq.iwater) then
               pl2(i,j)=rmsg
            else
               pl2(i,j)=1.609*pl2(i,j)  ! miles to km
            endif
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='FT-fcst. visibility'
         unwk(ipl)='km'
      elseif (cfeld(ifld,ipl)(1:3).eq.'clg'.and.
     &        cfeld(ifld,ipl)(4:4).ne.'o'.and.
     &        cfeld(ifld,ipl)(4:4).ne.'f') then! cloud ceiling, m
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
         if (cfeld(ifld,ipl)(4:4).eq.'b') then! Bocchieri distinction
            call wetbulbcalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
            call precprob(tmk,scr3a,ght,ter,scr3b,miy,mjx,mkzh)
            call addorfill(scr3b(1,1,3),scr2a,miy,mjx,mkzh,2,1,1.,0.)
         endif
         call readdat(iudata,fname,iendf1,'qcw       ',
     &      miy,mjx,mkzh,3,1,scr3a,istat)
         call readdat(iudata,fname,iendf1,'qra       ',
     &      miy,mjx,mkzh,3,1,scr3b,istat)
         if (iice.eq.1) then
            call readdat(iudata,fname,iendf1,'qci       ',
     &         miy,mjx,mkzh,3,1,scr3c,istat)
            call readdat(iudata,fname,iendf1,'qsn       ',
     &         miy,mjx,mkzh,3,1,scr3d,istat)
         endif
         call ceilingcalc(qvp,scr3a,scr3b,scr3c,scr3d,tmk,sigh,sigf,
     &      pstx,prs,scr2a,cfeld(ifld,ipl)(4:4),clgfac,pl2,
     &      miy,mjx,mkzh)
         do j=1,mjx-1
         do i=1,miy-1
            if (ter(i,j).le.0.0.and.nint(xlus(i,j)).eq.iwater) then
               pl2(i,j)=rmsg
            endif
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Cloud ceiling'
         unwk(ipl)='km'
      elseif (cfeld(ifld,ipl)(1:5).eq.'clgo ') then! anal obs clg., km
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call readdat(iudata,fname,iendf1,'clgo      ',
     &      miy,mjx,mkzh,2,1,pl2,istat)
         do j=1,mjx-1
         do i=1,miy-1
            if (ter(i,j).le.0.0.and.nint(xlus(i,j)).eq.iwater) then
               pl2(i,j)=rmsg
            else
               pl2(i,j)=.3048*pl2(i,j)  ! feet to meters
            endif
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Obs. cloud ceiling'
         unwk(ipl)='km'
      elseif (cfeld(ifld,ipl)(1:5).eq.'clgf ') then! anal FT clg., km
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call readdat(iudata,fname,iendf1,'clgf      ',
     &      miy,mjx,mkzh,2,1,pl2,istat)
         do j=1,mjx-1
         do i=1,miy-1
            if (ter(i,j).le.0.0.and.nint(xlus(i,j)).eq.iwater) then
               pl2(i,j)=rmsg
            else
               pl2(i,j)=.3048*pl2(i,j)  ! feet to meters
            endif
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='FT-fcst. cloud ceiling'
         unwk(ipl)='km'
      elseif (cfeld(ifld,ipl)(1:3).eq.'ext') then! extinc. coef, per km
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
         if (cfeld(ifld,ipl)(4:4).eq.'b') then! Bocchieri distinction
            call wetbulbcalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
            call precprob(tmk,scr3a,ght,ter,scr3b,miy,mjx,mkzh)
            call addorfill(scr3b(1,1,3),scr2a,miy,mjx,mkzh,2,1,1.,0.)
         endif
         call readdat(iudata,fname,iendf1,'qcw       ',
     &      miy,mjx,mkzh,3,1,scr3a,istat)
         call readdat(iudata,fname,iendf1,'qra       ',
     &      miy,mjx,mkzh,3,1,scr3b,istat)
         if (iice.eq.1) then
            call readdat(iudata,fname,iendf1,'qci       ',
     &         miy,mjx,mkzh,3,1,scr3c,istat)
            call readdat(iudata,fname,iendf1,'qsn       ',
     &         miy,mjx,mkzh,3,1,scr3d,istat)
         endif
         if (cfeld(ifld,ipl)(5:5).eq.'c') then! due to cloud water
            call extingcalc(qvp,scr3a,scr3b,scr3c,scr3d,tmk,
     &         prs,scr2a,cfeld(ifld,ipl)(4:4),1,pl2,
     &         miy,mjx,mkzh)
         elseif (cfeld(ifld,ipl)(5:5).eq.'r') then! due to rain
            call extingcalc(qvp,scr3a,scr3b,scr3c,scr3d,tmk,
     &         prs,scr2a,cfeld(ifld,ipl)(4:4),2,pl2,
     &         miy,mjx,mkzh)
         elseif (cfeld(ifld,ipl)(5:5).eq.'i') then! due to cloud ice
            call extingcalc(qvp,scr3a,scr3b,scr3c,scr3d,tmk,
     &         prs,scr2a,cfeld(ifld,ipl)(4:4),3,pl2,
     &         miy,mjx,mkzh)
         elseif (cfeld(ifld,ipl)(5:5).eq.'s') then! due to snow
            call extingcalc(qvp,scr3a,scr3b,scr3c,scr3d,tmk,
     &         prs,scr2a,cfeld(ifld,ipl)(4:4),4,pl2,
     &         miy,mjx,mkzh)
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Extinction coefficient due to '//
     &      cfeld(ifld,ipl)(5:5)
         unwk(ipl)='km~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:3).eq.'boc') then! Bocchieri prec prob
c
c      Fourth character determines desired prob:
c        'l':liquid; 'f':freezing; 'i':ice
c
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call wetbulbcalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
         call precprob(tmk,scr3a,ght,ter,scr3b,miy,mjx,mkzh)
         if (cfeld(ifld,ipl)(4:4).eq.'l') then
            n=1
            engplttl(ipl)='Prob. of liquid precip.'
         elseif (cfeld(ifld,ipl)(4:4).eq.'f') then
            n=2
            engplttl(ipl)='Prob. of freezing/mixed precip.'
         elseif (cfeld(ifld,ipl)(4:4).eq.'i') then
            n=3
            engplttl(ipl)='Prob. of frozen precip.'
         endif
         do j=1,mjx-1
         do i=1,miy-1
            pl2(i,j)=scr3b(i,j,n)
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='%'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tgc ') then! ground temp., deg.C
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call readdat(iudata,fname,iendf1,'tgk       ',
     &      miy,mjx,mkzh,2,1,pl2,istat)
         do j=1,mjx-1
         do i=1,miy-1
            pl2(i,j)=pl2(i,j)-celkel
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Ground/sea-surface temperature'
         unwk(ipl)='~S~o~N~C'
      elseif (cfeld(ifld,ipl)(1:4).eq.'pvo ') then! pot. vorticity, PVU
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call thecalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
         call pvocalc(sigh,xmap,uuu,vvv,cor,scr3a,prs,
     &      pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Potential vorticity'
         unwk(ipl)='PVU'
      elseif (cfeld(ifld,ipl)(1:4).eq.'pvm ') then! moist PV, PVU
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call eqthecalc(qvp,tmk,prs,scr3a,miy,mjx,mkzh)
         call pvocalc(sigh,xmap,uuu,vvv,cor,scr3a,prs,
     &      pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Moist potential vorticity'
         unwk(ipl)='PVU'
      elseif (cfeld(ifld,ipl)(1:4).eq.'vor '.or.! relative vort., per s
     &    cfeld(ifld,ipl)(1:4).eq.'avo ') then! absolute vort., per s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call derivcz(vvv,0,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,pl3,1,'x',miy,mjx,mkzh)
         call derivcz(uuu,0,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,scr3a,1,'y',miy,mjx,mkzh)
         if (cfeld(ifld,ipl)(1:1).eq.'v') then
            icor=0
            engplttl(ipl)='Relative vorticity'
         elseif (cfeld(ifld,ipl)(1:1).eq.'a') then
            icor=1
            engplttl(ipl)='Absolute vorticity'
         endif
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=pl3(i,j,k)-scr3a(i,j,k)+icor*cor(i,j)
c
c         Scale the vorticity by 1.e5
c
            pl3(i,j,k)=  pl3(i,j,k) * 1.e5
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='10~S~-5~N~ s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'div ') then! divergence, per s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call derivcz(uuu,0,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,pl3,1,'x',miy,mjx,mkzh)
         call derivcz(vvv,0,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,scr3a,1,'y',miy,mjx,mkzh)
         call addorfill(scr3a,pl3,miy,mjx,mkzh,3,1,1.,1.)
         call addorfill(pl3,pl3,miy,mjx,mkzh,3,1,1.e5,0.)  ! scale by 1.e5
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Divergence'
         unwk(ipl)='10~S~-5~N~ s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'rhu ') then! RH wrt liq, percent
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call rhucalc(qvp,tmk,prs,0,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Relative humidity (w.r.t. water)'
         unwk(ipl)='%'
      elseif (cfeld(ifld,ipl)(1:4).eq.'rhi ') then! RH wrt ice, percent
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call rhucalc(qvp,tmk,prs,1,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Relative humidity (w.r.t. ice)'
         unwk(ipl)='%'
      elseif (cfeld(ifld,ipl)(1:4).eq.'www '.or.
     &        cfeld(ifld,ipl)(1:4).eq.'omg '.or.
     &        cfeld(ifld,ipl)(1:4).eq.'sgd ') then
c
c      vert velocity (w in cm/s, omega in dPa/s, or d(sigma)/dt
c         in s**-1)
c
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         if (cfeld(ifld,ipl)(1:4).eq.'www ') then
            ind=1
            engplttl(ipl)='Vertical velocity'
            unwk(ipl)='cm s~S~-1~N~'
         elseif (cfeld(ifld,ipl)(1:4).eq.'omg ') then
            ind=2
            engplttl(ipl)='Omega'
            unwk(ipl)='dPa s~S~-1~N~'
         elseif (cfeld(ifld,ipl)(1:4).eq.'sgd ') then
            ind=3
            engplttl(ipl)='Sigma-dot'
            unwk(ipl)='s~S~-1~N~'
         endif
         call vervcalc(pstd,dmap,sigf,xmap,pstx,sigh,
     &      uuu,vvv,qvp,tmk,www,prs,pl3,ind,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
      elseif (cfeld(ifld,ipl)(1:7).eq.'netasc ') then
c
c      Net ascent of trajectories (hPa).  Must have a trajectory
c      position file from a properly organized trajectory run.
c
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call netasc(scr3a,scr3b,sigf,sigh,casename,iendc,
     &      ctjfl(ipl),rtjst(ipl),rtjen(ipl),pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Net ascent of trajectories'
         unwk(ipl)='hPa'
      elseif (cfeld(ifld,ipl)(1:4).eq.'qvx ') then! Q_vector_x,
c                                                    10**-6 m/s**3/hPa
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call vervcalc(pstd,dmap,sigf,xmap,pstx,sigh,
     &      uuu,vvv,qvp,tmk,www,prs,pl3,2,miy,mjx,mkzh)
         call bvfricalc(sigh,pstx,prs,tmk,scr3b,qvp,scr3c,
     &      uuu,vvv,0.,0.,0.,'bvfsqd',0,scr3a,
     &      miy,mjx,mkzh)
         numpas=iqgsm(ipl)
         call qgomg(prs,pl3,tmk,qvp,ght,scr3a,scr3b,scr3c,cor,xmap,ter,
     &      numpas,1,1,0,0,0,0,100.,
     &      fname,iendf1,ihrip,rhrip,fullsigma,halfsigma,chrip,
     &      vardesc,plchun,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Q (x-comp.)'
         unwk(ipl)='m s~S~-3~N~ hPa~S~-1~N~' !not enough room for 10**-6
      elseif (cfeld(ifld,ipl)(1:4).eq.'qvy ') then! Q_vector_y,
c                                                    10**-6 m/s**3/hPa
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call vervcalc(pstd,dmap,sigf,xmap,pstx,sigh,
     &      uuu,vvv,qvp,tmk,www,prs,pl3,2,miy,mjx,mkzh)
         call bvfricalc(sigh,pstx,prs,tmk,scr3b,qvp,scr3c,
     &      uuu,vvv,0.,0.,0.,'bvfsqd',0,scr3a,
     &      miy,mjx,mkzh)
         numpas=iqgsm(ipl)
         call qgomg(prs,pl3,tmk,qvp,ght,scr3a,scr3b,scr3c,cor,xmap,ter,
     &      numpas,2,1,0,0,0,0,100.,
     &      fname,iendf1,ihrip,rhrip,fullsigma,halfsigma,chrip,
     &      vardesc,plchun,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Q (y-comp.)'
         unwk(ipl)='m s~S~-3~N~ hPa~S~-1~N~' !not enough room for 10**-6
      elseif (cfeld(ifld,ipl)(1:6).eq.'qvdiv ') then !div-Q,
c                                                  10**-12(s**3 hPa)**-1
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call vervcalc(pstd,dmap,sigf,xmap,pstx,sigh,
     &      uuu,vvv,qvp,tmk,www,prs,pl3,2,miy,mjx,mkzh)
         call bvfricalc(sigh,pstx,prs,tmk,scr3b,qvp,scr3c,
     &      uuu,vvv,0.,0.,0.,'bvfsqd',0,scr3a,
     &      miy,mjx,mkzh)
         numpas=iqgsm(ipl)
         ivar2=igdir(ipl)
         call qgomg(prs,pl3,tmk,qvp,ght,scr3a,scr3b,scr3c,cor,xmap,ter,
     &      numpas,3,ivar2,0,0,0,0,100.,
     &      fname,iendf1,ihrip,rhrip,fullsigma,halfsigma,
     &      chrip,vardesc,plchun,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Q-vector divergence'
         unwk(ipl)='s~S~-3~N~ hPa~S~-1~N~' !not enough room for 10**-12
      elseif (cfeld(ifld,ipl)(1:5).eq.'qgomf') then !full omega, ubar/s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call vervcalc(pstd,dmap,sigf,xmap,pstx,sigh,
     &      uuu,vvv,qvp,tmk,www,prs,pl3,2,miy,mjx,mkzh)
         call bvfricalc(sigh,pstx,prs,tmk,scr3b,qvp,scr3c,
     &      uuu,vvv,0.,0.,0.,'bvfsqd',0,scr3a,
     &      miy,mjx,mkzh)
         numpas=iqgsm(ipl)
         ivar2=igdir(ipl)
         call qgomg(prs,pl3,tmk,qvp,ght,scr3a,scr3b,scr3c,cor,xmap,ter,
     &      numpas,0,ivar2,0,0,0,0,100.,
     &      fname,iendf1,ihrip,rhrip,fullsigma,halfsigma,
     &      chrip,vardesc,plchun,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Full omega'
         unwk(ipl)='dPa s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:5).eq.'qgomg'.or.
     &        cfeld(ifld,ipl)(1:5).eq.'qmomg') then! QG dp/dt, ubar/s
         read(cfeld(ifld,ipl)(6:6),'(i1)') iqvecforc
         read(cfeld(ifld,ipl)(7:7),'(i1)') itopobc
         read(cfeld(ifld,ipl)(8:8),'(i1)') iekmnbc
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call bvfricalc(sigh,pstx,prs,tmk,scr3b,qvp,scr3c,
     &      uuu,vvv,0.,0.,0.,'bvfsqd',0,scr3a,
     &      miy,mjx,mkzh)
         imo=0
         rhithresh=100.
         engplttl(ipl)='QG omega (dry)'
         if (cfeld(ifld,ipl)(2:2).eq.'m') then
            engplttl(ipl)='QG omega (moist)'
            imo=1
            call readdat(iudata,fname,iendf1,'qcw       ',
     &         miy,mjx,mkzh,3,1,pl3,istat)
            call readdat(iudata,fname,iendf1,'qra       ',
     &         miy,mjx,mkzh,3,1,scr3b,istat)
            call addorfill(scr3b,pl3,miy,mjx,mkzh,3,1,.001,1.)
            if (iice.eq.1) then
               call readdat(iudata,fname,iendf1,'qci       ',
     &            miy,mjx,mkzh,3,1,scr3b,istat)
               call addorfill(scr3b,pl3,miy,mjx,mkzh,3,1,.001,1.)
               call readdat(iudata,fname,iendf1,'qsn       ',
     &            miy,mjx,mkzh,3,1,scr3b,istat)
               call addorfill(scr3b,pl3,miy,mjx,mkzh,3,1,.001,1.)
               call readdat(iudata,fname,iendf1,'qgr       ',
     &            miy,mjx,mkzh,3,0,scr3b,istat)
               if (istat.ne.-1) call addorfill(scr3b,pl3,
     &            miy,mjx,mkzh,3,1,.001,1.)
            endif
            call bvfricalc(sigh,pstx,prs,tmk,scr3c,qvp,pl3,
     &         uuu,vvv,0.,0.,0.,'bvfsqi',1,scr3b,miy,mjx,mkzh)
            read(cfeld(ifld,ipl)(9:10),'(i2)') irhithresh
            rhithresh=float(irhithresh)
            call rhucalc(qvp,tmk,prs,1,scr3c,miy,mjx,mkzh)
         endif
         call vervcalc(pstd,dmap,sigf,xmap,pstx,sigh,
     &      uuu,vvv,qvp,tmk,www,prs,pl3,2,miy,mjx,mkzh)
         numpas=iqgsm(ipl)
         call qgomg(prs,pl3,tmk,qvp,ght,scr3a,scr3b,scr3c,cor,xmap,ter,
     &      numpas,4,1,iqvecforc,itopobc,iekmnbc,imo,rhithresh,
     &      fname,iendf1,ihrip,rhrip,fullsigma,halfsigma,
     &      chrip,vardesc,plchun,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='dPa s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:2).eq.'se'.or.
     &        cfeld(ifld,ipl)(1:2).eq.'sm') then ! Sawyer-Eliassen
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3d,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3e,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3f,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3g,wk,maxslab)
c
c   Put v-geos. in 3a, v-ageos. in 3b, u on cross points in 3c,
c   omega in 3d, and d(ter)/dx in 2a.
c   Note, y axis is assumed to be toward left of cross sec.,
c   x axis is assumed to be into page.
c
         call derivcz(prs,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,scr3a,1,'y',miy,mjx,mkzh)
         call derivcz(prs,1,ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,scr3b,1,'x',miy,mjx,mkzh)
         icen=nint(.5*(rcrag(2,ipl)+rcrbg(2,ipl)))
         jcen=nint(.5*(rcrag(1,ipl)+rcrbg(1,ipl)))
         fbar=cor(icen,jcen)
         cosa=xdist/xseclen
         sina=ydist/xseclen
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            tv=virtual(tmk(i,j,k),qvp(i,j,k))
            ugeo=-rgas*tv*scr3a(i,j,k)/(prs(i,j,k)*fbar)
            vgeo= rgas*tv*scr3b(i,j,k)/(prs(i,j,k)*fbar)
            scr3a(i,j,k)=-cosa*ugeo-sina*vgeo
            utot=.25*(uuu(i,j,k)+uuu(i+1,j,k)+uuu(i,j+1,k)+
     &         uuu(i+1,j+1,k))
            vtot=.25*(vvv(i,j,k)+vvv(i+1,j,k)+vvv(i,j+1,k)+
     &         vvv(i+1,j+1,k))
            scr3b(i,j,k)=-cosa*utot-sina*vtot-scr3a(i,j,k)
            scr3c(i,j,k)=-sina*utot+cosa*vtot
         enddo
         enddo
         enddo
         do j=1,mjx-1
            jp1=min(j+1,mjx-1)
            jm1=max(j-1,1)
         do i=1,miy-1
            ip1=min(i+1,miy-1)
            im1=max(i-1,1)
            dss=ds/xmap(i,j)
            dy=dss*(ip1-im1)
            dx=dss*(jp1-jm1)
            dterdx=(ter(i,jp1)-ter(i,jm1))/dx
            dterdy=(ter(ip1,j)-ter(im1,j))/dy
            scr2a(i,j)=-sina*dterdx+cosa*dterdy
         enddo
         enddo
         call vervcalc(pstd,dmap,sigf,xmap,pstx,sigh,
     &      uuu,vvv,qvp,tmk,www,prs,scr3d,2,miy,mjx,mkzh)
c
c      Get stability fields
c
         call bvfricalc(sigh,pstx,prs,tmk,scr3f,qvp,scr3g,
     &      uuu,vvv,0.,0.,0.,'bvfsqd    ',0,scr3e,
     &      miy,mjx,mkzh)
c
c   Set Saw.-El. parameters
c
         imaxlit=nint(rsepa(1,ipl))
         errmin=rsepa(2,ipl)
         alphor=rsepa(3,ipl)
         ixaverage=nint(rsepa(4,ipl))
         smfac=rsepa(5,ipl)
         imaxbig=nint(rsepa(6,ipl))
         rhithresh=rsepa(7,ipl)
         smfstb=rsepa(8,ipl)
         ilhs=nint(rsepa(9,ipl))
         if (cfeld(ifld,ipl)(2:2).eq.'m') then
            imo=1
            call readdat(iudata,fname,iendf1,'qcw       ',
     &         miy,mjx,mkzh,3,1,pl3,istat)
            call readdat(iudata,fname,iendf1,'qra       ',
     &         miy,mjx,mkzh,3,1,scr3f,istat)
            call addorfill(scr3f,pl3,miy,mjx,mkzh,3,1,.001,1.)
            if (iice.eq.1) then
               call readdat(iudata,fname,iendf1,'qci       ',
     &            miy,mjx,mkzh,3,1,scr3f,istat)
               call addorfill(scr3f,pl3,miy,mjx,mkzh,3,1,.001,1.)
               call readdat(iudata,fname,iendf1,'qsn       ',
     &            miy,mjx,mkzh,3,1,scr3f,istat)
               call addorfill(scr3f,pl3,miy,mjx,mkzh,3,1,.001,1.)
               call readdat(iudata,fname,iendf1,'qgr       ',
     &            miy,mjx,mkzh,3,0,scr3f,istat)
               if (istat.ne.-1) call addorfill(scr3f,pl3,
     &            miy,mjx,mkzh,3,1,.001,1.)
            endif
            call bvfricalc(sigh,pstx,prs,tmk,scr3g,qvp,pl3,
     &         uuu,vvv,0.,0.,0.,'bvfsqi    ',1,scr3f,
     &         miy,mjx,mkzh)
            call rhucalc(qvp,tmk,prs,1,scr3g,miy,mjx,mkzh)
         else
            imo=0
            imaxbig=1
         endif
c
c   Calculate number of pressure levels.  We want grid
c   aspect ratio to be similar to typical frontal slope,
c   i.e. dz/dx ~ .01
c
         pbotse=-9e9
         ptopse=-9e9
         do j=1,mjx-1
         do i=1,miy-1
            pbotse=max(pbotse,prs(i,j,mkzh))
            ptopse=max(ptopse,prs(i,j,1))
         enddo
         enddo
         ptopse=ptopse+1.
         dpapprox=ds*xseclen/(nscrs-1.)*.01/13.
         mkp=nint((pbotse-ptopse)/dpapprox)+1
c
         call fillarray(pl3,miy*mjx*mkzh,0.)
         call saweli(ght,tmk,qvp,prs,scr3a,scr3b,scr3c,scr3d,
     &      scr3e,scr3f,scr3g,imo,ilhs,rhithresh,imaxlit,imaxbig,
     &      errmin,alphor,ter,scr2a,cor,xmap,
     &      pl3,rcrag(1,ipl),rcrbg(1,ipl),ixaverage,smfac,smfstb,
     &      nscrs,nscrs,xdist,ydist,xseclen,ptopse,pbotse,
     &      cfeld(ifld,ipl),mkp,miy,mjx,mkzh)
c
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='unknown'
         if (cfeld(ifld,ipl)(3:5).eq.'vvv'.or.
     &       cfeld(ifld,ipl)(3:5).eq.'vge'.or.
     &       cfeld(ifld,ipl)(3:5).eq.'vag'.or.
     &       cfeld(ifld,ipl)(3:5).eq.'vab'.or.
     &       cfeld(ifld,ipl)(3:5).eq.'vtb') then
            unwk(ipl)='m s~S~-1~N~'
         elseif (cfeld(ifld,ipl)(3:5).eq.'omf'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'omb') then
            unwk(ipl)='dPa s~S~-1~N~'
         elseif (cfeld(ifld,ipl)(3:5).eq.'pvo') then
            unwk(ipl)='PVU'
         elseif (cfeld(ifld,ipl)(3:5).eq.'con'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'she'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'tot') then
            unwk(ipl)='m s~S~-2~N~ hPa~S~-1~N~'
         elseif (cfeld(ifld,ipl)(3:5).eq.'std'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'stm'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'ste'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'saj') then
            unwk(ipl)='m~S~2~N~/(s hPa)~S~2~N~'
         elseif (cfeld(ifld,ipl)(3:5).eq.'rhi') then
            unwk(ipl)='%'
         elseif (cfeld(ifld,ipl)(3:5).eq.'ghp') then
            unwk(ipl)='m'
         elseif (cfeld(ifld,ipl)(3:5).eq.'alp') then
            unwk(ipl)='m~S~2~N~/(hPa s~S~2~N~)'
         elseif (cfeld(ifld,ipl)(3:5).eq.'rsw'.or.
     &           (cfeld(ifld,ipl)(3:4).eq.'fr'.and.
     &            cfeld(ifld,ipl)(5:5).ne.'6')) then
            unwk(ipl)='none'
         elseif (cfeld(ifld,ipl)(3:5).eq.'vor'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'vaj') then
            unwk(ipl)='s~S~-1~N~'
         elseif (cfeld(ifld,ipl)(3:5).eq.'bcl'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'baj') then
            unwk(ipl)='m s~S~-1~N~ hPa~S~-1~N~'
         elseif (cfeld(ifld,ipl)(3:5).eq.'fst'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'faj') then
            unwk(ipl)='m s~S~-1~N~ hPa~S~-2~N~'
         elseif (cfeld(ifld,ipl)(3:5).eq.'fr6'.or.
     &           cfeld(ifld,ipl)(3:5).eq.'psi') then
            unwk(ipl)='m hPa s~S~-1~N~'
         endif
      elseif (cfeld(ifld,ipl)(1:5).eq.'bvfsq'.or.
     &        cfeld(ifld,ipl)(1:5).eq.'richn'.or.
     &        cfeld(ifld,ipl)(1:4).eq.'spsq') then
c
c      Brunt-Vaisala frequency squared (per sec squared), Richardson
c      number (dim'less), or Scorer parameter squared (per km squared)
c
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call fillarray(scr3b,miy*mjx*mkzh,0.)
         call readdat(iudata,fname,iendf1,'qcw       ',
     &      miy,mjx,mkzh,3,0,scr3b,istat)
         call readdat(iudata,fname,iendf1,'qra       ',
     &      miy,mjx,mkzh,3,0,pl3,istat)
         if (istat.ne.-1) call addorfill(pl3,scr3b,miy,mjx,mkzh,
     &      3,1,.001,1.)
         if (iice.eq.1) then
            call readdat(iudata,fname,iendf1,'qci       ',
     &         miy,mjx,mkzh,3,0,pl3,istat)
            if (istat.ne.-1) call addorfill(pl3,scr3b,miy,mjx,mkzh,
     &         3,1,.001,1.)
            call readdat(iudata,fname,iendf1,'qsn       ',
     &         miy,mjx,mkzh,3,0,pl3,istat)
            if (istat.ne.-1) call addorfill(pl3,scr3b,miy,mjx,mkzh,
     &         3,1,.001,1.)
            call readdat(iudata,fname,iendf1,'qgr       ',
     &         miy,mjx,mkzh,3,0,pl3,istat)
            if (istat.ne.-1) call addorfill(pl3,scr3b,miy,mjx,mkzh,
     &         3,1,.001,1.)
         endif
         cphase=sqrt(rstmv(2,ipl)*rstmv(2,ipl)+
     &      rstmv(1,ipl)*rstmv(1,ipl))
         if (cphase.gt.0.) then
            cosa=rstmv(2,ipl)/cphase
            sina=rstmv(1,ipl)/cphase
         else
            cosa=1.
            sina=0.
         endif
         call bvfricalc(sigh,pstx,prs,tmk,scr3a,qvp,scr3b,
     &      uuu,vvv,cosa,sina,cphase,cfeld(ifld,ipl),0,pl3,
     &      miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         if (cfeld(ifld,ipl)(1:5).eq.'bvfsq') then
            engplttl(ipl)='Brunt-Vaisala freq. (squared)'
            unwk(ipl)='s~S~-2~N~'
         elseif (cfeld(ifld,ipl)(1:5).eq.'richn') then
            engplttl(ipl)='Richardson number'
            unwk(ipl)='none'
         elseif (cfeld(ifld,ipl)(1:4).eq.'spsq') then
            engplttl(ipl)='Scorer parameter (squared)'
            unwk(ipl)='km~S~-2~N~'
         endif
      elseif (cfeld(ifld,ipl)(1:4).eq.'rib ') then
c
c      Near-surface Richardson number (dim'less), used in HIRPBL scheme
c
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call readdat(iudata,fname,iendf1,'tgk       ',
     &      miy,mjx,mkzh,2,1,scr2a,istat)
         do j=1,mjx-1
         do i=1,miy-1
            za=ght(i,j,mkzh)-ter(i,j)
            gammam=gamma*(1.+gammamd*qvp(i,j,mkzh))
            pfaca=(1000./prs(i,j,mkzh))**gammam
            psfc=prs(i,j,mkzh)+(1.-sigh(mkzh))*pstx(i,j)
            pfacg=(1000./psfc)**gammam
            tha=tmk(i,j,mkzh)*pfaca
            thva=virtual(tmk(i,j,mkzh),qvp(i,j,mkzh))*pfaca
            thvg=virtual(scr2a(i,j),qvp(i,j,mkzh))*pfacg
            ucross=.25*(uuu(i,j,mkzh)+uuu(i+1,j,mkzh)+
     &                  uuu(i,j+1,mkzh)+uuu(i+1,j+1,mkzh))
            vcross=.25*(vvv(i,j,mkzh)+vvv(i+1,j,mkzh)+
     &                  vvv(i,j+1,mkzh)+vvv(i+1,j+1,mkzh))
            vsq=ucross*ucross+vcross*vcross
            pl2(i,j)=grav*za*(thva-thvg)/(tha*vsq)
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         unwk(ipl)='s~S~-2~N~'
      elseif (cfeld(ifld,ipl)(1:8).eq.'condheat') then !cond. htg., K/day
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call vervcalc(pstd,dmap,sigf,xmap,pstx,sigh,
     &      uuu,vvv,qvp,tmk,www,prs,scr3a,2,miy,mjx,mkzh)
         call condheat(tmk,qvp,scr3a,0,prs,miy,mjx,mkzh,pl3)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Condensational heating'
         unwk(ipl)='K h~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:9).eq.'condheati') then !cond. htg., K/day
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call vervcalc(pstd,dmap,sigf,xmap,pstx,sigh,
     &      uuu,vvv,qvp,tmk,www,prs,scr3a,2,miy,mjx,mkzh)
         call condheat(tmk,qvp,scr3a,1,prs,miy,mjx,mkzh,pl3)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Condensational/fusional heating'
         unwk(ipl)='K h~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:9).eq.'vadvtheta') then !v.adv. of theta, K/hr
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
         call thecalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
         call smoothcp(scr3a,1,scr3c,prs,pslab1,iqgsm(ipl),
     &      miy,mjx,mkzh,mabpl,morpl)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            scr3b(i,j,k)=www(i,j,k)
         enddo
         enddo
         enddo
         call smoothcp(scr3b,1,scr3c,prs,pslab1,iqgsm(ipl),
     &      miy,mjx,mkzh,mabpl,morpl)
         do k=1,mkzh
            kp1=min(mkzh,k+1)
            km1=max(1,k-1)
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=-.01*scr3b(i,j,k)*
     &         (scr3a(i,j,kp1)-scr3a(i,j,km1))/
     &         (ght(i,j,kp1)-ght(i,j,km1))*3600.  ! K/s to K/hr
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Vert. adv. of potential temperature'
         unwk(ipl)='K h~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'wsp ' .or. 
     &        cfeld(ifld,ipl)(1:5).eq.'wspk ') then! Hor wind speed, m/s or kt
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call wspcalc(uuu,vvv,0.,0.,rstmv(1,ipl),pl3,
     &      miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Horizontal wind speed'
         if (cfeld(ifld,ipl)(1:4).eq.'wsp ') then
           unwk(ipl)='m s~S~-1~N~'
         else
           unwk(ipl)='kt'
           call addorfill(pl3,pl3,miy,mjx,mkzh,3,0,1.94,0.)
         endif
      elseif (cfeld(ifld,ipl)(1:4).eq.'wdr ') then! Hor wind dir, deg.
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call wdircalc(uuu,vvv,unorth,vnorth,rstmv(1,ipl),
     &      pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Horizontal wind direction'
         unwk(ipl)='degrees'
      elseif (cfeld(ifld,ipl)(1:4).eq.'xnt ') then! X-nrm tot wind, m/s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         cosa=xdist/xseclen
         sina=ydist/xseclen
         call wspcalc(uuu,vvv,-sina,cosa,rstmv(1,ipl),pl3,
     &      miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Horizontal wind into cross section'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'amt ') then! Abs. mom., m/s
c      Here amt is defined as total wind velocity into the page,
c      plus f (at middle of x-sec.) times left-to-right distance
c      along cross-section.  If the x-axis is taken to be parallel
c      to the cross section and pointing toward the right, then this
c      definition of abs. mom. corresponds to v+fx.
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         cosa=xdist/xseclen
         sina=ydist/xseclen
         call wspcalc(uuu,vvv,-sina,cosa,rstmv(1,ipl),pl3,
     &      miy,mjx,mkzh)
         caxgn=1.+(rcrag(2,ipl)-xjcorn)*refrat
         caygn=1.+(rcrag(1,ipl)-yicorn)*refrat
         cbxgn=1.+(rcrbg(2,ipl)-xjcorn)*refrat
         cbygn=1.+(rcrbg(1,ipl)-yicorn)*refrat
         fzero=cor(nint(.5*(caygn+cbygn)),nint(.5*(caxgn+cbxgn)))
         cosa=xdist/xseclen
         sina=ydist/xseclen
         do j=1,mjx
         do i=1,miy
            dltr=ds*((j-caxgn)*cosa+(i-caygn)*sina)
            do k=1,mkzh
               pl3(i,j,k)=pl3(i,j,k)+fzero*dltr
            enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Absolute momentum'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'xpt ') then! X-paral. wind, m/s
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         cosa=xdist/xseclen
         sina=ydist/xseclen
         call wspcalc(uuu,vvv,cosa,sina,rstmv(1,ipl),pl3,
     &      miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Horizontal wind along cross section'
         unwk(ipl)='m s~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:4).eq.'stb ') then! -d(theta)/dp, K/hPa
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call thecalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
         do k=1,mkzh
         do j=1,mjx
         do i=1,miy
            scr3a(i,j,k)=-scr3a(i,j,k)
         enddo
         enddo
         enddo
         call ddpcalc(prs,scr3a,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='-d(theta)/dp '
         unwk(ipl)='K hPa~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:5).eq.'stbz ') then! d(theta)/dz, K/km
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call thecalc(prs,tmk,qvp,scr3a,miy,mjx,mkzh)
         do k=1,mkzh
            kp1=min(k+1,mkzh)
            km1=max(k-1,1)
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=1000.*(scr3a(i,j,kp1)-scr3a(i,j,km1))/
     &            (ght(i,j,kp1)-ght(i,j,km1))
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='d(theta)/dz'
         unwk(ipl)='K km~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:7).eq.'dthtedz') then! d(theta_e)/dz, K/km
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call eqthecalc(qvp,tmk,prs,scr3a,miy,mjx,mkzh)
         do k=1,mkzh
            kp1=min(k+1,mkzh)
            km1=max(k-1,1)
         do j=1,mjx-1
         do i=1,miy-1
            pl3(i,j,k)=1000.*(scr3a(i,j,kp1)-scr3a(i,j,km1))/
     &            (ght(i,j,kp1)-ght(i,j,km1))
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='d(theta-e)/dz'
         unwk(ipl)='K km~S~-1~N~'
      elseif (cfeld(ifld,ipl)(1:7).eq.'brnshr ') then
c        Bulk Richardson Number Shear (a la Stensrud et al. 1997), m**2/s**2
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call pfcalc(sigh,sigf,pstx,prs,scr3a,miy,mjx,mkzh)
         call brnshr(uuu,vvv,ght,scr3a,ter,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         unwk(ipl)='m~S~2~N~ s~S~-2~N~'

c
c	The following are new fields added by Dave Ahijevych Dec 2001.
c
c	beta thresholding added Feb 2002.
c

      elseif (cfeld(ifld,ipl)(1:3).eq.'ush') then!  shear, m/s
c	ushxMMMNNN: u-component of vertical wind shear, m/s. (2D)
c	 
c
c	MMM and NNN are the lower and upper levels, respectively,
c	expressed with leading zeros included,
c	so that the total number of digits is always 6.
c	Vertical coordinate system is indicated by the character
c	in the 'x' position.
c
c	With height coordinates (both MSL and AGL)
c	MMM and NNN are in hectometers (hm).
c	With pressure coordinates, MMM and NNN are in
c	kPa.	
c	In sigma coordinates, MMM and NNN are merely
c	indices.
c	Right now, only sigma level indices may be used
c	in ushsMMMNNN and usfsMMMNNN.
c	In the future, I might add 'from the bottom'
c	notation 'b' and allow the user to specify the 
c	second from the bottom sigma level as 'b02'.
c
c

         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
c
c	subroutine pfcalc (text lifted from pfcalc.f)
c
c	calculates the pressure at full sigma levels. It is
c	assumed that pressure at the model bottom [i.e., the
c	last full sigma level] is pstx + ptop + p-prime(mkzh),
c	or pr(mkzh) + (1-sigh(mkzh))*pstx, and that the
c	pressure at the model top [i.e., the 1st full sigma
c	level] is ptop + p_prime(1), or pr(1)-sigh(1)*pstx.
c	The array only contains mkzh levels, so the pressure at
c	the model top is excluded.  The kth value of dummy argument
c	pf (or actual argument scr3a) is the
c	pressure at the (k+1)th full sigma level, which is the
c	lower bounding level for the kth sigma layer.  Thus, pf(mkzh),
c	or index mkzh of actual argument scr3a,
c	would be surface pressure.
c	scr3a is in mb or hPa.
c

         call pfcalc(sigh,sigf,pstx,prs,scr3a,miy,mjx,mkzh)


c
c	Read the bottom and top of the layer for which
c	to calculate the wind shear.
c

         read(cfeld(ifld,ipl)(5:10),'(2f3.0)') vlev_bot,vlev_top

c
c	Calculate the vertical wind shear between 
c	vlev_bot and vlev_top,
c	storing it in 2-dimensional variable pl2.
c	Remember to convert vlev_bot, vlev_top
c	to the units that wshear expects.
c
c	
c
c	Set the value for the string holding the 
c	vertical coordinate units, cvcorunits,
c	and create the string holding the plot title,
c	engplttl(ipl).
c
c	Use '-comp.' in the English plot title string,
c	engplttl(ipl).	Note the period.
c	This character pattern, which indicates
c	a "component" field, triggers 
c	pltitle.f to change the plot labels. 
c	This is desired when plotting dual-component
c	fields such as wind barbs or vectors. 
c

         IF (cfeld(ifld,ipl)(4:4) .eq. 'a') THEN
            engplttl(ipl)=cfeld(ifld,ipl)(5:7)//''//' to '//
     &      cfeld(ifld,ipl)(8:10)//' hm AGL wind shr (u-comp.)'
            cvcorunits = 'm AGL'
c		hm to m
            vlev_bot = vlev_bot * 100.
            vlev_top = vlev_top * 100.
         ELSE IF (cfeld(ifld,ipl)(4:4).eq.'z') THEN
            engplttl(ipl)=cfeld(ifld,ipl)(5:7)//''//' to '//
     &      cfeld(ifld,ipl)(8:10)//' hm MSL wind shr (u-comp.)'
            cvcorunits = 'm'
c		hm to m
            vlev_bot = vlev_bot * 100.
            vlev_top = vlev_top * 100.
         ELSE IF (cfeld(ifld,ipl)(4:4).eq.'p') THEN
            engplttl(ipl)=cfeld(ifld,ipl)(5:7)//'0'//' to '//
     &      cfeld(ifld,ipl)(8:10)//'0 hPa wind shear (u-comp.)'
            cvcorunits = 'hPa'
c		kPa to hPa
            vlev_bot = vlev_bot * 10.
            vlev_top = vlev_top * 10.
         ELSE IF (cfeld(ifld,ipl)(4:4).eq.'s') THEN
c		Convert from sigma level indices to
c		actual sigma values.
            WRITE(engplttl(ipl),'(F4.3," to ",F4.3,A)')
     &      sigh(NINT(vlev_bot)), sigh(NINT(vlev_top)),
     &      ' sig wind shr (u-comp.)'
            cvcorunits = ''
         ELSE
            WRITE(iup,*)"In subroutine field.f"
            WRITE(iup,*)"Unrecognized vertical coordinate "//
     +      "for u-comp. wind shear layer:"//cfeld(ifld,ipl)(4:4)
            WRITE(iup,*)"Stopping."
            STOP
         ENDIF

c
c
c	Subroutine wshear now calculates
c	wind shear using one of several  
c	vertical coordinate systems.
c
c	I appended arguments to the end of the
c	argument list for subroutine wshear.
c	One is a character indicating the
c	vertical coordinate system used in the 
c	wind shear layer specification.
c	It uses the same naming convention as does vcor,
c	with the additional option of 'a' for
c	height above ground level.
c
c	Current choices:
c		'a' = height (m AGL)
c		'z' = geopotential height (m MSL)
c		'p' = pressure (mb or hPa)
c		's' = sigma-surface (index)
c
c	Future choices:
c		't' = potential temp surface (index)
c
c
c	Another appended argument is a character string
c	holding the units of the vertical coordinate
c	system used for the wind shear layer specification:
c	cvcorunits.
c	Right now, this variable is only used to produce 
c	nice-looking warning messages in wshear.
c
c
c	Similar to the call to bshear, the first argument
c	to wshear is an integer that specifies
c	whether the u or v component of the wind shear
c	is returned in pl2.
c		1	=> u-comp.
c		not 1	=> v-comp.
c 


         call wshear (1,
     &                     uuu,vvv,ght,scr3a,ter,pl2,miy,mjx,mkzh,
     &                     vlev_bot,vlev_top,
     &                     cfeld(ifld,ipl)(4:4), cvcorunits )


c
c	According to Jim Bresch (MMM/NCAR)...
c
c	MM5 uses a staggered Arakawa-B grid.
c	Winds are defined on the so-called
c	"cross points" and the mass variables 
c	on the "dot points".
c	So, icdwk is the integer, cross-dot, work-array.
c	As you concluded, icdwk = 0 for wind computations.
c
c	indwk is the index into the correct 
c	position of the work array.
c	Following the pattern should be all right.
c

         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0

         unwk(ipl)='m s~S~-1~N~'

      elseif (cfeld(ifld,ipl)(1:3).eq.'vsh') then!  shear, m/s
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call pfcalc(sigh,sigf,pstx,prs,scr3a,miy,mjx,mkzh)

         read(cfeld(ifld,ipl)(5:10),'(2f3.0)') vlev_bot,vlev_top

         IF (cfeld(ifld,ipl)(4:4) .eq. 'a') THEN
            engplttl(ipl)=cfeld(ifld,ipl)(5:7)//''//' to '//
     &      cfeld(ifld,ipl)(8:10)//' hm AGL wind shr (v-comp.)'
            cvcorunits = 'm AGL'
            vlev_bot = vlev_bot * 100.
            vlev_top = vlev_top * 100.
         ELSE IF (cfeld(ifld,ipl)(4:4).eq.'z') THEN
            engplttl(ipl)=cfeld(ifld,ipl)(5:7)//''//' to '//
     &      cfeld(ifld,ipl)(8:10)//' hm MSL wind shr (v-comp.)'
            cvcorunits = 'm'
            vlev_bot = vlev_bot * 100.
            vlev_top = vlev_top * 100.
         ELSE IF (cfeld(ifld,ipl)(4:4).eq.'p') THEN
            engplttl(ipl)=cfeld(ifld,ipl)(5:7)//'0'//' to '//
     &      cfeld(ifld,ipl)(8:10)//'0 hPa wind shear (v-comp.)'
            cvcorunits = 'hPa'
            vlev_bot = vlev_bot * 10.
            vlev_top = vlev_top * 10.
         ELSE IF (cfeld(ifld,ipl)(4:4).eq.'s') THEN
            WRITE(engplttl(ipl),'(F4.3," to ",F4.3,A)')
     &      sigh(NINT(vlev_bot)), sigh(NINT(vlev_top)),
     &      ' sig wind shr (v-comp.)'
            cvcorunits = ''
         ELSE
            WRITE(iup,*)"In subroutine field.f"
            WRITE(iup,*)"Unrecognized vertical coordinate "//
     +      "for v-comp. wind shear layer:"//cfeld(ifld,ipl)(4:4)
            WRITE(iup,*)"Stopping."
            STOP
         ENDIF



         call wshear (2,
     &                        uuu,vvv,ght,scr3a,ter,pl2,miy,mjx,mkzh,
     &                        vlev_bot,vlev_top,
     &                        cfeld(ifld,ipl)(4:4), cvcorunits )
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0

         unwk(ipl)='m s~S~-1~N~'





c
c	End of new fields added by Dave Ahijevych Dec 2001
c

c
c   Following are new fields added by J.F. Bresch between Dec/97 and Mar/00
c
      elseif (cfeld(ifld,ipl)(1:4).eq.'tsfc' .or.    ! surface temp, deg C
     &        cfeld(ifld,ipl)(1:4).eq.'tsff') then   ! surface temp, deg F
         idimn(ipl)=2
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call readdat(iudata,fname,iendf1,'tgk       ',
     &      miy,mjx,mkzh,2,1,pl2,istat)
         if (iprog.le.3) then
            call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
            call readdat(iudata,fname,iendf1,'tmk_sfan  ',
     &         miy,mjx,mkzh,2,1,scr2a,istat)
         endif
         do j=1,mjx-1
         do i=1,miy-1
            if (iprog.le.3) then
               pl2(i,j)=scr2a(i,j)-celkel
            else
c
c            If air is warmer than ground, use air temp.
c            If air is colder than ground, use avg. of air temp.
c               and ground. temp.
c
               if (pl2(i,j) .le. tmk(i,j,mkzh)) then
                  pl2(i,j)= (tmk(i,j,mkzh)-celkel)
               else
                  pl2(i,j)=(0.5*(tmk(i,j,mkzh)+pl2(i,j))-celkel)
               endif
            endif
         enddo
         enddo
         engplttl(ipl)='Surface air temperature'
         unwk(ipl)='~S~o~N~C'
         if (cfeld(ifld,ipl)(1:4).eq.'tsff') then
           do j=1,mjx-1
           do i=1,miy-1
               pl2(i,j) = pl2(i,j) * 9./5. +32.
           enddo
           enddo
           unwk(ipl)='~S~o~N~F'
         endif
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
      elseif (cfeld(ifld,ipl)(1:4).eq.'tdf ') then! dewpoint, deg F
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call tdpcalc(qvp,prs,pl3,miy,mjx,mkzh)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
           pl3(i,j,k) = pl3(i,j,k) * 9./5. + 32.
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Dewpoint temperature'
         unwk(ipl)='~S~o~N~F'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tdd ') then! dewpoint depression, deg C
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call tdpcalc(qvp,prs,pl3,miy,mjx,mkzh)
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
           pl3(i,j,k) = tmk(i,j,k) - 273.16 - pl3(i,j,k) 
           if (pl3(i,j,k) .lt. 0. ) pl3(i,j,k) = 0.
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Dewpoint depression'
         unwk(ipl)='~S~o~N~C'
      elseif (cfeld(ifld,ipl)(1:4).eq.'tfp ') then! frostpoint, deg C
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call tfpcalc(qvp,prs,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Frostpoint temperature'
         unwk(ipl)='~S~o~N~C'
      elseif (cfeld(ifld,ipl)(1:4).eq.'thv ') then! virtual pot. temp., K
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         do k = 1, mkzh
         do j = 1, mjx-1
         do i = 1, miy-1
           tv=virtual(tmk(i,j,k),qvp(i,j,k))
           gammam=gamma*(1.+gammamd*qvp(i,j,k))
           pl3(i,j,k)=tv*(1000./prs(i,j,k))**gammam
         enddo
         enddo
         enddo
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Virtual potential temperature'
         unwk(ipl)='K'
      elseif (cfeld(ifld,ipl)(1:4).eq.'sreh') then! storm-rel helicity
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call relhl(uuu,vvv,ght,ter,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Storm-relative helicity 75%:30R'
         unwk(ipl)='m~S~2~N~ s~S~-2~N~'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:4).eq.'vgp') then! vorticity-gen. potential
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3c,wk,maxslab)
c
c      Put CAPE in mkzh slab of scr3a.
c
         call pfcalc(sigh,sigf,pstx,prs,scr3c,miy,mjx,mkzh)
         call capecalc3d(prs,tmk,qvp,ght,ter,scr3c,scr3a,scr3b,
     &      miy,mjx,mkzh,0)
c
c      Now put it in scr2a
c
         do j = 1, mjx-1
         do i = 1, miy-1
            scr2a(i,j)=scr3a(i,j,mkzh)
         enddo
         enddo
         call vgp (uuu,vvv,ght,ter,scr2a,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Vorticity generation potential'
         unwk(ipl)='m s~S~-2~N~'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:4).eq.'srfl') then! storm-rel low-level inflow
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call srflo(uuu,vvv,ght,ter,0,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Low-level storm-relative flow'
         unwk(ipl)='m s~S~-1~N~'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:4).eq.'srfh') then! storm-rel mid-level flow
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call srflo(uuu,vvv,ght,ter,1,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Mid-level storm-relative flow'
         unwk(ipl)='m s~S~-1~N~'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:4).eq.'uubs') then! bulk shear, m/s
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call bshear (1,uuu,vvv,ght,ter,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='0-6 km shear (x-comp.)'
         unwk(ipl)='m s~S~-1~N~'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:4).eq.'vvbs') then! bulk shear, m/s
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call bshear (2,uuu,vvv,ght,ter,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='0-6 km shear (y-comp.)'
         unwk(ipl)='m s~S~-1~N~'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:4).eq.'uusr') then! suprcell mvmnt x-comp., m/s
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call srflo4 (uuu,vvv,ght,ter,0,0,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Supercell motion (x-comp.)'
         unwk(ipl)='m s~S~-1~N~'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:4).eq.'vvsr') then! suprcell mvmnt y-comp., m/s
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call srflo4 (uuu,vvv,ght,ter,0,1,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=0
         engplttl(ipl)='Supercell motion (y-comp.)'
         unwk(ipl)='m s~S~-1~N~'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:3).eq.'sr9') then! storm-rel. flow at 9km, m/s
         call getpt(miy,mjx,mkzh,ifree,2,i_pl2,wk,maxslab)
         call srflo4 (uuu,vvv,ght,ter,1,0,pl2,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Supercell type (9-10 km rel. flow)'
         unwk(ipl)='m s~S~-1~N~'
         idimn(ipl)=2
      elseif (cfeld(ifld,ipl)(1:4).eq.'cat ') then! Turbulence index
         call getpt(miy,mjx,mkzh,ifree,3,i_pl3,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call turb(sigh,pstx,prs,uuu,vvv,ght,xmap,dmap,tmk,qvp,
     &     scr2a,scr2b,pl3,miy,mjx,mkzh)
         indwk(ifld,ipl)=incwk
         icdwk(ipl)=1
         engplttl(ipl)='Clear air turbulence index'
         unwk(ipl)='s~S~-2~N~'
c
c   End of new fields added by J.F. Bresch between Dec '97 and Mar '00
c
      else
         write(iup,*)'   Stopping in FIELDS: I don''t recognize',
     &      ' this field:'
         write(iup,*)'   <',cfeld(ifld,ipl),'>'
         stop
      endif
c
c   End of big "if" block
c
 1970 continue
c
      if (idimn(ipl).eq.3) then
         incwk=incwk+mkzh
      elseif (idimn(ipl).eq.2) then
         incwk=incwk+1
      else
         write(iup,*)'For ipl = ',ipl,'  idimn not 2 or 3.  idimn= ',
     &      idimn(ipl)
         stop
      endif
c
c   Do horizontal gradient, if asked for, and if field is 3D
c
      if (lgrad(ipl).or.llapl(ipl)) then
         ifree=incwk
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         ipass=0
 543     ipass=ipass+1
         if (idimn(ipl).eq.3) then
            call derivcz(pl3,icdwk(ipl),ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &         scr2a,scr2b,scr3a,icdwk(ipl),'x',miy,mjx,mkzh)
            if (ipass.eq.2.and.igdir(ipl).eq.361) call addorfill
     &         (scr3b,pl3,miy,mjx,mkzh,3,icdwk(ipl),1.,0.)
            call derivcz(pl3,icdwk(ipl),ght,xmap,dmap,pstx,qvp,tmk,
     &         sigh,scr2a,scr2b,scr3b,icdwk(ipl),'y',miy,mjx,mkzh)
         else
            if (icdwk(ipl).eq.0) then
               call ddx(pl2,icdwk(ipl),scr2a,xmap,dmap,icdwk(ipl),
     &              miy,mjx)
               if (ipass.eq.2.and.igdir(ipl).eq.361) call addorfill
     &            (scr2b,pl2,miy,mjx,mkzh,2,icdwk(ipl),1.,0.)
               call ddy(pl2,icdwk(ipl),scr2b,xmap,dmap,icdwk(ipl),
     &              miy,mjx)
            elseif (icdwk(ipl).eq.1) then
               call ddx(pl2,icdwk(ipl),scr2a,xmap,dmap,icdwk(ipl),
     &              miy,mjx)
               if (ipass.eq.2.and.igdir(ipl).eq.361) call addorfill
     &            (scr2b,pl2,miy,mjx,mkzh,2,icdwk(ipl),1.,0.)
               call ddy(pl2,icdwk(ipl),scr2b,xmap,dmap,icdwk(ipl),
     &              miy,mjx)
            endif
         endif
         if (igdir(ipl).ge.0.and.igdir(ipl).le.360) then ! gradient or
c                                            2nd deriv. in specified dir
            cosa=cos(rpd*(90-igdir(ipl)))
            sina=sin(rpd*(90-igdir(ipl)))
         elseif (igdir(ipl).eq.361) then  ! mag. of grad., or total lapl.
            cosa=0.
            sina=0.
         elseif (igdir(ipl).eq.362) then  ! grad. or 2nd drv along cross sec.
            cosa=xdist/xseclen
            sina=ydist/xseclen
         elseif (igdir(ipl).eq.363) then  ! grad. or 2nd drv perp to cross sec.
            cosa=-ydist/xseclen
            sina=xdist/xseclen
         endif
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            if (idimn(ipl).eq.3) then
               if (cosa.eq.0..and.sina.eq.0.) then
                  if (llapl(ipl).and.ipass.eq.1) then
                     pl3(i,j,k)=scr3a(i,j,k)
                  else
                     pl3(i,j,k)=sqrt(scr3a(i,j,k)*scr3a(i,j,k)+
     &                  scr3b(i,j,k)*scr3b(i,j,k))
                  endif
               else
                  pl3(i,j,k)=cosa*scr3a(i,j,k)+sina*scr3b(i,j,k)
               endif
            elseif (idimn(ipl).eq.2) then
               if (cosa.eq.0..and.sina.eq.0.) then
                  if (llapl(ipl).and.ipass.eq.1) then
                     pl2(i,j)=scr2a(i,j)
                  else
                     pl2(i,j)=sqrt(scr2a(i,j)*scr2a(i,j)+
     &                  scr2b(i,j)*scr2b(i,j))
                  endif
               else
                  pl2(i,j)=cosa*scr2a(i,j)+sina*scr2b(i,j)
               endif
            endif
         enddo
         enddo
         enddo
         if (llapl(ipl).and.ipass.eq.1) goto 543
      endif
c
c   Do advection, if asked for, and if field is 3D
c
      if (lhadv(ipl)) then
         if (idimn(ipl).ne.3) then
            write(iup,*)'Can only do advection for 3-D fields.'
         endif
         ifree=incwk
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2b,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3b,wk,maxslab)
         call derivcz(pl3,icdwk(ipl),ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,scr3a,icdwk(ipl),'x',miy,mjx,mkzh)
         call derivcz(pl3,icdwk(ipl),ght,xmap,dmap,pstx,qvp,tmk,sigh,
     &      scr2a,scr2b,scr3b,icdwk(ipl),'y',miy,mjx,mkzh)
         if (igdir(ipl).ge.0.and.igdir(ipl).le.360) then ! adv. in 
c                                                          specified dir.
            cosa=cos(rpd*(90-igdir(ipl)))
            sina=sin(rpd*(90-igdir(ipl)))
         elseif (igdir(ipl).eq.361) then  ! total adv.
            cosa=0.
            sina=0.
         elseif (igdir(ipl).eq.362) then  ! adv. along cross sec.
            cosa=xdist/xseclen
            sina=ydist/xseclen
         elseif (igdir(ipl).eq.363) then  ! adv. perp. to cross sec.
            cosa=-ydist/xseclen
            sina=xdist/xseclen
         endif
         do k=1,mkzh
         do j=1,mjx-1
         do i=1,miy-1
            uuucross=.25*(uuu(i,j,k)+uuu(i+1,j,k)+uuu(i,j+1,k)+
     &         uuu(i+1,j+1,k))
            vvvcross=.25*(vvv(i,j,k)+vvv(i+1,j,k)+vvv(i,j+1,k)+
     &         vvv(i+1,j+1,k))
            if (cosa.eq.0..and.sina.eq.0.) then
                  pl3(i,j,k)=
     &               -uuucross*scr3a(i,j,k)-vvvcross*scr3b(i,j,k)
            else
               pl3(i,j,k)=-(cosa*uuucross+sina*vvvcross)*
     &            (cosa*scr3a(i,j,k)+sina*scr3b(i,j,k))
            endif
         enddo
         enddo
         enddo
      endif
c
c   Do constant-pressure smoothing if asked for.
c
      if (ismcp(ipl).gt.0) then
         if (idimn(ipl).ne.3) then
            write(iup,*)'You can only do constant-pres.',
     &         ' smoothing of a 3-d field.'
            write(iup,*)'  ipl,ifld,cfeld=',ipl,ifld,cfeld(ifld,ipl)
            stop
         endif
         ifree=incwk
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
         call smoothcp(pl3,icdwk(ipl),
     &      scr3a,prs,pslab1,ismcp(ipl),miy,mjx,mkzh,
     &      mabpl,morpl)
      endif

c
c   Do time/model run differencing, if requested.
c
      if (cdiff(ipl)(1:5).ne.'none '.or.rdiff(ipl).ne.rmsg) then
c
      if (idiffpass.eq.2) goto 779
c
      ldiffsuccess=.true.
c
c   Save various parameters from present data set and time.
c
      xtime_sv=xtime
      casename_sv=casename
      iendc_sv=iendc
      fname_sv=fname
      iendf1_sv=iendf1
      nxtavl_sv=nxtavl
      do i=1,maxtavl
         cxtimeavl_sv(i)=cxtimeavl(i)
         xtimeavl_sv(i)=xtimeavl(i)
      enddo
      nxt_sv=nxt
c
c  Set new xtime
c
      if (rdiff(ipl).ne.rmsg) then
         if (ldfrl(ipl)) then
            xtime=xtime+rdiff(ipl)
         else
            xtime=rdiff(ipl)
         endif
      endif
c
c   Set new case name, load new available times
c
      if (cdiff(ipl)(1:5).ne.'none ') then
         casename=cdiff(ipl)
         iendc=index(casename,' ')-1
c
c      Temporarily use fname for new ".xtimes" file.
c
         fname=casename(1:iendc)//'.xtimes'
         open(unit=iuxtavl,file=fname,form='formatted',status='old')
         read(iuxtavl,*)nxtavl
         if (nxtavl.gt.maxtavl) then
            write(iup,*)'There are ',nxtavl,' times in the ".xtime"'
            write(iup,*)'file for the differencing data set, but'
            write(iup,*)'maxtavl is only ',maxtavl,'.  Increase maxtavl'
            write(iup,*)'in rip.f to be at least as big as ',nxtavl
            write(iup,*)'and then recompile and re-run rip.'
            stop
         endif
         do i=1,nxtavl
            read(iuxtavl,'(a9)') cxtimeavl(i)
            read(cxtimeavl(i),'(f9.5)') xtimeavl(i)
         enddo
         close (iuxtavl)
      endif
c
c   Determine if requested time is available in requested data set.
c
      iavail=0
      do i=1,nxtavl
         if (abs(xtimeavl(i)-xtime).le.tacch) then
            iavail=1
            nxt=i
            goto 640
         endif
      enddo
 640  continue
      if (iavail.eq.0) then
         write(iup,*)'   Requested time ',xtime,' is not available'
         write(iup,*)'   in the requested difference data set.'
         write(iup,*)'   RIP will plot undifferenced field.'
         ldiffsuccess=.false.  ! This is so difference info won't be
c                              ! printed in routine pltitle
c
c      If time was not available, restore some things and get out
c
         xtime=xtime_sv
         nxt=nxt_sv
         if (cdiff(ipl)(1:5).ne.'none ') then
            casename=casename_sv
            iendc=iendc_sv
            nxtavl=nxtavl_sv
            do i=1,maxtavl
               cxtimeavl(i)=cxtimeavl_sv(i)
               xtimeavl(i)=xtimeavl_sv(i)
            enddo
         endif
         goto 781
      endif
c
c   Make new fname
c
      iendf1=iendc+11
      fname=casename(1:iendc)//'_'//cxtimeavl(nxt)//'_'
c
c   Get new header information and check domain compatibility
c      (only if doing case differencing)
c
      if (cdiff(ipl)(1:5).ne.'none ') then
c
c      First read in the terrain file, without calling readdat,
c      so we can access the header information.
c
         fname(iendf1+1:)='ter'
         open (unit=iudata,file=fname,form='unformatted',
     &      status='unknown')
         read (iudata) vardesc,plchun,ihrip,rhrip,
     &      chrip,fullsigma,halfsigma
         close (iudata)
c
c      As far as compatibility between the original domain and the domain to
c      be subtracted, the only thing that we will require is that the
c      dimensions be the same in all three directions.  All other aspects
c      of the subtraction data set can (in principal) be different from the
c      original data set.
c
         il=ihrip(4)
         jl=ihrip(5)
         kl=ihrip(9)
         if(miy.ne.il.or.mjx.ne.jl.or.mkzh.ne.kl) then
            write(iup,*)'   Subtraction domain dimensions are ',
     &            il,jl,kl,'.'
            write(iup,*)'   This is incompatible with original domain'
            write(iup,*)'   whose dimensions are ',miy,mjx,mkzh
            write(iup,*)'   RIP will plot undifferenced field.'
            ldiffsuccess=.false.  ! This is so difference info won't be
c                                 ! printed in routine pltitle
c
c         If domain was incompatible, restore some things and get out
c
            xtime=xtime_sv
            casename=casename_sv
            iendc=iendc_sv
            nxtavl=nxtavl_sv
            do i=1,maxtavl
               cxtimeavl(i)=cxtimeavl_sv(i)
               xtimeavl(i)=xtimeavl_sv(i)
            enddo
            nxt=nxt_sv
            iendf1=iendf1_sv
            fname=fname_sv
            fname(iendf1+1:)='ter'
            open (unit=iudata,file=fname,form='unformatted',
     &         status='unknown')
            read (iudata) vardesc,plchun,ihrip,rhrip,
     &         chrip,fullsigma,halfsigma
            close (iudata)
            goto 781
         endif
c
c      Get information from new header record
c
         call getheadinfo(ihrip,rhrip,fullsigma,
     &      halfsigma,chrip,vardesc,plchun,
     &      sigf,sigh,nproj,miycors,mjxcors,inhyd,
     &      mdateb,mhourb,iice,iprog,ilandset,iwater,truelat1,truelat2,
     &      xlatc,xlonc,dskmc,dskm,yicorn,xjcorn,ptop,
     &      refslp,refslt,reflaps,refstratt,rhourb,dsc,ds,refrat,mkzh)
c
      endif
c
c   Read the basic fields.
c
      call readbasic(uuu,vvv,tmk,qvp,www,prs,ght,prs_tsf,ter,
     &   dmap,xmap,cor,xlus,pstx,pstd,sigh,sigf,ptop,inhyd,iprog,
     &   fname,iendf1,iudata,miy,mjx,mkzh)
c
c   Load/calculate field from new time/model run
c
      goto 35
c
 779  continue
c
c   Subtract the fields
c
      if (idimn(ipl).eq.3) then
c
c      New field is at slab marker incwk-mkzh, and is already in
c      pointed array pl3.  Old field is at slab marker incwk-2*mkzh.
c      Set up pointer so that old field is in pointed array scr3a.
c
         ifree=incwk-2*mkzh
         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)
c
c      Put subtracted field in scr3a
c
         do k=1,mkzh
         do j=1,mjx-icdwk(ipl)
         do i=1,miy-icdwk(ipl)
            scr3a(i,j,k)=scr3a(i,j,k)-pl3(i,j,k)
         enddo
         enddo
         enddo
c
c      Adjust slab markers
c
         indwk(ifld,ipl)=incwk-2*mkzh  ! where subtracted data is now located
         incwk=incwk-mkzh ! Reset to end of subtracted data.
c      Ready for next field.  Data from new time/model run will be overwritten
c
      elseif (idimn(ipl).eq.2) then
c
c      New field is at slab marker incwk-1, and is already in
c      pointed array pl2.  Old field is at slab marker incwk-2.
c      Set up pointer so that old field is in pointed array scr2a.
c
         ifree=incwk-2
         call getpt(miy,mjx,mkzh,ifree,2,i_scr2a,wk,maxslab)
c
c      Put subtracted field in scr2a
c
         do j=1,mjx-icdwk(ipl)
         do i=1,miy-icdwk(ipl)
            scr2a(i,j)=scr2a(i,j)-pl2(i,j)
         enddo
         enddo
c
c      Adjust slab markers
c
         indwk(ifld,ipl)=incwk-2  ! where subtracted data is now located
         incwk=incwk-1 ! Reset to end of subtracted data.
c      Data from new time/model run will be overwritten by next field.
c
      endif
c
c   Restore original parameters and data.
c
      xtime=xtime_sv
      nxt=nxt_sv
      fname=fname_sv
      iendf1=iendf1_sv
      if (cdiff(ipl)(1:5).ne.'none ') then
         casename=casename_sv
         iendc=iendc_sv
         nxtavl=nxtavl_sv
         do i=1,maxtavl
            cxtimeavl(i)=cxtimeavl_sv(i)
            xtimeavl(i)=xtimeavl_sv(i)
         enddo
         fname(iendf1+1:)='ter'
         open (unit=iudata,file=fname,form='unformatted',
     &      status='unknown')
         read (iudata) vardesc,plchun,ihrip,rhrip,
     &      chrip,fullsigma,halfsigma
         close (iudata)
         call getheadinfo(ihrip,rhrip,fullsigma,
     &      halfsigma,chrip,vardesc,plchun,
     &      sigf,sigh,nproj,miycors,mjxcors,inhyd,
     &      mdateb,mhourb,iice,iprog,ilandset,iwater,truelat1,truelat2,
     &      xlatc,xlonc,dskmc,dskm,yicorn,xjcorn,ptop,
     &      refslp,refslt,reflaps,refstratt,rhourb,dsc,ds,refrat,mkzh)
      endif
      call readbasic(uuu,vvv,tmk,qvp,www,prs,ght,prs_tsf,ter,
     &   dmap,xmap,cor,xlus,pstx,pstd,sigh,sigf,ptop,inhyd,iprog,
     &   fname,iendf1,iudata,miy,mjx,mkzh)
c
 781  continue
c
      endif  ! end of time/model run differencing
c
 1980 continue
c
c   Add the previous field(s) to this field if:
c	you are past the initial feld= line
c	the current feld= line has no addf request
c	the previous feld= line has a nonzero addf value
c
      if (ipl.gt.iplstrt.and.raddf(ipl).eq.0.0) then
      if (raddf(ipl-1).ne.0.0) then
         if (cfeld(ifld,ipl)(1:2).ne.'se'.and.
     &       cfeld(ifld,ipl)(1:2).ne.'sm') then
            icda=icdwk(ipl)
         else
            icda=0
         endif
         do ipla=ipl-1,iplstrt,-1
            if (raddf(ipla).eq.0.0) then
               goto 1900
            elseif (raddf(ipla).eq.1.) then
               do k=1,mkzh
                  kpl=indwk(ifld,ipl)-1+k
                  kpla=indwk(ifld,ipla)-1+k
               do j=1,mjx-icda
               do i=1,miy-icda
                  wk(i,j,kpl)=wk(i,j,kpl)+
     &               wk(i,j,kpla)
               enddo
               enddo
               enddo
            elseif (raddf(ipla).eq.-1.) then
               do k=1,mkzh
                  kpl=indwk(ifld,ipl)-1+k
                  kpla=indwk(ifld,ipla)-1+k
               do j=1,mjx-icda
               do i=1,miy-icda
                  wk(i,j,kpl)=wk(i,j,kpl)-
     &               wk(i,j,kpla)
               enddo
               enddo
               enddo
            else
               do k=1,mkzh
                  kpl=indwk(ifld,ipl)-1+k
                  kpla=indwk(ifld,ipla)-1+k
               do j=1,mjx-icda
               do i=1,miy-icda
                  wk(i,j,kpl)=wk(i,j,kpl)+
     &               raddf(ipla)*wk(i,j,kpla)
               enddo
               enddo
               enddo
            endif
         enddo
 1900    continue
      endif
      endif
c
c   Multiply by the previous field if asked for.
c
      if (ipl.gt.iplstrt.and.rmply(ipl).eq.0.0) then
      if (rmply(ipl-1).ne.0.0) then
         if (cfeld(ifld,ipl)(1:2).ne.'se'.and.
     &       cfeld(ifld,ipl)(1:2).ne.'sm') then
            icda=icdwk(ipl)
         else
            icda=0
         endif
         do ipla=ipl-1,iplstrt,-1
            if (rmply(ipla).eq.0.0) then
               goto 1901
            elseif (rmply(ipla).eq.1.) then
               do k=1,mkzh
                  kpl=indwk(ifld,ipl)-1+k
                  kpla=indwk(ifld,ipla)-1+k
               do j=1,mjx-icda
               do i=1,miy-icda
                  wk(i,j,kpl)=wk(i,j,kpl)*
     &               wk(i,j,kpla)
               enddo
               enddo
               enddo
            elseif (rmply(ipla).eq.-1.) then
               do k=1,mkzh
                  kpl=indwk(ifld,ipl)-1+k
                  kpla=indwk(ifld,ipla)-1+k
               do j=1,mjx-icda
               do i=1,miy-icda
                  IF(wk(i,j,kpla).eq.0.) THEN
                     WRITE(iup,*)'set undef quotient to 0 '//
     &                           'at i,j,kpla',i,j,kpla

                     wk(i,j,kpl)=0.
                  ELSE
                     wk(i,j,kpl)=wk(i,j,kpl)/ wk(i,j,kpla)
                  END IF
               enddo
               enddo
               enddo
            else
               do k=1,mkzh
                  kpl=indwk(ifld,ipl)-1+k
                  kpla=indwk(ifld,ipla)-1+k
               do j=1,mjx-icda
               do i=1,miy-icda
                  IF(rmply(ipla).LT.0. .AND. wk(i,j,kpla).eq.0.)THEN
                     WRITE(iup,*)'cannot divide by zero'
                     WRITE(iup,*)'Stopping.'
                     STOP
                  END IF
                  wk(i,j,kpl)=wk(i,j,kpl)*
     &               wk(i,j,kpla)**rmply(ipla)
               enddo
               enddo
               enddo
            endif
         enddo
 1901    continue
      endif
      endif
c
c This thresholding section was added by David Ahijevych (MMM)
c Feb 2002
c
c Threshold the field values if fdmn or fdmx are present in the feld= line.
c


      IF(rfdmn(ipl).ne.rmsg .or. rfdmx(ipl).ne.rmsg) then
         WRITE(iup,*)'Thresholding '//cfeld(ifld,ipl)
         WRITE(iup,*)'rfdmn,rfdmx=',rfdmn(ipl),rfdmx(ipl)



         if (cfeld(ifld,ipl)(1:2).ne.'se'.and.
     &       cfeld(ifld,ipl)(1:2).ne.'sm') then
            icda=icdwk(ipl)
         else
            icda=0
         endif
         do k=1,mkzh
            kpl=indwk(ifld,ipl)-1+k
            do j=1,mjx-icda
            do i=1,miy-icda
               IF(rfdmn(ipl).ne.rmsg)
     &           wk(i,j,kpl) = AMAX1(wk(i,j,kpl),rfdmn(ipl))
               IF(rfdmx(ipl).ne.rmsg)
     &           wk(i,j,kpl) = AMIN1(wk(i,j,kpl),rfdmx(ipl))
            enddo
            enddo

         enddo
            
c
c Write field name with the min and max thresholded values. 
c	
         WRITE(engplttl(ipl),'(F6.2,"<=",A,"<=",F6.2)')
     &      rfdmn(ipl),cfeld(ifld,ipl),rfdmx(ipl)
      ENDIF


c
c   Save the field to a file, if asked for.
c
      if (csave(ipl).ne.'dontsave  ') then
         vardesc='unknown'
         if (csave(ipl).eq.'sameasfeld') then
            varname=cfeld(ifld,ipl)
            plchun=unwk(ipl)
         else
            varname=csave(ipl)
            plchun='??'
         endif

c
c Get a 3-d scratch array and copy wk array elements into it, one-by-one.
c scr3a is the array that will be saved to file.
c If necessary, smooth scr3a, but leave wk alone.

         call getpt(miy,mjx,mkzh,ifree,3,i_scr3a,wk,maxslab)

c I don't know what icda means, but it's used above when
c dealing with raddf.

         if (cfeld(ifld,ipl)(1:2).ne.'se'.and.
     &       cfeld(ifld,ipl)(1:2).ne.'sm') then
            icda=icdwk(ipl)
         else
            icda=0
         endif
         do k=1,mkzh
            kpl=indwk(ifld,ipl)-1+k
            kpla=indwk(ifld,ipla)-1+k
            do j=1,mjx-icda
            do i=1,miy-icda
               scr3a(i,j,kpl)=wk(i,j,kpl)
            enddo
            enddo
c
c If first 4 characters of csave(ipl) are 'beta',
c smooth level kpl of scr3a.
c
c This ability to smooth the beta field before saving it to a file
c was added by
c David Ahijevych (MMM) Jan 2002
c

            if(csave(ipl)(1:4).eq.'beta') then
               call smooth(scr3a(:,:,kpl),ismth(ipl),miy,
     &                     miy-1,mjx-1) 
            endif

         enddo
            

         call writefile(scr3a(1,1,indwk(ifld,ipl)),varname,
     &      0,idimn(ipl),icdwk(ipl),vardesc,plchun,fname,iendf1,
     &      ihrip,rhrip,chrip,fullsigma,halfsigma,miy,mjx,mkzh)
      endif
c
 2000 continue
c
c   Set all levels to 1 for 2-d variables.
c
      if (idimn(ipl).eq.2) then
         do ilev=1,maxlev
            rlevl(ilev,ipl)=1
            rlavl(ilev,ipl)=1
         enddo
      endif
c
      return
      end