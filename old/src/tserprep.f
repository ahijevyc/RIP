c                                                                     c
c*********************************************************************c
c                                                                     c
      subroutine tserprep(tseryi,tserxj,tser3lc,tserloc,
     &   ntsers,ntsert,maxtsers,rip_root)
c
c   This subroutine reads the time series station locations and
c      calculates the corresponding x/y locations on the grid.
c
      dimension tseryi(maxtsers),tserxj(maxtsers)
      character tser3lc(maxtsers)*4,tserloc(maxtsers)*21,
     &   rip_root*256
c
      dimension slat(3000),slon(3000)
      character threelc(3000)*4,loc(3000)*21,threelctser*4,fnm*256
c
      include 'comconst'
c
      i=0
      iendci=index(rip_root,' ')-1
      fnm=rip_root(1:iendci)//'/stationlist'
      open (unit=iustnlist,file=fnm,form='formatted',status='old')
      read(iustnlist,*)
      read(iustnlist,*)
  200 i=i+1
      read(iustnlist,'(a21,2x,a4,5x,f6.2,f8.2)',end=205)
     &   loc(i),threelc(i),slat(i),slon(i)
      if (threelc(i)(1:1) .eq. ' ') threelc(i)(1:1) = 'K'
      goto 200
  205 nstot=i-1
      close (iustnlist)
c
      i=0
      open (unit=iutserstn,file='tserstn.dat',form='formatted',
     &   status='old')
  300 read(iutserstn,'(a4)',end=350) threelctser
      if (threelctser(1:1) .eq. ' ') threelctser(1:1) = 'K'
      do 310 j=1,nstot
         if (threelctser.eq.threelc(j)) then
            i=i+1
            tser3lc(i)=threelctser
            tserloc(i)=loc(j)
            call maptform(syy,sxx,slat(j),slon(j),-1)
            tseryi(i)=refrat*(syy-yicorn)+1.
            tserxj(i)=refrat*(sxx-xjcorn)+1.
            goto 315
         endif
  310 continue
  315 goto 300
  350 ntsers=i
      ntsert=0
      close (iutserstn)
      return
      end
