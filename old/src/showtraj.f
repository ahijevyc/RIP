      program showtraj
c
      parameter (maxtraj=1000,maxtrajtime=200)
c
      dimension stortr(maxtrajtime,maxtraj,3)
      character fname*256
c
c   RIP header variables
c
      dimension ihrip(32),rhrip(32),fullsigma(128),halfsigma(128)
      character chrip(64)*64,vardesc*64,plchun*24
c
      print*,'Enter the trajectory file name:'
      read(*,'(a256)') fname
      open (unit=10,file=fname,form='unformatted',status='old')
      open (unit=11,file='trajprint',form='formatted',status='unknown')
      read(10) vardesc,plchun,ihrip,rhrip,chrip,
     &      fullsigma,halfsigma
      read (10) rtim,ctim,dttraj,ntraj
      print*,'rtim,ctim,dttraj,ntraj=',rtim,ctim,dttraj,ntraj
      ntrajtime=nint(abs(rtim-ctim)/dttraj*3600) + 1
      if (rtim.lt.ctim) then
         trendtime=ctim
         trbegtime=rtim
         itm1=1
         itm2=ntrajtime
         itmi=1
      else
         trendtime=rtim
         trbegtime=ctim
         itm1=ntrajtime
         itm2=1
         itmi=-1
      endif
      do itm=itm1,itm2,itmi
         read(10) (stortr(itm,itr,1),itr=1,ntraj),
     &       (stortr(itm,itr,2),itr=1,ntraj),
     &       (stortr(itm,itr,3),itr=1,ntraj)
      enddo
      do itr=1,ntraj
         print*
         print*, 'Trajectory number ',itr,':'
         do itm=1,ntrajtime
            time=trbegtime+(itm-1)*dttraj/3600.
            xnest=1.+(stortr(itm,itr,2)-rhrip(8))*rhrip(5)/rhrip(6)
            ynest=1.+(stortr(itm,itr,1)-rhrip(7))*rhrip(5)/rhrip(6)
            write(6,'(a5,f7.3,a5,f7.3,a5,f7.3,a5,f7.5)')'time=',time,
     &         '   x=',xnest,'   y=',ynest,'   s=',stortr(itm,itr,3)
         enddo
      enddo
      stop
      end

