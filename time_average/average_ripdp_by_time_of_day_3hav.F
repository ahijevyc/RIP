! pgf90 -Mfree -byteswapio \
!        kwm_date_utilities.F readrip.F average_ripdp_by_time_of_day_3hav.F

program average_ripdp_by_time_of_day
  use ripdp_module
  use kwm_date_utilities
  implicit none

  character(len=256),parameter :: ripcase = "/aquula_d5/kmanning/EDAS_WRFIN/edas"

  character(len=13), parameter :: startdate = "2002-06-09_12"
  character(len=13), parameter :: enddate   = "2002-06-21_12"

  character(len=13) :: nowdate
  character(len=8) :: hdate
  integer :: idate

  integer :: ix, jx, kx, idt, ihr, ierr, i, ndim
  real, allocatable, dimension(:,:,:) :: avg2d
  real, allocatable, dimension(:,:,:,:) :: avg3d
  real, allocatable, dimension(:,:) :: tmp2d
  real, allocatable, dimension(:,:,:) :: tmp3d
  real, dimension(0:23) :: acount

  character(len=256) :: fieldname
  character(len=256) :: subcase
  character(len=64) :: vardesc_hold

  ! Command-line argument:  The name of the field you want to average.
  call getarg(1, fieldname)

  ! Read a ripdp header to get dimensions
  call readrip_header(0, trim(fieldname), trim(ripcase))
  ix = ihrip(4)
  jx = ihrip(5)
  ndim = ihrip(6)
  if (ndim == 2) then
     kx = 1
  else if (ndim == 3) then
     kx = ihrip(9)
  else
     write(*, '("Unrecognized number of dimensions: ", I)') ndim
     stop
  endif
  print*, 'ndim, ix, jx, kx = ', ndim, ix, jx, kx

  do i = 1, 32
     if (ihrip(i) < 999999998) then
        write(*, '(I2, 1x, I8, 1x, A)') i, ihrip(i), trim(chrip(i))
     endif
  enddo
  do i = 1, 32
     if (rhrip(i) < 999999998) then
        write(*,'(I2,1x,G15.4, 1x, A)') i, rhrip(i), trim(chrip(i+32))
     endif
  enddo

  if (ndim == 2) then
     allocate(tmp2d(ix,jx))
     allocate(avg2d(ix,jx,0:23))
     avg2d = 0.0
  else if (ndim == 3) then
     allocate(tmp3d(ix,jx,kx))
     allocate(avg3d(ix,jx,kx,0:23))
     avg2d = 0.0
  endif
  acount = 0.0

  call geth_newdate(nowdate, startdate, -3)
  TIMELOOP : do while ( nowdate < enddate )
     call geth_newdate(nowdate, nowdate, 3)
     call geth_idts(nowdate, startdate, idt)

     ! Skip the first day:
     ! if (idt < 25) then
     !    cycle TIMELOOP
     !endif
     
     read(nowdate(12:13), *) ihr

     if (ndim == 2) then
        call readrip(idt, trim(fieldname), trim(ripcase), tmp2d, ix, jx, ierr)
     else if (ndim == 3) then
        call readrip(idt, trim(fieldname), trim(ripcase), tmp3d, ix, jx, kx, &
             ierr)
     endif
     if (ierr /= 0) then
        idt = idt - 1
        exit TIMELOOP
     endif
     
     if (ndim == 2) then
        avg2d(:,:,ihr) = avg2d(:,:,ihr) + tmp2d
     else if (ndim == 3) then
        avg3d(:,:,:,ihr) = avg3d(:,:,:,ihr) + tmp3d
     endif
     acount(ihr) = acount(ihr) + 1.0

  enddo TIMELOOP

  do ihr = 0, 23, 3

     call geth_newdate(nowdate, startdate, ihr)
     hdate = nowdate(3:4)//nowdate(6:7)//nowdate(9:10)//nowdate(12:13)
     read(hdate,*) ihrip(11)
     rhrip(15) = ihr

     vardesc_hold = vardesc
     vardesc = "Average "//trim(vardesc)//", by hour."

     subcase = ripcase(index(ripcase,"/",.TRUE.)+1:len_trim(ripcase))
     if (ndim == 2) then
        avg2d(:,:,ihr) = avg2d(:,:,ihr) / acount(ihr)
        call writerip(ihr, trim(fieldname)//"A3H", trim(ripcase), &
             avg2d(:,:,ihr), ix, jx, ierr)

!KWM        call writerip(ihr, trim(fieldname)//"AVG", trim(subcase), &
!KWM             avg2d(:,:,ihr), ix, jx, ierr)

     else if (ndim == 3) then
        avg3d(:,:,:,ihr) = avg3d(:,:,:,ihr) / acount(ihr)

        call writerip(ihr, trim(fieldname)//"A3H", trim(ripcase), &
             avg3d(:,:,:,ihr), ix, jx, kx, ierr)

!KWM        call writerip(ihr, trim(fieldname)//"_AVG", trim(subcase), &
!KWM             avg3d(:,:,:,ihr), ix, jx, kx, ierr)

     endif
     vardesc = vardesc_hold


  enddo

  



end program average_ripdp_by_time_of_day

















