! pgf90 -Mfree -byteswapio \
!        kwm_date_utilities.F readrip.F average_specified_ripdp_files.F

program average_specified_ripdp_files
  use ripdp_module
  use kwm_date_utilities
  implicit none

  character(len=256) :: ripcase
  character(len=256) :: flnm

  character(len=8) :: hdate
  integer :: idate

  integer :: ix, jx, kx, idt, ihr, ierr, i, ndim
  real :: xdt
  real, allocatable, dimension(:,:) :: avg2d
  real, allocatable, dimension(:,:,:) :: avg3d
  real, allocatable, dimension(:,:) :: tmp2d
  real, allocatable, dimension(:,:,:) :: tmp3d
  real :: acount

  character(len=256) :: fieldname
  character(len=256) :: subcase
  character(len=64) :: vardesc_hold

  integer, external :: iargc
  integer :: iarg, numarg
  real    :: min_hour, max_hour
  real    :: increment
  character(len=40) :: txt
  
  numarg = iargc()
  if (numarg /= 5) call helpstop

  call getarg(1, ripcase)
  if ((ripcase == "-h") .or. (ripcase == "-help") .or. (ripcase == "help") .or. (ripcase == "--help")) call helpstop
  call getarg(2, fieldname)

  call getarg(3, txt)
  read(txt, *, iostat=ierr) min_hour
  if (ierr /= 0) then
     write(*,'(/,5x, "***** Problem with command-line argument <min_hour> *****")') 
     call helpstop
  endif

  call getarg(4, txt)
  read(txt, *, iostat=ierr) max_hour
  if (ierr /= 0) then
     write(*,'(/,5x, "***** Problem with command-line argument <max_hour> *****")') 
     call helpstop
  endif

  call getarg(5, txt)
  read(txt, *, iostat=ierr) increment
  if (ierr /= 0) then
     write(*,'(/,5x, "***** Problem with command-line argument <increment> *****")') 
     call helpstop
  endif

  if (min_hour >= max_hour) then
     write(*,'(/,5x, "***** <min_hour> must be less than <max_hour> *****")') 
     write(*,'(5x, "*****   You have requested:  ")')
     write(*,'(5x, "*****       <min_hour> = ", F20.10)') min_hour
     write(*,'(5x, "*****       <max_hour> = ", F20.10)') max_hour
     call helpstop
  endif

!KWM  if ( mod((max_hour-min_hour), increment) /= 0 ) then
!KWM     write(*,'(/,5x, "***** <min_hour> to <max_hour> does not fit with <increment> *****")') 
!KWM     write(*,'(5x, "*****   You have requested:  ")')
!KWM     write(*,'(5x, "*****       <min_hour> = ", F20.10)') min_hour
!KWM     write(*,'(5x, "*****       <max_hour> = ", F20.10)') max_hour
!KWM     write(*,'(5x, "*****       <increment> = ", F20.10)') increment
!KWM     call helpstop
!KWM  endif

  write(txt, '("from ", F8.4, " to ", F8.4, " by ", F8.6)') min_hour, max_hour, increment

  ! Read a ripdp header to get dimensions
  call readrip_header(min_hour, trim(fieldname), trim(ripcase))
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

#ifdef _DEBUG_
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
#endif

  if (ndim == 2) then
     allocate(tmp2d(ix,jx))
     allocate(avg2d(ix,jx))
     avg2d = 0.0
  else if (ndim == 3) then
     allocate(tmp3d(ix,jx,kx))
     allocate(avg3d(ix,jx,kx))
     avg3d = 0.0
  endif
  acount = 0.0

  xdt = min_hour
  FILELOOP : do 
     if (xdt > max_hour+0.02) exit FILELOOP

     print*, 'Reading data for hour ', xdt
     if (ndim == 2) then
        call readrip(xdt, trim(fieldname), trim(ripcase), tmp2d, ix, jx, ierr)
     else if (ndim == 3) then
        call readrip(xdt, trim(fieldname), trim(ripcase), tmp3d, ix, jx, kx, ierr)
     endif

     if (ierr /= 0) then
        stop "Problem reading data."
     endif

     if (ndim == 2) then
        avg2d(:,:) = avg2d(:,:) + tmp2d
     else if (ndim == 3) then
        avg3d(:,:,:) = avg3d(:,:,:) + tmp3d
        print*, "(350, 210, 20):  ", tmp3d(350,210,20), avg3d(350,210,20)
     endif
     acount = acount + 1.0

     xdt = xdt + increment

  enddo FILELOOP

  vardesc_hold = vardesc
  vardesc = "Average "//trim(vardesc)//" "//trim(txt)
  print*, 'vardesc = ', vardesc
  

  subcase = ripcase(index(ripcase,"/",.TRUE.)+1:len_trim(ripcase))
  ihr = 0

  if (ndim == 2) then

     avg2d = avg2d / acount
     call writerip(ihr, trim(fieldname)//"5", trim(ripcase), avg2d, ix, jx, ierr)

  else if (ndim == 3) then

     avg3d = avg3d / acount
        print*, "(350, 210, 20):  ", avg3d(350,210,20)
     call writerip(ihr, trim(fieldname)//"5", trim(ripcase), avg3d, ix, jx, kx, ierr)

  endif

end program average_specified_ripdp_files

subroutine helpstop
  implicit none
  character(len=256) :: cmd
  call getarg(0, cmd)

  print*
  print*, 'usage:  '//trim(cmd)//' <ripcase> <fieldname> <min_hour> <max_hour> <increment>'
  print*
  stop
end subroutine helpstop















