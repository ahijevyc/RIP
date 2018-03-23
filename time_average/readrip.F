module ripdp_module
  implicit none

  ! RIPDP variables:
  character(len=64) :: vardesc
  character(len=24) :: plchun
  integer, dimension(32) :: ihrip
  real, dimension(32) :: rhrip
  character(len=64), dimension(64) :: chrip

  integer, parameter :: nxtimes = 100000

  interface readrip
     module procedure readrip_2d, readrip_3d, freadrip_2d, freadrip_3d
  end interface

  interface readrip_header
     module procedure ireadrip_header, freadrip_header
  end interface

  interface writerip
     module procedure writerip_2d, writerip_3d
  end interface

contains

  subroutine freadrip_2d(xhour, fldname, casename, fld, ix, jx, ierr)
    implicit none
    real,                               intent(in)  :: xhour
    character(len=*),                   intent(in)  :: fldname
    character(len=*),                   intent(in)  :: casename
    integer,                            intent(in)  :: ix, jx
    real,             dimension(ix,jx), intent(out) :: fld
    integer,                            intent(out) :: ierr

    character(len=256) :: flnm

    real, dimension(:,:,:), allocatable :: dum
    integer :: i,j
    integer :: rderr
    real, dimension(nxtimes) :: xtimes
    real, parameter :: epsilon = 0.02 ! hours

    ierr = 0

    ! Read the xtimes file to get a list of available xhours
    call read_xtimes(casename, xtimes, nxtimes)

    ! Find the xtimes that matches (within some <epsilon> hours) the requested xhour
    SEARCHLOOP : do i = 1, nxtimes
       if (xtimes(i) < -1.E25) then
          write(*,'(/,"Requested xhour:  ", F20.10, " not found in ",/,"    xtimes file ''",A,"''",/)') xhour, trim(casename)//".xtimes"
          stop "TIME NOT FOUND"
       endif
       if ( ABS(xtimes(i)-xhour) < epsilon ) then
          ! Build ripdp file names
          write(flnm, '(A,"_",F10.5,"_",A)') casename, xtimes(i), fldname
          ! Blanks may be created by the F10.5 specifier.  Fill replace those blanks with "0"
          do j = 1, len_trim(flnm)
             if (flnm(j:j) == " ") flnm(j:j) = "0"
          enddo
          exit SEARCHLOOP
       endif
    enddo SEARCHLOOP

    open(10, file=trim(flnm), status='old', form='unformatted', action='read', iostat=rderr)
    if (rderr /= 0) then
       print*, 'File not found: "'//trim(flnm)//'"'
       ! stop "READRIP: File not found."
       ierr = 1
       close(10)
       return
    endif

    read(10, iostat=rderr) vardesc, plchun, ihrip, rhrip, chrip
    if(rderr /= 0) stop "READRIP:  Read Problem 4."

    if (ihrip(6) == 2) then

       read(10, iostat=rderr) fld
       if (rderr /= 0) stop "READRIP:  Read Problem 5."

    elseif (ihrip(6) == 3) then

       stop "READRIP:  readrip_2d trying to read 3d field"

    else
       print*, 'READRIP:  Not a 2d or 3d field:  "'//trim(fldname)//'"'
       stop
    endif

    close(10)

  end subroutine freadrip_2d

  subroutine readrip_2d(hour, fldname, casename, fld, ix, jx, ierr)
    implicit none
    integer,                            intent(in)  :: hour
    character(len=*),                   intent(in)  :: fldname
    character(len=*),                   intent(in)  :: casename
    integer,                            intent(in)  :: ix, jx
    real,             dimension(ix,jx), intent(out) :: fld
    integer,                            intent(out) :: ierr

    character(len=256) :: flnm

    real, dimension(:,:,:), allocatable :: dum
    integer :: i
    integer :: rderr

    ierr = 0

    ! Build ripdp file names
    write(flnm, '(A,"_",I4.4,".00000_",A)') &
         casename, hour, fldname

    open(10, file=trim(flnm), status='old', form='unformatted', &
         action='read', iostat=rderr)
    if (rderr /= 0) then
       print*, 'File not found: "'//trim(flnm)//'"'
       ! stop "READRIP: File not found."
       ierr = 1
       return
    endif

    read(10, iostat=rderr) vardesc, plchun, ihrip, rhrip, chrip
    if(rderr /= 0) stop "READRIP:  Read Problem 4."

    if (ihrip(6) == 2) then

       read(10, iostat=rderr) fld
       if (rderr /= 0) stop "READRIP:  Read Problem 5."

    elseif (ihrip(6) == 3) then

       stop "READRIP:  readrip_2d trying to read 3d field"
!KWM       allocate(dum(ihrip(4),ihrip(5),ihrip(9)))
!KWM       read(10,iostat=rderr) dum
!KWM       if(rderr /= 0) stop "READRIP:  Read Problem 6."
!KWM
!KWM       ! The lowest model level is the last slab of the 3-d array
!KWM       fld = dum(:,:,ihrip(9))
!KWM
!KWM       deallocate(dum)

    else
       print*, 'READRIP:  Not a 2d or 3d field:  "'//trim(fldname)//'"'
       stop
    endif

    close(10)
  end subroutine readrip_2d

  subroutine freadrip_3d(xhour, fldname, casename, fld, ix, jx, kx, ierr)
    implicit none
    real,                                  intent(in)  :: xhour
    character(len=*),                      intent(in)  :: fldname
    character(len=*),                      intent(in)  :: casename
    integer,                               intent(in)  :: ix, jx, kx
    real,             dimension(ix,jx,kx), intent(out) :: fld
    integer,                               intent(out) :: ierr

    character(len=256) :: flnm

    integer :: i, j
    integer :: rderr
    real, dimension(nxtimes) :: xtimes
    real, parameter :: epsilon = 0.02 ! hours

    ierr = 0

    ! Read the xtimes file to get a list of available xhours
    call read_xtimes(casename, xtimes, nxtimes)

    ! Find the xtimes that matches (within some <epsilon> hours) the requested xhour
    SEARCHLOOP : do i = 1, nxtimes
       if (xtimes(i) < -1.E25) then
          write(*,'(/,"Requested xhour:  ", F20.10, " not found in ",/,"    xtimes file ''",A,"''",/)') xhour, trim(casename)//".xtimes"
          stop "TIME NOT FOUND"
       endif
       if ( ABS(xtimes(i)-xhour) < epsilon ) then
          ! Build ripdp file names
          write(flnm, '(A,"_",F10.5,"_",A)') casename, xtimes(i), fldname
          ! Blanks may be created by the F10.5 specifier.  Fill replace those blanks with "0"
          do j = 1, len_trim(flnm)
             if (flnm(j:j) == " ") flnm(j:j) = "0"
          enddo
          exit SEARCHLOOP
       endif
    enddo SEARCHLOOP

    open(10, file=trim(flnm), status='old', form='unformatted',  action='read', iostat=rderr)
    if (rderr /= 0) then
       print*, 'File not found: "'//trim(flnm)//'"'
       ! stop "READRIP: File not found."
       ierr = 1
       close(10)
       return
    endif

    read(10, iostat=rderr) vardesc, plchun, ihrip, rhrip, chrip
    if(rderr /= 0) stop "READRIP:  Read Problem 4."

    if (ihrip(6) == 2) then

       stop "READRIP:  readrip_3d trying to read 2d field"
       read(10, iostat=rderr) fld

    elseif (ihrip(6) == 3) then

       read(10,iostat=rderr) fld
       if(rderr /= 0) stop "READRIP_3D:  Read Problem 6."

    else
       print*, 'READRIP:  Not a 2d or 3d field:  "'//trim(fldname)//'"'
       stop
    endif

    close(10)
  end subroutine freadrip_3d

  subroutine readrip_3d(hour, fldname, casename, fld, ix, jx, kx, ierr)
    implicit none
    integer,                               intent(in)  :: hour
    character(len=*),                      intent(in)  :: fldname
    character(len=*),                      intent(in)  :: casename
    integer,                               intent(in)  :: ix, jx, kx
    real,             dimension(ix,jx,kx), intent(out) :: fld
    integer,                               intent(out) :: ierr

    character(len=256) :: flnm

    integer :: i
    integer :: rderr

    ierr = 0

    ! Build ripdp file names
    write(flnm, '(A,"_",I4.4,".00000_",A)') &
         casename, hour, fldname

    open(10, file=trim(flnm), status='old', form='unformatted', &
         action='read', iostat=rderr)
    if (rderr /= 0) then
       print*, 'File not found: "'//trim(flnm)//'"'
       ! stop "READRIP: File not found."
       ierr = 1
       return
    endif

    read(10, iostat=rderr) vardesc, plchun, ihrip, rhrip, chrip
    if(rderr /= 0) stop "READRIP:  Read Problem 4."

    if (ihrip(6) == 2) then

       stop "READRIP:  readrip_3d trying to read 2d field"
       read(10, iostat=rderr) fld

    elseif (ihrip(6) == 3) then

       read(10,iostat=rderr) fld
       if(rderr /= 0) stop "READRIP_3D:  Read Problem 6."

    else
       print*, 'READRIP:  Not a 2d or 3d field:  "'//trim(fldname)//'"'
       stop
    endif

    close(10)
  end subroutine readrip_3d

  subroutine ireadrip_header(hour, fldname, casename, &
       ihrip_return, rhrip_return, chrip_return)
    implicit none
    integer,                          intent(in)  :: hour
    character(len=*),                 intent(in)  :: fldname
    character(len=*),                 intent(in)  :: casename
    integer,           dimension(32), intent(out), optional :: ihrip_return
    real,              dimension(32), intent(out), optional :: rhrip_return
    character(len=64), dimension(64), intent(out), optional :: chrip_return

    character(len=256) :: flnm

    integer :: ierr

    real, dimension(:,:,:), allocatable :: dum


    ! Build ripdp file names
    write(flnm, '(A,"_",I4.4,".00000_",A)') &
         casename, hour, fldname

    open(10, file=trim(flnm), status='old', form='unformatted', &
         action='read', iostat=ierr)
    if (ierr /= 0) then
       print*, 'Looking for file:  ', trim(flnm)
       stop "File not found."
    endif

    read(10, iostat=ierr) vardesc, plchun, ihrip, rhrip, chrip
    if(ierr /= 0) stop "Read Problem 4."
    close(10)
    if (present(ihrip_return)) then
       ihrip_return = ihrip
    endif
    if (present(rhrip_return)) then
       rhrip_return = rhrip
    endif
    if (present(chrip_return)) then
       chrip_return = chrip
    endif

  end subroutine ireadrip_header

  subroutine freadrip_header(xhour, fldname, casename, &
       ihrip_return, rhrip_return, chrip_return)
    implicit none
    real,                             intent(in)  :: xhour
    character(len=*),                 intent(in)  :: fldname
    character(len=*),                 intent(in)  :: casename
    integer,           dimension(32), intent(out), optional :: ihrip_return
    real,              dimension(32), intent(out), optional :: rhrip_return
    character(len=64), dimension(64), intent(out), optional :: chrip_return

    character(len=256) :: flnm

    integer :: ierr, j

    real, dimension(:,:,:), allocatable :: dum

    ! Build ripdp file names
    write(flnm, '(A,"_",F10.5,"_",A)') casename, xhour, fldname
    ! Blanks may be created by the F10.5 specifier.  Fill replace those blanks with "0"
    do j = 1, len_trim(flnm)
       if (flnm(j:j) == " ") flnm(j:j) = "0"
    enddo

    open(10, file=trim(flnm), status='old', form='unformatted', &
         action='read', iostat=ierr)
    if (ierr /= 0) then
       print*, 'Looking for file:  ', trim(flnm)
       stop "File not found."
    endif

    read(10, iostat=ierr) vardesc, plchun, ihrip, rhrip, chrip
    if(ierr /= 0) stop "Read Problem 4."
    close(10)
    if (present(ihrip_return)) then
       ihrip_return = ihrip
    endif
    if (present(rhrip_return)) then
       rhrip_return = rhrip
    endif
    if (present(chrip_return)) then
       chrip_return = chrip
    endif

  end subroutine freadrip_header

  subroutine writerip_2d(hour, fldname, casename, fld, ix, jx, ierr)
    ! Uses the ihrip, rhrip, chrip already 
    ! in the module variables.  

    implicit none
    integer,                            intent(in)  :: hour
    character(len=*),                   intent(in)  :: fldname
    character(len=*),                   intent(in)  :: casename
    integer,                            intent(in)  :: ix, jx
    real,             dimension(ix,jx), intent(in)  :: fld
    integer,                            intent(out) :: ierr

    character(len=256) :: flnm
    integer :: rderr

    ierr = 0

    ! Build ripdp file names
    write(flnm, '(A,"_",I4.4,".00000_",A)') casename, hour, fldname

    print*, 'Writerip_2d:  Opening file: "', trim(flnm), '"'
    open(10, file=trim(flnm), status='new', form='unformatted', &
         action='write', iostat=rderr)
    if (rderr /= 0) then
       print*, 'Problem opening file as new file: "'//trim(flnm)//'"'
       ierr = 1
       stop
       return
    endif

!KWM    if (fldname == "STAGE4") then
!KWM       vardesc = "Stage 4 precip analysis"
!KWM       plchun  = "mm"
!KWM    else if (fldname == "QFX_AVG") then
!KWM       vardesc = "Average QFX for time of day"
!KWM    else
!KWM       print*, 'Fieldname not recognized:  "', fldname, '"'
!KWM       stop
!KWM    endif

    write(10) vardesc, plchun, ihrip, rhrip, chrip
    write(10) fld
    close(10)
    ! print*, 'Closing unit 10'

  end subroutine writerip_2d

  subroutine writerip_3d(hour, fldname, casename, fld, ix, jx, kx, ierr)
    ! Uses the ihrip, rhrip, chrip already 
    ! in the module variables.  

    implicit none
    integer,                               intent(in)  :: hour
    character(len=*),                      intent(in)  :: fldname
    character(len=*),                      intent(in)  :: casename
    integer,                               intent(in)  :: ix, jx, kx
    real,             dimension(ix,jx,kx), intent(in)  :: fld
    integer,                               intent(out) :: ierr

    character(len=256) :: flnm
    integer :: rderr

    ierr = 0

    ! Build ripdp file names
    write(flnm, '(A,"_",I4.4,".00000_",A)') casename, hour, fldname

    print*, 'Writerip_3d:  Opening file: "', trim(flnm), '"'
    open(10, file=trim(flnm), status='new', form='unformatted', &
         action='write', iostat=rderr)
    if (rderr /= 0) then
       print*, 'Problem opening file as new file: "'//trim(flnm)//'"'
       ierr = 1
       stop
       return
    endif

!KWM    if (fldname == "STAGE4") then
!KWM       vardesc = "Stage 4 precip analysis"
!KWM       plchun  = "mm"
!KWM    else if (fldname == "QFX_AVG") then
!KWM       vardesc = "Average QFX for time of day"
!KWM    else
!KWM       print*, 'Fieldname not recognized:  "', fldname, '"'
!KWM       stop
!KWM    endif

    write(10) vardesc, plchun, ihrip, rhrip, chrip
    write(10) fld
    close(10)
    ! print*, 'Closing unit 10'

  end subroutine writerip_3d

  subroutine read_xtimes(casename, xtimes, nx)
    implicit none
    character(len=*),                   intent(in)  :: casename
    integer,                            intent(in)  :: nx
    real,            dimension(nx),     intent(out) :: xtimes

    integer :: ntimes
    integer :: i
    integer :: ierr

    xtimes = -1.E36
    
    open(14, file=trim(casename)//".xtimes", status='old', form='formatted', action='read')
    read(14,*) ntimes
    
    do i = 1, ntimes
       read(14,*,iostat=ierr) xtimes(i)
       if (ierr /= 0) then
          write(*,'("Problem reading xtimes file ''",A,"''")') trim(casename)//".xtimes"
          stop "ERROR EXIT"
       endif
    enddo

    close(14)

  end subroutine read_xtimes


end module ripdp_module

