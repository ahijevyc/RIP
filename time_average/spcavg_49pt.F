program timeavg
  implicit none

!
!   RIP header variables
!
  integer, dimension(32) :: ihrip
  real, dimension(32) :: rhrip
  real, dimension(128) :: fullsigma, halfsigma
  character(len=64), dimension(64) :: chrip
  character(len=64) :: vardesc
  character(len=24) :: plchun

  real, allocatable, dimension(:,:,:) :: arr
  real, allocatable, dimension(:,:,:) :: avg

  character(len=256) :: ripname, flnm
  integer :: iunit=15, ounit=16, ierr

  integer :: i, j, k, ix, jx, kx, filecount, numarg
  integer, external :: iargc

  numarg = iargc()
  do filecount = 1, numarg
     call getarg(filecount, ripname)

     ! File
     open(iunit, file=(trim(ripname)), status='old', form='unformatted', &
          action='read', iostat=ierr)
     read (iunit) vardesc,plchun,ihrip,rhrip,chrip

     ix = ihrip(4)
     jx = ihrip(5)
     if (ihrip(6) == 3) then
        kx = ihrip(9)
     else if ( ihrip(6) == 2) then
        kx = 1
     else
        print*, 'ndim = ', ihrip(6)
        print*,' This is a problem.'
        stop
     endif
     allocate(arr(ix, jx, kx))
     allocate(avg(ix, jx, kx))
     print*, 'shape arr = ', shape(arr)

     read (iunit, iostat=ierr) arr
     print*, 'ierr = ', ierr
     close(iunit)

     do k = 1, kx
        do i = 1, ix
           do j = 1, jx
              avg(i,j,k) = arr(i,j,k)
           enddo
        enddo
        do i = 4, ix-3
           do j = 4, jx-3
              avg(i,j,k) = (arr(i,j,k)+arr(i+1,j,k)+arr(i-1,j,k)+   &
                 arr(i,j+1,k)+arr(i,j-1,k)+arr(i+1,j+1,k)+          &
                 arr(i+1,j-1,k)+arr(i-1,j-1,k)+arr(i-1,j+1,k)+      &
                 arr(i,j+2,k)+arr(i+1,j+2,k)+arr(i+2,j+2,k)+        &
                 arr(i-2,j+2,k)+arr(i-1,j+2,k)+arr(i,j-2,k)+        &
                 arr(i+1,j-2,k)+arr(i+2,j-2,k)+arr(i-2,j-2,k)+      &
                 arr(i-1,j-2,k)+arr(i-2,j-1,k)+arr(i-2,j,k)+        &
                 arr(i-2,j+1,k)+arr(i+2,j-1,k)+arr(i+2,j,k)+        &
                 arr(i+2,j+1,k)+arr(i,j+3,k)+arr(i+1,j+3,k)+        &
                 arr(i+2,j+3,k)+arr(i+3,j+3,k)+arr(i-3,j+3,k)+      &
                 arr(i-2,j+3,k)+arr(i-1,j+3,k)+arr(i,j-3,k)+        &
                 arr(i+1,j-3,k)+arr(i+2,j-3,k)+arr(i+3,j-3,k)+      &
                 arr(i-3,j-3,k)+arr(i-3,j-2,k)+arr(i-3,j-1,k)+      &
                 arr(i-3,j-2,k)+arr(i-3,j-1,k)+arr(i-3,j,k)+        &
                 arr(i-3,j+1,k)+arr(i-3,j+2,k)+arr(i+3,j-2,k)+      &
                 arr(i+3,j-1,k)+arr(i+3,j,k)+arr(i+3,j+1,k)+        &
                 arr(i+3,j+2,k))/49.
           enddo
        enddo 
     enddo
        

     flnm = trim(ripname)//"s"
     open(ounit, file=(trim(flnm)), form='unformatted', action='write')
     write (ounit) &
          vardesc,plchun,ihrip,rhrip,chrip
     write (ounit) avg
     close(ounit)

     deallocate(arr, avg)
  enddo
        
end program timeavg
