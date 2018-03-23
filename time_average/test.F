program test
   implicit none 
   character(len=256) :: h
   call getarg(1, h)
   print*, 'h = ', trim(h)
end program test
