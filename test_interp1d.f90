program test_interp1d
use bspline_module
implicit none

integer :: nr
double precision, dimension(:), allocatable :: r
double precision, dimension(:), allocatable :: energy

double precision, dimension(:), allocatable :: tr

integer,parameter :: kr = 4     !! order in r
integer,parameter :: iknot = 0  !! automatically select the knots

double precision :: tol
logical :: fail

integer :: inbvx
integer :: iflag
integer :: idr

double precision :: newr, val

integer :: i, j ,k

! Reads the data

read(*,*)nr

allocate(r(nr))
allocate(energy(nr))

do i = 1, nr
      read(*,*)r(i),energy(i)
enddo 

! Set-up the interpolation
inbvx = 1
idr = 0

fail = .false.
tol = 1.0e-14

allocate(tr(nr+kr))

call db1ink(r,nr,energy,kr,iknot,tr,energy,iflag)

! Performs the interpolation

do i = 1, 580
   newr = 1.6 + (i-1)*0.005
   call db1val(newr,idr,tr,nr,kr,energy,val,iflag,inbvx)
   write(*,*) newr,val
enddo


deallocate(r,energy)
deallocate(tr)

end program test_interp1d
