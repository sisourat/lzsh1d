program test_interp3d
use bspline_module
implicit none

integer :: nr1, nr2, nthe
double precision, dimension(:), allocatable :: r1, r2, the
double precision, dimension(:,:,:), allocatable :: energy

double precision, dimension(:), allocatable :: tr1, tr2, tthe

integer,parameter :: kr1 = 4     !! order in r1
integer,parameter :: kr2 = 4     !! order in r2
integer,parameter :: kt = 4     !! order in theta
integer,parameter :: iknot = 0  !! automatically select the knots

double precision :: tol
logical :: fail

integer :: inbvx,inbvy,inbvz
integer :: iloy,iloz,iflag
integer :: idr1, idr2, idthe

double precision :: newr1, newr2, newthe, val

integer :: i, j ,k

! Reads the data

read(*,*)nr1, nr2, nthe

allocate(r1(nr1),r2(nr2),the(nthe))
allocate(energy(nr1,nr2,nthe))

do i = 1, nr1
  do j = 1, nr2
    do k = 1, nthe
      read(*,*)r1(i),r2(j),the(k),energy(i,j,k)
!      write(*,*)r1(i),r2(j),the(k),energy(i,j,k)
    enddo
  enddo
enddo 

! Set-up the interpolation
inbvx = 1
inbvy = 1
inbvz = 1
iloy  = 1
iloz  = 1

idr1 = 0
idr2 = 0
idthe = 0

fail = .false.
tol = 1.0e-14

allocate(tr1(nr1+kr1),tr2(nr2+kr2),tthe(nthe+kt))

call db3ink(r1,nr1,r2,nr2,the,nthe,energy,kr1,kr2,kt,iknot,tr1,tr2,tthe,energy,iflag)

! Performs the interpolation

do i = 1, 580
  write(*,*)
  do j = 1, 580
    do k = 1, 1

      newr1 = 1.6 + (i-1)*0.005
      newr2 = 1.6 + (j-1)*0.005
      newthe = 103.0 !+ (k-1) 

      call db3val(newr1,newr2,newthe,idr1,idr2,idthe,tr1,tr2,tthe,nr1,nr2,nthe,kr1,kr2,kt,energy,val,iflag,inbvx,inbvy,inbvz,iloy,iloz)
      write(*,*) newr1,newr2,val

    enddo
  enddo
enddo


deallocate(r1,r2,the,energy)
deallocate(tr1,tr2,tthe)

end program test_interp3d
