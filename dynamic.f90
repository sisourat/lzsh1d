program dynamic_1D
implicit none

integer, parameter :: npart = 2
double precision, dimension(3,npart) :: xyz, vxyz, xyzt, xyzm, xyznew
double precision, dimension(npart) :: mass

double precision :: r, rmax, ecoll
double precision :: ti

integer :: fsta
character(64) :: fpot

integer :: i, j, l

! reads the initial conditions
  read(*,*)fpot, rmax, fsta
do i = 1, npart
  read(*,*)mass(i),xyz(1,i),xyz(2,i),xyz(3,i),vxyz(1,i),vxyz(2,i),vxyz(3,i)
!  write(100,'(6(f20.15,1X))')xyz(1,i),xyz(2,i),xyz(3,i),vxyz(1,i),vxyz(2,i),vxyz(3,i)
  mass(i)=mass(i)*1836.15d0
enddo

  if(vxyz(1,1) /= 0d0 .or. vxyz(2,1) /= 0d0 .or. vxyz(3,1) /= 0d0) then
    write(*,*)"First atom should be fixed"
    stop
  endif
  if(vxyz(1,2) /= 0d0 .or. vxyz(2,2) /= 0d0) then
    write(*,*)"Velocity of the second atom should be along z"
    stop
  endif

! add up 1/R0 to last particule along z
 ecoll = 0.5d0*mass(2)*vxyz(3,2)**2 
! write(*,*)rmax,ecoll*27.211d0
 ecoll = ecoll + 1d0/rmax
! write(*,*)rmax,ecoll*27.211d0
 vxyz(3,2) = sqrt(2d0*ecoll/mass(2))
 ecoll = 0.5d0*mass(2)*vxyz(3,2)**2 
! write(*,*)rmax,ecoll*27.211d0

! starts with the dynamics

 ti = 0d0
 call dyn(npart,mass,xyz,vxyz,ti,rmax,fpot,fsta)

 call dist(npart,xyz,r)
! write(*,*)r
! write(*,*)
! do i = 1, npart
!   write(100,'(6(f20.15,1X))')xyz(:,i),vxyz(:,i)
! enddo
!  write(100,*)

end program dynamic_1D

