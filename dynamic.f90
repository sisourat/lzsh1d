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
 open(unit=10,file="traj")
  read(10,*)fpot, rmax, fsta
 do i = 1, npart
  read(10,*)mass(i),xyz(1,i),xyz(2,i),xyz(3,i),vxyz(1,i),vxyz(2,i),vxyz(3,i)
!  write(100,'(6(f20.15,1X))')xyz(1,i),xyz(2,i),xyz(3,i),vxyz(1,i),vxyz(2,i),vxyz(3,i)
  mass(i)=mass(i)*1836.15d0
 enddo
 close(10)

  if(vxyz(1,1) /= 0d0 .or. vxyz(2,1) /= 0d0 .or. vxyz(3,1) /= 0d0) then
    write(*,*)"First atom should be fixed"
    stop
  endif
  if(vxyz(1,2) /= 0d0 .or. vxyz(2,2) /= 0d0) then
    write(*,*)"Velocity of the second atom should be along z"
    stop
  endif

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

