subroutine dyn(npart,mass,xyz,vxyz,ti,rmax,fpot,fsta)
use bspline_module
use interpolation
use RDistributions
implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Nico 24.05.2016 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! performs the dynamics until r=rf                        !!
! uses Verlet algorithom                                !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! for the dynamics
double precision, parameter :: dt = 5d0, dr = 0.05d0  ! in atomic units
double precision, intent(inout) :: ti, rmax
double precision, dimension(3,npart) :: xyzm, xyzt, xyznew, xyzdr
double precision :: time, grade, ecoll, ekin1, ekin2, egap, vunscal, vscal, scal
double precision, dimension(:), allocatable :: valdrp, valdrm, val

! for the system
integer, intent(in) :: npart
double precision, dimension(npart), intent(in) :: mass
double precision, dimension(3,npart), intent(inout) :: xyz, vxyz

!for the electronic state
character(64) :: fpot

!for Landau-Zener 
double precision, dimension(:), allocatable :: epair, epair_t1, d_epair, d_epair_t1, d2_epair
double precision :: plz

!for the interpolation
integer :: nsta
integer :: nr
double precision, dimension(:), allocatable :: r
double precision, dimension(:,:), allocatable :: energy

integer, intent(inout) :: fsta
double precision, dimension(:,:), allocatable :: tr

integer,parameter :: kr = 4     !! order in r
integer,parameter :: iknot = 0  !! automatically select the knots

double precision :: tol
logical :: fail, file_e

integer :: inbvx
integer :: iflag
integer :: idr

double precision :: newr, rt, oldr, ddr, ddrold, r0

integer :: i, j, k, ista
character(3) :: ur, ue

call init_random_seed()
!write(*,*)rand_uniform(0d0,1d0)
! Reads the data

inquire( file=trim(fpot), exist=file_e )
if ( file_e .eqv. .false. ) then
 write(*,*) trim(fpot), " does not exist"
 stop
endif

!write(*,'(A,A)')'Read PECs in ',trim(fpot)
open(unit=10,file=trim(fpot))
 read(10,*)ur, ue, nr, nsta
! write(*,'(A,I4,A)')'There are',nr,' points'
! write(*,'(A,I4,A)')'There are', nsta, ' states'

allocate(r(nr))
allocate(energy(nsta,nr))
!
do i = 1, nr
  read(10,*)r(i),(energy(ista,i),ista=1,nsta)
  if(ur=='ang' .or. ur=='Ang' .or. ur=='ANG' .or. ur=='A') then
   r(i) = r(i)*1.889725989d0
  endif
  if(ue=='eV' .or. ue=='ev' .or. ue=='EV') then
   write(*,*)ue,ur
   energy(:,i) = energy(:,i)/27.211385d0
  endif
enddo 
close(10)
!write(*,'(A)')'Reading PECs done '

allocate(epair(nsta),epair_t1(nsta),d_epair(nsta),d_epair_t1(nsta),d2_epair(nsta))
d_epair_t1(:) = 0d0
epair_t1(:) = 0d0

!
! Set-up the interpolation
inbvx = 1
idr = 0
fail = .false.
tol = 1.0e-14

allocate(tr(nsta,nr+kr))
allocate(val(nsta),valdrp(nsta),valdrm(nsta))

do ista = 1, nsta
 call db1ink(r,nr,energy(ista,:),kr,iknot,tr(ista,:),energy(ista,:),iflag)
 if(iflag/=0) then
  write(*,*)"error in db1ink"
  return
 endif
enddo

! landau-zenner surf. hopp. stuff
 call dist(npart,xyz,newr)
 rt = newr
 r0 = newr
do ista=1,nsta
  call db1val(newr,idr,tr(ista,:),nr,kr,energy(ista,:),val(ista),iflag,inbvx)
 if(iflag/=0) then
  write(*,*)"error in db1val at time",time, newr
  stop
 endif
enddo

do ista=1,nsta
  epair(ista) = abs(val(ista)-val(fsta))
enddo

! first step of the Verlet algorithm
do i = 1, npart
  do j = 1, 3
    xyzm(j,i) = xyz(j,i)
    xyzt(j,i) = xyz(j,i) + vxyz(j,i)*dt
  enddo
enddo

time=ti
do while(rt<rmax)

! landau-zenner surf. hopp. stuff
 epair_t1(:) = epair(:)
 d_epair_t1(:) = d_epair(:)
 epair(:)=0d0
 
! computes the energy at position at time t
 oldr=rt
 call dist(npart,xyzt,newr)
 ddrold=ddr
 ddr = (newr-oldr)/dt
 rt=newr
! write(*,'(8(f20.10,1X))')time,newr,ecoll*27.211d0,val(fsta)*27.211d0,ecoll*27.211d0+val(fsta)*27.211d0!0.5d0*(mass(1)*mass(2)/(mass(1)+mass(2)))*ddr**2,ddr
! write(200,'(4(f20.10,1X),i3)')time,newr!,ecoll*27.211d0,val(fsta)*27.211d0,fsta!ddr,ddrold,fsta
 do ista=1,nsta
  call db1val(newr,idr,tr(ista,:),nr,kr,energy(ista,:),val(ista),iflag,inbvx)

!  if(newr>rmax) then
!   write(*,*)fsta,newr
!   stop
!  endif

  if(iflag/=0) then
   write(*,*)"error in db1val at time",time
   stop
  endif
 enddo

! apply LZ surface hopping here
 do ista=1,nsta
  epair(ista) = abs(val(ista)-val(fsta))
!  write(100+ista,*)newr,epair(ista)
 enddo

 d_epair(:) = (epair(:)-epair_t1(:))/dt
 d2_epair(:) = (d_epair(:)-d_epair_t1(:))/dt

 do ista=1,nsta 
  if(d_epair(ista)*d_epair_t1(ista) < 0d0 .and. d2_epair(ista)>0d0 .and. ddr*ddrold>0d0) then  ! ddr*ddrold>0d0  <=> no hopping at the turning point
   plz = exp(-0.5d0*pi*sqrt(abs(epair(ista))**3/d2_epair(ista)))
!  write(*,'(200(f20.10,1X))')time,newr,epair(ista),plz,ddr,ddrold
!  write(*,'(200(f20.10,1X))')time,newr,plz,float(ista)!, epair(ista), epair_t1(ista), float(fsta), d2_epair(ista)
!  plz = 0d0

! is there enough kinetic energy to fill the energy gap, if not => frustrated hop
      ekin1 = 0d0
        do i = 1, npart
!nico yz          do j = 1, 3
         do j = 2, 3
             vunscal = (xyzt(j,i)-xyzm(j,i))/dt
             ekin1 = ekin1 + 0.5d0*mass(i)*vunscal**2
          enddo
        enddo

!      write(*,'(2(f20.15,1X),2i5)')newr, plz,ista,fsta
!      plz=0d0
      if(plz>rand_uniform(0d0,1d0) .and. epair(ista)<ekin1) then

! rescale the energy              
          egap = val(fsta)-val(ista)
!          write(*,*)"Egap,Ekin",egap*27.211d0,ekin1*27.211d0
          ekin1 = 0d0
          ekin2 = 0d0
          scal = egap/6d0 ! the velocities scaled evenly in the collision plane (i.e. energy gap is spread over all 6 coord. e.g.  xOz)
          do i = 1, npart
!nico yz          do j = 1, 3
         do j = 2, 3
             vunscal = (xyzt(j,i)-xyzm(j,i))/dt
!nico yz             if(vunscal/=0d0) then
               ekin1 = ekin1 + 0.5d0*mass(i)*vunscal**2
               if(1d0+2d0*scal/(mass(i)*vunscal**2)>0d0) then
                  vscal = vunscal*sqrt(1d0+2d0*scal/(mass(i)*vunscal**2))
                  xyzt(j,i) = xyzt(j,i) + (vscal-vunscal)*dt
               else !! in this case the evenly distributed correction is larger than the veloc. in this coord=> we put the velocity to zero (energy is not strictly conserved)
                  xyzt(j,i) = xyzm(j,i)
               endif
               vunscal = (xyzt(j,i)-xyzm(j,i))/dt
!               write(*,*)vunscal
               ekin2 = ekin2 + 0.5d0*mass(i)*vunscal**2
!nico yz             endif
           enddo
          enddo
!        write(*,*)"Ekin_unscal,Ekinscal,Energy conservation",ekin1*27.211d0, ekin2*27.211d0, (ekin2-egap-ekin1)*27.211d0
!        write(*,*)"HOP",time,newr,ista
        fsta=ista
        epair_t1(:) = 0d0
        epair(:) = 0d0
        d_epair(:) = 0d0
        d_epair_t1(:) = 0d0
       exit ! only one hop allowed
      endif
  endif
 enddo

!! end surface hopping

! computes the new position
 ecoll = 0d0
 do i = 1, npart
!nico yz          do j = 1, 3
         do j = 2, 3
  
    xyzdr(:,:) =  xyzt(:,:) 
    xyzdr(j,i) = xyzt(j,i) + dr
!    write(*,'(i4,i4,20(f20.15))')j,i,xyzdr(:,:)
    call dist(npart,xyzdr,newr)
 do ista=1,nsta
    call db1val(newr,idr,tr(ista,:),nr,kr,energy(ista,:),valdrp(ista),iflag,inbvx)
 enddo
!    write(*,'(i4,i4,2(f20.15))')j,i,newr,valdrp(fsta)

    xyzdr(:,:) =  xyzt(:,:) 
    xyzdr(j,i) = xyzt(j,i) - dr
    call dist(npart,xyzdr,newr)
 do ista=1,nsta
    call db1val(newr,idr,tr(ista,:),nr,kr,energy(ista,:),valdrm(ista),iflag,inbvx)
 enddo
!    write(*,'(i4,i4,2(f20.15))')j,i,newr,valdrm(fsta)

    grade = 0.5*( (valdrp(fsta)-val(fsta))/dr - (valdrm(fsta)-val(fsta))/dr )
!    write(*,'(i4,i4,2(f20.15))')j,i,grade,newr
    xyznew(j,i) = 2*xyzt(j,i) - xyzm(j,i) - grade*dt**2/mass(i) 
    vxyz(j,i) = (xyzt(j,i) - xyzm(j,i))/dt
    ecoll = ecoll + 0.5d0*mass(i)*((xyzt(j,i) - xyzm(j,i))/dt)**2

   enddo
 enddo

 xyzm(:,:) = xyzt(:,:)
 xyzt(:,:) = xyznew(:,:)
 time = time + dt

enddo
!   write(*,*)fsta,newr
!stop

xyz(:,:) = xyznew(:,:)
vxyz(:,:) = (xyzt(:,:)-xyzm(:,:))/dt

ekin1 = 0d0
do i = 1, npart
!nico yz          do j = 1, 3
  do j = 2, 3
     vunscal = (xyzt(j,i)-xyzm(j,i))/dt
     ekin1 = ekin1 + 0.5d0*mass(i)*vunscal**2
  enddo
enddo
write(*,'(i4,1X,100(f20.16,1X))')fsta,ekin1*27.211d0,newr,xyz(:,:),vxyz(:,:)

deallocate(r,energy)
deallocate(tr)
deallocate(val,valdrp,valdrm)
deallocate(epair,epair_t1,d_epair,d2_epair)

end subroutine dyn


subroutine dist(npart,xyz,r)
implicit none

integer, intent(in) :: npart
double precision, dimension(3,npart), intent(in) :: xyz

double precision :: r

r = sqrt( (xyz(1,2)-xyz(1,1))**2 + (xyz(2,2)-xyz(2,1))**2 + (xyz(3,2)-xyz(3,1))**2 )

end subroutine dist

