
program  main
use quad3param

!*****************************************************************************80
!
!! P1_RECEIVE_INPUT receives input from process 0.
!  Parameters:
!
!    Output, integer ( kind = 4 ) ::  INPUT1, the value of the parameter.
!
  use mpi

  implicit none

  integer ::  rank, nprocs,  ierr,i
  real(qw) ::  in1
  integer, dimension(MPI_STATUS_SIZE) :: status1
! integer  ::  status(MPI_STATUS_SIZE) !source variable in subroutine 3/9 p0_receiveoutput
  integer ::  tag, seed, dest, ki
  integer :: seeds_recvd, seed_in
  !real(qw) :: start, terminate
!  integer :: spin_a, spin_b, spins_recvd
  character(len=1024) :: seedchar, inputseed
  real*8 :: begin, end
 Call MPI_Init(ierr)
 
 call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
 call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)

!goto 13
call getarg(1, inputseed)
read(inputseed, '(i8)') ki 

seed=ki+rank

 write(seedchar, '(i8)') seed

!if(rank==3)
print*,"seed=",seed
print*,"seedchar",seedchar
 
print*,"before dynamik"

! seeds_recvd = seeds_recvd + 1
call Dynamik(seed)!(seed_in, spin_a, spin_b)

print*,"after dynamik"
!goto 13

begin = MPI_Wtime()

do i=0,nprocs-1
  if(rank==i)then
!   start = MPI_Wtime()
   open(176,file='./qwa2b0/2emin3/spin_a_'//trim(adjustl(seedchar))//'.dat',status='new')
   write(176,*)spin_a
   close(176)

!goto 16
   open(189,file='./qwa2b0/2emin3/spin_b_'//trim(adjustl(seedchar))//'.dat',status='new')
   write(189, *) spin_b
   close(189)
16 continue

   print*,"rank=",rank
!   terminate = MPI_Wtime()
!   print *, terminate - start 
  endif
enddo


13 continue

end = MPI_Wtime()
print *, end - begin
 Call MPI_Finalize(ierr)

end


!program for Heisenberg and or driven dynamik

!define module

subroutine Dynamik(diff)

use quad3param
use quadmt19937
use quad3spindynamics
!use totalHxd
!use sitemagsqxd

Implicit None
!Integer, parameter :: xd = selected_real_kind(15,307)
character flindx*20			!harsh
CHARACTER*15:: Buf,foldernm, seedinput
integer :: ti
integer :: x
!integer :: conf
real(qw) :: epsil 
integer :: diff
real(qw):: sz_cos, sz_sin
real(qw), allocatable :: dummyzero(:)


allocate(dummyzero(3*L))
allocate(planar(L)); allocate(phi(L)); allocate(xyplanar(L))

allocate(k1(3*L)); allocate(k2(3*L)); allocate(k3(3*L)); allocate(k4(3*L))
allocate(s_a(3*L)); allocate(s_b(3*L)); allocate(s0_a(3*L)); allocate(s0_b(3*L))

!goto 99
allocate(spin_a(3*L, 0:steps/100)); allocate(spin_b(3*L, 0:steps/100))
!due to Fortran's column-based storage, use transpose of the intended array

k1 = 0.0_qw; k2 = 0.0_qw; k3 = 0.0_qw; k4 = 0.0_qw
s0_a = 0.0_qw; s0_b = 0.0_qw; s_a = 0.0_qw; s_b = 0.0_qw


!CALL GETARG(1,Buf) 
!CALL Getarg(2, seedinput)
conf=diff
!READ(Buf,*)conf
write(foldernm,'(i8)')conf		!harsh

!READ(seedinput,*)seedz

!print all the input params
print *, 'Pi = ',  Pi
print *, 'L = ', L
print *, 'steps = ', steps
print *, 'dt = ', dt
print *, 'lambda = ', lambda
print *, 'mu = ', mu
print *, 'a = dt*(Lambda+Mu) = ', a
print *, 'b = dt*(Lambda-Mu) = ', b
!print *, 'seed input = ', seedz

!open(72, file = 'config_seedxdd2048.txt', status = 'old', position = 'append')
!write(72,* ) conf, char(9),  seedz
! close(72)



!Initialization
!keeping epsil as a real16 variable instead of storing as a parameter, to avoid the error: Type mismatch in argument ‘eps’ at (1); passed REAL(8) to REAL(16)
epsil = 0.0_qw
print *, 'epsilon_a = ', epsil
call InitRandom(s0_a, epsil)


!s0_a(2) = s0_a(2)*(1 + epsil*sz_cos/sz_sin)
!s0_a(1) = s0_a(1)*(1 + epsil*sz_cos/sz_sin)

s_a = s0_a

ti=0; !write(flindx,'(i8)') ti
spin_a(:,ti) = s_a
! open(1, file = './spin_a_'//trim(adjustl(foldernm))//'.dat', status = 'new')
 !open(61, file = './L2048/xddrvn/2emin3/sitemagsqa_'//trim(adjustl(foldernm))//'.dat', status = 'new')
 !open(51, file = './L2048/xddrvn/2emin3/totalHa_'//trim(adjustl(foldernm))//'.dat', status = 'new')
 !do x = 1,3*L
! write(1,*) s_a !(x) 
 !end do 
!  call totalenergy(s_a, totalH)
 ! call sitemagzsq(s_a, sitemagsq)
 ! write(51,*) totalH
  !write (61, *) sitemagsq
! close(1)
 !close(61)
 !close(51)
!Performing the numerical dynamics

!open(1, file = './spin_a_'//trim(adjustl(foldernm))//'.dat', status = 'old', position = 'append')
!open(61, file = './L2048/xddrvn/2emin3/sitemagsqa_'//trim(adjustl(foldernm))//'.dat', status = 'old', position = 'append')
!open(51, file = './L2048/xddrvn/2emin3/totalHa_'//trim(adjustl(foldernm))//'.dat', !status = 'old', position = 'append')
do ti = 1,steps
	!write(flindx,'(i8)') ti     
    Call RK4dynamik(s_a)
    !Call sitemagzsq(s_a, sitemagsq)
    !Call totalenergy(s_a, totalH)
    !write(51,*) totalH
    !close(51)

    	if (mod(ti,100)==0) then
        spin_a(:,ti/100) = s_a
	!do x  = 1,3*L
!	write(1,*) s_a !(x) 
        !write(61,*) sitemagsq
	!end do 
	end if
end do
! close(1)
! close(61)
!spina = reshape(spin_a, shape(s_a)*(steps/100 + 1))

!clear(epsil)

epsil = 0.001_qw
print *, 'epsilon_b = ', epsil
call InitRandom(s0_b, epsil)

!s0_b = s0_a

s_b = s0_b
ti=0; !write(flindx,'(i8)') ti
! open(2, file = './spin_b_'//trim(adjustl(foldernm))//'.dat', status = 'new')
! open(62, file = './L2048/xddrvn/2emin3/sitemagsqb_'//trim(adjustl(foldernm))//'.dat', status = 'new')
 
spin_b(:,ti) = s_b
! write(2,*) s_b
  !call sitemagzsq(s_b, sitemagsq)
  ! write(62, *) sitemagsq
! close(2)
! close(62)

!must reset the RK4 subroutine to get spinb correct
dummyzero = 0.0_qw
Call RK4dynamik(dummyzero)

!Same dynamics for the initially perturbed system
!open(2, file = './spin_b_'//trim(adjustl(foldernm))//'.dat', status = 'old', position = 'append')
!open(62, file = './L2048/xddrvn/2emin3/sitemagsqb_'//trim(adjustl(foldernm))//'.dat', status = 'old',position = 'append')

do ti = 1,steps
	!write(flindx,'(i8)') ti
   Call RK4dynamik(s_b)
   !call sitemagzsq(s_b, sitemagsq)
   !   call totalenergy(s_b, totalH)

   	if (mod(ti,100)==0) then
        spin_b(:,ti/100) = s_b
	!do x = 1,3*L
!	write(2,*) s_b !(x)
    !write(62,*) sitemagsq
        !end do
	end if
	end do
! close(2)
 !close(62)
!spinb = reshape(spin_b, shape(s_b)*(steps/100 +1))

End subroutine Dynamik




