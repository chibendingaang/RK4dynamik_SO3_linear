program  main
use xdalphaparam

  implicit none

  integer ::  rank, nprocs,  ierr,i
  real(xd) ::  in1
  integer ::  seed
! real(xd) :: start, terminate
! integer :: spin_a, spin_b, spins_recvd
  character(len=1024) :: seedchar, inputseed
  real*8 :: begin, end


call getarg(1, inputseed)
read(inputseed, '(i8)') seed

write(seedchar, '(i8)') int(seed)
! single config, so whatever the seed only first digit of thousand order number (keeping the number small is arbitrary)

print*,"seed = ",seed
print*,"seedchar",seedchar
 

print*,"before dynamik"
call Dynamik(seed)
print*,"after dynamik"
   

print*,"rank=",rank

 open(177,file='./xpa2b0/L128/alpha_0pt8/2emin3/spin_a_'//trim(adjustl(seedchar))//'.dat',status='new')
 write(177, *) spin_a
 close(177)
 
 open(187,file='./xpa2b0/L128/alpha_0pt8/2emin3/spin_b_'//trim(adjustl(seedchar))//'.dat',status='new')
 write(187, *) spin_b
 close(187)

! print *, stoptime - starttime

end


!program for Heisenberg and or driven dynamik
subroutine Dynamik(diff)

    use xdalphaparam
    use mt19937xd
    use alphaspindynamics
    ! use totalHxd
    ! use sitemagsqxd

    Implicit None
    ! Integer, parameter :: xd = selected_real_kind(15,307)
    character flindx*20 
    CHARACTER*15:: Buf,foldernm, seedinput
    integer :: ti
    integer :: x
    !integer :: conf
    real(xd) :: epsil 
    integer :: diff
    real(xd):: sz_cos, sz_sin
    real(xd), allocatable :: dummyzero(:)


    allocate(dummyzero(3*L))
    allocate(planar(L)); allocate(phi(L)); allocate(xyplanar(L))

    allocate(k1(3*L)); allocate(k2(3*L)); allocate(k3(3*L)); allocate(k4(3*L))
    allocate(s_a(3*L)); allocate(s_b(3*L)); allocate(s0_a(3*L)); allocate(s0_b(3*L))


    allocate(spin_a(3*L, 0:steps/100)); allocate(spin_b(3*L, 0:steps/100))
    ! due to Fortran's column-based storage, use transpose of the intended array
    ! upon reshaping it in fortran: this will be of the dimensions: (steps/100, L, 3)

    k1 = 0.0_xd; k2 = 0.0_xd; k3 = 0.0_xd; k4 = 0.0_xd
    s0_a = 0.0_xd; s0_b = 0.0_xd; s_a = 0.0_xd; s_b = 0.0_xd

    conf=diff
    write(foldernm,'(i8)')conf !harsh

    ! READ(seedinput,*)seedz

    ! print all the input params
    print *, 'Pi = ',  Pi
    print *, 'L = ', L
    print *, 'steps = ', steps
    print *, 'dt = ', dt
    print *, 'lambda = ', lambda
    print *, 'mu = ', mu
    print *, 'a = dt*(Lambda+Mu) = ', a
    print *, 'b = dt*(Lambda-Mu) = ', b

    ! Initialization
    ! keeping epsil as a real16 variable instead of storing as a parameter, to avoid the error: Type mismatch in argument ‘eps’ at (1); passed REAL(8) to REAL(16)

    ! epsil is suppressed for InitPerturb funtion, but needed to call InitRandom
    epsil = 0.0_xd
    call InitRandom(s0_a, epsil)
    s_a = s0_a
    print *, 'epsilon_a = ', epsil
    deallocate(s0_a)

    ti=0; !write(flindx,'(i8)') ti
    spin_a(:,ti) = s_a

    do ti = 1,steps
    ! write(flindx,'(i8)') ti     
        Call RK4dynamik(s_a)
        if (mod(ti,100)==0) then
        spin_a(:,ti/100) = s_a
        end if

    end do

    epsil = 0.001_xd
    call InitRandom(s0_b, epsil)
    s_b = s0_b
    print *, 'epsilon_b = ', epsil
    deallocate(s0_b)

    ti=0;
    spin_b(:,ti) = s_b

    do ti = 1,steps
        Call RK4dynamik(s_b)
        if (mod(ti,100)==0) then
            spin_b(:,ti/100) = s_b
        end if

    end do

End subroutine Dynamik
