program  main
    use nnnparam_vec
    
      implicit none
      real :: start_time, stop_time 
      integer ::  rank, nprocs,  ierr,i
      real(real128) ::  in1
      integer ::  seed
    ! real(real128) :: start, terminate
    ! integer :: spin_a, spin_b, spins_recvd
      character(len=1024) :: seedchar, inputseed, command1, command2, command1pt0, command2pt0
      real*8 :: begin, end
    
    
    call cpu_time(start_time) 
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
    
     open(177,file='./xpdrvn_NNN/L512/2emin3/J1only/spin_a_'//trim(adjustl(seedchar))//'.dat',status='new')
     write(177, *) spin_a
     close(177)

    write(command1, '(A)') 'python save_dynamikvec_to_npz.py spin_a_'//trim(adjustl(seedchar))//''
    call execute_command_line(command1, wait=.true.)

    write(command1pt0, '(A)') 'rm ./xpdrvn_NNN/L512/2emin3/J1only/spin_a_'//trim(adjustl(seedchar))//'.dat'
    call execute_command_line(command1pt0, wait=.true.)

    open(187,file='./xpdrvn_NNN/L512/2emin3/J1only/spin_b_'//trim(adjustl(seedchar))//'.dat',status='new')
    write(187, *) spin_b
    close(187)


    write(command2, '(A)') 'python save_dynamikvec_to_npz.py spin_b_'//trim(adjustl(seedchar))//''
    call execute_command_line(command2, wait=.true.)

    write(command2pt0, '(A)') 'rm ./xpdrvn_NNN/L512/2emin3/J1only/spin_b_'//trim(adjustl(seedchar))//'.dat'
    call execute_command_line(command2pt0, wait=.true.)
 
    call cpu_time(stop_time)
    print '("Program executed in ",f10.3," seconds.")', stop_time - start_time 
    
    end
    
    
    !program for Heisenberg and or driven dynamik
    subroutine Dynamik(diff)
    
        use nnnparam_vec
        use mt19937_128
        use nnnspindynamik_vec
        ! use totalHxd
        ! use sitemagsqxd
    
        Implicit None
        character flindx*20 
        CHARACTER*15:: Buf,foldernm, seedinput
        integer :: ti
        integer :: x
        !integer :: conf
        real(real128) :: epsil 
        integer :: diff
        real(real128):: sz_cos, sz_sin
        real(real128), allocatable :: dummyzero(:,:)
    
    
        allocate(dummyzero(3,L))
        allocate(planar(L)); allocate(phi(L)); allocate(xyplanar(L))
    
        allocate(k1(3,L)); allocate(k2(3,L)); allocate(k3(3,L)); allocate(k4(3,L))
        allocate(s_a(3,L)); allocate(s_b(3,L)); allocate(s0_a(3,L)); allocate(s0_b(3,L))
    
    
        allocate(spin_a(3,L, 0:steps/100)); 
        allocate(spin_b(3,L, 0:steps/100))
        ! due to Fortran's column-based storage, use transpose of the intended array
        ! upon reshaping it in fortran: this will be of the dimensions: (steps/100, L, 3)
    
        k1 = 0.0_real128; k2 = 0.0_real128; k3 = 0.0_real128; k4 = 0.0_real128
        s0_a = 0.0_real128; s0_b = 0.0_real128; s_a = 0.0_real128; s_b = 0.0_real128
    
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
        !print *, 'a = dt*(Lambda+Mu) = ', a
        !print *, 'b = dt*(Lambda-Mu) = ', b
    
        ! Initialization
        ! keeping epsil as a real16 variable instead of storing as a parameter, to avoid the error: Type mismatch in argument ‘eps’ at (1); passed REAL(8) to REAL(16)
    
        ! epsil is suppressed for InitPerturb funtion, but needed to call InitRandom
        
        epsil = 0.0_real128
        print *, 'epsilon_a = ', epsil
        call InitRandom(s0_a, epsil)
        s_a = s0_a
        deallocate(s0_a)
    
        ti=0; !write(flindx,'(i8)') ti
        spin_a(:,:,ti) = s_a
        do ti = 1,steps
        ! write(flindx,'(i8)') ti     
            Call RK4dynamik_nnn(s_a)
            if (mod(ti,100)==0) then
            spin_a(:,:,ti/100) = s_a
            end if
        end do
    
    
        epsil = 0.001_real128
        print *, 'epsilon_b = ', epsil
        call InitRandom(s0_b, epsil)
        s_b = s0_b
        deallocate(s0_b)
        
        ti=0;
        spin_b(:,:,ti) = s_b
        do ti = 1,steps
            Call RK4dynamik_nnn(s_b)
            if (mod(ti,100)==0) then
                spin_b(:,:,ti/100) = s_b
            end if
        end do
    
    End subroutine Dynamik
    
