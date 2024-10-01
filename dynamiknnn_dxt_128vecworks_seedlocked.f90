program  main
    use nnnparam_128vec
    
      implicit none
      real :: start_time, stop_time 
      integer ::  rank, nprocs,  ierr,i 
      real(real128) ::  in1
      integer ::  seed, serial_number, io_status, line_number
    ! real(real128) :: start, terminate
    ! integer :: spin_a, spin_b, spins_recvd
      character(len=1024) :: seedkey, input_serial !, command, seed_lookup
      real*8 :: begin, end
    
    
    call cpu_time(start_time) 
    call getarg(1, input_serial) ! from 1 to 50,000 ideally (capacity of 1,000,000,000 for MT19937)
    read(input_serial, '(i8)') serial_number
    
    ! Construct and execute the command to run the Python script and capture the output seed
    !WRITE(command, '(A,I5)') "python3 lookuptable_script.py ", serial_number
    ! creating a temp seed_output.txt for contextual/input serial numbers may cause confusion later on
    ! better to rely on a lookup_table with all the serial numbers in it
    !CALL EXECUTE_COMMAND_LINE(command // " >> seed_output.txt")
    ! Open the output file
    !OPEN(UNIT=10, FILE="seed_output.txt", STATUS="OLD", ACTION="READ")

    ! Open the lookup table file
    OPEN(UNIT=10, FILE="lookup_table.txt", STATUS="OLD", ACTION="READ")

    ! Loop through the file to find the matching serial number
    DO
      ! Read the serial number and seed from the file
      READ(10, '(I5, 1X, I10)', IOSTAT=io_status) line_number, seed
      IF (io_status /= 0) THEN
        PRINT*, "Error: Serial number not found or end of file reached."
        EXIT
      END IF

      ! If the serial number matches, print the seed and exit
      IF (line_number == serial_number) THEN
        PRINT*, "The seed for serial number (seedkey) ", serial_number, " is ", seed
        EXIT
      END IF
    END DO

    ! Close the file
    CLOSE(10)

 
    
    write(seedkey, '(i8)') int(serial_number)
    ! single config, so whatever the seed only first digit of thousand order number (keeping the number small is arbitrary)
    
    !print*,"seed, seedkey = ",seed, seedkey
    
    
    print*,"before dynamik"
    call Dynamik(seed)
    print*,"after dynamik"
       
    
    ! print*,"rank=",rank
    
    ! file_a
    !open(177,file='./qpdrvn_NNN/L2048/1emin3/J2byJ1_0pt625/spin_a_'//trim(adjustl(seedkey))//'.dat',status='new')
    !write(177, *) spin_a
    !close(177)

    
    ! file_b
    !open(187,file='./qpdrvn_NNN/L2048/1emin3/J2byJ1_0pt625/spin_b_'//trim(adjustl(seedkey))//'.dat',status='new')
    !write(187, *) spin_b
    !close(187)
    
      ! Save Dxt
       open(197, file='./qpdrvn_NNN/L128/1emin3/J2byJ1_0pt250/Dxt_'//trim(adjustl(seedkey))//'.dat', status='new')
       write(197, *) Dxt
       close(197) 


    call cpu_time(stop_time)
    print '("Program executed in ",f10.3," seconds.")', stop_time - start_time 
    
    end
        
    !program for Heisenberg and or driven dynamik
    subroutine Dynamik(diff)
    
        use nnnparam_128vec
        use mt19937_128
        use nnnspindynamik_128vec
        ! use totalH_128
        ! use sitemagsq_128
    
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
        allocate(Dxt(L,0:steps/100))
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
      
        epsil = 0.001_real128
        print *, 'epsilon_b = ', epsil
        call InitRandom(s0_b, epsil)
        s_b = s0_b
        deallocate(s0_b)
        
        ti=0; !write(flindx,'(i8)') ti
        spin_a(:,:,ti) = s_a
        spin_b(:,:,ti) = s_b
        
        do ti = 1,steps
        ! write(flindx,'(i8)') ti     
            Call RK4dynamik_nnn(s_a)
            Call RK4dynamik_nnn(s_b)
            
            if (mod(ti,100)==0) then
                spin_a(:,:,ti/100) = s_a
                spin_b(:,:,ti/100) = s_b
                
                ! Compute Dxt for the current time step
                do x = 1, L
                    !dot_product = sum(spin_a(:,x,ti/100) * spin_b(:,x,ti/100))
                    Dxt(x, ti/100) = 1.0_real128 - dot_product(s_a(:,x), s_b(:,x))
                end do
                
            end if
        end do
    End subroutine Dynamik
    
