program  dynamiknnn_mpi_128vectorized
    use nnnparam_128vec
    use mpi
    
      implicit none
      
      real :: start_time, stop_time 
      integer ::  rank, nprocs,  ierr,i 
      integer, dimension(MPI_STATUS_SIZE) :: status1
      integer ::  tag, dest, ki
      integer :: seeds_recvd, serialkey_
      character(len=1024) :: seedchar, seed_index, command1, command2,command1pt0, command2pt0
      real ::  in1
      real*8 :: begin, end


     
    

    !call cpu_time(start_time) 
    call MPI_Init(ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierr)
    
    call getarg(1, seed_index)
    read(seed_index, '(i8)') ki 
    
    serialkey_ = ki+rank
    write(seedchar, '(i8)') int(serialkey_)
    ! single config, so whatever the seed only first digit of thousand order number (keeping the number small is arbitrary)
    
    print*,"serial key = ",serialkey_
    print*,"seedchar",seedchar
     
    
    print*,"before dynamik"
    call Dynamik(serialkey_)
    print*,"after dynamik"
       
    start_time = MPI_Wtime()

    do i=0,nprocs-1
      if(rank==i)then
        print*,"rank=",rank
    
    ! file_a
    !   open(177,file='./qpdrvn_NNN/L2048/1emin3/J2byJ1_0pt625/spin_a_'//trim(adjustl(seedchar))//'.dat',status='new')
    !    write(177, *) spin_a
    !   close(177)

    
        ! file_b
    !    open(187,file='./qpdrvn_NNN/L2048/1emin3/J2byJ1_0pt625/spin_b_'//trim(adjustl(seedchar))//'.dat',status='new')
    !    write(187, *) spin_b
    !    close(187)
    
    ! Save Dxt
        open(197, file='./qpdrvn_NNN/L512/1emin3/J1only/Dxt_'//trim(adjustl(seedchar))//'.dat', status='new')
        write(197, *) Dxt
        close(197)
      endif
    enddo


    stop_time = MPI_Wtime()
    !call cpu_time(stop_time)
    print '("Program executed in ",f10.3," seconds.")', stop_time - start_time 
     
    Call MPI_Finalize(ierr)
    
    end
    
    
    !program for Heisenberg and or driven dynamik
    subroutine Dynamik(initseed)
    
        use nnnparam_128vec
        use mt19937_128
        use nnnspindynamik_128vec
        ! use totalH_64
        ! use sitemagsq_64
    
        Implicit None
        character flindx*20 
        CHARACTER*15:: Buf,foldernm, seedinput
        integer :: ti
        integer :: x
        !integer :: conf
        real(real128) :: epsil 
        integer :: initseed, line_number, seed_val, io_status
        real(real128):: sz_cos, sz_sin
        real(real128), allocatable :: dummyzero(:,:)
        
        !real(real128) :: dot_product

      ! Open the lookup table file
      OPEN(UNIT=10, FILE="lookup_table.txt", STATUS="OLD", ACTION="READ")

      ! Loop through the file to find the matching serial number
      DO
        ! Read the serial number and seed from the file
        READ(10, '(I5, 1X, I10)', IOSTAT=io_status) line_number, seed_val
        IF (io_status /= 0) THEN
          PRINT*, "Error: Serial number not found or end of file reached."
          EXIT
        END IF

        ! If the serial number matches, print the seed and exit
        IF (line_number == initseed) THEN
          PRINT*, "The seed for serial number (serialkey_) ", initseed, " is ", seed_val
          EXIT
        END IF
      END DO

      ! Close the file
      CLOSE(10)

      write(seedinput, '(i8)') int(initseed)

   
        allocate(dummyzero(3,L))
        allocate(planar(L)); allocate(phi(L)); allocate(xyplanar(L))
    
        allocate(k1(3,L)); allocate(k2(3,L)); allocate(k3(3,L)); allocate(k4(3,L))
        allocate(s_a(3,L)); allocate(s_b(3,L)); allocate(s0_a(3,L)); allocate(s0_b(3,L))
    
    
        allocate(spin_a(3,L, 0:steps/100)); 
        allocate(spin_b(3,L, 0:steps/100))
        allocate(Dxt(L, 0:steps/100))

        ! due to Fortran's column-based storage, use transpose of the intended array
        ! upon reshaping it in fortran: this will be of the dimensions: (steps/100, L, 3)
    
        k1 = 0.0_real128; k2 = 0.0_real128; k3 = 0.0_real128; k4 = 0.0_real128
        s0_a = 0.0_real128; s0_b = 0.0_real128; s_a = 0.0_real128; s_b = 0.0_real128
    
        conf=initseed
        write(foldernm,'(i8)')conf 
    
        ! READ(seedinput,*)seedz
    
        ! print all the input params
        print *, 'Pi = ',  Pi
        print *, 'L = ', L
        print *, 'steps = ', steps
        print *, 'dt = ', dt
        print *, 'lambda = ', lambda
        print *, 'mu = ', mu
        print *, 'J1 = ', J1 
        print *, 'J2 = ', J2
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
    
