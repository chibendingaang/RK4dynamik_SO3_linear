!program for Heisenberg and or driven dynamik
!define module
module param!harsh
	implicit none
	save
	!List of Input paramters and constants
real, parameter :: PI =4.d0*atan(1.d0)
integer, parameter :: L = 128
integer, parameter :: steps = 10000
double precision, parameter :: dt = 0.01d0

double precision, parameter :: lambda = 1.0d0
double precision, parameter :: mu = 0.0d0
double precision, parameter :: epsil = 0.01d0
double precision, parameter :: a = dt*(lambda+mu), b = dt*(lambda-mu)

double precision, allocatable :: s_a(:), s0_a(:)
double precision, allocatable :: s_b(:), s0_b(:)
Double precision, allocatable :: k1(:), k2(:), k3(:), k4(:)
end module param

Program Dynamik
use param!harsh


Implicit None


character flindx*20!harsh
integer :: ti
integer :: x, xi
integer :: conf
CHARACTER*15:: Buf,foldernm
double precision :: sz_cos, sz_sin

allocate(k1(3*L)); allocate(k2(3*L)); allocate(k3(3*L)); allocate(k4(3*L))
allocate(s_a(3*L)); allocate(s_b(3*L)); allocate(s0_a(3*L)); allocate(s0_b(3*L))
k1 = 0.0d0; k2 = 0.0d0; k3 = 0.0d0; k4 = 0.0d0
s0_a = 0.0d0; s0_b = 0.0d0; s_a = 0.0d0; s_b = 0.0d0

!print all the input params
print *, PI
print *, 'L = ', L
print *, 'steps = ', steps
print *, 'dt = ', dt
print *, 'lambda = ', lambda
print *, 'mu = ', mu
print *, 'epsilon = ', epsil
print *, 'a = dt*(Lambda+Mu) = ', a
print *, 'b = dt*(Lambda-Mu) = ', b

!return

CALL GETARG(1,Buf) 
READ(Buf,*)conf
write(foldernm,'(i8)')conf!harsh

!Initialization
call InitRandom(s0_a)
s_a = s0_a

ti=0; write(flindx,'(i8)') ti
 open(111, file = './L128/init'//trim(adjustl(foldernm))//'/spin_a_'//trim(adjustl(flindx))//'.dat', status = 'new')!harsh
 do x = 1,3*L
 write(111,*) s0_a(x) !Sa(x), Sb(x)
 end do 
 close(111)

!Performing the numerical dynamics
do ti = 1,steps
	
	write(flindx,'(i8)') ti  
     Call RK4dynamik(s_a)
    	if (mod(ti,50)==0) then
		open(112, file = './L128/init'//trim(adjustl(foldernm))//'/spin_a_'//trim(adjustl(flindx))//'.dat', status = 'new')!harsh
		do x = 1,3*L
		write(112,*) s_a(x) !Sa(x), Sb(x)
		end do 
    	close(112)
		end if
end do

s0_b = s0_a
sz_cos = s0_b(3)
sz_sin = sqrt(1-sz_cos**2)

s0_b(3) = sz_cos -epsil*sz_sin
s0_b(2) = s0_b(2)*(1 + epsil*sz_cos/sz_sin)
s0_b(1) = s0_b(1)*(1 + epsil*sz_cos/sz_sin)

!using the perturbation cos(theta+epsil) = cos(theta)cos(epsil) - sin(theta)*sin(epsil)
!and... sin(theta+epsil) = sin(theta)cos(epsiL) + cos(theta)sin(epsil)
!sin(theta+epsil)cos(phi) = cos(phi)sin(theta) + cos(phi)*cos(theta)*epsil

s_b = s0_b
ti=0; write(flindx,'(i8)') ti
 open(113, file = './L128/init'//trim(adjustl(foldernm))//'/spin_b_'//trim(adjustl(flindx))//'.dat', status = 'new')!harsh
 do x = 1,3*L
 write(113,*) s0_b(x)  
 end do 
 close(113)

!Same dynamics for the initially perturbed system
do ti = 1,steps
	write(flindx,'(i8)') ti
	Call RK4dynamik(s_b)
    	if (mod(ti,50)==0) then
		open(114, file = './L128/init'//trim(adjustl(foldernm))//'/spin_b_'//trim(adjustl(flindx))//'.dat', status = 'new')!harsh
		do x = 1,3*L
		write(114,*) s_b(x) !Sa(x), Sb(x)
		end do 
    	close(114)
		end if 
end do

End Program Dynamik

!Contains
!FUNCTION to define periodic boundary conditions
Integer function bc(x,L)
    Integer, intent(in) :: x, L
    
    if (x==L+1) then
        bc =1
    else if (x==0) then
        bc = L
    else
        bc = x
    end if
End function bc



!SUBROUTINE#1: to initialize spins from a random distribution

Subroutine InitRandom(spin) 
use param
Implicit NOne

!---random seed setup: from ph.surrey.uk---
 ! ----- variables for portable seed setting -----
  INTEGER :: i_seed
  INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
  INTEGER, DIMENSION(1:8) :: dt_seed
  ! ----- end of variables for seed setting -----


Double precision, allocatable :: planar(:)
Double precision, allocatable :: phi(:)
Double precision, intent(out) :: spin(3*L)
Integer :: xi
Double precision, allocatable :: xyplanar(:)


!random_seed date_time idea from physics.surrey.uk
  ! ----- Set up random seed portably -----
!  CALL RANDOM_SEED(size=i_seed)
!  ALLOCATE(a_seed(1:i_seed))
!  CALL RANDOM_SEED(get=a_seed)
!  CALL DATE_AND_TIME(values=dt_seed)
!  a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
!  CALL RANDOM_SEED(put=a_seed)
!  DEALLOCATE(a_seed)



allocate(planar(L), phi(L), xyplanar(L))
!initialize immediately after allocation
planar = 0.0d0; phi = 0.0d0; xyplanar = 0.0d0

!CALL RANDOM_SEED (PUT=seed)
!call random_number(planar(1:L))
!call random_number(phi(1:L))

!phi = pi*(2.*phi -1.)
do xi=1,L
	call random_number(planar(xi))
	xyplanar(xi) = sqrt(1. - planar(xi)**2)
	call random_number(phi(xi))
	!planar(xi) = 2.*planar(xi) - 1.
	
	!phi(xi) = pi*(2.*(phi(xi)) -1.)
    spin(3*xi-2) = xyplanar(xi)*cos(phi(xi))
    spin(3*xi-1) = xyplanar(xi)*sin(phi(xi))
    spin(3*xi)   = planar(xi)
end do

!s_a(3) = s_a(3) + eps
!intended pertrubation, added externally

end Subroutine InitRandom

!SUBROUTINE#2: for the dynamics
Subroutine RK4dynamik(spin) 
use param!harsh
Implicit None

Double precision :: spin(3*L) !intent(inout)

!Other required variables
Integer, external :: bc   
Integer :: xi
    
    

    do xi = 1, L
        k1(3*xi-2) = a*(spin(3*xi-1) * spin(3*bc(xi+1,L))    - spin(3*xi) * spin(3*bc(xi+1,L)-1)) &
                    & + b*(spin(3*xi-1) * spin(3*bc(xi-1,L)) - spin(3*xi) * spin(3*bc(xi-1,L)-1)) 

        k1(3*xi-1) = a*(spin(3*xi) * spin(3*bc(xi+1,L)-2)    - spin(3*xi-2) * spin(3*bc(xi+1,L))) &
                    & + b*(spin(3*xi) * spin(3*bc(xi-1,L)-2) - spin(3*xi-2) * spin(3*bc(xi-1,L)))

        k1(3*xi) = a*(spin(3*xi-2) * spin(3*bc(xi+1,L) -1)    - spin(3*xi-1) * spin(3*bc(xi+1,L) -2)) &
                    & +b*(spin(3*xi-2) * spin(3*bc(xi-1,L)-1) - spin(3*xi-1) * spin(3*bc(xi-1,L) -2))
    end do
    
    do xi = 1, L
        k2(3*xi-2) = a*((spin(3*xi-1) +  k1(3*xi-1)*0.5) * (spin(3*bc(xi+1,L))   + k1(3*bc(xi+1,L))*0.5)   &
                    & -(spin(3*xi)  +    k1(3*xi)*0.5)   * (spin(3*bc(xi+1,L)-1) + k1(3*bc(xi+1,L)-1)*0.5))&
                    & +b*((spin(3*xi-1)+ k1(3*xi-1)*0.5)*(spin(3*bc(xi-1,L))     + k1(3*bc(xi-1,L))*0.5)   &
                    & -(spin(3*xi)   +   k1(3*xi)*0.5)  *(spin(3*bc(xi-1,L)-1)   + k1(3*bc(xi-1,L)-1)*0.5)) 

        k2(3*xi-1) = a*((spin(3*xi)  +  k1(3*xi)*0.5) * (spin(3*bc(xi+1,L)-2) + k1(3*bc(xi+1,L)-2)*0.5) &
                    & - (spin(3*xi-2) + k1(3*xi-2)*0.5) * (spin(3*bc(xi+1,L)) + k1(3*bc(xi+1,L))*0.5) ) &
                    & +b*((spin(3*xi) + k1(3*xi)*0.5) * (spin(3*bc(xi-1,L)-2) + k1(3*bc(xi-1,L)-2)*0.5) &
                    & - (spin(3*xi-2) + k1(3*xi-2)*0.5) * (spin(3*bc(xi-1,L)) + k1(3*bc(xi-1,L))*0.5) )

        k2(3*xi) = a*((spin(3*xi-2)    + k1(3*xi-2)*0.5)*(spin(3*bc(xi+1,L)-1) + k1(3*bc(xi+1,L)-1)*0.5) &
                    & -(spin(3*xi-1)   + k1(3*xi-1)*0.5)*(spin(3*bc(xi+1,L)-2) + k1(3*bc(xi+1,L)-2)*0.5))&
                    & +b*((spin(3*xi-2)+ k1(3*xi-2)*0.5)*(spin(3*bc(xi-1,L)-1) + k1(3*bc(xi-1,L)-1)*0.5) &
                    & - (spin(3*xi-1)  + k1(3*xi-1)*0.5)*(spin(3*bc(xi-1,L)-2) + k1(3*bc(xi-1,L)-2)*0.5))
     end do
     
     do xi = 1, L
        k3(3*xi-2) = a*((spin(3*xi-1)  +  k2(3*xi-1)*0.5)* (spin(3*bc(xi+1,L))   + k2(3*bc(xi+1,L))*0.5) &
                    & -(spin(3*xi)     +  k2(3*xi)*0.5)  * (spin(3*bc(xi+1,L)-1) + k2(3*bc(xi+1,L)-1)*0.5)) &
                    & +b*((spin(3*xi-1)+  k2(3*xi-1)*0.5)* (spin(3*bc(xi-1,L))   + k2(3*bc(xi-1,L))*0.5) &
                    & -(spin(3*xi)     +  k2(3*xi)*0.5)  * (spin(3*bc(xi-1,L)-1) + k2(3*bc(xi-1,L)-1)*0.5))

        k3(3*xi-1) = a*((spin(3*xi)   + k2(3*xi)*0.5) *(spin(3*bc(xi+1,L)-2) + k2(3*bc(xi+1,L)-2)*0.5) &
                    & - (spin(3*xi-2) + k2(3*xi-2)*0.5) *(spin(3*bc(xi+1,L)) + k2(3*bc(xi+1,L))*0.5)) &
                    & +b*((spin(3*xi) + k2(3*xi)*0.5) *(spin(3*bc(xi-1,L)-2) + k2(3*bc(xi-1,L)-2)*0.5) &
                    & - (spin(3*xi-2) + k2(3*xi-2)*0.5) *(spin(3*bc(xi-1,L)) + k2(3*bc(xi-1,L))*0.5))

        k3(3*xi) = a*((spin(3*xi-2)    + k2(3*xi-2)*0.5)*(spin(3*bc(xi+1,L)-1) + k2(3*bc(xi+1,L)-1)*0.5) &
                    & - (spin(3*xi-1)  + k2(3*xi-1)*0.5)*(spin(3*bc(xi+1,L)-2) + k2(3*bc(xi+1,L)-2)*0.5))&
                    & +b*((spin(3*xi-2)+ k2(3*xi-2)*0.5)*(spin(3*bc(xi-1,L)-1) + k2(3*bc(xi-1,L)-1)*0.5) &
                    & - (spin(3*xi-1)  + k2(3*xi-1)*0.5)*(spin(3*bc(xi-1,L)-2) + k2(3*bc(xi-1,L)-2)*0.5)) 
    end do     
    
    do xi = 1, L   
        k4(3*xi-2) = a*((spin(3*xi-1)  + k3(3*xi-1))*(spin(3*bc(xi+1,L))   + k3(3*bc(xi+1,L))) &
                    & -(spin(3*xi)     + k3(3*xi)) * (spin(3*bc(xi+1,L)-1) + k3(3*bc(xi+1,L)-1))) &
                    & +b*((spin(3*xi-1)+ k3(3*xi-1))*(spin(3*bc(xi-1,L))   + k3(3*bc(xi-1,L))) &
                    & -(spin(3*xi)     + k3(3*xi)) * (spin(3*bc(xi-1,L)-1) + k3(3*bc(xi-1,L)-1)))

        k4(3*xi-1) = a*((spin(3*xi)   + k3(3*xi))  *(spin(3*bc(xi+1,L)-2) + k3(3*bc(xi+1,L)-2)) &
                    & - (spin(3*xi-2) + k3(3*xi-2))*(spin(3*bc(xi+1,L))   + k3(3*bc(xi+1,L)))) &
                    & +b*((spin(3*xi) + k3(3*xi))  *(spin(3*bc(xi-1,L)-2) + k3(3*bc(xi-1,L)-2)) &
                    & - (spin(3*xi-2) + k3(3*xi-2))*(spin(3*bc(xi-1,L))   + k3(3*bc(xi-1,L))))

        k4(3*xi) = a*((spin(3*xi-2)    + k3(3*xi-2))*(spin(3*bc(xi+1,L)-1) + k3(3*bc(xi+1,L)-1)) &
                    & - (spin(3*xi-1)  + k3(3*xi-1))*(spin(3*bc(xi+1,L)-2) + k3(3*bc(xi+1,L)-2)))&
                    & +b*((spin(3*xi-2)+ k3(3*xi-2))*(spin(3*bc(xi-1,L)-1) + k3(3*bc(xi-1,L)-1)) &
                    & - (spin(3*xi-1)  + k3(3*xi-1))*(spin(3*bc(xi-1,L)-2) + k3(3*bc(xi-1,L)-2)))
    end do        
    
    do xi = 1, L
        spin(3*xi-2) = spin(3*xi-2) + (k1(3*xi-2) + 2*k2(3*xi-2) + 2*k3(3*xi-2) + k4(3*xi-2))/6.0
        spin(3*xi-1) = spin(3*xi-1) + (k1(3*xi-1) + 2*k2(3*xi-1) + 2*k3(3*xi-1) + k4(3*xi-1))/6.0
        spin(3*xi)   = spin(3*xi)   + (k1(3*xi)   + 2*k2(3*xi)   + 2*k3(3*xi)   + k4(3*xi))/6.0
!redundant: S_next instead of S_
    end do


End subroutine RK4dynamik


    





# RK4dynamik_SO3_linear
