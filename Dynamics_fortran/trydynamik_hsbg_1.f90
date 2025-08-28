!program for Heisenberg and or driven dynamik


MODULE mt19937
IMPLICIT NONE
INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12, 60)

! Period parameters
INTEGER, PARAMETER :: n = 624, n1 = n+1, m = 397, mata = -1727483681
!                                    constant vector a
INTEGER, PARAMETER :: umask = -2147483647 - 1
!                                    most significant w-r bits
INTEGER, PARAMETER :: lmask =  2147483647
!                                    least significant r bits
! Tempering parameters
INTEGER, PARAMETER :: tmaskb= -1658038656, tmaskc= -272236544

!                     the array for the state vector
INTEGER, SAVE      :: mt(0:n-1), mti = n1
!                     mti==N+1 means mt[N] is not initialized

PRIVATE
PUBLIC :: dp, sgrnd, grnd, init_genrand

CONTAINS


SUBROUTINE sgrnd(seed)
! This is the original version of the seeding routine.
! It was replaced in the Japanese version in C on 26 January 2002
! It is recommended that routine init_genrand is used instead.

INTEGER, INTENT(IN)   :: seed

!    setting initial seeds to mt[N] using the generator Line 25 of Table 1 in
!    [KNUTH 1981, The Art of Computer Programming Vol. 2 (2nd Ed.), pp102]

mt(0)= IAND(seed, -1)
DO  mti=1,n-1
  mt(mti) = IAND(69069 * mt(mti-1), -1)
END DO

RETURN
END SUBROUTINE sgrnd
!***********************************************************************

SUBROUTINE init_genrand(seed)
! This initialization is based upon the multiplier given on p.106 of the
! 3rd edition of Knuth, The Art of Computer Programming Vol. 2.

! This version assumes that integer overflow does NOT cause a crash.

INTEGER, INTENT(IN)  :: seed

INTEGER  :: latest

mt(0) = seed
latest = seed
DO mti = 1, n-1
  latest = IEOR( latest, ISHFT( latest, -30 ) )
  latest = latest * 1812433253 + mti
  mt(mti) = latest
END DO

RETURN
END SUBROUTINE init_genrand
!***********************************************************************

FUNCTION grnd() RESULT(fn_val)
REAL (dp) :: fn_val

INTEGER, SAVE :: mag01(0:1) = (/ 0, mata /)
!                        mag01(x) = x * MATA for x=0,1
INTEGER       :: kk, y

! These statement functions have been replaced with separate functions
! tshftu(y) = ISHFT(y,-11)
! tshfts(y) = ISHFT(y,7)
! tshftt(y) = ISHFT(y,15)
! tshftl(y) = ISHFT(y,-18)

IF(mti >= n) THEN
!                       generate N words at one time
  IF(mti == n+1) THEN
!                            if sgrnd() has not been called,
    CALL sgrnd(4357)
!                              a default initial seed is used
  END IF
  
  DO  kk = 0, n-m-1
    y = IOR(IAND(mt(kk),umask), IAND(mt(kk+1),lmask))
    mt(kk) = IEOR(IEOR(mt(kk+m), ISHFT(y,-1)),mag01(IAND(y,1)))
  END DO
  DO  kk = n-m, n-2
    y = IOR(IAND(mt(kk),umask), IAND(mt(kk+1),lmask))
    mt(kk) = IEOR(IEOR(mt(kk+(m-n)), ISHFT(y,-1)),mag01(IAND(y,1)))
  END DO
  y = IOR(IAND(mt(n-1),umask), IAND(mt(0),lmask))
  mt(n-1) = IEOR(IEOR(mt(m-1), ISHFT(y,-1)),mag01(IAND(y,1)))
  mti = 0
END IF

y = mt(mti)
mti = mti + 1
y = IEOR(y, tshftu(y))
y = IEOR(y, IAND(tshfts(y),tmaskb))
y = IEOR(y, IAND(tshftt(y),tmaskc))
y = IEOR(y, tshftl(y))

IF(y < 0) THEN
  fn_val = (DBLE(y) + 2.0D0**32) / (2.0D0**32 - 1.0D0)
ELSE
  fn_val = DBLE(y) / (2.0D0**32 - 1.0D0)
END IF

RETURN
END FUNCTION grnd


FUNCTION tshftu(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,-11)
RETURN
END FUNCTION tshftu


FUNCTION tshfts(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,7)
RETURN
END FUNCTION tshfts


FUNCTION tshftt(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,15)
RETURN
END FUNCTION tshftt


FUNCTION tshftl(y) RESULT(fn_val)
INTEGER, INTENT(IN) :: y
INTEGER             :: fn_val

fn_val = ISHFT(y,-18)
RETURN
END FUNCTION tshftl

END MODULE mt19937



!define module
module param!harsh
	implicit none
	save
	!List of Input paramters and constants
double precision, parameter :: PI =4.d0*atan(1.d0)!harsh
integer, parameter :: L = 128!harsh
integer, parameter :: steps = 4000!harsh
double precision, parameter :: dt = 0.001d0!harsh

double precision, parameter :: lambda = 1.0d0!harsh
double precision, parameter :: mu = 0.0d0!harsh
double precision, parameter :: epsil = 0.001d0
double precision, parameter :: a = dt*(lambda+mu), b = dt*(lambda-mu)

Double precision, allocatable :: zaxial(:)
Double precision, allocatable :: phi(:)
Double precision, allocatable :: xyplanar(:)

integer :: conf
double precision, allocatable :: s_a(:), sxt_a(:), S_at(:,:)
double precision, allocatable :: s_b(:), sxt_b(:), S_bt(:,:)
double precision, allocatable :: decorrx(:), Dxt(:,:)
Double precision, allocatable :: k1(:), k2(:), k3(:), k4(:)

end module param


Program Dynamik
use param!harsh
Implicit None

character flindx*20!harsh
integer :: ti, j
integer :: x
!integer :: conf
CHARACTER*15:: Buf,foldernm
!double precision :: sz_cos, sz_sin
double precision :: holdprod
double precision, allocatable :: dummyzero(:)


allocate(dummyzero(3*L))
allocate(zaxial(L)); allocate(phi(L)); allocate(xyplanar(L))
allocate(k1(3*L)); allocate(k2(3*L)); allocate(k3(3*L)); allocate(k4(3*L))
allocate(s_a(3*L)); allocate(s_b(3*L)); 
allocate(sxt_a(3*L*(int(steps/100)+1))); allocate(sxt_b(3*L*(int(steps/100)+1)))



allocate(S_at(1:3*L, 0:int(steps/100))); allocate(S_bt(1:3*L,0:int(steps/100)))
!the above have to be transposed and dot producted to fit the dimensions of Dxt
allocate(Dxt(0:int(steps/100), 1:L))
allocate(decorrx(1:L))

k1 = 0.0d0; k2 = 0.0d0; k3 = 0.0d0; k4 = 0.0d0
!sxt_a = 0.0d0; sxt_b = 0.0d0; 
s_a = 0.0d0; s_b = 0.0d0
!Dxt = 0.0d0
decorrx = 1.0d0

CALL GETARG(1,Buf) 
READ(Buf,*) conf
WRITE(foldernm,'(i8)')conf!harsh


!print all the input params
print *, 'Pi = ',  PI
print *, 'L = ', L
print *, 'steps = ', steps
print *, 'dt = ', dt
print *, 'lambda = ', lambda
print *, 'mu = ', mu
print *, 'epsilon = ', epsil
print *, 'a = dt*(Lambda+Mu) = ', a
print *, 'b = dt*(Lambda-Mu) = ', b


!Initialization
!s0_a is drawn from initrandom subroutine
!s0_b is the same configuration as s0_a with a perturbation added at one site
call InitRandom(s_a, 0.0d0)
!s_a = s0_a

ti=0; write(flindx,'(i8)') ti
 !open(11, file = './L128/hsbg/init'//trim(adjustl(foldernm))//'/spin_a_'//trim(adjustl(flindx))//'.dat', status = 'new')!harsh
 open(12, file = './L128/hsbg/init'//trim(adjustl(foldernm))//'/Spin_a.dat', status = 'new')!nsb
 do x = 1,3*L
 !write(11,*) s_a(x) 
 write(12,*) s_a(x)
 end do 
 !close(11)
 close(12)

!Performing the numerical dynamics for configuration_a
do ti = 1,steps
	write(flindx,'(i8)') ti
	open(12, file = './L128/hsbg/init'//trim(adjustl(foldernm))//'/Spin_a.dat', status = 'old', position = 'append')!nsb
    Call RK4dynamik(s_a)
    	if (mod(ti,100)==0) then
		do x = 1,3*L
		write(12,*) s_a(x)
		end do 
		end if
	close(12)
end do


deallocate(s_a); !allocate(s_a(3*L)); s_a = 0.0d0


call InitRandom(s_b, epsil)
!s_b = s0_b

!a primer on storing the data to the files correctly
ti=0; write(flindx,'(i8)') ti
 !open(21, file = './L128/hsbg/init'//trim(adjustl(foldernm))//'/spin_b_'//trim(adjustl(flindx))//'.dat', status = 'new')!harsh
 open(22, file = './L128/hsbg/init'//trim(adjustl(foldernm))//'/Spin_b.dat', status = 'new')
 
 !read(12,*) S_at
 !transposing a dynamical array will change its shape, hopefully!

 do x = 1,3*L

 !write(21,*) s_b(x) 
 write(22,*) s_b(x)

 end do 

 !close(21)
 close(22)


!print *, shape(S_at)
!go to 324
!must reset the RK4 subroutine to get spinb correct
dummyzero = 0.0d0
Call RK4dynamik(dummyzero)


!Performing the numerical dynamics
do ti = 1,steps
	write(flindx,'(i8)') ti
	
	open(22, file = './L128/hsbg/init'//trim(adjustl(foldernm))//'/Spin_b.dat', status = 'old', position = 'append')!nsb
    Call RK4dynamik(s_b)
    	if (mod(ti,100)==0) then
    	
		do x = 1,3*L
		 !this is the crucial step: read the correct position-time step from the 1-D array representation of the 2-D S(x,t) array
		!read(12,*) S_at(3*L*int(ti/100) + x)
		!print *, S_at(int(ti/100),x), s_b(x)
		write(22,*) s_b(x)
		end do 
		end if
	
	close(22)
	
end do

deallocate(s_b)

!324 continue

 open(12, file = './L128/hsbg/init'//trim(adjustl(foldernm))//'/Spin_a.dat', status = 'old')!nsb
 open(22, file = './L128/hsbg/init'//trim(adjustl(foldernm))//'/Spin_b.dat', status = 'old')!nsb
 open(31, file = './L128/hsbg/init'//trim(adjustl(foldernm))//'/Dxt.dat', status = 'new', position='append')

read(12,*) sxt_a
read(22,*) sxt_b

 !!!
print *
print *, shape(sxt_a)
print *, shape(sxt_b)

S_at = reshape(sxt_a, (/3*L,int(steps/100)+1/))
S_at = transpose(S_at)

S_bt = reshape(sxt_b, (/3*L,int(steps/100)+1/))
S_bt = transpose(S_bt)

print *, shape(S_at)
print *, shape(S_bt)
 !!!

do ti = 0, int(steps/100)
	do x = 1,L
		holdprod = 0.0d0
		do j = 1,3
			holdprod = holdprod + S_at(ti,3*x-3+j)*S_bt(ti,3*x-3+j)
		end do 
		
		Dxt(ti,x) = decorrx(x) - holdprod 
		print *, Dxt(ti,x)
		write(31,*) Dxt(ti,x)
	end do
	!print *, Dxt(ti)
end do 


 close(12)
 close(22)
 close(31)



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

Subroutine InitRandom(spin, eps) 
use param
USE mt19937
Implicit NOne


!double precision, intent(in) :: eps 
Double precision, intent(out) :: spin(3*L)
double precision, intent(in) :: eps
Integer :: xi
integer :: seed1 !, seed2
seed1 = conf;

!CALL RANDOM_SEED ()
!call random_number(zaxial)
!call random_number(phi)

CALL init_genrand(seed1)	

do xi=1,L
	zaxial(xi) = grnd()
	phi(xi) = grnd()
	if (xi==1) then
	zaxial(xi) = zaxial(xi) + eps
	end if
	
	phi(xi) = pi*(2.*(phi(xi)) -1.)
	zaxial(xi) = 2.*zaxial(xi) - 1.
	!call random_number(zaxial(xi))
	!call random_number(phi(xi))
	xyplanar(xi) = sqrt(1. - zaxial(xi)**2)
	!if (mod(xi,2)==1) then 
	!
	!xyplanar(xi) = 1
    !else 
    !phi(xi) = 0
    !xyplanar(xi) = 0
    !zaxial(xi) = 1.0
    !end if
    spin(3*xi-2) = xyplanar(xi)*cos(phi(xi))
    spin(3*xi-1) = xyplanar(xi)*sin(phi(xi))
    spin(3*xi)   = zaxial(xi)
    
end do

!spin(3) = zaxial(3) + eps
!spin(3) = spin(3) + eps
!intended pertrubation, added externally

end Subroutine InitRandom



!SUBROUTINE#2: for the dynamics
Subroutine RK4dynamik(spin)!harsh, replace S_ by spin_
use param
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
