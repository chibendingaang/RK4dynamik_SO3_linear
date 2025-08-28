module alphaspindynamics

use xdalphaparam
use mt19937xd
! use periodic_bc

implicit none 
! Integer, parameter :: dp = selected_real_kind(15,307)

Contains

! FUNCTION to define periodic boundary conditions
! This is to be excluded in alphaspindynamics
! due to fixed but time dependent boundary conditions 
!goto 45
! Integer function bc(x,L)
! Integer, intent(in) :: x, L
    
!    if (x==L+1) then
!        bc =1
!    else if (x==0) then
!        bc = L
!    else
!        bc = x
!    end if

! End function bc
!45 continue


! SUBROUTINE#1: to initialize spins from a random distribution
Subroutine InitRandom(spin, eps) 

!Implicit None
! double precision, intent(in) :: eps 
real(xd), intent(out) :: spin(3*L)
real(xd), intent(in) :: eps
Integer :: xi
integer :: seed1 !, seed2
seed1 = conf;

! CALL RANDOM_SEED ()
! call random_number(planar)
! call random_number(phi)
print*,"seed1=",seed1
CALL init_genrand(seed1)

do xi=1,L
        planar(xi) = grnd()
        phi(xi) = grnd()
        if (xi==1) then
        planar(xi) = planar(xi) + eps
        end if
        phi(xi) = pi*(2*(phi(xi)) -1.0_xd)
        planar(xi) = 2*planar(xi) - 1.0_xd
        xyplanar(xi) = sqrt(1.0_xd - planar(xi)**2)
    spin(3*xi-2) = xyplanar(xi)*cos(phi(xi))
    spin(3*xi-1) = xyplanar(xi)*sin(phi(xi))
    spin(3*xi)   = planar(xi)
    
end do

! if x = 2, L-1; then the following spin values are required:
! spin(1) = spin(4); spin(2) = spin(5); spin(3) = spin(6)
! spin(3*L-2) = spin(3*L-5); spin(3*L-1) = spin(3*L-4); spin(3*L) = spin(3*L-3)

end Subroutine InitRandom



! SUBROUTINE#2: for the dynamics
Subroutine RK4dynamik(spin)

real(xd) :: spin(3*L) !intent(inout)

! Other required variables
! a = lam + mu --> 1 ; b = lam - mu --> alpha
Integer :: xi
    k1(1) = a*(spin(2) * spin(6)    - spin(3) * spin(5))
    k1(2) = a*(spin(3) * spin(4)    - spin(1) * spin(6)) 
    k1(3) = a*(spin(1) * spin(5)    - spin(2) * spin(4)) 

    k1(3*L-2) = b*(spin(3*L-1) * spin(3*(L-1))    - spin(3*L) * spin(3*(L-1)-1))
    k1(3*L-1) = b*(spin(3*L) * spin(3*(L-1)-2)    - spin(3*L-2) * spin(3*(L-1))) 
    k1(3*L)   = b*(spin(3*L-2) * spin(3*(L-1)-1)    - spin(3*L-1) * spin(3*(L-1)-2)) 
    
    
    do xi = 2, L-1
        k1(3*xi-2) = a*(spin(3*xi-1) * spin(3*(xi+1))    - spin(3*xi) * spin(3*(xi+1)-1)) &
                    & + b*(spin(3*xi-1) * spin(3*(xi-1)) - spin(3*xi) * spin(3*(xi-1)-1)) 

        k1(3*xi-1) = a*(spin(3*xi) * spin(3*(xi+1)-2)    - spin(3*xi-2) * spin(3*(xi+1))) &
                    & + b*(spin(3*xi) * spin(3*(xi-1)-2) - spin(3*xi-2) * spin(3*(xi-1)))

        k1(3*xi) = a*(spin(3*xi-2) * spin(3*(xi+1) -1)    - spin(3*xi-1) * spin(3*(xi+1) - 2)) &
                    & +b*(spin(3*xi-2) * spin(3*(xi-1)-1) - spin(3*xi-1) * spin(3*(xi-1) - 2))
    end do
    
    
    k2(1) = a*((spin(2) + 0.5*k1(2)) * (spin(6) + 0.5*k1(6))   - (spin(3) + 0.5*k1(3)) * (spin(5) + 0.5*k1(5)))
    k2(2) = a*((spin(3) + 0.5*k1(3)) * (spin(4) + 0.5*k1(4))   - (spin(1) + 0.5*k1(1)) * (spin(6) + 0.5*k1(6))) 
    k2(3) = a*((spin(1) + 0.5*k1(1)) * (spin(5) + 0.5*k1(5))   - (spin(2) + 0.5*k1(2)) * (spin(4) + 0.5*k1(4))) 
    
    k2(3*L-2) = b*((spin(3*L-1) + 0.5*k1(3*L-1)) * (spin(3*(L-1)) + 0.5*k1(3*(L-1)))   &
                 - (spin(3*L) + 0.5*k1(3*L)) * (spin(3*(L-1)-1) + 0.5*k1(3*(L-1)-1)))
                 
    k2(3*L-1) = b*((spin(3*L) + 0.5*k1(3*L)) * (spin(3*(L-1)-2) + 0.5*k1(3*(L-1)-2))  &
                 - (spin(3*L-2) + 0.5*k1(3*L-2)) * (spin(3*(L-1)) + 0.5*k1(3*(L-1)))) 
                 
    k2(3*L)   = b*((spin(3*L-2) + 0.5*k1(3*L-2)) * (spin(3*(L-1)-1) + 0.5* k1(3*(L-1)-1)) &
                 - (spin(3*L-1) + 0.5*k1(3*L-1)) * (spin(3*(L-1)-2) + 0.5*k1(3*(L-1)-2))) 
    
    do xi = 2, L-1
        k2(3*xi-2) = a*((spin(3*xi-1) +  k1(3*xi-1)*0.5) * (spin(3*(xi+1))   + k1(3*(xi+1))*0.5)   &
                    & -(spin(3*xi)  +    k1(3*xi)*0.5)   * (spin(3*(xi+1)-1) + k1(3*(xi+1)-1)*0.5))&
                    & +b*((spin(3*xi-1)+ k1(3*xi-1)*0.5)*(spin(3*(xi-1))     + k1(3*(xi-1))*0.5)   &
                    & -(spin(3*xi)   +   k1(3*xi)*0.5)  *(spin(3*(xi-1)-1)   + k1(3*(xi-1)-1)*0.5)) 

        k2(3*xi-1) = a*((spin(3*xi)  +  k1(3*xi)*0.5) * (spin(3*(xi+1)-2) + k1(3*(xi+1)-2)*0.5) &
                    & - (spin(3*xi-2) + k1(3*xi-2)*0.5) * (spin(3*(xi+1)) + k1(3*(xi+1))*0.5) ) &
                    & +b*((spin(3*xi) + k1(3*xi)*0.5) * (spin(3*(xi-1)-2) + k1(3*(xi-1)-2)*0.5) &
                    & - (spin(3*xi-2) + k1(3*xi-2)*0.5) * (spin(3*(xi-1)) + k1(3*(xi-1))*0.5) )

        k2(3*xi) = a*((spin(3*xi-2)    + k1(3*xi-2)*0.5)*(spin(3*(xi+1)-1) + k1(3*(xi+1)-1)*0.5) &
                    & -(spin(3*xi-1)   + k1(3*xi-1)*0.5)*(spin(3*(xi+1)-2) + k1(3*(xi+1)-2)*0.5))&
                    & +b*((spin(3*xi-2)+ k1(3*xi-2)*0.5)*(spin(3*(xi-1)-1) + k1(3*(xi-1)-1)*0.5) &
                    & - (spin(3*xi-1)  + k1(3*xi-1)*0.5)*(spin(3*(xi-1)-2) + k1(3*(xi-1)-2)*0.5))
    end do
     
    k3(1) = a*((spin(2) + 0.5*k2(2)) * (spin(6) + 0.5*k2(6))   - (spin(3) + 0.5*k2(3)) * (spin(5) + 0.5*k2(5)))
    k3(2) = a*((spin(3) + 0.5*k2(3)) * (spin(4) + 0.5*k2(4))   - (spin(1) + 0.5*k2(1)) * (spin(6) + 0.5*k2(6))) 
    k3(3) = a*((spin(1) + 0.5*k2(1)) * (spin(5) + 0.5*k2(5))   - (spin(2) + 0.5*k2(2)) * (spin(4) + 0.5*k2(4))) 
    
    k3(3*L-2) = b*((spin(3*L-1) + 0.5*k2(3*L-1)) * (spin(3*(L-1)) + 0.5*k2(3*(L-1)))   &
                 - (spin(3*L) + 0.5*k2(3*L)) * (spin(3*(L-1)-1) + 0.5*k2(3*(L-1)-1)))
                 
    k3(3*L-1) = b*((spin(3*L) + 0.5*k2(3*L)) * (spin(3*(L-1)-2) + 0.5*k2(3*(L-1)-2))  &
                 - (spin(3*L-2) + 0.5*k2(3*L-2)) * (spin(3*(L-1)) + 0.5*k2(3*(L-1)))) 
                 
    k3(3*L)   = b*((spin(3*L-2) + 0.5*k2(3*L-2)) * (spin(3*(L-1)-1) + 0.5* k2(3*(L-1)-1)) &
                 - (spin(3*L-1) + 0.5*k2(3*L-1)) * (spin(3*(L-1)-2) + 0.5*k2(3*(L-1)-2))) 
    
    do xi = 2, L-1
        k3(3*xi-2) = a*((spin(3*xi-1)  +  k2(3*xi-1)*0.5)* (spin(3*(xi+1))   + k2(3*(xi+1))*0.5) &
                    & -(spin(3*xi)     +  k2(3*xi)*0.5)  * (spin(3*(xi+1)-1) + k2(3*(xi+1)-1)*0.5)) &
                    & +b*((spin(3*xi-1)+  k2(3*xi-1)*0.5)* (spin(3*(xi-1))   + k2(3*(xi-1))*0.5) &
                    & -(spin(3*xi)     +  k2(3*xi)*0.5)  * (spin(3*(xi-1)-1) + k2(3*(xi-1)-1)*0.5))

        k3(3*xi-1) = a*((spin(3*xi)   + k2(3*xi)*0.5) *(spin(3*(xi+1)-2) + k2(3*(xi+1)-2)*0.5) &
                    & - (spin(3*xi-2) + k2(3*xi-2)*0.5) *(spin(3*(xi+1)) + k2(3*(xi+1))*0.5)) &
                    & +b*((spin(3*xi) + k2(3*xi)*0.5) *(spin(3*(xi-1)-2) + k2(3*(xi-1)-2)*0.5) &
                    & - (spin(3*xi-2) + k2(3*xi-2)*0.5) *(spin(3*(xi-1)) + k2(3*(xi-1))*0.5))

        k3(3*xi) = a*((spin(3*xi-2)    + k2(3*xi-2)*0.5)*(spin(3*(xi+1)-1) + k2(3*(xi+1)-1)*0.5) &
                    & - (spin(3*xi-1)  + k2(3*xi-1)*0.5)*(spin(3*(xi+1)-2) + k2(3*(xi+1)-2)*0.5))&
                    & +b*((spin(3*xi-2)+ k2(3*xi-2)*0.5)*(spin(3*(xi-1)-1) + k2(3*(xi-1)-1)*0.5) &
                    & - (spin(3*xi-1)  + k2(3*xi-1)*0.5)*(spin(3*(xi-1)-2) + k2(3*(xi-1)-2)*0.5)) 
    end do     
    
    k4(1) = a*((spin(2) + k3(2)) * (spin(6) + k3(6))   - (spin(3) + k3(3)) * (spin(5) + k3(5)))
    k4(2) = a*((spin(3) + k3(3)) * (spin(4) + k3(4))   - (spin(1) + k3(1)) * (spin(6) + k3(6))) 
    k4(3) = a*((spin(1) + k3(1)) * (spin(5) + k3(5))   - (spin(2) + k3(2)) * (spin(4) + k3(4))) 
    
    k4(3*L-2) = b*((spin(3*L-1) + k3(3*L-1)) * (spin(3*(L-1)) + k3(3*(L-1)))   &
                 - (spin(3*L) +  k3(3*L)) * (spin(3*(L-1)-1) + k3(3*(L-1)-1)))
                 
    k4(3*L-1) = b*((spin(3*L) + k3(3*L)) * (spin(3*(L-1)-2) + k3(3*(L-1)-2))  &
                 - (spin(3*L-2) + k3(3*L-2)) * (spin(3*(L-1)) + k3(3*(L-1)))) 
                 
    k4(3*L)   = b*((spin(3*L-2) + k3(3*L-2)) * (spin(3*(L-1)-1) + k3(3*(L-1)-1)) &
                 - (spin(3*L-1) + k3(3*L-1)) * (spin(3*(L-1)-2) + k3(3*(L-1)-2))) 
    
    do xi = 2, L-1   
        k4(3*xi-2) = a*((spin(3*xi-1)  + k3(3*xi-1))*(spin(3*(xi+1))   + k3(3*(xi+1))) &
                    & -(spin(3*xi)     + k3(3*xi)) * (spin(3*(xi+1)-1) + k3(3*(xi+1)-1))) &
                    & +b*((spin(3*xi-1)+ k3(3*xi-1))*(spin(3*(xi-1))   + k3(3*(xi-1))) &
                    & -(spin(3*xi)     + k3(3*xi)) * (spin(3*(xi-1)-1) + k3(3*(xi-1)-1)))

        k4(3*xi-1) = a*((spin(3*xi)   + k3(3*xi))  *(spin(3*(xi+1)-2) + k3(3*(xi+1)-2)) &
                    & - (spin(3*xi-2) + k3(3*xi-2))*(spin(3*(xi+1))   + k3(3*(xi+1)))) &
                    & +b*((spin(3*xi) + k3(3*xi))  *(spin(3*(xi-1)-2) + k3(3*(xi-1)-2)) &
                    & - (spin(3*xi-2) + k3(3*xi-2))*(spin(3*(xi-1))   + k3(3*(xi-1))))

        k4(3*xi) = a*((spin(3*xi-2)    + k3(3*xi-2))*(spin(3*(xi+1)-1) + k3(3*(xi+1)-1)) &
                    & - (spin(3*xi-1)  + k3(3*xi-1))*(spin(3*(xi+1)-2) + k3(3*(xi+1)-2)))&
                    & +b*((spin(3*xi-2)+ k3(3*xi-2))*(spin(3*(xi-1)-1) + k3(3*(xi-1)-1)) &
                    & - (spin(3*xi-1)  + k3(3*xi-1))*(spin(3*(xi-1)-2) + k3(3*(xi-1)-2)))
    end do        
    
    do xi = 1, L
        spin(3*xi-2) = spin(3*xi-2) + (k1(3*xi-2) + 2*k2(3*xi-2) + 2*k3(3*xi-2) + k4(3*xi-2))/6.0
        spin(3*xi-1) = spin(3*xi-1) + (k1(3*xi-1) + 2*k2(3*xi-1) + 2*k3(3*xi-1) + k4(3*xi-1))/6.0
        spin(3*xi)   = spin(3*xi)   + (k1(3*xi)   + 2*k2(3*xi)   + 2*k3(3*xi)   + k4(3*xi))/6.0
    end do


End subroutine RK4dynamik
End module alphaspindynamics
