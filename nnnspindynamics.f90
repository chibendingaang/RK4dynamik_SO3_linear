module nnnspindynamics

use xdnnnparam
use mt19937xd
! use periodic_bc

implicit none 
! Integer, parameter :: dp = selected_real_kind(15,307)

Contains

!FUNCTION to define periodic boundary conditions

!FUNCTION to define periodic boundary conditions
Integer function bc(x,L)
    Integer, intent(in) :: x, L
    
    if (x<=0) then
        bc = x+L
    else if (x>0 .and. x<=L) then
        bc = x
    else
        bc = mod(x,L)
    end if

End function bc


!SUBROUTINE#1: to initialize spins from a random distribution

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
        phi(xi) = pi*(2*(phi(xi)) - 1.0_xd)
        planar(xi) = 2*planar(xi) - 1.0_xd
        if (xi==L/2) then
        planar(xi) = planar(xi) + eps
        end if
        ! important difference in initialization
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
Subroutine RK4dynamik_nnn(spin)

real(xd) :: spin(3*L) !intent(inout)

real(xd) :: sx_1, sx_2, sx_3
real(xd) :: sxleft_1, sxleft_2, sxleft_3, sxrght_1, sxrght_2, sxrght_3
real(xd) :: sxlleftt_1, sxlleftt_2, sxlleftt_3, sxrrghtt_1, sxrrghtt_2, sxrrghtt_3

real(xd) :: k1_1, k1_2, k1_3, k2_1, k2_2, k2_3, k3_1, k3_2, k3_3
real(xd) :: k1left_1, k1left_2, k1left_3, k1rght_1, k1rght_2, k1Rght_3
real(xd) :: k1LLeftt_1, k1LLeftt_2, k1LLeftt_3, k1Rrghtt_1, k1Rrghtt_2, k1Rrghtt_3
real(xd) :: k2left_1, k2left_2, k2left_3, k2rght_1, k2rght_2, k2Rght_3
real(xd) :: k2LLeftt_1, k2LLeftt_2, k2LLeftt_3, k2Rrghtt_1, k2Rrghtt_2, k2Rrghtt_3
real(xd) :: k3left_1, k3left_2, k3left_3, k3rght_1, k3rght_2, k3Rght_3
real(xd) :: k3LLeftt_1, k3LLeftt_2, k3LLeftt_3, k3Rrghtt_1, k3Rrghtt_2, k3Rrghtt_3

! Other required variables
! a = lam + mu --> 1 ; b = lam - mu --> alpha
Integer :: xi
    !k1(1) = a*(spin(2) * spin(6)    - spin(3) * spin(5))
    !k1(2) = a*(spin(3) * spin(4)    - spin(1) * spin(6)) 
    !k1(3) = a*(spin(1) * spin(5)    - spin(2) * spin(4)) 

    !k1(3*L-2) = b*(spin(3*L-1) * spin(3*(L-1))    - spin(3*L) * spin(3*(L-1)-1))
    !k1(3*L-1) = b*(spin(3*L) * spin(3*(L-1)-2)    - spin(3*L-2) * spin(3*(L-1))) 
    !k1(3*L)   = b*(spin(3*L-2) * spin(3*(L-1)-1)    - spin(3*L-1) * spin(3*(L-1)-2)) 
    
    
    do xi = 1, L
        sx_1 = spin(3*xi-2); sx_2 = spin(3*xi-1); sx_3 = spin(3*xi)

        sxleft_1 = J1*spin(3*bc(xi-1,L)-2); sxleft_2 = J1*spin(3*bc(xi-1,L)-1); sxleft_3 = J1*spin(3*bc(xi-1,L))
        sxrght_1 = J1*spin(3*bc(xi+1,L)-2); sxrght_2 = J1*spin(3*bc(xi+1,L)-1); sxrght_3 = J1*spin(3*bc(xi+1,l))

        sxlleftt_1 = J2*spin(3*bc(xi-2,L)-2); sxlleftt_2 = J2*spin(3*bc(xi-2,L)-1); sxlleftt_3 = J2*spin(3*bc(xi-2,L))
        sxrrghtt_1 = J2*spin(3*bc(xi+2,L)-2); sxrrghtt_2 = J2*spin(3*bc(xi+2,L)-1); sxrrghtt_3 = J2*spin(3*bc(xi+2,L))

!something's wrong here;? rectify it

        k1(3*xi-2) = sx_2*(a*(sxrght_3 + sxrrghtt_3) + b*(sxleft_3 + sxlleftt_3)) & 
                    - sx_3*(a*(sxrght_2 + sxrrghtt_2) + b*(sxleft_2 + sxlleftt_2)) 

        k1(3*xi-1) = sx_3*(a*(sxrght_1 + sxrrghtt_1) + b*(sxleft_1 + sxlleftt_1)) &
                    - sx_1*(a*(sxrght_3 + sxrrghtt_3) + b*(sxleft_3 + sxlleftt_3))

        k1(3*xi)   = sx_1*(a*(sxrght_2 + sxrrghtt_2) + b*(sxleft_2+ sxlleftt_2)) &
                    - sx_2*(a*(sxrght_1+ sxrrghtt_1) + b*(sxleft_1+ sxlleftt_1))
    end do
    
    
    !k2(1) = a*((spin(2) + 0.5*k1(2)) * (spin(6) + 0.5*k1(6))   - (spin(3) + 0.5*k1(3)) * (spin(5) + 0.5*k1(5)))
    !k2(2) = a*((spin(3) + 0.5*k1(3)) * (spin(4) + 0.5*k1(4))   - (spin(1) + 0.5*k1(1)) * (spin(6) + 0.5*k1(6))) 
    !k2(3) = a*((spin(1) + 0.5*k1(1)) * (spin(5) + 0.5*k1(5))   - (spin(2) + 0.5*k1(2)) * (spin(4) + 0.5*k1(4))) 
    
    !k2(3*L-2) = b*((spin(3*L-1) + 0.5*k1(3*L-1)) * (spin(3*(L-1)) + 0.5*k1(3*(L-1)))   &
    !             - (spin(3*L) + 0.5*k1(3*L)) * (spin(3*(L-1)-1) + 0.5*k1(3*(L-1)-1)))
                 
    !k2(3*L-1) = b*((spin(3*L) + 0.5*k1(3*L)) * (spin(3*(L-1)-2) + 0.5*k1(3*(L-1)-2))  &
    !             - (spin(3*L-2) + 0.5*k1(3*L-2)) * (spin(3*(L-1)) + 0.5*k1(3*(L-1)))) 
                 
    !k2(3*L)   = b*((spin(3*L-2) + 0.5*k1(3*L-2)) * (spin(3*(L-1)-1) + 0.5* k1(3*(L-1)-1)) &
    !             - (spin(3*L-1) + 0.5*k1(3*L-1)) * (spin(3*(L-1)-2) + 0.5*k1(3*(L-1)-2))) 
    
    do xi = 1, L
        sx_1 = spin(3*xi-2); sx_2 = spin(3*xi-1); sx_3 = spin(3*xi)

        sxleft_1 = J1*spin(3*bc(xi-1,L)-2); sxleft_2 = J1*spin(3*bc(xi-1,L)-1); sxleft_3 = J1*spin(3*bc(xi-1,L))
        sxrght_1 = J1*spin(3*bc(xi+1,L)-2); sxrght_2 = J1*spin(3*bc(xi+1,L)-1); sxrght_3 = J1*spin(3*bc(xi+1,l))

        sxlleftt_1 = J2*spin(3*bc(xi-2,L)-2); sxlleftt_2 = J2*spin(3*bc(xi-2,L)-1); sxlleftt_3 = J2*spin(3*bc(xi-2,L))
        sxrrghtt_1 = J2*spin(3*bc(xi+2,L)-2); sxrrghtt_2 = J2*spin(3*bc(xi+2,L)-1); sxrrghtt_3 = J2*spin(3*bc(xi+2,l))

        k1_1 = k1(3*xi-2); k1_2 = k1(3*xi-1); k1_3 = k1(3*xi)

        k1Left_1 = J1*k1(3*bc(xi-1,L)-2); k1Left_2 = J1*k1(3*bc(xi-1,L)-1); k1Left_3 = J1*k1(3*bc(xi-1,L))
        k1Rght_1 = J1*k1(3*bc(xi+1,L)-2); k1Rght_2 = J1*k1(3*bc(xi+1,L)-1); k1Rght_3 = J1*k1(3*bc(xi+1,L))

        k1LLeftt_1 = J2*k1(3*bc(xi-2,L)-2); k1LLeftt_2 = J2*k1(3*bc(xi-2,L)-1); k1LLeftt_3 = J2*k1(3*bc(xi-2,L))
        k1Rrghtt_1 = J2*k1(3*bc(xi+2,l)-2); k1Rrghtt_2 = J2*k1(3*bc(xi+2,l)-1); k1Rrghtt_3 = J2*k1(3*bc(xi+2,l))

        k2(3*xi-2) = (sx_2 + 0.5*k1_2) * (a*(sxrght_3 + 0.5*k1rght_3) + b*(sxleft_3 + 0.5*k1Left_3))   &
                    - (sx_3 + 0.5*k1_3)  * (a*(sxrght_2 + 0.5*k1rght_2) + b*(sxleft_2 + 0.5*k1left_2)) 

        k2(3*xi-1) = (sx_3 + 0.5*k1_3) * (a*(sxrght_1 + 0.5*k1rght_1)  + b*(sxleft_1 + 0.5*k1left_1))   &
                    -(sx_1 + 0.5*k1_1) * (a*(sxrght_3 + 0.5*k1rght_3) + b*(sxleft_3 + 0.5*k1left_3))

        k2(3*xi) = (sx_1 + 0.5*k1_1) * (a*(sxrght_2 + 0.5*k1rght_2) + b*(sxleft_2 + 0.5*k1left_2))   &
                  -(sx_2 + 0.5*k1_2) * (a*(sxrght_1 + 0.5*k1rght_1) + b*(sxleft_1 + 0.5*k1left_1))
    end do
     
    !k3(1) = a*((spin(2) + 0.5*k2(2)) * (spin(6) + 0.5*k2(6))   - (spin(3) + 0.5*k2(3)) * (spin(5) + 0.5*k2(5)))
    !k3(2) = a*((spin(3) + 0.5*k2(3)) * (spin(4) + 0.5*k2(4))   - (spin(1) + 0.5*k2(1)) * (spin(6) + 0.5*k2(6))) 
    !k3(3) = a*((spin(1) + 0.5*k2(1)) * (spin(5) + 0.5*k2(5))   - (spin(2) + 0.5*k2(2)) * (spin(4) + 0.5*k2(4))) 
    
    !k3(3*L-2) = b*((spin(3*L-1) + 0.5*k2(3*L-1)) * (spin(3*(L-1)) + 0.5*k2(3*(L-1)))   &
    !             - (spin(3*L) + 0.5*k2(3*L)) * (spin(3*(L-1)-1) + 0.5*k2(3*(L-1)-1)))
                 
    !k3(3*L-1) = b*((spin(3*L) + 0.5*k2(3*L)) * (spin(3*(L-1)-2) + 0.5*k2(3*(L-1)-2))  &
    !             - (spin(3*L-2) + 0.5*k2(3*L-2)) * (spin(3*(L-1)) + 0.5*k2(3*(L-1)))) 
                 
    !k3(3*L)   = b*((spin(3*L-2) + 0.5*k2(3*L-2)) * (spin(3*(L-1)-1) + 0.5* k2(3*(L-1)-1)) &
    !             - (spin(3*L-1) + 0.5*k2(3*L-1)) * (spin(3*(L-1)-2) + 0.5*k2(3*(L-1)-2))) 
    
    do xi = 1, L
        sx_1 = spin(3*xi-2); sx_2 = spin(3*xi-1); sx_3 = spin(3*xi)

        sxleft_1 = J1*spin(3*bc(xi-1,L)-2); sxleft_2 = J1*spin(3*bc(xi-1,L)-1); sxleft_3 = J1*spin(3*bc(xi-1,L))
        sxrght_1 = J1*spin(3*bc(xi+1,L)-2); sxrght_2 = J1*spin(3*bc(xi+1,L)-1); sxrght_3 = J1*spin(3*bc(xi+1,l))

        sxlleftt_1 = J2*spin(3*bc(xi-2,L)-2); sxlleftt_2 = J2*spin(3*bc(xi-2,L)-1); sxlleftt_3 = J2*spin(3*bc(xi-2,L))
        sxrrghtt_1 = J2*spin(3*bc(xi+2,L)-2); sxrrghtt_2 = J2*spin(3*bc(xi+2,L)-1); sxrrghtt_3 = J2*spin(3*bc(xi+2,l))

        k2_1 = k2(3*xi-2); k2_2 = k2(3*xi-1); k2_3 = k2(3*xi)

        k2Left_1 = J1*k2(3*bc(xi-1,L)-2); k2Left_2 = J1*k2(3*bc(xi-1,L)-1); k2Left_3 = J1*k2(3*bc(xi-1,L))
        k2Rght_1 = J1*k2(3*bc(xi+1,L)-2); k2Rght_2 = J1*k2(3*bc(xi+1,L)-1); k2Rght_3 = J1*k2(3*bc(xi+1,L))

        k2LLeftt_1 = J2*k2(3*bc(xi-2,L)-2); k2LLeftt_2 = J2*k2(3*bc(xi-2,L)-1); k2LLeftt_3 = J2*k2(3*bc(xi-2,L))
        k2Rrghtt_1 = J2*k2(3*bc(xi+2,l)-2); k2Rrghtt_2 = J2*k2(3*bc(xi+2,l)-1); k2Rrghtt_3 = J2*k2(3*bc(xi+2,l))

        k3(3*xi-2) = (sx_2 + 0.5*k2_2) * (a*(sxrght_3 + 0.5*k2rght_3) + b*(sxleft_3 + 0.5*k2Left_3))   &
                    - (sx_3 + 0.5*k2_3)  * (a*(sxrght_2 + 0.5*k2rght_2) + b*(sxleft_2 + 0.5*k2left_2)) 

        k3(3*xi-1) = (sx_3 + 0.5*k2_3) * (a*(sxrght_1 + 0.5*k2rght_1)  + b*(sxleft_1 + 0.5*k2left_1))   &
                    -(sx_1 + 0.5*k2_1) * (a*(sxrght_3 + 0.5*k2rght_3) + b*(sxleft_3 + 0.5*k2left_3))

        k3(3*xi) = (sx_1 + 0.5*k2_1) * (a*(sxrght_2 + 0.5*k2rght_2) + b*(sxleft_2 + 0.5*k2left_2))   &
                  -(sx_2 + 0.5*k2_2) * (a*(sxrght_1 + 0.5*k2rght_1) + b*(sxleft_1 + 0.5*k2left_1))
    end do     
    
    !k4(1) = a*((spin(2) + k3(2)) * (spin(6) + k3(6))   - (spin(3) + k3(3)) * (spin(5) + k3(5)))
    !k4(2) = a*((spin(3) + k3(3)) * (spin(4) + k3(4))   - (spin(1) + k3(1)) * (spin(6) + k3(6))) 
    !k4(3) = a*((spin(1) + k3(1)) * (spin(5) + k3(5))   - (spin(2) + k3(2)) * (spin(4) + k3(4))) 
    
    !k4(3*L-2) = b*((spin(3*L-1) + k3(3*L-1)) * (spin(3*(L-1)) + k3(3*(L-1)))   &
    !             - (spin(3*L) +  k3(3*L)) * (spin(3*(L-1)-1) + k3(3*(L-1)-1)))
                 
    !k4(3*L-1) = b*((spin(3*L) + k3(3*L)) * (spin(3*(L-1)-2) + k3(3*(L-1)-2))  &
    !             - (spin(3*L-2) + k3(3*L-2)) * (spin(3*(L-1)) + k3(3*(L-1)))) 
                 
    !k4(3*L)   = b*((spin(3*L-2) + k3(3*L-2)) * (spin(3*(L-1)-1) + k3(3*(L-1)-1)) &
    !             - (spin(3*L-1) + k3(3*L-1)) * (spin(3*(L-1)-2) + k3(3*(L-1)-2))) 
    
    do xi = 1, L  
        sx_1 = spin(3*xi-2); sx_2 = spin(3*xi-1); sx_3 = spin(3*xi)

        sxleft_1 = J1*spin(3*bc(xi-1,L)-2); sxleft_2 = J1*spin(3*bc(xi-1,L)-1); sxleft_3 = J1*spin(3*bc(xi-1,L))
        sxrght_1 = J1*spin(3*bc(xi+1,L)-2); sxrght_2 = J1*spin(3*bc(xi+1,L)-1); sxrght_3 = J1*spin(3*bc(xi+1,l))

        sxlleftt_1 = J2*spin(3*bc(xi-2,L)-2); sxlleftt_2 = J2*spin(3*bc(xi-2,L)-1); sxlleftt_3 = J2*spin(3*bc(xi-2,L))
        sxrrghtt_1 = J2*spin(3*bc(xi+2,L)-2); sxrrghtt_2 = J2*spin(3*bc(xi+2,L)-1); sxrrghtt_3 = J2*spin(3*bc(xi+2,l))

        k3_1 = k1(3*xi-2); k3_2 = k3(3*xi-1); k3_3 = k3(3*xi)

        k3Left_1 = J1*k3(3*bc(xi-1,L)-2); k3Left_2 = J1*k3(3*bc(xi-1,L)-1); k3Left_3 = J1*k3(3*bc(xi-1,L))
        k3Rght_1 = J1*k3(3*bc(xi+1,L)-2); k3Rght_2 = J1*k3(3*bc(xi+1,L)-1); k3Rght_3 = J1*k3(3*bc(xi+1,L))

        k3LLeftt_1 = J2*k3(3*bc(xi-2,L)-2); k3LLeftt_2 = J2*k3(3*bc(xi-2,L)-1); k3LLeftt_3 = J2*k3(3*bc(xi-2,L))
        k3Rrghtt_1 = J2*k3(3*bc(xi+2,l)-2); k3Rrghtt_2 = J2*k3(3*bc(xi+2,l)-1); k3Rrghtt_3 = J2*k3(3*bc(xi+2,l))

        k4(3*xi-2) = (sx_2 + k3_2) * (a*(sxrght_3 + k3rght_3) + b*(sxleft_3 + k3Left_3))   &
                    - (sx_3 + k3_3)  * (a*(sxrght_2 + k3rght_2) + b*(sxleft_2 + k3left_2)) 

        k4(3*xi-1) = (sx_3 + k3_3) * (a*(sxrght_1 + k3rght_1) + b*(sxleft_1 + k3left_1))   &
                    -(sx_1 + k3_1) * (a*(sxrght_3 + k3rght_3) + b*(sxleft_3 + k3left_3))

        k4(3*xi) = (sx_1 + k3_1) * (a*(sxrght_2 + k3rght_2) + b*(sxleft_2 + k3left_2))   &
                  -(sx_2 + k3_2) * (a*(sxrght_1 + k3rght_1) + b*(sxleft_1 + k3left_1))
    end do        
    
    do xi = 1, L
        spin(3*xi-2) = spin(3*xi-2) + (k1(3*xi-2) + 2*k2(3*xi-2) + 2*k3(3*xi-2) + k4(3*xi-2))/6.0
        spin(3*xi-1) = spin(3*xi-1) + (k1(3*xi-1) + 2*k2(3*xi-1) + 2*k3(3*xi-1) + k4(3*xi-1))/6.0
        spin(3*xi)   = spin(3*xi)   + (k1(3*xi)   + 2*k2(3*xi)   + 2*k3(3*xi)   + k4(3*xi))/6.0
    end do


End subroutine RK4dynamik_nnn
End module nnnspindynamics
