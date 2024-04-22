module spindynamicsxd

use paramxd
use mt19937xd
!use periodic_bc

implicit none 
!Integer, parameter :: dp = selected_real_kind(15,307)

Contains

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

!SUBROUTINE#1: to initializpe spins from a random distribution

Subroutine Initperturb(spin) !, eps
!Implicit None

!real(xd), intent(in) :: eps 
real(xd), intent(out) :: spin(3*L)
real(xd) :: eps = 0.01 !intent(in)
real(xd) :: modspin
Integer :: xi
integer :: n,m
integer :: seed1 !, seed2
seed1 = conf;
n = 8 ; m = 8

!CALL RANDOM_SEED ()
!call random_number(planar)
!call random_number(phi)
print*,"seed1=",seed1
CALL init_genrand(seed1)

do xi=1,L
        phi(xi) = 2*pi*n*xi/L !0.25*pi !grnd()
        !theta(xi) = 2*pi*m*xi/L !grnd()
        !if (xi==1) then
        !planar(xi) = planar(xi) + eps
        !end if
        !phi(xi) = pi*(2*(phi(xi)) -1.0_xd)
        !planar(xi) = 2*planar(xi) - 1.0_xd
        !call random_number(planar(xi))
        !call random_number(phi(xi))
        !xyplanar(xi) = sqrt(1.0_xd - planar(xi)**2)
    spin(3*xi-2) = sin(phi(xi)) !*cos(theta(xi))
    spin(3*xi-1) = 0.0_xd !sin(phi(xi))*sin(theta(xi)) !0.0_xd 
    spin(3*xi)   = cos(phi(xi)) !1.0_xd !sqrt(1-(spin(3*xi-2))**2); 
    
    !modspin = sqrt(spin(3*xi-2)**2 + spin(3*xi-1)**2  + spin(3*xi)**2)
    !remember the sqrt!!
    !spin(3*xi-2) = spin(3*xi-2)/modspin
    !spin(3*xi-1) = spin(3*xi-1)/modspin
    !spin(3*xi) =  spin(3*xi)/modspin

end do

!spin(3) = planar(3) + eps
!spin(3) = spin(3) + eps
!intended pertrubation, added externally
end Subroutine Initperturb



!SUBROUTINE#2: for the dynamics
Subroutine RK4dynamik(spin)!harsh, replace S_ by spin_

real(xd) :: spin(3*L) !intent(inout)
real(xd) :: zp(3*L)
!Other required variables

!this could very well be the error
!Integer, external :: bc   
!f(S_x, t) = S_x \cross (a*S_{x+1} + b*S_{x-1}); a = lambda + mu; b = lambda - mu
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


!deallocate(k1(3*L)); deallocate(k2(3*L)); deallocate(k3(3*L)); deallocate(k4(3*L))

!f(S_x,zp_x,  t) = S_x \cross (a*zp_{x+1} + b*zp_{x-1}) + zp_x \cross (a*S_{x+1} + b*S_{x-1});
!g(S_x,t)  = S_x \cross (a*S_{x+1} + b*S_{x-1})
!a = dt*(lambda + mu); b = dt*(lambda - mu)

    do xi = 1, L
        k1(3*xi-2) = a*(zp(3*xi-1)  * spin(3*bc(xi+1,L))    - zp(3*xi) * spin(3*bc(xi+1,L)-1) &
        	    & + spin(3*xi-1) * zp(3*bc(xi+1,L))    - spin(3*xi) * zp(3*bc(xi+1,L)-1)) &
                    & + b*(zp(3*xi-1) * spin(3*bc(xi-1,L)) - zp(3*xi) * spin(3*bc(xi-1,L)-1) &
                    & + spin(3*xi-1) * zp(3*bc(xi-1,L)) - spin(3*xi) * zp(3*bc(xi-1,L)-1)) 

        k1(3*xi-1) = a*(zp(3*xi) * spin(3*bc(xi+1,L)-2)    - zp(3*xi-2) * spin(3*bc(xi+1,L)) &
        	    & + spin(3*xi) * zp(3*bc(xi+1,L)-2)    - spin(3*xi-2) * zp(3*bc(xi+1,L))) &
                    & + b*(zp(3*xi) * spin(3*bc(xi-1,L)-2) - zp(3*xi-2) * spin(3*bc(xi-1,L)) &
                    & + spin(3*xi) * zp(3*bc(xi-1,L)-2) - spin(3*xi-2) * zp(3*bc(xi-1,L)))

        k1(3*xi) =  a*(zp(3*xi-2) * spin(3*bc(xi+1,L) -1)    - zp(3*xi-1) * spin(3*bc(xi+1,L) -2) &
        	    & + spin(3*xi-2) * zp(3*bc(xi+1,L) -1)    - spin(3*xi-1) * zp(3*bc(xi+1,L) -2)) &
                    & +b*(zp(3*xi-2) * spin(3*bc(xi-1,L)-1) - zp(3*xi-1) * spin(3*bc(xi-1,L) -2) &
                    & + spin(3*xi-2) * zp(3*bc(xi-1,L)-1) - spin(3*xi-1) * zp(3*bc(xi-1,L) -2))
    end do
    
    do xi = 1, L
        k2(3*xi-2) = a*((zp(3*xi-1) +  k1(3*xi-1)*0.5) * spin(3*bc(xi+1,L)) &
        	    & - (zp(3*xi) + k1(3*xi)*0.5)*(spin(3*bc(xi+1,L)-1)) & 
        	    & + spin(3*xi-1)*(zp(3*bc(xi+1,L)) + 0.5*k1(3*bc(xi+1,L)))  &
        	    & - spin(3*xi) *(zp(3*bc(xi+1,L)-1) + 0.5*k1(3*bc(xi+1,L)-1))) &
                    & + b*((zp(3*xi-1) + k1(3*xi-1)*0.5)* spin(3*bc(xi-1,L)) &
                    & -(zp(3*xi) + k1(3*xi)*0.5)*(spin(3*bc(xi-1,L)-1)) &
                    & + spin(3*xi- 1)*(zp(3*bc(xi-1,L))+ 0.5*k1(3*bc(xi-1,L))) & 
                    & - spin(3*xi) *(zp(3*bc(xi-1,L)-1) + 0.5*k1(3*bc(xi-1,L)-1))) 

        k2(3*xi-1) = a*((zp(3*xi) +  k1(3*xi)*0.5) * spin(3*bc(xi+1,L)-2) & 
        	    & - (zp(3*xi-2) + k1(3*xi-2)*0.5) *(spin(3*bc(xi+1,L))) & 
        	    & + spin(3*xi)*(zp(3*bc(xi+1,L)-2) + 0.5*k1(3*bc(xi+1,L)-2)) &
        	    & - spin(3*xi-2) * (zp(3*bc(xi+1,L)) + 0.5*k1(3*bc(xi+1,L))))&
                    & + b*((zp(3*xi) + k1(3*xi)*0.5)* spin(3*bc(xi-1,L)-2) &
                    & - (zp(3*xi-2) + k1(3*xi-2)*0.5) *(spin(3*bc(xi-1,L)))  &
                    & + spin(3*xi)*(zp(3*bc(xi-1,L)-2) + 0.5*k1(3*bc(xi-1,L)-2)) &
                    & - spin(3*xi-2) * (zp(3*bc(xi-1,L)) + 0.5*k1(3*bc(xi-1,L)))) 

        k2(3*xi) = a*((zp(3*xi-2) +  k1(3*xi-2)*0.5) * spin(3*bc(xi+1,L)-1) & 
        	    & - (zp(3*xi-1)  + k1(3*xi-1)*0.5) *(spin(3*bc(xi+1,L)-2)) & 
        	    & + spin(3*xi-2)*(zp(3*bc(xi+1,L)-1) + 0.5*k1(3*bc(xi+1,L)-1)) &
        	    & - spin(3*xi-1) * (zp(3*bc(xi+1,L)-2) + 0.5*k1(3*bc(xi+1,L)-2))) &
                    & + b*((zp(3*xi-2) + k1(3*xi-2)*0.5)* spin(3*bc(xi-1,L)-1) & 
                    & - (zp(3*xi-1)+ k1(3*xi-1)*0.5) *(spin(3*bc(xi-1,L)-2))  &
                    & + spin(3*xi- 2)*(zp(3*bc(xi-1,L)-1)+ 0.5*k1(3*bc(xi-1,L)-1)) &
                    & - spin(3*xi-1) * (zp(3*bc(xi-1,L)-2) + 0.5*k1(3*bc(xi-1,L)-2))) 
     end do
     
     do xi = 1, L
        k3(3*xi-2) = a*((zp(3*xi-1) +  k2(3*xi-1)*0.5) * spin(3*bc(xi+1,L)) & 
        	    & - (zp(3*xi) +  k2(3*xi)*0.5) *(spin(3*bc(xi+1,L)-1)) & 
        	    & + spin(3*xi-1)*(zp(3*bc(xi+1,L)) + 0.5*k2(3*bc(xi+1,L))) & 
        	    & - spin(3*xi) * (zp(3*bc(xi+1,L)-1) + 0.5*k2(3*bc(xi+1,L)-1))) &
                    & + b*((zp(3*xi-1) + k2(3*xi-1)*0.5)* spin(3*bc(xi-1,L)) & 
                    & - (zp(3*xi) +  k2(3*xi)*0.5) *(spin(3*bc(xi-1,L)-1))  &
                    & + spin(3*xi- 1)*(zp(3*bc(xi-1,L))+ 0.5*k2(3*bc(xi-1,L))) &
                    & - spin(3*xi) * (zp(3*bc(xi-1,L)-1) + 0.5*k2(3*bc(xi-1,L)-1))) 

        k3(3*xi-1) = a*((zp(3*xi) +  k2(3*xi)*0.5) * spin(3*bc(xi+1,L)-2) &
        	    & - (zp(3*xi-2) + k2(3*xi-2)*0.5) *(spin(3*bc(xi+1,L))) & 
        	    & + spin(3*xi)*(zp(3*bc(xi+1,L)-2) + 0.5*k2(3*bc(xi+1,L)-2)) &
        	    & - spin(3*xi-2) * (zp(3*bc(xi+1,L)) + 0.5*k2(3*bc(xi+1,L)))) &
                    & + b*((zp(3*xi) + k2(3*xi)*0.5)* spin(3*bc(xi-1,L)-2) & 
                    & - (zp(3*xi-2) + k2(3*xi-2)*0.5) *(spin(3*bc(xi-1,L)))  &
                    & + spin(3*xi)*(zp(3*bc(xi-1,L)-2) + 0.5*k2(3*bc(xi-1,L)-2)) &
                    & - spin(3*xi-2) * (zp(3*bc(xi-1,L)) + 0.5*k2(3*bc(xi-1,L)))) 

        k3(3*xi) = a*((zp(3*xi-2) +  k2(3*xi-2)*0.5) * spin(3*bc(xi+1,L)-1) & 
        	    & - (zp(3*xi-1) + k2(3*xi-1)*0.5) *(spin(3*bc(xi+1,L)-2)) & 
        	    & + spin(3*xi-2)*(zp(3*bc(xi+1,L)-1) + 0.5*k2(3*bc(xi+1,L)-1)) &
        	    & - spin(3*xi-1) * (zp(3*bc(xi+1,L)-2) + 0.5*k2(3*bc(xi+1,L)-2))) &
                    & + b*((zp(3*xi-2) + k2(3*xi-2)*0.5)* spin(3*bc(xi-1,L)-1) & 
                    & - (zp(3*xi-1) + k2(3*xi-1)*0.5) *(spin(3*bc(xi-1,L)-2))  &
                    & + spin(3*xi- 2)*(zp(3*bc(xi-1,L)-1)+ 0.5*k2(3*bc(xi-1,L)-1)) &
                    & - spin(3*xi-1) * (zp(3*bc(xi-1,L)-2) + 0.5*k2(3*bc(xi-1,L)-2)))  
    end do     
    
    do xi = 1, L   
        k4(3*xi-2) = a*((zp(3*xi-1)  + k3(3*xi-1))*spin(3*bc(xi+1,L)) &
        	    & -  (zp(3*xi) + k3(3*xi)) * spin(3*bc(xi+1,L)-1) &
                    & + spin(3*xi-1) * (zp(3*bc(xi+1,L)) + k3(3*bc(xi+1,L))) &
                    & - spin(3*xi) * (zp(3*bc(xi+1,L)-1) + k3(3*bc(xi+1,L)-1))) &
                    & + b*((zp(3*xi-1)+ k3(3*xi-1))*spin(3*bc(xi-1,L)) &
                    & - (zp(3*xi) + k3(3*xi)) * spin(3*bc(xi-1,L)-1) &
                    & + spin(3*xi-1) * (zp(3*bc(xi-1,L)) + k3(3*bc(xi-1,L))) &
                    & - spin(3*xi) * (zp(3*bc(xi-1,L)-1) + k3(3*bc(xi-1,L)-1)))
                    
        k4(3*xi-1) = a*((zp(3*xi)  + k3(3*xi))*spin(3*bc(xi+1,L)-2) &
        	    & -  (zp(3*xi-2) + k3(3*xi-2)) * spin(3*bc(xi+1,L)) &
                    & + spin(3*xi) * (zp(3*bc(xi+1,L)-2) + k3(3*bc(xi+1,L)-2)) &
                    & - spin(3*xi-2) * (zp(3*bc(xi+1,L)) + k3(3*bc(xi+1,L)))) &
                    & + b*((zp(3*xi)+ k3(3*xi))*spin(3*bc(xi-1,L)-2) & 
                    & - (zp(3*xi-2) + k3(3*xi-2)) * spin(3*bc(xi-1,L)) &
                    & + spin(3*xi) * (zp(3*bc(xi-1,L)-2) + k3(3*bc(xi-1,L)-2)) &
                    & - spin(3*xi-2) * (zp(3*bc(xi-1,L)) + k3(3*bc(xi-1,L))))

        k4(3*xi) = a*((zp(3*xi-2)  + k3(3*xi-2))*spin(3*bc(xi+1,L)-1) &
        	    & - (zp(3*xi-1) + k3(3*xi-1)) * spin(3*bc(xi+1,L)-2) &
                    & + spin(3*xi-2) * (zp(3*bc(xi+1,L)-1) + k3(3*bc(xi+1,L)-1)) &
                    & - spin(3*xi-2) * (zp(3*bc(xi+1,L)-1) + k3(3*bc(xi+1,L)-1))) &
                    & + b*((zp(3*xi-2) + k3(3*xi-2)) * spin(3*bc(xi-1,L)-1) &
                    & -  (zp(3*xi-1) +  k3(3*xi-1))  * spin(3*bc(xi-1,L)-2) &
                    & + spin(3*xi-2) * (zp(3*bc(xi-1,L)-1) + k3(3*bc(xi-1,L)-1)) &
                    & - spin(3*xi-2) * (zp(3*bc(xi-1,L)-1) + k3(3*bc(xi-1,L)-1)))
    end do        
    
    do xi = 1, L
        zp(3*xi-2) = zp(3*xi-2) + (k1(3*xi-2) + 2*k2(3*xi-2) + 2*k3(3*xi-2) + k4(3*xi-2))/6.0
        zp(3*xi-1) = zp(3*xi-1) + (k1(3*xi-1) + 2*k2(3*xi-1) + 2*k3(3*xi-1) + k4(3*xi-1))/6.0
        zp(3*xi)   = zp(3*xi)   + (k1(3*xi)   + 2*k2(3*xi)   + 2*k3(3*xi)   + k4(3*xi))/6.0
    end do 

End subroutine RK4dynamik

End module spindynamicsxd
