module nnnspindynamik_128vec

use nnnparam_128vec
use mt19937_128
! use periodic_bc

implicit none 
! Integer, parameter :: dp = selected_real_kind(15,307)

Contains

function compute_f(S) result(f)
    real(real128), intent(in) :: S(3,L)
    real(real128) :: f(3,L), H(3,L)
    integer :: i

    ! Computes effective field H from neighbouring contributions
    H = 0.0_real128
    !do i = 1, L
     H = mu * (J1 * (cshift(S, 1, 2) - cshift(S, -1, 2)) + &
                          J2 * (cshift(S, 2, 2) - cshift(S, -2, 2)))
    !end do
    ! here taking either mu or lambda is irrelevant as long as the sign tacked
    ! on to the left/right neighbours dictates the dyanmics

    ! Compute cross product S × H
    f(1,:) = S(2,:) * H(3,:) - S(3,:) * H(2,:)
    f(2,:) = S(3,:) * H(1,:) - S(1,:) * H(3,:)
    f(3,:) = S(1,:) * H(2,:) - S(2,:) * H(1,:)
  end function compute_f

!FUNCTION to define periodic boundary conditions
!Integer function bc(x,L)
!    Integer, intent(in) :: x, L
    
!    if (x==L+1) then
!        bc =1
!    else if (x==0) then
!        bc = L
!    else
!        bc = x
!    end if

!End function bc

!45 continue


! SUBROUTINE#1: to initialize spins from a random distribution
Subroutine InitRandom(spin, eps) 

!Implicit None
! double precision, intent(in) :: eps 
real(real128), intent(out) :: spin(3,L)
real(real128), intent(in) :: eps
Integer :: xi
integer :: seed_val !, seed2
seed_val = serialkey_value;

! CALL RANDOM_SEED ()
! call random_number(planar)
! call random_number(phi)
print*,"seed key =",seed_val
CALL init_genrand(seed_val)

do xi=1,L
        planar(xi) = grnd()
        phi(xi) = grnd()
        if (xi==L/2) then
        planar(xi) = planar(xi) + eps
        end if
        phi(xi) = pi*(2*(phi(xi)) - 1.0_real128)
        planar(xi) = 2*planar(xi) - 1.0_real128
        xyplanar(xi) = sqrt(1.0_real128 - planar(xi)**2)
    spin(1,xi) = xyplanar(xi)*cos(phi(xi))
    spin(2,xi) = xyplanar(xi)*sin(phi(xi))
    spin(3,xi) = planar(xi)
    
end do

! if x = 2, L-1; then the following spin values are required:
! spin(1) = spin(4); spin(2) = spin(5); spin(3) = spin(6)
! spin(3*L-2) = spin(3*L-5); spin(3*L-1) = spin(3*L-4); spin(3*L) = spin(3*L-3)

end Subroutine InitRandom



! SUBROUTINE#2: for the dynamics
Subroutine RK4dynamik_nnn(spin)

real(real128), intent(inout) :: spin(3,L) !intent(inout)
real(real128) :: k1(3,L), k2(3,L), k3(3,L), k4(3,L), S_temp(3,L)


k1 = dt * compute_f(spin)
S_temp = spin + 0.5_real128 * k1
k2 = dt * compute_f(S_temp)
S_temp = spin + 0.5_real128 * k2
k3 = dt * compute_f(S_temp)
S_temp = spin + k3
k4 = dt * compute_f(S_temp)
  
spin = spin + (k1 + 2.0_real128*k2 + 2.0_real128*k3 + k4) / 6.0_real128

End subroutine RK4dynamik_nnn
End module nnnspindynamik_128vec
