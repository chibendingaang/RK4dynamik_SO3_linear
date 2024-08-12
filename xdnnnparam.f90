module xdnnnparam 
implicit none
INTEGER, PARAMETER :: xd = SELECTED_REAL_KIND(15, 307)
save
!PUBLIC :: xd

!List of Input paramters and constants
real(xd), parameter :: PI =4.0_8*atan(1.0_xd)
integer, parameter :: L = 514
integer, parameter :: steps = 200000
real(xd), parameter :: dt = 0.005_xd

real(xd), parameter :: lambda = 0.0_xd
real(xd), parameter :: mu = 1.0_xd
!real(xd), parameter :: epsil = 0.001_xd
real(xd), parameter :: a = dt*(lambda+mu), b = dt*(lambda-mu)
real(xd), parameter :: J1 = 0.0_xd, J2 = 1.0_xd

real(xd), allocatable :: planar(:)
real(xd), allocatable :: phi(:)
real(xd), allocatable :: xyplanar(:)
real(xd) :: totalH, sitemagsq(L) !else it throws an error: sitemag has no IMPLICIT TYPE

integer :: conf !, seedz
real(xd), allocatable :: s_a(:), s0_a(:)
real(xd), allocatable :: s_b(:), s0_b(:)
real(xd), allocatable :: spin_a(:, :), spin_b(:, :)
real(xd), allocatable :: spina(:), spinb(:)
real(xd), allocatable :: k1(:), k2(:), k3(:), k4(:)

end module xdnnnparam

! https://stackoverflow.com/questions/16088892/fortran-compilation-error-undefined-reference
! https://stackoverflow.com/questions/11512897/linking-fortran-module-undefined-reference
