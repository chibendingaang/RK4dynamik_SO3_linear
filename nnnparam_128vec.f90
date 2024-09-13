module nnnparam_128vec 
use, intrinsic :: iso_fortran_env, only: real128

implicit none
!INTEGER, PARAMETER :: xd = SELECTED_REAL_KIND(15, 307)
!save
!PUBLIC :: xd

!List of Input paramters and constants
real(real128), parameter :: PI =4.0_real128*atan(1.0_real128)
integer, parameter :: L = 2048
integer, parameter :: steps = 450000
real(real128), parameter :: dt = 0.001_real128

real(real128), parameter :: lambda = 0.0_real128
real(real128), parameter :: mu = 1.0_real128
!real(real128), parameter :: epsil = 0.001_real128
real(real128), parameter :: a = dt*(lambda+mu), b = dt*(lambda-mu)
real(real128), parameter :: J1 = 1.0_real128, J2 = 0.625_real128

real(real128), allocatable :: planar(:)
real(real128), allocatable :: phi(:)
real(real128), allocatable :: xyplanar(:)
real(real128) :: totalH, sitemagsq(L) !else it throws an error: sitemag has no IMPLICIT TYPE

integer :: conf !, seedz
real(real128), allocatable :: s_a(:,:), s0_a(:,:)
real(real128), allocatable :: s_b(:,:), s0_b(:,:)
real(real128), allocatable :: spin_a(:, :, :), spin_b(:, :, :)
real(real128), allocatable :: spina(:,:), spinb(:,:)
real(real128), allocatable :: k1(:,:), k2(:,:), k3(:,:), k4(:,:)
real(real128), allocatable :: H(:,:)

end module nnnparam_128vec

! https://stackoverflow.com/questions/16088892/fortran-compilation-error-undefined-reference
! https://stackoverflow.com/questions/11512897/linking-fortran-module-undefined-reference
