! Useful constants and conversion factors
module definitions
implicit none

! technical
public :: wp
integer,parameter :: wp = 8


! FOR PERIODIC BOUNDARY CONDITIONS
public :: md_pbc, pbc_debug, md_boxlength
logical :: md_pbc = .true., pbc_debug = .true.
real(kind=wp) :: md_boxlength = 10

! scientific constants
public ::  pi, kcal_to_kJ, safeguard, boltzmann, kB, proton_mass,avogad
real(kind = wp), parameter ::   pi = 3.14159265358979
real(kind = wp), parameter ::   kB = 1.38e-23               !m^2 kg / (s^2 K)
real(kind = wp), parameter ::   proton_mass = 1.6605e-27     !kg
real(kind = wp), parameter ::   avogad = 6.02214076 * 1e23  ! 1/mol
real(kind=wp), parameter :: kcal_to_kJ = 4.184
real(kind=wp), parameter :: safeguard = 0.000000001 ! 1e-9
real(kind=wp), parameter :: boltzmann = 8.314462618e-3 ! kJ/(mol K)




end module
