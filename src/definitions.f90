! Useful constants and conversion factors
module definitions
implicit none

! technical
public :: wp
integer,parameter :: wp = 8 

! scientific constants
real(kind = wp), parameter ::   pi = 3.14159265358979, &
                                kB = 1.38e-23,&               !m^2 kg / (s^2 K)
                                proton_mass = 1.6605e-27     !kg
public ::  pi, kcal_to_kJ, safeguard, boltzmann
real(kind=wp), parameter :: kcal_to_kJ = 4.184
real(kind=wp), parameter :: safeguard = 0.000000001 ! 1e-9
real(kind=wp), parameter :: boltzmann = 8.314462618e-3 ! kJ/(mol K)




end module
