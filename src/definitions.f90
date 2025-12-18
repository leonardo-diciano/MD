! Useful constants and conversion factors
module definitions
implicit none

public :: wp, pi, kcal_to_kJ, safeguard, boltzmann
integer,parameter :: wp = 8 
real(kind=wp), parameter :: pi = 3.14159265358979
real(kind=wp), parameter :: kcal_to_kJ = 4.184
real(kind=wp), parameter :: safeguard = 0.000000001 ! 1e-9
real(kind=wp), parameter :: boltzmann = 8.314462618e-3 ! kJ/(mol K)



end module

