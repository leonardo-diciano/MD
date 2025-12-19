module definitions

public :: wp

! technical
integer,parameter :: wp = 8


! scientific constants
real(kind = wp), parameter ::   pi = 3.14159265358979, &
                                kB = 1.38e-23,&               !m^2 kg / (s^2 K)
                                proton_mass = 1.6605e-27     !kg



end module
