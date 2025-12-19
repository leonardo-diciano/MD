module user_settings
use definitions, only: wp
implicit none
public:: temperature, etol, ftol
real(kind = wp),parameter :: temperature=298, etol=1e-6, ftol=1e-4
  


end module user_settings
