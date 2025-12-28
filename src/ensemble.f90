! Thermostat and Barostat module with subroutines for:
!   - Bussi/CSVR thermostat
!   - Berendsen barostat
! Author: Leonardo Di Ciano (2025)

module ensemble_mod

contains

subroutine bussi_thermostat(new_K,ndeg,alpha)
use definitions, only: wp,boltzmann
use lin_alg, only: gauss_distrib, sumnoises
use parser_mod, only: md_temp, bus_tau, md_ts

implicit none
    
real(kind=wp), intent(in) :: new_K
integer, intent(in) :: ndeg
real(kind=wp), intent(out) :: alpha
real(kind=wp) :: rr, factor, targ_K

! Calculation of the target kinetic energy by equipartion theorem 
targ_K = (ndeg / 2) * boltzmann * md_temp

! check that thermostat relaxation time is larger than timestep
if (bus_tau > md_ts) then
    factor=exp(- md_ts / bus_tau)
else
    factor = 0
end if
! Extraction of a random number from a gaussian distribution 
rr = gauss_distrib()

! Calculation of alpha (rescaling factor)
! sumnoises function calculates the sum of ndeg-1 independent gaussian numbers squared
alpha = SQRT( factor + (targ_K/(ndeg*new_K)) * (1-factor) * (rr**2 + sumnoises(ndeg-1)) +&
         2 * SQRT(factor) * SQRT((targ_K / (ndeg *new_K)) * (1-factor)*rr) )
end subroutine



subroutine berendsen_barostat(positions,new_P)
use definitions, only: wp
use force_field_mod, only: n_atoms
use parser_mod, only: md_ts, md_press,ber_tau,ber_k

implicit none
    
real(kind=wp), intent(in) :: new_P
real(kind=wp), intent(inout) :: positions (n_atoms,3)
real(kind=wp) :: lambda

! Calculate the rescaling of the cell axis given by different P
lambda = ( 1 + ( ((ber_k * md_ts)/ber_tau) * (new_P - md_press)))**(1/3)

! Rescale coordinates to a new "box volume"
positions = lambda * positions

end subroutine



end module