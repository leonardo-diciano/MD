! Thermostat and Barostat module with subroutines for:
!   - Bussi/CSVR thermostat
!   - Berendsen barostat
! Author: Leonardo Di Ciano (2025)

module ensemble_mod

contains

subroutine bussi_thermostat(new_K,ndeg, tau, alpha, Dt)
use definitions, only: wp,boltzmann
use lin_alg, only: gauss_distrib, sumnoises
use parser_mod, only: md_temp

implicit none
    
real(kind=wp), intent(in) :: new_K, tau, Dt
integer, intent(in) :: ndeg
real(kind=wp), intent(out) :: alpha
real(kind=wp) :: rr, factor, targ_K

targ_K = (ndeg / 2) * boltzmann * md_temp

if (tau > Dt) then
    factor=exp(- Dt / tau)
else
    factor = 0
end if
rr = gauss_distrib()
alpha = SQRT( factor + (targ_K/(ndeg*new_K)) * (1-factor) * (rr**2 + sumnoises(ndeg-1)) +&
         2 * SQRT(factor) * SQRT((targ_K / (ndeg *new_K)) * (1-factor)*rr) )
end subroutine



subroutine berendsen_barostat(positions,new_P,targ_P,tau,Dt,k)
use definitions, only: wp

implicit none
    
real(kind=wp), intent(in) :: new_P, targ_P, tau, Dt
real(kind=wp), allocatable, intent(inout) :: positions (:,:)
real(kind=wp) :: lambda, k

! Just calculate the rescaling of the cell given by different P
lambda = ( 1 - ( ((k * Dt)/tau) * (new_P - targ_P)))**(-3)

! Rescale coordinates to a new "box volume"
positions = lambda * positions

end subroutine



end module