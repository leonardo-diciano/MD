! Module to handle metadynamics
!


module metadynamics_mod
implicit none

contains

subroutine deposit_bias_potential(meta_counter,cv_value,init_cv_value,bias_param,tot_bias_pot)
use definitions, only: wp
use force_field_mod, only: n_atoms
use parser, only: meta_cv, meta_dT, meta_omega, meta_sigma, meta_nsteps 

implicit none 

real(kind=wp), intent(inout) :: tot_bias_pot, bias_param(meta_nsteps,2) 
integer, intent(inout) :: meta_counter
real(kind=wp), intent(in) :: cv_value,init_cv_value, positions(n_atoms,3)
real(kind=wp) :: gaussian_potential, height

! calculate the height based on the previous bias potentials and the constants
height = meta_omega * meta_tau * exp(-(tot_bias_pot/meta_dT))

! Calculate the new gaussian potential 
if (meta_counter > 1) then
    gaussian_potential = exp(-0.5*((cv_value-bias_param(meta_counter - 1,2)) / meta_sigma)**2)
else
    gaussian_potential = exp(-0.5*((cv_value-init_cv_value) / meta_sigma)**2)
end if

! Calculate the total bias potential applied at this deposition time
tot_bias_pot = tot_bias_pot + (height * gaussian_potential)

! Store params for evaluation of bias at each MD timestep
bias_param(meta_counter,1) = height
bias_param(meta_counter,2) = cv_value

! Update the counter
meta_counter = meta_counter + 1
end subroutine


subroutine metadynamics_propagation()
use definitions, only: wp
use parser, only: meta_cv, meta_tau, meta_nsteps, meta_cv_type

implicit none
integer :: meta_counter = 1
real(kind=wp) :: bias_param(meta_nsteps,2), tot_bias_pot, cv_value

bias_param(:,:) = 0.0
tot_bias_pot = 0.0
init_cv_value = 


! Velocity Verlet Propagation scheme with integrated well-tempered metadynamics
do i=1, md_nstepsnstep
    
    !calculate forces like usual
    call get_energy_gradient(positions,tot_pot,forces, gradnorm, suppress_flag)

    if (meta_cv_type == "distance") then
        CALL distance_bias(positions,tot_pot,forces,cv_value)
    elseif (meta_cv_type == "angle") then
        CALL angle_bias(positions,tot_pot,forces,cv_value)
    elseif (meta_cv_type == "dihedral") then
        CALL dihedral_bias(positions,tot_pot,forces,cv_value)
    else
        write(*,*) "Error in CV type"
        stop
    end if
    
    CALL velocity_verlet(positions,velocities,acceleration)

    ! Deposit a new gaussian bias potential every meta_tau fs
    if (MOD((md_ts * i),meta_tau) == 0) then
        CALL deposit_bias_potential(meta_counter,cv_value,init_cv_value,bias_param,tot_bias_pot)
    end if



end do
end subroutine


subroutine distance_bias(bias_param,positions,tot_pot,forces,cv_value)
use definitions, only: wp
use force_field_mod, only: n_atoms
use parser, only: meta_cv, meta_nsteps, meta_cv_type, meta_sigma

implicit none

real(kind=wp), intent(in) :: bias_param(meta_nsteps,2)
real(kind=wp), intent(inout) :: positions(n_atoms,3),tot_pot,forces(n_atoms,3)
real(kind=wp), intent(out) :: cv_value 
real(kind=wp) :: derivative(3)

distance = positions(meta_cv(2),:) - positions(meta_cv(1),:)
cv_value = SQRT(dot_product(distance,distance))

do i=1, meta_nsteps
    tot_pot = tot_pot + bias_param(i,1) * ( exp(-0.5*((cv_value-bias_param(meta_counter,2)) / meta_sigma)**2))

    der_1= - bias_param(i,1) * ((cv_value - bias_param(meta_counter,2))/meta_sigma) * &
            exp(-0.5*((cv_value-bias_param(meta_counter,2)) / meta_sigma)**2) * &
            (positions(meta_cv(2),:) - positions(meta_cv(1),:)) / distance

    der_2 = bias_param(i,1) * ((cv_value - bias_param(meta_counter,2))/meta_sigma) * &
            exp(-0.5*((cv_value-bias_param(meta_counter,2)) / meta_sigma)**2) * &
            (positions(meta_cv(2),:) - positions(meta_cv(1),:)) / distance

    force(meta_cv(1),:) = force(meta_cv(1),:) + der_1
    force(meta_cv(2),:) = force(meta_cv(2),:) + der_2
end do

end subroutine

subroutine angle_bias()
implicit none

end subroutine

subroutine dihedral_bias()
implicit none

end subroutine

function distance_CV(positions) result(cv_value)

end function

end module