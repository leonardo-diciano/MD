! Module to handle metadynamics
!


module metadynamics_mod
implicit none

contains

subroutine deposit_bias_potential(meta_counter,cv_value,bias_param,tot_bias_pot)
use definitions, only: wp
use force_field_mod, only: n_atoms
use parser, only: meta_cv, meta_dT, meta_omega, meta_sigma, meta_nsteps 

implicit none 

real(kind=wp), intent(inout) :: tot_bias_pot, bias_param(meta_nsteps,2) 
integer, intent(inout) :: meta_counter
real(kind=wp), intent(in) :: cv_value, positions(n_atoms,3)
real(kind=wp) :: gaussian_potential, height

! calculate the height based on the previous bias potentials and the constants
height = meta_omega * meta_tau * exp(-(tot_bias_pot/meta_dT))

! Calculate the new gaussian potential 
gaussian_potential = exp(-0.5*((cv_value-bias_param(meta_counter - 1,2)) / meta_sigma)**2)

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
cv_value = 0.0


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
        CALL deposit_bias_potential(meta_counter,cv_value,bias_param,tot_bias_pot)
    end if



end do
end subroutine


subroutine distance_bias(positions,tot_pot,forces)
use definitions, only: wp
use force_field_mod, only: n_atoms
use parser, only: meta_cv, meta_tau, meta_nsteps, meta_cv_type

implicit none

real(kind=wp), intent(inout) :: positions(n_atoms,3),tot_pot,forces(n_atoms,3)
real(kind=wp), intent()


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