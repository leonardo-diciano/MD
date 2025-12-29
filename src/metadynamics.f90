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
    
    CALL velocity_verlet_position(positions_current, velocities, acceleration, positions_new)
    CALL velocity_verlet_velocity(old_acceleration, new_acceleration, velocities)

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
real(kind=wp) :: der_1(3), der_2(3), distance

distance = positions(meta_cv(2),:) - positions(meta_cv(1),:)
cv_value = SQRT(dot_product(distance,distance))

do i=1, meta_nsteps
    tot_pot = tot_pot + bias_param(i,1) * ( exp(-0.5*((cv_value-bias_param(i,2)) / meta_sigma)**2))

    der_1= - bias_param(i,1) * ((cv_value - bias_param(i,2))/meta_sigma) * &
            exp(-0.5*((cv_value-bias_param(i,2)) / meta_sigma)**2) * &
            (positions(meta_cv(2),:) - positions(meta_cv(1),:)) / distance

    der_2 = bias_param(i,1) * ((cv_value - bias_param(i,2))/meta_sigma) * &
            exp(-0.5*((cv_value-bias_param(i,2)) / meta_sigma)**2) * &
            (positions(meta_cv(2),:) - positions(meta_cv(1),:)) / distance

    force(meta_cv(1),:) = force(meta_cv(1),:) + der_1
    force(meta_cv(2),:) = force(meta_cv(2),:) + der_2
end do

end subroutine

subroutine angle_bias(bias_param,positions,tot_pot,forces,cv_value)
use definitions, only: wp
use force_field_mod, only: n_atoms
use parser, only: meta_cv, meta_nsteps, meta_cv_type, meta_sigma

implicit none

real(kind=wp), intent(in) :: bias_param(meta_nsteps,2)
real(kind=wp), intent(inout) :: positions(n_atoms,3),tot_pot,forces(n_atoms,3)
real(kind=wp), intent(out) :: cv_value 
real(kind=wp) :: d12(3), d23(3), d12_norm, d23_norm, der_magn, der_1(3), der_3(3)


d12 = positions(meta_cv(1),1:3) - positions(meta_cv(2),1:3)
d23 = positions(meta_cv(3),1:3) - positions(meta_cv(2),1:3)
d12_norm = SQRT(dot_product(d12,d12))
d23_norm = SQRT(dot_product(d23,d23))

! Calculate the angle (in radians)
cv_value = acos(dot_product(d12,d23) / (d12_norm * d23_norm))

do i=1, meta_nsteps
    tot_pot = tot_pot + bias_param(i,1) * ( exp(-0.5*((cv_value-bias_param(i,2)) / meta_sigma)**2))

    der_magn = bias_param(i,1) * ((cv_value - bias_param(i,2))/meta_sigma) * &
            exp(-0.5*((cv_value-bias_param(i,2)) / meta_sigma)**2)

    der_1 = - der_magn * (1 / ( sin(cv_value) * d12_norm + safeguard)) * &
                    ( cos(cv_value) * (d12 / d12_norm) - (d23 / d23_norm ))

    der_3 = - der_magn * (1 / ( sin(cv_value) * d23_norm + safeguard)) * &
                    (cos(cv_value) * (d23 / d23_norm) - (d12 / d12_norm))

    forces(meta_cv(1),1:3) = forces(meta_cv(1),1:3) + der_1
    forces(meta_cv(2),1:3) = forces(meta_cv(2),1:3) - (der_1 + der_3)
    forces(meta_cv(3),1:3) = forces(meta_cv(3),1:3) + der_3
end do

end subroutine

subroutine dihedral_bias(bias_param,positions,tot_pot,forces,cv_value)
use definitions, only: wp
use force_field_mod, only: n_atoms
use parser, only: meta_cv, meta_nsteps, meta_cv_type, meta_sigma

implicit none

real(kind=wp), intent(in) :: bias_param(meta_nsteps,2)
real(kind=wp), intent(inout) :: positions(n_atoms,3),tot_pot,forces(n_atoms,3)
real(kind=wp), intent(out) :: cv_value 
real(kind=wp) :: d12(3), d23(3), d34(3), a(3), b(3), a_norm, b_norm, der_magn, &
                ratio, cap_A, cap_B, der_1(3), der_2(3), der_3(3), der_4(3)

d12 = positions(meta_cv(2),1:3) - positions(meta_cv(1),1:3)
d23 = positions(meta_cv(3),1:3) - positions(meta_cv(2),1:3)
d34 = positions(meta_cv(4),1:3) - positions(meta_cv(3),1:3)

! Define a (and b) as cross product and its norm a_norm (and b_norm)
a = cross_product(d12,d23)
a_norm = SQRT(dot_product(a,a))
b = cross_product(d23,d34)
b_norm = SQRT(dot_product(b,b))

ratio= dot_product(a,b) / (a_norm * b_norm)
if (ratio > 1) then
    if (ratio - 1 > 0.1) then   ! check that the difference is not too high
        stop
    else
        ratio = 1.0
    endif
elseif (ratio < -1) then
    if (ratio + 1 < -0.1) then
        stop
    else
        ratio = -1.0
    endif
else
    ratio = ratio
endif
! Calculate the dihedral angle
cv_value = acos(ratio)

! Define A and B terms
cap_A = (b / (a_norm * b_norm)) - ((cos(cv_value) * a) / (a_norm**2))
cap_B = (a / (a_norm * b_norm)) - ((cos(cv_value) * b) / (b_norm**2 ))

do i=1, meta_nsteps
    tot_pot = tot_pot + bias_param(i,1) * ( exp(-0.5*((cv_value-bias_param(i,2)) / meta_sigma)**2))

    der_magn = bias_param(i,1) * ((cv_value - bias_param(i,2))/meta_sigma) * &
            exp(-0.5*((cv_value-bias_param(i,2)) / meta_sigma)**2) 

    der_1 = der_magn * (1/ (sin(cv_value) + safeguard)) * cross_product(cap_A,d23)
    der_4 = - der_magn * (1/ (sin(cv_value) + safeguard)) * cross_product(d23,cap_B)
    der_2 = - der_1 + (der_magn * (1/ (sin(cv_value) + safeguard)) * (cross_product(d12,cap_A) + cross_product(cap_B,d34)))
    der_3 = - der_4 - (der_magn * (1/ (sin(cv_value) + safeguard)) * (cross_product(cap_A,d12) + cross_product(cap_B,d34)))

    forces(meta_cv(1),1:3) = forces(meta_cv(1),1:3) + der_1
    forces(meta_cv(2),1:3) = forces(meta_cv(2),1:3) + der_2 
    forces(meta_cv(3),1:3) = forces(meta_cv(3),1:3) + der_3
    forces(meta_cv(4),1:3) = forces(meta_cv(4),1:3) + der_4
end do

end subroutine



end module