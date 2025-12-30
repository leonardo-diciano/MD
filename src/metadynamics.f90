! Module to handle metadynamics
!


module metadynamics_mod
implicit none

contains

subroutine metadynamics_propagation(positions,xyzfile)
! Velocity Verlet Propagation scheme with integrated well-tempered metadynamics
use definitions, only: wp
use simulation_subroutines, only: init_v
use force_field_mod, only: get_energy_gradient, mweights
use parser, only: meta_cv, meta_tau, meta_nsteps, meta_cv_type

implicit none
real(kind=wp),intent(inout) :: positions(n_atoms,3)
character(len=256), intent(in) :: xyzfile
integer :: meta_counter = 1
real(kind=wp) :: bias_param(meta_nsteps,2), tot_bias_pot, cv_value

real(kind=wp) :: displacement(n_atoms), positions_previous(n_atoms,3), input_positions(n_atoms,3), forces(n_atoms,3), &
                total_displacement(n_atoms), velocities(n_atoms,3)
real(kind=wp) :: positions_list(2,n_atoms,3), acceleration_list(2,n_atoms,3)
real(kind=wp) :: gradnorm, tot_pot, instant_temp, pressure, E_kin, tot_momentum(3), tot_momentum_norm
character(len=256) :: traj_xyzfile, properties_outfile
integer :: istep, icartesian,i, dot, current, previous, new
logical :: suppress_flag = .true.


! for intuitive storing of position data 
current = 1
new = 2

! INITIALIZATION
istep = 0
total_displacement(:) = 0
input_positions(:,:) = positions(:,:) !x(t=0)
positions_list(current,:,:) = positions(:,:) !x(t=0)
acceleration_list(:,:,:) = 0
bias_param(:,:) = 0.0
tot_bias_pot = 0.0
tot_pot = 0.0
forces(:,:) = 0.0

if (meta_cv_type == "distance") then
    CALL distance_bias(positions,tot_pot,forces,init_cv_value)
elseif (meta_cv_type == "angle") then
    CALL angle_bias(positions,tot_pot,forces,init_cv_value)
elseif (meta_cv_type == "dihedral") then
    CALL dihedral_bias(positions,tot_pot,forces,init_cv_value)
else
    write(*,*) "Error in CV type"
    stop
end if

call get_energy_gradient(positions_list(current,:,:),tot_pot,forces, gradnorm, suppress_flag)
call init_v(velocities(:,:)) 

do icartesian = 1,3
    acceleration_list(current,:,icartesian) = 1e-4 * forces(:,icartesian) / mweights(:)
    !acceleration in Å/fs^2                    in kJ/mol/Å           in g/mol;
end do

! PREPARE FILE THAT TRACKS PROPERTIES
dot = index(xyzfile, ".", back=.true.)     ! find last "." in xyzfile name
properties_outfile = xyzfile(:dot-1) // ".properties.txt"
open(97, file=properties_outfile, status='replace', action='write')
write(97,"(A10,9(A20))") "istep", "E_tot", "E_kin","E_pot", "F_norm", "Bias", "Free_energy" , "Temp", "Pressure", "COM_momentum"
write(97,"(A10,9(A20))") "none","kJ/mol", "kJ/mol","kJ/mol", "kJ/(Åmol)", "kJ/mol","kJ/mol","K", "Pa", "gÅ/fs"

! PREPARE TRAJECTORY FILE: traj_xyzfile
traj_xyzfile = xyzfile(:dot-1) // "_meta.traj" // xyzfile(dot:)
open(98, file=traj_xyzfile, status='replace', action='write')
write(98,*) n_atoms
write(98,"(A,F6.2,A)") "atomic positions at t = ",istep * md_ts, " fs"
do i=1, size(positions,1), 1
    write(98,FMT='(A3,3(2X,F15.8))') atomnames(i), positions(i,1:3)
end do

 do while (istep<md_nsteps)
        istep = istep +1

    ! update position with velocity verlet
    CALL velocity_verlet_position(positions_list(current,:,:), velocities(:,:),acceleration_list(current,:,:), &
                            positions_list(new,:,:))

    
    CALL get_energy_gradient(positions,tot_pot,forces, gradnorm, suppress_flag)

    ! Deposit a new gaussian bias potential every meta_tau fs
    if (MOD((md_ts * i),meta_tau) == 0) then
        CALL calc_cv(positions_list(current,:,:))
        CALL deposit_bias_potential(meta_counter,cv_value,init_cv_value,bias_param,tot_bias_pot)
    end if

    if (meta_cv_type == "distance") then
        CALL distance_bias(positions_list(),tot_pot,forces,cv_value)
    elseif (meta_cv_type == "angle") then
        CALL angle_bias(positions_list(),tot_pot,forces,cv_value)
    elseif (meta_cv_type == "dihedral") then
        CALL dihedral_bias(positions_list(),tot_pot,forces,cv_value)
    else
        write(*,*) "Error in CV type"
        stop
    end if
    

    do icartesian = 1,3
        acceleration_list(new,:,icartesian) = 1e-4 * forces(:,icartesian) / mweights(:)
        !acceleration in Å/fs^2                    in kJ/mol/Å           in g/mol;
    end do
    
    CALL velocity_verlet_velocity(old_acceleration, new_acceleration, velocities)




end do
end subroutine

subroutine deposit_bias_potential(meta_counter,cv_value,init_cv_value,bias_param,tot_bias_pot)
use definitions, only: wp
use force_field_mod, only: n_atoms
use parser, only: meta_cv, meta_dT, meta_omega, meta_sigma, meta_nsteps 

implicit none 

real(kind=wp), intent(inout) :: tot_bias_pot, bias_param(meta_nsteps,2),init_cv_value 
integer, intent(inout) :: meta_counter
real(kind=wp), intent(in) :: cv_value, positions(n_atoms,3)
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
init_cv_value = cv_value
end subroutine

function calc_cv(positions) result(cv_value)
use definitions, only: wp
use force_field_mod, only: n_atoms
use parser, only: meta_cv, meta_cv_type

implicit none
real(kind=wp), intent(in) :: positions(n_atoms,3)
real(kind=wp), intent(out) :: cv_value
real(kind=wp) :: distance, d12(3), d23(3), d12_norm, d23_norm, d34(3), a(3), b(3), a_norm, b_norm, ratio

if (meta_cv_type == "distance") then
    distance = positions(meta_cv(2),:) - positions(meta_cv(1),:)
    cv_value = SQRT(dot_product(distance,distance))
elseif (meta_cv_type == "angle") then
    d12 = positions(meta_cv(1),1:3) - positions(meta_cv(2),1:3)
    d23 = positions(meta_cv(3),1:3) - positions(meta_cv(2),1:3)
    d12_norm = SQRT(dot_product(d12,d12))
    d23_norm = SQRT(dot_product(d23,d23))
    cv_value = acos(dot_product(d12,d23) / (d12_norm * d23_norm))
elseif (meta_cv_type == "dihedral") then
    d12 = positions(meta_cv(2),1:3) - positions(meta_cv(1),1:3)
    d23 = positions(meta_cv(3),1:3) - positions(meta_cv(2),1:3)
    d34 = positions(meta_cv(4),1:3) - positions(meta_cv(3),1:3)
    a = cross_product(d12,d23)
    a_norm = SQRT(dot_product(a,a))
    b = cross_product(d23,d34)
    b_norm = SQRT(dot_product(b,b))
    ratio= dot_product(a,b) / (a_norm * b_norm)
    if (ratio > 1) then
        if (ratio - 1 > 0.1) then   
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
    cv_value = acos(ratio)
else
    stop
end if
end function

subroutine distance_bias(bias_param,positions,tot_pot,forces,cv_value)
use definitions, only: wp
use force_field_mod, only: n_atoms
use parser, only: meta_cv, meta_nsteps, meta_cv_type, meta_sigma

implicit none

real(kind=wp), intent(in) :: bias_param(meta_nsteps,2), positions(n_atoms,3)
real(kind=wp), intent(inout) :: tot_pot,forces(n_atoms,3),cv_value 
real(kind=wp) :: distance, der_1(3), der_2(3)

distance = positions(meta_cv(2),:) - positions(meta_cv(1),:)
cv_value = SQRT(dot_product(distance,distance))

do i=1, meta_nsteps
    tot_pot = tot_pot + bias_param(i,1) * ( exp(-0.5*((cv_value-bias_param(i,2)) / meta_sigma)**2))

    der_1= - bias_param(i,1) * ((cv_value - bias_param(i,2))/meta_sigma) * &
            exp(-0.5*((cv_value-bias_param(i,2)) / meta_sigma)**2) * &
            (positions(meta_cv(2),:) - positions(meta_cv(1),:)) / cv_value

    der_2 = bias_param(i,1) * ((cv_value - bias_param(i,2))/meta_sigma) * &
            exp(-0.5*((cv_value-bias_param(i,2)) / meta_sigma)**2) * &
            (positions(meta_cv(2),:) - positions(meta_cv(1),:)) / cv_value

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
real(kind=wp), intent(inout) :: positions(n_atoms,3),tot_pot,forces(n_atoms,3), cv_value 
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
real(kind=wp), intent(inout) :: positions(n_atoms,3),tot_pot,forces(n_atoms,3), cv_value
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

function calc_free_energy()

implicit none

! calculate probability histogram and reweight 


end function

end module