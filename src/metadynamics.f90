! Module for well-tempered metadynamics 
! simulations with subroutines to add 
! the bias potential and perform the 
! simulation with velocity Verlet
!
! Author: Leonardo Di Ciano (2025)

module metadynamics_mod
implicit none

contains

subroutine run_metadynamics(positions,xyzfile)
! The subroutine handles the initialization of a metadynamics simulation
use definitions, only: wp
use force_field_mod, only: n_atoms
use parser_mod, only: meta_cv, meta_tau, meta_nsteps, meta_cv_type, meta_dT, meta_omega, meta_sigma, &
                     md_nsteps, md_ts, md_ensemble, md_fix_com_mom, md_temp, md_press, bus_tau, ber_k,&
                     ber_tau
implicit none

real(kind=wp), intent(inout) :: positions(n_atoms,3)
character(len=256), intent(in) :: xyzfile

! Some security checks and adjustments before actually running the metadynamics
if (meta_tau < md_ts) then
    md_ts = meta_tau 
end if 

if ((meta_nsteps * meta_tau) > (md_nsteps * md_ts)) then
    md_nsteps = meta_nsteps + 10
end if

! Print settings of the metadynamics
write(*,*)" "
write(*,*) "METADYNAMICS module"
write(*,*)
write(*,*) "Settings:"
write(*,"(A22,I5)") "  MD number of steps: ", md_nsteps
write(*,"(A17,A16)") "  MD integrator: ", " Velocity Verlet"
write(*,"(A15,F10.3,A3)") "  MD timestep: ", md_ts, " fs"
write(*,"(A15,A3)") "  MD ensemble: ", md_ensemble
write(*,"(A18,F10.3,A2)") "  MD temperature: ", md_temp, " K"
write(*,"(A15,F10.3,A3)") "  MD pressure: ", md_press, " Pa"
if (md_ensemble == "NVT") then
    write(*,"(A20,A3)") "  Bussi thermostat: ", " ON"
    write(*,"(A31,F10.3)") "    Bussi time constant (tau): ", bus_tau
elseif (md_ensemble == "NPT") then
    write(*,"(A20,A3)") "  Bussi thermostat: ", " ON"
    write(*,"(A31,F10.3)") "    Bussi time constant (tau): ", bus_tau
    write(*,"(A22,A3)") "  Berendsen barostat: ", " ON"
    write(*,"(A35,F10.3)") "    Berendsen time constant (tau): ", ber_tau
    write(*,"(A28,F20.15)") "    Berendsen constant (k): ", ber_k
end if
if (md_fix_com_mom) then
    write(*,"(A23,2X,A4)") "  MD fix COM momentum: ", "True"
else
    write(*,"(A23,2X,A5)") "  MD fix COM momentum: ", "False"
end if
write(*,"(A19,A8,2X,4(I4,2X))") "  Metadynamics CV: ", meta_cv_type, meta_cv(:)
write(*,"(A32,I10)") "  Metadynamics number of steps: ", meta_nsteps
write(*,"(A31,F10.3,A3)") "  Metadynamics timestep (tau): ", meta_tau, " fs"
write(*,"(A48,F10.3,A7)") "  Metadynamics initial deposition rate (omega): ", meta_omega, " kJ/mol"
write(*,"(A36,F10.3,A2)") "  Metadynamics bias parameter (dT): ", meta_dT, " K"
write(*,"(A39,F10.3)") "  Metadynamics Gaussian width (sigma): ", meta_sigma
write(*,*) " "
write(*,*) "Starting the metadynamics run"

! Start the simulation subroutine
CALL metadynamics_propagation(positions,xyzfile)

end subroutine

subroutine metadynamics_propagation(positions,xyzfile)
! Well-Tempered metadynamics with velocity Verlet integration scheme
use definitions, only: wp, avogad
use print_mod, only: recprt2, recprt3
use lin_alg, only: displacement_vec
use simulation_subroutines, only: init_v, get_pressure, get_temperature, get_tot_momentum
use ensemble_mod, only: berendsen_barostat, bussi_thermostat
use force_field_mod, only: get_energy_gradient, mweights, n_atoms
use parser_mod, only: atomnames, meta_cv, meta_tau, meta_nsteps, meta_cv_type, &
                md_nsteps, md_ts, md_ensemble, md_fix_com_mom, md_debug, md_temp, md_press
use propagators, only: velocity_verlet_position, velocity_verlet_velocity

implicit none
real(kind=wp),intent(inout) :: positions(n_atoms,3)
character(len=256), intent(in) :: xyzfile
integer :: meta_counter = 1
real(kind=wp) :: bias_param(meta_nsteps,2), tot_bias_pot, cv_value,step_bias_pot, init_cv_value
real(kind=wp) :: displacement(n_atoms), positions_previous(n_atoms,3), input_positions(n_atoms,3), forces(n_atoms,3), &
                total_displacement(n_atoms), velocities(n_atoms,3)
real(kind=wp) :: positions_list(2,n_atoms,3), acceleration_list(2,n_atoms,3)
real(kind=wp) :: gradnorm, tot_pot, instant_temp, pressure, E_kin, tot_momentum(3), tot_momentum_norm
character(len=256) :: traj_xyzfile, properties_outfile
integer :: istep, icartesian,i, dot, current, previous, new
logical :: suppress_flag = .true., debug_print_all_matrices = .false.

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

! Evaluate the intial value of the CV, depending on its type
if (meta_cv_type == "distance") then
    CALL distance_bias(bias_param,positions_list(current,:,:),tot_pot,forces,init_cv_value,step_bias_pot)
elseif (meta_cv_type == "angle") then
    CALL angle_bias(bias_param,positions_list(current,:,:),tot_pot,forces,init_cv_value,step_bias_pot)
elseif (meta_cv_type == "dihedral") then
    CALL dihedral_bias(bias_param,positions_list(current,:,:),tot_pot,forces,init_cv_value,step_bias_pot)
else
    write(*,*) "Error in CV type"
    stop
end if

! Initial forces and velocities
call get_energy_gradient(positions_list(current,:,:),tot_pot,forces, gradnorm, suppress_flag)
call init_v(velocities(:,:)) 

! initial acceleration, no bias on the first step
do icartesian = 1,3
    acceleration_list(current,:,icartesian) = 1e-4 * forces(:,icartesian) / mweights(:)
    !acceleration in Å/fs^2                    in kJ/mol/Å           in g/mol;
end do

! PREPARE FILE THAT TRACKS PROPERTIES
dot = index(xyzfile, ".", back=.true.)     ! find last "." in xyzfile name
properties_outfile = xyzfile(:dot-1) // ".properties.txt"
open(34, file=properties_outfile, status='replace', action='write')
write(34,"(A10,10(A20))") "istep", "E_tot", "E_kin","E_pot", "F_norm", "Temp", "Pressure", "COM_momentum", "CV_value", &
                        "Inst_Bias" , "Tot_Bias"
write(34,"(A10,10(A20))") "none","kJ/mol", "kJ/mol","kJ/mol", "kJ/(Åmol)", "K", "Pa", "gÅ/fs", " ","kJ/mol","kJ/mol"

! PREPARE TRAJECTORY FILE: traj_xyzfile
traj_xyzfile = xyzfile(:dot-1) // "_meta.traj" // xyzfile(dot:)
open(35, file=traj_xyzfile, status='replace', action='write')
write(35,*) n_atoms
write(35,"(A,F6.2,A)") "atomic positions at t = ",istep * md_ts, " fs"
do i=1, size(positions,1), 1
    write(35,FMT='(A3,3(2X,F15.8))') atomnames(i), positions(i,1:3)
end do

write(*,"(A10,7(A20))") "istep", "E_tot", "E_kin","E_pot","Temp", "Pressure", "CV_value", "Tot_Bias"
write(*,"(A10,7(A20))") "none","kJ/mol", "kJ/mol","kJ/mol", "K", "Pa", "none","kJ/mol"
write(*,'(A)') repeat('-', 150)


 do while (istep<md_nsteps)
        istep = istep +1

    ! update position with velocity Verlet
    CALL velocity_verlet_position(positions_list(current,:,:), velocities(:,:),acceleration_list(current,:,:), &
                            positions_list(new,:,:))

    ! calculate the forces with new positions
    CALL get_energy_gradient(positions_list(current,:,:),tot_pot,forces, gradnorm, suppress_flag)

    ! Deposit a new gaussian bias potential every meta_tau fs
    if (MOD((md_ts * istep),meta_tau) == 0) then
        cv_value = calc_cv(positions_list(current,:,:))
        CALL deposit_bias_potential(meta_counter,cv_value,init_cv_value,bias_param,tot_bias_pot)
    end if

    ! Calculate the bias potential and the forces contribution with the instantenous CV value
    if (meta_cv_type == "distance") then
        CALL distance_bias(bias_param,positions_list(current,:,:),tot_pot,forces,cv_value,step_bias_pot)
    elseif (meta_cv_type == "angle") then
        CALL angle_bias(bias_param,positions_list(current,:,:),tot_pot,forces,cv_value,step_bias_pot)
    elseif (meta_cv_type == "dihedral") then
        CALL dihedral_bias(bias_param,positions_list(current,:,:),tot_pot,forces,cv_value,step_bias_pot)
    else
        write(*,*) "Error in CV type"
        stop
    end if
    

    do icartesian = 1,3
        acceleration_list(new,:,icartesian) = 1e-4 * forces(:,icartesian) / mweights(:)
        !acceleration in Å/fs^2                    in kJ/mol/Å           in g/mol;
    end do
    
    ! Update velocities with Velocity Verlet scheme
    CALL velocity_verlet_velocity(acceleration_list(current,:,:), acceleration_list(new,:,:), velocities(:,:))

    ! SCALE VELOCITIES TO ENSURE ZERO MOMENTUM OF THE CENTER OF MASS
    call get_tot_momentum(velocities, tot_momentum)
    tot_momentum_norm = 0
    do icartesian = 1,3
        tot_momentum_norm = tot_momentum_norm + tot_momentum(icartesian)**2
    end do
    tot_momentum_norm = SQRT(tot_momentum_norm)

    if (md_fix_com_mom) then
        do icartesian = 1,3
            velocities(:,icartesian) = velocities(:,icartesian) - tot_momentum(icartesian) / SUM(mweights)
        end do

        if (md_debug) then
            write(*,"(/A)") "After adapting the velocities with respect to the momentum and mass:"
            call get_tot_momentum(velocities, tot_momentum) !THIS SHOULD NOW BE ZERO
            if (debug_print_all_matrices) then
                call recprt3("v(t_0) = velocities(:,:) [Å/fs]",velocities(:,:),n_atoms)
            end if
        end if
    end if

    ! GET PROPERTIES
    call get_temperature(velocities, instant_temp, E_kin) ! use v(t)
    call get_pressure(positions_list(current,:,:), forces,instant_temp, pressure) ! use x(t) and v(t)

    ! Apply thermostat/barostat constraints and compute again the properties
    if (md_ensemble == "NVT") then
        CALL bussi_thermostat(E_kin,(3*n_atoms)-3,velocities)
        call get_temperature(velocities, instant_temp, E_kin) 
        call get_pressure(positions_list(current,:,:), forces,instant_temp, pressure)
    elseif (md_ensemble == "NPT") then
        CALL bussi_thermostat(E_kin,(3*n_atoms)-3,velocities)
        CALL berendsen_barostat(positions_list(current,:,:),pressure)
        call get_temperature(velocities, instant_temp, E_kin) 
        call get_pressure(positions_list(current,:,:), forces,instant_temp, pressure)
    end if

    ! WRITE QUANTITIES FILE
    ! this prints info from the previous step
    write(34,"(I8,2x,10(F20.8))") istep-1,E_kin+tot_pot,E_kin,tot_pot, gradnorm, instant_temp, pressure,& 
                                tot_momentum_norm/avogad, cv_value, step_bias_pot, tot_bias_pot

    if (debug_print_all_matrices) then
        write(*,"(/A,I5)") "New quantities at step ",istep
        call recprt2("r(t) = positions_list(current,:,:) [Å]",atomnames,positions_list(current,:,:),n_atoms)
        call recprt2("F(t) = forces(:,:) [kJ/mol]",atomnames,forces(:,:),n_atoms)
        call recprt2("a(t) = acceleration(:,:) [Å/(fs)^2]",atomnames,acceleration_list(new,:,:),n_atoms)
        call recprt2("v(t) = velocities(:,:) [Å/fs]",atomnames,velocities(:,:),n_atoms)
    end if

    if (mod(istep,100) == 1) then
        write(*,"(I8,2x,7(F20.8))") istep-1,E_kin+tot_pot,E_kin,tot_pot, instant_temp, pressure,& 
                                 cv_value, tot_bias_pot
    end if

    ! TRACK DISPLACEMENT OF THE ATOMS
    call displacement_vec(positions_list(new,:,:),positions_list(current,:,:),displacement,n_atoms,atomnames)
    total_displacement(:) = total_displacement(:) + displacement(:)

    ! WRITE TRAJECTORY FILE
    write(35,*) n_atoms
    write(35,"(A,F6.2,A)") "atomic positions at t = ",istep * md_ts, " fs"
    do i=1, size(positions,1), 1
        write(35,FMT='(A3,3(2X,F15.8))') atomnames(i), positions_list(new,i,:)
    end do

    ! PREPARE NEXT STEP
    positions_list(current,:,:) = positions_list(new,:,:)
    acceleration_list(current,:,:) = acceleration_list(new,:,:)

end do
write(*,*) "last step"
write(*,"(I8,2x,7(F20.8))") istep-1,E_kin+tot_pot,E_kin,tot_pot, instant_temp, pressure,& 
                                 cv_value, tot_bias_pot
write(*,*) " "
write(*,*) "Metadynamics succesfully completed"
end subroutine

subroutine deposit_bias_potential(meta_counter,cv_value,init_cv_value,bias_param,tot_bias_pot)
! Subroutine to deposit additional Gaussian function to bias potential
use definitions, only: wp
use force_field_mod, only: n_atoms
use parser_mod, only: meta_cv, meta_dT, meta_omega, meta_sigma, meta_nsteps, meta_tau 

implicit none 

real(kind=wp), intent(inout) :: tot_bias_pot, bias_param(meta_nsteps,2),init_cv_value 
integer, intent(inout) :: meta_counter
real(kind=wp), intent(in) :: cv_value
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
! Calculate the instaneous value of the CV
use definitions, only: wp
use force_field_mod, only: n_atoms, cross_product
use parser_mod, only: meta_cv, meta_cv_type

implicit none
real(kind=wp), intent(in) :: positions(n_atoms,3)
real(kind=wp) :: cv_value
real(kind=wp) :: distance(3), d12(3), d23(3), d12_norm, d23_norm, d34(3), a(3), b(3), a_norm, b_norm, ratio

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

subroutine distance_bias(bias_param,positions,tot_pot,forces,cv_value,bias_pot)
! Compute the bias potential and the force contribution for the distance CV type
use definitions, only: wp
use force_field_mod, only: n_atoms
use parser_mod, only: meta_cv, meta_nsteps, meta_cv_type, meta_sigma

implicit none

real(kind=wp), intent(in) :: bias_param(meta_nsteps,2), positions(n_atoms,3)
real(kind=wp), intent(inout) :: tot_pot,forces(n_atoms,3),cv_value 
real(kind=wp), intent(out) :: bias_pot
real(kind=wp) :: distance(3), der_1(3), der_2(3)
integer :: i
bias_pot = 0

! Distance CV
distance = positions(meta_cv(2),:) - positions(meta_cv(1),:)
cv_value = SQRT(dot_product(distance,distance))

! iterate over the sum terms of the bias potential function
do i=1, meta_nsteps
    ! Bias potential calculation
    bias_pot =bias_pot + bias_param(i,1) * ( exp(-0.5*((cv_value-bias_param(i,2)) / meta_sigma)**2))

    ! Derivative over atomic coordinates of atom 1
    der_1= - bias_param(i,1) * ((cv_value - bias_param(i,2))/meta_sigma) * &
            exp(-0.5*((cv_value-bias_param(i,2)) / meta_sigma)**2) * &
            (positions(meta_cv(2),:) - positions(meta_cv(1),:)) / cv_value

    ! Derivative over atomic coordinates of atom 2
    der_2 = bias_param(i,1) * ((cv_value - bias_param(i,2))/meta_sigma) * &
            exp(-0.5*((cv_value-bias_param(i,2)) / meta_sigma)**2) * &
            (positions(meta_cv(2),:) - positions(meta_cv(1),:)) / cv_value

    ! Update the forces with the bias contribution
    forces(meta_cv(1),:) = forces(meta_cv(1),:) + der_1
    forces(meta_cv(2),:) = forces(meta_cv(2),:) + der_2
end do

! Update the total potential with the bias contribution
tot_pot = tot_pot + bias_pot
end subroutine

subroutine angle_bias(bias_param,positions,tot_pot,forces,cv_value,bias_pot)
! Compute the bias potential and the force contribution for the angle CV type
use definitions, only: wp, safeguard
use force_field_mod, only: n_atoms
use parser_mod, only: meta_cv, meta_nsteps, meta_cv_type, meta_sigma

implicit none

real(kind=wp), intent(in) :: bias_param(meta_nsteps,2)
real(kind=wp), intent(inout) :: positions(n_atoms,3),tot_pot,forces(n_atoms,3), cv_value 
real(kind=wp), intent(out) :: bias_pot
integer :: i
real(kind=wp) :: d12(3), d23(3), d12_norm, d23_norm, der_magn, der_1(3), der_3(3)

bias_pot = 0
d12 = positions(meta_cv(1),1:3) - positions(meta_cv(2),1:3)
d23 = positions(meta_cv(3),1:3) - positions(meta_cv(2),1:3)
d12_norm = SQRT(dot_product(d12,d12))
d23_norm = SQRT(dot_product(d23,d23))

! Calculate the angle (in radians)
cv_value = acos(dot_product(d12,d23) / (d12_norm * d23_norm))

! iterate over the sum terms of the bias potential function
do i=1, meta_nsteps
    ! Bias potential calculation
    bias_pot = bias_pot + bias_param(i,1) * ( exp(-0.5*((cv_value-bias_param(i,2)) / meta_sigma)**2))

    ! Magnitude of the gradient
    der_magn = bias_param(i,1) * ((cv_value - bias_param(i,2))/meta_sigma) * &
            exp(-0.5*((cv_value-bias_param(i,2)) / meta_sigma)**2)

    ! Derivative over atomic coordinates of atom 1
    der_1 = - der_magn * (1 / ( sin(cv_value) * d12_norm + safeguard)) * &
                    ( cos(cv_value) * (d12 / d12_norm) - (d23 / d23_norm ))

    ! Derivative over atomic coordinates of atom 3
    der_3 = - der_magn * (1 / ( sin(cv_value) * d23_norm + safeguard)) * &
                    (cos(cv_value) * (d23 / d23_norm) - (d12 / d12_norm))

    ! Update the forces with the bias contribution                    
    forces(meta_cv(1),1:3) = forces(meta_cv(1),1:3) + der_1
    forces(meta_cv(2),1:3) = forces(meta_cv(2),1:3) - (der_1 + der_3)
    forces(meta_cv(3),1:3) = forces(meta_cv(3),1:3) + der_3
end do

! Update the total potential with the bias contribution
tot_pot = tot_pot + bias_pot
end subroutine

subroutine dihedral_bias(bias_param,positions,tot_pot,forces,cv_value,bias_pot)
! Compute the bias potential and the force contribution for the dihedral CV type
use definitions, only: wp, safeguard
use force_field_mod, only: n_atoms, cross_product
use parser_mod, only: meta_cv, meta_nsteps, meta_cv_type, meta_sigma

implicit none

real(kind=wp), intent(in) :: bias_param(meta_nsteps,2)
real(kind=wp), intent(inout) :: positions(n_atoms,3),tot_pot,forces(n_atoms,3), cv_value
real(kind=wp), intent(out) :: bias_pot
integer :: i
real(kind=wp) :: d12(3), d23(3), d34(3), a(3), b(3), a_norm, b_norm, der_magn, &
                ratio, cap_A(3), cap_B(3), der_1(3), der_2(3), der_3(3), der_4(3)
bias_pot = 0

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
! Calculate the dihedral angle for the CV
cv_value = acos(ratio)

! Define A and B terms
cap_A = (b / (a_norm * b_norm)) - ((cos(cv_value) * a) / (a_norm**2))
cap_B = (a / (a_norm * b_norm)) - ((cos(cv_value) * b) / (b_norm**2 ))

! iterate over the sum terms of the bias potential function
do i=1, meta_nsteps
    ! Bias potential calculation
    bias_pot = bias_pot + bias_param(i,1) * ( exp(-0.5*((cv_value-bias_param(i,2)) / meta_sigma)**2))

    ! Magnitude of the gradient
    der_magn = bias_param(i,1) * ((cv_value - bias_param(i,2))/meta_sigma) * &
            exp(-0.5*((cv_value-bias_param(i,2)) / meta_sigma)**2) 

    ! Derivative over atomic coordinates
    der_1 = der_magn * (1/ (sin(cv_value) + safeguard)) * cross_product(cap_A,d23)
    der_4 = - der_magn * (1/ (sin(cv_value) + safeguard)) * cross_product(d23,cap_B)
    der_2 = - der_1 + (der_magn * (1/ (sin(cv_value) + safeguard)) * (cross_product(d12,cap_A) + cross_product(cap_B,d34)))
    der_3 = - der_4 - (der_magn * (1/ (sin(cv_value) + safeguard)) * (cross_product(cap_A,d12) + cross_product(cap_B,d34)))

    ! Update the forces with the bias contribution 
    forces(meta_cv(1),1:3) = forces(meta_cv(1),1:3) + der_1
    forces(meta_cv(2),1:3) = forces(meta_cv(2),1:3) + der_2 
    forces(meta_cv(3),1:3) = forces(meta_cv(3),1:3) + der_3
    forces(meta_cv(4),1:3) = forces(meta_cv(4),1:3) + der_4
end do

! Update the total potential with the bias contribution
tot_pot = tot_pot + bias_pot
end subroutine

end module
