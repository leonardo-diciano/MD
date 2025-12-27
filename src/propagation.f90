module propagation
contains

subroutine propagator(positions,positions_previous,mweights,n_atoms, debug_flag, atomnames,xyzfile)!,md_ts,nsteps)
    use definitions, only: wp
    use print_mod, only: recprt2
    use lin_alg, only: displacement_vec
    use force_field_mod, only: get_energy_gradient
    use parser_mod, only: md_ts,md_nsteps,md_ensemble,md_barostat,md_thermostat,md_temp

    implicit none

    real(kind=wp),intent(in) ::  mweights(n_atoms)
    real(kind=wp),intent(inout) :: positions(n_atoms,3)
    integer, intent(in) :: n_atoms
    logical, intent(in) :: debug_flag
    character(len=2), intent(in) :: atomnames(:)
    character(len=256), intent(in) :: xyzfile


    real(kind=wp) :: displacement(n_atoms), positions_previous(n_atoms,3), input_positions(n_atoms,3), forces(n_atoms,3), &
                    acceleration(n_atoms,3), total_displacement(n_atoms), velocities(n_atoms,3)
    real(kind=wp) :: positions_list(3,n_atoms,3)
    real(kind=wp) :: gradnorm, tot_pot, instant_temp, pressure
    character(len=256) :: traj_xyzfile, properties_outfile
    integer :: istep, icartesian,i, dot, current, previous, new
    logical :: suppress_flag = .true., debug = .false.

    ! for intuitive storing of position data (since Verlet requires storing previous positions)
    previous = 1
    current = 2
    new = 3

    ! INITIALIZATION
    istep = 0
    total_displacement(:) = 0
    input_positions(:,:) = positions(:,:)

    call init_v(input_positions,velocities, n_atoms, mweights, debug_flag, md_temp)
    call get_pressure(positions, forces,instant_temp,n_atoms, pressure)

    positions_list(previous,:,:) = positions(:,:) - velocities(:,:) * md_ts
    positions_list(current,:,:) = positions(:,:)

    if (debug) then
        write(*,"(/A,/A,/A)") "In propagation","---------------------------------------","Initialization:"
        call recprt2("forces",atomnames,forces,n_atoms)
        write(*,"(A,*(/,F10.6))") "masses: [g/mol]", mweights
       call recprt2("acceleration = forces / masses",atomnames,acceleration,n_atoms)
        write(*,"(A)") "Start performing steps ..."
    end if

    ! PREPARE TRAJECTORY FILE: traj_xyzfile
    dot = index(xyzfile, ".", back=.true.)     ! find last "." in xyzfile name
    traj_xyzfile = xyzfile(:dot-1) // ".traj" // xyzfile(dot:)
    open(98, file=traj_xyzfile, status='replace', action='write')
    write(98,*) n_atoms
    write(98,"(A,F6.2,A)") "atomic positions at t = ",istep * md_ts, " fs"
    do i=1, size(positions,1), 1
        write(98,FMT='(A3,3(2X,F15.8))') atomnames(i), positions(i,1:3)
    end do

    ! PREPARE FILE THAT TRACKS PROPERTIES
    properties_outfile = xyzfile(:dot-1) // ".properties.txt"
    open(97, file=properties_outfile, status='replace', action='write')
    write(97,*) "properties"
    write(*,*) "E_tot", "E_kin","E_pot", "F_norm", "Temp"


    if (debug) then
        write(*,"(/A,I5)") "Initial quantities at step ",istep
        call recprt2("r(t-Δt)",atomnames,positions_list(previous,:,:),n_atoms)
        call recprt2("r(t)",atomnames,positions_list(current,:,:),n_atoms)
        !call recprt2("r(t+Δt)",atomnames,positions_list(new,:,:),n_atoms)
        !call recprt2("acceleration",atomnames,acceleration,n_atoms)
    end if

    ! HERE THE STEPS ARE TAKEN
    do while (istep<md_nsteps)
        istep = istep +1

        !CALCULATE FORCES / ACCELERATION AT CURRENT POSITION
        call get_energy_gradient(positions_list(current,:,:),tot_pot,forces, gradnorm, suppress_flag)
        do icartesian = 1,3
            acceleration(:,icartesian) = 1e-4 * forces(:,icartesian) / mweights(:)
            !acceleration in Å/fs^2                    in kJ/mol/Å           in g/mol;
        end do

        !UPDATE POSITIONS
        !positions_list(new,:,:) = 0

        ! WRITE QUANTITIES FILE
        open(97, file=properties_outfile, status='old', action='write')
        call update_pos_Verlet(positions_list(previous,:,:),positions_list(current,:,:),md_ts,acceleration(:,:), &
                            n_atoms, positions_list(new,:,:))

        call displacement_vec(positions_list(new,:,:),positions_list(current,:,:),n_atoms,atomnames,displacement)
        total_displacement(:) = total_displacement(:) + displacement(:)


        if (debug) then
            call recprt2("New forces",atomnames,forces,n_atoms)
            call recprt2("New accelerations",atomnames,acceleration,n_atoms)
            call displacement_vec(positions_list(new,:,:),positions_list(current,:,:),n_atoms,atomnames,displacement)
            total_displacement(:) = total_displacement(:) + displacement
            write(*,*) "Displacements"
            do i = 1, n_atoms
                write(*,"(I3,1x,A3,1x,F16.12,1x,A)") i,atomnames(i),displacement(i),"Å"
            end do
            write(*,"(/A,F12.8,A)") "  sum = ", sum(displacement(:)), " Å"
        end if
        positions_list(new,:,:) = 0
        positions_list(new,:,:) = 2 * positions_list(current,:,:) - positions_list(previous,:,:) &
                                    + md_ts**2 * acceleration(:,:)


        ! WRITE TRAJECTORY FILE
        open(98, file=traj_xyzfile, status='old', action='write')
        write(98,*) n_atoms
        write(98,"(A,F6.2,A)") "atomic positions at t = ",istep * md_ts, " fs"
        do i=1, size(positions,1), 1
            write(98,FMT='(A3,3(2X,F15.8))') atomnames(i), positions_list(new,i,:)
        end do


        if (debug_flag) then
            write(*,"(/A,I5)") "New quantities at step ",istep
            call recprt2("r(t-Δt) = positions_list(previous,:,:) [Å]",atomnames,positions_list(previous,:,:),n_atoms)
            call recprt2("r(t) = positions_list(current,:,:) [Å]",atomnames,positions_list(current,:,:),n_atoms)
            call recprt2("r(t+Δt) = positions_list(new,:,:) [Å]",atomnames,positions_list(new,:,:),n_atoms)
            call recprt2("F(t) = forces(:,:) [kJ/mol]",atomnames,forces(:,:),n_atoms)
            call recprt2("a(t) = acceleration(:,:) [Å/(fs)^2]",atomnames,acceleration(:,:),n_atoms)
        end if

        positions_list(previous,:,:) = positions_list(current,:,:)
        positions_list(current,:,:) = positions_list(new,:,:)

    end do

    displacement(:) = 0
    write(*,"(/A)") "Throughout the simulation, the atoms displaced: (no MSD, but initial vs final coords)"
    call displacement_vec(positions_list(current,:,:),input_positions,n_atoms,atomnames, displacement)
    write(*,*) "Displacements (summed all steps)"
    do i = 1, n_atoms
        write(*,"(I3,1x,A3,1x,  F16.12,1x,A)") i,atomnames(i),displacement(i),"Å"
    end do

    write(*,"(/A)") "Displacements (summed all steps)"
    do i = 1, n_atoms
        write(*,"(I3,1x,A3,1x,F16.12,1x,A)") i,atomnames(i),total_displacement(i),"Å"
    end do

    write(*,"(/A,A)") "Trajectory was written to: ", traj_xyzfile
    write(*,"(A,F10.2,A)") "Total simulation time", md_ts * istep, " fs"

end subroutine propagator

subroutine update_pos_Verlet(positions_previous,positions_current,md_ts,acceleration, n_atoms, positions_new)
    use definitions, only: wp
    implicit none
    real(kind=wp), intent(in) :: positions_previous(n_atoms,3), positions_current(n_atoms,3),md_ts,acceleration(n_atoms,3)
    integer, intent(in) :: n_atoms
    real(kind=wp), intent(out) :: positions_new(n_atoms,3)

    positions_new(:,:) = 2 * positions_current(:,:) - positions_previous(:,:) + md_ts**2 * acceleration(:,:)

    !positions_list(new,:,:) = 2 * positions_list(current,:,:) - positions_list(previous,:,:) &
    !                        + md_ts**2 * acceleration(:,:)

end subroutine update_pos_Verlet

subroutine init_v(positions,velocities, n_atoms, mweights, debug_flag, md_temp)
    use definitions, only: wp, pi, kB, proton_mass
    use print_mod, only: recprt3
    implicit none
    real(kind=wp), intent(in) :: mweights(n_atoms), positions(n_atoms,3), md_temp
    integer, intent(in) :: n_atoms
    logical, intent(in) :: debug_flag
    real(kind=wp), intent(out) :: velocities(n_atoms,3)

    real(kind=wp) :: rand1,rand2, rand_gaussian, v_ix, v_atoms(n_atoms), tot_momentum(3)
    real(kind=wp) :: E_kin, instant_temp
    integer :: iatom, icartesian

    write(*,"(/A,/A)") "in init_v", "-------------------------------"

    do iatom = 1,n_atoms
        do icartesian = 1, 3
            call random_number(rand1)
            call random_number(rand2)
            rand_gaussian = SQRT(-2*LOG(rand1)) * COS(2*pi*rand2)
            v_ix = SQRT(kB*md_temp/(proton_mass*mweights(iatom))) * rand_gaussian

            !write(*,"(/,2(A,1x,F12.8,1x),A,F16.8)") "rand1 =",rand1, "rand2 =", rand2,"rand_gaussian = ",rand_gaussian
            !write(*,"(A,F16.4,A)") "v_ix = ", v_ix, " m/s"

            ! convert to Å/fs -> 1m/s = 10^10/10^15 Å/fs = 10^-5 Å/fs

            v_ix = 1e-5 * v_ix
            !write(*,"(A,F14.8,A)") "v_ix = ", v_ix, " Å/fs"

            velocities(iatom,icartesian) = v_ix
        end do
    end do

    !call recprt3("v(t_0) = velocities(:,:) [Å/fs]",velocities(:,:),n_atoms)
    !call get_v_atoms(v_atoms,velocities,n_atoms,.true.)

    ! THIS IS LIKELY NOT YET ZERO
    call get_tot_momentum(velocities,mweights, n_atoms, tot_momentum)

    ! NOW WE RECALCULATE THE VELOCITIES TO MAKE THE TOTAL MOMENTUM ZERO (Leach, p.365)
    do icartesian = 1,3
        velocities(:,icartesian) = velocities(:,icartesian) - tot_momentum(icartesian) / SUM(mweights)
    end do
    write(*,"(/A)") "After adapting the velocities with respect to the momentum and mass:"
    call get_tot_momentum(velocities,mweights, n_atoms, tot_momentum) !THIS SHOULD NOW BE ZERO

    call recprt3("v(t_0) = velocities(:,:) [Å/fs]",velocities(:,:),n_atoms)
    call get_v_atoms(v_atoms,velocities,n_atoms,.true.)

    call get_temperature(velocities, mweights,n_atoms, instant_temp)

    end subroutine init_v

subroutine get_v_atoms(v_atoms,velocities,n_atoms,printopt)
    use definitions, only: wp
    implicit none
    real(kind=wp),intent(out) :: v_atoms(n_atoms)
    logical, intent(in) :: printopt
    integer, intent(in) :: n_atoms
    real(kind=wp), intent(in) :: velocities(n_atoms,3)

    integer :: iatom, icartesian
    real(kind=wp) :: v_atom

    v_atoms(:) = 0
    if (printopt) then
        write(*,*) "v_atoms = ["
    end if
    !calculate the norm of the velocity vector on each atom
    do iatom = 1,n_atoms
        v_atom = 0
        do icartesian = 1,3
            v_atom = v_atom + velocities(iatom,icartesian)**2
        end do
        v_atom = SQRT(v_atom)
        v_atoms(iatom) = v_atom
        !if (printopt) then
            write(*,"(F16.8,A)") v_atom, ","
            !write(*,"(A,I3,A,F16.8,A)") "v_atom of atom ", iatom, "= ", v_atom, " Å/fs"
        !end if
    end do
    if (printopt) then
        write(*,*) "]"
    end if
end subroutine get_v_atoms

subroutine get_tot_momentum(velocities,mweights, n_atoms, tot_momentum)
    use definitions, only: wp, proton_mass
    implicit none
    real(kind=wp), intent(in) :: velocities(n_atoms,3), mweights(n_atoms)
    integer, intent(in) :: n_atoms
    real(kind=wp), intent(out) :: tot_momentum(3)
    integer :: iatom, icartesian

    tot_momentum(:) = 0

    do iatom=1,n_atoms
        do icartesian=1,3
            tot_momentum(icartesian) = tot_momentum(icartesian) + mweights(iatom) * velocities(iatom,icartesian)
            !g Å/(mol fs)                                            !g/mol                            !Å/fs
        end do
    end do

    write(*,"(/A,3(F10.6),A)") "The total momentum of the system is ", tot_momentum, " g Å/(mol fs)"
end subroutine get_tot_momentum

subroutine get_E_kin(velocities, mweights, n_atoms, E_kin)
    use definitions, only: wp
    implicit none
    real(kind=wp), intent(in) :: velocities(n_atoms,3), mweights(n_atoms)
    integer, intent(in) :: n_atoms
    real(kind=wp), intent(out) :: E_kin

    integer :: icartesian, iatom

    E_kin = 0

    do iatom = 1,n_atoms
        do icartesian = 1, 3
            !E_kin = E_kin + 0.5 * mweights(iatom)*velocities(iatom,icartesian)**2
                                    !g/mol               !Å^2/fs^2
            !E_kin = E_kin + 0.5 * 1e-3 * mweights(iatom)* 1e10 * velocities(iatom,icartesian)**2
            !J/mol                       !kg/mol               ! m^2/s^2
            E_kin = E_kin + 0.5 * 1e4 * mweights(iatom) * velocities(iatom,icartesian)**2
            !kJ/mol                       !kg/mol               ! m^2/s^2
        end do
    end do

    write(*,*) "The kinetic energy is", E_kin, "kJ/mol"
end subroutine get_E_kin

subroutine get_temperature(velocities,mweights,n_atoms, instant_temp)
    use definitions, only: wp, boltzmann
    implicit none
    real(kind=wp), intent(in) ::velocities(n_atoms,3), mweights(n_atoms)
    integer, intent(in) :: n_atoms
    real(kind=wp), intent(out) :: instant_temp !/K

    real(kind=wp) :: E_kin !kJ/mol

    call get_E_kin(velocities, mweights, n_atoms, E_kin)
    instant_temp = 2 * E_kin /(boltzmann*3*n_atoms)
                       !kJ/mol  !kJ/(K mol)

    write(*,*) "The temperature is ", instant_temp, "K"
end subroutine get_temperature

subroutine get_pressure(positions, forces,instant_temp,n_atoms, pressure)
    use definitions, only: wp, boltzmann, avogad
    use parser_mod, only: md_boxlength
    implicit none
    real(kind=wp), intent(in) ::instant_temp, positions(n_atoms,3),forces(n_atoms,3)
    integer, intent(in) :: n_atoms
    real(kind=wp), intent(out) :: pressure

    real(kind=wp) :: volume, avg
    integer :: icartesian, iatom

    volume = md_boxlength * md_boxlength * md_boxlength ! Å^3

    avg = 0
    do iatom = 1, n_atoms
        do icartesian = 1, 3
            avg = avg + positions(iatom, icartesian) * forces(iatom, icartesian)
            !kJ/mol          !Å                             !kJ/mol/Å
        end do
    end do
    pressure = 1/(3*volume)*(avg+3*n_atoms*boltzmann*instant_temp)
    !kJ/mol/Å^3      !Å^3              !kJ/mol

    pressure = pressure * 1e30 / avogad

    write(*,*) "The pressure is ",pressure,"kJ/m^3 = kPa"

end subroutine get_pressure
end module propagation
