module simulation_mod
contains

subroutine simulation(positions,xyzfile)
    use definitions, only: wp
    use force_field_mod, only: n_atoms
    use parser_mod, only: md_int
    implicit none
    real(kind=wp),intent(inout) :: positions(n_atoms,3)
    character(len=256), intent(in) :: xyzfile

    if (md_int == "velocity_verlet" ) then
        call simulation_vel_verlet(positions,xyzfile)
    else
        call simulation_verlet(positions,xyzfile)
    end if

end subroutine simulation

subroutine simulation_verlet(positions,xyzfile)
    use definitions, only: wp, avogad
    use print_mod, only: recprt2, recprt3
    use lin_alg, only: displacement_vec
    use force_field_mod, only: get_energy_gradient, n_atoms, mweights
    use parser_mod, only: atomnames,md_ts,md_nsteps,md_ensemble,md_temp, md_press, md_debug, md_fix_com_mom
    use simulation_subroutines, only: init_v, get_pressure, get_temperature, get_tot_momentum
    use ensemble_mod, only: bussi_thermostat, berendsen_barostat
    use propagators, only: Verlet

    implicit none

    real(kind=wp),intent(inout) :: positions(n_atoms,3)
    character(len=256), intent(in) :: xyzfile


    real(kind=wp) :: displacement(n_atoms), positions_previous(n_atoms,3), input_positions(n_atoms,3), forces(n_atoms,3), &
                    acceleration(n_atoms,3), total_displacement(n_atoms), velocities(n_atoms,3)
    real(kind=wp) :: positions_list(3,n_atoms,3)
    real(kind=wp) :: gradnorm, tot_pot, instant_temp, pressure, E_kin, tot_momentum(3), tot_momentum_norm
    character(len=256) :: traj_xyzfile, properties_outfile
    integer :: istep, icartesian,i, dot, current, previous, new
    logical :: suppress_flag = .true., debug_print_all_matrices = .false.

    ! for intuitive storing of position data (since Verlet requires storing previous positions)
    previous = 1
    current = 2
    new = 3

    ! INITIALIZATION
    istep = 0
    total_displacement(:) = 0
    input_positions(:,:) = positions(:,:) !x(t=0)
    positions_list(current,:,:) = positions(:,:) !x(t=0)

    call get_energy_gradient(positions_list(current,:,:),tot_pot,forces, gradnorm, suppress_flag) !F(t=0), E_pot(t=0)
    call init_v(velocities(:,:)) !v(t=-1)
    do icartesian = 1,3
        acceleration(:,icartesian) = 1e-4 * forces(:,icartesian) / mweights(:)
        !acceleration in Å/fs^2                    in kJ/mol/Å           in g/mol;
    end do

    ! PREPARE FILE THAT TRACKS PROPERTIES
    dot = index(xyzfile, ".", back=.true.)     ! find last "." in xyzfile name
    properties_outfile = xyzfile(:dot-1) // ".properties.txt"
    open(97, file=properties_outfile, status='replace', action='write')
    write(97,"(A10,7(A20))") "istep", "E_tot", "E_kin","E_pot", "F_norm", "Temp", "Pressure", "COM_momentum"
    write(97,"(A10,7(A20))") "none","kJ/mol", "kJ/mol","kJ/mol", "kJ/(Åmol)", "K", "Pa", "gÅ/fs"

    positions_list(previous,:,:) = positions(:,:) - velocities(:,:) * md_ts

    if (md_debug) then
        write(*,"(/A,/A,/A)") "In propagation","---------------------------------------","Initialization:"
        call recprt2("forces",atomnames,forces,n_atoms)
        write(*,"(A,*(/,F10.6))") "masses: [g/mol]", mweights
       call recprt2("acceleration = forces / masses",atomnames,acceleration,n_atoms)
        write(*,"(A)") "Start performing steps ..."
    end if

    ! PREPARE TRAJECTORY FILE: traj_xyzfile
    traj_xyzfile = xyzfile(:dot-1) // ".traj" // xyzfile(dot:)
    open(98, file=traj_xyzfile, status='replace', action='write')
    write(98,*) n_atoms
    write(98,"(A,F6.2,A)") "atomic positions at t = ",istep * md_ts, " fs"
    do i=1, size(positions,1), 1
        write(98,FMT='(A3,3(2X,F15.8))') atomnames(i), positions(i,1:3)
    end do


    if (md_debug) then
        write(*,"(/A,I5)") "Initial quantities at step ",istep
        call recprt2("r(t-Δt)",atomnames,positions_list(previous,:,:),n_atoms)
        call recprt2("r(t)",atomnames,positions_list(current,:,:),n_atoms)
        !call recprt2("r(t+Δt)",atomnames,positions_list(new,:,:),n_atoms)
        !call recprt2("acceleration",atomnames,acceleration,n_atoms)
    end if

    if (md_ensemble == "NVT") then
        write(*,*) "Turning on the Bussi thermostat"
    elseif (md_ensemble == "NPT") then
        write(*,*) "Turning on the Bussi thermostat"
        write(*,*) "Turning on the Berendsend barostat"
    end if
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! PERFORM THE STEPS ITERATIVELY
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do while (istep<md_nsteps)
        istep = istep +1

        ! UPDATE POSITIONS
        call Verlet(positions_list(previous,:,:),positions_list(current,:,:),acceleration(:,:), &
                            positions_list(new,:,:), velocities(:,:)) ! get x(t+1) and v(t)

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

        ! Apply thermostat/barostat constraints
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
        open(97, file=properties_outfile, status='old', action='write')
        ! this prints info from the previous step
        write(97,"(I8,2x,7(F20.8))") istep-1,E_kin+tot_pot,E_kin,tot_pot, gradnorm, instant_temp, pressure, &
                                            tot_momentum_norm/avogad
        if (debug_print_all_matrices) then
            write(*,"(/A,I5)") "New quantities at step ",istep
            call recprt2("r(t) = positions_list(current,:,:) [Å]",atomnames,positions_list(current,:,:),n_atoms)
            call recprt2("F(t) = forces(:,:) [kJ/mol]",atomnames,forces(:,:),n_atoms)
            call recprt2("a(t) = acceleration(:,:) [Å/(fs)^2]",atomnames,acceleration(:,:),n_atoms)
            call recprt2("v(t) = velocities(:,:) [Å/fs]",atomnames,velocities(:,:),n_atoms)
        end if

        ! TRACK DISPLACEMENT OF THE ATOMS
        call displacement_vec(positions_list(new,:,:),positions_list(current,:,:),displacement,n_atoms,atomnames)
        total_displacement(:) = total_displacement(:) + displacement(:)

        ! WRITE TRAJECTORY FILE
        open(98, file=traj_xyzfile, status='old', action='write')
        write(98,*) n_atoms
        write(98,"(A,F6.2,A)") "atomic positions at t = ",istep * md_ts, " fs"
        do i=1, size(positions,1), 1
            write(98,FMT='(A3,3(2X,F15.8))') atomnames(i), positions_list(new,i,:)
        end do

        ! PREPARE NEXT STEP
        positions_list(previous,:,:) = positions_list(current,:,:)
        positions_list(current,:,:) = positions_list(new,:,:)

        ! CALCULATE NEW FORCES / ACCELERATION AT NEW POSITION
        call get_energy_gradient(positions_list(new,:,:),tot_pot,forces, gradnorm, suppress_flag)
        do icartesian = 1,3
            acceleration(:,icartesian) = 1e-4 * forces(:,icartesian) / mweights(:)
            !acceleration in Å/fs^2                    in kJ/mol/Å           in g/mol;
        end do

    end do

    displacement(:) = 0

    if (md_debug) then
        write(*,"(/A)") "Throughout the simulation, the atoms displaced: (no MSD, but initial vs final coords)"
        call displacement_vec(positions_list(current,:,:),input_positions, displacement,n_atoms,atomnames)
        write(*,*) "Displacements (summed all steps)"
        do i = 1, n_atoms
            write(*,"(I3,1x,A3,1x,  F16.12,1x,A)") i,atomnames(i),displacement(i),"Å"
        end do

        write(*,"(/A)") "Displacements (summed all steps)"
        do i = 1, n_atoms
            write(*,"(I3,1x,A3,1x,F16.12,1x,A)") i,atomnames(i),total_displacement(i),"Å"
        end do
    end if

    write(*,"(/A,//A)") "================================================================","Simulation finished"
    write(*,"(/A,F10.2,A)") "Total simulation time", md_ts * istep, " fs"
    write(*,"(/A,A)") "Wrote properties to ", properties_outfile
    write(*,"(A,A)") "Wrote trajectory to ", traj_xyzfile
    write(*,"(     A)") "================================================================"

end subroutine simulation_verlet


subroutine simulation_vel_verlet(positions,xyzfile)
    use definitions, only: wp, avogad
    use print_mod, only: recprt2, recprt3
    use lin_alg, only: displacement_vec
    use force_field_mod, only: get_energy_gradient, n_atoms, mweights
    use parser_mod, only: atomnames,md_ts,md_nsteps,md_ensemble,md_temp, md_press, md_debug, md_fix_com_mom
    use simulation_subroutines, only: init_v, get_pressure, get_temperature, get_tot_momentum
    use propagators, only: velocity_verlet_position, velocity_verlet_velocity

    implicit none

    real(kind=wp),intent(inout) :: positions(n_atoms,3)
    character(len=256), intent(in) :: xyzfile


    real(kind=wp) :: displacement(n_atoms), positions_previous(n_atoms,3), input_positions(n_atoms,3), forces(n_atoms,3), &
                    total_displacement(n_atoms), velocities(n_atoms,3)
    real(kind=wp) :: positions_list(2,n_atoms,3), acceleration_list(2,n_atoms,3)
    real(kind=wp) :: gradnorm, tot_pot, instant_temp, pressure, E_kin, tot_momentum(3), tot_momentum_norm
    character(len=256) :: traj_xyzfile, properties_outfile
    integer :: istep, icartesian,i, dot, current, previous, new
    logical :: suppress_flag = .true., debug_print_all_matrices = .false.

    ! for intuitive storing of position data (since Verlet requires storing previous positions)
    current = 1
    new = 2

    ! INITIALIZATION
    istep = 0
    total_displacement(:) = 0
    input_positions(:,:) = positions(:,:) !x(t=0)
    positions_list(current,:,:) = positions(:,:) !x(t=0)
    acceleration_list(:,:,:) = 0

    call get_energy_gradient(positions_list(current,:,:),tot_pot,forces, gradnorm, suppress_flag) !F(t=0), E_pot(t=0)
    
    call init_v(velocities(:,:)) !v(t=-1)
    do icartesian = 1,3
        acceleration_list(current,:,icartesian) = 1e-4 * forces(:,icartesian) / mweights(:)
        !acceleration in Å/fs^2                    in kJ/mol/Å           in g/mol;
    end do

    ! PREPARE FILE THAT TRACKS PROPERTIES
    dot = index(xyzfile, ".", back=.true.)     ! find last "." in xyzfile name
    properties_outfile = xyzfile(:dot-1) // ".properties.txt"
    open(97, file=properties_outfile, status='replace', action='write')
    write(97,"(A10,7(A20))") "istep", "E_tot", "E_kin","E_pot", "F_norm", "Temp", "Pressure", "COM_momentum"
    write(97,"(A10,7(A20))") "none","kJ/mol", "kJ/mol","kJ/mol", "kJ/(Åmol)", "K", "Pa", "gÅ/fs"

    if (md_debug) then
        write(*,"(/A,/A,/A)") "In propagation","---------------------------------------","Initialization:"
        call recprt2("forces",atomnames,forces,n_atoms)
        write(*,"(A,*(/,F10.6))") "masses: [g/mol]", mweights
       !call recprt2("acceleration = forces / masses",atomnames,acceleration,n_atoms)
        write(*,"(A)") "Start performing steps ..."
    end if

    ! PREPARE TRAJECTORY FILE: traj_xyzfile
    traj_xyzfile = xyzfile(:dot-1) // ".traj" // xyzfile(dot:)
    open(98, file=traj_xyzfile, status='replace', action='write')
    write(98,*) n_atoms
    write(98,"(A,F6.2,A)") "atomic positions at t = ",istep * md_ts, " fs"
    do i=1, size(positions,1), 1
        write(98,FMT='(A3,3(2X,F15.8))') atomnames(i), positions(i,1:3)
    end do

    if (md_ensemble == "NVT") then
        write(*,*) "Turning on the Bussi thermostat"
    elseif (md_ensemble == "NPT") then
        write(*,*) "Turning on the Bussi thermostat"
        write(*,*) "Turning on the Berendsend barostat"
    end if

    if (md_debug) then
        write(*,"(/A,I5)") "Initial quantities at step ",istep
        call recprt2("r(t)",atomnames,positions_list(current,:,:),n_atoms)
    end if

    do while (istep<md_nsteps)
        istep = istep +1

        ! UPDATE POSITIONS
        call velocity_verlet_position(positions_list(current,:,:), velocities(:,:),acceleration_list(current,:,:), &
                            positions_list(new,:,:))

        ! CALCULATE NEW FORCES / ACCELERATION AT NEW POSITION and update velocities 
        call get_energy_gradient(positions_list(new,:,:),tot_pot,forces, gradnorm, suppress_flag)
        do icartesian = 1,3
            acceleration_list(new,:,icartesian) = 1e-4 * forces(:,icartesian) / mweights(:)
            !acceleration in Å/fs^2                    in kJ/mol/Å           in g/mol;
        end do

        call velocity_verlet_velocity(acceleration_list(current,:,:), acceleration_list(new,:,:), velocities)

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

        ! Apply thermostat/barostat constraints
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
        open(97, file=properties_outfile, status='old', action='write')
        ! this prints info from the previous step
        write(97,"(I8,2x,7(F20.8))") istep-1,E_kin+tot_pot,E_kin,tot_pot, gradnorm, instant_temp, pressure, &
                                            tot_momentum_norm/avogad
        if (debug_print_all_matrices) then
            write(*,"(/A,I5)") "New quantities at step ",istep
            call recprt2("r(t) = positions_list(current,:,:) [Å]",atomnames,positions_list(current,:,:),n_atoms)
            call recprt2("F(t) = forces(:,:) [kJ/mol]",atomnames,forces(:,:),n_atoms)
            call recprt2("a(t) = acceleration(:,:) [Å/(fs)^2]",atomnames,acceleration_list(new,:,:),n_atoms)
            call recprt2("v(t) = velocities(:,:) [Å/fs]",atomnames,velocities(:,:),n_atoms)
        end if

        ! TRACK DISPLACEMENT OF THE ATOMS
        call displacement_vec(positions_list(new,:,:),positions_list(current,:,:),displacement,n_atoms,atomnames)
        total_displacement(:) = total_displacement(:) + displacement(:)

        ! WRITE TRAJECTORY FILE
        open(98, file=traj_xyzfile, status='old', action='write')
        write(98,*) n_atoms
        write(98,"(A,F6.2,A)") "atomic positions at t = ",istep * md_ts, " fs"
        do i=1, size(positions,1), 1
            write(98,FMT='(A3,3(2X,F15.8))') atomnames(i), positions_list(new,i,:)
        end do

        ! PREPARE NEXT STEP
        positions_list(current,:,:) = positions_list(new,:,:)
        acceleration_list(current,:,:) = acceleration_list(new,:,:)

    end do

    displacement(:) = 0

    if (md_debug) then
        write(*,"(/A)") "Throughout the simulation, the atoms displaced: (no MSD, but initial vs final coords)"
        call displacement_vec(positions_list(current,:,:),input_positions, displacement,n_atoms,atomnames)
        write(*,*) "Displacements (summed all steps)"
        do i = 1, n_atoms
            write(*,"(I3,1x,A3,1x,  F16.12,1x,A)") i,atomnames(i),displacement(i),"Å"
        end do

        write(*,"(/A)") "Displacements (summed all steps)"
        do i = 1, n_atoms
            write(*,"(I3,1x,A3,1x,F16.12,1x,A)") i,atomnames(i),total_displacement(i),"Å"
        end do
    end if

    write(*,"(/A,//A)") "================================================================","Simulation finished"
    write(*,"(/A,F10.2,A)") "Total simulation time", md_ts * istep, " fs"
    write(*,"(/A,A)") "Wrote properties to ", properties_outfile
    write(*,"(A,A)") "Wrote trajectory to ", traj_xyzfile
    write(*,"(     A)") "================================================================"

end subroutine simulation_vel_verlet

end module simulation_mod
