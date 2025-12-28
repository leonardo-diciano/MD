module propagation
contains

subroutine Verlet_propagator(positions,mweights,n_atoms, debug_flag, atomnames,xyzfile)!,timestep,nsteps)
    use definitions, only: wp
    use print_mod, only: recprt2
    use lin_alg, only: displacement_vec
    use force_field_mod, only: get_energy_gradient
    use parser_mod, only: md_ts,md_nsteps,md_ensemble,md_barostat,md_thermostat,md_temp
    implicit none

    real(kind=wp),intent(in) ::  mweights(n_atoms) !, timestep
    real(kind=wp),intent(inout) :: positions(n_atoms,3)
    integer, intent(in) :: n_atoms !,nsteps
    logical, intent(in) :: debug_flag
    character(len=2), intent(in) :: atomnames(:)
    character(len=256), intent(in) :: xyzfile


    real(kind=wp) :: displacement(n_atoms), input_positions(n_atoms,3), forces(n_atoms,3), &
                    acceleration(n_atoms,3), total_displacement(n_atoms)
    real(kind=wp) :: positions_list(3,n_atoms,3)
    integer :: istep, icartesian,i, dot, current, previous, new ! ,nsteps
    real(kind=wp) ::  gradnorm, tot_pot, kin_en, v(n_atoms,3) !,timestep 
    logical :: suppress_flag = .true.
    character(len=256) :: traj_xyzfile

    !nsteps = 1000
    !timestep = 1.0 !in fs

    ! for intuitive storing of position data (since Verlet requires storing previous positions)
    previous = 1
    current = 2
    new = 3

    ! INITIALIZATION
    istep = 0
    total_displacement(:) = 0
    input_positions(:,:) = positions(:,:)
    
    positions_list(previous,:,:) = positions(:,:)
    positions_list(current,:,:) = positions(:,:)
    !positions_list(previous,:,:) = 0

    if (debug_flag) then
        write(*,"(/A,/A,/A)") "In propagation","---------------------------------------","Initialization:"
    !    call recprt2("forces",atomnames,forces,n_atoms)
        write(*,"(A,*(/,F10.6))") "masses: [g/mol]", mweights
    !   call recprt2("acceleration = forces / masses",atomnames,acceleration,n_atoms)
        write(*,"(A)") "Start performing steps ..."
    end if

    ! PREPARE TRAJECTORY FILE: traj_xyzfile
    dot = index(xyzfile, ".", back=.true.)     ! find last "." in xyzfile name
    traj_xyzfile = xyzfile(:dot-1) // ".traj" // xyzfile(dot:)
    open(98, file=traj_xyzfile, status='replace', action='write')

    write(98,*) n_atoms
    write(98,"(A,F6.2,A)") "atomic positions at t = ",istep * md_ts, " fs"
    do i=1, size(positions,1), 1
        write(98,FMT='(A3,3(2X,F15.8))') atomnames(i), positions(i,1), positions(i,2),positions(i,3)
    end do


    if (debug_flag) then
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
            acceleration(:,icartesian) = 0.0001 *  forces(:,icartesian) / mweights(:)
            !acceleration in Å/fs^2                    in kJ/mol/Å           in g/mol;
        end do

        !UPDATE POSITIONS
        positions_list(new,:,:) = 0
        positions_list(new,:,:) = 2 * positions_list(current,:,:) - positions_list(previous,:,:) &
                                    + md_ts**2 * acceleration(:,:)

        if (debug_flag) then
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

    write(*,"(/A)") "Throughout the simulation, the atoms displaced: (no MSD, but initial vs final coords)"
    call displacement_vec(positions_list(current,:,:),input_positions,n_atoms,atomnames, displacement)
    write(*,*) "Displacements (summed all steps)"
    do i = 1, n_atoms
        write(*,"(I3,1x,A3,1x,F16.12,1x,A)") i,atomnames(i),displacement(i),"Å"
    end do

    write(*,"(/A)") "Displacements (summed all steps)"
    do i = 1, n_atoms
        write(*,"(I3,1x,A3,1x,F16.12,1x,A)") i,atomnames(i),total_displacement(i),"Å"
    end do

    write(*,"(/A,A)") "Trajectory was written to: ", traj_xyzfile
    write(*,"(A,F10.2,A)") "Total simulation time", md_ts * istep, " fs"

end subroutine Verlet_propagator




end module propagation
