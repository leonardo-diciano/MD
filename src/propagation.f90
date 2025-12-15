module propagation
contains

subroutine Verlet_propagator(positions,positions_previous,mweights,n_atoms, debug_flag, atomnames,xyzfile)!,timestep,nsteps)
    use definitions, only: wp
    use print_mod, only: recprt2
    use lin_alg, only: displacement_vec
    use force_field_mod, only: get_energy_gradient
    implicit none

    real(kind=wp),intent(in) ::  mweights(n_atoms) !, timestep
    real(kind=wp),intent(inout) :: positions(n_atoms,3)
    integer, intent(in) :: n_atoms !,nsteps
    logical, intent(in) :: debug_flag
    character(len=2), intent(in) :: atomnames(:)
    character(len=256), intent(in) :: xyzfile


    real(kind=wp) :: displacement(n_atoms), positions_previous(n_atoms,3), input_positions(n_atoms,3), forces(n_atoms,3), &
                    acceleration(n_atoms,3), total_displacement(n_atoms)
    integer :: istep, icartesian,i,nsteps, dot
    real(kind=wp) :: timestep, gradnorm, tot_pot
    logical :: suppress_flag = .true.
    character(len=256) :: traj_xyzfile


    nsteps = 10
    timestep =1e-5
    istep = 0

    ! INITIALIZATION
    input_positions(:,:) = positions(:,:)

    call get_energy_gradient(positions,tot_pot,forces, gradnorm, suppress_flag)
    total_displacement(:) = 0
    acceleration(:,:) = 0
    do icartesian = 1,3
        acceleration(:,icartesian) = forces(:,icartesian) / mweights(:)
    end do

    if (debug_flag) then
        write(*,"(/A,/A,/A)") "In propagation","---------------------------------------","Initialization:"
        call recprt2("forces",atomnames,forces,n_atoms)
        write(*,"(A,*(/,F10.6))") "masses:", mweights
        call recprt2("acceleration = forces / masses",atomnames,acceleration,n_atoms)
        write(*,"(A)") "Start performing steps ..."
    end if


    ! PREPARE TRAJECTORY FILE: traj_xyzfile
    dot = index(xyzfile, ".", back=.true.)     ! find last "." in xyzfile name
    traj_xyzfile = xyzfile(:dot-1) // ".traj" // xyzfile(dot:)
    open(98, file=traj_xyzfile, status='replace', action='write')

    write(98,*) n_atoms
    write(98,"(A,ES10.8)") "atomic positions at t = ",istep*timestep
    do i=1, size(positions,1), 1
        write(98,FMT='(A3,3(2X,F15.8))') atomnames(i), positions(i,1), positions(i,2),positions(i,3)
    end do


    ! HERE THE STEPS ARE TAKEN
    do while (istep<nsteps)
        istep = istep +1

        !UPDATE POSITIONS
        positions_previous(:,:) = positions(:,:)
        positions(:,:) = 2 * positions(:,:) - positions_previous(:,:) + timestep**2 * acceleration


        if (debug_flag) then
            write(*,"(/A,I5)") "New quantities at step ",istep
            call recprt2("positions",atomnames,positions,n_atoms)
        end if

        !UPDATE FORCES / ACCELERATION
        call get_energy_gradient(positions,tot_pot,forces, gradnorm, suppress_flag)
        do icartesian = 1,3
            acceleration(:,icartesian) = forces(:,icartesian) / mweights(:)
        end do

        if (debug_flag) then
            call recprt2("New forces",atomnames,forces,n_atoms)
            call recprt2("New accelerations",atomnames,acceleration,n_atoms)
            call displacement_vec(positions,positions_previous,n_atoms,atomnames,displacement)
            total_displacement(:) = total_displacement(:) + displacement

        end if

        ! WRITE TRAJECTORY FILE
        open(98, file=traj_xyzfile, status='old', action='write')
        write(98,*) n_atoms
        write(98,*) "atomic positions: at t = ",istep*timestep
        do i=1, size(positions,1), 1
            write(98,FMT='(A3,3(2X,F15.8))') atomnames(i), positions(i,1), positions(i,2),positions(i,3)
        end do

    end do

    write(*,"(/A)") "Throughout the simulation, the atoms displaced: (no MSD, but initial vs final coords)"
    call displacement_vec(positions,input_positions,n_atoms,atomnames, displacement)
    write(*,*) "Displacements (summed all steps)"
    do i = 1, n_atoms
        write(*,"(I3,1x,A3,1x,F16.12,1x,A)") i,atomnames(i),displacement(i),"Å"
    end do

    write(*,"(/A)") "Displacements (summed all steps)"
    do i = 1, n_atoms
        write(*,"(I3,1x,A3,1x,F16.12,1x,A)") i,atomnames(i),total_displacement(i),"Å"
    end do

    write(*,*) "Trajectory was written to: ", traj_xyzfile

end subroutine Verlet_propagator




end module propagation
