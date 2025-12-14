module propagation
contains

subroutine Verlet_propagator(positions,positions_previous,forces,n_atoms, timestep,nsteps)
    use definitions, only: wp
    implicit none

    real(kind=wp),intent(in) :: forces(n_atoms,3), timestep
    real(kind=wp),intent(inout) :: positions(n_atoms,3),positions_previous(n_atoms,3)
    integer, intent(in) :: n_atoms,nsteps

    integer :: istep
    real(kind=wp) :: acceleration(n_atoms,3), masses(n_atoms,3)


    masses =
    acceleration(:,:) = forces(:,:) / masses(:,:)

    do while (istep<nsteps)
        positions = 2 * positions(:,:) - positions_previous(:,:) + timestep**2 * acceleration



    positions_previous(:,:) = positions(:,:)

end subroutine Verlet_propagator




end module propagation
