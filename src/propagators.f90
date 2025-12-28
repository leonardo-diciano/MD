module propagators
    contains

    subroutine Verlet(positions_previous,positions_current,acceleration, n_atoms, positions_new, velocities)
        use definitions, only: wp
        use parser_mod, only: md_ts
        implicit none
        real(kind=wp), intent(in) :: positions_previous(n_atoms,3), positions_current(n_atoms,3),acceleration(n_atoms,3)
        integer, intent(in) :: n_atoms
        real(kind=wp), intent(out) :: positions_new(n_atoms,3), velocities(n_atoms,3)

        ! UPDATE POSITIONS
        positions_new(:,:) = 2 * positions_current(:,:) - positions_previous(:,:) + md_ts**2 * acceleration(:,:) !x(t+1)

        !UPDATE VELOCITIES
        velocities(:,:) = (positions_new(:,:) - positions_previous(:,:))/(2*md_ts) !v(t)
    end subroutine Verlet


end module propagators
