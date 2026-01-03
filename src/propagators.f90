module propagators
    contains

    subroutine Verlet(positions_previous,positions_current,acceleration, positions_new, velocities)
        use definitions, only: wp
        use parser_mod, only: md_ts
        use force_field_mod, only: n_atoms

        implicit none
        real(kind=wp), intent(in) :: positions_previous(n_atoms,3), positions_current(n_atoms,3),acceleration(n_atoms,3)
        real(kind=wp), intent(out) :: positions_new(n_atoms,3), velocities(n_atoms,3)

        ! UPDATE POSITIONS
        positions_new(:,:) = 2 * positions_current(:,:) - positions_previous(:,:) + md_ts**2 * acceleration(:,:) !x(t+1)

        !UPDATE VELOCITIES
        velocities(:,:) = (positions_new(:,:) - positions_previous(:,:))/(2*md_ts) !v(t)
    end subroutine Verlet

    subroutine velocity_verlet_position(positions_current, velocities, acceleration, positions_new)
        use definitions, only: wp
        use parser_mod, only: md_ts
        use force_field_mod, only: n_atoms

        implicit none
        real(kind=wp), intent(in) :: positions_current(n_atoms,3),acceleration(n_atoms,3), velocities(n_atoms,3)
        real(kind=wp), intent(out) :: positions_new(n_atoms,3)

        ! UPDATE POSITIONS
        positions_new(:,:) = positions_current(:,:)+ md_ts * velocities(:,:) + md_ts**2 * acceleration(:,:)

    end subroutine

    subroutine velocity_verlet_velocity(old_acceleration, new_acceleration, velocities)
        use definitions, only: wp
        use parser_mod, only: md_ts
        use force_field_mod, only: n_atoms

        implicit none
        real(kind=wp), intent(in) :: old_acceleration(n_atoms,3), new_acceleration(n_atoms,3)
        real(kind=wp), intent(inout) :: velocities(n_atoms,3)

        !UPDATE VELOCITIES
        velocities(:,:) = velocities(:,:) + (md_ts / 2 ) * (old_acceleration(:,:) + new_acceleration(:,:))
    end subroutine

end module propagators
