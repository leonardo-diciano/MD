module propagators
    contains

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


end module propagators
