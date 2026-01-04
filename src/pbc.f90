! module contains
! box definition
! setting periodic boundary conditions
! filling the box with particles

module pbc_mod
contains
subroutine pbc_ctrl_positions(positions)
    ! for now just a 3D cube, where (0,0,0) is in the center; try allowing different shapes later
    use definitions, only: wp, md_boxlength
    use force_field_mod, only: n_atoms
    implicit none
    real(kind=wp), intent(inout) :: positions(n_atoms,3)
    real(kind=wp) :: displacement, origin(3)
    integer :: icartesian, iatom

    origin(:) = 0.0

    do iatom = 1, n_atoms
        do icartesian = 1, 3
            if (positions(iatom,icartesian) >= 0.5*md_boxlength) then
                !write(*,*) "atom outside of box in direction ", icartesian, positions(iatom,icartesian)
                positions(iatom,icartesian) = positions(iatom,icartesian) - md_boxlength
                !write(*,*) "corrected to: ",positions(iatom,icartesian)
            else if (positions(iatom,icartesian) <= -0.5*md_boxlength) then
                !write(*,*) "atom outside of box in direction ", icartesian, positions(iatom,icartesian)
                positions(iatom,icartesian) = positions(iatom,icartesian) + md_boxlength
                !write(*,*) "corrected to: ",positions(iatom,icartesian)
            end if
        end do
    end do

end subroutine pbc_ctrl_positions

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



end module pbc_mod
