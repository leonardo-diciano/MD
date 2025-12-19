module lin_alg
contains

subroutine mat_norm(mat,n_atoms,norm)
    ! subroutine to get the norm of a vector that is stored as (n,3) matrix
    use definitions, only: wp
    implicit none
    real(kind=wp),intent(in) :: mat(n_atoms,3)
    integer, intent(in) :: n_atoms
    real(kind=wp),intent(out) :: norm

    integer :: iatom, icartesian

    norm = 0
    do iatom=1,n_atoms
        do icartesian=1,3
        norm = norm + mat(iatom,icartesian)**2
        end do
    end do
    norm = SQRT(norm)


end subroutine mat_norm

subroutine F_tot_norm(forces,n_atoms,norm)
    ! subroutine that first sums up all atomic force vectors to
    ! a whole force vector and then calculates the norm of the latter
    use definitions, only: wp
    use print_mod, only: recprt
    implicit none
    real(kind=wp),intent(in) :: forces(n_atoms,3)
    integer, intent(in) :: n_atoms
    real(kind=wp),intent(out) :: norm

    integer :: iatom, icartesian
    real(kind=wp) :: F_tot(3)

    !call recprt("Forces","(*(F12.6))",forces,n_atoms,3)
    F_tot(:) = 0
    do iatom=1,n_atoms
        F_tot(:) = F_tot(:) + forces(iatom,:)
    end do

    !write(*,*) "F_tot = ",F_tot(:)
    norm = 0
    do icartesian=1,3
        norm = norm + F_tot(icartesian)**2
    end do
    norm = SQRT(norm)

    !write(*,*) "F_norm = ",norm

end subroutine F_tot_norm


subroutine displacement_vec(pos1,pos2,n_atoms,atomnames, displacement)
use definitions, only: wp
implicit none
real(kind=wp),intent(in) :: pos1(n_atoms,3), pos2(n_atoms,3)
character(len=2), intent(in) :: atomnames(:)
integer, intent(in) :: n_atoms
real(kind=wp), intent(out) ::  displacement(n_atoms)

integer :: icartesian
real(kind=wp) :: diff(n_atoms,3)

diff(:,:) = pos2(:,:)-pos1(:,:)
displacement(:) = 0

do icartesian=1,3
    displacement(:) = displacement(:) + diff(:,icartesian)**2
end do

displacement(:) = SQRT(displacement)




end subroutine displacement_vec


end module lin_alg
