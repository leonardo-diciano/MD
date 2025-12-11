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


end module lin_alg
