module lin_alg
contains

subroutine mat_norm(mat,norm,n_atoms)
    ! subroutine to get the norm of a vector that is stored as (n,3) matrix
    use definitions, only: wp
    !use force_field_mod, only: n_atoms
    implicit none
    integer, intent(in) :: n_atoms
    real(kind=wp),intent(in) :: mat(n_atoms,3)
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

subroutine F_tot_norm(forces,norm,n_atoms)
    ! subroutine that first sums up all atomic force vectors to
    ! a whole force vector and then calculates the norm of the latter
    use definitions, only: wp
    use print_mod, only: recprt
    !use force_field_mod, only: n_atoms

    implicit none
    integer, intent(in) :: n_atoms
    real(kind=wp),intent(in) :: forces(n_atoms,3)
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


subroutine displacement_vec(pos1,pos2, displacement, n_atoms)
use definitions, only: wp
implicit none
real(kind=wp),intent(in) :: pos1(n_atoms,3), pos2(n_atoms,3)
real(kind=wp), intent(out) ::  displacement(n_atoms)
integer, intent(in) :: n_atoms
integer :: icartesian
real(kind=wp) :: diff(n_atoms,3)

diff(:,:) = pos2(:,:)-pos1(:,:)
displacement(:) = 0

do icartesian=1,3
    displacement(:) = displacement(:) + diff(:,icartesian)**2
end do

displacement(:) = SQRT(displacement)

end subroutine displacement_vec

function gamma_distrib(ia) result(gamma_dist)
use definitions, only: wp

implicit none
integer, intent(in) :: ia
real(kind=wp) :: gamma_dist
real(kind=wp) :: rand1, rand2, x, am, s, v1, v2, acc_prob, y
integer :: j

! Check if ia is a good input
if (ia < 1 ) then
    write(*,*) "Error in the input parameter of the gamma distribution"
    stop
! If ia is small, follow uniform distribution
elseif (ia < 6) then
    x=1.
    do j=1, ia
        call random_number(rand1)
        x = x*rand1
    end do
    gamma_dist = -log(x)
else
! Calculate the Gamma distribution function
    rand1=0.1
    acc_prob=0
    do while (rand1>acc_prob)
        x=-0.1
        do while(x<0)
            v1 = 1.1
            v2 = 1.1
            do while ((v1**2 + v2**2)> 1) ! v1 and v2 need to be in the unit circle
                call random_number(rand1)
                call random_number(rand2)
                v1 = 2.* rand1 - 1.
                v2 = 2.* rand2 - 1.
            end do

            y = v2/v1
            am = ia-1
            s=SQRT(2.*am+1.)
            x = s*y+am
        end do
        acc_prob=(1.+y**2)*exp(am*log(x/am)-s*y)
    end do
gamma_dist=x
end if
end function

function gauss_distrib() result(gau_distrib)
use definitions, only: wp

implicit none
real(kind=wp) :: gau_distrib
real(kind=wp), save :: gset
real(kind=wp) :: rand1, rand2, v1, v2, fac, rsq
logical, save :: iset = .false.

if (iset .eqv. .false.) then
    rsq=0.
    do while (rsq > 1 .or. rsq == 0)
        call random_number(rand1)
        call random_number(rand2)
        v1 = 2.* rand1 - 1.
        v2 = 2.* rand2 - 1.
        rsq = v1**2 + v2**2
    end do
    fac=SQRT(-2.*log(rsq)/rsq)
    gset = v1*fac
    gau_distrib = v2*fac
    iset = .true.
else
    gau_distrib = gset
    iset = .false.
endif
end function

function sumnoises(num) result(sumnoise)
use definitions, only: wp

implicit none
integer, intent(in) :: num
real(kind=wp) :: sumnoise
real(kind=wp) :: gamdev, gasdev

if (num == 0) then
    sumnoise = 0.0
else if (num == 1) then
    sumnoise = gauss_distrib()**2
else if (MOD(num,2) == 0) then
    sumnoise = 2.0*gamma_distrib(num/2)
else
    sumnoise = 2.0*gamma_distrib((num-1)/2) + gauss_distrib()**2
end if
end function



end module lin_alg
