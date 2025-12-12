module ensemble_mod

contains

subroutine bussi_thermostat(new_K, targ_K, ndeg, tau, alpha)
use definitions, only: wp
    
implicit none
    
real(kind=wp), intent(in) :: new_K, targ_K, tau, ndeg
real(kind=wp), intent(out) :: alpha


end subroutine

function gamma_distrib(ndeg) result(distrib)

implicit none
real(kind=wp), intent(in) :: ndeg
real(kind=wp), intent(out) :: distrib
real(kind=wp) :: rand1, rand2

do
    call random_number(rand1)
    rand1 = 1 - rand1
    y = -log(rand1)
    t = (y/exp(y-1))**(ndeg-1)
    call random_stduniform(rand2)
    rand2 = 1 - rand2
    if(rand2 <= t) then
        distrib = ndeg * y
        exit
    end if
end do
end function

end module