! Thermostat and Barostat module with several subroutines:
!   - Bussi/CSVR thermostat
!   - Nosé-Hoover thermostat
! Author: Leonardo Di Ciano (2025)

module ensemble_mod

contains

subroutine bussi_thermostat(new_K, targ_K, ndeg, tau, alpha, Dt)
use definitions, only: wp

implicit none
    
real(kind=wp), intent(in) :: new_K, targ_K, tau, Dt
integer, intent(in) :: ndeg
real(kind=wp), intent(out) :: alpha
real(kind=wp) :: rr, factor

if (tau > Dt) then
    factor=exp(- Dt / tau)
else
    factor = 0
end if
rr = gauss_distrib()
alpha = SQRT( factor + (targ_K/(ndeg*new_K)) * (1-factor) * (rr**2 + sumnoises(ndeg-1)) +&
         2 * SQRT(factor) * SQRT((targ_K / (ndeg *new_K)) * (1-factor)*rr) )
end subroutine

subroutine nosehoover_thermostat(n_atoms,cur_K,targ_T,Dt,Q,friction)
use definitions, only: wp, boltzmann

implicit none
real(kind=wp), intent(in) :: cur_K, targ_T, Dt, Q
integer, intent(in) :: n_atoms
real(kind=wp), intent(inout) :: friction

friction = friction + (Dt/(2*Q))*(cur_K - (((3*n_atoms + 1)/2) * boltzmann * targ_T))

end subroutine

subroutine par_rah_barostat()


end subroutine


function gamma_distrib(ia) result(gamma_dist)
use definitions, only: wp

implicit none
integer, intent(in) :: ia
real(kind=wp) :: gamma_dist
real(kind=wp) :: rand1, rand2, x, am, s, v1, v2, acc_prob
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
            am = ndeg-1
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

end module