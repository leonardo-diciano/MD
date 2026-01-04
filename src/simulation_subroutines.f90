! Complementary subroutines for
! molecular dynamics simulations
!
! Authors: Lila Zapp (2025)

module simulation_subroutines
contains


subroutine init_v(velocities)
    use definitions, only: wp, pi, kB, proton_mass
    use print_mod, only: recprt3
    use parser_mod, only: md_debug, md_temp
    use force_field_mod, only: n_atoms, mweights

    implicit none
    real(kind=wp), intent(out) :: velocities(n_atoms,3)

    real(kind=wp) :: rand1,rand2, rand_gaussian, v_ix, v_atoms(n_atoms), tot_momentum(3)
    integer :: iatom, icartesian

    if (md_debug) then
        write(*,"(/A,/A)") "in init_v", "-------------------------------"
    end if

    do iatom = 1,n_atoms
        do icartesian = 1, 3
            call random_number(rand1)
            call random_number(rand2)
            rand_gaussian = SQRT(-2*LOG(rand1)) * COS(2*pi*rand2)
            v_ix = SQRT(kB*md_temp/(proton_mass*mweights(iatom))) * rand_gaussian

            !write(*,"(/,2(A,1x,F12.8,1x),A,F16.8)") "rand1 =",rand1, "rand2 =", rand2,"rand_gaussian = ",rand_gaussian
            !write(*,"(A,F16.4,A)") "v_ix = ", v_ix, " m/s"

            ! convert to Å/fs -> 1m/s = 10^10/10^15 Å/fs = 10^-5 Å/fs

            v_ix = 1e-5 * v_ix
            !write(*,"(A,F14.8,A)") "v_ix = ", v_ix, " Å/fs"

            velocities(iatom,icartesian) = v_ix
        end do
    end do

    !call recprt3("v(t_0) = velocities(:,:) [Å/fs]",velocities(:,:),n_atoms)
    !call get_v_atoms(v_atoms,velocities,n_atoms,.true.)

    ! THIS IS LIKELY NOT YET ZERO
    call get_tot_momentum(velocities,tot_momentum)

    ! NOW WE RECALCULATE THE VELOCITIES TO MAKE THE TOTAL MOMENTUM ZERO (Leach, p.365)
    do icartesian = 1,3
        velocities(:,icartesian) = velocities(:,icartesian) - tot_momentum(icartesian) / SUM(mweights)
    end do

    if (md_debug) then
        write(*,"(/A)") "After adapting the velocities with respect to the momentum and mass:"
        call get_tot_momentum(velocities, tot_momentum) !THIS SHOULD NOW BE ZERO
        call recprt3("v(t_0) = velocities(:,:) [Å/fs]",velocities(:,:),n_atoms)
    end if
    call get_v_atoms(v_atoms,velocities,md_debug)

    end subroutine init_v

subroutine get_v_atoms(v_atoms,velocities,printopt)
    use definitions, only: wp
    use force_field_mod, only: n_atoms
    
    implicit none
    real(kind=wp),intent(out) :: v_atoms(n_atoms)
    logical, intent(in) :: printopt
    real(kind=wp), intent(in) :: velocities(n_atoms,3)

    integer :: iatom, icartesian
    real(kind=wp) :: v_atom

    v_atoms(:) = 0
    if (printopt) then
        write(*,*) "v_atoms = ["
    end if
    !calculate the norm of the velocity vector on each atom
    do iatom = 1,n_atoms
        v_atom = 0
        do icartesian = 1,3
            v_atom = v_atom + velocities(iatom,icartesian)**2
        end do
        v_atom = SQRT(v_atom)
        v_atoms(iatom) = v_atom
        if (printopt) then
            write(*,"(F16.8,A)") v_atom, ","
            !write(*,"(A,I3,A,F16.8,A)") "v_atom of atom ", iatom, "= ", v_atom, " Å/fs"
        end if
    end do
    if (printopt) then
        write(*,*) "]"
    end if
end subroutine get_v_atoms

subroutine get_tot_momentum(velocities, tot_momentum)
    use definitions, only: wp, proton_mass
    use parser_mod, only: md_debug
    use force_field_mod, only: n_atoms, mweights
    
    implicit none
    real(kind=wp), intent(in) :: velocities(n_atoms,3)
    real(kind=wp), intent(out) :: tot_momentum(3)
    integer :: iatom, icartesian

    tot_momentum(:) = 0

    do iatom=1,n_atoms
        do icartesian=1,3
            tot_momentum(icartesian) = tot_momentum(icartesian) + mweights(iatom) * velocities(iatom,icartesian)
            !g Å/(mol fs)                                            !g/mol                            !Å/fs
        end do
    end do

end subroutine get_tot_momentum

subroutine get_E_kin(velocities, E_kin)
    use definitions, only: wp
    use force_field_mod, only: n_atoms, mweights

    implicit none
    real(kind=wp), intent(in) :: velocities(n_atoms,3)
    real(kind=wp), intent(out) :: E_kin

    integer :: icartesian, iatom

    E_kin = 0

    do iatom = 1,n_atoms
        do icartesian = 1, 3
            E_kin = E_kin + 0.5 * 1e4 * mweights(iatom) * velocities(iatom,icartesian)**2
            !kJ/mol                       !kg/mol               ! m^2/s^2
        end do
    end do

end subroutine get_E_kin

subroutine get_temperature(velocities, instant_temp, E_kin)
    use definitions, only: wp, boltzmann
    use force_field_mod, only: n_atoms, mweights

    implicit none
    real(kind=wp), intent(in) ::velocities(n_atoms,3)
    real(kind=wp), intent(out) :: instant_temp

    real(kind=wp) :: E_kin !kJ/mol

    call get_E_kin(velocities, E_kin)
    instant_temp = 2 * E_kin /(boltzmann*3*n_atoms)
                       !kJ/mol  !kJ/(K mol)

end subroutine get_temperature

subroutine get_pressure(positions, forces,instant_temp, pressure)
    use definitions, only: wp, boltzmann, avogad
    use parser_mod, only: md_boxlength, md_pbc
    use force_field_mod, only: n_atoms
    
    implicit none
    real(kind=wp), intent(in) ::instant_temp, positions(n_atoms,3),forces(n_atoms,3)
    real(kind=wp), intent(out) :: pressure

    real(kind=wp) :: volume, avg, max_x, min_x, max_y, min_y, max_z, min_z
    integer :: icartesian, iatom

    if (md_pbc) then
        volume = md_boxlength * md_boxlength * md_boxlength ! Å^3
    else
        max_x = maxval(positions(:,1))
        min_x = minval(positions(:,1))
        max_y = maxval(positions(:,2))
        min_y = minval(positions(:,2))
        max_z = maxval(positions(:,3))
        min_z = minval(positions(:,3))
        volume = (max_x-min_x) * (max_y-min_y) * (max_z-min_z) ! Å^3
    end if

    avg = 0
    do iatom = 1, n_atoms
        do icartesian = 1, 3
            avg = avg + positions(iatom, icartesian) * forces(iatom, icartesian)
            !kJ/mol          !Å                             !kJ/mol/Å
        end do
    end do
    pressure = 1/(3*volume)*(avg+3*n_atoms*boltzmann*instant_temp)
    !kJ/mol/Å^3      !Å^3              !kJ/mol

    pressure = pressure * 1e30 / avogad

end subroutine get_pressure




end module simulation_subroutines
