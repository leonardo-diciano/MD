module simulation_subroutines
contains


subroutine init_v(positions,velocities, n_atoms, mweights, debug_flag, md_temp)
    use definitions, only: wp, pi, kB, proton_mass
    use print_mod, only: recprt3
    implicit none
    real(kind=wp), intent(in) :: mweights(n_atoms), positions(n_atoms,3), md_temp
    integer, intent(in) :: n_atoms
    logical, intent(in) :: debug_flag
    real(kind=wp), intent(out) :: velocities(n_atoms,3)

    real(kind=wp) :: rand1,rand2, rand_gaussian, v_ix, v_atoms(n_atoms), tot_momentum(3)
    real(kind=wp) :: E_kin, instant_temp
    integer :: iatom, icartesian

    write(*,"(/A,/A)") "in init_v", "-------------------------------"

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
    call get_tot_momentum(velocities,mweights, n_atoms, tot_momentum)

    ! NOW WE RECALCULATE THE VELOCITIES TO MAKE THE TOTAL MOMENTUM ZERO (Leach, p.365)
    do icartesian = 1,3
        velocities(:,icartesian) = velocities(:,icartesian) - tot_momentum(icartesian) / SUM(mweights)
    end do
    write(*,"(/A)") "After adapting the velocities with respect to the momentum and mass:"
    call get_tot_momentum(velocities,mweights, n_atoms, tot_momentum) !THIS SHOULD NOW BE ZERO

    call recprt3("v(t_0) = velocities(:,:) [Å/fs]",velocities(:,:),n_atoms)
    call get_v_atoms(v_atoms,velocities,n_atoms,.true.)

    call get_temperature(velocities, mweights,n_atoms, instant_temp)

    end subroutine init_v

subroutine get_v_atoms(v_atoms,velocities,n_atoms,printopt)
    use definitions, only: wp
    implicit none
    real(kind=wp),intent(out) :: v_atoms(n_atoms)
    logical, intent(in) :: printopt
    integer, intent(in) :: n_atoms
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
        !if (printopt) then
            write(*,"(F16.8,A)") v_atom, ","
            !write(*,"(A,I3,A,F16.8,A)") "v_atom of atom ", iatom, "= ", v_atom, " Å/fs"
        !end if
    end do
    if (printopt) then
        write(*,*) "]"
    end if
end subroutine get_v_atoms

subroutine get_tot_momentum(velocities,mweights, n_atoms, tot_momentum)
    use definitions, only: wp, proton_mass
    implicit none
    real(kind=wp), intent(in) :: velocities(n_atoms,3), mweights(n_atoms)
    integer, intent(in) :: n_atoms
    real(kind=wp), intent(out) :: tot_momentum(3)
    integer :: iatom, icartesian

    tot_momentum(:) = 0

    do iatom=1,n_atoms
        do icartesian=1,3
            tot_momentum(icartesian) = tot_momentum(icartesian) + mweights(iatom) * velocities(iatom,icartesian)
            !g Å/(mol fs)                                            !g/mol                            !Å/fs
        end do
    end do

    write(*,"(/A,3(F10.6),A)") "The total momentum of the system is ", tot_momentum, " g Å/(mol fs)"
end subroutine get_tot_momentum

subroutine get_E_kin(velocities, mweights, n_atoms, E_kin)
    use definitions, only: wp
    implicit none
    real(kind=wp), intent(in) :: velocities(n_atoms,3), mweights(n_atoms)
    integer, intent(in) :: n_atoms
    real(kind=wp), intent(out) :: E_kin

    integer :: icartesian, iatom

    E_kin = 0

    do iatom = 1,n_atoms
        do icartesian = 1, 3
            !E_kin = E_kin + 0.5 * mweights(iatom)*velocities(iatom,icartesian)**2
                                    !g/mol               !Å^2/fs^2
            !E_kin = E_kin + 0.5 * 1e-3 * mweights(iatom)* 1e10 * velocities(iatom,icartesian)**2
            !J/mol                       !kg/mol               ! m^2/s^2
            E_kin = E_kin + 0.5 * 1e4 * mweights(iatom) * velocities(iatom,icartesian)**2
            !kJ/mol                       !kg/mol               ! m^2/s^2
        end do
    end do

    write(*,*) "The kinetic energy is", E_kin, "kJ/mol"
end subroutine get_E_kin

subroutine get_temperature(velocities,mweights,n_atoms, instant_temp)
    use definitions, only: wp, boltzmann
    implicit none
    real(kind=wp), intent(in) ::velocities(n_atoms,3), mweights(n_atoms)
    integer, intent(in) :: n_atoms
    real(kind=wp), intent(out) :: instant_temp !/K

    real(kind=wp) :: E_kin !kJ/mol

    call get_E_kin(velocities, mweights, n_atoms, E_kin)
    instant_temp = 2 * E_kin /(boltzmann*3*n_atoms)
                       !kJ/mol  !kJ/(K mol)

    write(*,*) "The temperature is ", instant_temp, "K"
end subroutine get_temperature

subroutine get_pressure(positions, forces,instant_temp,n_atoms, pressure)
    use definitions, only: wp, boltzmann, avogad
    use parser_mod, only: md_boxlength
    implicit none
    real(kind=wp), intent(in) ::instant_temp, positions(n_atoms,3),forces(n_atoms,3)
    integer, intent(in) :: n_atoms
    real(kind=wp), intent(out) :: pressure

    real(kind=wp) :: volume, avg
    integer :: icartesian, iatom

    volume = md_boxlength * md_boxlength * md_boxlength ! Å^3

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

    write(*,*) "The pressure is ",pressure,"kJ/m^3 = kPa"

end subroutine get_pressure




end module simulation_subroutines
