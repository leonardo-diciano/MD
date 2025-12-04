! Forces and Potential Energy calculator 
! ***NOT TESTED, FINAL CHECK NEEDED***
! Author: Leonardo Di Ciano (2025)

module force_field_mod

contains

function cross_product(a, b) result(cross)
  real, dimension(3) :: cross
  real, dimension(3), intent(in) :: a, b

  cross(1) = a(2) * b(3) - a(3) * b(2)
  cross(2) = a(3) * b(1) - a(1) * b(3)
  cross(3) = a(1) * b(2) - a(2) * b(1)
end function

subroutine force_field_calc(n_atoms,n_bonds,n_angles,n_impdie,n_torsions,positions,bond_params,angle_params,&
            impdihedrals_params,tors_params,lj_params,resp_charges,tot_pot,forces)

implicit none

integer, intent(in) :: n_atoms,n_bonds,n_angles,n_torsions,n_impdie
real, intent(in) :: positions(n_atoms,3),bond_params(n_bonds,4),angle_params(n_angles,5),&
                    impdihedrals_params(n_impdie,8),tors_params(n_torsions,8),lj_params(n_atoms,3),resp_charges(n_atoms,2)

real, intent(out) :: tot_pot
real, allocatable, intent(out) :: forces(:,:)
real :: pi, kcal_to_kJ, charge_to_kJ_mol, bond_pot, distance, angle_pot, angle , die_pot, imp_die_pot, coulomb_pot, lj_pot, pot_14
real :: d12(3), d23(3), d34(3), f_magnitude, epsilon, sigma, f1(3), f3(3), f2(3), f4(3), d12_norm, d23_norm
integer :: i, j, k, row, a1, a2, a3, a4, pair(2),one_bond_list(n_bonds,2), two_bonds_list(n_angles,2)  
integer, allocatable :: dummy_pairs(:,:), non_bonded_pairs(:,:), three_bonds_list(:,:) 
real :: a(3),b(3),a_norm,b_norm, dihedral, cap_A(3), cap_B(3), cap_C(3)
logical :: is_bonded

write(*,*) "Entering FF "

! Useful costants and conversion factors
pi = 3.14159265
kcal_to_kJ= 4.184
charge_to_kJ_mol =  332.05 / kcal_to_kJ



allocate(forces(n_atoms,3))
forces(:,:) = 0         ! initialize the force vector

! Bonds term
bond_pot=0  ! initialize the potential
do i=1, n_bonds, 1

    ! Get atoms' indexes 
    a1 = int(bond_params(i,1))
    a2 = int(bond_params(i,2))
    one_bond_list(i,:)=[a1, a2]
    ! Calculate distance between a1 and a2 
    distance = SQRT(dot_product(positions(a1,1:3) - positions(a2,1:3), positions(a1,1:3) - positions(a2,1:3)))

    ! Calculate the force magnitude 
                    !force constant                            equilibrium distance
    f_magnitude = - bond_params(i,3) * kcal_to_kJ * (distance - bond_params(i,4))  

    ! Calculate the force on each atom and update the force vector
    forces(a1,1:3) = forces(a1,1:3) + f_magnitude * ( (positions(a1,1:3) - positions(a2,1:3)) / distance)
    forces(a2,1:3) = forces(a2,1:3) - f_magnitude * ( (positions(a1,1:3) - positions(a2,1:3)) / distance)
    
    ! Calculate the bond contributions to potential energy
    !                           force constant                      equilibrium distance
    bond_pot = bond_pot + bond_params(i,3) * kcal_to_kJ * (distance - bond_params(i,4))**2
end do


! Angles term
angle_pot=0     ! initialize the potential
do i=1, n_angles, 1

    ! Get atoms' indexes 
    a1 = int(angle_params(i,1))
    a2 = int(angle_params(i,2))
    a3 = int(angle_params(i,3))
    two_bonds_list(i,:)=[a1, a3]
    ! Calculate the bond/distance vectors and their norm
    d12 = positions(a1,1:3) - positions(a2,1:3)
    d23 = positions(a3,1:3) - positions(a2,1:3)
    d12_norm = SQRT(dot_product(d12,d12))
    d23_norm = SQRT(dot_product(d23,d23))

    ! Calculate the angle (in radians)
    angle= acos(dot_product(d12,d23) / (d12_norm * d23_norm)) 

    ! Calculate force magnitude
    !                    angle force constant                   equilibrium angle in radians
    f_magnitude = -2 * angle_params(i,4) * kcal_to_kJ * (angle - (angle_params(i,5) * pi / 180))

    ! Calculate the force acting on atom 1 and 3
    f1 = f_magnitude * (1 / ( sin(angle) * d12_norm )) * &
                    (cos(angle) * (d12 / d12_norm) - (d23 / d23_norm ))
                    
    f3 = f_magnitude * (1 / ( sin(angle) * d23_norm )) * &
                    (cos(angle) * (d23 / d23_norm) - (d12 / d12_norm))

    ! Update the force vector
    forces(a1,1:3) = forces(a1,1:3) + f1
    forces(a2,1:3) = forces(a2,1:3) - (f1+f3)
    forces(a3,1:3) = forces(a3,1:3) + f3

    ! Calculate the angle contributions to potential energy
    !                          angle force constant                     equilibrium angle in radians
    angle_pot = angle_pot + ( angle_params(i,4) * kcal_to_kJ * (angle - (angle_params(i,5) * pi / 180 ))**2 )
end do


! Dihedrals term
k=1
die_pot = 0     ! initialize the potential
allocate(three_bonds_list(2 * n_torsions,2))
three_bonds_list(:,:) = 0   ! Initialize the 1-4 interaction atoms list
do i=1, n_torsions, 1

    ! Get atoms' indexes
    a1 = int(tors_params(i,1))
    a2 = int(tors_params(i,2))
    a3 = int(tors_params(i,3))
    a4 = int(tors_params(i,4))

    ! Add atoms 1 and 4 in the list for 1-4 non-bonded interactions
    if (tors_params(i,8) > 0) then  ! negative multiplicity are given to repeated terms by the FF parameters
        three_bonds_list(k,:) = [a1, a4]
        k = k + 1
    end if

    ! Calculate bond/distance vectors
    d12 = positions(a1,1:3) - positions(a2,1:3)
    d23 = positions(a2,1:3) - positions(a3,1:3)
    d34 = positions(a4,1:3) - positions(a3,1:3)

    ! Define a (and b) as cross product and its norm a_norm (and b_norm)
    a = cross_product(d12,d23)          
    a_norm = SQRT(dot_product(a,a))
    b = cross_product(d23,d34)
    b_norm = SQRT(dot_product(b,b))

    ! Calculate the dihedral angle
    dihedral= acos( - dot_product(a,b) / (a_norm * b_norm))

    ! Define A, B and C terms
    cap_A = (b / (a_norm * b_norm)) - ((dot_product(a,b) * a) / (a_norm**3 * b_norm))
    cap_B = (a / (a_norm * b_norm)) - ((dot_product(a,b) * b) / (a_norm * b_norm**3 ))
    cap_C = cross_product(cap_A,d12) + cross_product(cap_B,d34)

    ! Calculate the force magnitude * ( 1 / sin(phi))
                ! torsional barrier                 divider            periodicity
    f_magnitude = tors_params(i,6) * kcal_to_kJ / tors_params(i,5) * tors_params(i,8) * &
                       !  periodicity                  phase in radians
                    sin(tors_params(i,8) * dihedral - tors_params(i,7) * pi / 180 ) / sin(dihedral)

    ! Calculate the force on each atom
    f1 = - f_magnitude * cross_product(cap_A,d23)
    f4 = f_magnitude * cross_product(cap_B,d23)
    f2 = - f1 - (f_magnitude * cap_C)
    f3 = - f4 + (f_magnitude * cap_C)

    ! Update the force vector
    forces(a1,1:3 )= forces(a1,1:3) + f1
    forces(a2,1:3) = forces(a2,1:3) + f2
    forces(a3,1:3) = forces(a3,1:3) + f3
    forces(a4,1:3) = forces(a4,1:3) + f4

    ! Calculate the dihedral contributions to potential energy
                    ! torsional barrier                     divider                
    die_pot = die_pot + tors_params(i,6) * kcal_to_kJ / tors_params(i,5) * &
                !            periodicity                   phase in radians
                (1 + cos(tors_params(i,8) * dihedral - ( tors_params(i,7) * pi / 180)))

end do
if (k > 1) then
    allocate(dummy_pairs(k-1,2))
    dummy_pairs = three_bonds_list(1:k-1, :)  ! copy only used rows
    deallocate( three_bonds_list)
    allocate( three_bonds_list(k-1, 2))
    three_bonds_list = dummy_pairs
    deallocate(dummy_pairs)
end if

! Improper dihedrals term
imp_die_pot=0   ! initialize the potential

do i=1, n_impdie, 1

    ! Get atoms' indexes
    a1 = int(impdihedrals_params(i,1))
    a2 = int(impdihedrals_params(i,2))
    a3 = int(impdihedrals_params(i,3))
    a4 = int(impdihedrals_params(i,4))

    ! Calculate bond/distance vectors
    d12 = positions(a1,1:3) - positions(a2,1:3)
    d23 = positions(a2,1:3) - positions(a3,1:3)
    d34 = positions(a4,1:3) - positions(a3,1:3)

    ! Define a (and b) as cross product and its norm a_norm (and b_norm)
    a = cross_product(d12,d23)          
    a_norm = SQRT(dot_product(a,a))
    b = cross_product(d23,d34)
    b_norm = SQRT(dot_product(b,b))

    ! Calculate the improper dihedral angle
    dihedral= acos( - dot_product(a,b) / (a_norm * b_norm))

    ! Define A, B and C terms
    cap_A = (b / (a_norm * b_norm)) - ((dot_product(a,b) * a) / (a_norm**3 * b_norm))
    cap_B = (a / (a_norm * b_norm)) - ((dot_product(a,b) * b) / (a_norm * b_norm**3 ))
    cap_C = cross_product(cap_A,d12) + cross_product(cap_B,d34)

    ! Calculate the force magnitude * ( 1 / sin(phi))
                  ! torsional barrier                           divider            periodicity
    f_magnitude = impdihedrals_params(i,6) * kcal_to_kJ / impdihedrals_params(i,5) * impdihedrals_params(i,8) * &
                         !  periodicity                          phase 
                    sin(impdihedrals_params(i,8) * dihedral - impdihedrals_params(i,7) * pi / 180) / sin(dihedral)

    ! Calculate the force on each atom
    f1 = - f_magnitude * cross_product(cap_A,d23)
    f4 = f_magnitude * cross_product(cap_B,d23)
    f2 = - f1 - (f_magnitude * cap_C)
    f3 = - f4 + (f_magnitude * cap_C)

    ! Update the force vector
    forces(a1,1:3) = forces(a1,1:3) + f1
    forces(a2,1:3) = forces(a2,1:3) + f2
    forces(a3,1:3) = forces(a3,1:3) + f3
    forces(a4,1:3) = forces(a4,1:3) + f4

    ! Calculate the improper dihedral contributions to potential energy
                                     ! torsional barrier                    divider                                            
    imp_die_pot = imp_die_pot + impdihedrals_params(i,6) * kcal_to_kJ / impdihedrals_params(i,5) * &
                !          periodicity                              phase
                (1 + cos(impdihedrals_params(i,8) * dihedral - impdihedrals_params(i,7) * pi / 180))

end do
write(*,*) "Improper dihedrals term done: ", imp_die_pot


! Create list of atoms for non-bonded interactions (>3 bonds distance)
k=1
allocate(non_bonded_pairs(n_atoms*(n_atoms-1)/2, 2))
do i=1, n_atoms, 1
    do j=i+1, n_atoms, 1
        pair=[i,j]
        
        is_bonded = .false.
        do row = 1, size(one_bond_list,1)   ! Exclude if directly bonded
            if ( (one_bond_list(row,1) == i .and. one_bond_list(row,2) == j) .or. &
                (one_bond_list(row,1) == j .and. one_bond_list(row,2) == i) ) then
                is_bonded = .true.
            endif
        end do
        if (is_bonded) cycle
        
        do row = 1, size(two_bonds_list,1)   ! Exclude two bonds distance
            if ( (two_bonds_list(row,1) == i .and. two_bonds_list(row,2) == j) .or. &
                (two_bonds_list(row,1) == j .and. two_bonds_list(row,2) == i) ) then
                is_bonded = .true.
            endif
        enddo
        if (is_bonded) cycle

        do row = 1, size(three_bonds_list,1)    ! Exclude three bonds distance atoms (1-4 separately)
            if ( (three_bonds_list(row,1) == i .and. three_bonds_list(row,2) == j) .or. &
                (three_bonds_list(row,1) == j .and. three_bonds_list(row,2) == i) ) then
                is_bonded = .true.
            endif
        enddo

        if (is_bonded) then
            cycle
        else                        ! If all checks are passed, assign to non bonded list
            non_bonded_pairs(k,1) = i
            non_bonded_pairs(k,2) = j
            k = k + 1
        end if
    end do
end do
write(*,*) "Non bonded pairs prepared"
if (k > 1) then
    allocate(dummy_pairs(k-1,2))
    dummy_pairs = non_bonded_pairs(1:k-1, :)  ! copy only used rows
    deallocate(non_bonded_pairs)
    allocate(non_bonded_pairs(k-1, 2))
    non_bonded_pairs = dummy_pairs
    deallocate(dummy_pairs)
end if
write(*,*) "Non bonded pairs reshaped"

!LJ term (> 1-4 interactions)
lj_pot =  0    ! Initialize the potential
do i=1, size(non_bonded_pairs,dim=1), 1      ! Iterate over the number of non-bonded pairs

    ! Get atoms' indexes
    a1 = non_bonded_pairs(i,1)
    a2 = non_bonded_pairs(i,2)
    write(*,*) a1 , a2

    ! Calculate the distance between the atoms
    distance = SQRT(dot_product(positions(a1,1:3)-positions(a2,1:3),positions(a1,1:3)-positions(a2,1:3)))

    if (distance < 6 ) then
        ! Calculate the parameters following Lorentz/Berthelot mixing rules
                        ! epsilon a1       epsilon a2
        epsilon = SQRT(lj_params(a1,3) * lj_params(a2,3)) * kcal_to_kJ
                    ! sigma a1           sigma a2
        sigma = 0.5 * (lj_params(a1,2) + lj_params(a2,2))

        ! Calculate the force magnitude
        f_magnitude = 24 * epsilon * ( 2 * (sigma / distance)**12 - (sigma / distance)**6 )

        ! Calculate the force on each atom and update the force vector
        forces(a1,1:3) = forces(a1,1:3) + f_magnitude * ( (positions(a1,1:3) - positions(a1,1:3)) / distance**2 )
        forces(a2,1:3) = forces(a2,1:3) - f_magnitude * ( (positions(a2,1:3) - positions(a2,1:3)) / distance**2 )

        ! Calculate LJ contributions to potential energy
        lj_pot = lj_pot + (4 * epsilon * ( (sigma / distance)**12 - (sigma / distance)**6 ))
        write(*,*) "LJ", i, distance, epsilon, sigma, 4 * epsilon * ( (sigma / distance)**12 - (sigma / distance)**6 ), lj_pot
    else
        cycle
    endif

end do
write(*,*) "LJ done", lj_pot

!Coulomb term 
coulomb_pot = 0   ! initialize the potential
do i=1, k-1, 1

    ! Get atoms' indexes
    a1=non_bonded_pairs(i,1)
    a2=non_bonded_pairs(i,2)

    ! Calculate the distance between the atoms
    distance = SQRT(dot_product(positions(a1,1:3)-positions(a2,1:3),positions(a1,1:3)-positions(a2,1:3)))
    
    ! Calculate the force magnitude
                !       charge on a1        charge on a2
    f_magnitude = - ((resp_charges(a1,2) * resp_charges(a2,2) * charge_to_kJ_mol ) / (distance**2))  

    ! Calculate the force on each atom and update the force vector
    forces(a1,1:3) = forces(a1,1:3) + f_magnitude * ( (positions(a1,1:3) - positions(a1,1:3)) / distance)
    forces(a2,1:3) = forces(a2,1:3) - f_magnitude * ( (positions(a2,1:3) - positions(a2,1:3)) / distance)
    
    ! Calculate Coulombic contributions to the potential energy
    coulomb_pot = coulomb_pot + ((resp_charges(a1,2) * resp_charges(a2,2) *  charge_to_kJ_mol ) / ( distance))

end do
write(*,*) "Coulomb done", coulomb_pot

! 1 - 4 interactions term
pot_14 = 0      ! initialize the potential
do i=1, size(three_bonds_list,dim=1) , 1

    ! Get atoms' indexes
    a1=three_bonds_list(i,1)
    a2=three_bonds_list(i,2)

    ! Calculate the distance between the atoms
    distance = SQRT(dot_product(positions(a1,1:3)-positions(a2,1:3),positions(a1,1:3)-positions(a2,1:3)))

    ! LJ term  (scaled down by 2)

    if (distance < 6 ) then
        ! Calculate the parameters following Lorentz/Berthelot mixing rules
                        ! epsilon a1       epsilon a2
        epsilon = SQRT(lj_params(a1,3) * lj_params(a2,3)) * kcal_to_kJ
                    ! sigma a1           sigma a2
        sigma = 0.5 * (lj_params(a1,2) + lj_params(a2,2))

        ! Calculate the force magnitude
        f_magnitude = 12 * epsilon * ( 2 * (sigma / distance)**12 - (sigma / distance)**6 )

        ! Calculate the force on each atom and update the force vector
        forces(a1,1:3) = forces(a1,1:3) + f_magnitude * ( (positions(a1,1:3) - positions(a1,1:3)) / distance**2 )
        forces(a2,1:3) = forces(a2,1:3) - f_magnitude * ( (positions(a2,1:3) - positions(a2,1:3)) / distance**2 )

        ! Calculate LJ contributions to potential energy
        lj_pot = lj_pot + 2 * epsilon * ( (sigma / distance)**12 - (sigma / distance)**6 )
        ! Calculate 1-4 interaction contributions to potential energy
        pot_14 = pot_14 + 2 * epsilon * ( (sigma / distance)**12 - (sigma / distance)**6 )
    end if
    write(*,*) "1-4", distance, pot_14
    ! Coulomb term (scaled down by 1.2)

    ! Calculate the force magnitude
                    !   charge on a1        charge on a2
    f_magnitude = - ((resp_charges(a1,2) * resp_charges(a2,2) * charge_to_kJ_mol ) / ( distance**2))  / 1.2 

    ! Calculate the force on each atom and update the force vector
    forces(a1,1:3) = forces(a1,1:3) + f_magnitude * ( (positions(a1,1:3) - positions(a1,1:3)) / distance) 
    forces(a2,1:3) = forces(a2,1:3) - f_magnitude * ( (positions(a2,1:3) - positions(a2,1:3)) / distance)
    
    ! Calculate Coulombic contributions to the potential energy
    coulomb_pot = coulomb_pot + ((resp_charges(a1,2) * resp_charges(a2,2) *  charge_to_kJ_mol ) / ( distance)) /1.2 

    ! Calculate 1-4 interaction contributions to potential energy
    pot_14 = pot_14 + ((resp_charges(a1,2) * resp_charges(a2,2) *  charge_to_kJ_mol ) / ( distance)) /1.2 
    write(*,*) "1-4", distance, pot_14
end do
write(*,*) "1-4 done", pot_14 

! Calculate the total potential energy
tot_pot = bond_pot + angle_pot + die_pot + imp_die_pot + lj_pot + coulomb_pot
end subroutine

end module 
