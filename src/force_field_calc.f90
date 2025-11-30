! Forces calculator subroutine
!

subroutine force_field_calc(n_atoms,n_bonds,n_angles,n_impdie,n_torsions,mweights,positions,bond_params,angle_params,&
            impdihedrals_params,tors_params,lj_params,resp_charges,tot_pot,forces)

implicit none

integer, intent(in) :: n_atoms,n_bonds,n_angles,n_torsions,n_impdie
real, intent(in) :: mweights(n_atoms),positions(n_atoms,4),bond_params(n_bonds,4),angle_params(n_angles,5),&
                    impdihedrals_params(n_impdie,8),tors_params(n_torsions,8),lj_params(n_atoms,3),resp_charges(n_atoms,2)

real, intent(out) :: tot_pot, forces(n_atoms,3)
real :: pi, epsilon_0, bond_pot, distance, angle_pot, angle , dihedral, die_pot, imp_die_pot, coulomb_pot, lj_pot
real :: AB(3), BC(3), CD(3), f_magnitude, epsilon, sigma,  angle_0, f1, f3
integer :: i, j, a1, a2, a3, a4

forces(:,:) = 0.0

! Bonds/Stretching term
bond_pot=0
do i=1, n_bonds, 1

    a1 = int(bond_params(i,1))
    a2 = int(bond_params(i,2))

    distance = SQRT(dot_product(positions(a1,2:4) - positions(a2,2:4), positions(a1,2:4) - positions(a2,2:4)))
    f_magnitude = - bond_params(i,4) * (distance - bond_params(i,3))  

    forces(a1,1:3) = forces(a1,1:3) + f_magnitude * ( (positions(a1,2:4) - positions(a2,2:4)) / distance)
    forces(a2,1:3) = forces(a2,1:3) - f_magnitude * ( (positions(a1,2:4) - positions(a2,2:4)) / distance)
    
    bond_pot = bond_pot + 0.5 * bond_params(i,4) * (distance - bond_params(i,3))**2
end do

! Angles term
pi = 3.14
angle_pot=0
do i=1, n_angles, 1

    a1 = int(angle_params(i,1))
    a2 = int(angle_params(i,2))
    a3 = int(angle_params(i,3))

    AB = positions(a1,2:4) - positions(a2,2:4)
    BC = positions(a3,2:4) - positions(a2,2:4)

    !  careful on radians conversion
    angle= acos(dot_product(AB,BC) / (SQRT(dot_product(AB,AB)) * SQRT(dot_product(BC,BC)))) * 180 / pi
    angle_0=angle_params(i,4)*180/pi
    f_magnitude=-2*angle_params(i,5)*(angle - angle_0)
    f1 = f_magnitude * (1 / ( sin(angle) * SQRT(dot_product(AB,AB)) )) * &
                    (cos(angle) * (AB/SQRT(dot_product(AB,AB))) - (BC/SQRT(dot_product(BC,BC))))
    f3 = f_magnitude * (1 / ( sin(angle) * SQRT(dot_product(BC,BC)) )) * &
                    (cos(angle) * (BC/SQRT(dot_product(BC,BC))) - (AB/SQRT(dot_product(AB,AB))))
    forces(a1,1:3)= forces(a1,1:3) + f1
    forces(a2,1:3)= forces(a2,1:3) - (f1+f3)
    forces(a3,1:3)= forces(a3,1:3) + f3

    angle_pot = angle_pot + 0.5 * angle_params(i,5) * (angle - angle_0)**2
end do

! Dihedrals term

do i=1, n_torsions, 1
    a1 = int(tors_params(i,1))
    a2 = int(tors_params(i,2))
    a3 = int(tors_params(i,3))
    a4 = int(tors_params(i,4))

    AB = positions(a1,2:4) - positions(a2,2:4)
    BC = positions(a2,2:4) - positions(a3,2:4)
    CD = positions(a4,2:4) - positions(a3,2:4)

    !  careful on radians conversion
    ! add cross product function
    angle= acos(- (dot_product((AB,BC)))/(SQRT()*SQRT()))

end do

! Improper dihedrals term

!LJ term
lj_pot=0
do i=1, n_atoms, 1
    do j=i+1, n_atoms, 1

        if ( i .NE. j) then

            distance = SQRT(dot_product(positions(i,2:4)-positions(j,2:4),positions(i,2:4)-positions(j,2:4)))
            epsilon= 0.5 * (lj_params(i,3) + lj_params(j,3))
            sigma= SQRT(lj_params(i,2) * lj_params(j,2))
            f_magnitude = 24 * epsilon * ( 2 * (sigma / distance)**12 - (sigma / distance)**6 )

            forces(a1,1:3) = forces(i,1:3) + f_magnitude * ( (positions(i,2:4) - positions(i,2:4)) / distance**2 )
            forces(a2,1:3) = forces(j,1:3) - f_magnitude * ( (positions(j,2:4) - positions(j,2:4)) / distance**2 )
            
            lj_pot = lj_pot + 4 * epsilon * ( (sigma / distance)**12 - (sigma / distance)**6 )
        else
            cycle
        end if
    end do
end do


!Coulomb term 
coulomb_pot=0
epsilon_0=1234
do i=1, n_atoms, 1
    do j=i+1, n_atoms, 1

        if ( i .NE. j) then

            distance = SQRT(dot_product(positions(i,2:4)-positions(j,2:4),positions(i,2:4)-positions(j,2:4)))
            
            f_magnitude = - ((resp_charges(i,2) * resp_charges(j,2)) / (4 * pi * epsilon_0 * distance**2))  

            forces(a1,1:3) = forces(i,1:3) + f_magnitude * ( (positions(i,2:4) - positions(j,2:4)) / distance)
            forces(a2,1:3) = forces(j,1:3) - f_magnitude * ( (positions(i,2:4) - positions(j,2:4)) / distance)
            
            coulomb_pot = coulomb_pot + ((resp_charges(i,2) * resp_charges(j,2)) / (4 * pi * epsilon_0 * distance))
        else
            cycle
        end if
    end do
end do



! Sum all terms to get the total potential
tot_pot = bond_pot + angle_pot + die_pot + imp_die_pot + lj_pot + coulomb_pot
end subroutine

function cross_product()

end function