! module contains
! box definition
! setting periodic boundary conditions
! filling the box with particles

module pbc_mod
contains
subroutine pbc_ctrl_positions(positions)
    ! for now just a 3D cube, where (0,0,0) is in the center; try allowing different shapes later
    use definitions, only: wp
    use parser_mod, only: md_boxlength
    use force_field_mod, only: n_atoms
    implicit none
    real(kind=wp), intent(inout) :: positions(n_atoms,3)
    real(kind=wp) :: displacement, origin(3)
    integer :: icartesian, iatom

    origin(:) = 0.0

    do iatom = 1, n_atoms
        do icartesian = 1, 3
            if (positions(iatom,icartesian) >= 0.5*md_boxlength) then
                write(*,*) "atom outside of box in direction ", icartesian, positions(iatom,icartesian)
                positions(iatom,icartesian) = positions(iatom,icartesian) - md_boxlength
                write(*,*) "corrected to: ",positions(iatom,icartesian)
            else if (positions(iatom,icartesian) <= -0.5*md_boxlength) then
                write(*,*) "atom outside of box in direction ", icartesian, positions(iatom,icartesian)
                positions(iatom,icartesian) = positions(iatom,icartesian) + md_boxlength
                write(*,*) "corrected to: ",positions(iatom,icartesian)
            end if
        end do
    end do
    
end subroutine pbc_ctrl_positions
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine define_box()
    ! for now just a 3D cube; try allowing different shapes later
    use definitions, only: wp
    use parser_mod, only: md_boxlength
    implicit none

    logical:: within, outside
    real(kind = wp) :: centerpoint(3), vertices(8,3), faces(6,4,3), measures(3)
    integer :: sign1, sign2, sign3, vertnr, i

    measures = [md_boxlength, md_boxlength, md_boxlength]
    centerpoint(:) = 0



    write(*,"(/A,/A,/)") "DEFINE THE BOX", "---------------------------"
    ! find the positions of the vertices by expanding around the center using the input measurements
    ! the following code results in 8 vertex definitions, where the first one has all negative coords if the center is at 0

    vertnr = 0
    do sign1 = -1,1,2
        do sign2 = -1,1,2
            do sign3 = -1,1,2
                vertnr = vertnr + 1
                !write(*,*) sign1, sign2, sign3
                vertices(vertnr,1) = centerpoint(1) + sign1 * 0.5 * measures(1)
                vertices(vertnr,2) = centerpoint(2) + sign2 * 0.5 * measures(2)
                vertices(vertnr,3) = centerpoint(3) + sign3 * 0.5 * measures(3)
                write(*,"(A,I3,A,3(F15.3))")"vertex ",vertnr," lies at ", vertices(vertnr,:)
            end do
        end do
    end do

    ! example vertices with measures = 4,4,4
    ! vertex   1 lies at   -2.000  -2.000  -2.000
    ! vertex   2 lies at   -2.000  -2.000   2.000
    ! vertex   3 lies at   -2.000   2.000  -2.000
    ! vertex   4 lies at   -2.000   2.000   2.000
    ! vertex   5 lies at    2.000  -2.000  -2.000
    ! vertex   6 lies at    2.000  -2.000   2.000
    ! vertex   7 lies at    2.000   2.000  -2.000
    ! vertex   8 lies at    2.000   2.000   2.000

    ! facing the xy plane with the z axis going through the desktop plane in standard avogadro
    ! define the cube faces by connecting the correct vertices (we don't need to define edges explicitly here)
    ! left and right (fixed x): vertices 1234 and 5678 respectively
    ! top and bottom (fixed y): 3478 and 1256
    ! front and back (fixed z): 2468 and 1357

    ! left
    faces(1,1,:) = vertices(1,:)
    faces(1,2,:) = vertices(2,:)
    faces(1,3,:) = vertices(3,:)
    faces(1,4,:) = vertices(4,:)

    ! right
    faces(2,1,:) = vertices(5,:)
    faces(2,2,:) = vertices(6,:)
    faces(2,3,:) = vertices(7,:)
    faces(2,4,:) = vertices(8,:)

    ! top
    faces(3,1,:) = vertices(3,:)
    faces(3,2,:) = vertices(4,:)
    faces(3,3,:) = vertices(7,:)
    faces(3,4,:) = vertices(8,:)

    ! bottom
    faces(4,1,:) = vertices(1,:)
    faces(4,2,:) = vertices(2,:)
    faces(4,3,:) = vertices(5,:)
    faces(4,4,:) = vertices(6,:)

    ! front
    faces(5,1,:) = vertices(2,:)
    faces(5,2,:) = vertices(4,:)
    faces(5,3,:) = vertices(6,:)
    faces(5,4,:) = vertices(8,:)

    ! back
    faces(6,1,:) = vertices(1,:)
    faces(6,2,:) = vertices(3,:)
    faces(6,3,:) = vertices(5,:)
    faces(6,4,:) = vertices(7,:)

    do i = 1, 6
        write(*,"(/A,I2)") "face nr ", i
        do vertnr = 1,4
            write(*,"(3(F15.3))") faces(i,vertnr,:)
        end do
    end do

    end subroutine define_box


end module pbc_mod
