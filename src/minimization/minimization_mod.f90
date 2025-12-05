module minimization_mod
    implicit none

contains

subroutine get_forces(coordlist,natoms,topofile, F_bonds, U_bonds,gradnorm)
    implicit none

    integer, intent(in) :: natoms
    double precision, intent(in) :: coordlist(3*natoms)
    character(len=200), intent(in) :: topofile
    double precision, intent(out) :: F_bonds(3*natoms)
    double precision, intent(out) :: U_bonds, gradnorm

    character(len=200) :: line, term
    double precision :: x1,x2,y1,y2,z1,z2, b12, b0,k0
    integer :: n1,n2,atom1,atom2
    integer :: i,n, n_term, n_bonds

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! READ THE BONDS FROM THE TOPOLOGY FILE
    ! later: exclude this from the subroutine to not read topofile over and over again;
    ! rather save bond info in an array
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(unit=11, file=topofile)

    ![AtomTypes]
    read(11, *) term, n_term ![AtomTypes]
    !print "(a20,i5)", term, n_term
    n=0
    do while (n<natoms)
        n= n+1
        read(11,*) line
    end do

    ![Bonds]
    read (11,*) term, n_bonds
    !print "(a20,i5)", term, n_bonds

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! from current coordlist(3*natoms) calculate potential and gradient
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    U_bonds = 0
    F_bonds(:) = 0
    do i = 1, n_bonds ! for each bond in the molecule
        ! get info from topo file
        read(11,*) atom1,atom2,b0,k0

        ! index of x coordinate of atom1/2 in coordlist (y coord is at n+1, z is at n+2)
        n1 = 3*(atom1-1)+1
        n2 = 3*(atom2-1)+1

        x1 = coordlist(n1)
        x2 = coordlist(n2)
        y1 = coordlist(n1+1)
        y2 = coordlist(n2+1)
        z1 = coordlist(n1+2)
        z2 = coordlist(n2+2)
        !print *, "x1,y1,z1,x2,y2,z2",x1,y1,z1,x2,y2,z2

        ! calculate bond length from xyz
        b12 = SQRT((x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2)
        print "(a5,f8.1,a7,f6.4,a9,f6.4)", "k0 = ",k0,", b0 = ",b0,", b12 = ",b12

        ! calculate energy
        U_bonds = U_bonds + 0.5 * k0*(b12-b0)**2

        ! calculate force (3N vector)
        F_bonds(n1) = F_bonds(n1) - k0*(b12-b0) * (x1-x2) !-dU/dx1
        F_bonds(n1+1) = F_bonds(n1) - k0*(b12-b0) * (y1-y2) !-dU/dy1
        F_bonds(n1+2) = F_bonds(n1) - k0*(b12-b0) * (z1-z2) !-dU/dz1
        F_bonds(n2) = F_bonds(n2) -k0*(b12-b0) * (x2-x1) !-dU/dx2
        F_bonds(n2+1) = F_bonds(n2) -k0*(b12-b0) * (y2-y1) !-dU/dy2
        F_bonds(n2+2) = F_bonds(n2) -k0*(b12-b0) * (z2-z1) !-dU/dz2

    end do

    gradnorm = SQRT(dot_product(F_bonds(:),F_bonds(:)))
    print*, "b12=",b12
    !print*, "U_bonds = ", U_bonds
    !print*,"F_bonds = ", F_bonds(:)
    print*, "gradnorm = ",gradnorm
    !print*, "coordlist = ", coordlist(:)
    close(11)
end subroutine get_forces


end module minimization_mod
