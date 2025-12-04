!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program to find the next local minimum
!
! Author: Lila Zapp (2025)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program minimization
    implicit none

    integer, parameter :: maxatoms = 1000
    character(len=2), allocatable :: symbol(:)
    real, allocatable  :: x(:), y(:), z(:), coordlist(:)
    character(len=200) :: line, term,topofile,xyzfile
    integer            :: i,n, natoms,n_term, n_bonds

    !for reading topofile
    real :: b0,k0
    integer :: atom1,atom2, n1,n2

    ! for the minimization
    real,allocatable :: F_bonds(:)
    real :: b12, U_bonds
    real :: x1,x2,y1,y2,z1,z2
    real :: alpha, gradnorm, gradnorm_previous
    real, parameter :: conv_gradnorm=1e-6
    alpha = 0.001

    topofile = "H2_distorted.top"
    xyzfile = "H2_distorted.xyz"


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! READ THE COORDINATES INTO A LIST (x1,y1,z1,x2,y2,z2,...)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(unit=10, file=xyzfile)
    
    ! read first line as natoms
    read(10, *) natoms
    allocate(x(natoms))
    allocate(y(natoms))
    allocate(z(natoms))
    allocate(symbol(natoms))
    allocate(coordlist(3*natoms))
    allocate(F_bonds(3*natoms))

    ! skip second (comment line)
    read(10, *) line

    print "(a20,i3)", "natoms = ", natoms
    n=0
    do while (n<natoms)
        n = n + 1
        read(10,*) symbol(n),x(n),y(n),z(n)
        !print '(I3,2X,A2,3F12.6)', n, symbol(n), x(n), y(n), z(n)
        i = 3*(n-1)+1
        !print *, i,i+1,i+2
        coordlist(i)=x(n)
        coordlist(i+1)=y(n)
        coordlist(i+2)=z(n)
    end do
    print *,"coordinates in a list (x1,y1,z1,x2,y2,z2,...) = ",coordlist

    
    ! xyz file in angstrom, convert to nm:
    coordlist = coordlist*0.1
    print *, "coordinates transformed into nm", coordlist


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! READ THE BONDS FROM THE TOPOLOGY FILE
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    open(unit=11, file=topofile)

    ![AtomTypes]
    read(11, *) term, n_term ![AtomTypes]
    print "(a20,i5)", term, n_term
    n=0
    do while (n<natoms)
        n= n+1
        read(11,*) line
    end do

    ![Bonds]
    read (11,*) term, n_bonds
    print "(a20,i5)", term, n_bonds

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! from current coordlist(3*natoms) calculate potential and gradient
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    U_bonds = 0
    F_bonds(:) = 0
    do i = 1, n_bonds ! for each bond in the molecule
        ! get info from topo file
        read(11,*) atom1,atom2,b0,k0
        print *, "atom",atom1, symbol(atom1), "is bonded to atom", atom2, symbol(atom2)
        
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
        print*,"0.5 * k0*(b12-b0)**2 = ",0.5 * k0*(b12-b0)**2
        
        ! calculate force (3N vector)
        F_bonds(n1) = F_bonds(n1) - k0*(b12-b0) * (x1-x2) !-dU/dx1
        F_bonds(n1+1) = F_bonds(n1) - k0*(b12-b0) * (y1-y2) !-dU/dy1
        F_bonds(n1+2) = F_bonds(n1) - k0*(b12-b0) * (z1-z2) !-dU/dz1
        F_bonds(n2) = F_bonds(n2) -k0*(b12-b0) * (x2-x1) !-dU/dx2
        F_bonds(n2+1) = F_bonds(n2) -k0*(b12-b0) * (y2-y1) !-dU/dy2
        F_bonds(n2+2) = F_bonds(n2) -k0*(b12-b0) * (z2-z1) !-dU/dz2

    end do
    print*, "U_bonds = ", U_bonds
    print*,"F_bonds = ", F_bonds

    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! now that we have the gradient (force vector), 
    ! we can follow it (gradient descent) to the minimum
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    gradnorm_previous = SQRT(dot_product(F_bonds,F_bonds))
    print*,"gradnorm_previous = ", gradnorm_previous
    !do while ABS(gradnorm-gradnorm_previous)>conv_grad
    
    coordlist(:) = coordlist(:) + alpha*F_bonds(:)
    print *, coordlist(:)
    
end program minimization

