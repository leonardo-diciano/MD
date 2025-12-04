!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program to find the next local minimum
!
! Author: Lila Zapp (2025)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


program minimization
    use minimization_mod
    implicit none

    character(len=2), allocatable :: symbol(:)
    double precision, allocatable  :: x(:), y(:), z(:), coordlist(:)
    character(len=200) :: line, term,topofile,xyzfile
    integer            :: iter, i,n, natoms,n_term, n_bonds

    !for reading topofile
    double precision :: b0,k0
    integer :: atom1,atom2, n1,n2

    ! for the minimization
    double precision,allocatable :: F(:),F_P1(:), F_P2(:), F_P3(:), point1(:),point2(:),point3(:)
    double precision :: b12, U,U_P1, U_P2, U_P3
    double precision :: x1,x2,y1,y2,z1,z2
    double precision :: alpha, gradnorm,gradnorm_P1, gradnorm_P2, gradnorm_P3, gradnorm_previous
    double precision, parameter :: conv_gradnorm=1e-6
    double precision :: a,b,best_step

    alpha = 0.01

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
    allocate(F(3*natoms))

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


    call get_forces(coordlist,natoms,topofile, F, U,gradnorm)
    gradnorm_previous = 0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! now that we have the gradient (force vector),
    ! we can follow it (gradient descent) to the minimum
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


    ! line search: calculate energy at three points; first one is already done
    allocate(point1(3*natoms))
    allocate(point2(3*natoms))
    allocate(point3(3*natoms))
    allocate(F_P1(3*natoms))
    allocate(F_P2(3*natoms))
    allocate(F_P3(3*natoms))
    print*, "delta = ", ABS(gradnorm-gradnorm_previous)
    iter = 1
    do while (ABS(gradnorm-gradnorm_previous)>conv_gradnorm)
        print*, ""
        print*, "Iteration =", iter
        print*, "delta = ", ABS(gradnorm-gradnorm_previous)
        point1(:) = coordlist(:)
        point2(:) = coordlist(:) + 0.5 * alpha * F(:) / gradnorm
        point3(:) = coordlist(:) + alpha * F(:) / gradnorm

        print*, " "
        print*, "Point1 (eta = 0)"
        !print*, point1(:)
        !print*, "coordlist = [", coordlist(:)
        call get_forces(point1,natoms,topofile, F_P1, U_P1,gradnorm_P1)
        print*, " "
        print*, "Point2 (eta =",0.5*alpha,")"
        !print*, point2(:)
        !print*, "coordlist = [", coordlist(:)
        !print*, "alpha * 0.5 * F(:) / gradnorm", alpha * 0.5 * F(:) / gradnorm
        call get_forces(point2,natoms,topofile, F_P2, U_P2,gradnorm_P2)
        print*, " "
        print*, "Point3 (eta =",alpha,")"
        !print*, point3(:)
        !print*, "coordlist = [", coordlist(:)
        !print*, "alpha * F(:) / gradnorm", alpha * F(:) / gradnorm
        call get_forces(point3,natoms,topofile, F_P3, U_P3,gradnorm_P3)

        !extrapolation
        b = ((U_P3-U_P1) * (0.5*alpha)**2 - (U_P2 - U_P1)*alpha**2) / (-0.5*alpha**2 *0.5*alpha)
        a = (U_P3-U_P1) / alpha**2 - b / alpha
        best_step = - b / (2*a)

        print*,"a = ",a,", b = ",b, ", best_step = ", best_step

        coordlist(:) = coordlist(:) + best_step * F(:) / gradnorm

        print*, coordlist(:)

        gradnorm_previous = gradnorm
        call get_forces(coordlist,natoms,topofile, F, U,gradnorm)
        iter = iter + 1

    end do

print*, "==================================================================="
print*, " HAPPY LANDING "
print "(a32,i3,a12)","Gradient Descent converged in ", iter,"iterations"
print*, "==================================================================="
    
end program minimization
