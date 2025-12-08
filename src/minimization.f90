!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program to find the next local minimum
!
! Author: Lila Zapp (2025)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module minimization_mod

contains

    subroutine minimization(positions,n_atoms,tot_pot,forces, debug_flag)
        implicit none

        logical, intent(in) :: debug_flag
        integer, intent(in) :: n_atoms
        double precision, intent(inout) :: tot_pot
        double precision, intent(inout) :: positions(n_atoms,3), forces(n_atoms,3)

        ! No use for these outside of this subroutine
        double precision :: forces_P1(n_atoms,3), forces_P2(n_atoms,3), forces_P3(n_atoms,3), point1(n_atoms,3),&
                            point2(n_atoms,3),point3(n_atoms,3)
        integer :: iter
        double precision :: tot_pot_P1, tot_pot_P2, tot_pot_P3, gradnorm, gradnorm_previous,a,b,best_step
        double precision, parameter :: conv_gradnorm=1e-6,alpha = 0.01

        call get_energy_gradient(positions,tot_pot,forces)


        gradnorm_previous = SQRT(dot_product(forces,forces))


        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! now that we have the gradient (force vector),
        ! we can follow it (gradient descent) to the minimum
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! line search: calculate energy at three points; first one is already done
        print*, "delta = ", ABS(gradnorm-gradnorm_previous)
        
        iter = 1
        do while (ABS(gradnorm-gradnorm_previous)>conv_gradnorm)
            print*, ""
            print*, "Iteration =", iter
            print*, "delta = ", ABS(gradnorm-gradnorm_previous)
            point1(:,:) = positions(:,:)
            point2(:,:) = positions(:,:) + 0.5 * alpha * forces(:,:) / gradnorm
            point3(:,:) = positions(:,:) + alpha * forces(:,:) / gradnorm

            print*, " "
            print*, "Point1 (eta = 0)"
            !print*, point1(:)
            !print*, "positions = [", positions(:,:)
            CALL get_energy_gradient(point1,tot_pot,forces)
            print*, " "
            print*, "Point2 (eta =",0.5*alpha,")"
            !print*, point2(:)
            !print*, "positions = [", positions(:,:)
            !print*, "alpha * 0.5 * forces(:,:) / gradnorm", alpha * 0.5 * forces(:,:) / gradnorm
            CALL get_energy_gradient(point2,tot_pot,forces)
            print*, " "
            print*, "Point3 (eta =",alpha,")"
            !print*, point3(:)
            !print*, "positions = [", positions(:,:)
            !print*, "alpha * forces(:,:) / gradnorm", alpha * forces(:,:) / gradnorm
            CALL get_energy_gradient(point3,tot_pot,forces)

            !extrapolation
            b = ((tot_pot_P3-tot_pot_P1) * (0.5*alpha)**2 - (tot_pot_P2 - tot_pot_P1)*alpha**2) / (-0.5*alpha**2 *0.5*alpha)
            a = (tot_pot_P3-tot_pot_P1) / alpha**2 - b / alpha
            best_step = - b / (2*a)

            print*,"a = ",a,", b = ",b, ", best_step = ", best_step

            positions(:,:) = positions(:,:) + best_step * forces(:,:) / gradnorm

            print*, positions(:,:)

            gradnorm_previous = gradnorm
            CALL get_energy_gradient(positions,tot_pot,forces)
            iter = iter + 1
            gradnorm = SQRT(dot_product(forces,forces))

        end do

    print*, "==================================================================="
    print*, " HAPPY LANDING "
    print "(a32,i3,a12)","Gradient Descent converged in ", iter,"iterations"
    print*, "==================================================================="

    end subroutine minimization

end module minimization_mod
