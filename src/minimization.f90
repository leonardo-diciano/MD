!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Program to find the next local minimum
!
! Author: Lila Zapp (2025)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module minimization_mod

contains

    subroutine minimization(positions,n_atoms,tot_pot,forces, debug_flag)
        use definitions, only: wp
        use force_field_mod, only: get_energy_gradient
        implicit none

        logical, intent(in) :: debug_flag
        integer, intent(in) :: n_atoms
        real(kind=wp), intent(inout) :: tot_pot
        !real(4), intent(inout) :: positions(n_atoms,3)
        real(4), allocatable, intent(inout) :: positions(:,:),forces(:,:)

        ! No use for these outside of this subroutine
        !real(4) :: forces_P1(n_atoms,3), forces_P2(n_atoms,3), forces_P3(n_atoms,3), positions_P1(n_atoms,3),&
        !                    positions_P2(n_atoms,3),positions_P3(n_atoms,3)
        real(4), allocatable :: forces_P1(:,:), forces_P2(:,:), forces_P3(:,:), positions_P1(:,:),&
                            positions_P2(:,:),positions_P3(:,:)
        integer :: iter
        real(4) :: tot_pot_P1, tot_pot_P2, tot_pot_P3, gradnorm_P1, gradnorm_P2, gradnorm_P3, &
                            gradnorm, gradnorm_previous,a,b,best_step
        real(4), parameter :: conv_gradnorm=1e-6,alpha = 0.01

        call get_energy_gradient(positions,tot_pot,forces, gradnorm)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! now that we have the gradient (force vector),
        ! we can follow it (gradient descent) to the minimum
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! line search: calculate energy at three positions_Ps; first one is already done
        if (debug_flag) then
            write(*,"(A/A)") "In minimization:","----------------"
            write(*,*) "delta = ", ABS(gradnorm-gradnorm_previous)
        end if

        allocate(positions_P1(n_atoms,3))
        allocate(positions_P2(n_atoms,3))
        allocate(positions_P3(n_atoms,3))

        iter = 1
        do while (ABS(gradnorm-gradnorm_previous)>conv_gradnorm)
            positions_P1(:,:) = positions(:,:)
            positions_P2(:,:) = positions(:,:) + 0.5 * alpha * forces(:,:) / gradnorm
            positions_P3(:,:) = positions(:,:) + alpha * forces(:,:) / gradnorm

            CALL get_energy_gradient(positions_P1,tot_pot_P1,forces_P1,gradnorm_P1)
            deallocate(forces_P1)


            CALL get_energy_gradient(positions_P2,tot_pot_P2,forces_P2,gradnorm_P2)
            deallocate(forces_P2)

            CALL get_energy_gradient(positions_P3,tot_pot_P3,forces_P3,gradnorm_P3)
            deallocate(forces_P3)

            !extrapolation
            b = ((tot_pot_P3-tot_pot_P1) * (0.5*alpha)**2 - (tot_pot_P2 - tot_pot_P1)*alpha**2) / (-0.5*alpha**2 *0.5*alpha)
            a = (tot_pot_P3-tot_pot_P1) / alpha**2 - b / alpha
            best_step = - b / (2*a)

            write(*,"(/A,I5,/A)") "In iteration ",iter,"----------------"
            write(*,*) "delta = ", ABS(gradnorm-gradnorm_previous)
            if (debug_flag) then
                write(*,*) ""
                write(*,*) "Iteration =", iter
                write(*,*) "delta = ", ABS(gradnorm-gradnorm_previous)

                !write(*,*) " "
                !write(*,*) "positions_P1 (eta = 0)"
                !write(*,*) positions_P1(:)
                !write(*,*) "positions = [", positions(:,:)

                !write(*,*) " "
                !write(*,*) "positions_P2 (eta =",0.5*alpha,")"
                !write(*,*) positions_P2(:)
                !write(*,*) "positions = [", positions(:,:)
                !write(*,*) "alpha * 0.5 * forces(:,:) / gradnorm", alpha * 0.5 * forces(:,:) / gradnorm

                !write(*,*) " "
                !write(*,*) "positions_P3 (eta =",alpha,")"
                !write(*,*) positions_P3(:)
                !write(*,*) "positions = [", positions(:,:)
                !write(*,*) "alpha * forces(:,:) / gradnorm", alpha * forces(:,:) / gradnorm

                !write(*,*)"a = ",a,", b = ",b, ", best_step = ", best_step
            end if

            positions(:,:) = positions(:,:) + best_step * forces(:,:) / gradnorm

            write(*,*) positions(:,:)

            gradnorm_previous = gradnorm
            CALL get_energy_gradient(positions,tot_pot,forces,gradnorm)
            iter = iter + 1

        end do

        deallocate(positions_P1)
        deallocate(positions_P2)
        deallocate(positions_P3)


        write(*,*) "==================================================================="
        write(*,*) " HAPPY LANDING "
        print "(a32,i3,a12)","Gradient Descent converged in ", iter,"iterations"
        write(*,*) "==================================================================="

    end subroutine minimization

end module minimization_mod
