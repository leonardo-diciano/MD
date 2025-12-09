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
        !real(kind=wp), intent(inout) :: positions(n_atoms,3)
        real(kind=wp), allocatable, intent(inout) :: positions(:,:),forces(:,:)

        ! No use for these outside of this subroutine
        !real(kind=wp) :: forces_P1(n_atoms,3), forces_P2(n_atoms,3), forces_P3(n_atoms,3), positions_P1(n_atoms,3),&
        !                    positions_P2(n_atoms,3),positions_P3(n_atoms,3)
        real(kind=wp), allocatable :: forces_P1(:,:), forces_P2(:,:), forces_P3(:,:), positions_P1(:,:),&
                            positions_P2(:,:),positions_P3(:,:)
        integer :: iter,i
        real(kind=wp) :: tot_pot_P1, tot_pot_P2, tot_pot_P3, gradnorm_P1, gradnorm_P2, gradnorm_P3, &
                            gradnorm, gradnorm_previous,tot_pot_previous,a,b,best_step
        real(kind=wp), parameter :: conv_pot=1e-8, conv_gradnorm=1e-6,alpha = 0.01
        integer, parameter :: maxiter = 300

        logical :: suppress_flag

        call get_energy_gradient(positions,tot_pot,forces, gradnorm, suppress_flag = .true.)

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


        write(*,"(/A/A)") "Energy minimization:","-------------------"
        write(*,"(/A)") "iteration          U_tot              delta U            F_tot              delta F"
        write(*,"(A)") "----------------------------------------------------------------------------------------"

        do while ((ABS(tot_pot-tot_pot_previous)>conv_pot) .and. (iter < maxiter))
        !do while ((ABS(gradnorm-gradnorm_previous)>conv_gradnorm) .or. (ABS(tot_pot-tot_pot_previous)>conv_pot)&
            !.and. (iter < maxiter)) ! if we want both conv criteria to be met
            write(*,"(i5,10x,4(F16.8,3x))") iter, tot_pot, tot_pot-tot_pot_previous, gradnorm,&
                                                                gradnorm-gradnorm_previous

            positions_P1(:,:) = positions(:,:)
            positions_P2(:,:) = positions(:,:) + 0.5 * alpha * forces(:,:) / gradnorm
            positions_P3(:,:) = positions(:,:) + alpha * forces(:,:) / gradnorm

            CALL get_energy_gradient(positions_P1,tot_pot_P1,forces_P1,gradnorm_P1,suppress_flag = .true.)
            deallocate(forces_P1)


            CALL get_energy_gradient(positions_P2,tot_pot_P2,forces_P2,gradnorm_P2,suppress_flag = .true.)
            deallocate(forces_P2)

            CALL get_energy_gradient(positions_P3,tot_pot_P3,forces_P3,gradnorm_P3,suppress_flag = .true.)
            deallocate(forces_P3)

            !extrapolation
            b = ((tot_pot_P3-tot_pot_P1) * (0.5*alpha)**2 - (tot_pot_P2 - tot_pot_P1)*alpha**2) / (-0.5*alpha**2 *0.5*alpha)
            a = (tot_pot_P3-tot_pot_P1) / alpha**2 - b / alpha
            best_step = - b / (2*a)

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

            gradnorm_previous = gradnorm
            tot_pot_previous = tot_pot
            CALL get_energy_gradient(positions,tot_pot,forces,gradnorm, suppress_flag = .true.)
            iter = iter + 1
            !write(*,"(/A)") "iteration           U_tot          F_norm           delta"

        end do

        deallocate(positions_P1)
        deallocate(positions_P2)
        deallocate(positions_P3)


        write(*,"(/A)") "==================================================================="
        write(*,"(a32,i3,a12/)")"Gradient Descent converged in ", iter,"iterations"
        write(*,*) "Cartesian forces over the atoms:"
        do i=1, size(forces,1), 1
            write(*,FMT='(I3,5X,F15.6,2X,F15.6,2X,F15.6)') i, forces(i,1), forces(i,2),forces(i,3)
        end do
        write(*,"(/A,F16.8)") "Resulting potential energy:", tot_pot
        write(*,"(/A)") "New atomic positions:"
        do i=1, size(positions,1), 1
            write(*,FMT='(I3,5X,F15.6,2X,F15.6,2X,F15.6)') i, positions(i,1), positions(i,2),positions(i,3)
        end do
        write(*,*) "==================================================================="

    end subroutine minimization

end module minimization_mod
