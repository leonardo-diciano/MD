!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subroutine that performs energy minimization of the structure stored in the "positions(n_atoms,3)" array
! the algorithm employed here is a simple steepest descent
! to find the optimal step size in the gradient direction, a line search is done by fitting a parabola
! through the points P1,P2,P3. The points are defined using a constant alpha:
!   P1 = positions; P2 = positions + 0.5 * alpha * grad; P3 = positions + min_alpha * grad
!
! Author: Lila Zapp (2025)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module minimization_mod

contains

    subroutine minimization(positions,n_atoms,tot_pot,forces, debug_flag, xyzfile, atomnames,opt_method)
        use definitions, only: wp
        use print_mod, only: recprt,recprt2,recprt3
        use force_field_mod, only: get_energy_gradient
        use lin_alg, only: mat_norm, displacement_vec
        use parser_mod, only: min_max_iter, min_etol, min_ftol,min_alpha

        implicit none
        logical, intent(in) :: debug_flag
        integer, intent(in) :: n_atoms,opt_method
        character(len=256), intent(in) :: xyzfile
        character(len=2), intent(in) :: atomnames(:)
        real(kind=wp), intent(inout) :: tot_pot
        real(kind=wp), intent(inout) :: forces(n_atoms,3),positions(n_atoms,3)

        real(kind=wp) :: forces_P1(n_atoms,3), forces_P2(n_atoms,3), forces_P3(n_atoms,3)
        real(kind=wp), allocatable :: positions_P1(:,:),positions_P2(:,:),positions_P3(:,:)
        real(kind=wp) :: input_positions(n_atoms,3), forces_previous(n_atoms,3), search_dir(n_atoms,3), &
                        search_dir_previous(n_atoms,3)
        integer :: iter,i, dot, icartesian, iatom
        real(kind=wp) :: tot_pot_P1, tot_pot_P2, tot_pot_P3, gradnorm_P1, gradnorm_P2, gradnorm_P3, &
                            gradnorm, gradnorm_previous,tot_pot_previous,a,b,best_step, dummy_real,&
                            beta_numerator,beta_denominator, beta
        logical :: suppress_flag, converged_pot =.false., converged_grad = .false., converged = .false., conj_grad = .false.
        character(len=256) :: minimized_xyzfile, traj_xyzfile

        real(kind=wp) :: displacement(n_atoms)

        write(*,"(/A/A)") "Energy minimization:","-------------------"
        if (.not. debug_flag) then
            write(*,"(/A)") "iteration          U_tot              delta U            F_tot              delta F"
            write(*,"(A)") "----------------------------------------------------------------------------------------"
        else
            write(*,*) "Gradient Descent: positions(:,:) = positions(:,:) - eta * gradient(:,:)/gradnorm &
                                                         = positions(:,:) + eta * forces(:,:)/F_tot"
            write(*,*) "where the optimal eta in the current gradient vector direction is obtained in a line search:"
            write(*,"(A/A/A)") "Three points (U_P1, U_P2, U_P3) in the direction of the gradient vector are calculated,",&
                       " a 2nd order rational function is fitted; of which the minimum is at the position of   ",&
                      " the optimal eta (best_step). U_P1 is the energhy at initial positions(:,:);",&
                      " U_P2 and P_3 are the energy at positions(:,:) + 0.5 min_alpha * forces(:,:)/F_tot and 1 min_alpha", &
                      "respectively."
            write(*,"(3(A,F10.8,2x),A)") "min_alpha = ", min_alpha, "Å;     0.5*min_alpha = ", 0.5*min_alpha, "Å"
            write(*,"(/A)") "iteration  U_tot (U_P1)  delta U   F_tot        delta F       &
                U_P1   U_P3  [  function to obtain optimal step size    ]  best_step   max displacement"
            write(*,"(A)") "-------------------------------------------------------------------------------------------&
                            ------------------------------------------------"
        end if

        input_positions(:,:) = positions(:,:)

        call get_energy_gradient(positions,tot_pot,forces, gradnorm, suppress_flag = .true.)
        forces_previous = forces(:,:)
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! now that we have the gradient (force vector),
        ! we can follow it (gradient descent) to the minimum
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! line search: calculate energy at three positions_Ps; first one is already done
        allocate(positions_P1(n_atoms,3))
        allocate(positions_P2(n_atoms,3))
        allocate(positions_P3(n_atoms,3))

        iter = 0

        ! Print initial data into iteration table
        if (.not. debug_flag) then
            write(*,"(i5,10x,2(F16.8,13x,A,8x))") iter, tot_pot,"/" , gradnorm,"/"
        else
            write(*,"(i5,2x,2(F12.6,8x,A,5x))") &
                iter, tot_pot, "/", gradnorm,"/"
        end if


        !prepare traj file
        dot = index(xyzfile, ".", back=.true.)     ! find last "." in xyzfile name
        traj_xyzfile = xyzfile(:dot-1) // ".minimization.traj" // xyzfile(dot:)
        open(97, file=traj_xyzfile, status='replace', action='write')

        write(97,*) n_atoms
        write(97,*) "atomic positions at iteration = ",iter
        do i=1, size(positions,1), 1
            write(97,FMT='(A3,3(2X,F15.8))') atomnames(i), positions(i,1), positions(i,2),positions(i,3)
        end do


        ! START ITERATING
        do while (.not. converged)
            iter = iter + 1

            if (opt_method==2 .and. iter==1) then
                write(*,*) "Activating conjugate gradient procedure"
                conj_grad = .false.
            end if

            if (.not. conj_grad) then
                search_dir(:,:) = forces(:,:) / gradnorm !unit vector that is the opposite of gradient vector
            else if (conj_grad) then
                search_dir_previous(:,:) = search_dir(:,:)

                beta = 0
                beta_numerator = 0
                beta_denominator = 0
                ! Polak-Ribiere
                do icartesian=1,3
                    do iatom =1,n_atoms
                        beta_numerator = beta_numerator + forces(iatom,icartesian)*(forces(iatom,icartesian)&
                                        -forces_previous(iatom,icartesian))
                        beta_denominator = beta_denominator + forces_previous(iatom,icartesian)&
                                        *forces_previous(iatom,icartesian)
                    end do
                end do
                beta=beta_numerator/beta_denominator

                search_dir(:,:) = forces(:,:) / gradnorm + beta * search_dir_previous(:,:)
                if (debug_flag) then
                    write(*,*) "beta_denominator =", beta_denominator, "beta_numerator =", beta_numerator, "beta = ",beta
                    call recprt2("forces(:,:)/gradnorm",atomnames,forces(:,:)/gradnorm,n_atoms)
                    call recprt2("search direction",atomnames,search_dir,n_atoms)
                end if
            end if

            ! get coordinates for the three points used in the line search:
            positions_P1(:,:) = positions(:,:)
            positions_P2(:,:) = positions(:,:) + 0.5 * min_alpha * search_dir(:,:)
            positions_P3(:,:) = positions(:,:) + min_alpha * search_dir(:,:)

            ! get energies for the three points used in the line search:
            CALL get_energy_gradient(positions_P1,tot_pot_P1,forces_P1,gradnorm_P1,suppress_flag = .true.)
            CALL get_energy_gradient(positions_P2,tot_pot_P2,forces_P2,gradnorm_P2,suppress_flag = .true.)
            CALL get_energy_gradient(positions_P3,tot_pot_P3,forces_P3,gradnorm_P3,suppress_flag = .true.)

            ! fitting a parabola through the three points with respect to the step size; find x=best_step of minimum
            b = ((tot_pot_P3-tot_pot_P1) * (0.5*min_alpha)**2 - (tot_pot_P2 - tot_pot_P1)*min_alpha**2) / &
                                                                        (-0.5*min_alpha**2 *0.5*min_alpha)
            a = (tot_pot_P3-tot_pot_P1) / min_alpha**2 - b / min_alpha
            best_step = - b / (2*a)

            ! Perform the step along the gradient vector
            positions(:,:) = positions(:,:) + best_step * search_dir(:,:)



            ! uncomment this to see the displacements of all atoms (default: only max value is printed (max displacement))
            !if (debug_flag) then
            !    call recprt2("Forces",atomnames,forces,n_atoms)
            !    call recprt2("Step taken",atomnames,best_step*forces/gradnorm,n_atoms)
            !    call mat_norm(forces/gradnorm, n_atoms,dummy_real)
            !    write(*,*) "norm =", dummy_real !forces/gradnorm is a unit vector
            !end if

            ! save previous values to track progress and convergence
            gradnorm_previous = gradnorm
            tot_pot_previous = tot_pot

            if (conj_grad) then
                forces_previous(:,:) = forces(:,:)
            end if

            ! get energies and forces at the new coordinates
            CALL get_energy_gradient(positions,tot_pot,forces,gradnorm, suppress_flag = .true.)


            ! Print iteration data
            if (.not. debug_flag) then
                write(*,"(i5,10x,4(F18.8,3x))") iter, tot_pot, tot_pot-tot_pot_previous, gradnorm,&
                                                gradnorm-gradnorm_previous
            else
                write(*,"(i5,2x,4(F12.6,1x),2(f7.1),A,2(F10.2,A),4x,F6.3,2x, F10.6)") &
                    iter, tot_pot, tot_pot-tot_pot_previous, gradnorm, gradnorm-gradnorm_previous,&
                    tot_pot_P1, tot_pot_P3, "  [ f(x) = (", a,") x^2 + (",b,") x ]", &
                    best_step, maxval(best_step*forces(:,:)/gradnorm)
            end if

            ! Write trajectory file for the minimization
            open(97, file=traj_xyzfile, status='old', action='write')
            write(97,*) n_atoms
            write(97,*) "atomic positions at iteration = ",iter
            do i=1, size(positions,1), 1
                write(97,FMT='(A3,3(2X,F15.8))') atomnames(i), positions(i,1:3)
            end do

            ! Check for convergence
            !if (ABS(gradnorm-gradnorm_previous)<conv_gradnorm .and. .not. converged_grad) then
            if (ABS(gradnorm-gradnorm_previous)<min_ftol .and. .not. converged_grad) then
                write(*,*) "gradient converged"
                converged_grad = .true.
            end if
            if (ABS(tot_pot-tot_pot_previous)<min_etol .and. .not. converged_pot) then
                write(*,*) "potential converged"
                converged_pot = .true.
            end if
            if (converged_grad .and. converged_pot) then
                converged = .true.
            end if
            if (iter == min_max_iter) then
                exit
            end if

            if (opt_method == 2) then
                conj_grad = .true.
            end if
        end do

        deallocate(positions_P1)
        deallocate(positions_P2)
        deallocate(positions_P3)

        ! Print convergence message
        write(*,"(/A)") "==================================================================="
        if (converged) then
            if (opt_method == 2) then
                write(*,"(A,i4,A/)")"Conjugate gradient minimization converged in", iter," iterations"
            else
                write(*,"(A,i4,A/)")"Gradient Descent converged in", iter," iterations"
            end if
        else
            write(*,"(A/)") "WARNING: energy minimization did not converge in the max number of iterations"
        end if

        if (debug_flag) then
            call recprt2("forces on the atoms after minimization",atomnames(:),forces(:,:),n_atoms)
            write(*,"(A,F16.8)") "Resulting potential energy:", tot_pot
            call recprt2("Atomic coordinates BEFORE minimization",atomnames(:),input_positions(:,:),n_atoms)
            call recprt2("Atomic coordinates AFTER minimization",atomnames(:),positions(:,:),n_atoms)
            call recprt2("Change in atomic coordinates THROUGH minimization",atomnames(:),&
                                                            positions(:,:)-input_positions(:,:),n_atoms)
        end if
        call displacement_vec(input_positions,positions,n_atoms,atomnames, displacement)
        write(*,*) "Displacements"
        do i = 1, n_atoms
            write(*,"(I3,1x,A3,1x,F16.12,1x,A)") i,atomnames(i),displacement(i),"Å"
        end do
        write(*,"(/A,F12.8,A)") "  sum = ", sum(displacement(:)), " Å"

        write(*,*) "==================================================================="

        ! Writing the updated coordinates to an xyzfile
        dot = index(xyzfile, ".", back=.true.)     ! find last "." in xyzfile name
        minimized_xyzfile = xyzfile(:dot-1) // ".minimized" // xyzfile(dot:)
        open(99, file=minimized_xyzfile, status='replace', action='write')
        write(99,*) n_atoms
        write(99,"(A)") "New atomic positions:"
        do i=1, size(positions,1), 1
            write(99,FMT='(A3,5X,F15.6,2X,F15.6,2X,F15.6)') atomnames(i), positions(i,1), positions(i,2),positions(i,3)
        end do


        write(*,*) "Wrote updated coordinates to ", minimized_xyzfile
        write(*,*) "Wrote minimization trajectory to ", traj_xyzfile


    end subroutine minimization



end module minimization_mod
