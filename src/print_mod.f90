module print_mod

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Module containing three subroutines for rectangular matrix printing:
    !
    ! all of the recprt subroutines specify title string (printed as header for the matrix); the matrix to be printed (mat); and its dimensions
    !
    ! recprt1: general, requires input of number printing format, e.g. FmtIn="(*(F10.6,1X))"
    !
    ! recprt2: for printing atom coordinate dependent properties -> matrices of dimensionality (n_atoms,3);
    !          format of printed numbers is fixed in the subroutine; this subroutines requires specifying atomnames(:);
    !
    ! recprt3: like recprt2, but without atomnames(:) -> instead, only line /atom numbers are printed
    !
    ! Author: Lila Zapp (2025)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


contains
    subroutine recprt(title,fmtin,mat,nrow,ncol)
        ! for printing a matrix along with its name
        use definitions, only: wp

        implicit none
        integer, intent(in) :: nrow,ncol
        character(len=*), intent(in) :: Title,FmtIn
        real(kind=wp), intent(in) :: mat(:,:)
        integer :: i

        !FmtIn is a format statement for the number print; e.g. "(*(F10.6,1X))"
        write(*,"(/AA)") Title, ":"
        do i=1,nrow
            write(*,FmtIn) mat(i,1:ncol)
        end do
        write(*,*) ""

    end subroutine recprt
!=====================================================================================================================
    subroutine recprt2(title,atomnames,mat,n_atoms)
        ! for printing a matrix along with its name; this routine is for (n_atoms,3) matrices specifically;
        !   it prints also the atomnames for each row
        use definitions, only: wp
        implicit none
        integer, intent(in) :: n_atoms
        character(len=*), intent(in) :: Title
        character(len=2), intent(in) :: atomnames(:)
        real(kind=wp), intent(in) :: mat(n_atoms,3)
        integer :: i

        write(*,"(/AA)") Title, ":"
        do i=1,n_atoms
            write(*,"(I3,2X,A2,2X,3(F18.8,1X))") i,atomnames(i),mat(i,1:3)
        end do
        write(*,*) ""

    end subroutine recprt2

!=====================================================================================================================

    subroutine recprt3(title,mat,n_atoms)
        ! for printing a matrix along with its number (doesnt require to define atomnames);
        !this routine is for (n_atoms,3) matrices specifically;
        use definitions, only: wp
        implicit none
        integer, intent(in) :: n_atoms
        character(len=*), intent(in) :: Title
        real(kind=wp), intent(in) :: mat(n_atoms,3)
        integer :: i

        write(*,"(/AA)") Title, ":"
        do i=1,n_atoms
            write(*,"(I3,3X,3(F16.8,1X))") i,mat(i,1:3)
        end do
        write(*,*) ""

    end subroutine recprt3

!=====================================================================================================================




end module print_mod
