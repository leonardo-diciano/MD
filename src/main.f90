
program MD 
use f90getopt
use header_mod
use parser_mod
use force_field_mod

implicit none
character(len=256) :: xyzfile, topofile
integer :: n_atoms,n_bonds,n_angles,n_torsions,n_impdie
character(len=2), allocatable :: atomtypes(:)
real, allocatable :: mweights(:),positions(:,:),bond_params(:,:),angle_params(:,:),impdihedrals_params(:,:),&
                                    tors_params(:,:),lj_params(:,:), resp_charges(:,:), forces(:,:)
real :: tot_pot
character(len=1) :: short
logical :: t_present = .false. , c_present = .false., debug_flag = .false.



! Parse the command line arguments with f90getopt library

! Declaration of options
type(option_s) :: opts(5)
opts(1) = option_s("top",.true.,  "t") ! has argument
opts(2) = option_s("coord",.true.,  "c") ! has argument
opts(3) = option_s("debug",.false., "d") ! has no argument
opts(4) = option_s("help",.false., "h") ! has no argument
opts(5) = option_s("minimize",.false., "m") ! has no argument

! Handle absence of options
if (command_argument_count() .eq. 0) then
    write(*,*) "ERROR: Input data are missing. Use options -h or --help for details"
    stop
end if

CALL pine_tree()

do
    short = getopt("t:c:dhm", opts) 
    select case(short)
        case(char(0))
            exit
        case("t") ! option -t --topology
            t_present=.true.
            if (trim(optarg) > "") then 
                topofile = trim(optarg)
            else
                stop
            end if
        case("c") ! option -c --coord
            c_present=.true.
            if (trim(optarg) > "") then 
                xyzfile = trim(optarg)
            else
                stop
            end if
        case("d")   ! option -d --debug 
            debug_flag = .true.
        case("h") ! help output
            write(*, '(6(A/),/,4(A/))')&
                "Usage: md_program [options] …",&
                "Options:",&
                "  -t file.top  --top=file.top      Topology file - Required",&
                "  -c coord.xyz --coord=coord.xyz   XYZ coordinate file - Required",&
                "  -m           --minimize          Require minimization",&
                "  -h           --help              Print this help screen",&
                "  -d           --debug             Print extended output for debug",&
                "Examples:",&
                "  main -t file.top -c coord.xyz",&
                "  main -top=file.top -coord=coord.xyz"
            stop
    end select
end do

if (t_present .and. c_present) then ! both mandatory options are present
    write(*,FMT = '( "The topology file is: ", A50 )') topofile
    write(*,FMT = '( "The coordinate file is: ", A50 )') xyzfile
else if (t_present .and. (.not. c_present)) then
    write(*,*) "ERROR: Options -c (--coord) with argument is missing"
    stop
else if (c_present .and. (.not. t_present)) then
    write(*,*) "ERROR: Options -t (--top) with argument is mising" 
    stop
end if


CALL parser(xyzfile,topofile,n_atoms,n_bonds,n_angles,n_impdie,n_torsions,mweights,positions,atomtypes,bond_params,&
                angle_params,impdihedrals_params,tors_params,lj_params,resp_charges,debug_flag)


CALL force_field_calc(n_atoms,n_bonds,n_angles,n_impdie,n_torsions,positions,bond_params,angle_params,&
            impdihedrals_params,tors_params,lj_params,resp_charges,tot_pot,forces,debug_flag)


CALL final_phrase()

end program
