
program MD
use definitions, only: wp
use f90getopt
use header_mod, only: pine_tree, final_phrase
use parser_mod, only: parser
use force_field_mod
use minimization_mod, only: minimization
use propagation, only: Verlet_propagator
use init_v_mod, only: init_v

implicit none
character(len=256) :: xyzfile, topofile
character(len=2), allocatable :: atomtypes(:),atomnames(:)
real(kind=wp), allocatable :: mweights(:),positions_previous(:,:),positions(:,:), forces(:,:)
real(kind=wp) :: start_time, end_time, tot_pot, gradnorm
character(len=1) :: short
logical :: t_present = .false. , c_present = .false., m_present = .false., m1_present =.false., p_present = .false.

! Parse the command line arguments with f90getopt library

! Declaration of options
type(option_s) :: opts(6)
opts(1) = option_s("top",.true.,  "t") ! has argument
opts(2) = option_s("coord",.true.,  "c") ! has argument
opts(3) = option_s("debug",.false., "d") ! has no argument
opts(4) = option_s("help",.false., "h") ! has no argument
opts(5) = option_s("minimize",.false., "m") ! has no argument
opts(6) = option_s("propagate",.false., "p") ! has no argument

! Handle absence of options
if (command_argument_count() .eq. 0) then
    write(*,*) "ERROR: Input data are missing. Use options -h or --help for details"
    stop
end if

do
    short = getopt("t:c:dhmnp", opts)
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
        case("n") !minimize with gradient descent if -m1 flag
                m1_present = .true.
        case("m") !minimize with conjugate gradient if -m flag
                m_present = .true.
        case("p") !propagate if -p flag present
                p_present = .true.
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

CALL CPU_TIME(start_time)
CALL pine_tree()

CALL parser(xyzfile,topofile,n_atoms,n_bonds,n_angles,n_impdie,n_torsions,mweights,positions,atomtypes,bond_params,&
                angle_params,impdihedrals_params,tors_params,lj_params,resp_charges,debug_flag,atomnames)



allocate(forces(n_atoms,3))
CALL force_field_calc(n_atoms,n_bonds,n_angles,n_impdie,n_torsions,positions,bond_params,angle_params,&
            impdihedrals_params,tors_params,lj_params,resp_charges,tot_pot,forces,debug_flag, suppress_flag = .false.)

! do minimization if -m flag active
if (m_present) then
    CALL minimization(positions,n_atoms,tot_pot,forces, debug_flag,xyzfile,atomnames,2) ! conjugate gradient
    !
else if (m1_present) then
    CALL minimization(positions,n_atoms,tot_pot,forces, debug_flag,xyzfile,atomnames,1) !steepest descent
end if

if (p_present) then
    allocate(positions_previous(n_atoms,3))
    !call init_v()
    call Verlet_propagator(positions,positions_previous,mweights,n_atoms,debug_flag,atomnames,xyzfile)!,timestep,nsteps)
end if

deallocate(forces)

CALL CPU_TIME(end_time)
write(*,*) "Total CPU time: ", end_time - start_time, " seconds"

CALL final_phrase()

end program
