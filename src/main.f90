
program MD
use definitions, only: wp
use f90getopt
use header_mod, only: pine_tree, final_phrase
use parser_mod, only: parser_top, parser_input, atomnames
use force_field_mod, only: n_atoms,n_bonds,n_angles,n_torsions,n_impdie, bond_params, angle_params,& 
                    impdihedrals_params, tors_params,lj_params,resp_charges,mweights, force_field_calc
use minimization_mod, only: minimization
use pbc_mod, only: define_box
use simulation_mod, only: simulation
use metadynamics_mod, only: run_metadynamics

implicit none
character(len=256) :: xyzfile, topofile, inputfile
character(len=2), allocatable :: atomtypes(:)
real(kind=wp), allocatable :: positions(:,:), forces(:,:)
real(kind=wp) :: start_time, end_time, tot_pot, gradnorm
character(len=1) :: short
logical :: t_present = .false. , c_present = .false., m_present = .false., m1_present =.false.,&
         p_present = .false., i_present = .false., suppress_flag = .true., debug_flag = .false., meta_present = .false.

! Parse the command line arguments with f90getopt library

! Declaration of options
type(option_s) :: opts(6)
opts(1) = option_s("top",.true.,  "t") ! has argument
opts(2) = option_s("coord",.true.,  "c") ! has argument
opts(3) = option_s("debug",.false., "d") ! has no argument
opts(4) = option_s("help",.false., "h") ! has no argument
opts(5) = option_s("minimize",.true., "m") ! has argument (SD or CG method)
opts(6) = option_s("propagate",.false., "p") ! has no argument

! Handle absence of options
if (command_argument_count() .eq. 0) then
    write(*,*) "ERROR: Input data are missing. Use options -h or --help for details"
    stop
end if

do
    short = getopt("i:t:c:m:dhp", opts)
    select case(short)
        case(char(0))
            exit
        case ("i") ! if input file given, we skip all command line parsing
            i_present = .true.
            if (trim(optarg) > "") then
                inputfile = trim(optarg)
            else
                stop
            end if
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
                "  -m [sd,cg]   --minimize          Require minimization with steepest descent(sd) or conjugate gradient(cg)",&
                "  -p           --minimize          Propagate the system using the Verlet integrator and default settings",&
                "  -h           --help              Print this help screen",&
                "  -d           --debug             Print extended output for debug",&
                "Examples:",&
                "  main -t file.top -c coord.xyz",&
                "  main -top=file.top -coord=coord.xyz"
            stop
        case("m") !minimize
            if ((trim(optarg) > "") .and. (trim(optarg)) == 'sd') then
                m1_present = .true. ! if explicitly requested for sd
            else
                m_present = .true. ! default choice is cg
            end if
        case("p") !propagate if -p flag present
                p_present = .true.
    end select
end do

if (i_present) then
    CALL parser_input(inputfile,xyzfile, topofile, t_present, c_present, m_present, m1_present,&
        p_present, meta_present)
else
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
end if

CALL CPU_TIME(start_time)
CALL pine_tree()

CALL parser_top(xyzfile,topofile,positions,atomtypes,debug_flag)


allocate(forces(n_atoms,3))
CALL force_field_calc(positions,tot_pot,forces,debug_flag, suppress_flag = .false.)

! do minimization if -m flag active
if (m_present) then
    CALL minimization(positions,tot_pot,forces,xyzfile,2) ! conjugate gradient
    !
else if (m1_present) then
    CALL minimization(positions,tot_pot,forces,xyzfile,1) !steepest descent
end if

if (meta_present) then
    call run_metadynamics(positions,xyzfile) ! run metadynamics
elseif (p_present) then
    call simulation(positions,xyzfile)      ! run MD
end if


!call define_box()


deallocate(forces)

CALL CPU_TIME(end_time)
write(*,"(/A,ES18.8,A)") "Total CPU time: ", end_time - start_time, " seconds"

CALL final_phrase()

end program
