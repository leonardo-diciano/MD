! Module with subroutines to parse
! topology files, input files and store data
!
! Author: Leonardo Di Ciano (2025)

module parser_mod
use definitions, only: wp
implicit none

! General
public :: atomnames
character(len=2), allocatable :: atomnames(:)
logical :: debug_flag = .false.

! minimization params
public :: min_max_iter, min_etol, min_ftol, min_alpha, min_debug
integer :: min_max_iter=10000
real(kind=wp) :: min_etol=1.0e-6, min_ftol=1.0e-4, min_alpha = 1e-3
logical :: min_debug = .false.

! MD params
public :: md_ts,md_nsteps,md_ensemble,md_temp, md_press, md_boxlength, md_debug, md_int
integer :: md_nsteps=1000
real(kind=wp) :: md_ts=1.0 ,md_boxlength = 15, md_temp=300.0, md_press=100000.0! in Pa
character(len=32) :: md_ensemble="NVE", md_int="verlet"
logical :: md_debug = .false., debug_print_all_matrices = .false., md_fix_com_mom = .false., md_pbc = .false.

! Bussi thermostat params
public :: bus_tau
real(kind=wp) :: bus_tau = 50.0 ! fs

! Berendsen barostat params
public :: ber_tau, ber_k
real(kind=wp) :: ber_tau = 5000.0 !fs, following GROMACS default  
real(kind=wp) :: ber_k = 4.6e-10 ! Pa^{-1} for water at ~ 300K and 1 atm

! Metadynamics
public :: meta_cv, meta_tau, meta_nsteps, meta_cv_type, meta_dT, meta_omega, meta_sigma
integer(kind=wp), allocatable :: meta_cv(:)
real(kind=wp) :: meta_tau = 10.0 ,meta_dT = 2700, meta_omega = 1, meta_sigma = 0.2
integer(kind=wp) ::  meta_nsteps = 100
character(len=32) :: meta_cv_type

contains

subroutine parser_top(xyzfile,topofile,positions,atomtypes,debug_flag)
use definitions, only: wp
use force_field_mod, only: n_atoms,n_bonds,n_angles,n_torsions,n_impdie, bond_params, angle_params,& 
                    impdihedrals_params, tors_params,lj_params,resp_charges,mweights
implicit none

character(len=*), intent(in) :: xyzfile, topofile
logical,intent(in) :: debug_flag
character(len=2), allocatable, intent(out) :: atomtypes(:)
real(kind=wp), allocatable, intent(out) :: positions(:,:)
character(len=256) :: line
logical :: inAtomBlock,inBondBlock,inAngleBlock,inImpDieBlock, inDieBlock, inLJBlock, inChgBlock
integer :: count, dummy_idx, io, i
character(len=256) :: dummy_symb

! Initialize the block flags
inAtomBlock = .false.
inBondBlock = .false.
inAngleBlock = .false.
inImpDieBlock = .false.
inDieBlock = .false.
inLJBlock = .false.
inChgBlock = .false.

! Open topology file
open(unit=10,file=topofile,status='old',access='sequential',action='read')

if (debug_flag) then
    write(*,*) ""
    write(*,*) "Parsing the topology file"
end if

do
    read(10,'(A)',iostat=io) line
    if (io /= 0) exit

    if (index(line,'[AtomTypes]') == 1) then !Starting Block [Atom]
        inAtomBlock = .true.
        read(line, *) dummy_symb, n_atoms
        allocate(atomtypes(n_atoms))    !Allocating the arrays for storing data
        allocate(mweights(n_atoms))
        count=0   !the count is resetted in each section, for correct array indexing

        if (debug_flag) then
            write(*,*) "Reading ", n_atoms ," Atom Types"
        end if

        cycle

    elseif (index(line,'[Bonds]') == 1) then !Starting Block [Bonds]
        inBondBlock = .true.
        read(line, *) dummy_symb, n_bonds
        allocate(bond_params(n_bonds,4))
        count=0

        if (debug_flag) then
            write(*,*) "Reading ", n_bonds ," bonds terms"
        end if

        cycle

    elseif (index(line,'[Angles]') == 1) then !Starting Block [Angles]
        inAngleBlock = .true.
        read(line, *) dummy_symb, n_angles
        allocate(angle_params(n_angles,5))
        count=0

        if (debug_flag) then
            write(*,*) "Reading ", n_angles ," angle terms"
        end if

        cycle

    elseif (index(line,'[ImproperDihedrals]') == 1) then !Starting Block [ImproperDihedrals]
        inImpDieBlock = .true.
        read(line, *) dummy_symb, n_impdie
        if (n_impdie == 0 ) then
            allocate(impdihedrals_params(0,7))
        else
            allocate(impdihedrals_params(n_impdie,7))
        end if
        count=0

        if (debug_flag) then
            write(*,*) "Reading ", n_impdie ," improper dihedral terms"
        end if

        cycle

    elseif (index(line,'[Dihedrals]') == 1) then !Starting Block [Dihedrals]
        inDieBlock = .true.
        read(line, *) dummy_symb, n_torsions
        if (n_torsions == 0 ) then
            allocate(tors_params(0,8))
        else
            allocate(tors_params(n_torsions,8))
        end if
        count=0

        if (debug_flag) then
            write(*,*) "Reading ", n_torsions ," dihedral terms"
        end if

        cycle

    elseif (index(line, '[LJ]') == 1) then !Starting Block [LJ]
        inLJBlock = .true.
        count=0
        allocate(lj_params(n_atoms,3))

        if (debug_flag) then
            write(*,*) "Reading ", n_atoms ," LJ parameters"
        end if

        cycle

    elseif (index(line, '[Charges]') == 1) then !Starting Block [Charges]
        inChgBlock = .true.
        count=0
        allocate(resp_charges(n_atoms,2))

        if (debug_flag) then
            write(*,*) "Reading ", n_atoms ," partial charges"
        end if

        cycle

    end if

    if (inChgBlock .and. trim(line) == "" ) then
        exit
    elseif (trim(line) == "" ) then !If blank line, new block incoming, reset all status
        inAtomBlock = .false.
        inBondBlock = .false.
        inAngleBlock = .false.
        inImpDieBlock = .false.
        inDieBlock = .false.
        inLJBlock = .false.
        inChgBlock = .false.
    end if

    if (inAtomBlock) then
        count = count + 1
        read(line, *) dummy_idx, dummy_symb, atomtypes(count), mweights(count)
    elseif (inBondBlock) then
        count = count + 1
        read(line, *) bond_params(count,1), bond_params(count,2), bond_params(count,3), bond_params(count,4)
    elseif (inAngleBlock) then
        count = count + 1
        read(line, *) dummy_symb, angle_params(count,1), angle_params(count,2), angle_params(count,3), angle_params(count,4), &
                                 angle_params(count,5)
    elseif (inImpDieBlock) then
        count = count + 1
        read(line, *) dummy_symb, impdihedrals_params(count,1),impdihedrals_params(count,2),impdihedrals_params(count,3), &
                                impdihedrals_params(count,4),impdihedrals_params(count,5),impdihedrals_params(count,6), &
                                impdihedrals_params(count,7)
    elseif (inDieBlock) then
        count = count + 1
        read(line, *) dummy_symb, tors_params(count,1),tors_params(count,2),tors_params(count,3),tors_params(count,4), &
                    tors_params(count,5),tors_params(count,6),tors_params(count,7),tors_params(count,8)
    elseif (inLJBlock) then
        count = count + 1
        read(line, *) lj_params(count,1), lj_params(count,2), lj_params(count,3)
    elseif (inChgBlock) then
        count = count + 1
        read(line, *) resp_charges(count,1), resp_charges(count,2)
    end if

end do
close(10)

if (debug_flag) then
    write(*,*) ""
    write(*,*) "Parsing the topology file - Done"
    write(*,*) "Summary of stored data sizes: "
    write(*,FMT='( "Bond parameters: ", I4," rows ",I4," cols")') size(bond_params,1),size(bond_params,2)
    write(*,FMT='( "Angle parameters: ", I4," rows ",I4," cols")') size(angle_params,1),size(angle_params,2)
    write(*,FMT='( "Dihedrals parameters: ", I4," rows ",I4," cols")') size(tors_params,1),size(tors_params,2)
    write(*,FMT='( "Improper dihedrals parameters: ", I4," rows ",I4," cols")') size(impdihedrals_params,1),&
                    size(impdihedrals_params,2)
    write(*,FMT='( "LJ parameters: ", I4," rows ",I4," cols")') size(lj_params,1),size(lj_params,2)
    write(*,FMT='( "Partial charges: ", I4," rows ",I4," cols")') size(resp_charges),size(resp_charges,2)
end if

! now parsing the XYZ file
open(unit=11,file=xyzfile,status='old',access='sequential',action='read')
allocate(positions(n_atoms,3))

if (debug_flag) then
    write(*,*) ""
    write(*,*) "Reading XYZ coordinates"
end if


allocate(atomnames(n_atoms))    !Allocating the arrays for storing data
count = 0
do
    read(11,'(A)',iostat=io) line
    if (io /= 0) then
        exit
    elseif ( count < 2 ) then !Skip first two lines of XYZ file
        count = count + 1
        cycle
    elseif (count > n_atoms + 1) then
        exit
    else
        !use count-2 to match the atom index when writing the array
        read(line, *) atomnames(count-1), positions(count-1,1), positions(count-1,2), positions(count-1,3)
        count = count + 1
    end if
end do
close(11)
if (debug_flag) then
    write(*,*) "Reading XYZ coordinates - Done"
end if

! Write output
write(*,*)""
write(*,*)"Parsing of input files completed"
write(*,*)""
write(*,*) "Atom Types and XYZ coordinates (Å)"
write(*,FMT='("Index",2X,"Atom Type",8X,"X",15X,"Y",15X,"Z")')
do i=1, size(positions,1), 1
    write(*,FMT='(I3,5X,A2,2X,F15.6,2X,F15.6,2X,F15.6)') i, atomtypes(i),positions(i,1), positions(i,2),positions(i,3)
end do
end subroutine



subroutine parser_input(inputfile,xyzfile, topofile, t_present, c_present, m_present, m1_present,&
         p_present, meta_present)
! Subroutine to parse input files 

character(len=256), intent(in) :: inputfile
character(len=256), intent(out) :: xyzfile, topofile
logical, intent(inout) :: t_present, c_present, m_present, m1_present,p_present, meta_present
integer :: io, dummy_idx
logical :: mini_block = .false., md_block = .false., meta_block = .false.
character(len=256) :: line
character(len=32) :: dummy_symb


open(unit=12,file=inputfile,status='old',access='sequential',action='read')

do
    read(12,'(A)',iostat=io) line

    if (io /= 0) exit

    if (index(line,'#') == 1) cycle ! Use # for input file comments

    if (index(trim(line), "topology") == 1) then
        read(line, *) dummy_symb, topofile
        if (len_trim(topofile) > 0) then
            t_present = .true.
        end if
    end if

    if (index(trim(line), "coords") == 1) then
        read(line, *) dummy_symb, xyzfile
        if (len_trim(xyzfile) > 0) then
            c_present = .true.
        end if
    end if

    if (index(trim(line), "debug") == 1) then
        read(line, *) dummy_symb, dummy_idx
        if (dummy_idx .le. 0 ) then
            debug_flag = .false.
            md_debug = .false.
            min_debug = .false.
            debug_print_all_matrices = .false.
        elseif (dummy_idx == 1) then
            debug_flag = .true.
            md_debug = .true.
            min_debug = .true.
        elseif (dummy_idx > 1) then
            debug_flag = .true.
            md_debug = .true.
            min_debug = .true.
            debug_print_all_matrices = .true.
        end if
    end if

    if (index(trim(line),"[minimize]") == 1) then
        mini_block = .true.
        cycle
    end if

    if (index(trim(line),"[run_md]") == 1) then
        mini_block = .false.
        md_block = .true.
        p_present = .true.
        cycle
    end if

    if (index(trim(line),"[run_meta]") == 1) then
        mini_block = .false.
        md_block = .false.
        meta_block = .true.
        p_present = .false.
        meta_present = .true.
        cycle
    end if


    if (mini_block) then
        if (index(trim(line),"conj_grad") == 1) then
            m_present = .true.
        elseif (index(trim(line),"steep_desc") == 1) then
            m1_present = .true.
        elseif (index(trim(line), "etol") == 1) then
            read(line,*) dummy_symb, min_etol
        elseif (index(trim(line),"ftol") == 1) then
            read(line,*) dummy_symb, min_ftol
        elseif (index(trim(line),"maxiter") == 1) then
            read(line,*) dummy_symb, min_max_iter
        endif
    end if

    if (md_block) then
        if (index(trim(line),"ts") == 1) then
            read(line,*) dummy_symb, md_ts
        elseif (index(trim(line),"nsteps") == 1) then
            read(line,*) dummy_symb, md_nsteps
        elseif (index(trim(line),"integrator") == 1) then
            read(line,*) dummy_symb, md_int
        elseif (index(trim(line),"temp") == 1) then
            read(line,*) dummy_symb, md_temp
        elseif (index(trim(line), "ensemble") == 1) then
            read(line,*) dummy_symb, md_ensemble
        elseif (index(trim(line),"press") == 1) then
            read(line,*) dummy_symb, md_press
        elseif (index(trim(line),"fix_com") == 1) then
            read(line,*) dummy_symb, dummy_symb
            if (dummy_symb == "true") then
                md_fix_com_mom = .true.
            endif
        elseif (index(trim(line),"debug") == 1) then
            md_debug = .true.
        elseif (index(trim(line),"berendsen_tau") == 1) then
            read(line,*) dummy_symb, ber_tau
        elseif (index(trim(line),"berendsen_k") == 1) then
            read(line,*) dummy_symb, ber_k
        elseif (index(trim(line),"bussi_tau") == 1) then
            read(line,*) dummy_symb, bus_tau
        end if
    end if

    if (meta_block) then
        if (index(trim(line),"cv") == 1) then
            read(line,*) dummy_symb, meta_cv_type
            if (meta_cv_type == "distance") then
                allocate(meta_cv(2))
                read(line,*) dummy_symb, meta_cv_type, meta_cv(1), meta_cv(2)
            elseif (meta_cv_type == "angle") then
                allocate(meta_cv(3))
                read(line,*) dummy_symb, meta_cv_type, meta_cv(1), meta_cv(2), meta_cv(3)
            elseif (meta_cv_type == "dihedral") then
                allocate(meta_cv(4))
                read(line,*) dummy_symb, meta_cv_type, meta_cv(1), meta_cv(2), meta_cv(3), meta_cv(4)
            else
                write(*,*) "Unrecognized CV type for metadynamics: ", meta_cv_type
                write(*,*) "The allowed ones are: distance, angle, dihedral"
                stop
            end if 
        elseif (index(trim(line),"nsteps") == 1) then
            read(line,*) dummy_symb, meta_nsteps
        elseif (index(trim(line),"tau") == 1) then
            read(line,*) dummy_symb, meta_tau
        elseif (index(trim(line),"dT") == 1) then
            read(line,*) dummy_symb, meta_dT
        elseif (index(trim(line),"omega") == 1) then
            read(line,*) dummy_symb, meta_omega
        elseif (index(trim(line),"sigma") == 1) then
            read(line,*) dummy_symb, meta_sigma
        end if
    end if

end do

close(12)
end subroutine

end module
