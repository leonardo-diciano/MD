
program MD 

implicit none
character(len=256) :: xyzfile, topofile 
character, allocatable :: atomtypes(:)
integer :: n_atoms,n_bonds,n_angles,n_torsions,n_impdie
real, allocatable :: mweights(:),positions(:,:),bond_params(:,:),angle_params(:,:),impdihedrals_params(:,:),&
                                    tors_params(:,:),lj_params(:,:), resp_charges(:,:)

xyzfile = "../topology_code/benzene.xyz"
topofile = "../topology_code/benzene.top"

CALL parser(xyzfile,topofile,n_atoms,n_bonds,n_angles,n_impdie,n_torsions,mweights,positions,atomtypes,bond_params,angle_params,&
            impdihedrals_params,tors_params,lj_params,resp_charges)




contains
! Subroutine to parse topology files and store data 
! ***NOT TESTED version***
! Author: Leonardo Di Ciano (2025)
subroutine parser(xyzfile,topofile,n_atoms,n_bonds,n_angles,n_impdie,n_torsions,mweights,positions,atomtypes,bond_params,&
                angle_params,impdihedrals_params,tors_params,lj_params,resp_charges)

implicit none

character(len=*), intent(in) :: xyzfile, topofile 
integer :: n_atoms,n_bonds,n_angles,n_torsions,n_impdie
character, allocatable, intent(out) :: atomtypes(:)
real, allocatable, intent(out) :: mweights(:),positions(:,:),bond_params(:,:),angle_params(:,:),impdihedrals_params(:,:),&
                                    tors_params(:,:),lj_params(:,:), resp_charges(:,:)
character(len=256) :: line
logical :: inAtomBlock,inBondBlock,inAngleBlock,inImpDieBlock, inDieBlock, inLJBlock, inChgBlock 
integer :: count, dummy_idx, io
character(len=20) :: dummy_symb

inAtomBlock = .false.
inBondBlock = .false.
inAngleBlock = .false.
inImpDieBlock = .false.
inDieBlock = .false. 
inLJBlock = .false. 
inChgBlock = .false.

open(unit=10,file=topofile,status='old',access='sequential',action='read')


do
    read(10,'(A)',iostat=io) line
    if (io /= 0) exit  

    if (index(line,'[AtomTypes]') == 1) then !Starting Block [Atom]
        inAtomBlock = .true.
        read(line, *) dummy_symb, n_atoms   
        allocate(atomtypes(n_atoms))    !Allocating the arrays for storing data
        allocate(mweights(n_atoms))
        count=0   !the count is resetted in each section, for correct array indexing
        cycle
    elseif (index(line,'[Bonds]') == 1) then !Starting Block [Bonds]
        inBondBlock = .true.
        read(line, *) dummy_symb, n_bonds
        allocate(bond_params(n_bonds,4))
        count=0
        cycle
    elseif (index(line,'[Angles]') == 1) then !Starting Block [Angles]
        inAngleBlock = .true.
        read(line, *) dummy_symb, n_angles
        allocate(angle_params(n_angles,5))
        count=0
        cycle
    elseif (trim(line) == '[ImproperDihedrals]') then !Starting Block [ImproperDihedrals]
        inImpDieBlock = .true.
        read(line, *) dummy_symb, n_impdie
        if (n_impdie == 0 ) then
            allocate(impdihedrals_params(1,1))
            impdihedrals_params(1,1)=0
        else
            allocate(impdihedrals_params(n_impdie,8))
        end if
        count=0
        cycle
    elseif (trim(line) == '[Dihedrals]') then !Starting Block [Dihedrals]
        inDieBlock = .true.
        read(line, *) dummy_symb, n_torsions
        if (n_torsions == 0 ) then
            allocate(tors_params(1,1))
            tors_params(1,1)=0
        else
            allocate(tors_params(n_torsions,8))
        end if
        count=0
        cycle
    elseif (trim(line) == '[LJ]') then !Starting Block [LJ]
        inLJBlock = .true.
        count=0
        allocate(lj_params(n_atoms,3))
        cycle
    elseif (trim(line) == '[Charges]') then !Starting Block [Charges]
        inChgBlock = .true.
        count=0
        allocate(resp_charges(n_atoms,2))
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
                                impdihedrals_params(count,7),impdihedrals_params(count,8)            
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

! now parsing the XYZ file
open(unit=11,file=xyzfile,status='old',access='sequential',action='read')
allocate(positions(n_atoms,4))
count=0
do
    read(11,'(A)') line
    if ( count < 2 ) then !Skip first two lines of XYZ file
        count = count + 1
        cycle
    else
        !use count-1 to match the atom index
        read(line, *) positions(count-1,1), positions(count - 1,2), positions(count-1,3), positions(count-1,4)
        count = count + 1
    end if
end do
close(11)

end subroutine

end program
