! Subroutine to parse topology files
!
! Author: Leonardo Di Ciano (2025)

subroutine parser(xyzfile,topofile,n_atoms,n_bonds,n_angles,n_impdie,n_torsions,mweights,positions,atomtypes,bond_params,angle_params,impdihedrals_params, &
                tors_params,lj_params,resp_charges)

implicit none

character, intent(in) :: xyzfile, topofile 
integer :: n_atoms,n_bonds,n_angles,n_torsions,n_impdie
character, allocatable, intent(out) :: atomtypes(:)
real, allocatable, intent(out) :: mweights(:),positions(:,:),bond_params(:,:),angle_params(:,:),impdihedrals_params(:,:),tors_params(:,:),&
                                  lj_params(:,:), resp_charges(:,:)
character(len=256) :: line
logical :: inAtomBlock,inBondBlock,inAngleBlock,inImpDieBlock, inDieBlock, inLJBlock, inChgBlock 
integer :: count
integer :: dummy_idx
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
    read(10,'(A)') line

    if (trim(line) == '[AtomTypes]') then !Starting Block [Atom]
        inAtomBlock = .true.
        read(line, *) dummy_symb, n_atoms
        allocate(atomtypes(n_atoms))
        count=0
        cycle
    elseif (trim(line) == '[Bonds]') then !Starting Block [Bonds]
        inBondBlock = .true.
        read(line, *) dummy_symb, n_bonds
        allocate(bond_params(n_bonds,4))
        count=0
        cycle
    elseif (trim(line) == '[Angles]') then !Starting Block [Angles]
        inAngleBlock = .true.
        read(line, *) dummy_symb, n_angles
        allocate(angle_params(n_angles,5))
        count=0
        cycle
    elseif (trim(line) == '[ImproperDihedrals]') then !Starting Block [ImproperDihedrals]
        inImpDieBlock = .true.
        read(line, *) dummy_symb, n_impdie
        allocate(impdihedrals_params(n_impdie,7))
        count=0
        cycle
    elseif (trim(line) == '[Dihedrals]') then !Starting Block [Dihedrals]
        inDieBlock = .true.
        read(line, *) dummy_symb, n_torsions
        allocate(tors_params(n_torsions,7))
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
    endif

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
    endif

    if (inAtomBlock) then
        count = count + 1 
        read(line, *) dummy_idx, dummy_symb, atomtypes(count), mweights(count)
    elseif (inBondBlock) then
        count = count + 1 
        read(line, *) 
    end if
    
    

enddo
close(10)


end subroutine
