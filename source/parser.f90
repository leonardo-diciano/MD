subroutine parser(xyzfile,topofile,n_atoms,n_bonds,n_angles,n_torsions,positions,atomtypes,bond_const,bond_dist,angle_const,tors_const,&
    n_value,tors_phase,lj_eps,lj_sigma,charges)

implicit none

character, intent(in) :: xyzfile 
character, intent(in) :: topofile 
integer :: n_atoms
integer :: n_bonds
integer :: n_angles
integer :: n_torsions
real, allocatable, intent(out) :: positions(:,:)
character, allocatable, intent(out) :: atomtypes(:)
real, allocatable, intent(out) :: bond_const(:)
real, allocatable, intent (out) :: bond_dist(:)
real, allocatable, intent(out) :: angle_const(:)
real, allocatable, intent(out) :: tors_const(:)
real, allocatable, intent(out) :: tors_phase(:)
real, allocatable, intent (out) :: n_value(:)
real, allocatable, intent(out) :: lj_eps(:)
real, allocatable, intent(out) :: lj_sigma(:)
real, allocatable, intent(out) :: charges(:)



    
end subroutine