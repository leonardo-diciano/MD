subroutine parser(xyzfile,topofile,n_atoms,n_bonds,n_angles,n_torsions,positions,atomtypes,bond_const,bond_dist,angle_const,tors_const,&
    n_value,tors_phase,lj_eps,lj_sigma,charges)

implicit none

character, intent(in) :: xyzfile, topofile 
integer :: n_atoms,n_bonds,n_angles,n_torsions
character, allocatable, intent(out) :: atomtypes(:)
real, allocatable, intent(out) :: positions(:,:),bond_const(:),bond_dist(:),angle_const(:),tors_const(:),tors_phase(:),&
                                  n_value(:), lj_eps(:),lj_sigma(:), charges(:)
end subroutine
