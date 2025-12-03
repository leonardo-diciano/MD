
program MD 
use parser_mod
use force_field_mod

implicit none
character(len=256) :: xyzfile, topofile
integer :: n_atoms,n_bonds,n_angles,n_torsions,n_impdie
character, allocatable :: atomtypes(:)
real, allocatable :: mweights(:),positions(:,:),bond_params(:,:),angle_params(:,:),impdihedrals_params(:,:),&
                                    tors_params(:,:),lj_params(:,:), resp_charges(:,:)




CALL parser(xyzfile,topofile,n_atoms,n_bonds,n_angles,n_impdie,n_torsions,mweights,positions,atomtypes,bond_params,&
                angle_params,impdihedrals_params,tors_params,lj_params,resp_charges)



end program
