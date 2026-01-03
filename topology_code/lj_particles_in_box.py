# script to generate a .top and a .xyz file as input for the MDProgram
# AIM: with knowledge of the (cubic) box size, LJ particles of different species can be created to get an initial system in a box

import sys
import random
main_name = "LJ_particles"

helpmode = False

for flag in sys.argv:
    if flag == "-h":
        helpmode = True
        print("usage: python3 lj_particles_in_box.py [-h] 'properties_file' [-] [-e_temp] [...] [-savefig]")
        print("\nyou may use different flags to track different properties:")
        print("     -e_tot          plot the total energy")
        print("     -e_pot          plot the potential energy")


#input paramters
n_atoms = 100
atomtype_short = "Ar"
atomtype_long = "Ar"
weight = 40
boxlength = 20
sigma = 2
eps = 0.1

#n_species

#usage: get_topology.py [-h] (-s SMI | -p PDB) [-o filename] [-m | -r | -a ASSIGN] [-c CONSTRAINTS] [-d]


topfile = str(main_name)+".top"

with open(topfile, "w") as f:
    f.write("[AtomTypes] %i"%n_atoms)
    for atom in range(n_atoms):
        f.write("\n%-4i%-4s%-4s%f"%(atom+1,atomtype_short,atomtype_long,weight))
    f.write("\n\n[Bonds] 0")
    f.write("\n\n[Angles] 0")
    f.write("\n\n[Dihedrals] 0")
    f.write("\n\n[ImproperDihedrals] 0")
    f.write("\n\n[LJ]")
    for atom in range(n_atoms):
        f.write("\n%-4i%10.5f%10.5f"%(atom+1,sigma,eps))
    f.write("\n\n[Charges]")


xyzfile = str(main_name) + ".xyz"

with open(xyzfile, "w") as d:
    d.write(str(n_atoms))
    d.write("\nxyz file created with lj_particles_in_box.py")
    for atom in range(n_atoms):
        x = boxlength * random.random()
        y = boxlength * random.random()
        z = boxlength * random.random()
        #d.write("\n%-4s"%atomtype_short+3*"%12.6f"%(x,y,z))
        d.write("\n%-4s"%weight+3*"%12.6f"%(x,y,z))
    d.write("\n")
