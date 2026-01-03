# script to generate a .top and a .xyz file as input for the MDProgram
# AIM: with knowledge of the (cubic) box size, LJ particles of ONE species can be created to get an initial system in a box
# Author: Lila Zapp (2026)
import sys
import random

main_name = "LJ_particles"

helpmode = False

#input paramters
n_atoms = 1000
atomtype_short = "Ar"
atomtype_long = "Ar"
weight = 40
boxlength = 40
sigma = 2
eps = 0.1


topfile = str(main_name)+".top"

with open(topfile, "w") as f:
    f.write("[AtomTypes] %i"%n_atoms)
    for atom in range(n_atoms):
        f.write("\n%-6i%-4s%-4s%f"%(atom+1,atomtype_short,atomtype_long,weight))
    f.write("\n\n[Bonds] 0")
    f.write("\n\n[Angles] 0")
    f.write("\n\n[Dihedrals] 0")
    f.write("\n\n[ImproperDihedrals] 0")
    f.write("\n\n[LJ]")
    for atom in range(n_atoms):
        f.write("\n%-6i%10.5f%10.5f"%(atom+1,sigma,eps))
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
