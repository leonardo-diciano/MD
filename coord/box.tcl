# Load pbctools (needed for drawing the box)
package require pbctools

# Set periodic box dimensions: a b c alpha beta gamma
pbc set {40 40 40 90 90 90}

# Draw the periodic boundary condition box
#pbc box -center origin -color black -width 2
pbc box

# Set atom representation to VDW spheres
mol delrep 0 top
#mol representation VDW 1.0 12.0
mol representation CPK
mol color Name
mol selection all
mol material Opaque
mol addrep top

# Optional: make the view nicer
display projection Orthographic
axes location Off

