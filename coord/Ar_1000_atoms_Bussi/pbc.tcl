package require pbctools

# Box parameters
set box {40 40 40 90 90 90}

# Procedure to update the PBC box
proc update_pbc_box {args} {
    global box
    pbc set $box
    pbc box
}

# Remove any existing callbacks to avoid duplicates
trace remove variable ::vmd_frame write update_pbc_box

# Register callback: called every time the frame changes
trace add variable ::vmd_frame write update_pbc_box

# Set atom representation to VDW spheres
mol delrep 0 top
mol representation VDW 0.5 12.0
mol color Name
mol selection all
mol material Opaque
mol addrep top

# Optional: make the view nicer
display projection Orthographic
axes location Off

