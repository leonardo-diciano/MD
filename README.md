# MDProgram

MDProgram is a molecular dynamics code developed by Leonardo Di Ciano and Lila Zapp as part of the Statistical Thermodynamics: Theory and Simulation course at Uppsala University (2025-2026). The program natively implements the General Amber Force Field (GAFF v1.4) and performs: single-point calculation, minimizations, molecular dynamics simulation with NVE, NVT and NPT ensembles and well-tempered metadynamics. The program is supported by a python-based module to generate topology files necessary to its functioning and by a second python module focused on the analysis of the simulated MD trajectories.
 
## Installation

To install MDProgram, clone the repository and execute the installer _install.sh_ in the MD directory.
```
>>  git clone https://github.com/lilazapp/MD.git
>>  cd MD
>>  ./install.sh
```
The compilation uses gfortran (>=13) and produces the md\_program binary.
The topology and trajectory analysis codes are available in their directories. The necessary python packages can be installed using conda and the dependencies file.
```
>>  conda create --file python_dependencies.yml
```

## Use of MDProgram 
The topology code can be used from command line.
```
>>  python topology_code/get_topology.py -h
 
usage: get_topology.py [-h] (-s SMI | -p PDB) [-o filename] [-m | -r | -a ASSIGN] [-c CONSTRAINTS] [-d]

options:
  -h, --help            show this help message and exit
  -s, --smi SMI         Input SMILES string
  -p, --pdb PDB         Input PDB file path
  -o, --output filename
                        Give the name for topology and XYZ coordinate file
  -m, --mulliken        Select Mulliken Charges calculation to obtain partial charges
  -r, --resp            Select RESP charges calculation to obtain partial charges
  -a, --assign ASSIGN   Give a list of atom types and corresponding charge to assign
  -c, --constraints CONSTRAINTS
                        Give a list of constrained groups of atoms for RESP charges calculation, i.e.
                        '[[1,2,][3,4,5,6,7,8]]' for ethane
  -d, --debug           Activate debug level of printing
```  
MDProgram can be run from command line inputs, with limited possibility of customization, or by using an input file.
```
>>  ./md_program -h
Usage: md_program [options] …
Options:
  -t file.top  --top file.top      Topology file - Required
  -c coord.xyz --coord coord.xyz   XYZ coordinate file - Required
  -i file.inp  --input file.inp    Input file with all information
  -m [sd,cg]   --minimize          Require minimization with steepest descent(sd) or conjugate gradient(cg)

  -p           --propagate         Propagate the system using the Verlet integrator and default settings
  -h           --help              Print this help screen
  -d           --debug             Print extended output for debug
Examples:

  main -i file.inp
  main -t file.top -c coord.xyz
  main -top=file.top -coord=coord.xyz

```
The structure of the input file with all possible parameters is provided in the MDProgram directory.
Finally, the trajectory analysis script can be run from command line only and produces the requested plots either as popup windows or as saved figures.
```
>>  python trajectory_analysis/plot.py -h
usage: python3 plot.py [-h] 'properties_file' [-energies] [-e_temp] [...] [-savefig]

you may use different flags to track different properties:
     -e_tot          plot the total energy
     -e_pot          plot the potential energy
     -e_kin          plot the kinetic energy
     -energies       plot all three energy components in one window
     -f_norm         plot the norm of the force
     -temp           plot the temperature
     -pressure       plot the pressure
     -com_mom        plot the norm of the momentum of the center of mass (to check for translational motion)
     -bias           plot the total bias potential evolution
     -metaG          plot the reconstructed free energy surface from metadynamics

     -savefig        instead of pop up windows, save a png of the plot
```
