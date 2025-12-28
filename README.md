# MD

## example use: simple MD simulation starting from a distorted butane structure 

### minimization (requires gfortran compiler):
cd build && ./compile.sh && cd .. && ./md_program -c coord/butane_distorted.xyz -t coord/butane.top -mcg

### md simulation (requires gfortran compiler):
cd build && ./compile.sh && cd .. && ./md_program -c coord/butane_distorted.minimized.xyz -t coord/butane.top -p

### analysis (requires NumPy and matplotlib.pyplot):
cd trajectory_analysis/ && python3 plot.py -h && cd .. <p>
cd trajectory_analysis/ && python3 plot.py ../coord/butane_distorted.minimized.properties.txt -energies && cd ..
