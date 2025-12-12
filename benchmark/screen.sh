#!/bin/bash

list=$(ls *.pdb)

for i in $list; do
    name=${i//.pdb/}
    echo $name
#    python ../topology_code/get_topology.py -p $i -o $name -m 
#    mkdir $name
    cd $name
#    mkdir md_prog
    cd obabel
    python ../../pybel_calc.py $i >> output.txt
    cd ../..
    echo "${name} done"
done

    


