#!/bin/bash

list=$(ls *.pdb)

for i in $list; do
    name=${i//.pdb/}
    echo $name
    python ../topology_code/get_topology.py -p $i -o $name -m 
    mkdir $name
    cd $name
    mkdir md_prog
    cd md_prog
    mv ../../${name}.xyz ../../${name}.top .
    ../../../md_program -c ${name}.xyz -t ${name}.top -d >> output.txt
    cd ..
    mkdir obabel
    cp ../$i obabel
    cd ..
    echo "${name} done"
done

    


