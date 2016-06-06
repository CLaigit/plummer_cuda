#!/bin/bash

Tempurature=(1.9 2.0 2.1 2.2 2.3 2.4 2.5 2.6 2.7 2.8)
LATTICE_LENGTH=(30)
# Input Tempurature you want to test, use space to separate different value

for l in ${LATTICE_LENGTH[@]}; do
    for t in ${Tempurature[@]}; do
        echo Tempurature: $t lattice_length $l   Outputfile: ising3dT${t}L${l}.txt
        mkdir data/3d/lattice$l
        ./src/ising3d.o $t  $l >   data/3d/lattice$l/ising3dT${t}L${l}.txt
    done
done
