#!/bin/bash

Tempurature=(2.3)
LATTICE_LENGTH=(256)
# Input Tempurature you want to test, use space to separate different value

for l in ${LATTICE_LENGTH[@]}; do
    for t in ${Tempurature[@]}; do
        echo Tempurature: $t lattice_length $l   Outputfile: ising3dT${t}L${l}.txt
        mkdir data/3d/lattice$l
        ./src/ising3d.o $t  $l >   data/3d/lattice$l/ising3dT${t}L${l}.txt
    done
done
