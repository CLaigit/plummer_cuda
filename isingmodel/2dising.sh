#!/bin/bash

Tempurature=(1.6 1.7 1.8 1.9 1.95 2.0 2.05 2.1 2.15 2.2 2.25 2.3 2.35 2.4 2.5 2.6 2.7 2.8 2.9 3.0 3.1 3.2 3.3 3.4 3.5)
LATTICE_LENGTH=(20 30)
# Input Tempurature you want to test, use space to separate different value

for l in ${LATTICE_LENGTH[@]}; do
    for t in ${Tempurature[@]}; do
        echo Tempurature: $t lattice_length $l   Outputfile: isingT${t}L${l}.dat
        mkdir -p data/2d/lattice$l
        ./src/ising.o $t  $l >   data/2d/lattice$l/isingT${t}L${l}.dat
    done
done
