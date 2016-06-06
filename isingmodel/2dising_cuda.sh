#!/bin/bash

Tempurature=(2)
# Input Tempurature you want to test, use space to separate different value
for a in ${Tempurature[@]}; do
    echo Tempurature: $a    Outputfile: isingT$a.dat
    ./src/ising2d.o $a  >   data/ising_cuda_T$a.dat
done
