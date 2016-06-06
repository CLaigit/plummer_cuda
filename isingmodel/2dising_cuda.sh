#!/bin/bash

Tempurature=(2)
# Input Tempurature you want to test, use space to separate different value
for a in ${Tempurature[@]}; do
    echo Tempurature: $a    Outputfile: isingT$a.dat
    ./ising2d.o $a  >   ising_cuda_T$a.dat
done
