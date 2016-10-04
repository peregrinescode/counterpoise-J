#! /bin/bash

for Z in 3.0 4.0 5.0 6.0 7.0 8.0 9.0 # be careful here; treated as strings!
do
    cat etheneB.xyz | sed s/3.500/${Z}/ > etheneZ.xyz

    . counterpoise.sh etheneA.xyz etheneZ.xyz ethene

    g09_batch_laptop.sh */*.com

    python ProJ.py part1/ethenepart1.log  part2/ethenepart2.log dim/ethenedim.log > ${Z}.dat
done


