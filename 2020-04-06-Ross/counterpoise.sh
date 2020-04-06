#! /bin/bash

method="b3lyp/6-31g*"
namemol1=$1
namemol2=$2
namecalc=$3

echo "usage: counterpoise.sh molA.xyz molB.xyz CalculationPrefix" 
echo "I will make directories, and compose the ghost atom lines, and Gaussian jobs to run."
echo "Currently: molA=$namemol1 molB=$namemol2 CalculationPrefix=$namecalc" 
echo "Internally specified:  method=$method"
echo "NB: at the moment you cannot specify chkpoint change by hand if needed, also _you_ gotta run the Gaussian jobs"

echo "%mem=14GB
%nprocshared=16
#p SCF(Tight,Conver=8) Integral(Grid=UltraFine) IOp(6/7=3) METHOD nosymm

autogen

0 1" > header

mkdir part1 
mkdir part2
mkdir dim

#### molA ####
sed "s:METHOD:${method}:" header > part1/${namecalc}-part1.com
cat $namemol1 >> part1/${namecalc}-part1.com
awk '{printf "%s-Bq %f %f %f \n", $1, $2, $3,$4}' $namemol2 >> part1/${namecalc}-part1.com
echo >>  part1/${namecalc}-part1.com

#### molB ####
sed "s:METHOD:${method}:" header > part2/${namecalc}-part2.com
awk '{printf "%s-Bq %f %f %f \n", $1, $2, $3,$4}' $namemol1 >> part2/${namecalc}-part2.com
cat $namemol2 >> part2/${namecalc}-part2.com
echo >>  part2/${namecalc}-part2.com

#### molecular pair / dimer ####
echo "%chk=dimer.chk" > dim/${namecalc}-dim.com
sed "s:METHOD:${method}:" header >> dim/${namecalc}-dim.com
cat $namemol1 >> dim/${namecalc}-dim.com
cat $namemol2 >> dim/${namecalc}-dim.com
echo >>  dim/${namecalc}-dim.com

# The molecular pair / dimer calc has a Link1 (Gaussian restart) to separately perform IOp(3/33=1)
echo "--Link1--
%chk=dimer.chk
%nprocshared=16
%mem=14Gb
#p geom(allcheck) guess(read,only) IOp(3/33=1) ${method} nosymm

" >> dim/${namecalc}-dim.com

