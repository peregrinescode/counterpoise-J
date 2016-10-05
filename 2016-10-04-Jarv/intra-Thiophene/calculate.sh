# Jarvist Moore Frost 2016-10-05
. ../counterpoise.sh thiopheneA.xyz thiopheneB.xyz thiophene
# NB: part1 and part2 fragment jobs were made into open shell doublets '0 2' calculations; by hand
vim part?/*.com
g09_batch_laptop.sh */*.com
python ../ProJ.py part1/thiophenepart1.log part2/thiophenepart2.log dim/thiophenedim.log
