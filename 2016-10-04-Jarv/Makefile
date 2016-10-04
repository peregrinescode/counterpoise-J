# A crummy Makefile just to avoid making bash horror scripts - JMF 2016-10-04

Counterpoise_method:
	python Counterpoise_method.py  part1/ethenepart1.log  part2/ethenepart2.log dim/ethenedim.log

ProJ:
	python ProJ.py part1/ethenepart1.log  part2/ethenepart2.log dim/ethenedim.log

counterpoise:
	sh counterpoise.sh etheneA.xyz etheneB.xyz ethene

run_gaussian:
	g09_batch_laptop.sh */*.com
