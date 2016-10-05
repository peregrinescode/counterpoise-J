# Jarvist Moore Frost 2016-10-04
# I took the opportunity of having a working codebase again to look at
# Intra-PFO transfer integrals.  
# This is for a relaxed (~23 degree) dimer. The code produces two values,
# I assume for the two symmetries of the open-shell calculation. I'm not sure
# how trustworthy these are, I think you should probably do something a bit
# more sophisticated than King Solomon in splitting the dimer!  It seems to
# agree with the 'energy splitting in dimer method' (Esplit=2*J).
#
# (( Polyfluorene Dimer; open shell doublet calculation for the fragments. Structures from:
# jmf02@login-2:/work/jmf02/absorption_pfo/absorption_pfo2/
# some Dec 2008 work I did on PFO absorption))

. ../counterpoise.sh pfoA.xyz pfoB.xyz pfo-intra
# NB: part1 and part2 fragment jobs were made into open shell doublets '0 2' calculations; by hand
g09_batch_laptop.sh */*.com
python ../ProJ.py part1/pfo-intrapart1.log part2/pfo-intrapart2.log dim/pfo-intradim.log
