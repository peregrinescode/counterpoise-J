import numpy as np
import scipy as sp
from scipy import linalg
import sys
from cclib.parser import ccopen


def CountPJ():
# Read in molecule log files for counterpoise method. Requires IOp(6/7=3) in Gaussian header + ghost atoms to make up basis sets in individual com files  
	MOLA_CP=sys.argv[1]
	MOLB_CP=sys.argv[2]
	MOLAB_CP=sys.argv[3]

# Open the log files
	molA_parser=ccopen("%s"%(MOLA_CP))
	molB_parser=ccopen("%s"%(MOLB_CP))
	molAB_parser=ccopen("%s"%(MOLAB_CP))

# Parse the relevant data
	molA=molA_parser.parse()
	molB=molB_parser.parse()
	molAB=molAB_parser.parse()

	print ("Parsed...")

# Size of basis sets
	Nbasis=molAB.nbasis

# HOMO index
	nhomoA=molA.homos
	nhomoB=molB.homos
	nhomoAB=molAB.homos

	print "HOMO A: ", nhomoA
	print "HOMO B: ", nhomoB
	print "HOMO AB: ", nhomoAB

# Every basis set should have the same size (the size of the pair) 
	if molA.nbasis!=molB.nbasis:
		print("Count of basis functions doesn't match. Failing.")
		return False

# Get molecular orbitals
	MOsA=(molA.mocoeffs[0])
	MOsB=(molB.mocoeffs[0])
	MOsAB=(molAB.mocoeffs[0])

	#print "MOsA: ", MOsA
	#print "MOsB: ", MOsB
	#print "MOsAB: ", MOsAB

# Get overlaps
	SAB=molAB.aooverlaps

	#print "Overlaps: ", SAB


# Get eigenvalues of pair
	EvalsAB=molAB.moenergies[0]

	print "Energies: ", molAB.moenergies
	print "Energies: ", EvalsAB
# Find HOMO and LUMO from energy splitting in dimer

	print "ESID HOMO-HOMO coupling", 0.5*(EvalsAB[nhomoAB]-EvalsAB[nhomoAB-1])
	print "ESID LUMO-LUMO coupling", 0.5*(EvalsAB[nhomoAB+2]-EvalsAB[nhomoAB+1])

# Calculate the molecular orbitals of A and B in the AB basis set

	MolAB_Pro = (np.dot(MOsAB,SAB)).T
	PsiA_AB_BS = np.dot(MOsA, MolAB_Pro)
	PsiB_AB_BS = np.dot(MOsB, MolAB_Pro)

	#MolAB_Pro = (np.dot(MOsAB,SAB)).T
	#PsiA_AB_BS = np.dot(MOsA, MolAB_Pro)
	#PsiB_AB_BS = np.dot(MOsB, MolAB_Pro)

	#print "MolA in basis of dimer: ", PsiA_AB_BS

	#print "MolB in basis of dimer: ", PsiB_AB_BS

# Calculate the matrix of transfer integrals
	JAB=np.dot(np.dot(np.diagflat(EvalsAB),PsiA_AB_BS),PsiB_AB_BS.T)
	print "J matrix: ", JAB

# Print the HOMO-HOMO and LUMO-LUMO coupling
	print "HOMO-HOMO coupling: ", JAB[nhomoA,nhomoB]
	print "LUMO-LUMO coupling: ", JAB[nhomoA+1,nhomoB+1]


CountPJ()



