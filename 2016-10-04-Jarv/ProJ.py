import numpy as np
from cclib.parser import ccopen
import sys

def CalcJ():
    #molA_parser=ccopen("part1/"+namejob+"part1.log")
    #molB_parser=ccopen("part2/"+namejob+"part2.log")
    #molAB_parser=ccopen("dim/"+namejob+"dim.log")


	MOLA_CP=sys.argv[1]
	MOLB_CP=sys.argv[2]
	MOLAB_CP=sys.argv[3]

# Open the log files
	molA_parser=ccopen("%s"%(MOLA_CP))
	molB_parser=ccopen("%s"%(MOLB_CP))
	molAB_parser=ccopen("%s"%(MOLAB_CP))


	molA=molA_parser.parse()
	molB=molB_parser.parse()
	molAB=molAB_parser.parse()
	print ("Parsed...")

	nbs=molAB.nbasis
	if molA.nbasis!=molB.nbasis:
		print("Count of basis functions doesn't match. Failing.")
		return False

	for mole in molA,molB,molAB:
		if len(mole.atomcoords)!=1:
			print(mole," calculation appears to be an optimisation! Failing.")

	nhomoA=molA.homos
	nhomoB=molB.homos
	nhomoAB=molAB.homos

#Take mocoeffs[0] - the alpha electrons to get matrix in correct order
	MolAB_Pro = np.transpose(np.dot(molAB.mocoeffs[0],molAB.aooverlaps))
#    print "DimPro via JKP:", DimPro, "DimPro2 via JMF", DimPro2

#    print "Psi1DimBS via JKP: ", Psi1DimBS
	PsiA_DimBS = np.dot(molA.mocoeffs[0], MolAB_Pro)
	PsiB_DimBS = np.dot(molB.mocoeffs[0], MolAB_Pro)
#    print "Dim Mocoeffs: ", molAB.moenergies[0]/27.211


#Note: moenergies in eV, so converted to Hartree for checking with JKP code
	JAB=np.dot(np.dot(np.diagflat(molAB.moenergies[0]),PsiA_DimBS), np.transpose(PsiB_DimBS) )
	JAA=np.dot(np.dot(np.diagflat(molAB.moenergies[0]),PsiA_DimBS), np.transpose(PsiA_DimBS) )
	JBB=np.dot(np.dot(np.diagflat(molAB.moenergies[0]),PsiB_DimBS), np.transpose(PsiB_DimBS) )

	#print "JAB", JAB
	#print "JAA", JAA
	#print "JBB", JBB

	print "HOMO-HOMO coupling: ", JAB[nhomoA,nhomoB]
	print "LUMO-LUMO coupling: ", JAB[nhomoA+1,nhomoB+1]


	return [JAB, JAA, JBB]

CalcJ()
