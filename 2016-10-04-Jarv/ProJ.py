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


    JAB=np.dot(np.dot(np.diagflat(molAB.moenergies[0]),PsiA_DimBS), np.transpose(PsiB_DimBS) )
    JAA=np.dot(np.dot(np.diagflat(molAB.moenergies[0]),PsiA_DimBS), np.transpose(PsiA_DimBS) )
    JBB=np.dot(np.dot(np.diagflat(molAB.moenergies[0]),PsiB_DimBS), np.transpose(PsiB_DimBS) )

    #print "JAB", JAB
    #print "JAA", JAA
    #print "JBB", JBB

    print "HOMO-HOMO coupling: ", JAB[nhomoA,nhomoB], " eV"
    print "LUMO-LUMO coupling: ", JAB[nhomoA+1,nhomoB+1], " eV"
#Note: CCLIB values in eV, so converted to Hartree for checking with JKP original code
    print "HOMO-HOMO coupling: ", JAB[nhomoA,nhomoB]/27.211, " Ha"
    print "LUMO-LUMO coupling: ", JAB[nhomoA+1,nhomoB+1]/27.211, " Ha"

# JKP - 2009-09-10 update (script pasted into email)
    print("James KP 2009-09-10 re-order matrix multiplication...")
# this is the bit I changed: note the different order of multiplication
# I suspect that most mistakes will tend to be of that kind... x.y != y.x alas...

    JAB = np.dot(np.dot( PsiB_DimBS, np.diagflat(molAB.moenergies[0])) , np.transpose(PsiA_DimBS) )
    JAA = np.dot(np.dot( PsiA_DimBS, np.diagflat(molAB.moenergies[0])) , np.transpose(PsiA_DimBS) )
    JBB = np.dot(np.dot( PsiB_DimBS, np.diagflat(molAB.moenergies[0])) , np.transpose(PsiB_DimBS) )

# I think these are hardwired by James to look at ethene molecular orbital overlaps
   # print "PSI1 in dimer", PsiA_DimBS[3][0], PsiA_DimBS[3][1], PsiA_DimBS[3][2], PsiA_DimBS[3][3] , PsiA_DimBS[3][4], PsiA_DimBS[3][5], PsiA_DimBS[3][6], PsiA_DimBS[3][7]
   # print "PSI1 in dimer", PsiB_DimBS[3][0], PsiB_DimBS[3][1], PsiB_DimBS[3][2], PsiB_DimBS[3][3] , PsiB_DimBS[3][4], PsiB_DimBS[3][5], PsiB_DimBS[3][6], PsiB_DimBS[3][7]

    print "HOMO-HOMO coupling: ", JAB[nhomoA,nhomoB], " eV"
    print "LUMO-LUMO coupling: ", JAB[nhomoA+1,nhomoB+1], " eV"
#Note: CCLIB values in eV, so converted to Hartree for checking with JKP original code
    print "HOMO-HOMO coupling: ", JAB[nhomoA,nhomoB]/27.211, " Ha"
    print "LUMO-LUMO coupling: ", JAB[nhomoA+1,nhomoB+1]/27.211, " Ha"

# These are I think due to Bjorn Baumeir 2009-09-17

# Determine overlap analogous to JAB
    SAB = np.dot(PsiB_DimBS , np.transpose(PsiA_DimBS) )
    # Calculate JAB_eff according to Eq.10 in JACS 128, 9884 (2006)
    # !only for the desired orbitals!

    print "Bjorn Baumeir JAB_eff (Eq.10 in JACS 128, 9884 (2006)) "
    orbA=nhomoA
    orbB=nhomoB
    JAB_eff = (JAB[orbA,orbB] - 0.5*(JAA[orbA,orbA]+JBB[orbB,orbB])*SAB[orbA,orbB])/(1.0 - SAB[orbA,orbB]*SAB[orbA,orbB])
    print "HOMO-HOMO - Jeff", JAB_eff

    orbA=nhomoA+1
    orbB=nhomoB+1
    JAB_eff = (JAB[orbA,orbB] - 0.5*(JAA[orbA,orbA]+JBB[orbB,orbB])*SAB[orbA,orbB])/(1.0 - SAB[orbA,orbB]*SAB[orbA,orbB])
    print "LUMO-LUMO - Jeff", JAB_eff

    return [JAB, JAA, JBB]

CalcJ()
