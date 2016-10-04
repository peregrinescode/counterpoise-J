# Code to calculate transfer integrals via a Counterpoise (Gaussian 'ghost' atoms) method
# Original code by James Kirkpatrick 2009
# Reimplemented with CCLIB by Jarvist Moore frost 2009
# Addition by Bjorn Baumeier to calculate JAB_eff (Eq.10 in JACS 128, 9884 (2006)), 2009
# Aggregation of various code bases and attempt to figure out what on earth is going on; Beth Rice and Jarvist Moore Frost 2016.

import numpy as np
from cclib.parser import ccopen
import sys

def CalcJ():
    #
    #molA_parser=ccopen("part1/"+namejob+"part1.log")
    #molB_parser=ccopen("part2/"+namejob+"part2.log")
    #molAB_parser=ccopen("dim/"+namejob+"dim.log")

    MOLA_CP=sys.argv[1]
    MOLB_CP=sys.argv[2]
    MOLAB_CP=sys.argv[3]
# Open the Gaussian 'log' files with CCLIB
    molA_parser=ccopen("%s"%(MOLA_CP))
    molB_parser=ccopen("%s"%(MOLB_CP))
    molAB_parser=ccopen("%s"%(MOLAB_CP))

# Parse the logfiles into mol{A,B,AB} objects
    molA=molA_parser.parse()
    molB=molB_parser.parse()
    molAB=molAB_parser.parse()
    print ("Gaussian logfiles Loaded and Parsed...")

    nbs=molAB.nbasis
    if molA.nbasis!=molB.nbasis:
        print("Number of basis functions in molA and molB. Failing.")
        return False

    for mole in molA,molB,molAB:
        if len(mole.atomcoords)!=1:
            print(mole," calculation appears to be an optimisation (multiple steps)! Warning!")

    nhomoA=molA.homos
    nhomoB=molB.homos
    nhomoAB=molAB.homos

    print "Energy splitting in dimer"
    print "HOMO-HOMO coupling: ",(molAB.moenergies[0][nhomoAB]-molAB.moenergies[0][nhomoAB-1])/2.0," eV"
    print "LUMO-LUMO coupling: ",(molAB.moenergies[0][nhomoAB+2]-molAB.moenergies[0][nhomoAB+1])/2.0," eV"

#Take mocoeffs[0] - the alpha electrons, to get matrix in correct order
    MolAB_Pro = np.transpose(np.dot(molAB.mocoeffs[0],molAB.aooverlaps))
   
    PsiA_DimBS = np.dot(molA.mocoeffs[0], MolAB_Pro)
    PsiB_DimBS = np.dot(molB.mocoeffs[0], MolAB_Pro)
#    print "Dim Mocoeffs: ", molAB.moenergies[0]/27.211

# JMF 2009-08 CCLIB code to reproduce JKP's original values for Ethene. BELIEVED INCORRECT.
    print "Original James KP method; matrix multiplication believed out of order."
    JAB=np.dot(np.dot(np.diagflat(molAB.moenergies[0]),PsiA_DimBS), np.transpose(PsiB_DimBS) )
    JAA=np.dot(np.dot(np.diagflat(molAB.moenergies[0]),PsiA_DimBS), np.transpose(PsiA_DimBS) )
    JBB=np.dot(np.dot(np.diagflat(molAB.moenergies[0]),PsiB_DimBS), np.transpose(PsiB_DimBS) )

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

    #print "JAB", JAB
    #print "JAA", JAA
    #print "JBB", JBB

# I think these are hardwired by James to look at ethene molecular orbital overlaps
    # print "PSI1 in dimer", PsiA_DimBS[3][0], PsiA_DimBS[3][1], PsiA_DimBS[3][2], PsiA_DimBS[3][3] , PsiA_DimBS[3][4], PsiA_DimBS[3][5], PsiA_DimBS[3][6], PsiA_DimBS[3][7]
    # print "PSI1 in dimer", PsiB_DimBS[3][0], PsiB_DimBS[3][1], PsiB_DimBS[3][2], PsiB_DimBS[3][3] , PsiB_DimBS[3][4], PsiB_DimBS[3][5], PsiB_DimBS[3][6], PsiB_DimBS[3][7]

    print "HOMO-HOMO coupling: ", JAB[nhomoA,nhomoB], " eV"
    print "LUMO-LUMO coupling: ", JAB[nhomoA+1,nhomoB+1], " eV"
#Note: CCLIB values in eV, so converted to Hartree for checking with JKP original code
    print "HOMO-HOMO coupling: ", JAB[nhomoA,nhomoB]/27.211, " Ha"
    print "LUMO-LUMO coupling: ", JAB[nhomoA+1,nhomoB+1]/27.211, " Ha"

# JMF 2016-10-04: These additions are I think due to Bjorn Baumeier 2009-09-17
# Determine overlap analogous to JAB
    SAB = np.dot(PsiB_DimBS , np.transpose(PsiA_DimBS) )
    # Calculate JAB_eff according to Eq.10 in JACS 128, 9884 (2006)
    # !only for the desired orbitals!

    print "Bjorn Baumeir JAB_eff (Eq.10 in JACS 128, 9884 (2006)) "
    orbA=nhomoA
    orbB=nhomoB
    JAB_eff = (JAB[orbA,orbB] - 0.5*(JAA[orbA,orbA]+JBB[orbB,orbB])*SAB[orbA,orbB])/(1.0 - SAB[orbA,orbB]*SAB[orbA,orbB])
    print "HOMO-HOMO Jeff-coupling", JAB_eff," eV"

    orbA=nhomoA+1
    orbB=nhomoB+1
    JAB_eff = (JAB[orbA,orbB] - 0.5*(JAA[orbA,orbA]+JBB[orbB,orbB])*SAB[orbA,orbB])/(1.0 - SAB[orbA,orbB]*SAB[orbA,orbB])
    print "LUMO-LUMO Jeff-coupling", JAB_eff, " eV"

    return [JAB, JAA, JBB]

CalcJ()
