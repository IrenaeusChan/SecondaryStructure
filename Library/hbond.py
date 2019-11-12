import sys
import glob
import os
import string
import numpy as np
sys.path.append(os.path.realpath("Library"))
import protein
import sheet
from vector import *

hbondList = []

#Example
print np.std([1.5, 2.5, 2.5, 2.75, 3.25, 4.75], ddof=1)

#Checks normal formatting
def format(filename):
	if os.path.isfile(filename) == True and filename.endswith(".pdb"):
		return True
	else:
		return False

if format(sys.argv[1]):
	filename = sys.argv[1]
	print "\nComputing Using File: {0}".format(sys.argv[1])

p = protein.buildProtein(filename)
sheetList = sheet.buildSheet(filename, p)

for theSheet in sheetList:
	for aStrand in theSheet.strandList:
		
		if aStrand.thisStrand != None:
			print ""
			print prevStrand
			print aStrand
			print ""
			print aStrand.otherStrand
			print aStrand.thisStrand
			for aa in prevStrand.aminoAcidList:
				if aStrand.otherStrand == aa.position:
					prevAA = aa
			for aa in aStrand.aminoAcidList:
				if aStrand.thisStrand == aa.position:
					thisAA = aa
			NFromPrev = (prevAA.backboneAtoms[0].x, prevAA.backboneAtoms[0].y, prevAA.backboneAtoms[0].z)
			OFromPrev = (prevAA.sidechainAtoms[0].x, prevAA.sidechainAtoms[0].y, prevAA.sidechainAtoms[0].z)
			OFromCurrent = (thisAA.sidechainAtoms[0].x, thisAA.sidechainAtoms[0].y, thisAA.sidechainAtoms[0].z)
			NFromCurrent = (thisAA.backboneAtoms[0].x, thisAA.backboneAtoms[0].y, thisAA.backboneAtoms[0].z)
			print ""
			print vectorMagnitude(vectorCalculation(NFromPrev, OFromCurrent))
			print vectorMagnitude(vectorCalculation(OFromPrev, NFromCurrent))
			print ""
			hbond = min(vectorMagnitude(vectorCalculation(NFromPrev, OFromCurrent)), vectorMagnitude(vectorCalculation(OFromPrev, NFromCurrent)))
			hbondList.append(hbond)
			print hbond
			print ""
		prevStrand = aStrand
	print hbondList
	print ""
	print np.mean(hbondList)
	print np.mean(hbondList) - np.std(hbondList, ddof = 1)
	print np.mean(hbondList) + np.std(hbondList, ddof = 1)

	hbondList = []