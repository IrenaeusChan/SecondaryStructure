"""
Irenaeus Chan
5/11/2016

Sheet Angles

Checklist:
- Evaluate a systematic approach to geometrically organizing the information based on sequences
"""
import math
import vector

def printSheetBackbone(sheetBackbone):
	with open('output.txt', "w") as out:
		out.write("x,y,z\n")
		for atom in sheetBackbone:
			out.write(str(atom[0])+","+str(atom[1])+","+str(atom[2])+"\n")

#Code to convert the SheetList into a single List of Coordinates
def organizeSheet (sheet):
	sheetBackbone = []
	centerAtoms = []
	listOfCoords = []
	totalx, totaly, totalz = 0, 0, 0
	#position = []					Not sure if we need these yet...
	#seqres = []
	for AA in sheet.amino_acids:
		for atom in AA.backbone:
			totalx += atom.x
		 	totaly += atom.y
		 	totalz += atom.z
		 	sheetBackbone.append(atom)
		 	listOfCoords.append([atom.x, atom.y, atom.z])
		totalx /= 3
		totaly /= 3
		totalz /= 3
		centerAtoms.append([totalx, totaly, totalz])
		totalx, totaly, totalz = 0, 0, 0
	return sheetBackbone, centerAtoms, listOfCoords

def angleCalculation(vectorList, atomList, listType, filename):
	filePrep = filename.split('.')
	with open('{0}sheet.txt'.format(filePrep[0]), "a") as out:
		out.write(listType + "\n")
		for i in range(0,len(vectorList)-1):
			angle = vector.dihedralAngle(vectorList[i], vectorList[i+1])
			out.write(str(angle) + " " + str(atomList[i].atom) + " " + str(atomList[i].position) + " " + str(atomList[i+1].atom) + " " + str(atomList[i+1].position) + "\n")
		out.write("\n")

def evaluateAngles (sheet, filename):
	topAtoms, bottomAtoms = [], []
	sheetBackbone, centerAtoms, listOfCoords = organizeSheet(sheet)
	#printSheetBackbone(sheetBackbone)
	#printSheetBackbone(centerAtoms)
	#regressionVector, regressionPoint = vector.orthogonalDistanceRegression(centerAtoms)
	regressionVector, regressionPoint = vector.orthogonalDistanceRegression(listOfCoords)
	for pos, atom in enumerate(sheetBackbone):
		orthogonalVector = vector.orthogonalVectorCalculation(regressionPoint, regressionVector, [atom.x, atom.y, atom.z])
		if pos % 2 == 0:			#if the number is even
			topAtoms.append([orthogonalVector, atom])
		else:
			bottomAtoms.append([orthogonalVector, atom])

	angleCalculation([i[0] for i in topAtoms], [a[1] for a in topAtoms], "top", filename)
	angleCalculation([i[0] for i in bottomAtoms], [a[1] for a in bottomAtoms], "bottom", filename)