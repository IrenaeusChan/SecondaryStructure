import sys
import glob
import os
import re
import math
import numpy as np
from itertools import tee, islice, chain, izip		#May not be needed
sys.path.append(os.path.realpath("Library"))
import protein
import sheet
import vector

MAXVALUE = 500				#Used for cases where there are... like... 1 million files
ELEMENTS = {'N':14, 'C':12, 'O':16, 'S':32, 'H':1, 'P':31, 'D':1}

def detectBetaBarrel(listOfSheets):
	listOfBetaBarrel = []
	for everySheet in listOfSheets:
		center, centroids = centerOfSheet(everySheet)		#Grab the center and centroids
		shearNumber = shearNumberCalculation(everySheet)
		calcRadius = radiusEquation(shearNumber, everySheet.totalStrand)
		roughRadius = roughRadiusCalculations(center, centroids)

		print shearNumber
		print calcRadius
		print roughRadius
		
		if abs(calcRadius - roughRadius) < 3:
			#listOfBetaBarrel.append((str(everySheet), str(calcRadius), str(roughRadius)))
			listOfBetaBarrel.append(everySheet)
	return listOfBetaBarrel

def radiusEquation(shearNumber, numberOfStrands):
	#Adapted from Principles determining the structures of B-sheet barrels in proteins I. A theoretical analysis
	return (math.sqrt((shearNumber*3.3)**2 + (numberOfStrands*4.4)**2))/(2*numberOfStrands*math.sin(math.pi/numberOfStrands))
	#return (4.4)/(2*math.sin(math.pi/numberOfStrands)*math.cos(angle))

#This has too many problems with it. We need something else
def shearNumberCalculation(sheet):
	#Adapted from Shear number of protein B-barrels; Definition refinements and statistics by Wei-Min Liu
	closedPath = []
	shearNumber, direction = 0, 0
	for i, everyStrand in enumerate(sheet.strandList):
		if everyStrand.strandType == 0: 
			direction = 1 	#The direction is used to determine what the next strand direction is
			continue
		else:																#ThisSheet PreviousSheet Result
			if everyStrand.strandType == -1 and direction == 1: direction = -1 		#Anti 		-->			 <--
			elif everyStrand.strandType == -1 and direction == -1: direction = 1 	#Anti 		<-- 		 -->
			elif everyStrand.strandType == 1 and direction == 1: direction = 1 		#Parallel	--> 		 -->
			elif everyStrand.strandType == 1 and direction == -1: direction = -1 	#Parallel	<-- 		 <--
			closedPath.append([everyStrand.otherStrand, everyStrand.thisStrand, direction])
			#With this we are creating a list of how the strands connect to one another and the direction of the strands
	#print closedPath
	if not detectBridge(sheet):
		firstAA, lastAA = hbond(sheet)
		if firstAA is not None and lastAA is not None:
			closedPath.append([lastAA.position, firstAA.position, direction])
	for i in range(0, len(closedPath)):
		#According to the paper, this would be the first calculation, but we are making it last due to how
		# our list is currently ordered. For example if the series of connections was
		# (118, 217; 216, 192; 194, 137; 130, 239; 244, 119)
		# then the generated list would look like:
		# 118, 217, -1
		# 216, 192, -1
		# 194, 137, -1
		# 130, 239, 1
		# 244, 119, 1
		if i == len(closedPath)-1: 
			shearNumber += closedPath[0][0] - closedPath[-1][1]		#Using our example, this would be +(118 - 119)
		else:
			if closedPath[i][2] == -1: 
				shearNumber -= (closedPath[i+1][0] - closedPath[i][1])
				#These would be for:
				# -(216 - 217)
				# -(194 - 192)
				# -(130 - 137)
			elif closedPath[i][2] == 1: 
				shearNumber += (closedPath[i+1][0] - closedPath[i][1])
				#These would be for:
				# +(244 - 239)
		#Put together it would be (118 - 119) - (216 - 217) - (194 - 192) - (130 - 137) + (244 - 239)
		# Which is exactly how it should be according to the paper resulting in 10.
	return shearNumber

def detectBridge(wholeSheet):
	firstStrand = wholeSheet.strandList[0]
	lastStrand = wholeSheet.strandList[-1]
	if (((firstStrand.start-1) == lastStrand.stop) or ((firstStrand.stop+1) == lastStrand.start)):
		return True
	elif ((firstStrand.start >= lastStrand.start and firstStrand.start <= lastStrand.stop) \
		or (firstStrand.stop >= lastStrand.start and firstStrand.stop <= lastStrand.stop)):
		return True
	elif ((lastStrand.start >= firstStrand.start and lastStrand.start <= firstStrand.stop) \
		or (lastStrand.stop >= firstStrand.start and lastStrand.stop <= firstStrand.stop)):
		return True

def centerOfSheet(sheet):
	#Adapted from Shear number of protein B-barrels; Definition refinements and statistics by Wei-Min Liu
	totalX, totalY, totalZ, totalMass = 0, 0, 0, 0
	listOfStrandCentroids = []

	for everyStrand in sheet.strandList:
		for everyAminoAcid in everyStrand.aminoAcidList:
			for everyAtom in everyAminoAcid.backboneAtoms:
				totalX += everyAtom.x * ELEMENTS[everyAtom.element]
				totalY += everyAtom.y * ELEMENTS[everyAtom.element]
				totalZ += everyAtom.z * ELEMENTS[everyAtom.element]
				totalMass += ELEMENTS[everyAtom.element]
		if totalMass != 0:
			totalX = totalX/totalMass
			totalY = totalY/totalMass
			totalZ = totalZ/totalMass
			listOfStrandCentroids.append((totalX, totalY, totalZ))
		totalX, totalY, totalZ, totalMass = 0, 0, 0, 0					#Forgot to make TotalMass 0...
	for centroids in listOfStrandCentroids:
		totalX += centroids[0]
		totalY += centroids[1]
		totalZ += centroids[2]
	totalX = totalX/len(listOfStrandCentroids)
	totalY = totalY/len(listOfStrandCentroids)
	totalZ = totalZ/len(listOfStrandCentroids)
	center = (totalX, totalY, totalZ)
	return center, listOfStrandCentroids

def roughRadiusCalculations(center, centroids):
	sumOfMagnitudes = 0
	for centroid in centroids:
		sumOfMagnitudes += vector.vectorMagnitude(vector.vectorCalculation(center, centroid))
	return sumOfMagnitudes/len(centroids)

def hbond(sheet):
	#Used to detect the presence of Hydrogen Bonds that may exist between amino acids between two ajoining B-Strands
	hbondList = []
	firstStrandAA, lastStrandAA = None, None
	print sheet.PDBFile
	for aStrand in sheet.strandList:
		#Basically if it's the first Strand, there won't be any information regarding the bonds
		if aStrand.strandType != 0 and (aStrand.otherStrand != 0 or aStrand.thisStrand != 0):
			for aa in prevStrand.aminoAcidList:
				if aStrand.otherStrand == aa.position: prevAA = aa 		#We are trying to save the AA
			for aa in aStrand.aminoAcidList:
				if aStrand.thisStrand == aa.position: thisAA = aa
			#In BetaSheets, the hydrogen bond typically goes from the Nitrogen to the Oxygen on the Carboxyl Group
			# But the hydrogen bond can be in both directions, so we just want to check all possible bond configurations
			NFromPrev = (prevAA.backboneAtoms[0].x, prevAA.backboneAtoms[0].y, prevAA.backboneAtoms[0].z)
			OFromPrev = (prevAA.sidechainAtoms[0].x, prevAA.sidechainAtoms[0].y, prevAA.sidechainAtoms[0].z)
			OFromCurrent = (thisAA.sidechainAtoms[0].x, thisAA.sidechainAtoms[0].y, thisAA.sidechainAtoms[0].z)
			NFromCurrent = (thisAA.backboneAtoms[0].x, thisAA.backboneAtoms[0].y, thisAA.backboneAtoms[0].z)
			#Which is why we want to look for the smallest bond distance which would indicate is the "actual" hydrogen bond
			hbond = min(vector.vectorMagnitude(vector.vectorCalculation(NFromPrev, OFromCurrent)), \
				vector.vectorMagnitude(vector.vectorCalculation(OFromPrev, NFromCurrent)))
			hbondList.append(hbond)		#Grab all the hydrogen bonds throughout the entire Strands
		prevStrand = aStrand

	if len(hbondList) > 1:	#If the sheet is only 2 Strands, it won't have enough values to give proper std
		lowerLimit = np.mean(hbondList) - np.std(hbondList, ddof = 1)
		upperLimit = np.mean(hbondList) + np.std(hbondList, ddof = 1)
	else: upperLimit, lowerLimit = 0,0

	print lowerLimit
	print upperLimit
	print ""
	firstStrand = sheet.strandList[0]
	lastStrand = sheet.strandList[-1]
	#Now we possess the average Hbond that occurs witin the structure as well as the standard deviations
	# we will use this information to attempt to determine if there is a hydrogen bond between the first and last
	# strand. If a hydrogen bond exists between the first and last strand, then we can assume that a complete
	# circle is created and then we can safely assume it is a beta barrel
	for aa in firstStrand.aminoAcidList:		#Look at the First Strand
		NFromFirst = (aa.backboneAtoms[0].x, aa.backboneAtoms[0].y, aa.backboneAtoms[0].z)
		OFromFirst = (aa.sidechainAtoms[0].x, aa.sidechainAtoms[0].y, aa.sidechainAtoms[0].z)
		for bb in lastStrand.aminoAcidList:		#Compare it with the AA in Last Strand
			NFromLast = (bb.backboneAtoms[0].x, bb.backboneAtoms[0].y, bb.backboneAtoms[0].z)
			OFromLast = (bb.sidechainAtoms[0].x, bb.sidechainAtoms[0].y, bb.sidechainAtoms[0].z)
			
			#If the bond distances between the FirstAA and LastAA are within the parameters, then we can assume
			# that it is a possible hydrogen bond as determined by our whole structure
			hbond = min(vector.vectorMagnitude(vector.vectorCalculation(NFromFirst, OFromLast)), \
				vector.vectorMagnitude(vector.vectorCalculation(OFromFirst, NFromLast)))
			print hbond
			print hbond < upperLimit
			print hbond > lowerLimit
			print ""
			if (hbond < upperLimit) and (hbond > lowerLimit):
				lastStrandAA = bb
				break
		if hbond < upperLimit and hbond > lowerLimit:
			firstStrandAA = aa
			break
	return firstStrandAA, lastStrandAA		
	#Return the AA that we suspect the bond exists, which will help us complete the "Path"

def computeAll(path):
	only100 = 0
	count = 0
	for filename in glob.glob(os.path.join(path, '*.pdb')):
		drive, pathAndFile = os.path.splitdrive(filename)			#http://stackoverflow.com/questions/3167154/how-to-split-a-dos-path-into-its-components-in-python
		filePath, file = os.path.split(pathAndFile)
		if (only100 == MAXVALUE):
			break
		else:
			only100+=1
			p = protein.buildProtein(filename)
			sheetList = sheet.buildSheet(filename, p)
			print "{0} files completed...".format(only100)
			listOfBetaBarrel = detectBetaBarrel(sheetList)
			if listOfBetaBarrel:
				print "{0} has a Beta Barrel".format(file)
				toHTML(listOfBetaBarrel, file)
				count+=1
			else:
				print "{0} is NOT a Beta Barrel".format(file)
	with open("twist.html", 'a') as output:						#My quick solution to ending HTML file
		output.write("</body>\n</html>")
	print count

def computeOne(filename):
	p = protein.buildProtein(filename)
	sheetList = sheet.buildSheet(filename, p)
	listOfBetaBarrel = detectBetaBarrel(sheetList)
	if listOfBetaBarrel:
		print "{0} has a Beta Barrel".format(filename)
		toHTML(listOfBetaBarrel, filename)
	else:
		print "{0} is NOT a Beta Barrel".format(filename)
	with open("twist.html", 'a') as output:						#My quick solution to ending HTML file
		output.write("</body>\n</html>")

def toHTML(listOfBetaBarrel, filename):
	with open('betabarrel.html', 'a') as output:
		if os.stat("betabarrel.html").st_size == 0:
			output.write("<html>\n<head>\n</head>\n<body>")
		output.write("<pre>" + filename + "</pre>")
		for theSheet in listOfBetaBarrel:
			output.write("<pre>" + "Theoretical: " + theSheet[1] + "</br>" + "Empirical: " + theSheet[2] + "</pre>")
			output.write("<pre>")
			for i, line in enumerate(theSheet[0].splitlines()):
				output.write(line + "</br>") 
			output.write("</pre>")