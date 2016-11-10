import sys
import glob
import os
import string
sys.path.append(os.path.realpath("Library"))
import protein
import sheet

MAXVALUE = 500				#Used for cases where there are... like... 1 million files
path = os.getcwd()

#Checks normal formatting
def format(filename):
	if os.path.isfile(filename) == True and filename.endswith(".pdb"):
		return True
	else:
		return False

def detectBridge(listOfSheets):
	for i, currentSheet in enumerate(listOfSheets):
		for otherSheet in listOfSheets[i+1:len(listOfSheets)]:
			if otherSheet.seqres != currentSheet.seqres:
				continue
			for currentStrand in currentSheet.strandList:
				for otherStrand in otherSheet.strandList:
					if (((currentStrand.start-1) == otherStrand.stop) or ((currentStrand.stop+1) == otherStrand.start)):
						print "First"
						print "{0} {1}{2}".format(currentSheet.seqres, currentSheet.sheetIden, currentStrand.strandNum)
						print "{0} {1}{2}".format(otherSheet.seqres, otherSheet.sheetIden, otherStrand.strandNum)
						return True
					elif ((currentStrand.start >= otherStrand.start and currentStrand.start <= otherStrand.stop) 
						or (currentStrand.stop >= otherStrand.start and currentStrand.stop <= otherStrand.stop)):
						print "Second"
						print "{0} {1}{2}".format(currentSheet.seqres, currentSheet.sheetIden, currentStrand.strandNum)
						print "{0} {1}{2}".format(otherSheet.seqres, otherSheet.sheetIden, otherStrand.strandNum)
						return True
					elif ((otherStrand.start >= currentStrand.start and otherStrand.start <= currentStrand.stop) 
						or (otherStrand.stop >= currentStrand.start and otherStrand.stop <= currentStrand.stop)):
						print "Third"
						print "{0} {1}{2}".format(currentSheet.seqres, currentSheet.sheetIden, currentStrand.strandNum)
						print "{0} {1}{2}".format(otherSheet.seqres, otherSheet.sheetIden, otherStrand.strandNum)
						return True
	return False

def computeAll(filename):
	only100 = 0
	count = 0
	for filename in glob.glob(os.path.join(path, '*.pdb')):
		drive, pathAndFile = os.path.splitdrive(filename)			#http://stackoverflow.com/questions/3167154/how-to-split-a-dos-path-into-its-components-in-python
		filePath, file = os.path.split(pathAndFile)
		if (only100 == MAXVALUE):
			break
		else:
			only100+=1
			#print filename
			p = protein.buildProtein(filename)
			sheetList = sheet.buildSheet(filename, p)
			print "{0} files completed...".format(only100)
			if detectBridge(sheetList) == True:
				print "{0} is BRIDGE".format(file)
				count+=1
			else:
				print "{0} is NOT BRIDGE".format(file)
	return count
				
def computeOne(filename):
	p = protein.buildProtein(filename)
	sheetList = sheet.buildSheet(filename, p)
	if detectBridge(sheetList) == True:
		print "{0} is BRIDGE".format(filename)
	else:
		print "{0} is NOT BRIDGE".format(filename)
	for sheet1 in sheetList:
		print sheet1

if __name__ == '__main__':
	if (os.path.isdir(path) == True):
		path += sys.argv[1]
		print "\nComputing All Files."
		amount = computeAll(sys.argv[1])
		#computeOne(sys.argv[1])
		print amount
	else:
		print "\nERROR: This folder does not exist"
		print "\nERROR: Please make sure directory format is correct"
		print "\nFORMAT: python script.py \Directory"