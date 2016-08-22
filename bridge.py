import sys
import glob
import os
sys.path.append(os.path.realpath("Library"))
import protein
import sheet

MAXVALUE = 100				#Used for cases where there are... like... 1 million files
path = os.getcwd()

#Checks normal formatting
def format(filename):
	if os.path.isfile(filename) == True and filename.endswith(".pdb"):
		return True
	else:
		return False

def computeAll(filename):
	only100 = 0
	for filename in glob.glob(os.path.join(path, '*.pdb')):
		if (only100 == MAXVALUE):
			break
		else:
			only100+=1
			p = protein.buildProtein(filename)
			sheetList = sheet.buildSheet(filename, p)
			print "{0} files completed...".format(only100)
			for sheet in sheetList:
				print "===> Do Calculations Here <==="

if __name__ == '__main__':
	if (os.path.isdir(path) == True):
		path += sys.argv[1]
		print "\nComputing All Files."
		#computeAll(sys.argv[1])
		p = protein.buildProtein(sys.argv[1])
		sheetList = sheet.buildSheet(sys.argv[1], p)
		for sheet in sheetList:
			print sheet
	else:
		print "\nERROR: This folder does not exist"
		print "\nERROR: Please make sure directory format is correct"
		print "\nFORMAT: python script.py \Directory"