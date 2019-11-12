import sys
import glob
import os
import itertools
sys.path.append(os.path.realpath("Library"))
import bridge, twotwo, twist, betabarrel, protein, sheet, cylindrical, small, gradual, bulge, generic

MAXVALUE = 500				#Used for cases where there are... like... 1 million files

#Creates a dictionary of function pointers. This will be used to call the functions found in the separate modules
STRUCTURES = {'Cylindrical':cylindrical.detectCylindrical, 'BetaBarrel':betabarrel.detectBetaBarrel, 'BBarrelExtension':3, \
'UBarrel':4, 'UBarrelExtension':5, 'Spiral':6, 'Flat':7, 'Small':small.detectSmall, 'Twotwo':twotwo.detectTwoTwo, 'Large':10, 'Twist':twist.detectTwist, \
'Bulge':bulge.detectBulge, 'Gradual':gradual.detectGradual, 'Bridge':bridge.detectBridge, 'Generic':generic.detectGeneric}

#Checks normal formatting
def format(filename):
	if os.path.isfile(filename) == True and filename.endswith(".pdb"):
		return True
	else:
		print "\nERROR: Your file <{0}> either is incorrect or does not exist".format(filename)
		return False

def compute(setChoice):
	sheetList = []
	filename = sys.argv[2]
	if setChoice.lower() in list(k.lower() for k in STRUCTURES.keys()):
		functionToBeUsed = dict((key.lower(), value) for key, value in STRUCTURES.iteritems())[setChoice.lower()]
	else: 
		print "\nERROR: Choice of SET was unavailable"
		return
	if (os.path.isdir(path) == True):
		only_100 = 0
		for filename in glob.glob(os.path.join(path, '*.pdb')):
			if (only_100 == MAXVALUE): break
			else:
				only_100+=1
				print filename
				p = protein.buildProtein(filename)
				sheetList += sheet.buildSheet(filename, p)
	else:
		p = protein.buildProtein(filename)
		sheetList = sheet.buildSheet(filename, p)
	return functionToBeUsed(sheetList)
		
if __name__ == '__main__':
	path = ''
	if len(sys.argv) < 2:
		print "\nERROR: No file was provided"
		print "\nFORMAT: File[.pdb]"
		print "\nHELP: Type -h"
		sys.exit(1)
	elif (sys.argv[1] == "-c" and len(sys.argv) > 2 and format(sys.argv[2])):
		if len(sys.argv) > 3: sheetList = compute(sys.argv[3])
		if not sheetList:
			print "\nThe set <{0}> was not detected in file <{1}>".format(sys.argv[3], sys.argv[2])
		else: 
			for sheet in sheetList: print sheet
	elif (sys.argv[1] == "-a" and len(sys.argv) > 2):
		path = os.getcwd() + sys.argv[2]
		if len(sys.argv) > 3: sheetList = compute(sys.argv[3])
		if not sheetList:
			print "\nThe set <{0}> was not detected in directory <{1}>".format(sys.argv[3], sys.argv[2])
		else: 
			for sheet in sheetList: 
				print "{0} {1}".format(sheet.PDBFile, sheet.sheetIden)
			print len(sheetList)
	elif (sys.argv[1] == "-p" and len(sys.argv) > 2):
		p = protein.buildProtein(sys.argv[2])
		sheetList = sheet.buildSheet(sys.argv[2], p)
		with open('print.html', 'a') as output:
			if os.stat("print.html").st_size == 0:
				output.write("<html>\n<head>\n</head>\n<body>")
			output.write("<pre>" + sys.argv[2] + "</pre>")
			for theSheet in sheetList:
				output.write("<pre>")
				for i, line in enumerate(str(theSheet).splitlines()):
					output.write(line + "</br>") 
				output.write("</pre>")
			output.write("</body>\n</html>")
	elif (sys.argv[1] == "-h") or len(sys.argv) >= 2 and len(sys.argv) <= 3:
		print "\nPossible Commands"
		print "------------------"
		print "1. Compute Single File\t\tpython structures.py -c File[.pdb]"
		print "2. Compute Multiple Files\tpython structures.py -a PathToFiles"
		print "3. Print Sheets in File\t\tpython structures.py -p File[.pdb]"
