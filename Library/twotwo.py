import sys
import glob
import os
import re
sys.path.append(os.path.realpath("Library"))
import protein
import sheet

MAXVALUE = 500				#Used for cases where there are... like... 1 million files

def detectTwoTwo(listOfSheets):
	listOfTwoTwo = []
	for everySheet in listOfSheets:
		if everySheet.totalStrand != 4:		#Because 2 2 can only appear in Sheets with 4 Strands
			continue
		else:
			#As long as the Strands are similar in length this is possible
			if similarInLength(everySheet.strandList[0], everySheet.strandList[1]) and \
				similarInLength(everySheet.strandList[2], everySheet.strandList[3]):
				current = everySheet.strandList[2].otherStrand		#This finds where they line up according to PDB
				#We want to find the true ending of the Strand, in this case, the "true" ending is the end
				# which faces the direction of the other 2 strands, see Write-Up for more information
				if abs(current - everySheet.strandList[1].aminoAcidList[-1].position) <= 2:
					end = everySheet.strandList[1].aminoAcidList[-1].position
					if checkTwoTwo(everySheet.strandList[2], everySheet.strandList[1], current, end) == True:
						#listOfTwoTwo.append(str(everySheet))
						listOfTwoTwo.append(everySheet)
				#The first case is if the end is towards the smaller position, this is because Strands may
				# have numerical representations for their Amino Acids, but that number may not be representitive
				# of what is the head or tail or the strand itself.
				elif abs(current - everySheet.strandList[1].aminoAcidList[0].position) <= 2:
					end = everySheet.strandList[1].aminoAcidList[0].position
					if checkTwoTwo(everySheet.strandList[2], everySheet.strandList[1], current, end) == True:
						#listOfTwoTwo.append(str(everySheet))
						listOfTwoTwo.append(everySheet)
	return listOfTwoTwo

#This makes sense if you look at the written notes
def checkTwoTwo(sheetThree, sheetTwo, current, end):
	#We use this to check for directionality, based on where the larger position is, we can infer where the end
	# of the next strand will be. If the current - end is positive, then we know that the current is larger
	# which means if the next strand is antiparallel, the end of the next strand, which is opposite to the current
	# strands directionality would be in the opposite direction, meaning we are moving towards the smaller position.

	#Assume 102 is the current, and 101 is the end
	# ... 103 <-- 102 <-- 101
	#     070 --> 071 --> 072 ...
	#In this case, the end of the 2nd Strand heads towards the smaller position at 070

	#However, in this case:
	# ...101 --> 102 --> 103
	#    072 <-- 071 <-- 070 ...
	#We see that the end of the 2nd strand heads towards the larger position at 072

	#Hence, this is why these series of If functions exist, to manage this propery
	if (current - end) > 0:
		if sheetThree.strandType == -1:
			if (sheetThree.thisStrand - 1) == sheetThree.aminoAcidList[-1].position \
			or (sheetThree.thisStrand - 1) == sheetThree.aminoAcidList[0].position:
				return True
		elif sheetThree.strandType == 1:
			if (sheetThree.thisStrand + 1) == sheetThree.aminoAcidList[-1].position \
			or (sheetThree.thisStrand + 1) == sheetThree.aminoAcidList[0].position:
				return True
	#However, if it is less than, we know the current position is the smaller one, which means that if we wanted
	# to move to the end of the next strand, we would have to flip it around, and for antiparallel, instead of
	# moving towards the smaller position, we want to move towards the larger position
	elif (current - end) < 0:
		if sheetThree.strandType == -1:
			if (sheetThree.thisStrand + 1) == sheetThree.aminoAcidList[-1].position \
			or (sheetThree.thisStrand + 1) == sheetThree.aminoAcidList[0].position:
				return True
		elif sheetThree.strandType == 1:
			if (sheetThree.thisStrand - 1) == sheetThree.aminoAcidList[-1].position \
			or (sheetThree.thisStrand - 1) == sheetThree.aminoAcidList[0].position:
				return True
	#This is in case the PDB file connects the two strands via the "last" strand
	#If this is the case, then we have to assume the other end of SheetTwo will be our new "current"
	# in order for our directionality to exist
	elif (current - end) == 0:
		if sheetTwo.aminoAcidList[-1].position - end > 0:
			if sheetThree.strandType == -1:
				if (sheetThree.thisStrand - 2) == sheetThree.aminoAcidList[-1].position \
				or (sheetThree.thisStrand - 2) == sheetThree.aminoAcidList[0].position:
					return True
			elif sheetThree.strandType == 1:
				if (sheetThree.thisStrand + 2) == sheetThree.aminoAcidList[-1].position \
				or (sheetThree.thisStrand + 2) == sheetThree.aminoAcidList[0].position:
					return True
		elif sheetTwo.aminoAcidList[0].position - end < 0:
			if sheetThree.strandType == -1:
				if (sheetThree.thisStrand + 2) == sheetThree.aminoAcidList[-1].position \
				or (sheetThree.thisStrand + 2) == sheetThree.aminoAcidList[0].position:
					return True
			elif sheetThree.strandType == 1:
				if (sheetThree.thisStrand - 2) == sheetThree.aminoAcidList[-1].position \
				or (sheetThree.thisStrand - 2) == sheetThree.aminoAcidList[0].position:
					return True
	#However, if it is connected 2 away from the end of SheetTwo, we want to see if that position on SheetThree
	# is the end. If it is the end then we don't have to worry that much, however if isn't the end, then we have
	# to check if the distance between the current and the true end is within the third sheet range

	#This is explained badly, but I can't remember why I wrote this part...
	if (sheetThree.thisStrand == sheetThree.aminoAcidList[-1].position \
		and current - end != 0) or (sheetThree.thisStrand == sheetThree.aminoAcidList[0].position \
		and current - end != 0):
		if current - end < 0:
			value = abs(current - end)
			if sheetThree.thisStrand-value in range(sheetThree.start, sheetThree.stop):
				return True
		elif current - end > 0:
			value = current - end
			if sheetThree.thisStrand+value in range(sheetThree.start, sheetThree.stop):
				return True
	return False

#We want to check if the strands are "similar" in length, in this case, if they are off by less than 2 amino acids
#We also want to ensure that the lengths aren't shorter than 3, because if it is, it wouldn't satisfy the property of
# 2 2's properly. This is arbitrary from visual inspection. Not that great. But meh.
def similarInLength(strandOne, strandTwo):
	if ((abs(strandOne.strandLength() - strandTwo.strandLength()) <= 2) and \
		(strandOne.strandLength() > 3 or strandTwo.strandLength() > 3)):		#I'm retarded... see Facebook post.
		return True
	return False

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
			listOfTwoTwo = detectTwoTwo(sheetList)
			if listOfTwoTwo:
				print "{0} is 2 2".format(file)
				toHTML(listOfTwoTwo, file)
				count+=1
			else:
				print "{0} is NOT 2 2".format(file)
	with open("twotwo.html", 'a') as output:						#My quick solution to ending HTML file
		output.write("</body>\n</html>")
	print count
				
def computeOne(filename):
	p = protein.buildProtein(filename)
	sheetList = sheet.buildSheet(filename, p)
	listOfTwoTwo = detectTwoTwo(sheetList)
	if listOfTwoTwo:
		print "{0} is 2 2".format(filename)
		toHTML(listOfTwoTwo, filename)
	else:
		print "{0} is NOT 2 2".format(filename)
	with open("twotwo.html", 'a') as output:						#My quick solution to ending HTML file
		output.write("</body>\n</html>")

def toHTML(listOfTwoTwo, filename):
	with open('twotwo.html', 'a') as output:
		if os.stat("twotwo.html").st_size == 0:
			output.write("<html>\n<head>\n</head>\n<body>")
		output.write("<pre>" + filename + "</pre>")
		for theSheet in listOfTwoTwo:
			output.write("<pre>")
			for i, line in enumerate(str(theSheet).splitlines()):
				output.write(line + "</br>") 
			output.write("</pre>")