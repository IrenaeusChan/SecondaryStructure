def detectBulge(listOfSheets):
	listOfBulge = []
	for everySheet in listOfSheets:
		if everySheet.totalStrand >=5:
		#Anything smaller than 5 risks being small without much of a gradual bulge
			if everySheet.totalStrand % 2 != 0:	#The Sheet is odd
				#We need to understand here that the Middle Position is not necessarily the Middle Number
				# in the sense of the list, the Middle Position is correct. However in the sense of
				# realistic representation the Middle Position is actually ONE less than the proper value
				# If we are trying to model starting at 0, middlePos is correct, but if we are modelling
				# based on a starting 1 strand, then we have to remember to add 1 at all times.
				middlePos = (everySheet.totalStrand/2)
				middle = everySheet.strandList[middlePos]		#We find the "middle" of the even sheet
				#Example of a 7 Strand Sheet
				#    | |   | |
				#    | | | | |
				#| | | | |
				#| | | |
				listOfBulge.extend(goingLeftAndRight(middlePos, middle, everySheet))
			else: #The sheet is even, meaning there won't be a single bulge
				middlePos = len(everySheet.strandList)/2
				#Basically we want to know if we found the correct middle. Since this is even there are 2 possible "middles"
				# which means that we need to know which middle is the longest of the two. However, if they are the same
				# then the result doesn't matter. Only if the middle we have is shorter than the other middle.
				if (everySheet.strandList[middlePos].strandLength()) < (everySheet.strandList[middlePos-1].strandLength()):
					middle = everySheet.strandList[middlePos-1]
					middlePos -= 1
				else: middle = everySheet.strandList[middlePos]
				#Example of a 6 Strand Sheet
				#    | 
				#| | | | | |
				#| | | | |
				#    | |
				listOfBulge.extend(goingLeftAndRight(middlePos, middle, everySheet))
	return listOfBulge

def goingLeftAndRight(middlePos, middle, sheet):
	listOfBulge = []
	isBulge = True
	maximumStrands = sheet.totalStrand/2	#This is the maximum strands that are allowed to have maximum length
	for i in range(middlePos-1,-1,-1):
		#print i+1
		#print sheet.strandList[i].strandLength()
		#print middle.strandLength()
		middleLengthLeft = middle.strandLength()
		#Detects the moment when the strands are less than the middle strand length
		#if sheet.strandList[i].strandLength() < middleLengthLeft and middleLengthLeft == middle.strandLength():
		#	middleLengthLeft -= 1
		if sheet.strandList[i].strandLength() == middleLengthLeft and maximumStrands > 0: maximumStrands -=1
		elif maximumStrands <= 0: isBulge = False
		#print sheet.strandList[i].strandLength() <= middleLengthLeft
		if sheet.strandList[i].strandLength() <= middleLengthLeft and isBulge:
			middleLengthLeft = sheet.strandList[i].strandLength()
			isBulge = True
		else:
			isBulge = False
			break
	for i in range(middlePos+1,len(sheet.strandList)):
		#print i+1
		#print sheet.strandList[i].strandLength()
		#print middle.strandLength()
		middleLengthRight = middle.strandLength()
		#if sheet.strandList[i].strandLength() < middleLengthRight and middleLengthRight == middle.strandLength():
		#	middleLengthRight -= 1
		if sheet.strandList[i].strandLength() == middleLengthRight and maximumStrands > 0: maximumStrands -=1
		elif maximumStrands <= 0: isBulge = False
		#print sheet.strandList[i].strandLength() <= middleLengthRight
		if sheet.strandList[i].strandLength() <= middleLengthRight and isBulge: 
			middleLengthLeft = sheet.strandList[i].strandLength()
			isBulge = True
		else:
			isBulge = False
			break
	if isBulge: listOfBulge.append(sheet)
	return listOfBulge

#Algorithm currently has the problem of potentially detecting a square.
# Should consider using a length relative measurement for what constitutes as, "middle"
# This way, after said limit has been reached, the resultant strands of the same length, will
# no longer be considered the "middle" of the Beta-Sheet, therefore resulting in the overall
# Sheet to be void in the case of this definition.

#Half?
# 6 would be 3
# 7 would be 3... should it be 4?
# 8 would be 4
# 9 would be 4...
# 10 would be 5
# |     | | | |  
# | | | | | | | |  
# | | | | | | | | | |
#   | | | | | | | | |
#     | | | | | |  
# I think half is reasonable.