def detectSmall(listOfSheets):
	listOfSmall = []
	isSmall = False
	for everySheet in listOfSheets:
		if everySheet.totalStrand > 4 and everySheet.totalStrand <= 10:
			# for everyStrand in everySheet.strandList:
			# 	if everyStrand.strandLength() <= 4:
			# 		isSmall = True
			# 	else:
			# 		isSmall = False
			# 		break
			# if isSmall == True:
			listOfSmall.append(everySheet)
			isSmall = False

	return listOfSmall