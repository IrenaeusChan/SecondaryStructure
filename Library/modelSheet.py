"""
Irenaeus Chan
12/01/2016

Sheets and Strands (that make up Sheets)
Not the same as the other Sheet.py, this one is less detailed
"""

from itertools import tee, islice, chain, izip
import re
import string

class Strand(object):
	def __init__(self, strandNum, start, stop, seqres, strandType, aminoAcidNumbers, thisStrand, otherStrand):
		"""Creates a new STRAND made up of Amino Acids

		Arguments:
			strandNum: This will be a LETTER representing which SHEET it belongs to and a NUMBER 
						which will correspond to the order of the STRAND in the SHEET
			start: The starting position of the STRAND
			stop: The ending position of the STRAND
			seqres: Which chain does the STRAND belong to
			strandType: The number which determines Parallel (1) or Anti-Parallel (-1) or Starting (0)
			aminoAcidNumbers: The NUMERICAL positions of the AMINO ACIDS that build up the STRAND

			thisStrand: RESIDUE NUMBER of THIS STRAND that connects to otherStrand
			otherStrand: RESIDUE NUMBER of OTHER STRAND that connects to thisStrand

		Exceptions:
			ValuError: If given invalid start, stop, seqres, sheetType, or amino_acids
		"""

		if isinstance(strandNum, basestring):
			self.strandNum = strandNum
		else:
			raise ValueError('Invalid SEQRES {0}'.format(strandNum))

		if isinstance(start, int):
			self.start = start
		else: 
			raise ValueError('Invalid Start Position {0}'.format(start))

		if isinstance(stop, int):
			self.stop = stop
		else:
			raise ValueError('Invalid Stop Position {0}'.format(stop))

		if isinstance(seqres, basestring):
			self.seqres = seqres
		else:
			raise ValueError('Invalid SEQRES {0}'.format(seqres))

		if isinstance(strandType, int):
			self.strandType = strandType
		else:
			raise ValueError('Invalid Type {0}'.format(strandType))
			
		self.aminoAcidNumbers = aminoAcidNumbers

		if isinstance(thisStrand, int):
			self.thisStrand = thisStrand
		elif thisStrand == None:
			self.thisStrand = None
		else:
			raise ValueError('Invalid Type {0}'.format(otherStrand))

		if isinstance(otherStrand, int):
			self.otherStrand = otherStrand
		elif otherStrand == None:
			self.otherStrand = None
		else:
			raise ValueError('Invalid Type {0}'.format(otherStrand))

	def __hash__(self):
		return hash(self.__repr__())

	def __eq__(self,other):
		return self.aminoAcidNumbers == other.aminoAcidNumbers

	def __ne__(self, other):
		return not self.__eq__(other)

	"""
	Assumption: We must assume the FIRST STRAND (0) in the SHEET goes from LEFT to RIGHT
				From here, ANY successive "Parallel (1)" STRAND will also go LEFT to RIGHT
				AND ANY successive "Anti-Parallel (-1)" STRAND will therefore have to go LEFT to RIGHT
	"""
	def __repr__(self):
		strand_sequence = ""
		for number in self.aminoAcidNumbers:
			if number != self.aminoAcidNumbers[-1]:
				strand_sequence += "{0} --- ".format(number)
			else:
				strand_sequence += "{0}".format(number)
		return strand_sequence

class Sheet(object):
	def __init__(self, sheetIden, seqres, totalStrand, strandList):
		"""Creates a new Sheet made up of Amino Acids

		Arguments:
			sheetIden: The identification character of the SHEET
			seqres: Which chain does the SHEET belong to
			totalStrand: The total number of STRANDS that make up the SHEET
			strandList: The list of STRANDS that make up the SHEET

		Exceptions:
			ValuError: If given invalid start, stop, seqres, sheetType, or amino_acids
		"""

		if isinstance(sheetIden, basestring):
			self.sheetIden = sheetIden
		else:
			raise ValueError('Invalid SEQRES {0}'.format(sheetIden))

		if isinstance(seqres, basestring):
			self.seqres = seqres
		else:
			raise ValueError('Invalid SEQRES {0}'.format(seqres))

		if isinstance(totalStrand, int):
			self.totalStrand = totalStrand
		else:
			raise ValueError('Invalid SEQRES {0}'.format(totalStrand))
			
		self.strandList = strandList

	def breakApart(self):
		broken = []
		for strand in self.strandList:
			if strand.strandType == 0:
				beforeEditStrand = strand.strandNum + " " + "{0}".format(strand).replace("---", "-->")
			elif re.search('-->', beforeEditStrand) and strand.strandType != 0:
				if strand.strandType == 1: #-->
					beforeEditStrand = strand.strandNum + " " + "{0}".format(strand).replace("---", "-->")
				elif strand.strandType == -1: #<--
					beforeEditStrand = strand.strandNum + " " + reverseDirection("{0}".format(strand).replace("---", "<--"))
			elif re.search('<--', beforeEditStrand) and strand.strandType != 0:
				if strand.strandType == 1:
					beforeEditStrand = strand.strandNum + " " + reverseDirection("{0}".format(strand).replace("---", "<--"))
				elif strand.strandType == -1:
					beforeEditStrand = strand.strandNum + " " + "{0}".format(strand).replace("---", "-->")
			broken.append(beforeEditStrand.split())
		return broken

	def addCase2and3(self, otherSheet):
		for selfStrand in self.strandList:
			for otherStrand in otherSheet.strandList:
				if ((selfStrand.start >= otherStrand.start and selfStrand.start <= otherStrand.stop) 
						or (selfStrand.stop >= otherStrand.start and selfStrand.stop <= otherStrand.stop)):
					currentStrand.strandNum
					otherStrand.strandNum
				elif ((otherStrand.start >= selfStrand.start and otherStrand.start <= selfStrand.stop) 
						or (otherStrand.stop >= selfStrand.start and otherStrand.stop <= selfStrand.stop)):
					currentStrand.strandNum
					otherStrand.strandNum


	def __hash__(self):
		return hash(self.__repr__())

	def __eq__(self,other):
		return self.strandList == other.strandList

	def __ne__(self, other):
		return not self.__eq__(other)

	def __repr__(self):
		#sheetRepresent = "Sheet Identifier: {0}, SeqRes: {1}, Total Strands: {2}\n".format(self.sheetIden.split()[0], self.seqres, self.totalStrand)
		sheetRepresent = ""
		listOfStrands, templistOfStrands = [], []
		nextStrandSpacing = 0
		for prev, currentStrand, nxt in previousAndNext(self.strandList):
			if currentStrand.strandType == 0:
				beforeEditStrand = "{0}\n".format(currentStrand).replace("---", "-->")
				
				if nxt.strandType == 1:
					differenceValue = abs(nxt.thisStrand - nxt.start) - abs(nxt.otherStrand - currentStrand.start)
				elif nxt.strandType == -1:
					differenceValue = abs(nxt.thisStrand - nxt.stop) - abs(nxt.otherStrand - currentStrand.start)
				
				if differenceValue > 0:
					beforeEditStrand = " "*8*differenceValue + beforeEditStrand
				elif differenceValue < 0:
					nextStrandSpacing += abs(differenceValue)
			elif re.search('-->', beforeEditStrand) and currentStrand.strandType != 0:
				if currentStrand.strandType == 1: #-->
					beforeEditStrand = " "*8*nextStrandSpacing + "{0}".format(currentStrand).replace("---", "-->") + "\n"
					if nxt is not None:
						if nxt.strandType == 1:	#Next Strand is --> Meaning Start
							differenceValue = abs(nxt.thisStrand - nxt.start) - abs(nxt.otherStrand - currentStrand.start)
						elif nxt.strandType == -1:	#Next Strand is <-- Meaning Stop
							differenceValue = abs(nxt.thisStrand - nxt.stop) - abs(nxt.otherStrand - currentStrand.start)
						
						if differenceValue > 0:
							beforeEditStrand = " "*8*differenceValue + beforeEditStrand
							for strandString in listOfStrands:
								templistOfStrands.append([strandString[0], " "*8*differenceValue + strandString[1]])
							listOfStrands = []
							listOfStrands = templistOfStrands
						elif differenceValue < 0:
							nextStrandSpacing += abs(differenceValue)
				elif currentStrand.strandType == -1: #<--
					beforeEditStrand = " "*8*nextStrandSpacing + reverseDirection("{0}".format(currentStrand).replace("---", "<--")) + "\n"
					if nxt is not None:
						if nxt.strandType == 1:	#Next Strand is <-- Meaning Stop
							differenceValue = abs(nxt.thisStrand - nxt.stop) - abs(nxt.otherStrand - currentStrand.stop)
						elif nxt.strandType == -1:	#Next Strand is --> Meaning Start
							differenceValue = abs(nxt.thisStrand - nxt.start) - abs(nxt.otherStrand - currentStrand.stop)
						
						if differenceValue > 0:
							beforeEditStrand = " "*8*differenceValue + beforeEditStrand
							for strandString in listOfStrands:
								templistOfStrands.append([strandString[0], " "*8*differenceValue + strandString[1]])
							listOfStrands = []
							listOfStrands = templistOfStrands
						elif differenceValue < 0:
							nextStrandSpacing += abs(differenceValue)
			elif re.search('<--', beforeEditStrand) and currentStrand.strandType != 0:
				if currentStrand.strandType == 1:
					beforeEditStrand = " "*8*nextStrandSpacing + reverseDirection("{0}".format(currentStrand).replace("---", "<--")) + "\n"
					if nxt is not None:
						if nxt.strandType == 1:
							differenceValue = abs(nxt.thisStrand - nxt.stop) - abs(nxt.otherStrand - currentStrand.stop)
						elif nxt.strandType == -1:
							differenceValue = abs(nxt.thisStrand - nxt.start) - abs(nxt.otherStrand - currentStrand.stop)
						
						if differenceValue > 0:
							beforeEditStrand = " "*8*differenceValue + beforeEditStrand
							for strandString in listOfStrands:
								templistOfStrands.append([strandString[0], " "*8*differenceValue + strandString[1]])
							listOfStrands = []
							listOfStrands = templistOfStrands
						elif differenceValue < 0:
							nextStrandSpacing += abs(differenceValue)
				elif currentStrand.strandType == -1:
					beforeEditStrand = " "*8*nextStrandSpacing + "{0}".format(currentStrand).replace("---", "-->") + "\n"
					if nxt is not None:
						if nxt.strandType == 1:
							differenceValue = abs(nxt.thisStrand - nxt.start) - abs(nxt.otherStrand - currentStrand.start)
						elif nxt.strandType == -1:
							differenceValue = abs(nxt.thisStrand - nxt.stop) - abs(nxt.otherStrand - currentStrand.start)
						
						if differenceValue > 0:
							beforeEditStrand = " "*8*differenceValue + beforeEditStrand
							for strandString in listOfStrands:
								templistOfStrands.append([strandString[0], " "*8*differenceValue + strandString[1]])
							listOfStrands = []
							listOfStrands = templistOfStrands
						elif differenceValue < 0:
							nextStrandSpacing += abs(differenceValue)
			listOfStrands.append([currentStrand.strandNum, beforeEditStrand])
			templistOfStrands = []
		sheetRepresent += ''.join(word[0] + " " + word[1] for word in listOfStrands)
		return sheetRepresent

def previousAndNext(some_iterable):
	#http://stackoverflow.com/questions/1011938/python-previous-and-next-values-inside-a-loop
	prevs, items, nexts = tee(some_iterable, 3)
	prevs = chain([None], prevs)
	nexts = chain(islice(nexts, 1, None), [None])
	return izip(prevs, items, nexts)

def reverseDirection(stringOfNumbers):
	toFlip = stringOfNumbers.split(" ")
	return ' '.join(word for word in toFlip[::-1])

def buildSheet(file_name):
	aminoAcidNumbers, strandList, sheetList = [], [], []
	thisStrand, otherStrand = None, None
	with open(file_name, "r") as stream:
		for line in stream:
			#Reads from the file the information regarding the SHEET structures
			if (line[0:5] == "SHEET"):
				sheetIden = str(line[11:14])
				totalStrand = int(line[14:16])		#The number of Strands in the Sheet
				strandNum = str(line[8:10])
				start = int(line[22:26])
				stop = int(line[33:37])
				seqres = str(line[21:22])
				strandType = int(line[38:40])
				if strandType != 0:
					thisStrand = int(line[51:55])
					otherStrand = int(line[66:69])

				for i in range(start, stop+1):
					if len(str(i)) == 2:
						aminoAcidNumbers.append("0" + str(i))
					elif len(str(i)) == 1:
						aminoAcidNumbers.append("00" + str(i))
					else:
						aminoAcidNumbers.append(str(i))

				#Appends the STRAND sequence to a list of other sequences
				strandList.append(Strand(sheetIden.strip() + strandNum.strip(), start, stop, seqres, +
					+ strandType, aminoAcidNumbers, thisStrand, otherStrand))
				aminoAcidNumbers = []

				if len(strandList) == totalStrand:
					sheetList.append(Sheet(sheetIden, seqres, totalStrand, strandList))
					strandList = []
	return sheetList