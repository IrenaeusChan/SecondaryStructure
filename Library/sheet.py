"""
Irenaeus Chan
11/27/2015

Sheets and Strands (that make up Sheets)
"""

from atom import Atom
from aminoacid import AminoAcid

class Strand(object):
	def __init__(self, strandNum, start, stop, seqres, strandType, amino_acids):
		"""Creates a new STRAND made up of Amino Acids

		Arguments:
			strandNum: The order of the STRAND in the SHEET
			start: The starting position of the STRAND
			stop: The ending position of the STRAND
			seqres: Which chain does the STRAND belong to
			amino_acids: All the Amino Acids that build up the STRAND

		Exceptions:
			ValuError: If given invalid start, stop, seqres, sheetType, or amino_acids
		"""

		if isinstance(strandNum, int):
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
			
		self.amino_acids = amino_acids

	def __hash__(self):
		return hash(self.__repr__())

	def __eq__(self,other):
		return self.amino_acids == other.amino_acids

	def __ne__(self, other):
		return not self.__eq__(other)

	def __repr__(self):
		strand_sequence = ""
		for AA in self.amino_acids:
			strand_sequence += "{0} {1}\n".format(AA.amino_acid, AA.position)
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

	def __hash__(self):
		return hash(self.__repr__())

	def __eq__(self,other):
		return self.strandList == other.strandList

	def __ne__(self, other):
		return not self.__eq__(other)

	def __repr__(self):
		sheet_sequence = "Sheet Identifier: {0}\nSeqRes: {1}\nTotal Strands: {2}\n".format(self.sheetIden, self.seqres, self.totalStrand)
		for strand in self.strandList:
			sheet_sequence += "{0}\n".format(strand)
		return sheet_sequence

def buildSheet(file_name, protein):
	sequence, strandList, sheetList = [], [], []
	with open(file_name, "r") as stream:
		for line in stream:
			#Reads from the file the information regarding the SHEET structures
			if (line[0:5] == "SHEET"):
				sheetIden = str(line[11:14])
				totalStrand = int(line[14:16])		#The number of Strands in the Sheet
				strandNum = int(line[8:10])
				start = int(line[22:26])
				stop = int(line[33:37])
				seqres = str(line[21:22])
				strandType = int(line[38:40])

				#Looks for the position of the Amino Acid using the already parsed sequence and copies the sequence
				for AA in protein.amino_acids:
					if (AA.position >= start and AA.position <= stop and AA.seqres == seqres):
						sequence.append(AA)
					if (AA.position == stop and AA.seqres == seqres):
						break

				#Appends the STRAND sequence to a list of other sequences
				strandList.append(Strand(strandNum, start, stop, seqres, strandType, sequence))
				sequence = []

				if len(strandList) == totalStrand:
					sheetList.append(Sheet(sheetIden, seqres, totalStrand, strandList))
					strandList = []
	return sheetList