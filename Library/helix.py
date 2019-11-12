"""
Irenaeus Chan
12/13/2016

Helix Class
Used for the BINF6210 Final Project
"""

#Configuration for a Helix Object
class Helix(object):
	def __init__(self, helixIden, start, stop, seqres, helixType, aminoAcidList):
		"""Creates a new Helix made up of Amino Acids

		Arguments:
			start: The starting position of the Helix
			stop: The ending position of the Helix
			seqres: Which chain does the Helix belong to
			aminoAcidList: All the Amino Acids that build up the Helix

		Exceptions:
			ValuError: If given any invalid parameters
		"""
		if isinstance(helixIden, basestring): self.helixIden = helixIden
		else: raise ValueError('Invalid SEQRES {0}'.format(helixIden))

		if isinstance(start, int): self.start = start
		else: raise ValueError('Invalid Start Position {0}'.format(start))

		if isinstance(stop, int): self.stop = stop
		else: raise ValueError('Invalid Stop Position {0}'.format(stop))

		if isinstance(seqres, basestring): self.seqres = seqres
		else: raise ValueError('Invalid SEQRES {0}'.format(seqres))

		if isinstance(helixType, int): self.helixType = helixType
		else: raise ValueError('Invalid Type {0}'.format(helixType))
			
		self.aminoAcidList = aminoAcidList

	def __eq__(self, other): return self.__dict__ == other.__dict__
	def __ne__(self, other): return not self.__eq__(other)
	def __repr__(self):
		helixSequence = ""
		aarepr = ""
		for aa in self.amino_acids:
			if aa.position != self.aminoAcidList[-1].position:
				aarepr = str(aa.position)
				if len(aarepr) == 2: aarepr = "0" + aarepr
				elif len(aarepr) == 1: aarepr = "00" + aarepr
				helixSequence += aarepr + " --- "
			else: helixSequence += aarepr
		return helixSequence

def buildHelix(filename, protein):
	helixAminoAcidList, helixList = [], []

	with open(filename, "r") as stream:
		for line in stream:
			#Reads from the file the information regarding the HELIX structures
			if (line[0:5] == "HELIX"):
				helixIden = str(line[11:14])		#The PDB Identifier Number
				start = int(line[21:25])
				stop = int(line[33:37])
				seqres = str(line[19:20])
				helixType = int(line[39:40])
				#Looks for the position of the Amino Acid using the already parsed sequence and copies the sequence
				for aa in protein.aminoAcidList:
					if (aa.position >= start and aa.position <= stop and aa.seqres == seqres): helixAminoAcidList.append(aa)
					if (aa.position == stop and aa.seqres == seqres): break

				#Appends the HELIX sequence to a list of other sequences
				helixList.append(Helix(helixIden, start, stop, seqres, helixType, helixAminoAcidList))
				helixAminoAcidList = []
	return helixList