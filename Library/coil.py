"""
Irenaeus Chan
12/13/2016

Coil Class
Used for the BINF6210 Final Project
"""

#Configuration for a Coil Object
class Coil(object):
	def __init__(self, start, stop, aminoAcidList):
		"""Creates a new Coil made up of Amino Acids

		Arguments:
			start: The starting position of the Coil
			stop: The ending position of the Coil
			seqres: Which chain does the Coil belong to
			amino_acids: All the Amino Acids that build up the Coil

		Exceptions:
			ValuError: If given any invalid parameters
		"""

		if isinstance(start, int): self.start = start
		else: raise ValueError('Invalid Start Position {0}'.format(start))

		if isinstance(stop, int): self.stop = stop
		else: raise ValueError('Invalid Stop Position {0}'.format(stop))
	
		self.aminoAcidList = aminoAcidList

	def __eq__(self, other): return self.__dict__ == other.__dict__
	def __ne__(self, other): return not self.__eq__(other)
	def __repr__(self):
		coilSequence = ""
		aarepr = ""
		for aa in self.aminoAcidList:
			if aa.position != self.aminoAcidList[-1].position:
				aarepr = str(aa.position)
				if len(aarepr) == 2: aarepr = "0" + aarepr
				elif len(aarepr) == 1: aarepr = "00" + aarepr
				coilSequence += aarepr + " --- "
			else: coilSequence += aarepr
		return coilSequence

def buildCoil(filename, protein):
	coilAminoAcidList, coilList, startList, stopList = [], [], [0], []

	with open(filename, "r") as stream:
		for line in stream:
			#Reads from the file the information regarding the HELIX structures
			if (line[0:5] == "HELIX"):
				start = int(line[21:25])
				stop = int(line[33:37])
				stopList.append(start)
				startList.append(stop)
			
			#Reads from the file the information regarding the SHEET structures
			if (line[0:5] == "SHEET"):
				start = int(line[22:26])
				stop = int(line[33:37])
				stopList.append(start)
				startList.append(stop)
				#There's a problem here where the beta turns are not "coils"

		stopList.append(protein.aminoAcidList[-1].position)

		#There is a logic problem here with how it partitions the coils that won't work
		# properly with the Phi and Psi calculator
		for start, stop in zip(startList, stopList):
			for aa in protein.aminoAcidList:
				if (aa.position > start and aa.position < stop): 
					coilAminoAcidList.append(aa)
			tcoil = Coil(start, stop, coilAminoAcidList)
			print tcoil
			coilList.append(tcoil)
			coilAminoAcidList = []
	return coilList