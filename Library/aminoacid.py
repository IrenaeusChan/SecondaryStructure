"""
Irenaeus Chan
11/27/2015

Amino Acid Class w/ Functions
"""

from atom import Atom

AMINO_ACIDS = {'GLY', 'ALA', 'SER', 'THR', 'CYS', 'VAL', 'LEU', 'ILE', 'MET', 'PRO', 'PHE', 'TYR', 'TRP', 'ASP', 'GLU', 'ASN', 'GLN', 'HIS', 'LYS', 'ARG', 'TER'}
ELEMENTS = {'N':14, 'C':12, 'O':16, 'S':32, 'H':1, 'P':31}

class AminoAcid (object):
	"""A configuration for a single Amino Acid"""

	def __init__(self, amino_acid, seqres, position, backbone, sidechain):
		"""Creates a new Amino Acid

		Arguments:
			amino_acid: The specific Amino Acid
			seqres: Which chain does the Amino Acid belong to
			position: The exact position of the Amino Acid
			backbone: The backbone atoms
			sidechain: The sidechain atoms

		Exceptions:
			ValuError: If given invalid amino_acid, seqres, position, or backbone
		"""

		if amino_acid in AMINO_ACIDS:
			self.amino_acid = amino_acid
		else:
			raise ValueError('Invalid Amino Acid {0}'.format(amino_acid))

		if isinstance(seqres, basestring):
			self.seqres = seqres
		else:
			raise ValueError('Invalid SEQRES {0}'.format(seqres))

		if isinstance(position, int):
			self.position = position
		else:
			raise ValueError('Invalid Position {0}'.format(position))

		self.backbone = backbone
		self.sidechain = sidechain

		self.avgx, self.avgy, self.avgz = weightedAverage(self)

	def __hash__(self):
		return hash(self.__repr__())

	def __eq__(self, other):
		return self.__dict__ == other.__dict__

	def __ne__(self, other):
		return not self.__eq__(other)

	def __repr__(self):
		AA = "Amino Acid: {0}\nSEQRES: {1}\nPosition: {2}\n".format(self.amino_acid, self.seqres, self.position)
		for a in self.backbone:
			AA += "{0}".format(a)
		for b in self.sidechain:
			AA += "{0}".format(b)
		AA += "\nWeighted Means: {0}, {1}, {2}\n".format(self.avgx, self.avgy, self.avgz)
		return AA

def weightedAverage(self):
	totalX = 0
	totalY = 0
	totalZ = 0
	totalMass = 0

	for atom in self.backbone:
		totalX += atom.x * ELEMENTS[atom.element]
		totalY += atom.y * ELEMENTS[atom.element]
		totalZ += atom.z * ELEMENTS[atom.element]
		totalMass += ELEMENTS[atom.element]
	for atom in self.sidechain:
		totalX += atom.x * ELEMENTS[atom.element]
		totalY += atom.y * ELEMENTS[atom.element]
		totalZ += atom.z * ELEMENTS[atom.element]
		totalMass += ELEMENTS[atom.element]

	totalX = totalX/totalMass
	totalY = totalY/totalMass
	totalZ = totalZ/totalMass

	return totalX, totalY, totalZ