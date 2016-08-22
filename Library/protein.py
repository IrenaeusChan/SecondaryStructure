"""
Irenaeus Chan
11/27/2015

Protein Class w/ Functions
"""

backboneAtoms = (" N  ", " CA ", " C  ")
AMINO_ACIDS = {'GLY':1, 'ALA':2, 'SER':3, 'THR':4, 'CYS':3, 'VAL':4, 'LEU':5, 'ILE':5, 'MET':5, 'PRO':4, 'PHE':8, 'TYR':9, 'TRP':11, 'ASP':5, 'GLU':6, 'ASN':5, 'GLN':6, 'HIS':7, 'LYS':6, 'ARG':8, 'TER':0}
ELEMENTS = {'N':14, 'C':12, 'O':16, 'S':32, 'H':1, 'P':31}

import sys
import vector
from atom import Atom
from aminoacid import AminoAcid

class Protein(object):
	def __init__(self, amino_acids):
		"""Creates a Protein built of Amino Acids

		Arguments:
			amino_acids: A list of Amino Acids
		"""
		self.amino_acids = amino_acids
		self.length = amino_acids[-1].position

	def __hash__(self):
		return hash(self.__repr__())

	def __eq__(self,other):
		return self.amino_acids == other.amino_acids

	def __ne__(self, other):
		return not self.__eq__(other)

	def __repr__(self):
		protein_sequence = ""
		for AA in self.amino_acids:
			protein_sequence += "{0}\n".format(AA)
		return protein_sequence

def buildProtein(file_name):
	backbone = []
	sidechain = []
	protein = []
	currentAminoAcid = ""
	currentSeq = ""
	currentPos = 0

	with open(file_name, "r") as stream:
		for line in stream:
			if(line[0:4] == "ATOM"):
				if ((currentPos != int(line[22:26])) and currentPos != 0):
					#Amino Acid, SEQRES, Position, Backbone[N][Ca][C], SideChain[1]...[*]
					aminoacid = AminoAcid(currentAminoAcid, currentSeq, currentPos, list(backbone), list(sidechain))
					#print aminoacid
					protein.append(aminoacid)
					backbone = []
					sidechain = []
				currentAminoAcid = str(line[17:20])
				currentSeq = str(line[21:22])
				currentPos = int(line[22:26])
				if (line[12:16] in backboneAtoms):
					backbone.append(Atom(line[12:16].replace(" ", ""),float(line[31:38]), float(line[39:46]), float(line[47:54]), currentPos, str(line[77:78])))
				else:
					sidechain.append(Atom(line[12:16].replace(" ", ""), float(line[31:38]), float(line[39:46]), float(line[47:54]), currentPos, str(line[77:78])))
		aminoacid = AminoAcid(currentAminoAcid, currentSeq, currentPos, list(backbone), list(sidechain))
		protein.append(aminoacid)
	return Protein(list(protein))

def buildProteinOriginal(file_name): #This is the proper one, do not change
	backbone = []
	sidechain = []
	protein = []

	with open(file_name, "r") as stream:
		for line in stream:
			if(line[0:4] == "ATOM") and line[12:16] in backboneAtoms:			#or line[0:6] == "HETATM	<- Not sure what to do with this"
				#Once we know which Amino Acid we are reading, we can instantiate the number of atoms in the sidechain
				sidechain_count = AMINO_ACIDS[line[17:20]]
				backbone_count = 3
				#Creates the Atom: Atom, X, Y, Z, Element
				if (backbone_count > 0):
					backbone.append(Atom(line[12:16].replace(" ", ""),float(line[31:38]), float(line[39:46]), float(line[47:54]), line[77:78]))
				backbone_count-=1
			elif(line[0:4] == "ATOM") and line[77:78] != 'H': 	#For Side Chains...
			#I am deliberately ignoring H here because the mass of Hydrogen is very insignificant, also some PDBs have it, some don't
				if (sidechain_count > 0):
					sidechain.append(Atom(line[12:16].replace(" ", ""), float(line[31:38]), float(line[39:46]), float(line[47:54]), line[77:78]))
				sidechain_count-=1

				#We put this here because the PDB file has the backbone atoms first, then the side chain, so once side chains are
				# finished, then it implies that all the information required for the Amino Acid has been stored
				if (sidechain_count == 0):
					#Amino Acid, SEQRES, Position, Backbone[N][Ca][C], SideChain[1]...[*]
					aminoacid = AminoAcid(line[17:20], line[21:22], int(line[22:26]), list(backbone), list(sidechain))
					#print aminoacid
					protein.append(aminoacid)
					backbone = []
					sidechain = []
	
	return Protein(list(protein))

def weightedAverage(protein):
	totalX = 0
	totalY = 0
	totalZ = 0
	totalMass = 0

	for AA in protein.amino_acids:
		for atom in AA.backbone:
			totalX += atom.x * ELEMENTS[atom.element]
			totalY += atom.y * ELEMENTS[atom.element]
			totalZ += atom.z * ELEMENTS[atom.element]
			totalMass += ELEMENTS[atom.element]
		for atom in AA.sidechain:
			totalX += atom.x * ELEMENTS[atom.element]
			totalY += atom.y * ELEMENTS[atom.element]
			totalZ += atom.z * ELEMENTS[atom.element]
			totalMass += ELEMENTS[atom.element]

	totalX = totalX/totalMass
	totalY = totalY/totalMass
	totalZ = totalZ/totalMass

	return totalX, totalY, totalZ

def relativeToCenter(protein, center):
	write = 'w'
	if (sys.argv[1] == "all" and len(sys.argv) > 2):
		write = 'a'
	with open("distances.txt", write) as output:
		for AA in protein.amino_acids:
			aminoacid = [AA.avgx, AA.avgy, AA.avgz]
			extraInfo = AA.amino_acid
			d = vector.vectorMagnitude(vector.vectorCalculation(center, aminoacid))
			output.write(extraInfo + ' ' + str(d) + '\n')