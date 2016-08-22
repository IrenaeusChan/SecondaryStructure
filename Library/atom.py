"""
Irenaeus Chan
11/27/2015

Atom Class
"""

BACKBONE_ATOMS = {'N', 'Ca', 'C'}
ELEMENTS = {'N', 'C', 'O', 'S', 'H', 'P'}

import sys

class Atom (object):
	"""A configuration for a single Backbone Atom"""

	def __init__(self, atom, x, y, z, position, element):
		"""Creates a new Atom

		Arguments:
			atom: The specific Backbone Atom
			x: X position of the atom
			y: Y position of the atom
			z: Z position of the atom
			position: The Amino Acid which the Atom belongs to
			element: Which element the atom is made of (mainly important for side chains)

		Exceptions:
			ValuError: If given invalid atom, x, y, or z
		"""

		#if atom in BACKBONE_ATOMS:
		self.atom = atom
		#else:
			#raise ValueError('Invalid Atom {0}'.format(atom))

		if isinstance(x, float):
			self.x = x
		else:
			raise ValueError('Invalid X {0}'.format(x))

		if isinstance(y, float):
			self.y = y
		else:
			raise ValueError('Invalid Y {0}'.format(y))

		if isinstance(z, float):
			self.z = z
		else:
			raise ValueError('Invalid Z {0}'.format(z))

		if isinstance(position, int):
			self.position = position
		else:
			raise ValueError('Invalid Position {0}'.format(position))

		if element in ELEMENTS:
			self.element = element
		else:
			raise ValueError('Invalid Element {0}'.format(element))

	def __hash__(self):
		return hash(self.__repr__())

	def __eq__(self, other):
		return self.__dict__ == other.__dict__

	def __ne__(self, other):
		return not self.__eq__(other)

	def __repr__(self):
		return "\nAtom: {0} at ({1}, {2}, {3}) and {4}".format(self.atom, self.x, self.y, self.z, self.element)