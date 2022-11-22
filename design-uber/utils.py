"""
Utilities for uber array design
"""

def ReverseComplement(seq):
	"""
	Return the reverse complement of a string

	Arguments
	---------
	seq : str
	   Input sequence

	Returns
	-------
	revcomp : str
	   Reverse complement of the input sequence

	Example
	-------
	ReverseComplement("ATTC")
	> "GAAT"
	"""
	return "NNNNN" # Susan!

def Canonicalize(seq):
	"""
	Return the canonical version of a repeat unit

	Arguments
	---------
	seq : str
	   Input repeat unit sequence

	Returns
	-------
	canon_seq : str
	   Canonicalized version of the repeat unit

	Example
	-------
	Canonicalize("GTT")
	> "AAC"
	"""
	return seq # TODO !!! Not implemented Susan will do this

if __name__ == "__main__":
	print("Testing utils functions")
	print(ReverseComplement("AG"))