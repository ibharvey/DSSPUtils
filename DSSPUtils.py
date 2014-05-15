#!/usr/bin/python

##=========================================================================================================
##						DSSPUtils: 
##
##				For use with Define Secondary Structure Prediction
##
##	    Specifically this file provides utilities to easily work with standard *.dssp files
##
##
##
##	Ian Harvey
##	Pappu Lab
##	May 2014
##=========================================================================================================


class DsspHBond:
	# argument is a string array of:
	# (1) The partner residue index
	# (2) The energy of the interaction
	def __init__(self,arrPE):
		if len(arrPE) == 2:
			self.partner = int(arrPE[0])
			self.energy = float(arrPE[1])
		else:
			raise ValueError("The argument must list the partner and energy")
		
class DsspResidue:
	def __init__(self, line):
		if len(line) == 137:
			self.dsspIndex = int(line[:5])
			self.pdbIndex = int(line[6:10])
			self.chainId = line[11]
			self.resId = line[13]
			self.structure = line[15:25]
			self.bp1 = int(line[26:29])
			self.bp2 = int(line[30:33])
			self.bp2C = line[33]
			self.sas = int(line[34:38])
			self.nh_o_1 = DsspHBond(line[40:50].split(','))
			self.o_nh_1 = DsspHBond(line[51:61].split(','))
			self.nh_o_2 = DsspHBond(line[62:72].split(','))
			self.o_nh_2 = DsspHBond(line[73:83].split(','))
			self.tco = float(line[85:91])
			self.kappa = float(line[91:97])
			self.alpha = float(line[97:103])
			self.phi = float(line[103:109])
			self.psi = float(line[109:115])
			self.xca = float(line[116:122])
			self.yca = float(line[122:129])
			self.zca = float(line[129:])
		else:
			raise ValueError("The input DSSP residue format was not 137 characters, can't parse.")
		
class DsspChain:
	def __init__(self, lines):
		self.residues = []
		for index in range(0,len(lines)):
			aRes = DsspResidue(lines[index])
			self.residues.append(aRes)

class DsspStructure:
	def __init__(self, DsspFile):
		# Open the DSSP file to be parsed
		oneDSSP = open(DsspFile,'r').readlines()
		
		# Find index of the start of DSSP per residue data
		startIndex = next(i for i in range(0,len(oneDSSP)) if '#' in oneDSSP[i]) + 1
		
		# DSSP uses '!*' to demarcate chains
		endChainIndeces = [ i for i in range(0,len(oneDSSP)) if '!*' in oneDSSP[i] ]
		
		# Array of chains
		self.chains = []
		
		# Iterate through chains to get DSSP info per residue
		for i in range(0,len(endChainIndeces)+1):
			# if i is the length of the endChainIndeces, the next end is the true endChainIndeces
			if i == len(endChainIndeces):
				self.chains.append(DsspChain(oneDSSP[startIndex:]))
			else:
				self.chains.append(DsspChain(oneDSSP[startIndex:endChainIndeces[i]]))
				startIndex = endChainIndeces[i] + 1
			setattr(self, 'Chain'+chr(65+i), self.chains[i])
