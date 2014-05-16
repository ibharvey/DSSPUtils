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
			# If the residue is a lowercase letter, it's a cysteine
			self.resId = line[13] if not line[13].islower() else 'C'
			self.ss = line[16]
			# parse this better
			self.structure = line[17:25]
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
		self.name = self.residues[0].chainId

class DsspStructure:
	def __init__(self, DsspFile):
		# Open the DSSP file to be parsed
		oneDSSP = open(DsspFile,'r').readlines()
		
		# Find index of the start of DSSP per residue data
		header =   '#  RESIDUE AA STRUCTURE BP1 BP2  ACC     N-H-->O    O-->H-N    N-H-->O    O-->H-N    TCO  KAPPA ALPHA  PHI   PSI    X-CA   Y-CA   Z-CA'
		startIndex = next(i for i in range(0,len(oneDSSP)) if header in oneDSSP[i]) + 1
		
		# DSSP uses '!*' to demarcate chains
		# DSSP uses '!' to signify breaks within chains
		endChainIndeces = [ i for i in range(0,len(oneDSSP)) if '!' in oneDSSP[i] ]
		
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
			
			# Setting individual chain variables
			
			# If this chain type is already split into atleast 2, find next number to use
			if 'Chain'+self.chains[i].name + '1' in self.__dict__.keys():
				iter = 2
				while 'Chain'+self.chains[i].name + str(iter) in self.__dict__.keys():
					iter += 1
				setattr(self, 'Chain'+self.chains[i].name + str(iter), self.chains[i])
			
			# If this is the second occurence of this chain type, split the variable
			elif 'Chain' + self.chains[i].name in self.__dict__.keys():
				# Modify the variable named 'ChainX' to 'ChainX1'
				setattr(self, 'Chain'+self.chains[i].name + '1', getattr(self, 'Chain' + self.chains[i].name))
				delattr(self, 'Chain'+self.chains[i].name)
				
				# Add the new variable as 'ChainX2'
				setattr(self, 'Chain'+self.chains[i].name + '2', self.chains[i])
			
			else:
				# Don't assume there is a break in a chain
				setattr(self, 'Chain'+self.chains[i].name, self.chains[i])



