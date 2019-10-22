import numpy as np
import warnings

warnings.simplefilter(action='ignore', category=FutureWarning)


class DescendantProgenitor:
	'A wrapper/class (structure-like) to save a few properties of the main branch halos. Useful for dictionaries.'

	progID = '123123123123123123'
	descID = '123123123123123123'
	progNPart = 0
	descNPart = 0
	descNProg = 0
	
	def __init__(self, progID, progNPart, descID, descNPart, descNProg):
		self.progID = progID
		self.descID = descID
		self.progNPart = progNPart
		self.descNPart = descNPart
		self.descNProg = descNProg


class MTree:
        'Basic Merger tree class for the database'

        # How many snapshots is our merger tree composed of (upper limit)
        nSteps = 0

        # Main halo ID at z = 0. We use strings as we will be using dictionaries all along to connect these objects
        mainID = '123123123123123123123'

        # Time, expansion factor and redshift
        t = np.zeros((0))
        a = np.zeros((0))
        z = np.zeros((0))

        # Main branch halo: number of particles, main halo ID and number of descendants
        mainBranchID = []
        mainBranchNPart = np.full((0), 0)
        mainBranchNProg = np.full((0), 0)

        # Keep track of the steps at which the halo was lost and replaced by a token instead
        isToken = np.full((0), False)

        # Initialize the Merger Tree by its (expected) number of steps
        def __init__(self, n, ID):
                self.mainID = ID
                self.nSteps = n
                self.z = np.zeros((n))
                self.a = np.zeros((n))
                self.t = np.zeros((n))
                self.mainBranchID = [''] * n
                self.mainBranchNPart = np.full((n), 0)
                self.mainBranchNProg = np.full((n), 0)
                self.isToken = np.full((n), False)

       	def fill_mass_id(self, nps, ids):
                for iM in range(0 , self.nSteps):
                        self.mainBranchNPart[iM] = nps[iM]
                        self.mainBranchID[iM] = ids[iM]

        def get_mass_id(self):
                return [self.mainBranchNPart, self.mainBranchID]

        # Print the number of particles and mass ID corresponding 
        def print_mass_id(self):
                for iM in range(0, self.nSteps):
                        print("%s %d" % (self.mainBranchID[iM], self.mainBranchNPart[iM]))

        # Return the normalized (z=0) mass across all history
        def norm_mass(self):
                norm_mass = []
                m_zero = self.mainBranchNPart[0]

                return self.mainBranchNPart[:]/m_zero

        # Add the descendant progenitor pairs to the next step
        def update_desc_prog(self, iStep, dP):
                self.mainBranchID[iStep] = dP.descID
                self.mainBranchNPart[iStep] = dP.descNPart
                self.mainBranchNProg[iStep] = dP.descNProg

        # Read & implement the steps from an expansion factor file
        def init_a(self, file_name):
                self.a = read_a(file_name)

        # Read & implement the steps from a redshift file
        def init_z(self, file_name):
                self.z = read_z(file_name)

        # Print the number of particles and mass ID corresponding 
        def dump_to_mass_id(self):
                f_name = 'halo_' + self.mainID + '.full_tree'
                f_out = open(f_name, 'w')

                for iM in range(0, self.nSteps):
                        line = "%s %d" % (self.mainBranchID[iM], self.mainBranchNPart[iM]) + '\n'
                        #print(line, file=f_out)
                        #print(line)
                        f_out.write(line)

        def dump_to_file_mass_id(self, f_name):
                #f_name = 'halo_' + self.mainID + '.full_tree'
                f_out = open(f_name, 'w')

                for iM in range(0, self.nSteps):
                        line = "%s %d" % (self.mainBranchID[iM], self.mainBranchNPart[iM]) + '\n'
                        f_out.write(line)




class MergerTree:
	# Major merger
	MM = 0

	# Formation time
	FT = 0

	# Number of steps to track the tree
	nSteps = 0

	# Number of particles per step
	nPart = np.zeros((0))

	# Normalized to z=0
	nPartNorm = np.zeros((0))

	# Normalized to z=0
	nPartSmooth = np.zeros((0))

	def __init__(self, nSteps, readPart):
		self.MM = 0
		self.FT = 0
		self.nSteps = nSteps
		self.nPart = readPart 
		self.nPartNorm = np.zeros((nSteps))
		self.nPartSmooth = np.zeros((nSteps))

		print(readPart)

		nNorm = float(self.nPart[0])

		iN = 0
		for thisNPart in self.nPart:
			self.nPartNorm[iN] = (1.0 * thisNPart) / nNorm
			iN += 1

		# Automatically smooth the tree when reading in
		self.smooth_tree()

	# We are assuming that N_Particles(z = 0) is the first element, nPart[0]
	def last_major_merger(self, smooth):
		MajorMerg = 0.1 

		if smooth == True:
			thesePart = self.nPartSmooth
		else:
			thesePart = self.nPartNorm
			

		for iMM in range(0, self.nSteps-1):
			nNP0 = thesePart[iMM]
			nNP1 = thesePart[iMM+1]
			
			dNP = abs(nNP1 - nNP0) / nNP1

			if dNP > MajorMerg:
				#print(self.nPartNorm[iMM-2:iMM+3])
				return iMM
	
	def formation_time(self, smooth):
		mHalf = 0.5
		mZero = 2.0

		if smooth == True:
			thesePart = self.nPartSmooth
		else:
			thesePart = self.nPartNorm
			

		indexForm = np.where(thesePart < 0.5)
		thisIndex = indexForm[0]

		return thisIndex[0]

	def smooth_tree(self):
		nPtsAvg = 2
		
		for iPts in range(0, nPtsAvg):
			self.nPartSmooth[iPts] = self.nPartNorm[iPts]

		for iPts in range(0, nPtsAvg):
			self.nPartSmooth[self.nSteps-iPts-1] = self.nPartNorm[self.nSteps-iPts-1]

		for iPts in range(nPtsAvg, self.nSteps-nPtsAvg):
			
			for iSmooth in range(iPts-nPtsAvg, iPts+nPtsAvg+1):
				if iSmooth < iPts-1:
					thisNorm = self.nPartSmooth[iSmooth] 
				else:
					thisNorm = self.nPartNorm[iSmooth] 

				if thisNorm > 1.0:
					self.nPartSmooth[iPts] += thisNorm * abs(1-0.75*(thisNorm - 1.0))
				else:
					self.nPartSmooth[iPts] += thisNorm
	
			self.nPartSmooth[iPts] /= float(2*nPtsAvg+1)
			#print(iPts, iSmooth, self.nPartSmooth[iPts], self.nPartNorm[iPts])

		return self.nPartSmooth

	def print_tree(self):
		
		tree_out = 'halo_' + str(ID)

#		for i_step in range(0, self.nSteps):


	def info(self):
		print('MTree with: %d steps' % self.nSteps)

