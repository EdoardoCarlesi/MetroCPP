'''
	TODO: this is all work in progess!
	Anyway - this script should read the halo ids from a db, find the orphan ones, identify their particle content in the snapshots
	and trace it back. Ideally it would compute center of mass, position, velocity etc. and use this to replace the "original" real halo 
'''

import numpy as np
# import readsnap 

# List of orphan trees IDs
haloIDs = []

# We trace all the particles of each orphan halo at each step 
partIDs = np.zeros((nOrphans, nParticles))




