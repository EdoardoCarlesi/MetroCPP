#!/usr/bin/python3

import numpy as np



class MTree:
	'The Merger tree class'
	
	nSteps = 0
	
	# Time, expansion factor and redshift
	t = np.zeros((0))
	a = np.zeros((0))
	z = np.zeros((0))
	
	# Main branch halo: number of particles and halo ID
	mainBranchID = np.zeros((0))
	mainBranchNP = np.zeros((0))
	 
	# Keep track of the steps at which the halo was lost and replaced by a token instead
	isToken = np.full((1), False)

	def __init__(self, n):
		self.nSteps = n
		self.z = np.zeros((n))
		self.a = np.zeros((n))
		self.t = np.zeros((n))
		self.isToken = np.full((n), False)

	# When did the halo experience its last major merger
	def last_major_merger(self):
		lmm = 0.0
		return lmm

	# When did the halo form
	def formation_time(self):
		ft = 0.0
		return ft

