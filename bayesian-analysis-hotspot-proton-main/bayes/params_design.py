import logging
import sys
import numpy as np
import pandas as pd
from scipy.stats import qmc
import subprocess
from pyDOE import *
from global_variables import *

class design:
	"""
	Latin-hypercube model design.
    	Public attributes:

	- ``scenario``: The scenario string to indicate the parameter schemes
    	- ``keys``: list of parameter keys
    	- ``range``: list of parameter (min, max) tuples
    	- ``min``, ``max``: numpy arrays of parameter min and max
    	- ``ndim``: number of parameters (i.e. dimensions)
    	- ``n_train``: number of training points to generate
    	- ``points``: list of design point names (formatted numbers) (for saving)
    	- ``array``: the actual design array
    	
	"""
	def __init__(self,scenario,scheme_Nh,scheme_Qs,scheme_wf,note,n_train,index=0):
		self.scenario = scenario
		self.scheme_Nh = scheme_Nh
		self.scheme_Qs = scheme_Qs
		self.scheme_wf = scheme_wf
		self.n_train = n_train
		# read prior ranges of parameters
		prior_infile = PRIOR_DIR_ + scenario + "_" + scheme_Nh + "_" + scheme_Qs + "_" + scheme_wf + "_" + note + ".txt"
		prior_df = pd.read_csv(prior_infile,sep='\t')
		self.key = prior_df["PARAMS"].tolist()
		self.min = prior_df["MIN"].to_numpy()
		self.max = prior_df["MAX"].to_numpy()			
		self.ndim = len(self.key)
		## for file naming 
		#fmt = '{:03d' + str(len(str(n_train - 1))) + 'd}'
		fmt = '{:04d}'
		self.points = [fmt.format(i) for i in range(n_train)]
		## Latin hypercube sampling for the above parameters and prior ranges 
		#sampler = qmc.LatinHypercube(d=self.ndim) # in [0,1]^ndim
		sample = lhs(n=self.ndim, samples=n_train, criterion="maximin")#lhs_sampler.random(n_train)	
		self.array = qmc.scale(sample, self.min, self.max) # scaled to the prior ranges
	def __array__(self):
		return self.array
	def write_to_file(self, destination, root_name):
		for i in range(self.n_train):
			name = "%s_%s.txt" % (root_name,self.points[i])
			f = open(destination + name,"w")
			for j in range(len(self.key)):
				if (self.scenario=="no_evolution") and (self.scheme_Nh == "without_Nh_fluct") and (self.key[j]=="p0"):
					f.write("%s\t%.6lf\n" % ( self.key[j], self.array[i,j] ))
				else:
					f.write("%s\t%.6lf\n" % ( self.key[j], self.array[i,j] ))
			f.close()

def run_Cpp(command, command_loc, inputfile, inputfile_loc):
	subprocess.run([command_loc + command, inputfile],cwd = inputfile_loc)	

def main():
	to_dir = "/home/anh-dung/projects/test/"
	#s = design(scenario="no_evolution")
	s = design(scenario="no_evolution")
	s.write_to_file(to_dir)
	
	command = "./test_wo"
	command_loc = "/home/anh-dung/projects/current_projects/hot_spot_NR/bin/"
	infile = "no_evolution_params_train_487.txt"
	infile_loc = "/home/anh-dung/projects/test/"
	subprocess.run([command_loc + command, infile],cwd = infile_loc)
	
if __name__ == '__main__':
	main()
