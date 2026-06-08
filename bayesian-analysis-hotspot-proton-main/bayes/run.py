import logging
import sys
import os
import glob
import numpy as np
import pandas as pd
import emcee
import h5py
import corner
import matplotlib.pyplot as plt
from scipy.stats import qmc
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from global_variables import *
from params_design import design
from emulator import GP_Emulator as GPE
import subprocess
from pyDOE import *
from mcmc import *



def main():
	"""
	Generating training data points
	- scenario: 
		-> no_evolution_W78
		-> no_evolution_W55
		-> Nh_evolution
		-> rH_evolution
	- scheme_wf:
		-> rel
		-> nonrel
	- scheme_Nh:
		-> with_Nh_fluct
		-> without_Nh_fluct
	- scheme_Qs:
		-> with_Qs_fluct
		-> without_Qs_fluct
	"""
	#files = glob.glob(TRAINING_DIR_+'*.txt')
	#for f in files:
	#	os.remove(f)
	print ('argument list', sys.argv)
	scenario = sys.argv[1]
	scenario_prior = sys.argv[2]
	scheme_wf = sys.argv[3]
	scheme_Nh = sys.argv[4]
	scheme_Qs = sys.argv[5]
	note = sys.argv[6]
	
	## designed parameter space
	#scenario = "no_evolution_W78"
	#scenario_prior = "no_evolution"	
	#scheme_wf = "rel"
	#scheme_Nh = "without_Nh_fluct"
	#scheme_Qs = "with_Qs_fluct"
	if scheme_wf == "rel":
		n_trains = NUM_TRAINED_POINTS_REL
	else:
		n_trains = NUM_TRAINED_POINTS_NONREL
	root_name_params = "training_params_%s_%s_%s_%s_%s" % (scenario,scheme_Nh,scheme_Qs,scheme_wf,note)
	root_name_data = "training_data_%s_%s_%s_%s_%s" % (scenario,scheme_Nh,scheme_Qs,scheme_wf,note)
	s = design(scenario=scenario_prior, scheme_Nh=scheme_Nh, scheme_Qs=scheme_Qs, scheme_wf=scheme_wf, note=note, n_train=n_trains)
	s.write_to_file(TRAINING_DIR_,root_name_params)
	print("Training parameter generated successfully!")	
	
	# implementing C++ code for the model outputs at design points 
	command = "./run_prior"
	for i in range(n_trains):
		print("Training point k = %04d" % i)
		infile = "%s_%04d.txt" % (root_name_params,i)
		outfile = "%s_%04d.txt" % (root_name_data,i)
		subprocess.run([CXX_DIR_ + command, infile, outfile, scenario, scheme_wf, scheme_Nh, scheme_Qs],cwd = TRAINING_DIR_)
	print("Training model outputs generated successfully!")
	print("%d generated training (couple) files are stored in %s"% (n_trains,TRAINING_DIR_))
	## 

if __name__ == '__main__':
	main()
