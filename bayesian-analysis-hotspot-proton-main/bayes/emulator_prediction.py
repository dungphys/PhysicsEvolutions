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
from scipy.linalg import lapack
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from global_variables import *
from params_design import design
from emulator import GP_Emulator as GPE
import subprocess
from scipy import stats
from pyDOE import *
from mcmc import *

def main():
	print ('argument list', sys.argv)
	scenario = sys.argv[1]
	scenario_prior = sys.argv[2]
	scheme_wf = sys.argv[3]
	scheme_Nh = sys.argv[4]
	scheme_Qs = sys.argv[5]
	note = sys.argv[6]
	
	X = [2.857, 0.076, 2.31, 0.001, 1.608]
	X_mod = [2.857, 0.076, 2.31, 0.2, 1.608]
	X = np.array(X, ndmin=2)
	X_mod = np.array(X_mod, ndmin=2)
	
	if scheme_wf == "rel":
			n_trains = NUM_TRAINED_POINTS_REL
	else:
			n_trains = NUM_TRAINED_POINTS_NONREL
	emu = GPE(scenario,scenario_prior,scheme_Nh,scheme_Qs,scheme_wf,note=note,n_train=n_trains)
	y_emul, cov = emu.predict(X, return_cov = True, extra_std = 0.0)
	y_emul_mod, cov_mod = emu.predict(X_mod, return_cov = True, extra_std = 0.0)
	name = "emulator_predict_%s_%s_%s_%s_%s.txt" % (scenario,scheme_Nh,scheme_Qs,scheme_wf,note) 
	with open(name, "w") as f:
		f.write("MAP\tmMAP\n")
		for i in range(len(y_emul[0])):
			f.write("%.4e\t%.4e\n" % (y_emul[0][i],y_emul_mod[0][i]))

if __name__ == '__main__':
	main()
