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
	#scenario = "no_evolution_W78"
	#scenario_prior = "no_evolution"
	#scheme_wf = "rel"
	#scheme_Nh = "with_Nh_fluct"
	#scheme_Qs = "without_Qs_fluct"
	while True:
		chain = mcmc_chain(scenario=scenario,scenario_prior=scenario_prior,scheme_Nh=scheme_Nh,scheme_Qs=scheme_Qs,scheme_wf=scheme_wf,note=note)
		samples, proba = chain.run_mcmc()
		print ("Acceptance ratios:")
		print("\tmean: ", chain.acceptance_ratio.mean())
		print("\tmin: ", chain.acceptance_ratio.min())
		print("\tmax: ",chain.acceptance_ratio.max())
		print("\tstd: ",chain.acceptance_ratio.std())
		if np.mean(chain.acceptance_ratio) > 0.2:
			break
	#for i in range(len(chain.label)):
		#if chain.label[i]=="p0" and scenario_prior == "no_evolution" and scheme_Nh =="without_Nh_fluct":
			#samples[:,i] = np.round(samples[:,i])
	criteria = 'MAP' 
	MED, MAP, hi, lo = credible_interval(samples, chain.label, scenario_prior, scheme_Nh, ci=.95)
	#for i in range(len(chain.label)):
		#if chain.label[i]=="p0" and scenario_prior == "no_evolution" and scheme_Nh =="without_Nh_fluct":
			#MAP[i] = stats.mode(samples[:,i], keepdims=True, axis=None)[0].item()
	center = [MAP, MED]
	name = "params_bayes_output_%s_%s_%s_%s_%s.txt" % (scenario,scheme_Nh,scheme_Qs,scheme_wf,note) 
	with open(name, "w") as f:
		f.write("Params\tMAP\terr_MAP_m\terr_MAP_p\tMED\terr_MED_m\terr_MED_p\n")
		for i in range(len(MAP)):
			f.write("%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (chain.label[i],MAP[i],MAP[i]-lo[i],hi[i]-MAP[i],MED[i],MED[i]-lo[i],hi[i]-MED[i]))
	name = "samples_bayes_output_%s_%s_%s_%s_%s.txt" % (scenario,scheme_Nh,scheme_Qs,scheme_wf,note)	
	label = chain.labeltex	
	with open(name, "w") as f:
		f.write("Type")
		for i in range(len(label)):
			f.write("\t%s"%(label[i]))
		f.write("\tProba")
		f.write("\n")
		for i in range(len(samples)):
			f.write("%s"%(scheme_wf))
			for j in range(len(samples[i])):
				f.write("\t%.6f"%(samples[i][j]))
			f.write("\t%.6f"%(proba[i]))
			f.write("\n")
	name = "CL_%s_%s_%s_%s_%s.txt" % (scenario,scheme_Nh,scheme_Qs,scheme_wf,note) 
	with open(name, "w") as f:
		f.write("Params\tlow\thight\n")
		for i in range(len(MAP)):
			f.write("%s\t%.3f\t%.3f\n" % (chain.label[i],lo[i],hi[i]))
	name = "corner_%s_%s_%s_%s_%s.pdf" % (scenario,scheme_Nh,scheme_Qs,scheme_wf,note) 
	print ("Acceptance ratios:")
	print("\tmean: ", np.mean(chain.acceptance_ratio))
	print("\tmin: ", min(chain.acceptance_ratio))
	print("\tmax: ",max(chain.acceptance_ratio))
	makeplot(samples, chain.labeltex, center, hi,lo, name)

if __name__ == '__main__':
	main()
