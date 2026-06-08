import logging
import sys
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
from pyDOE import *
from emulator import GP_Emulator as GPE


def check_symmetric(a, rtol=1e-06, atol=1e-09):
	return np.allclose(a, a.T, rtol=rtol, atol=atol)

def mvn_loglikelihood(y, cov):
	"""
	Log-likelihood functions
	L = -0.5*y.T*(cov^-1)*y - 0.5*log(det(cov)) - (n/2)* log(2*pi)
	Form MCMC, the last normalization is not important
	 	
	"""
	# check symmetry
	#print(cov.shape)
	#print(y.shape)
	symmetry = check_symmetric(cov)
	if not symmetry:
		print("cov = ")
		print(cov)
		print("Covariance matrix is not symmetric")
		sys.exit(1)
	# Cholesky decomp.
	L, info = lapack.dpotrf(cov, clean=False)
	if info < 0:
		raise ValueError('lapack dpotrf err: the {}-th argument had an illegal value'.format(-info))
	elif info > 0:
		raise np.linalg.LinAlgError('lapack dpotrf error: the leading minor of order {} is not positive definite'.format(info))
	
	# Solve for alpha = cov^-1.y using the Cholesky decomp.
	alpha, info = lapack.dpotrs(L, y)
	if info != 0:
		raise ValueError('lapack dpotrs err: the {}-th argument had an illegal value'.format(-info))
	
	return -.5*np.dot(y, alpha) - np.log(L.diagonal()).sum()

class mcmc_chain:
	def __init__(self, scenario, scenario_prior, scheme_Nh,scheme_Qs,scheme_wf,note,nwalkers=1000, nsteps = 3800, index=0):
		self.nsteps = nsteps	# number of steps for each walker
		self.nwalkers = nwalkers	# number of walkers
		self.scenario = scenario
		self.scenario_prior = scenario_prior
		self.scheme_Nh = scheme_Nh
		self.scheme_Qs = scheme_Qs
		self.scheme_wf = scheme_wf
		self.note = note  
		prior_infile = PRIOR_DIR_ + "%s_%s_%s_%s_%s.txt" % (self.scenario_prior, self.scheme_Nh, self.scheme_Qs, self.scheme_wf,self.note)
		prior_df = pd.read_csv(prior_infile,sep='\t')
		self.label = prior_df["PARAMS"].tolist() 
		self.min = prior_df["MIN"].to_numpy()
		self.max = prior_df["MAX"].to_numpy()
		self.ndim = len(self.label)
		self.labeltex = [0 for i in range(self.ndim)]
		#To latex
		for i in range(self.ndim):
			match self.label[i]:
				case "sigma_Qs":
					self.labeltex[i] = r"$\sigma_s$"
				case "R0":
					if "no_evolution" in self.scenario:
						self.labeltex[i] = r"$R_p ~\rm[GeV^{-1}]$"
					else:
						self.labeltex[i] = r"$R_0 ~\rm[GeV^{-1}]$"
				case "rH0":
					if ("no_evolution" in self.scenario) or (self.scenario=="Nh_evolution"):
						self.labeltex[i] = r"$r_H ~\rm[GeV^{-1}]$"
					else:
						self.labeltex[i] = r"$r_{H0} ~\rm[GeV^{-1}]$"
				case "gmu0":
					self.labeltex[i] = r"$g\mu_0$"
				case "m2":
					self.labeltex[i] = r'$m^2 ~\rm[GeV^2]$'
				case "lambda_p":
					self.labeltex[i] = r'$\lambda_{p}$'
				case "lambda_H":
					self.labeltex[i] = r'$\lambda_{H}$'
				case "p0":
					if (self.scenario=="Nh_evolution"):
						self.labeltex[i] = r'$p_0$'
					else:
						if (self.scheme_Nh=="with_Nh_fluct"):
							self.labeltex[i] = r'$\bar{N}_h$'
						else:
							self.labeltex[i] = r'$N_h$'
				case "p1":
					self.labeltex[i] = r'$p_1$'
				case "p2":
					self.labeltex[i] = r'$p_2$'
 
		self._y_exp, self._cov_exp = self._load_exp_data("HERA_coherent_data.txt","HERA_incoherent_data.txt")
		self.nobs = len(self._y_exp)
		print("Run emulator")
		if scheme_wf == "rel":
			n_trains = NUM_TRAINED_POINTS_REL
		else:
			n_trains = NUM_TRAINED_POINTS_NONREL
		self.emu = GPE(scenario,scenario_prior,scheme_Nh,scheme_Qs,scheme_wf,note=note,n_train=n_trains)
		self.mcsampler = emcee.EnsembleSampler(
			self.nwalkers, 
			self.ndim, 
			self.log_posterior
		)
		lsum = 0
		#fmt = '{:0' + str(len(str(n_train - 1))) + 'd}'
		par_infile = TRAINING_DIR_ + "training_params_%s_%s_%s_%s_%s_%04d.txt" % (scenario,scheme_Nh,scheme_Qs,scheme_wf,note,0)
		par_df = pd.read_csv(
					par_infile, 
					names=['par_name', 'value'], 
					header=None, 
					sep='\t'
				)
		key_values = par_df["value"].tolist()
		model_output_filename = TRAINING_DIR_ + "training_data_%s_%s_%s_%s_%s_%04d.txt" % (scenario,scheme_Nh,scheme_Qs,scheme_wf,note,0)
		output_df = pd.read_csv(
					model_output_filename, 
					sep='\t'
				)
		model_output = output_df["dsigma/dt"].tolist()
		l = len(model_output)
		for j in range(l):
			lsum += (-np.power(self._y_exp[j]-model_output[j],2)/self._cov_exp[j,j])
		max_sum = lsum
		self.init_param_val = key_values
		print(lsum)
		for i in range(1,n_trains):
			lsum = 0
			#fmt = '{:0' + str(len(str(n_train - 1))) + 'd}'
			par_infile = TRAINING_DIR_ + "training_params_%s_%s_%s_%s_%s_%04d.txt" % (scenario,scheme_Nh,scheme_Qs,scheme_wf,note,i)
			par_df = pd.read_csv(
						par_infile, 
						names=['par_name', 'value'], 
						header=None, 
						sep='\t'
					)
			key_values = par_df["value"].tolist()
			model_output_filename = TRAINING_DIR_ + "training_data_%s_%s_%s_%s_%s_%04d.txt" % (scenario,scheme_Nh,scheme_Qs,scheme_wf,note,i)
			output_df = pd.read_csv(
						model_output_filename, 
						sep='\t'
					)
			model_output = output_df["dsigma/dt"].tolist()
			l = len(model_output)
			for j in range(l):
				lsum += (-np.power(self._y_exp[j]-model_output[j],2)/self._cov_exp[j,j])
			print(lsum)
			if lsum>max_sum:
				self.init_param_val = key_values
				max_sum = lsum
			self.init_param_val = np.array(self.init_param_val)
		
	def _load_exp_data(self,filename_coh,filename_incoh):
		print('load exp data!')
		filepath_coh = EXP_DIR_ + filename_coh
		filepath_incoh = EXP_DIR_ + filename_incoh
		df_coh = pd.read_csv(filepath_coh,sep='\t')
		df_incoh = pd.read_csv(filepath_incoh,sep='\t')
		
		Q2_coh = df_coh["Q2"].tolist()
		W_coh = df_coh["W"].tolist()
		t_coh = df_coh["t"].tolist()
		ds_coh = df_coh["dsigma/dt"].tolist()
		stat_coh = df_coh["stat"].tolist()
		sys_coh =  df_coh["sys"].tolist()
		
		Q2_incoh = df_incoh["Q2"].tolist()
		W_incoh = df_incoh["W"].tolist()
		t_incoh = df_incoh["t"].tolist()
		ds_incoh = df_incoh["dsigma/dt"].tolist()
		stat_incoh = df_incoh["stat"].tolist()
		sys_incoh =  df_incoh["sys"].tolist()
		
		ndata_incoh = len(Q2_incoh) 
		ndata_coh = len(Q2_coh) 
		ndata = ndata_incoh + ndata_coh
		
		dsigma = [] 
		err2 = []
		
		if self.scenario == "no_evolution_W78":
			for i in range(ndata_coh):
				if (Q2_coh[i]==0.1) and (W_coh[i]==78):
					dsigma.append(ds_coh[i])
					err2.append(stat_coh[i]*stat_coh[i] + sys_coh[i]*sys_coh[i])
			for i in range(ndata_incoh):
				if (Q2_incoh[i]==0.1) and (W_incoh[i]==78):
					dsigma.append(ds_incoh[i])
					err2.append(stat_incoh[i]*stat_incoh[i] + sys_incoh[i]*sys_incoh[i])
		elif self.scenario == "no_evolution_W55":	
			for i in range(ndata_coh):
				if (Q2_coh[i]==0.1) and (W_coh[i]==55):
					dsigma.append(ds_coh[i])
					err2.append(stat_coh[i]*stat_coh[i] + sys_coh[i]*sys_coh[i])
			for i in range(ndata_incoh):
				if (Q2_incoh[i]==0.1) and (W_incoh[i]==55):
					dsigma.append(ds_incoh[i])
					err2.append(stat_incoh[i]*stat_incoh[i] + sys_incoh[i]*sys_incoh[i])
		else:
			dsigma = ds_coh + ds_incoh
			for i in range(ndata_coh):
				err2.append(stat_coh[i]*stat_coh[i] + sys_coh[i]*sys_coh[i])
			for i in range(ndata_incoh):
				err2.append(stat_incoh[i]*stat_incoh[i] + sys_incoh[i]*sys_incoh[i])
					
		nobs = len(dsigma)
		cov_exp = np.zeros([nobs, nobs])
		for i in range(nobs):
			cov_exp[i,i] = err2[i]
		print('data load successfully!')
		return dsigma, cov_exp
		
	def log_posterior(self, X, extra_std_prior_scale=.05):
		X = np.array(X, ndmin=2)
		in_hypercube = np.all((X > self.min) & (X < self.max), axis=1)
		lp = np.zeros(X.shape[0])
		# zero proba. for points outside the hypercube
		lp[~in_hypercube] = -np.inf 
		nsamples = np.count_nonzero(in_hypercube) # count how many point inside the hypercube
		if nsamples>0:
			extra_std = 0.0 # no addition sys. uncertainty for simplicity
			# model prediction using emulator at point(s) X inside the hypercube
			y_model, cov_model = self.emu.predict(
						X[in_hypercube], 
						return_cov = True, 
						extra_std = extra_std
					)		
			dY = np.empty([nsamples, self.nobs])	
			cov = np.empty([nsamples, self.nobs, self.nobs])
			for i in range(nsamples):
				dY[i] = y_model[i] - self._y_exp
				cov[i] = cov_model[i] + self._cov_exp 
			#print(cov.shape)
			#print(dY.shape)
			#self.cov = cov
			lp[in_hypercube] += list(map(mvn_loglikelihood, dY, cov))
			# add model sys error in form of gamma distribution (by assumption)
			lp[in_hypercube] += (2*np.log(extra_std + 1e-15) - (extra_std/extra_std_prior_scale))
		return lp
	#	#
	def run_mcmc(self,ratio_burn_in=2.0,ratio_thin=0.2):
		nburnin = 400 ##int(autocorr_len*ratio_burn_in)
		nthin = 1 #max(int(autocorr_len*ratio_thin),1)
		# creat starting points for nwalkers in the parameter space
		print("creat starting points for nwalkers in the parameter space")
		#X0 = np.array([self.init_param_val + 0.01*np.random.randn(self.ndim) for i in range(self.nwalkers)])
		lhs_sampler = qmc.LatinHypercube(d=self.ndim)
		X0_unscaled = lhs(n=self.ndim, samples=self.nwalkers, criterion = 'maximin')
		X0 = qmc.scale(X0_unscaled, self.min, self.max)
		#X0 = np.random.uniform(self.min, self.max, (self.nwalkers, self.ndim))
		print("X0 = ",X0)
		# run and compute the autocorrelation length
		print("run and compute the autocorrelation length")
		nburnin0 = nburnin//2
		self.mcsampler.run_mcmc(X0, nburnin0, progress=False)
		X0 = self.mcsampler.flatchain[
			np.unique(self.mcsampler.flatlnprobability, return_index=True)[1][-self.nwalkers:]
		]
		self.mcsampler.reset()
		#X0 = self.mcsampler.run_mcmc(X0, nburnin-nburnin0, progress=False)
				#autocorr_len = max(self.mcsampler.get_autocorr_time())
		# extract uncorrelated chains for each parameter
		#self.mcsampler.reset()
		self.mcsampler.run_mcmc(X0, self.nsteps, progress=False)
		flat_samples = self.mcsampler.get_chain(flat=True,discard=nburnin//10,thin=nthin)
		flat_prob = self.mcsampler.get_log_prob(flat=True,discard=nburnin//10,thin=nthin) 
		self.acceptance_ratio = self.mcsampler.acceptance_fraction
		return flat_samples, flat_prob
	#	#

def credible_interval(samples, label, scenario, scheme_Nh, ci=.95):
	print(samples)
	n0 = samples.shape[0]
	n1 = samples.shape[1]
	nci = int(n0*ci)
	print(nci)
	print("Samples: ", samples[:,0])
	
	samples_cp = np.zeros(samples.shape)
	for i in range(n1):
		samples_cp[:,i] = np.sort(samples[:,i])
	intervals = np.zeros((n0-nci+1,n1)) 
	for i in range(n0-nci+1): 
		for j in range(n1):
			intervals[i,j] = samples_cp[nci-1+i,j] - samples_cp[i,j] 
			
	ihpd = np.zeros(n1,dtype=int)
	MED = np.zeros(n1)
	MAP = np.zeros(n1)
	hi = np.zeros(n1)
	lo = np.zeros(n1)
	for i in range(n1):
		if (label[i] == "p0") and (scenario=="no_evolution") and (scheme_Nh=="without_Nh_fluct"):
			bins=100
		else:
			bins=100
		ihpd[i] = np.argmin(intervals[:,i])
		hi[i] = samples_cp[nci-1+ihpd[i],i]
		lo[i] = samples_cp[ihpd[i],i]
		MED[i] = np.median(samples_cp[:,i])
		hist, edges = np.histogram(samples[:,i],bins=bins)
		MAP[i] = (edges[np.argmax(hist)] + edges[np.argmax(hist)+1])/2.		

	return MED, MAP, hi, lo

def makeplot(samples, params_labels, center, hi, lo, name):

	from matplotlib import rcParams

	rcParams["savefig.dpi"] = 300
	rcParams["figure.dpi"] = 300
	rcParams["font.size"] = 16
	rcParams['font.family'] = 'sans-serif'
	rcParams["font.sans-serif"] = ["Nimbus Sans"]
	rcParams['text.usetex'] == True
		# corner plot
	ndim = len(params_labels)
	bins = [40 for i in range(ndim)]
	for i in range(ndim):
		if params_labels[i] == r"$N_h$":
			bins[i] = 8
		else:
			bins[i] = 40
	figure = corner.corner(
		samples, 
		bins=bins,
		labels = params_labels,
		color = "blue"
		#truths = center[0]
	)
	axes = np.array(figure.axes).reshape((ndim, ndim))
	# show the center value with uncertainty for each density plot in diagonal
	#for i in range(ndim):
		#title = r"MAP: $%.3f^{+%.3f}_{-%.3f}$" % (center[0][i], hi[i]-center[0][i], center[0][i]-lo[i])
	#	ax = axes[i, i]
	#	ax.set_title(title)
	plt.savefig(name,dpi=300)
	#plt.show()




		
			
			
		
		
	 
		
		
		
		
			

