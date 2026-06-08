import logging
import sys
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
from sklearn.gaussian_process import GaussianProcessRegressor
from sklearn.gaussian_process.kernels import RBF, WhiteKernel
from global_variables import *

class GP_Emulator:
	"""
	Gaussian emulator for posterior distribution
	Public attributes:
	- ``scenario``		: scenario string for parameters
	- ``n_train``		: number of training points
	- ``key_names``	: parameters' names
	- ``key_values[]``	: values of training parameters
	- ``model_output[]``: values of training model outputs
	- ``min`` 		: prior min values of parameters
	- ``max`` 		: prior max values of parameters
	- ``scaler``		: StandardScaler for data (preprocessing)
	- ``pca``			: PCA machine
	- ``gps``			: Gaussian Process sampler
	- ``npc``			: number of PC
	- ``nobs``		: number of 'observables'
	"""
	def __init__(self,
			scenario,
			scenario_prior,
			scheme_Nh,
			scheme_Qs,
			scheme_wf,
			note,
			n_train,
			threshold_explvar=THRESHOLD_RATIO,
			rbf_coeff=1., 
			white_coeff=1.,
			index=0
		):
		self.scenario = scenario
		self.scenario_prior = scenario_prior
		self.n_train = n_train
		self.key_values = [0 for i in range(n_train)]
		self.model_output = [0 for i in range(n_train)]
		print("load training data")
		for i in range(n_train):
			# read training params for each training point
			fmt = '{:0' + str(len(str(n_train - 1))) + 'd}'
			par_infile = TRAINING_DIR_ + "training_params_%s_%s_%s_%s_%s_%04d.txt" % (scenario,scheme_Nh,scheme_Qs,scheme_wf,note,i)
			par_df = pd.read_csv(
						par_infile, 
						names=['par_name', 'value'], 
						header=None, 
						sep='\t'
					)
			self.key_names = par_df["par_name"].tolist()
			self.key_values[i] = par_df["value"].tolist()
			
			# load training model outputs for each training point
			model_output_filename = TRAINING_DIR_ + "training_data_%s_%s_%s_%s_%s_%04d.txt" % (scenario,scheme_Nh,scheme_Qs,scheme_wf,note,i)
			output_df = pd.read_csv(
						model_output_filename, 
						sep='\t'
					)
			self.model_output[i] = output_df["dsigma/dt"].tolist()
		## Before dim. reduction
		## 		dim(params) = ntrain x nparams
		## 		dim(output) = ntrain x nobs
		## -------------------------------------------
		self.nobs = len(self.model_output[0])	# number of observables
		self.key_values = np.array(self.key_values) # 2D array storing training params values
		self.model_output = np.array(self.model_output) # 2D array storing training model output values
		# PCA processing
		print("PCA processing")
		self.scaler = StandardScaler(copy=False)
		self.pca = PCA(copy=False, whiten=True, svd_solver='full')
		#
		scaled_data = self.scaler.fit_transform(self.model_output)
		PC_data = self.pca.fit_transform(scaled_data)
		# extract the number of PCs from the threshold for explained variance ratio
		self.npc = np.where(np.cumsum(self.pca.explained_variance_ratio_)>=threshold_explvar)[0][0] + 1
		# extract npc components
		print("nPC = %d"%self.npc)
		emulate_data = PC_data[:,:self.npc]
		emulate_data_variance = self.pca.explained_variance_[:self.npc]
		## After dim. reduction
		## 		dim(params) = ntrain x nparams
		## 		dim(output) = ntrain x npc
		## -------------------------------------------
		# kernel to used: RBF + white noise
		prior_infile = PRIOR_DIR_ + "%s_%s_%s_%s_%s.txt" % (scenario_prior, scheme_Nh, scheme_Qs, scheme_wf, note)
		prior_df = pd.read_csv(prior_infile,sep='\t')
		self.min = prior_df["MIN"].to_numpy()
		self.max = prior_df["MAX"].to_numpy()
		length_scale = self.max-self.min
		print(length_scale)
		kernel = (
			rbf_coeff * RBF(
				length_scale=length_scale, 
				length_scale_bounds=np.outer(length_scale,(0.1, 10)) 
			) +
			white_coeff * WhiteKernel(
				noise_level = 0.1**2,
				noise_level_bounds = (1e-10, 10)
			)
		)				
		# Gaussian emulator for each PC
		print("Run GPS for each PC")
		self.gps = [
			GaussianProcessRegressor(kernel=kernel, 
					alpha=1e-36, 
					n_restarts_optimizer = 25
			).fit(self.key_values, z)
			for z in emulate_data.T 
		]
		# Construct the full linear transformation matrix
		print("linear transformation matrix")
		self._trans_matrix = (
			self.pca.components_
			* np.sqrt(self.pca.explained_variance_[:, np.newaxis])
			* self.scaler.scale_
		)
		## Full transformation:
		##	dim(pca.component_) = nobs x nobs
		##	dim(pca.explained_variance_[:, np.newaxis]) = nobs x 1
		##   dim(scaler.scale_) = 1 x nobs
		## -> dim(#1*#2*#3) = nobs x nobs
		# Pre-calculate some arrays for inverse transforming the predictive
        	# variance (from PC space to physical space).
        	# Dim(A) = npc x nobs
		A = self._trans_matrix[:self.npc]
		self._var_trans = np.einsum(
			'ki,kj->kij', A, A, optimize=False
		).reshape(self.npc,self.nobs**2)
		# Compute the covariance matrix for the remaining neglected PCs
		# (truncation error).  These components always have variance == 1.		
		B = self._trans_matrix[self.npc:]
		self._cov_trunc = np.dot(B.T, B)
		# Add small term to diagonal for numerical stability.
		self._cov_trunc.flat[::self.nobs + 1] += 1e-4 * self.scaler.var_
	def _inverse_transform (self, Z):
		'''
		Transform principal components Z back to observable Y
		'''
		Y = np.dot(Z,self._trans_matrix[:Z.shape[-1]])
		Y += self.scaler.mean_
		return Y
	def predict(self, X, return_cov=False, extra_std=0):
		"""
		Predict model output at X: dim(X) = n_predict x nparams
		"""
		# for each pricipal component, make prediction at the point X 
		# return mean value and covariant matrix for each PC (so npc members)
		gp_mean = [gp.predict(X, return_cov=return_cov) for gp in self.gps]
		if return_cov:
			gp_mean, gp_cov = zip(*gp_mean)
		# 
		# dim(mean) = n_predict x nobs
		mean = self._inverse_transform(
			np.concatenate([np.array(m)[:, np.newaxis] for m in gp_mean], axis=1)
		)
		
		if return_cov:
			# collecting predictive variances (diagonal) for PC 
			gp_var = np.concatenate(
				[c.diagonal()[:,np.newaxis] for c in gp_cov],
				axis = 1
			)
			# add systematic uncertainty given by extra_std
			extra_std = np.array(extra_std).reshape(-1, 1)
			gp_var += extra_std**2
			# transform back to physical space
			cov = np.dot(gp_var, self._var_trans).reshape(
				X.shape[0], self.nobs, self.nobs
			)
			# add truncation error
			cov += self._cov_trunc
			## dim(cov) = n_predict arrays of nobs x nobs 
			return mean, cov
		else:
			return mean

def main():
	link = "/home/anh-dung/projects/test/"
	scenario = "no_evolution"
	
	E = Emulation(link,scenario)
	print(E.key_values)
	print(E.scaled_data)
	print(E.PC_data)
	print(E.npc)

if __name__ == '__main__':
	main()		
