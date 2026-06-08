import sys
import os
import glob
import numpy as np
import pandas as pd

scenario = "no_evolution_W78"
scheme_Nh = "with_Nh_fluct"
scheme_Qs = "with_Qs_fluct"
scheme_wf = "rel"
note = "none"
name = "samples_bayes_output_%s_%s_%s_%s_%s.txt" % (scenario,scheme_Nh,scheme_Qs,scheme_wf,note)
df = pd.read_csv(name,sep='\t')
data = df.iloc[::,1::].values
proba = df.loc[::,"Proba"].values
l = len(proba)

data_filter = []
proba_filter = []

for i in range(l):
	if (data[i,4]>0.0) and (data[i,2]<6.08) :
		data_filter.append(data[i])
		proba_filter.append(proba[i])
data_filter = np.array(data_filter)
proba_filter = np.array(proba_filter)
proba_max = np.max(proba_filter)
imax = np.argmax(proba_filter) 

for i in range(len(proba_filter)):
	if proba_filter[i] >= (proba_max-0.1):
		print(data_filter[i])
	
print(proba_max)
print(data_filter[imax])	


