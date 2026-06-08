import sys
import os
import glob
import numpy as np
import pandas as pd
import h5py
import corner
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
from matplotlib import rcParams
import seaborn as sns

CORNER_KWARGS = dict(
    smooth=0.95,
    label_kwargs=dict(fontsize=16),
    title_kwargs=dict(fontsize=16),
    plot_density=True,
    plot_datapoints=False,
    fill_contours=True,
    show_titles=False,
    labelpad=0.1,
    max_n_ticks=5,
    top_ticks = False,
    title_fmt=".2E"
)

def get_normalisation_weight(len_current_samples, len_of_longest_samples):
    return np.ones(len_current_samples) * (len_of_longest_samples / len_current_samples)



def overlaid_corner(samples_list, sample_labels, name, labels, bins, ranges):
    """Plots multiple corners on top of each other"""
    # get some constants
    n = len(samples_list)
    _, ndim = samples_list[0].shape
    #print(samples_list[0])
    #print(ndim)
    max_len = max([len(s) for s in samples_list])
    print(len(samples_list[0]))
    print("max len, ",max_len)
    #cmap = plt.get_cmap('gist_rainbow', n)
    colors = ["red", "blue", "green"]

    plot_range = []
    for dim in range(ndim):
        plot_range.append(
            [
                min([min(samples_list[i].T[dim]) for i in range(n)]),
                max([max(samples_list[i].T[dim]) for i in range(n)]),
            ]
        )
	
    #CORNER_KWARGS.update(range=ranges[0])

    fig = corner.corner(
        samples_list[0],
        bins=bins,
        labels=labels,
        color=colors[0],
        #range=ranges[0],
        **CORNER_KWARGS
    )
    
    

    for idx in range(1, n):
        fig = corner.corner(
            samples_list[idx],
            fig=fig,
            bins=bins,
            labels=labels,
            #range=ranges[idx],
            reverse=False,
            color=colors[idx],
            **CORNER_KWARGS
        )

    plt.legend(
        handles=[
            mlines.Line2D([], [], color=colors[i], label=sample_labels[i])
            for i in range(n)
        ],
        fontsize=25, frameon=False,
        bbox_to_anchor=(1, ndim), loc="upper right"
    )
    plt.savefig(name,bbox_inches='tight')
    plt.show()
    plt.close()

def sample_maxcount(samples):
	sample_unique = np.unique(samples)
	cnt = [0 for i in range(len(sample_unique))]
	for i in range(len(sample_unique)):
		for j in range(len(samples)):
			if sample_unique[i] == samples[j]:
				cnt[i] += 1
	return sample_unique[np.argmax(cnt)]  

def credible_interval(samples, bins, ci=.95):
	#print(len(samples))
	#print(samples[0],samples[5])
	nsample = len(samples)
	nci = int(ci*nsample)
	samples_cp = np.sort(samples)
	interval = np.zeros(nsample-nci+1)
	for i in range(nsample-nci+1):
		interval[i] = samples_cp[nci-1+i] - samples_cp[i]
	ihpd = np.argmin(interval)
	hi = samples_cp[nci-1+ihpd]
	lo = samples_cp[ihpd]
	MED = np.median(samples_cp)
	hist, edges = np.histogram(samples,bins=bins)
	print("hist ok!")
	if bins == 80:
		#MAP = (edges[np.argmax(hist)] + edges[np.argmax(hist)+1])/2.
		MAP = np.mean(samples_cp)
	else:
		#MAP = (edges[np.argmax(hist)+1])#sample_maxcount(samples_cp)
		MAP = np.mean(samples_cp)
	#print(MAP)
	return hi, lo, MED, MAP

def main():
	scenario = "no_evolution_W78"
	scheme_Nh = "without_Nh_fluct"
	scheme_Qs = "without_Qs_fluct"
	note = "Nvaried"
	name_rel = "samples_bayes_output_%s_%s_%s_rel_%s.txt" % (scenario,scheme_Nh,scheme_Qs,note)
	df_rel = pd.read_csv(name_rel,sep='\t')
	name_nonrel = "samples_bayes_output_%s_%s_%s_nonrel_%s.txt" % (scenario,scheme_Nh,scheme_Qs,note)
	df_nonrel = pd.read_csv(name_nonrel,sep='\t')
	print("read ok:\n")
	row, col = df_rel.shape
	label_rel = list(df_rel.columns[1:col-1])
	label_nonrel = list(df_nonrel.columns[1:col-1])
	
	
	
	df = pd.concat([df_rel,df_nonrel])
	
	ucut = [0 for i in range(len(label_rel))]
	lcut = [0 for i in range(len(label_rel))]
	
	data_rel = df_rel.iloc[:,1:col-1].values
	data_nonrel = df_nonrel.iloc[:,1:col-1].values
	print("len", len(data_rel))
	#rnd = np.unique(np.random.randint(1,len(data_rel),500000))
	#data_rel = data_rel[rnd,:]
	#data_nonrel = data_nonrel[rnd,:]
	#print("len", len(data_rel))
	#rnd = np.unique(np.random.randint(1,len(data_rel),50000))
		
	for i in range(len(label_rel)):
		ucut[i] = max(float(max(data_rel[:,i])),float(max(data_nonrel[:,i])))
		lcut[i] = min(float(min(data_rel[:,i])),float(min(data_nonrel[:,i])))
	# for w Nh w Qs
	#ucut[2] = 0.21	
	#ucut[1] = 0.118
	#
	#lcut[3] = 0.08
	# for wo Nh w Qs
	#lcut[5] = 0.13
	#ucut[1] = 0.12
	#ucut[5] = 7.0
	# for wo Nh w Qs
	#ucut[4] = 4.7
	#ucut[2] = 0.51
	#lcut[4] = 2.1
	#ucut[4] = 4.6 
		
	#print(ucut)
	data_rel_renorm = []
	data_nonrel_renorm = []
	nsamples_rel = len(data_rel)
	nsamples_nonrel = len(data_nonrel)
	nlabel = len(label_rel)
	for i in range(nsamples_rel):
		count_true = 0
		for j in range(nlabel):
			if data_rel[i,j]<=ucut[j] and data_rel[i,j]>=lcut[j]:
				count_true +=1
		if count_true ==nlabel:
			data_rel_renorm.append(data_rel[i])
	for i in range(nsamples_nonrel):
		count_true = 0
		for j in range(nlabel):
			if data_nonrel[i,j]<=ucut[j] and data_nonrel[i,j]>=lcut[j]:
				count_true +=1
		if count_true ==nlabel:
			data_nonrel_renorm.append(data_nonrel[i])
	
	data_rel_renorm = np.array(data_rel_renorm)
	data_nonrel_renorm = np.array(data_nonrel_renorm)	
	ranges = [[0 for i in range(len(label_rel))],[0 for i in range(len(label_rel))]]
	#ranges[0] = [(2.2,4.2), (0.005,0.1), (1.4,3.2), (0, 0.7), (0.5,8.5)]
	#ranges[1] = [(2.2,4.2), (0.005,0.2), (1.4,3.2), (0, 0.7), (0.5,8.5)]
	label = [0 for i in range(len(label_rel))]
	
	bins = [80 for i in range(len(label_rel))]
	for i in range(len(label_rel)):
		ranges[0][i] = (float(min(data_rel_renorm[:,i])),float(max(data_rel_renorm[:,i])))
		ranges[1][i] = (float(min(data_nonrel_renorm[:,i])),float(max(data_nonrel_renorm[:,i])))
		if label_rel[i] == '$N_h$':
			tmp1 = int(np.round(float(max(data_rel_renorm[:,i])) - float(min(data_rel_renorm[:,i]))))
			tmp2 = int(np.round(float(max(data_nonrel_renorm[:,i])) - float(min(data_nonrel_renorm[:,i]))))
			bins[i] = 80#8*max(tmp1,tmp2) 
			#ranges[0][i] = (float(min(data_rel_renorm[:,i]))-0.5,float(max(data_rel_renorm[:,i]))+0.5)
			#ranges[1][i] = (float(min(data_nonrel_renorm[:,i]))-0.5,float(max(data_nonrel_renorm[:,i]))+0.5)
			#bins[i] = 20
			label[i] = "p0"
		if label_rel[i] == r"$g\mu_0$":
			label[i] = "gmu0"
		if label_rel[i] == r"$m^2 ~\rm[GeV^2]$":
			label[i] = "m2"
		if label_rel[i] == r"$R_p ~\rm[GeV^{-1}]$":
			label[i] = "R0"
		if label_rel[i] == r"$r_H ~\rm[GeV^{-1}]$":
			label[i] = "rH0"
		if label_rel[i] == r"$\bar{N}_h$":
			label[i] = "p0"
		if label_rel[i] == r"$\sigma_s$":
			label[i] = "sigma_Qs"
	
	cre_in_hi_rel = np.zeros(len(label_rel))
	cre_in_hi_nonrel = np.zeros(len(label_rel))
	MED_rel = np.zeros(len(label_rel))
	MAP_rel = np.zeros(len(label_rel))
	cre_in_lo_rel = np.zeros(len(label_rel))
	cre_in_lo_nonrel = np.zeros(len(label_rel))
	MED_nonrel = np.zeros(len(label_rel))
	MAP_nonrel = np.zeros(len(label_rel))
	for i in range(len(label_rel)):
		cre_in_hi_rel[i], cre_in_lo_rel[i], MED_rel[i], MAP_rel[i] = credible_interval(data_rel_renorm.T[i], bins[i])
		cre_in_hi_nonrel[i], cre_in_lo_nonrel[i], MED_nonrel[i], MAP_nonrel[i] = credible_interval(data_nonrel_renorm.T[i], bins[i])
		
	print(MED_rel, MAP_rel)
	print(MED_nonrel, MAP_nonrel)
	
	#Nhcut_rel = 4.2
	#Nhcut_nonrel = 6.1
	#Nh_nonrel = data_nonrel_renorm.T[5][data_nonrel_renorm.T[5]<=Nhcut_nonrel]
	#Nh_rel = data_rel_renorm.T[5][data_rel_renorm.T[5]<=Nhcut_rel]
	#cre_in_hi_nonrel[5], cre_in_lo_nonrel[5], MED_nonrel[5], MAP_nonrel[5] = credible_interval(Nh_nonrel, bins[5])
	#cre_in_hi_rel[5], cre_in_lo_rel[5], MED_rel[5], MAP_rel[5] = credible_interval(Nh_rel, bins[5])
	
	## save to file
	name = "params_bayes_output_%s_%s_%s_rel_%s.txt" % (scenario,scheme_Nh,scheme_Qs,note) 
	with open(name, "w") as f:
		f.write("Params\tMAP\terr_MAP_m\terr_MAP_p\tMED\terr_MED_m\terr_MED_p\n")
		for i in range(len(MAP_rel)):
			f.write("%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (label[i],MAP_rel[i],MAP_rel[i]-cre_in_lo_rel[i],cre_in_hi_rel[i]-MAP_rel[i],MED_rel[i],MED_rel[i]-cre_in_lo_rel[i],cre_in_hi_rel[i]-MED_rel[i]))
	name = "params_bayes_output_%s_%s_%s_nonrel_%s.txt" % (scenario,scheme_Nh,scheme_Qs,note) 
	with open(name, "w") as f:
		f.write("Params\tMAP\terr_MAP_m\terr_MAP_p\tMED\terr_MED_m\terr_MED_p\n")
		for i in range(len(MAP_nonrel)):
			f.write("%s\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\t%.3f\n" % (label[i], MAP_nonrel[i], MAP_nonrel[i]-cre_in_lo_nonrel[i], cre_in_hi_nonrel[i]-MAP_nonrel[i], MED_nonrel[i], MED_nonrel[i]-cre_in_lo_nonrel[i], cre_in_hi_nonrel[i]-MED_nonrel[i]))
	 
	rcParams["savefig.dpi"] = 300
	rcParams["figure.dpi"] = 300
	rcParams["font.size"] = 20
	rcParams['font.family'] = 'sans-serif'
	rcParams["font.sans-serif"] = ["Nimbus Sans"]
	rcParams['text.usetex'] == True
		# corner plot
	#g = sns.PairGrid(df, hue='Type',diag_sharey=False)
	#g.map_upper(sns.kdeplot)
	#g.map_lower(sns.kdeplot)	
	#g.map_diag(sns.kdeplot)
	#g.add_legend()
	name = "posterior_%s_%s_%s_%s.pdf" % (scenario,scheme_Nh,scheme_Qs,note) 
	leg = [r"$J/\psi$ wave function with relativistic correction", r"Non-relativistic $J/\psi$ wave function"]
	overlaid_corner([data_rel_renorm, data_nonrel_renorm], leg, name, label_rel, bins, ranges)
	#figure = corner.corner(
	#	df_rel[1:], 
	#	bins=40,
	#	labels = label_rel,
	#	color = "blue"
	#	#truths = center[0]
	#)
	
	plt.show()
	
if __name__ == '__main__':
	main()
