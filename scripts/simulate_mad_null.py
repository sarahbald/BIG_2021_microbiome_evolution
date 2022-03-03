import pickle
import numpy as np
#import matplotlib.pyplot as plt
#import math


data_directory = "C:/Users/sarah/Garud Lab/"
filename = data_directory+"lognormal_params_dict.pickle"
with open(filename, 'rb') as handle:
        lognormal_params_dict = pickle.load(handle)

means = [lognormal_params_dict['Africa']['mu'], lognormal_params_dict['North America']['mu']]
sigmas = [lognormal_params_dict['Africa']['sigma'], lognormal_params_dict['North America']['sigma']]
cutoffs = [lognormal_params_dict['Africa']['cutoff'], lognormal_params_dict['North America']['cutoff']]
mads = [lognormal_params_dict['Africa']['mad'], lognormal_params_dict['North America']['mad']]
#cov = lognormal_params_dict['Log-Cov Matrix']

#cov = lognormal_params_dict['Log-Cov Matrix'][0,1]
#cov_matrix = np.asarray([ [lognormal_params_dict['Africa']['sigma']**2, cov], [cov, lognormal_params_dict['North America']['sigma']**2 ]  ])
#multi_norm_sample = np.random.multivariate_normal(means, cov_matrix, 100000)
#multi_lognorm_sample = np.exp(multi_norm_sample)
#ratio_mad = (multi_lognorm_sample[:,0]) / (multi_lognorm_sample[:,1])
#ratio_mad_log10 = np.log10(ratio_mad)
#​fig, ax = plt.subplots(figsize=(4,4))
#​ax.hist(ratio_mad_log10, alpha=0.8, bins= 20, density=True)
##ax.set_xscale('log', basex=10)
​
#ax.set_xlabel('Log of ratio of mean abundances, '  + r'$\frac{\bar{x}_{\mathrm{Africa}}}{\bar{x}_{\mathrm{America}}}$', fontsize=12)
#ax.set_ylabel('Probability density', fontsize=12)
​
#fig_name = data_directory+ 'test.png'
#fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
#plt.close()


#exp_multinorm_sample = []
#for sample in multi_norm_sample:
#    exp_multinorm_sample.append(np.exp(sample))
    
#af_sim_sample = [val[0] for val in exp_multinorm_sample]
#na_sim_sample = [val[1] for val in exp_multinorm_sample]



