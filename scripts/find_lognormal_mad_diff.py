import pickle
import numpy as np
#import matplotlib.pyplot as plt
#import math
#import scipy.special as sp
import matplotlib.pyplot as plt
import statsmodels.stats.multitest
import gene_prev_utils


data_directory = "C:/Users/sarah/Garud Lab/"
#data_directory = "/Users/williamrshoemaker/Desktop/"
filename = data_directory+"lognormal_params_dict.pickle"
with open(filename, 'rb') as handle:
    lognormal_params_dict = pickle.load(handle)

n_iter = 100000

means = [lognormal_params_dict['Africa']['mu'], lognormal_params_dict['North America']['mu']]
sigmas = [lognormal_params_dict['Africa']['sigma'], lognormal_params_dict['North America']['sigma']]
cutoffs = [lognormal_params_dict['Africa']['cutoff'], lognormal_params_dict['North America']['cutoff']]
mads = [lognormal_params_dict['Africa']['mad'], lognormal_params_dict['North America']['mad']]
cov = lognormal_params_dict['Log-Cov Matrix'][0,1]


cov_matrix = np.asarray([ [lognormal_params_dict['Africa']['sigma']**2, cov], [cov, lognormal_params_dict['North America']['sigma']**2 ]  ])

multi_norm_sample = np.random.multivariate_normal(means, cov_matrix, n_iter)
multi_lognorm_sample = np.exp(multi_norm_sample)
ratio_mad = multi_lognorm_sample[:,0] / multi_lognorm_sample[:,1]
ratio_mad_log10 = np.log10(ratio_mad)
upper_ci = np.percentile(ratio_mad_log10, 97.5)
lower_ci = np.percentile(ratio_mad_log10, 2.5)


fig, ax = plt.subplots(figsize=(4,4))
ax.hist(ratio_mad_log10, alpha=0.8, bins= 20, density=True)
#ax.set_xscale('log', basex=10)
ax.vlines([upper_ci, lower_ci], ymin=[0,0], ymax=[0.28, 0.28], color ='red')
ax.set_xlabel('Log of ratio of mean abundances, '  + r'$\frac{\bar{x}_{\mathrm{Africa}}}{\bar{x}_{\mathrm{America}}}$', fontsize=12)
ax.set_ylabel('Probability density', fontsize=12)
fig_name = data_directory+ 'sim_log_ratio_dist.png'
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
#plt.show()
plt.close()

#plot hist of emperical data log ratio
species_ordered = []
ratio_emp_mad = []
for species in lognormal_params_dict['Species Abundance Dict']:
    species_ordered.append(species)
    ratio_emp_mad.append(lognormal_params_dict['Species Abundance Dict'][species]["Africa"] / 
                         lognormal_params_dict['Species Abundance Dict'][species]["North America"])
    
#ratio_emp_mad = mads[0] / mads[1]
ratio_emp_mad_log10 = np.log10(ratio_emp_mad)

fig, ax = plt.subplots(figsize=(4,4))
ax.hist(ratio_emp_mad_log10, alpha=0.8, bins= 20, histtype='barstacked', rwidth=0.9, density=True)
#ax.set_xscale('log', basex=10)
#ax.vlines([upper_ci, lower_ci], ymin=[0,0], ymax=[0.45, 0.45], color ='red')
ax.vlines([0.0, np.mean(ratio_mad_log10)], ymin=[0,0], ymax=[0.45, 0.45], colors = ['black', 'black'], linestyles = ['dotted', 'dashed'], linewidth = 2)

###add tick marks where evolutionarily relevent species are
evo_species = gene_prev_utils.load_predicted_prevalence_subsample_dict(want_names=True)
tick_values = []
tick_labels = []
for i, name in enumerate(species_ordered):
    if name in evo_species:
        tick_values.append(ratio_emp_mad_log10[i])
        tick_labels.append(name)

regular_ticks = list(ax.get_xticks())
all_ticks = regular_ticks + tick_values
tick_colors = ['black','black','black','black','black','black','red','red','red','red','red','red','red','red']
ax.vlines(tick_values, ymin=0, ymax=0.1, colors = 'red', linewidth = 0.8) 
ax.set_xlabel('Log of ratio of mean abundances, '  + r'$\frac{\bar{x}_{\mathrm{Africa}}}{\bar{x}_{\mathrm{America}}}$', fontsize=12)
ax.set_ylabel('Probability density', fontsize=12)
fig_name = data_directory+ 'emp_log_ratio_dist.png'
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.show()
plt.close()


###calculate corrected p-values for each mean abundance observed

p_vals = []
for val in ratio_emp_mad_log10:
    
    p_val = (ratio_mad_log10 > val).sum()
    if p_val > (n_iter/2.0):
        p_val = n_iter - p_val

    p_vals.append(p_val/(n_iter*1.0))
    
#correct p-values
reject_bool, corrected_p_vals, sidak_alpha, bonf_alpha = statsmodels.stats.multitest.multipletests(p_vals, alpha=0.05,method = 'fdr_bh')

num_sig_enriched_species = reject_bool.sum()
#print(corrected_p_vals)
print(str(num_sig_enriched_species) + " enriched species")
        
