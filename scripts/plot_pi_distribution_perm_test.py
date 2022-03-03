import numpy as np
import random
import pickle
import matplotlib.pyplot as plt

data_directory = "C:/Users/sarah/Garud Lab/BIG_2021_microbiome_evolution/"
output_directory = "C:/Users/sarah/Garud Lab/BIG_2021_microbiome_evolution/plots/"
#data_directory = "/Users/williamrshoemaker/Desktop/"
filename = data_directory+"pi_dict.pickle"
with open(filename, 'rb') as handle:
    pi_dict = pickle.load(handle)
    
n_iter = 10000
    
obs_log10_ratio_mean_pi = {}
for species in pi_dict.keys():
    mean_na = np.mean(pi_dict[species]['North America'].values())
    mean_af = np.mean(pi_dict[species]['Africa'].values())
    obs_log10_ratio_mean_pi[species] = np.log10(mean_af / mean_na)
    
    #shuffle labels in groups to get null distribution of log ratio mean pi
    n_na = len(pi_dict[species]['North America'].values())
    n_af = len(pi_dict[species]['Africa'].values())
    all_pi = pi_dict[species]['North America'].values() + pi_dict[species]['Africa'].values()
    perm_log_ratio_pi = []
    for i in range(n_iter):
        af_idx = random.sample(range(len(all_pi)), n_af)
        perm_af_mean = 0.0
        perm_na_mean = 0.0
        for i, pi in enumerate(all_pi):
            if i in af_idx:
                perm_af_mean += pi
            else:
                perm_na_mean += pi     
        perm_af_mean = perm_af_mean/n_af
        perm_na_mean = perm_na_mean/n_na
        perm_log_ratio_pi.append(np.log10(perm_af_mean / perm_na_mean))
    
    p_val = float((perm_log_ratio_pi > obs_log10_ratio_mean_pi[species]).sum())
    p_val = p_val/n_iter
    #print("p-val for %s is %s" % (species, p_val))
    
    fig, ax = plt.subplots(figsize=(4,4))
    ax.hist(perm_log_ratio_pi, bins = 20, histtype = 'barstacked', rwidth = 0.85, alpha = 0.8, density=True)
    axes = plt.gca()
    y_min, y_max = axes.get_ylim()
    #if p_val < 0.05:
    #    ax.text(0.98*obs_log10_ratio_mean_pi[species], (y_max*0.5), '*', size = 'medium')
    ax.vlines(obs_log10_ratio_mean_pi[species], ymin = y_min, ymax = y_max, color = 'red')
    ax.vlines(0, ymin = y_min, ymax = y_max, color = 'black', linestyle = '--')
    ax.set_xlabel('Log ratio of mean synonymous diversity, ' + r'$\frac{\bar{x}_{\mathrm{Africa}}}{\bar{x}_{\mathrm{America}}}$', fontsize=12)
    ax.set_ylabel('Probability density', fontsize=12)
    plt.title(str(species))
    filename = output_directory + species + '_pi_distribution.png'
    plt.savefig(filename, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
    plt.show()
    plt.clf()
    plt.close()
    
        

    

