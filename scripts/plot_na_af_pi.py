import numpy as np
import pickle
import seaborn
from scipy.stats import sem
import matplotlib.pyplot as plt

data_directory = "C:/Users/sarah/Garud Lab/BIG_2021_microbiome_evolution/"
output_directory = "C:/Users/sarah/Garud Lab/BIG_2021_microbiome_evolution/plots/"
#data_directory = "/Users/williamrshoemaker/Desktop/"
filename = data_directory+"pi_dict.pickle"
with open(filename, 'rb') as handle:
    pi_dict = pickle.load(handle)
 
#setup color dictionary
good_species = pi_dict.keys()   
species_color_dict = {}
colors = seaborn.color_palette('muted', n_colors=len(good_species))
i = 0
good_species.sort()
for species in good_species:
    species_color_dict[species] = colors[i]
    i+=1
    
#get standard error bars
mean_pi_dict = {}
for species in pi_dict.keys():
    na_pi = pi_dict[species]['North America'].values()
    log10_na_pi = np.log10(pi_dict[species]['North America'].values())
    af_pi = pi_dict[species]['Africa'].values()
    log10_af_pi = np.log10(pi_dict[species]['Africa'].values())
    print("Africa: min pi for %s is %s and max is %s" % (species, min(af_pi), max(af_pi)))
    
    mean_na = np.mean(na_pi)
    mean_af = np.mean(af_pi)
    mean_pi_dict[species] = {}
    mean_pi_dict[species]['North America'] = {}
    mean_pi_dict[species]['North America']['mean'] = mean_na
    mean_pi_dict[species]['North America']['std error'] = sem(na_pi)
    mean_pi_dict[species]['Africa'] = {}
    mean_pi_dict[species]['Africa']['mean'] = mean_af
    mean_pi_dict[species]['Africa']['std error'] = sem(af_pi)
 

fig, ax = plt.subplots(figsize = (6,6))
for species in pi_dict.keys():
    ax.scatter(mean_pi_dict[species]['Africa']['mean'], mean_pi_dict[species]['North America']['mean'], color = species_color_dict[species], label = species)
    ax.errorbar(mean_pi_dict[species]['Africa']['mean'], mean_pi_dict[species]['North America']['mean'], yerr=(mean_pi_dict[species]['North America']['std error']), xerr=(mean_pi_dict[species]['Africa']['std error']), color = species_color_dict[species])

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Synonymous Diversity in Non-Western Hosts')
ax.set_ylabel('Synonymous Diversity in Western Hosts')
ax.plot([10**(-2.5), 10**-1], [10**(-2.5), 10**-1], color = "grey", alpha = 0.5)
ax.legend(bbox_to_anchor = (1.05,1), loc='upper left')
filename = output_directory + "na_versus_af_pi.png"
plt.savefig(filename, format='png', bbox_inches = 'tight', dpi = 600)

  