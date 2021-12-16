import pickle
import matplotlib.pyplot as plt
import figure_utils
from operator import itemgetter
from scipy.stats import sem
import seaborn
import numpy as np

filename = 'C:\Users\sarah\Garud Lab\BIG_2021_microbiome_evolution\strain_dict.pickle'
with open(filename, 'rb') as handle:
    strain_dict = pickle.load(handle)

good_species = []    
for species_name in strain_dict.keys():
    if (len(strain_dict[species_name]['Africa'].keys()) >= 20):
        if (len(strain_dict[species_name]['North America'].keys()) >= 20):
            good_species.append(species_name)
            
#clean up dict, random 0.0 freq
for species in good_species:
    for sample in strain_dict[species]['Africa']:
        for freq in strain_dict[species]['Africa'][sample]:
            if float(freq) == 0.0:
                strain_dict[species]['Africa'][sample].remove(freq)
    for sample in strain_dict[species]['North America']:
        for freq in strain_dict[species]['North America'][sample]:
            if float(freq) == 0.0:
                strain_dict[species]['North America'][sample].remove(freq)
            
species_color_dict = {}
colors = seaborn.color_palette('muted', n_colors=len(good_species))
i = 0
good_species.sort()
for species in good_species:
    species_color_dict[species] = colors[i]
    i+=1
    

def plot_strain_samples(strain_dict):
    num_samples_na = {}
    num_samples_af = {}
    for species_name in strain_dict.keys():
        num_samples_af[species_name] = len(strain_dict[species_name]['Africa'])
        num_samples_na[species_name] = len(strain_dict[species_name]['North America'])
    sorted_num_samples_na = sorted(num_samples_na.items(), key=itemgetter(1))
    species_names = [species[0] for species in sorted_num_samples_na]
    plot_num_samples_na = [species[1] for species in sorted_num_samples_na]
    plot_num_samples_af = []
    for species in species_names:
        plot_num_samples_af.append(num_samples_af[species])
    
    plt.figure(figsize=(12,8))
    plt.scatter(species_names, plot_num_samples_na, label = 'North America')
    plt.scatter(species_names, plot_num_samples_af, label = 'Africa') 
    plt.axhline(y=20, linestyle='--', color = 'gray')
    plt.title("Distribution of Samples for Estimating Strain Structure")
    plt.ylabel("Number of Samples")
    plt.xlabel("Ranked Species")
    plt.text(65, 22, 'minimum samples', color = 'grey')
    plt.legend()
    plt.xticks([])
    plt.savefig('C:/Users/sarah/Garud Lab/plots/num_strain_samples.png', dpi=600)
    
    return

def plot_mean_strains(strain_dict, good_species, plot_evenness = False):
    mean_strains_na = {}
    mean_strains_af = {}
    std_err_na = {}
    std_err_af = {}
    if plot_evenness == True:
        mean_evenness_na = {}
        mean_evenness_af = {}
        std_err_evenness_na = {}
        std_err_evenness_af = {}
    for species_name in strain_dict.keys():
        num_strains_na = []
        num_strains_af = []
        if plot_evenness == True:
            evenness_af = []
            evenness_na = []            
        if species_name in good_species:
            for sample in strain_dict[species_name]['Africa'].keys():
                num_strains_af.append(len(strain_dict[species_name]['Africa'][sample]))
                if plot_evenness == True:
                    shan_diversity = 0.0
                    for freq in strain_dict[species_name]['Africa'][sample]:
                        shan_diversity += -1*(float(freq)*np.log(float(freq)))
                    evenness_af.append(shan_diversity/np.log(len(strain_dict[species_name]['Africa'][sample])))
            for sample in strain_dict[species_name]['North America'].keys():
                num_strains_na.append(len(strain_dict[species_name]['North America'][sample]))
                if plot_evenness == True:
                    shan_diversity = 0.0
                    for freq in strain_dict[species_name]['North America'][sample]:
                        shan_diversity += -1*(float(freq)*np.log(float(freq)))
                    evenness_na.append(shan_diversity/np.log(len(strain_dict[species_name]['North America'][sample])))
            mean_strains_af[species_name] = float(sum(num_strains_af))/float(len(num_strains_af))
            std_err_af[species_name] = sem(num_strains_af)
            mean_strains_na[species_name] = float(sum(num_strains_na))/float(len(num_strains_na))
            std_err_na[species_name] = sem(num_strains_na)
            if plot_evenness == True:
                mean_evenness_af[species_name] = float(sum(evenness_af))/float(len(evenness_af))
                std_err_evenness_af[species_name] = sem(evenness_af)
                mean_evenness_na[species_name] = float(sum(evenness_na))/float(len(evenness_na))
                std_err_evenness_na[species_name] = sem(evenness_na)
    sorted_mean_strains_na = sorted(mean_strains_na.items(), key=itemgetter(1))
    species_names = [species[0] for species in sorted_mean_strains_na]
    if plot_evenness == True:
        plot_mean_evenness_na = []
        plot_mean_evenness_af = []
        plot_std_err_evenness_na = []
        plot_std_err_evenness_af = []
        legend = []
        for species in species_names:
            plot_mean_evenness_af.append(mean_evenness_af[species])
            plot_mean_evenness_na.append(mean_evenness_na[species])
            plot_std_err_evenness_na.append(std_err_evenness_na[species])
            plot_std_err_evenness_af.append(std_err_evenness_af[species])
            legend.append(species)
        for i in range(len(legend)):
            if legend[i].startswith('Faecalibacterium_prau'):
                legend[i] = figure_utils.get_pretty_species_name(legend[i], include_number=True)
            else:
                legend[i] = figure_utils.get_pretty_species_name(legend[i])
        plt.clf()    
        plt.figure()  
        for i in range(len(species_names)):
            plt.scatter(plot_mean_evenness_af[i], plot_mean_evenness_na[i], color=species_color_dict[species_names[i]], label = legend[i])
        plt.errorbar(plot_mean_evenness_af, plot_mean_evenness_na, yerr=plot_std_err_evenness_na, xerr=plot_std_err_evenness_af, ecolor = [species_color_dict[i] for i in species_names], ls='none')
        plt.plot(np.linspace(.95,.99,10), np.linspace(.94,.99,10), ls='--', color='grey')
        plt.legend(bbox_to_anchor = (1.02,1), loc='upper left')
        plt.title("Strain Evenness")
        plt.ylabel("Mean Strain Evenness in Western Hosts")
        plt.xlabel("Mean Strain Evenness in non-Western Hosts")
        plt.savefig('C:/Users/sarah/Garud Lab/plots/strain_evenness.png', dpi=600, bbox_inches = 'tight')    
    else:        
        plot_mean_strains_na = [species[1] for species in sorted_mean_strains_na]
        plot_mean_strains_af = []
        plot_std_err_na = []
        plot_std_err_af = []
        legend = []
        for species in species_names:
            plot_mean_strains_af.append(mean_strains_af[species])
            plot_std_err_na.append(std_err_na[species])
            plot_std_err_af.append(std_err_af[species])
            legend.append(species)
        for i in range(len(legend)):
            if legend[i].startswith('Faecalibacterium_prau'):
                legend[i] = figure_utils.get_pretty_species_name(legend[i], include_number=True)
            else:
                legend[i] = figure_utils.get_pretty_species_name(legend[i])
        plt.clf()    
        plt.figure()  
        for i in range(len(species_names)):
            plt.scatter(plot_mean_strains_af[i], plot_mean_strains_na[i], color=species_color_dict[species_names[i]], label = legend[i])
        plt.errorbar(plot_mean_strains_af, plot_mean_strains_na, yerr=plot_std_err_na, xerr=plot_std_err_af, ecolor = [species_color_dict[i] for i in species_names], ls='none')
        plt.plot(np.linspace(2.4,3.2,10), np.linspace(2.4,3.2,10), ls='--', color='grey')
        plt.legend(bbox_to_anchor = (1.05,1), loc='upper left')
        plt.title("Strain Richness")
        plt.ylabel("Mean Richness in Western Hosts")
        plt.xlabel("Mean Richness in non-Western Hosts")
        plt.savefig('C:/Users/sarah/Garud Lab/plots/strain_richness.png', dpi=600, bbox_inches = 'tight')    
       
    return

def plot_num_strain_dist(strain_dict, good_species):
    fig, axs = plt.subplots(1, 2, sharex=True, sharey=True)
    af_num_strain_dict = {}
    na_num_strain_dict = {}
    for species in good_species:
        af_num_strains = []
        na_num_strains = []
        for sample in strain_dict[species]['Africa'].keys():
            af_num_strains.append(len(strain_dict[species]['Africa'][sample]))
        for sample in strain_dict[species]['North America'].keys():
            na_num_strains.append(len(strain_dict[species]['North America'][sample]))
        af_num_strain_dict[species] = np.sort(af_num_strains)
        na_num_strain_dict[species] = np.sort(na_num_strains)
    for species in good_species:
        na_survival_plot = 1 - np.arange(len(na_num_strain_dict[species]))/float(len(na_num_strain_dict[species]))
        af_survival_plot = 1 - np.arange(len(af_num_strain_dict[species]))/float(len(af_num_strain_dict[species]))
        axs[0].plot(na_num_strain_dict[species], na_survival_plot)
        axs[1].plot(af_num_strain_dict[species], af_survival_plot)
    axs[0].set_title('Non-Western Hosts') 
    axs[1].set_title('Western Hosts')
    axs[0].set_xlabel('Number of Strains')
    axs[1].set_xlabel('Number of Strains')
    axs[0].set_ylabel('Proportion of Hosts')
    legend = good_species
    for i in range(len(legend)):
        if legend[i].startswith('Faecalibacterium_prau'):
            legend[i] = figure_utils.get_pretty_species_name(legend[i], include_number=True)
        else:
            legend[i] = figure_utils.get_pretty_species_name(legend[i])
    plt.legend(legend, bbox_to_anchor = (1.02,1), loc='upper left')
    plt.savefig('C:/Users/sarah/Garud Lab/plots/strain_dist.png', dpi=600)
    return

#plot_num_strain_dist(strain_dict, good_species)
#plot_strain_samples(strain_dict)
#plot_mean_strains(strain_dict, good_species, plot_evenness = True)
plot_num_strain_dist(strain_dict, good_species)