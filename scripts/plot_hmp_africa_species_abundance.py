import calculate_species_relative_abundance
import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats


data_directory = "C:/Users/sarah/Garud Lab/"

prev_dict = calculate_species_relative_abundance.load_species_prev_dict()
africa_prev_dict = {}
hmp_prev_dict = {}
for species_prev in prev_dict.keys():
    africa_prev_dict[species_prev] = prev_dict[species_prev]['Africa']
    hmp_prev_dict[species_prev] = prev_dict[species_prev]['North America']
    
africa_abundances = africa_prev_dict.values()
hmp_abundances = hmp_prev_dict.values()

species_names = [species for species in africa_prev_dict.keys()]

select_species = ["Eubacterium_rectale_56927", "Faecalibacterium_cf_62236", "Faecalibacterium_prausnitzii_57453", 
                  "Faecalibacterium_prausnitzii_61481", "Faecalibacterium_prausnitzii_62201", "Oscillibacter_sp_60799", 
                  "Prevotella_copri_61740", "Roseburia_inulinivorans_61943"]

col = []
for i in range(len(species_names)):
    if species_names[i] in select_species:
        col.append('red')
    else:
        col.append('blue')
    

plt.figure(figsize=(20,15))
for i in range(len(africa_abundances)):
    if col[i] == 'red':
        plt.scatter(africa_abundances[i], hmp_abundances[i], c = col[i], zorder=2)
    if col[i] == 'blue':
        plt.scatter(africa_abundances[i], hmp_abundances[i], c = col[i], zorder=1)

slope, intercept, r_value, p_value, se = stats.linregress(np.log10(africa_abundances), np.log10(hmp_abundances))
x_log10_range =  np.linspace(min(np.log10(africa_abundances)) , max(np.log10(africa_abundances)) , 10000)
y_log10_fit_range = slope*x_log10_range + intercept
plt.plot(10**x_log10_range, 10**y_log10_fit_range, c='k', lw=2, linestyle='-', zorder=2, label=str("Slope: "+str(slope)+" Intercept: "+str(intercept)+" R-value: "+str(r_value)))
plt.legend()
plt.title("", fontsize=30)
plt.ylabel("North America Species Relative Abundance", fontsize=20)
plt.yscale('log')
plt.yticks([10e-3,10e-2, 10e-1], fontsize=10)
plt.xlabel("Africa Species Relative Abundance", fontsize=20)
plt.xscale('log')
plt.xticks([10e-3,10e-2, 10e-1], fontsize=10)
plt.xlim([(min(africa_abundances) - .0005),(max(africa_abundances)+ .1)])
plt.ylim([(min(hmp_abundances) - .0005),(max(hmp_abundances)+ .1)])
#plt.margins(x=0,y=0)
plt.savefig('C:/Users/sarah/Garud Lab/plots/africa_hmp_species_abundance.png', bbox_inches='tight', dpi=600)