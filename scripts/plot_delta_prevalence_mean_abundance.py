import calculate_species_relative_abundance
import matplotlib.pyplot as plt


data_directory = "C:/Users/sarah/Garud Lab/"

prev_dict = calculate_species_relative_abundance.load_species_prev_dict()
delta_prev_dict = {}
for species_prev in prev_dict.keys():
    delta_prev = prev_dict[species_prev]['Africa'] - prev_dict[species_prev]['North America']
    delta_prev_dict[species_prev] = delta_prev

species_names = [species for species in delta_prev_dict.keys()]
delta_prevalences = [species for species in delta_prev_dict.values()]

mean_abundance_dict = calculate_species_relative_abundance.load_species_mean_abun_dict()
delta_abundance_dict = {}
for species_abundance in mean_abundance_dict.keys():
    delta_abun = mean_abundance_dict[species_abundance]['Africa'] - mean_abundance_dict[species_abundance]['North America']
    delta_abundance_dict[species_abundance] = delta_abun

abundances = [species for species in delta_abundance_dict.values()]

select_species = ["Eubacterium_rectale_56927", "Faecalibacterium_cf_62236", "Faecalibacterium_prausnitzii_57453", 
                  "Faecalibacterium_prausnitzii_61481", "Faecalibacterium_prausnitzii_62201", "Oscillibacter_sp_60799", 
                  "Prevotella_copri_61740", "Roseburia_inulinivorans_61943"]

col = []
for i in range(len(species_names)):
    if species_names[i] in select_species:
        col.append('red')
    else:
        col.append('blue')
    

plt.figure(figsize=(25,20))
for i in range(len(delta_prevalences)):
    if col[i] == 'red':
        plt.scatter(delta_prevalences[i], abundances[i], c = col[i], zorder=2)
    if col[i] == 'blue':
        plt.scatter(delta_prevalences[i], abundances[i], c = col[i], zorder=1)
    
plt.axvline(x=0, linestyle='--', color = 'gray')
plt.title("Difference in Species Prevalence\n Between African and North American Hosts", fontsize=50)
plt.ylabel("Difference in Mean Relative Abundance(Africa - NA)", fontsize=35)
plt.yscale('symlog', linthreshy = 1e-6)
plt.yticks([-10e-1, -10e-3,-10e-5,0,10e-5,10e-3, 10e-1], fontsize=25)
#plt.yscale('symlog', linthreshy = 1e-2)
#plt.yticks([-10e-1, -10e-2,-10e-3,0,10e-3,10e-2, 10e-1], fontsize=25)
plt.xlabel("Difference in Prevalence(Africa - NA)", fontsize=35)
#plt.xscale('symlog', linthreshx = 1e-5)
#plt.xticks([-10e-1, -10e-3,-10e-5,0,10e-5,10e-3, 10e-1], fontsize=25)
plt.xscale('symlog', linthreshx = 1e-2)
plt.xticks([-10e-1, -10e-2,-10e-3,0,10e-3,10e-2, 10e-1], fontsize=25)
plt.savefig('C:/Users/sarah/Garud Lab/plots/species_prevalence_abundance.png', dpi=600)