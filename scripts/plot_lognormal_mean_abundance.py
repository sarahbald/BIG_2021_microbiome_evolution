import calculate_species_relative_abundance
import numpy as np
import matplotlib.pyplot as plt

mean_abundance_dict = calculate_species_relative_abundance.load_species_mean_abun_dict()
africa_mean_abun = []
north_america_mean_abun = []
for species_abundance in mean_abundance_dict.keys():
    africa_mean_abun.append(np.log10(mean_abundance_dict[species_abundance]['Africa']))
    north_america_mean_abun.append(np.log10(mean_abundance_dict[species_abundance]['North America']))
    
africa_mean_abun = np.array(africa_mean_abun)
standardized_africa_abun = (africa_mean_abun - africa_mean_abun.mean())/africa_mean_abun.std()
north_america_mean_abun = np.array(north_america_mean_abun)
standardized_na_abun = (north_america_mean_abun - north_america_mean_abun.mean())/north_america_mean_abun.std()
    
plt.figure()
plt.hist(standardized_africa_abun, bins=30, alpha=0.5, label='Africa')
plt.hist(standardized_na_abun, bins=30, alpha=0.5, label='North America')
plt.xlim(-3,3)
plt.xlabel("Rescaled Log of Abundances")
plt.legend(loc='upper right')