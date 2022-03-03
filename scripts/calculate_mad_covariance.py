import numpy as np
import pickle

data_directory = "C:/Users/sarah/Garud Lab/"

filename = data_directory+"species_abun_dict.dat"
with open(filename, 'rb') as handle:
        species_abun_dict_zeros = pickle.load(handle)
na_mean_species_abun_zeros = []
af_mean_species_abun_zeros = []
for species in species_abun_dict_zeros.keys():
    mean_abun_af = [float(i) for i in species_abun_dict_zeros[species]['Africa']]
    mean_abun_af = np.mean(mean_abun_af)
    af_mean_species_abun_zeros.append(mean_abun_af)
    mean_abun_na = [float(i) for i in species_abun_dict_zeros[species]['North America']]
    mean_abun_na = np.mean(mean_abun_na)
    na_mean_species_abun_zeros.append(mean_abun_na)
    
cov_mat = np.cov(na_mean_species_abun_zeros, af_mean_species_abun_zeros)


