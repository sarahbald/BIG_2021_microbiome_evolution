import bz2
import pickle

data_directory = "C:/Users/sarah/Garud Lab/"

rel_abundance_file = bz2.BZ2File("C:/Users/sarah/Garud Lab/relative_abundance.txt.bz2" , "r")
samples = rel_abundance_file.readline() # header
samples = samples.decode().split('\t')
    
#get idxs for africa and HMP samples
africa_idx = []
north_america_idx = []
    
for i in range(len(samples)):
        
    if samples[i].startswith('S'):
        africa_idx.append(i)
        
    elif samples[i].startswith('7'):
        north_america_idx.append(i)
            
    else:
        continue
        
species_mean_abun_dict = {}
species_in_neither = []
species_prev_dict = {}
species_abun_dict = {}

#loop through all species
for line in rel_abundance_file:
#for i in range(1,20):
    
    line = line.strip()
    line = line.split('\t')
    species_name = line[0]

        
    africa_mean_abun = 0.0
    north_america_mean_abun = 0.0
    africa_prev = 0.0
    north_america_prev = 0.0
    
    species_abun_dict[species_name] = {}
    species_abun_dict[species_name]['Africa']=[]
    species_abun_dict[species_name]['North America']=[]
        
    for idx in africa_idx:
        africa_mean_abun += float(line[idx])
        if float(line[idx]) > 0.0:
            africa_prev += 1.0   
        species_abun_dict[species_name]['Africa'].append(line[idx])
            
    africa_mean_abun = africa_mean_abun/len(africa_idx)
    africa_prev = africa_prev/len(africa_idx)
    for idx in north_america_idx:
        north_america_mean_abun += float(line[idx])
        if float(line[idx]) > 0.0:
            north_america_prev += 1.0
        if species_abun_dict.has_key(species_name):
            species_abun_dict[species_name]['North America'].append(line[idx])
    north_america_mean_abun = north_america_mean_abun/len(north_america_idx)
    north_america_prev = north_america_prev/len(north_america_idx)
    
    if africa_mean_abun == 0.0 or north_america_mean_abun == 0.0:
        species_abun_dict.pop(species_name)     
        continue
    
    species_mean_abun_dict[species_name] = {}
    species_mean_abun_dict[species_name]['Africa'] = africa_mean_abun
    species_mean_abun_dict[species_name]['North America'] = north_america_mean_abun
    
    species_prev_dict[species_name] = {}
    species_prev_dict[species_name]['Africa'] = africa_prev
    species_prev_dict[species_name]['North America'] = north_america_prev


#remove species with no North America Abundances from species_abun_dict
for species in species_abun_dict.keys():
    if species_abun_dict[species]['North America']==[] or species_abun_dict[species]['Africa'] == []:
        species_abun_dict.pop(species)
    
intermediate_filename_template = data_directory+"species_mean_relative_abundance.dat"
intermediate_filename = intermediate_filename_template

with open(intermediate_filename, 'wb') as handle:
    pickle.dump(species_mean_abun_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
intermediate_filename_template = data_directory+"species_prevalence.dat"
intermediate_filename = intermediate_filename_template

with open(intermediate_filename, 'wb') as handle:
    pickle.dump(species_prev_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
intermediate_filename_template = data_directory+"species_abun_dict.dat"
intermediate_filename = intermediate_filename_template

with open(intermediate_filename, 'wb') as handle:
    pickle.dump(species_abun_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)


def load_species_mean_abun_dict():

    intermediate_filename = data_directory+"species_mean_relative_abundance.dat"

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)
    return b

def load_species_prev_dict():

    intermediate_filename = data_directory+"species_prevalence.dat"

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)
    return b

def load_species_abun_dict():

    intermediate_filename = data_directory+"species_abun_dict.dat"

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)
    return b
