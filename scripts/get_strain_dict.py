import os
import pickle

strain_dict = {}

strainfinder_dir = "/u/project/ngarud/sarahbal/BIG_2021_microbiome_evolution/strainfinder_output/"
species = os.listdir(strainfinder_dir)

for name in species:

    strain_dict[name]={}
    strain_dict[name]['Africa'] = {}
    strain_dict[name]['North America'] = {}
    output_file = strainfinder_dir + "%s/final_strain_freq.txt" % name
#output_file = "C:/Users/sarah/Garud Lab/final_strain_freq.txt"
    with open(output_file) as f:
        for line in f:
            line = line.strip()
            line = line.split(' ')
            if line[0].startswith('B'):
                freq = []
                freq.append(line[10].strip('['))
                for i in range(11,len(line)):
                    if line[i].strip(']') == '':
                        continue
                    freq.append(line[i].strip(']'))
                if line[6].startswith('7'):
                    strain_dict[name]['North America'][line[6].strip(':')] = freq
                if line[6].startswith('S'):
                    strain_dict[name]['Africa'][line[6].strip(':')] = freq    

pickle.dump(strain_dict, strainfinder_dir + "strain_dict.pickle")
                    