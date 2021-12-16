import bz2
import sys
import os

species_names_file =  "/u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/snps/species_snps.txt"
species_names = open(species_names_file, 'r')
names = []
for line in species_names:
    line = line.strip('\n')
    names.append(line)


snp_directory = "/u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/snps"

for species_name in names:
    if os.path.isfile(os.path.join(snp_directory, species_name, "annotated_snps.txt.bz2")):
        snps_file = bz2.BZ2File(snp_directory+"/"+species_name+"/annotated_snps.txt.bz2" , "r")
        samples = snps_file.readline() # header
        samples = samples.decode().split('\t')
        output_directory = "/u/project/ngarud/sarahbal/BIG_2021_microbiome_evolution/scripts/sample_names/%s.txt" % species_name
        f = open(output_directory, 'w')
        for i in range(1,len(samples)):
            if samples[i].endswith('c'):
                f.write(samples[i][:-1]+'\n') 
            else:
                f.write(samples[i]+'\n')
        f.close()



