import parse_HMP_data
#import parse_africa_data

import sample_utils
import diversity_utils




subject_sample_map = parse_HMP_data.parse_subject_sample_map()

good_species_list = ['Bacteroides_vulgatus_57955']

#for species_name in good_species_list:

#    snp_samples = diversity_utils.calculate_haploid_samples(species_name)

#    snp_samples = snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]

#    print(snp_samples)



#Traceback (most recent call last):
 # File "/u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/scripts/calculate_singletons.py", line 163, in <module>
#    snp_samples = snp_samples[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]
 # File "/u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/scripts/sample_utils.py", line 550, in calculate_unique_samples
#    subject = sample_subject_map[sample]
#KeyError: '700016456'


# second item is sample name
#'700175181': ('764851951', '700175181', 'SRS149879', 'United States', 'North America', 1L)
#'V50_02_1FE': ('PRJNA485056_subject_47', 'V50_02_1FE', 'SRR7658625', 'Madagascar', 'Africa', 1)

#replace third with second for africa



def flatten_samples(subject_sample_map):

    grouping_replicate_map = {}
    for subject in sorted(subject_sample_map.keys()):
        for sample in sorted(subject_sample_map[subject].keys()):
            grouping_replicate_map[sample] = subject_sample_map[subject][sample]

    return grouping_replicate_map

def calculate_subject_pairs(subject_sample_map, sample_list=[]):

    sample_continent_map = parse_HMP_data.parse_sample_continent_map()

    if len(sample_list)==0:
        sample_list = list(sorted(flatten_samples(subject_sample_map).keys()))

    new_sample_list = []
    for sample in sample_list:
        if sample.endswith('c'):
            new_sample_list.append(sample[:-1])
        else:
            new_sample_list.append(sample)

    sample_list = new_sample_list

    # invert subject sample map
    sample_subject_map = {}
    for subject in subject_sample_map.keys():
        for sample in subject_sample_map[subject].keys():
            sample_subject_map[sample] = subject

    #sample_subject_map = {}
    #for subject in subject_sample_map.keys():
    #    for sample in list(subject_sample_map[subject].values()[0]):
    #        sample_subject_map[sample] = subject

    #print(sample_subject_map.keys())


    #if '700161629' in sample_subject_map:
    #    print(True)
    #else:
    #    print(False)


#calculate_subject_pairs(subject_sample_map_hmp)


#same_sample_idxs, same_subject_idxs, diff_subject_idxs = sample_utils.calculate_subject_pairs(subject_sample_map_hmp)

#print(diff_subject_idxs)


#print(sample_continent_map)

#if '700161629' in sample_continent_map:
#    print(True)

#else:
#    print(False)





x = sample_utils.parse_sample_metadata_map()

print(x)
