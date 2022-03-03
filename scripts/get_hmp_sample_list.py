import os
import sample_utils
import parse_HMP_data

stem = "/u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/midas_output/"
hmp_sample_names = os.listdir(stem)

subject_sample_map = parse_HMP_data.parse_subject_sample_map()
samples_all_unique = hmp_sample_names[sample_utils.calculate_unique_samples(subject_sample_map, sample_list=hmp_sample_names)]

for sample in samples_all_unique:
    with open("/u/project/ngarud/sarahbal/BIG_2021_microbiome_evolution/scripts/hmp_sample_names.txt", "wb") as f:
        f.write(str(sample)+ '\n')
        f.close