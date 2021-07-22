import os
import sys


rootdir_africa = '/u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/midas_output_v1.2.1'
rootdir_hmp = '/u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/midas_output'


file = open('/u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/scripts/midas_output_paths.txt', 'w')



for sample in os.listdir(rootdir_africa):

    if sample == 'nohup.out':
        continue

    out_line = '%s/%s\n' % (rootdir_africa, sample)
    file.write(out_line)


for sample in os.listdir(rootdir_hmp):

    if sample == 'nohup.out':
        continue

    out_line = '%s/%s\n' % (rootdir_hmp, sample)
    file.write(out_line)


file.close()
