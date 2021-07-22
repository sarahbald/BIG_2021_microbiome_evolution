#!/bin/bash
#$ -N postproc_midas
#$ -e /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/scripts/postproc_error
#$ -o /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/scripts/postproc_output
#$ -cwd
#$ -r y
#$ -j y
#$ -l time=120:00:00
#$ -l h_data=34G
#$ -l highp

####ttt##$-tc 100 # Throttle to max 100 tasks at a time
########$-t 1-158
#######$-ltime=23:00:00


. /u/local/Modules/default/init/modules.sh



module unload python
module load anaconda/python2-4.2

source activate midas

module load samtools


module unload gcc
module load gcc/9.3.0



export PYTHONPATH=$PYTHONPATH:/u/project/ngarud/Garud_lab/MIDAS_mod
export PATH=$PATH:/u/project/ngarud/Garud_lab/MIDAS_mod/scripts
export MIDAS_DB=/u/project/ngarud/Garud_lab/midas_db_v1.2


###########################
# Postprocessing of MIDAS #
###########################


#bzip2 /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/species/coverage.txt &
#bzip2 /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/species/count_reads.txt &
#bzip2 /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/species/relative_abundance.txt &

# bzip all other species before running folowup

#for dir in /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/genes/*/; do
#    bzip2 ${dir}/genes_copynum.txt &
#    bzip2 ${dir}/genes_depth.txt &
#    bzip2 ${dir}/genes_presabs.txt &
#    bzip2 ${dir}/genes_reads.txt &
#done


# SNPs
#for dir in /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/snps/*/; do
#    echo ${dir}
#    bzip2 ${dir}/snps_depth.txt &
#    bzip2 ${dir}/snps_info.txt &
#    bzip2 ${dir}/snps_ref_freq.txt &
#    bzip2 ${dir}/snps_alt_allele.txt &
#done



# make species_snps.txt
#ls /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/snps/ > /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/snps/species_snps.txt

#ls /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/genes/ > /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/genes/species_genes.txt



# The purpose of this step is to and assign a probability tvhat a site's genootype is real or an error.

# calculate_marker_gene_coverage.py -- (another possibility is to change to using the CNV output)
    #output: marker_coverage.txt.bz2

# calculate_coverage_distribution.py -- This filters which genes have less than or greater than 2-fold difference from the median coverage in the data set. We want to filter these genes out (i.e. not enough info to call a SNP or too much coverage means that there might be a CNV)
    #output: sample_coverage_distribution.txt.bz2
    #output: sample_gene_coverage.txt.bz2

# calculate_error_pvalues.py -- (this uses the cpp code Ben wrote)
    # output: sample_annotated_snps.txt.bz2

#compile the error_pvalues cpp code as follows:
#g++ -std=c++11 -O3 /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/scripts/error_filtering_cpp_code/*.cpp -o /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/scripts/annotate_pvalue



# run core genes first:
#python /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/scripts/core_gene_utils.py

# python script to run the above:
#python postprocess_midas_data.py $species



# qsub script to run the above:
#python /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/scripts/print_good_species_list.py name> /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/tmp_intermediate_files/tmp_species_list.txt

#qsub /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/scripts/scripts/qsub_scripts/run_postprocessing_scripts


python /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/scripts/postprocess_all_midas_data_serial.py
