#!/bin/bash
#$ -N gene_prev
#$ -e /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/scripts/gene_prev_error
#$ -o /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/scripts/gene_prev_output
#$ -cwd
#$ -r y
#$ -j y
#$ -l time=24:00:00
#$ -l h_data=12G
#$ -l highp

####ttt##$-tc 100 # Throttle to max 100 tasks at a time
########$-t 1-158
#######$-ltime=23:00:00


. /u/local/Modules/default/init/modules.sh



module unload python
module load anaconda/python2-4.2

source activate midas


python /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/scripts/calculate_delta_prevalence.py
