#!/bin/bash

#$ -N all_sp_sam
#$ -e /u/project/ngarud/sarahbal/BIG_2021_microbiome_evolution/strainfinder_input/$TASK_ID
#$ -o /u/project/ngarud/sarahbal/BIG_2021_microbiome_evolution/strainfinder_input/$TASK_ID
#$ -cwd
#$ -r y
#$ -j y
#$ -l time=72:00:00
#$ -l h_data=12G
#$ -l highp
#$ -t 1-249



#############ONLY CHANGED TEMPORARILY TO COMPARE W RICKYS DATA
SAMPLE=$(sed "${SGE_TASK_ID}q;d" "hmp_sample_names.txt")
#SAMPLE=$(sed "${SGE_TASK_ID}q;d" "sample_names/${SPECIES}.txt")
#DIR_PATH="/u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/midas_output_v1.2.1/"
DIR_PATH="/u/project/ngarud/Garud_lab/metagenomic_fastq_files/HMP1-2/midas_output/"

echo ${SAMPLE}
echo $DIR_PATH
#list of files in form {SPECIES}.snps.gz
SPECIES_LIST=$(ls $DIR_PATH${SAMPLE}/snps/output/)

echo ${SPECIES_LIST}
. /u/local/Modules/default/init/modules.sh

module unload python
module load python/2.7.15 ##--enable-unicode=ucs4
##./configure --enable-unicode=ucs4
####source activate midas
pip uninstall -y scipy
pip install scipy
pip install openopt
pip install FuncDesigner
pip install DerApproximator


for species in $SPECIES_LIST
do
    SPECIES=$species
    
    echo $SPECIES

    python /u/project/ngarud/sarahbal/BIG_2021_microbiome_evolution/scripts/test_ricky_get_strainfinder_input.py  ${SPECIES}  ${SAMPLE}  /u/project/ngarud/sarahbal/BIG_2021_microbiome_evolution/strainfinder_input/

    python /u/project/ngarud/sarahbal/BIG_2021_microbiome_evolution/scripts/test_ricky_strainfinder.py  ${SPECIES}  ${SAMPLE} 

    python /u/project/ngarud/sarahbal/BIG_2021_microbiome_evolution/scripts/test_ricky_postprocess_strainfinder_output.py  ${SPECIES}  ${SAMPLE}
done



