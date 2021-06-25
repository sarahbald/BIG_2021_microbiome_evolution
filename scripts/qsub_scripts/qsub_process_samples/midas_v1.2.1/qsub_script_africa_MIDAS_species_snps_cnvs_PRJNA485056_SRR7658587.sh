#!/bin/bash
#$ -N MIDAS_6
#$ -e /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/postproc_error
#$ -o /u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/postproc_output
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

#conda create -n midas python=2.7
#conda install biopython
#conda install numpy

source activate midas

module load samtools



export PYTHONPATH=$PYTHONPATH:/u/project/ngarud/Garud_lab/MIDAS_mod
export PATH=$PATH:/u/project/ngarud/Garud_lab/MIDAS_mod/scripts
export MIDAS_DB=/u/project/ngarud/Garud_lab/midas_db_v1.2




acc=SRR7658587



echo $acc

OUTDIR=/u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/midas_output_v1.2.1/${acc}_indiv

mkdir -p $OUTDIR

fastq1=/u/project/ngarud/wrshoema/madagascar_data/data/SRA_files/PRJNA485056/${acc}_1.fastq.gz
fastq2=/u/project/ngarud/wrshoema/madagascar_data/data/SRA_files/PRJNA485056/${acc}_2.fastq.gz
run_midas.py species $OUTDIR -1 $fastq1 -2 $fastq2 --remove_temp
run_midas.py genes $OUTDIR -1 $fastq1 -2 $fastq2 --remove_temp
run_midas.py snps $OUTDIR -1 $fastq1 -2 $fastq2 --remove_temp
