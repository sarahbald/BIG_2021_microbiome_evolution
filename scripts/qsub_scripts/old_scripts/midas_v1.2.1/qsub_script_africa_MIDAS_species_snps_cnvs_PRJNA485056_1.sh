#!/bin/bash
#$ -N MIDAS_1
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

#module unload python
#module load python/2.7
#module load python/2.7.13




module unload python
module load anaconda/python2-4.2

#conda create -n midas python=2.7
#conda install biopython
#conda install numpy

source activate midas

module load samtools
module load bowtie2



#export PYTHONPATH=$PYTHONPATH:/u/project/ngarud/Garud_lab/MIDAS
#export PATH=$PATH:/u/project/ngarud/Garud_lab/MIDAS/scripts
export PYTHONPATH=$PYTHONPATH:/u/project/ngarud/Garud_lab/MIDAS_mod
export PATH=$PATH:/u/project/ngarud/Garud_lab/MIDAS_mod/scripts
export MIDAS_DB=/u/project/ngarud/Garud_lab/midas_db_v1.2

#export LD_LIBRARY_PATH=../lib:${LD_LIBRARY_PATH}


samples_PRJNA485056=()

{
read
while IFS= read -r acc; do
    samples_PRJNA485056+=($acc)

done
} < /u/project/ngarud/wrshoema/negative_selection_microbiome/data/madagascar/PRJNA485056_run_accessions_1.txt



for acc in "${samples_PRJNA485056[@]}"
do
  echo $acc
  OUTDIR=/u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/midas_output_v1.2.1/${acc}

  mkdir -p $OUTDIR

  fastq1=/u/project/ngarud/wrshoema/madagascar_data/data/SRA_files/PRJNA485056/${acc}_1.fastq.gz
  fastq2=/u/project/ngarud/wrshoema/madagascar_data/data/SRA_files/PRJNA485056/${acc}_2.fastq.gz
  run_midas.py species $OUTDIR -1 $fastq1 -2 $fastq2 --remove_temp
  run_midas.py genes $OUTDIR -1 $fastq1 -2 $fastq2 --remove_temp
  run_midas.py snps $OUTDIR -1 $fastq1 -2 $fastq2 --remove_temp
done
