# BIG_2021_microbiome_evolution

Repository for code used by the 2021 BIG summer students in the Garud lab. This code builds off of the [microbiome_evolution](https://github.com/benjaminhgood/microbiome_evolution) repo and retains the same license.


### Setting up your working environment

For big students:

The code in this repo is primarily written in Python 2.7 and uses the following Python packages: numpy, matplotlib, and scipy. If you do not have Python 2.7 on your local machine, you should create an environment with 2.7 using conda. First, create the environment

`conda create -n microbiome_evolution python=2.7`

Then, enter the environment and install the necessary packages

`source activate microbiome_evolution`

`conda install numpy`

`conda install matplotlib`

`conda install scipy`

Then you will have to create a local copy of this repo. First, make sure you have a GitHub account. Then, fork `BIG_2021_microbiome_evolution` by clicking the "Fork" button and then clicking your account. Next, you'll need to create a copy of your forked repo on your computer. We can go over the intricacies of GitHub later in the summer, but for now the instructions to set up a GitHub repo on your computer can be found [in this PDF](https://github.com/QuantitativeBiodiversity/QB-2019/blob/master/1.HandOuts/2.Github/2.GitHandout.pdf )


### Microbiome data


The data is on Hoffman2 under the following directory

`/u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution`

It has been processed using the pipeline from Garud et al., 2019. If you want to see what commands were run for each data type, look into the commands that are called in `./scripts/postprocess_midas_data.py`


A brief overview of the data types

In the directory `genes` each species has the following files

- `genes_copynum.txt.bz2`: Copy numbers of each gene in each host

- `genes_depth.txt.bz2`: Depth of sequencing coverage at each gene in each sample

- `genes_presabs.txt.bz2`: Whether a gene is present or absent in a given host

- `genes_reads.txt.bz2`: Total number of reads that mapped to a gene in a given host


The directory `core_genes` contains `*_gene_freqs.txt.gz` files, which lists the proportion of hosts where a given gene is present

The directory `snps` contains the following files

- `snps_alt_allele.txt.bz2`: The alternative allelic state for each site.

- `snps_depth.txt.bz2`: The depth of coverage at a given site where a SNP was called

- `snps_info.txt.bz2`: A lot of information. The file contains the mean frequency, mean depth of coverage, the proportion of times where a site was found, the allelic state, the type of site (1D, 4D, etc), gene name, and amino acid.

- `snps_ref_freq.txt.bz2`: The frequency of the reference allele at a given site for each host.

- `snps_summary.txt`: The size of the genome, number of covered bases, fraction covered, and mean coverage for each sample

Each species has a file in `snp_prevalences` containing the frequency of the alternative allele and the frequency of a SNP for each site for each host.

Files in `singleton_rates` are the result of comparing pairs of hosts to determine what mutations are **unique** to a given host. It includes information on the type of sites that are compared (e.g., sites in core genes, 4D sites, etc), the number of mutations, number of reversions and number of opportunities for a mutation or reversion



The directory `snv_distances` contains the output of genetic distances calculated from phased haplotypes for each species. For each site within a species, the derived allele count, ancestral allele count, and minimum and maximum genetic distance are calculated.

The directory `private_snvs` contains information about what SNVs are unique to a given host.

The directory `substitution_rates` contains information about substitution rates calculated by comparing a given pair of samples. It has a similar structure to `singleton_rates`

`linkage_disequilibria` contains linkage disequilibria calculations for each species. 








- `core_genes`:   
