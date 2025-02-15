###############################################################################
#
# Set up default source and output directories
#
###############################################################################
import os.path
import os
from math import log10

if os.geteuid() == 501:

    data_directory = os.path.expanduser("~/GitHub/BIG_2021_microbiome_evolution/data/")
    analysis_directory = os.path.expanduser("~/GitHub/BIG_2021_microbiome_evolution/analysis/")
    scripts_directory = os.path.expanduser("~/GitHub/BIG_2021_microbiome_evolution/scripts/")

    midas_directory = os.path.expanduser("~/GitHub/BIG_2021_microbiome_evolution/data/midas_db_v1.2/")

    #patric_directory = os.path.expanduser("~/GitHub/BIG_2021_microbiome_evolution/data/patric_db/")
    patric_directory = os.path.expanduser("~/GitHub/BIG_2021_microbiome_evolution/data/patric_db/")

    metadata_directory = os.path.expanduser("~/GitHub/BIG_2021_microbiome_evolution/scripts/metadata/")



else:

    data_directory = os.path.expanduser("/u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/")
    analysis_directory = os.path.expanduser("/u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/analysis/")
    scripts_directory = os.path.expanduser("/u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/scripts/")

    midas_directory = os.path.expanduser("/u/project/ngarud/Garud_lab/midas_db_v1.2/")
    #patric_directory = os.path.expanduser("/u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/data/patric_db/")
    #patric_directory = os.path.expanduser("/u/project/ngarud/Garud_lab/PATRIC/")
    patric_directory = os.path.expanduser("/u/project/ngarud/Garud_lab/software/PATRIC/")
    

    metadata_directory = os.path.expanduser("/u/project/ngarud/Garud_lab/BIG_2021_microbiome_evolution/scripts/metadata/")

#patric_directory = os.path.expanduser("~/patric_db/")
#patric_directory = os.path.expanduser("~/patric_db/")

#midas_directory = os.path.expanduser("/u/project/ngarud/wrshoema/negative_selection_microbiome/data/midas_db_v1.2/")

# We use this one to debug because it was the first one we looked at
debug_species_name = 'Bacteroides_uniformis_57318'

good_species_min_coverage = 10
good_species_min_prevalence = 10

min_median_coverage = 20

consensus_lower_threshold = 0.2
consensus_upper_threshold = 0.8
fixation_min_change = consensus_upper_threshold-consensus_lower_threshold
fixation_log10_depth_ratio_threshold = log10(3)

threshold_within_between_fraction = 0.1
threshold_pi = 1e-03

min_opportunities = 100000

modification_difference_threshold = 20
replacement_difference_threshold = 500

twin_modification_difference_threshold = 1000
twin_replacement_difference_threshold = 1000

gainloss_max_absent_copynum = 0.05
gainloss_min_normal_copynum = 0.6
gainloss_max_normal_copynum = 1.2

core_genome_min_copynum = 0.3
core_genome_max_copynum = 3 # BG: should we use a maximum for "core genome"? I'm going to go w/ yes for now
core_genome_min_prevalence = 0.9
shared_genome_min_copynum = 3

# Default parameters for pipe snps
# (Initial filtering for snps, done during postprocessing)
pipe_snps_min_samples=4
pipe_snps_min_nonzero_median_coverage=5
pipe_snps_lower_depth_factor=0.3
pipe_snps_upper_depth_factor=3

parse_snps_min_freq = 0.05

between_host_min_sample_size = 33
between_host_ld_min_sample_size = 10
within_host_min_sample_size = 3
within_host_min_haploid_sample_size = 10

between_low_divergence_threshold = 2e-04



# Comment this out
#from parse_HMP_data import *
# and uncomment this
#from parse_simulated_data import *
# for isolate data
