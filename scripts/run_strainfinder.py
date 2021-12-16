import sys
import os

species_name = sys.argv[1]
sample_name = sys.argv[2]
Nmax = 4
Lmax = 20000


output_directory = os.path.expanduser("/u/project/ngarud/sarahbal/BIG_2021_microbiome_evolution/strainfinder_output/")

####################

## FOR SARAH BALD

## You'll need to modify this script so the input/output not only contains what species to run strainfinder on, but also what host/sample to run it on. 

###################

filename_prefix = "%s%s/%s" % (output_directory, species_name, sample_name)

os.system('python create_StrainFinderInput.py -o %s --species %s --sample %s -L %d' % (output_directory, species_name, sample_name, Lmax))
#sys.exit(0)
for N in xrange(1,Nmax+1):

    alignment_filename = filename_prefix+".strainfinder.p" 
    em_filename = filename_prefix+(".%d.em.cpickle" % N)
    log_filename = filename_prefix+(".%d.log.txt" % N)
    otu_filename = filename_prefix+(".%d.otu.txt" % N)

    os.system('python StrainFinder.py --aln %s -N %d --max_reps 10 --dtol 1 --ntol 2 --max_time 20000 --converge --em %s --em_out %s --otu_out %s --log %s --n_keep 2 --force_update --merge_out --msg' % (alignment_filename,N,em_filename, em_filename,otu_filename,log_filename))
    
    
