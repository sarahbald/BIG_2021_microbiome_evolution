
import sys
import numpy
import config

import diversity_utils
import parse_midas_data
import parse_HMP_data
import calculate_substitution_rates


import matplotlib.pyplot as plt




from scipy.cluster.hierarchy import dendrogram, linkage
from scipy.cluster.hierarchy import cophenet
from scipy.cluster.hierarchy import fcluster



#good_species_list = ['Escherichia_coli_58110', 'Eubacterium_rectale_56927']


species_name = 'Eubacterium_rectale_56927'



# Load subject and sample metadata
sys.stderr.write("Loading sample metadata...\n")
subject_sample_map = parse_HMP_data.parse_subject_sample_map()
sample_continent_map = parse_HMP_data.parse_sample_continent_map()



sys.stderr.write("Loading haploid samples...\n")
    # Only plot samples above a certain depth threshold that are "haploids"
snp_samples = diversity_utils.calculate_haploid_samples(species_name)

sys.stderr.write("Calculating unique samples...\n")
# Only consider one sample per person
snp_samples = snp_samples[parse_midas_data.calculate_unique_samples(subject_sample_map, sample_list=snp_samples)]



# Load divergence matrices
sys.stderr.write("Loading pre-computed substitution rates for %s...\n" % species_name)
substitution_rate_map = calculate_substitution_rates.load_substitution_rate_map(species_name)
sys.stderr.write("Calculating matrix...\n")
dummy_samples, snp_difference_matrix, snp_opportunity_matrix = calculate_substitution_rates.calculate_matrices_from_substitution_rate_map(substitution_rate_map, 'core', allowed_samples=snp_samples)
snp_samples = dummy_samples
sys.stderr.write("Done!\n")

snp_substitution_matrix = snp_difference_matrix*1.0/(snp_opportunity_matrix+(snp_opportunity_matrix==0))





#print(snp_substitution_matrix)



#snp_substitution_rate = divergence_matrices[species_name]
snp_substitution_rate = snp_substitution_matrix
snp_substitution_rate = numpy.clip(snp_substitution_rate,1e-11,10)

sys.stderr.write("Calculating UPGMA dendrogram...\n")
# calculate compressed distance matrix suitable for agglomerative clustering
Y = []
for i in xrange(0,snp_substitution_rate.shape[0]):
    for j in xrange(i+1,snp_substitution_rate.shape[1]):
        Y.append(snp_substitution_rate[i,j])
Y = numpy.array(Y)
Z = linkage(Y, method='average')
c, coph_dists = cophenet(Z, Y)
ddata = dendrogram(Z, no_plot=True)
sys.stderr.write("Done! cophenetic correlation: %g\n" % c)


#################################################
#
# Plot dendrogram figure
#
#######

# calculate second minimum y value
ys = []
xs = []
for i, d in zip(ddata['icoord'], ddata['dcoord']):
    ys.extend(d)
    xs.extend(i)

xs = list(set(xs))
xs.sort()
xs = numpy.array(xs)

dx = xs[-1]-xs[0]
xmin = xs[0]-dx*0.025
xmax = xs[-1]+dx*0.025

ys = list(set(ys))
ys.sort()
ys = numpy.array(ys)

if ys[0]<1e-10:
    y_penultimin = ys[1]/2
else:
    y_penultimin = ys[0]/2

y_penultimax = ys[-1]

ymin = 1e-06
#ymin=2e-10
ymax=1e-01

yplotmin = 1e-06
yplotmax = 1e-01


leaf_xs = []


fig, dendrogram_axis = plt.subplots(figsize=(4,4))



for icoord, dcoord in zip(ddata['icoord'], ddata['dcoord']):
    for idx in xrange(0,len(icoord)-1):
        x0 = icoord[idx]
        y0 = dcoord[idx]
        if y0<1e-10:
            y0 = ymin
        x1 = icoord[idx+1]
        y1 = dcoord[idx+1]
        if y1<1e-10:
            y1 = ymin

        if (y0==ymin):
            leaf_xs.append(x0)

        if (y1==ymin):
            leaf_xs.append(x1)

        if (y0<2e-04) and (y1<2e-04):
            linewidth=0.75
            color='0.4'
        else:
            linewidth=0.3
            color='0.6'

        #print x0, '->', x1, '\t',y0, '->', y1
        print([x0,x1],[y0,y1])
        dendrogram_axis.semilogy([x0,x1],[y0,y1],'-',color=color,linewidth=linewidth)

        if (y0==y_penultimax) and (y1==y_penultimax):
            # it's the cross bar that bridges the two most-diverged clades
            # so plot a root branch to the top of the plot
            xavg = (x0+x1)*0.5


            dendrogram_axis.semilogy([xavg,xavg],[y_penultimax, ymax],'-',color=color,linewidth=linewidth)


#([635.0, 635.0], [1e-06, 1e-09])
#([635.0, 645.0], [1e-09, 1e-09])
#([645.0, 645.0], [1e-09, 1e-06])



leaf_xs = list(sorted(set(leaf_xs)))

xticks = []
xticklabels = []
samples = []

print species_name
#outFile.write(species_name +'\n')

for i in xrange(0,len(ddata['ivl'])):

    idx = long(ddata['ivl'][i])
    x = leaf_xs[i]
    y = yplotmin

    sample = snp_samples[idx]
    xticks.append(x)
    xticklabels.append(str(i))
    samples.append(sample)

    #print i, sample
    #outFile.write(str(i) + ' ' + sample +'\n')

    if sample_continent_map[sample]=='North America':
        color = '#deebf7'
        #if sample_phenotype_map[sample]==0:
        #    color = '#9ecae1'
        #elif sample_phenotype_map[sample]==1:
        #    color = '#3182bd'
        #else:
        #    color = '#deebf7'
    elif sample_continent_map[sample]=='Africa':
        color = '#31a354'
    else:
        color = '#de2d26'

    dendrogram_axis.plot([x],[y],'o',color=color,markeredgewidth=0,markersize=2)


dendrogram_axis.set_xticklabels(xticklabels)


fig.tight_layout()
#fig.subplots_adjust(hspace=0.2)
fig.savefig("%sdendrogram_%s.png" % (config.analysis_directory, species_name), format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
