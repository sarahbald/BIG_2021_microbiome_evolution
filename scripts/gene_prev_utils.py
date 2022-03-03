import pickle
import numpy

def load_predicted_prevalence_subsample_dict(want_names=False):

    intermediate_filename = 'C:/Users/sarah/Garud Lab/delta_prevalence_test.dat'

    with open(intermediate_filename, 'rb') as handle:
#        b = pickle.load(handle, encoding='latin1')
        b=pickle.load(handle)
    if want_names == False:
        return b
    if want_names == True:
        return b.keys()



def get_kegg_matrix(species_name):
    
    prevalence_dict = load_predicted_prevalence_subsample_dict()
    gene_dict = prevalence_dict[species_name]['genes']

    kegg_count_dict = {}
    delta_prevalence = []
    pathway_annotation_dict = {}
    genes = []
    for gene in gene_dict.keys():

        kegg_id = gene_dict[gene]['kegg_id']


        if kegg_id == [['', '']]:
            continue

        delta_prevalence.append(gene_dict[gene]['delta_prevalence'])
        genes.append(gene)

        kegg_count_dict[gene] = {}

        for k in kegg_id:
            #if k == ['', '']:
            if k == [['', '']]:
                continue
            pathway = k[0]
            annotation = k[1]

            if pathway not in pathway_annotation_dict:
                pathway_annotation_dict[pathway] = annotation

                kegg_count_dict[gene][pathway] = 1

            kegg_count_dict[gene]['delta_prevalence'] = gene_dict[gene]['delta_prevalence']



    return genes, delta_prevalence, kegg_count_dict, pathway_annotation_dict


def calculate_r_squared(true_prev, predicted_prev):
    mean_y = numpy.mean(predicted_prev)
    tss = 0
    for y in predicted_prev:
        tss += (y - mean_y)**2
    rss = 0
    for i in range(len(true_prev)):
        rss += (true_prev[i] - predicted_prev[i])**2
    r_squared = (tss-rss)/tss
    return r_squared

