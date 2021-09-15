from __future__ import division
import numpy
import pickle
import config

import parse_midas_data
import midas_db_utils

import sample_utils
import parse_HMP_data
import parse_patric



permutations = 10000

n_min_gene_observations = 5

n_min_sampels = 20

def calculate_delta_prevalence():

    extracted_sample_metadata_map = sample_utils.parse_sample_continent_map()

    prevalence_dict = {}

    #good_species_list = ['Prevotella_copri_61740']
    good_species_list = parse_midas_data.parse_good_species_list()

    for species_name in good_species_list:


        midas_shared_genes = midas_db_utils.parse_midas_shared_genes(species_name)

        desired_samples, new_gene_names, gene_presence_matrix, gene_depth_matrix, marker_coverages, gene_reads_matrix = parse_midas_data.parse_pangenome_data(species_name)

        # first filter genes by removing the genes that are in multiple species
        # we can't say whether theyse genes are present or absent for a given species
        new_gene_names_idx = numpy.asarray([ numpy.where(new_gene_names == s)[0][0] for s in new_gene_names if s not in midas_shared_genes])

        if len(new_gene_names_idx) == 0:
            continue

        gene_presence_matrix = gene_presence_matrix[new_gene_names_idx,:]
        new_gene_names = [s for s in new_gene_names if s not in midas_shared_genes]
        new_gene_names = numpy.asarray(new_gene_names)


        # remove genes with few observations across all samples
        new_gene_names = new_gene_names[gene_presence_matrix.sum(axis=1) > n_min_gene_observations]
        gene_presence_matrix = gene_presence_matrix[gene_presence_matrix.sum(axis=1) > n_min_gene_observations,:]


        desired_samples_africa = []
        desired_samples_north_america = []

        for desired_sample in desired_samples:

            if desired_sample in extracted_sample_metadata_map:

                desired_sample_continent = extracted_sample_metadata_map[desired_sample]

                if desired_sample_continent == 'Africa':
                    desired_samples_africa.append(desired_sample)

                elif desired_sample_continent == 'North America':
                    desired_samples_north_america.append(desired_sample)

                else:
                    continue

        # only examine species with at least 20 hosts for each continent
        n_desired_samples_africa = len(desired_samples_africa)
        n_desired_samples_north_america = len(desired_samples_north_america)

        if (n_desired_samples_africa < n_min_sampels) or (n_desired_samples_north_america < n_min_sampels):
            continue

        print(species_name)

        desired_samples_africa_idx = numpy.asarray([ numpy.where(desired_samples == s)[0][0] for s in desired_samples_africa])
        desired_samples_north_america_idx = numpy.asarray([ numpy.where(desired_samples == s)[0][0] for s in desired_samples_north_america])

        gene_presence_matrix_africa = gene_presence_matrix[:,desired_samples_africa_idx]
        gene_presence_matrix_north_america = gene_presence_matrix[:,desired_samples_north_america_idx]

        # new matrix of all africa and n america samples
        gene_presence_matrix = numpy.concatenate((gene_presence_matrix_africa, gene_presence_matrix_north_america), axis=1)


        prevalence_africa = gene_presence_matrix_africa.sum(axis=1) / gene_presence_matrix_africa.shape[0]
        prevalence_north_america = gene_presence_matrix_north_america.sum(axis=1) / gene_presence_matrix_north_america.shape[0]


        delta_prevalence = prevalence_africa - prevalence_north_america

        delta_prevalence_null_arrays = []

        for p in range(permutations):

            numpy.random.shuffle(gene_presence_matrix.T)

            gene_presence_matrix_africa_null = gene_presence_matrix[:,:n_desired_samples_africa]
            gene_presence_matrix_north_america_null = gene_presence_matrix[:,n_desired_samples_africa:]

            prevalence_africa_null = gene_presence_matrix_africa_null.sum(axis=1) / gene_presence_matrix_africa_null.shape[0]
            prevalence_north_america_null = gene_presence_matrix_north_america_null.sum(axis=1) / gene_presence_matrix_north_america_null.shape[0]

            delta_prevalence_null = prevalence_africa_null - prevalence_north_america_null
            delta_prevalence_null_arrays.append(delta_prevalence_null)

        prevalence_dict[species_name] = {}
        prevalence_dict[species_name]['n_africa'] = n_desired_samples_africa
        prevalence_dict[species_name]['n_north_america'] = n_desired_samples_north_america
        prevalence_dict[species_name]['n_permutations'] = permutations


        delta_prevalence_null_matrix = numpy.array(delta_prevalence_null_arrays)
        #we're looking at the absolute value of delta prevalence
        mean_delta_prevalence_null = numpy.mean(numpy.absolute(delta_prevalence_null_matrix), axis=0)


        prevalence_dict[species_name]['delta_prevalence_observed'] = numpy.absolute(delta_prevalence).tolist()
        prevalence_dict[species_name]['delta_prevalence_null'] = mean_delta_prevalence_null.tolist()

        prevalence_dict[species_name]['genes'] = {}


        groups = new_gene_names[0].split('.')
        genome_id = '.'.join(groups[:2])#, '.'.join(groups[2:])

        kegg_ids = parse_patric.load_kegg_annotations([genome_id])

        #537011.5.peg.2684

        # loop through, get 95% CIs for each gene
        for gene_delta_prevalence_null_idx, gene_delta_prevalence_null in enumerate(delta_prevalence_null_matrix.T):

            gene_name = new_gene_names[gene_delta_prevalence_null_idx]

            gene_delta_prevalence_null = numpy.sort(gene_delta_prevalence_null)

            gene_delta_prevalence = delta_prevalence[gene_delta_prevalence_null_idx]

            gene_delta_prevalence_low_ci = gene_delta_prevalence_null[int(permutations*0.025)]
            gene_delta_prevalence_high_ci = gene_delta_prevalence_null[int(permutations*0.975)]

            p_value = (sum(numpy.absolute(gene_delta_prevalence_null) > numpy.absolute(gene_delta_prevalence))+1) / permutations

            prevalence_dict[species_name]['genes'][gene_name] = {}
            prevalence_dict[species_name]['genes'][gene_name]['delta_prevalence'] = gene_delta_prevalence
            prevalence_dict[species_name]['genes'][gene_name]['low_ci'] = gene_delta_prevalence_low_ci
            prevalence_dict[species_name]['genes'][gene_name]['high_ci'] = gene_delta_prevalence_high_ci
            prevalence_dict[species_name]['genes'][gene_name]['p_value'] = p_value

            if gene_name in kegg_ids:
                prevalence_dict[species_name]['genes'][gene_name]['kegg_id'] = kegg_ids[gene_name]

            else:
                prevalence_dict[species_name]['genes'][gene_name]['kegg_id'] = (gene_name, [['', '']])



    intermediate_filename_template = config.data_directory+"delta_prevalence_test.dat"
    intermediate_filename = intermediate_filename_template

    with open(intermediate_filename, 'wb') as handle:
        pickle.dump(prevalence_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)




def load_predicted_prevalence_subsample_dict():

    intermediate_filename = config.data_directory+"delta_prevalence_test.dat"

    with open(intermediate_filename, 'rb') as handle:
        b = pickle.load(handle)
    return b


def get_kegg_matrix(species_name='Prevotella_copri_61740'):


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
            if k == ['', '']:
                continue

            pathway = k[0]
            annotation = k[1]

            if pathway not in pathway_annotation_dict:
                pathway_annotation_dict[pathway] = annotation

                kegg_count_dict[gene][pathway] = 1

            kegg_count_dict[gene]['delta_prevalence'] = gene_dict[gene]['delta_prevalence']



    return genes, delta_prevalence, kegg_count_dict, pathway_annotation_dict




if __name__=='__main__':

    calculate_delta_prevalence()
