from __future__ import division
import numpy
import pickle
import config
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import figure_utils

#import calculate_delta_prevalence
import parse_midas_data
import parse_patric
import core_gene_utils

import pandas


from sklearn.linear_model import LassoCV
#from sklearn.linear_model import Lasso
#from sklearn.linear_model import LinearRegression
#from sklearn.model_selection import train_test_split



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




# load data
genes, delta_prevalence, kegg_count_dict, pathway_annotation_dict = get_kegg_matrix()


#make matrix of explanatory variables
kegg_count_df = pandas.DataFrame.from_dict(kegg_count_dict)
kegg_count_df = kegg_count_df.fillna(0)
kegg_count_df = kegg_count_df.T


kegg_ids = kegg_count_df.columns.values[:-1]

kegg_count_array = kegg_count_df.values
X, y = kegg_count_array[:, :-1], kegg_count_array[:, -1]
# X = annotation presence-absence matrix (independent variable)
# y= delta prevalence array (dependent variable)

# this is your range of penelaty terms  in log10 space
# -9 = log10(10**-9), 0 = log10(1)
alpha_array = numpy.logspace(-9, 0, num=50)
#smaller alpha = closer to OLS regression

# fitting and cross-validating the lasso regression
reg = LassoCV(cv=50, alphas = alpha_array, random_state=0).fit(X, y)

# get the predicted prevalence
predicted_delta_prev = reg.predict(X)

# get the slopes of each function
slopes = numpy.abs(reg.coef_)



# things for Sarah to do
# two figures
# 1) plot the predicted and the observed delta prevalence values as a 1:1 plot
#    i) print the coefficeint of determination (R^2) on the plot
#     (hint, calculate the residual sum of squares and total sum of squares )
#    ii) (optional) calculate cross-validated R^2

# 2) Plot the slopes of all selected KEGG IDs
