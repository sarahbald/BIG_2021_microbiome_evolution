from __future__ import division
import numpy
#import pickle
#import config
import matplotlib.pyplot as plt
#import matplotlib.gridspec as gridspec
import gene_prev_utils

import pandas


from sklearn.linear_model import LassoCV
#from sklearn.linear_model import LassoLarsCV
#from sklearn.linear_model import Lasso
#from sklearn.linear_model import LinearRegression
#from sklearn.model_selection import train_test_split



all_species_names = gene_prev_utils.load_predicted_prevalence_subsample_dict(want_names=True)
important_kegg_ids = []
nonzero_slope_values = []
important_gene_dict = {}
predicted_delta_prev_plot = []
true_delta_prev_plot = []
coef_function_dict = {}
pathway_annotation_dict = {}
kegg_count_dict = {}
genes = []
delta_prevalence = []
predicted_delta_prev = []

#run model for all species separately, and compile results in 1 plot
for species in all_species_names:
    #load data
    new_genes, new_delta_prevalence, new_kegg_count_dict, new_pathway_annotation_dict = gene_prev_utils.get_kegg_matrix(species)
    #accumulate the dictionaries and other outputs between species to have cummulative lists
    pathway_annotation_dict.update(new_pathway_annotation_dict) 
    kegg_count_dict.update(new_kegg_count_dict)
    genes.extend(new_genes)
    delta_prevalence.extend(new_delta_prevalence)
    
    #make matrix of explanatory variables for 1 species
    kegg_count_df = pandas.DataFrame.from_dict(new_kegg_count_dict)
    kegg_count_df = kegg_count_df.fillna(0)
    kegg_count_df = kegg_count_df.T
    kegg_count_df = kegg_count_df.drop('delta_prevalence', axis=1)

    
    kegg_ids = kegg_count_df.columns.values

    kegg_count_array = kegg_count_df.values
    X, y = kegg_count_array, new_delta_prevalence
    # X = annotation presence-absence matrix (independent variable)
    # y= delta prevalence array (dependent variable)

    # this is your range of penelaty terms  in log10 space
    # -9 = log10(10**-9), 0 = log10(1)
    alpha_array = numpy.logspace(-9, 0, num=50)
    #smaller alpha = closer to OLS regression

    # fitting and cross-validating the lasso regression
    reg = LassoCV(cv=50, alphas = alpha_array, random_state=0).fit(X, y)

    # get the predicted prevalence
    new_predicted_delta_prev = reg.predict(X)
    predicted_delta_prev.extend(new_predicted_delta_prev)

    # get the slopes of each function
    slopes = numpy.abs(reg.coef_)
    coefficients = reg.coef_

    #get list of kegg IDs that correspond to nonzero slopes
    nonzero_slope_idx = numpy.nonzero(slopes)[-1]
    important_kegg_ids.append(kegg_ids[nonzero_slope_idx])
    
    #find genes with these important kegg ids, and save in dictionary that maps gene to kegg ID and delta prevalence
    for value in coefficients:
        if not value == 0:
            nonzero_slope_values.append(value)
            for gene_dict in new_kegg_count_dict.keys():
                gene_added = False
                for kegg_id in important_kegg_ids:
                    if kegg_id.shape == (0,):
                        continue
                    if gene_added == True:
                        continue
                    if kegg_id[0] in new_kegg_count_dict[gene_dict].keys():
                        important_gene_dict[gene_dict] = new_kegg_count_dict[gene_dict]
                        gene_added = True

#construct coef_function_dict that maps function to nonzero slope value from model
k=0
for i in range(len(important_kegg_ids)):
    if important_kegg_ids[i].shape == (0,):
        continue
    for j in range(len(important_kegg_ids[i])):
        if pathway_annotation_dict[important_kegg_ids[i][j]] not in coef_function_dict:
            coef_function_dict[pathway_annotation_dict[important_kegg_ids[i][j]]] = []
        coef_function_dict[pathway_annotation_dict[important_kegg_ids[i][j]]].append(nonzero_slope_values[k])
        k+=1

#recover the true and predicted values for select subset of genes
idxs = []
predicted_delta_prev_plot = []
true_delta_prev_plot = []
for gene_id in important_gene_dict.keys():
    idxs.append(genes.index(gene_id))
for idx in idxs:
    predicted_delta_prev_plot.append(predicted_delta_prev[idx])
    true_delta_prev_plot.append(delta_prevalence[idx])
 
#plot observed delta prev vs predicted delta prev plot
plt.figure(1)
plt.scatter(true_delta_prev_plot, predicted_delta_prev_plot, label='gene')
xy_line = numpy.linspace(min(delta_prevalence)-0.01, max(delta_prevalence)+0.01)
plt.plot(xy_line,xy_line, color = 'red', label = '1:1', linewidth=0.5)
plt.xlim(min(delta_prevalence)-0.01,max(delta_prevalence)+0.01)
plt.ylim(min(delta_prevalence)-0.01,max(delta_prevalence)+0.01)
plt.xlabel('Observed ' + r'$\Delta p$')
plt.ylabel('Predicted ' + r'$\Delta p$')
plt.title(r'$\Delta p$' + ' ~ Kegg IDs, Composite Lasso Regression for 8 Prevalent Species', fontdict = {'fontsize' : 10})
r_squared = round(gene_prev_utils.calculate_r_squared(true_delta_prev_plot, predicted_delta_prev_plot), 3)
plt.text(-0.02,-0.05, r'$R^{2}=$' + str(r_squared))
plt.legend()
#plt.savefig('C:/Users/sarah/Garud Lab/plots/lasso_regression.png', dpi=600)

#drop non-zero coefficients with values less than 1e-15
function_names_plot = []
function_coef_plot = []
for function in coef_function_dict.keys():
    for coef in coef_function_dict[function]:
        if abs(coef) < 1E-15:
            continue
        function_names_plot.append(function)
        function_coef_plot.append(coef)

#plot histogram of nonzero 
plt.figure(2)
plt.bar(range(len(function_names_plot)),function_coef_plot)
plt.xticks(range(len(function_names_plot)), function_names_plot, rotation='vertical')
plt.title("Significant Functions for predicting Delta Prevalence")
#plt.savefig('C:/Users/sarah/Garud Lab/plots/function_coef.png', dpi=600, bbox_inches='tight')



#ONE MODEL RUN TOGETHER ACROSS ALL SPECIES INSTEAD OF SEPARATELY
#make matrix of explanatory variables for 1 species
kegg_count_df = pandas.DataFrame.from_dict(kegg_count_dict)
kegg_count_df = kegg_count_df.fillna(0)
kegg_count_df = kegg_count_df.T
kegg_count_df = kegg_count_df.drop('delta_prevalence', axis=1)

    
kegg_ids = kegg_count_df.columns.values
kegg_count_array = kegg_count_df.values
X, y = kegg_count_array, delta_prevalence
# X = annotation presence-absence matrix (independent variable)
# y= delta prevalence array (dependent variable)

# this is your range of penelaty terms  in log10 space
# -9 = log10(10**-9), 0 = log10(1)
alpha_array = numpy.logspace(-4, 0, num=50)
#smaller alpha = closer to OLS regression

# fitting and cross-validating the lasso regression
reg_combined = LassoCV(cv=50, alphas=alpha_array, random_state=0).fit(X, y)

# get the predicted prevalence
predicted_delta_prev_combined = reg_combined.predict(X)

# get the slopes of each function
slopes_combined = numpy.abs(reg_combined.coef_)
coefficients_combined = reg_combined.coef_

#get list of kegg IDs that correspond to nonzero slopes
nonzero_slope_idx_combined = numpy.nonzero(slopes_combined)[-1]
important_kegg_ids_combined = kegg_ids[nonzero_slope_idx]
    
nonzero_slope_values_combined = []
important_gene_dict_combined = {}
#find genes with these important kegg ids, and save in dictionary that maps gene to kegg ID and delta prevalence
for value in coefficients_combined:
    if not value == 0:
        nonzero_slope_values_combined.append(value)
        for gene_dict in kegg_count_dict.keys():
            gene_added = False
            for kegg_id in important_kegg_ids_combined:
                if kegg_id.shape == (0,):
                    continue
                if gene_added == True:
                    continue
                if kegg_id[0] in kegg_count_dict[gene_dict].keys():
                    print('kegg_id')
                    important_gene_dict_combined[gene_dict] = kegg_count_dict[gene_dict]
                    gene_added = True

#construct coef_function_dict that maps function to nonzero slope value from model
coef_function_dict_combined = {}
k=0
for i in range(len(important_kegg_ids_combined)):
    if important_kegg_ids_combined[i].shape == (0,):
        continue
    for j in range(len(important_kegg_ids_combined[i])):
        if pathway_annotation_dict[important_kegg_ids_combined[i][j]] not in coef_function_dict_combined:
            coef_function_dict_combined[pathway_annotation_dict[important_kegg_ids_combined[i][j]]] = []
        coef_function_dict_combined[pathway_annotation_dict[important_kegg_ids_combined[i][j]]].append(nonzero_slope_values_combined[k])
        k+=1

idxs = []
predicted_delta_prev_plot = []
true_delta_prev_plot = []
for gene_id in important_gene_dict_combined.keys():
    idxs.append(genes.index(gene_id))
for idx in idxs:
    predicted_delta_prev_plot.append(predicted_delta_prev[idx])
    true_delta_prev_plot.append(delta_prevalence[idx])
 
#plot observed delta prev vs predicted delta prev plot
plt.figure(3)
plt.scatter(true_delta_prev_plot, predicted_delta_prev_plot, label='gene')
xy_line = numpy.linspace(min(delta_prevalence)-0.01, max(delta_prevalence)+0.01)
plt.plot(xy_line,xy_line, color = 'red', label = '1:1', linewidth=0.5)
plt.xlim(min(delta_prevalence)-0.01,max(delta_prevalence)+0.01)
plt.ylim(min(delta_prevalence)-0.01,max(delta_prevalence)+0.01)
plt.xlabel('Observed ' + r'$\Delta p$')
plt.ylabel('Predicted ' + r'$\Delta p$')
plt.title(r'$\Delta p$' + ' ~ Kegg IDs, Combined Lasso Regression for 8 Prevalent Species', fontdict = {'fontsize' : 10})
r_squared = round(gene_prev_utils.calculate_r_squared(true_delta_prev_plot, predicted_delta_prev_plot), 3)
plt.text(-0.02,-0.05, r'$R^{2}=$' + str(r_squared))
plt.legend()
plt.savefig('C:/Users/sarah/Garud Lab/plots/lasso_regression_combined.png', dpi=600)

#drop non-zero coefficients with values less than 1e-15
function_names_plot = []
function_coef_plot = []
for function in coef_function_dict_combined.keys():
    for coef in coef_function_dict_combined[function]:
        if abs(coef) < 1E-15:
            continue
        function_names_plot.append(function)
        function_coef_plot.append(coef)

#plot histogram of nonzero 
plt.figure(4)
plt.bar(range(len(function_names_plot)),function_coef_plot)
plt.xticks(range(len(function_names_plot)), function_names_plot, rotation='vertical')
plt.title("Significant Functions for predicting Delta Prevalence")
plt.savefig('C:/Users/sarah/Garud Lab/plots/function_coef_combined.png', dpi=600, bbox_inches='tight')




















###BOXPLOT CODE FROM BIG SUMMER PLOT, VERY SPECIFIC TO THAT ONE GRAPH

#small_coef_values = []
#small_coef_functions = []
#big_coef_values = []
#big_coef_functions = []

#for value in coef_function_dict.keys():
#    if numpy.abs(value) > 10**-5:
#        big_coef_functions.append(coef_function_dict[value])
#        big_coef_values.append(value)
#    else:
#        small_coef_functions.append(coef_function_dict[value])
#        small_coef_values.append(value)

#plt.clf()
#x_pos = [i for i, _ in enumerate(small_coef_functions)]
#fig, axs = plt.subplots(2)
#plt.subplots_adjust(wspace=0.2, 
#                    hspace=0.4)
#axs[1].text(-1.5,4e-17,"Coefficient Value", rotation='vertical')
#fig.tight_layout()
#fig.suptitle('Feature Coefficients in Regression Model for Faecalibacterium cf')
#axs[0].bar(x_pos, big_coef_values)
#axs[0].set_yscale('symlog', linthresh=10e-3)
#axs[0].set_yticks([-10e-2,-10e-3, 0, 10e-3, 10e-2, 0.04])
#axs[0].set_yticklabels([-10e-2,-10e-3, 0, 10e-3, 10e-2, 0.04])
#axs[0].set_ylim(-10e-3,5e-2)
#axs[0].set_xticks(x_pos)
#axs[0].set_xticklabels(['Streptomycin \n biosynthesis', 'Folate \n biosynthesis', 'Galactose \n metabolism', 'Fructose and  \n mannose metabolism', 'Tetrachloroethene \n degradation'], fontsize = 6)

#axs[1].bar(x_pos, small_coef_values)
#axs[1].axhline(0.0,color='black',linewidth=0.5)
#axs[1].set_yscale('symlog', linthresh=10e-20) 
#axs[1].set_ylim(-10e-18, 5e-17)
#axs[1].set_yticks([-10e-18, 0, 10e-18, 4e-17])
#axs[1].set_yticklabels([-10e-17, 0, 10e-18, 4e-17])
#axs[1].set_xticks(x_pos)
#axs[1].set_xticklabels(['Polyketide sugar \n unit biosynthesis', 'Glycosaminoglycan \n degradation', 'O-Glycan \n biosynthesis', 'Retinol \n metabolism', 'Metabolism \n of xenobiotics'], fontsize = 6)

#plt.savefig('Coefficient_barplot.png', dpi=600)

