import pickle
#import numpy
import gene_prev_utils
import matplotlib.pyplot as plt

data_directory = "C:/Users/sarah/Garud Lab/"
#data_directory = "/Users/williamrshoemaker/Desktop/"
filename = data_directory+"lognormal_params_dict.pickle"
with open(filename, 'rb') as handle:
    lognormal_params_dict = pickle.load(handle)

species_list = []
af_mad = []
na_mad = []
for species in lognormal_params_dict["Species Abundance Dict"].keys():
    species_list.append(species)
    af_mad.append(lognormal_params_dict["Species Abundance Dict"][species]["Africa"])
    na_mad.append(lognormal_params_dict["Species Abundance Dict"][species]["North America"])
evo_species = gene_prev_utils.load_predicted_prevalence_subsample_dict(want_names=True)
fig, ax = plt.subplots(figsize=(4,4))
for i, val in enumerate(af_mad):
    if species_list[i] in evo_species:
        ax.scatter(na_mad[i], val, c='r', alpha = 0.4)
    else:
        ax.scatter(na_mad[i], val, c='blue', alpha = 0.4)
ax.plot([1e-8, 1e-2, 1e-1],[1e-8, 1e-2, 1e-1], linestyle="-")
ax.set_xlabel('North America Species Mean Abundance')
ax.set_ylabel('Africa Species Mean Abundance')
ax.set_ylim(min(af_mad)*.5, max(af_mad)*1.8)
ax.set_xlim(min(na_mad)*.5, max(na_mad)*1.8)
ax.set_xscale("log")
ax.set_yscale("log")

fig_name = data_directory+ 'mad_na_versus_af.png'
fig.savefig(fig_name, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.show()
plt.close()
    