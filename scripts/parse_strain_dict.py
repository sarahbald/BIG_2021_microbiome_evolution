import pickle

data_directory = "C:/Users/sarah/Garud Lab/BIG_2021_microbiome_evolution/"
#data_directory = "/Users/williamrshoemaker/Desktop/"
filename = data_directory+"strain_dict.pickle"
with open(filename, 'rb') as handle:
    strain_dict = pickle.load(handle)