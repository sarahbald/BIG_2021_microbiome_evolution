import bz2
import numpy as np
import matplotlib.pyplot as plt
import math
import sympy as sp
import scipy
#from scipy.special import erf
#import pickle

abun_dat_path = 'C:/Users/sarah/Garud Lab/BIG_2021_microbiome_evolution/abundance_stoolsubset.txt.bz2'
dat_file = bz2.BZ2File(abun_dat_path)
abundances = []
groups = dat_file.readline()
groups = groups.strip()
groups = groups.split('\t')
for line in dat_file:
    line = line.split('\t')
    info = [line[0]]
    #isolate hmp samples
    for i in range(317,469):
        info.append(line[i])
    abundances.append(info)

abun_data = [abundances[0]]
abun_data.append(abundances[1])
abun_data.append(abundances[77])
for i in range(210,len(abundances)):
    abun_data.append(abundances[i])

abun_data_transpose = [[row[i] for row in abun_data] for i in range(len(abun_data[0]))]
abun_matrix_raw = np.array(abun_data_transpose)
abun_matrix = []
for idx in range(abun_matrix_raw.shape[1]):
    label = abun_matrix_raw[0,idx]
    if 'k__' in label and 's__' in label and 't__' not in label:
        abun_matrix.append(abun_matrix_raw[:, idx])
abun_matrix = np.array(abun_matrix)
zero_abun_species = []
for idx in range(abun_matrix.shape[0]):
    abun = [float(i) for i in abun_matrix[idx,1:]]
    if np.sum(abun) == 0.0:
        zero_abun_species.append(idx)
abun_matrix = np.delete(abun_matrix, zero_abun_species, 0)

mean_abun = []
for idx in range(abun_matrix.shape[0]):
    abun = [float(i) for i in abun_matrix[idx,1:]]
    mean_abun.append(np.sum(abun)/(100*len(abun)))   




def Klogn(emp_mad, c, mu0=-19,s0=5):
    # This function estimates the parameters (mu, s) of the lognormal distribution of K
    m1 = np.mean(np.log(emp_mad[emp_mad>c]))
    m2 = np.mean(np.log(emp_mad[emp_mad>c])**2)
    xmu = sp.Symbol('xmu')
    xs = sp.Symbol('xs')
    eq1 = - m1 + xmu + np.sqrt(2/math.pi)*xs*(sp.exp(-((np.log(c)-xmu)**2)/2/(xs**2))/(sp.erfc((np.log(c)-xmu)/np.sqrt(2)/xs)))
    eq2 = - m2 + xs**2 + m1*xmu+np.log(c)*m1-xmu*np.log(c)
    sol = sp.nsolve([eq1,eq2],[xmu,xs],[mu0,s0])

    return(float(sol[0]),float(sol[1]))

mad_log10 = np.log10(mean_abun)
hist_mad, bin_edges_mad = np.histogram(mad_log10, density=True, bins=20)
bins_mean_mad = [0.5 * (bin_edges_mad[i] + bin_edges_mad[i+1]) for i in range(len(bin_edges_mad)-1 )]
prob_to_plot = [sum( (mad_log10>=bin_edges_mad[i]) & (mad_log10<=bin_edges_mad[i+1])  ) / (len(mad_log10)*1.0) for i in range(len(bin_edges_mad)-1 )]
bins_mean_plot = [10**(bins_mean_mad[i]) for i in range(len(bins_mean_mad))]

x=np.arange(-14,-1,step=0.1)

plt.figure()

plt.scatter(bins_mean_plot, prob_to_plot, alpha=0.5, s=30)
plt.xlabel('Mean Abundance Distribution')
plt.xscale('log')
plt.xlim(10e-9,0.5)
plt.ylim(10e-3, .2)
plt.savefig('C:/Users/sarah/Garud Lab/BIG_2021_microbiome_evolution/plots/metaphlan_mad.png', dpi=600)
plt.yscale('log')

c = 10e-6
(mu,sigma) = Klogn(np.array(mean_abun), c)
#plt.plot(np.exp(x), np.sqrt(2/math.pi)/sigma *np.exp(-(x-mu)**2 /2/(sigma**2))/scipy.special.erfc((np.log(c)-mu)/np.sqrt(2)/sigma))
c = 10e-7
(mu,sigma) = Klogn(np.array(mean_abun), c)
plt.plot(np.exp(x), np.sqrt(2/math.pi)/sigma *np.exp(-(x-mu)**2 /2/(sigma**2))/scipy.special.erfc((np.log(c)-mu)/np.sqrt(2)/sigma))  
c = 10e-8
(mu,sigma) = Klogn(np.array(mean_abun), c)
plt.plot(np.exp(x), np.sqrt(2/math.pi)/sigma *np.exp(-(x-mu)**2 /2/(sigma**2))/scipy.special.erfc((np.log(c)-mu)/np.sqrt(2)/sigma))
plt.savefig('C:/Users/sarah/Garud Lab/BIG_2021_microbiome_evolution/plots/metaphlan_mad_estimates.png', dpi=600)
        
    