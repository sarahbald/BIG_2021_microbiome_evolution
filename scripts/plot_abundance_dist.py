from __future__ import division
#import calculate_species_relative_abundance
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as stats
import pickle
#from scipy.special import erf

import math
import sympy as sp
import scipy

data_directory = "C:/Users/sarah/Garud Lab/"

filename = data_directory+"species_abun_dict.dat"
with open(filename, 'rb') as handle:
        species_abun_dict_zeros = pickle.load(handle)
species_abun_dict = {}
for species in species_abun_dict_zeros.keys():
    species_abun_dict[species] = {}
    na_nonzero_abun = []
    af_nonzero_abun = []
    for abun in species_abun_dict_zeros[species]['Africa']:
        if abun != '0.0':
            af_nonzero_abun.append(abun)
    for abun in species_abun_dict_zeros[species]['North America']:
        if abun != '0.0':
            na_nonzero_abun.append(abun)
    species_abun_dict[species]['Africa'] = af_nonzero_abun
    species_abun_dict[species]['North America'] = na_nonzero_abun
    
    
#####collect all species abundances and variance of abundance within group for each species
varplot_dict = {}
africa_abun = []
na_abun = []
for species in species_abun_dict.keys():
    africa_abun.append(species_abun_dict[species]['Africa'])
    na_abun.append(species_abun_dict[species]['North America'])
    if len(species_abun_dict[species]['Africa']) > 1 and len(species_abun_dict[species]['North America']) > 1:
        varplot_dict[species] = {}
        varplot_dict[species]['Africa'] = []
        varplot_dict[species]['North America'] = []
        varplot_dict[species]['Africa'].append(np.var(np.array(species_abun_dict_zeros[species]['Africa']).astype(np.float)))
        varplot_dict[species]['North America'].append(np.var(np.array(species_abun_dict_zeros[species]['North America']).astype(np.float)))

#flatten lists of lists into 1 list    
africa_abun = [float(item) for sublist in africa_abun for item in sublist]
na_abun = [float(item) for sublist in na_abun for item in sublist]

#####collect log mean abundance within group for each species
filename = data_directory+"species_mean_relative_abundance.dat"
with open(filename, 'rb') as handle:
        mean_abundance_dict = pickle.load(handle)
log_africa_mean_abun = []
log_north_america_mean_abun = []
for species_abundance in mean_abundance_dict.keys():
    log_africa_mean_abun.append(np.log10(mean_abundance_dict[species_abundance]['Africa']))
    log_north_america_mean_abun.append(np.log10(mean_abundance_dict[species_abundance]['North America'])) 

#get list of mean abundances for plotting and mean abundances for varplot_dict
africa_mean_abun = []
na_mean_abun = []
species_order_mean = []
for species in mean_abundance_dict.keys():
    species_order_mean.append(species)
    africa_mean_abun.append(mean_abundance_dict[species]['Africa'])
    na_mean_abun.append(mean_abundance_dict[species]['North America'])
    if len(species_abun_dict[species]['Africa']) > 1 and len(species_abun_dict[species]['North America']) > 1:
        varplot_dict[species]['Africa'].append(mean_abundance_dict[species]['Africa'])
        varplot_dict[species]['North America'].append(mean_abundance_dict[species]['North America'])

filename = data_directory+"species_prevalence.dat"
with open(filename, 'rb') as handle:
        species_prev_dict = pickle.load(handle)
africa_prev = []
na_prev = []
species_order_prev = []
for species in species_prev_dict.keys():
    species_order_prev.append(species)
    africa_prev.append(species_prev_dict[species]['Africa'])
    na_prev.append(species_prev_dict[species]['North America'])
    
#rescale abundances
africa_abun = np.array(africa_abun)
na_abun = np.array(na_abun)
rescaled_africa_abun = (np.log10(africa_abun) - np.log10(africa_abun).mean())/np.log10(africa_abun).std()
rescaled_na_abun = (np.log10(na_abun) - np.log10(na_abun).mean())/np.log10(na_abun).std()

fig, axs = plt.subplots(2, 2, figsize=(15,15))


hist_na, bin_edges_na = np.histogram(rescaled_na_abun, density=True, bins=30)
bins_mean_na = [0.5 * (bin_edges_na[i] + bin_edges_na[i+1]) for i in range(0, len(bin_edges_na)-1 )]

hist_af, bin_edges_af = np.histogram(rescaled_africa_abun, density=True, bins=30)
bins_mean_af = [0.5 * (bin_edges_af[i] + bin_edges_af[i+1]) for i in range(0, len(bin_edges_af)-1 )]


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

na_mad_log10 = np.log10(na_mean_abun)
na_hist_mad, na_bin_edges_mad = np.histogram(na_mad_log10, density=True, bins=20)
na_bins_mean_mad = [0.5 * (na_bin_edges_mad[i] + na_bin_edges_mad[i+1]) for i in range(0, len(na_bin_edges_mad)-1 )]
na_prob_to_plot = [sum( (na_mad_log10>=na_bin_edges_mad[i]) & (na_mad_log10<=na_bin_edges_mad[i+1])  ) / len(na_mad_log10) for i in range(0, len(na_bin_edges_mad)-1 )]
na_bins_mean_plot = [10**(na_bins_mean_mad[i]) for i in range(len(na_bins_mean_mad))]

af_mad_log10 = np.log10(africa_mean_abun)
af_hist_mad, af_bin_edges_mad = np.histogram(af_mad_log10, density=True, bins=20)
af_bins_mean_mad = [0.5 * (af_bin_edges_mad[i] + af_bin_edges_mad[i+1]) for i in range(0, len(af_bin_edges_mad)-1 )]
af_prob_to_plot = [sum( (af_mad_log10>=af_bin_edges_mad[i]) & (af_mad_log10<=af_bin_edges_mad[i+1])  ) / len(af_mad_log10) for i in range(0, len(af_bin_edges_mad)-1 )]
af_bins_mean_plot = [10**(af_bins_mean_mad[i]) for i in range(len(af_bins_mean_mad))]


plt.clf()
fig, axs = plt.subplots(2, 2, figsize=(15,15))

# range of abundances to plot
x=np.arange(-14,-1,step=0.1)

axs[0,0].scatter(na_bins_mean_plot, na_prob_to_plot, alpha=0.5, s=30, label = "North America")
axs[0,0].scatter(af_bins_mean_plot, af_prob_to_plot, alpha=0.5, s=30, label = "Africa")
#ax.scatter(bins_mean_af, hist_af, alpha=0.5, s=30, label = "Africa")
#c_values = [10e-8, 10e-7, 10e-6]
c = [1e-8, 1e-6]
#for value in numpy.linspace(numpy.quantile(emp_mad, 0.25), numpy.quantile(emp_mad, 0.5), 5):
(mu,sigma) = Klogn(np.array(na_mean_abun), c[0])
axs[0,0].plot(np.exp(x), np.sqrt(2/math.pi)/sigma *np.exp(-(x-mu)**2 /2/(sigma**2))/scipy.special.erfc((np.log(c[0])-mu)/np.sqrt(2)/sigma))
(mu,sigma) = Klogn(np.array(africa_mean_abun), c[1])
axs[0,0].plot(np.exp(x), np.sqrt(2/math.pi)/sigma *np.exp(-(x-mu)**2 /2/(sigma**2))/scipy.special.erfc((np.log(c[1])-mu)/np.sqrt(2)/sigma))
#plot truncated lognormal fitted for K>log(c)
#ax.plot(numpy.log10(numpy.exp(x)),numpy.sqrt(2/math.pi)/sigma *numpy.exp(-(x-mu)**2 /2/(sigma**2))/scipy.special.erfc((numpy.log(c)-mu)/numpy.sqrt(2)/sigma),'-k')
#ax.plot(numpy.exp(x), numpy.sqrt(2/math.pi)/sigma *numpy.exp(-(x-mu)**2 /2/(sigma**2))/scipy.special.erfc((numpy.log(c)-mu)/numpy.sqrt(2)/sigma),'-k')
axs[0,0].set_xscale('log', basex=10)
axs[0,0].set_yscale('log', basey=10)
axs[0,0].set_xlim(10**(-9), 1)
axs[0,0].legend()
axs[0,0].set_xlabel('Mean Relative Abundance Distribution')

axs[0, 1].scatter(bins_mean_na, hist_na, alpha=0.5, s=30, label = "North America")
axs[0, 1].scatter(bins_mean_af, hist_af, alpha=0.5, s=30, label = "Africa")
axs[0, 1].set_xlabel('Rescaled Log Abundance Distribution')
axs[0, 1].legend()



na_var = [varplot_dict[s]['North America'][0] for s in varplot_dict.keys()]
af_var = [varplot_dict[s]['Africa'][0] for s in varplot_dict.keys()]
na_mean_varplot = [varplot_dict[s]['North America'][1] for s in varplot_dict.keys()]
af_mean_varplot = [varplot_dict[s]['Africa'][1] for s in varplot_dict.keys()]

na_slope, intercept, r_value, p_value, se = stats.linregress(np.log10(na_mean_varplot), np.log10(na_var))
x_log10_range =  np.linspace(min(np.log10(na_mean_varplot)) , max(np.log10(na_mean_varplot)) , 1000)
y_log10_fit_range = na_slope*x_log10_range + intercept
axs[1, 0].plot(10**x_log10_range, 10**y_log10_fit_range, lw=2, linestyle='--', zorder=2)

af_slope, intercept, r_value, p_value, se = stats.linregress(np.log10(af_mean_varplot), np.log10(af_var))
x_log10_range =  np.linspace(min(np.log10(af_mean_varplot)) , max(np.log10(af_mean_varplot)) , 1000)
y_log10_fit_range = af_slope*x_log10_range + intercept
axs[1, 0].plot(10**x_log10_range, 10**y_log10_fit_range, lw=2, linestyle='--', zorder=2)

axs[1, 0].scatter(na_mean_varplot, na_var, alpha=0.5, label=str('North America (slope = %.2f)' % na_slope))
axs[1, 0].scatter(af_mean_varplot, af_var, alpha=0.5, label =str('Africa (slope = %.2f)' % af_slope))
axs[1, 0].set_xlabel('Mean Relative Abundance')
axs[1, 0].set_ylabel('Variance of Abundance')
axs[1, 0].set_xscale("log")
axs[1, 0].set_yscale("log")
axs[1, 0].legend(loc = 'upper left')

axs[1, 1].scatter(na_mean_abun, na_prev, alpha=0.5, label='North America')
axs[1, 1].scatter(africa_mean_abun, africa_prev, alpha=0.3, label='Africa')
axs[1, 1].set_xlabel('Mean Relative Abundance')
axs[1, 1].set_ylabel('Prevalence')
axs[1, 1].set_xlim(10e-9,20e-2)
axs[1, 1].set_ylim(0.002,1)
axs[1, 1].set_yscale("log")
axs[1, 1].set_xscale("log")
axs[1, 1].legend()
plt.savefig('C:/Users/sarah/Garud Lab/BIG_2021_microbiome_evolution/plots/macroecology_plot.png', bbox_inches='tight', dpi=600)
