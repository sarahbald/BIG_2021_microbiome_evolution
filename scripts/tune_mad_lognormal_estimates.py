import numpy as np
import math
import sympy as sp
import scipy.special
import pickle
import matplotlib.pyplot as plt


def Klogn(emp_mad, c, mu0, s0):
    # This function estimates the parameters (mu, s) of the lognormal distribution of K
    m1 = np.mean(np.log(emp_mad[emp_mad>c]))
    m2 = np.mean(np.log(emp_mad[emp_mad>c])**2)
    xmu = sp.Symbol('xmu')
    xs = sp.Symbol('xs')
    eq1 = - m1 + xmu + np.sqrt(2/math.pi)*xs*(sp.exp(-((np.log(c)-xmu)**2)/2/(xs**2))/(sp.erfc((np.log(c)-xmu)/np.sqrt(2)/xs)))
    eq2 = - m2 + xs**2 + m1*xmu+np.log(c)*m1-xmu*np.log(c)
    sol = sp.nsolve([eq1,eq2],[xmu,xs],[mu0,s0])

    return(float(sol[0]),float(sol[1]))


data_directory = "C:/Users/sarah/Garud Lab/"
filename = data_directory+"species_mean_relative_abundance.dat"
with open(filename, 'rb') as handle:
        mean_abundance_dict = pickle.load(handle)
        
af_abun_all = []
na_abun_all = []
for species_abundance in mean_abundance_dict.keys():
    af_abun_all.append(mean_abundance_dict[species_abundance]['Africa'])
    na_abun_all.append(mean_abundance_dict[species_abundance]['North America'])   
midas_af_mad = np.array(af_abun_all)
midas_hmp_mad = np.array(na_abun_all)


#test c-values that span from the minimum to the median of the aundance distribution 8.2/7.5
min_af_abun = np.log10(min(midas_af_mad))
min_na_abun = np.log10(min(midas_hmp_mad))
c_values_na = np.logspace(-8.2, -6.2, num=20, base=10)
c_values_af = np.logspace(-7.5, -5.5, num=20, base=10)
#initial_mu = [-10, -1, -1e-1, 1e-1, 1, 1e1, 1e2]
initial_mu = [-10]
#initial_s0 = [1e-1, 1, 1e1]
initial_s0 = [1]
mu_estimates_na = []
mu_estimates_af = []
sigma_estimates_na = []
sigma_estimates_af = []
initial_conditions = []
for j,s0 in enumerate(initial_s0):
    for i,m0 in enumerate(initial_mu):
        for k,val in enumerate(c_values_na):
            c_na = val
            c_af = c_values_af[k]
            (mu_na,sigma_na) = Klogn(np.array(midas_hmp_mad), c_na, m0, s0)
            mu_estimates_na.append(mu_na)
            sigma_estimates_na.append(sigma_na)
            (mu_af,sigma_af) = Klogn(np.array(midas_af_mad), c_af, m0, s0)
            mu_estimates_af.append(mu_af)
            sigma_estimates_af.append(sigma_af)


def calc_log_likelihood(mu, sigma, mad, c):
    #plug each observed abundance into EQ to get probability vector
    probabilities = []
    for abun in mad:
        if abun > c:
            prob = (1/(np.sqrt(2*np.pi*(sigma**2))*abun))*np.exp((-1*(np.log(abun) - mu)**2)/(2*(sigma**2)))
            probabilities.append(prob)
    #take ln of probability vector and sum to get log likelihood
    log_prob = [np.log(p) for p in probabilities]
    log_likelihood = sum(log_prob)/len(probabilities)
    return log_likelihood

log_l_na = []
log_l_af = []
for i, val in enumerate(c_values_na):
    log_l_na.append(calc_log_likelihood(mu_estimates_na[i], sigma_estimates_na[i], midas_hmp_mad, val))
    log_l_af.append(calc_log_likelihood(mu_estimates_na[i], sigma_estimates_na[i], midas_af_mad, val))
    
#select c values that give lowest likelihood
min_index_na = log_l_na.index(max(log_l_na))
min_index_af = log_l_af.index(max(log_l_af))
opt_c_na = c_values_na[min_index_na]
opt_c_af = c_values_af[min_index_af]
opt_mu_na = mu_estimates_na[min_index_na]
opt_mu_af = mu_estimates_af[min_index_af]
opt_sigma_na = sigma_estimates_na[min_index_na]
opt_sigma_af = sigma_estimates_af[min_index_af]

#plot best estimates
na_mad_log = np.log(midas_hmp_mad)
na_hist_mad, na_bin_edges_mad = np.histogram(na_mad_log, density=True, bins=20)
na_bins_mean_mad = [0.5 * (na_bin_edges_mad[i] + na_bin_edges_mad[i+1]) for i in range(0, len(na_bin_edges_mad)-1 )]
na_prob_to_plot = [sum( (na_mad_log>=na_bin_edges_mad[i]) & (na_mad_log<=na_bin_edges_mad[i+1])  ) / (len(na_mad_log)*1.0) for i in range(0, len(na_bin_edges_mad)-1 )]
na_bins_mean_plot = [np.exp(na_bins_mean_mad[i]) for i in range(len(na_bins_mean_mad))]

af_mad_log = np.log(midas_af_mad)
af_hist_mad, af_bin_edges_mad = np.histogram(af_mad_log, density=True, bins=20)
af_bins_mean_mad = [0.5 * (af_bin_edges_mad[i] + af_bin_edges_mad[i+1]) for i in range(0, len(af_bin_edges_mad)-1 )]
af_prob_to_plot = [sum( (af_mad_log>=af_bin_edges_mad[i]) & (af_mad_log<=af_bin_edges_mad[i+1])  ) / (len(af_mad_log)*1.0) for i in range(0, len(af_bin_edges_mad)-1 )]
af_bins_mean_plot = [np.exp(af_bins_mean_mad[i]) for i in range(len(af_bins_mean_mad))]


# range of abundances to plot
x=np.arange(-14,-1,step=0.1)
fig, ax = plt.subplots()

ax.scatter(na_bins_mean_plot, na_prob_to_plot, alpha=0.5, s=30)
ax.scatter(af_bins_mean_plot, af_prob_to_plot, alpha=0.5, s=30)

#for i in range(5,10):
#    opt_sigma_na = sigma_estimates_na[i]
#    opt_mu_na = mu_estimates_na[i]
#    opt_c_na = c_values_na[i]
#    opt_sigma_af = sigma_estimates_af[i]
#    opt_mu_af = mu_estimates_af[i]
#    opt_c_af = c_values_af[i]
    
ax.plot(np.exp(x), np.sqrt(2/math.pi)/opt_sigma_na *np.exp(-(x-opt_mu_na)**2 /2/(opt_sigma_na**2))/scipy.special.erfc((np.log(opt_c_na)-opt_mu_na)/np.sqrt(2)/opt_sigma_na), label = 'NA {:0.2e}'.format(opt_c_na))
ax.plot(np.exp(x), np.sqrt(2/math.pi)/opt_sigma_af *np.exp(-(x-opt_mu_af)**2 /2/(opt_sigma_af**2))/scipy.special.erfc((np.log(opt_c_af)-opt_mu_af)/np.sqrt(2)/opt_sigma_af), label = 'AF {:0.2e}'.format(opt_c_af))
#plot truncated lognormal fitted for K>log(c)
#ax.plot(numpy.log10(numpy.exp(x)),numpy.sqrt(2/math.pi)/sigma *numpy.exp(-(x-mu)**2 /2/(sigma**2))/scipy.special.erfc((numpy.log(c)-mu)/numpy.sqrt(2)/sigma),'-k')
#ax.plot(numpy.exp(x), numpy.sqrt(2/math.pi)/sigma *numpy.exp(-(x-mu)**2 /2/(sigma**2))/scipy.special.erfc((numpy.log(c)-mu)/numpy.sqrt(2)/sigma),'-k')
ax.set_xscale('log', basex=10)
ax.set_xlabel('Mean Abundance Distribution')
ax.set_yscale('log', basey=10)
ax.set_xlim(10**(-9), 0.5)
ax.set_ylim(10**(-3), 0.5)
ax.legend()
plt.savefig('C:/Users/sarah/Garud Lab/BIG_2021_microbiome_evolution/plots/mad_lognormal_estimates.png', dpi=600)

#Calculate covariance between African and North American mad including zeros
data_directory = "C:/Users/sarah/Garud Lab/"

filename = data_directory+"species_mean_relative_abundance_zeros.dat"
with open(filename, 'rb') as handle:
        species_abun_dict_zeros = pickle.load(handle)
na_mean_species_abun_zeros = []
af_mean_species_abun_zeros = []
for species in species_abun_dict_zeros.keys():
    if species_abun_dict_zeros[species]['Africa'] == 0.0 and species_abun_dict_zeros[species]['North America'] == 0.0:
        continue
    af_mean_species_abun_zeros.append(species_abun_dict_zeros[species]['Africa'])
    na_mean_species_abun_zeros.append(species_abun_dict_zeros[species]['North America'])
    
cov_mat = np.cov(na_mean_species_abun_zeros, af_mean_species_abun_zeros)
log_cov_mat = np.cov(na_mad_log, af_mad_log)

#save parameters into a dict
lognormal_fit_params_dict = {}
lognormal_fit_params_dict['Africa'] = {}
lognormal_fit_params_dict['North America'] = {}
lognormal_fit_params_dict['Africa']['mu'] = opt_mu_af
lognormal_fit_params_dict['Africa']['sigma'] = opt_sigma_af
lognormal_fit_params_dict['Africa']['cutoff'] = opt_c_af
lognormal_fit_params_dict['Africa']['mad'] = midas_af_mad
lognormal_fit_params_dict['North America']['mu'] = opt_mu_na
lognormal_fit_params_dict['North America']['sigma'] = opt_sigma_na
lognormal_fit_params_dict['North America']['cutoff'] = opt_c_na
lognormal_fit_params_dict['North America']['mad'] = midas_hmp_mad
lognormal_fit_params_dict['Species Abundance Dict'] = mean_abundance_dict
lognormal_fit_params_dict['Covariance Matrix'] = cov_mat
lognormal_fit_params_dict['Log-Cov Matrix'] = log_cov_mat


with open(data_directory+'lognormal_params_dict.pickle', 'wb') as handle:
    pickle.dump(lognormal_fit_params_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)