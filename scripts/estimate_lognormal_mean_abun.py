import calculate_species_relative_abundance
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
from scipy.special import erf
from seaborn import kdeplot

mean_abundance_dict = calculate_species_relative_abundance.load_species_mean_abun_dict()

af_abun_all = []
na_abun_all = []
for species_abundance in mean_abundance_dict.keys():
    af_abun_all.append(mean_abundance_dict[species_abundance]['Africa'])
    na_abun_all.append(mean_abundance_dict[species_abundance]['North America'])
    
af_abun_all = np.array(af_abun_all)
na_abun_all = np.array(na_abun_all)
    
af_log_abun = [np.log10(mean) for mean in af_abun_all]
na_log_abun = [np.log10(mean) for mean in na_abun_all]
af_log_abun = np.array(af_log_abun)
na_log_abun = np.array(na_log_abun)

standardized_africa_abun = (af_log_abun - af_log_abun.mean())/af_log_abun.std()
standardized_na_abun = (na_log_abun - na_log_abun.mean())/na_log_abun.std()

fig, axs = plt.subplots(1, 2, figsize=(15,7))
n_af, bins_af, patches = axs[1].hist(np.log10(af_abun_all), bins = 40)
n_na, bins_na, patches = axs[0].hist(np.log10(na_abun_all), bins = 40)
mode_na = (bins_na[n_na.argmax()] + bins_na[n_na.argmax()+1])/2
mode_af = (bins_af[n_af.argmax()] + bins_af[n_af.argmax()+1])/2

n = float(na_abun_all.shape[0])

na_m1 = (1/n)*na_log_abun.sum()
na_m2 = (1/n)*((na_log_abun)**2).sum()

af_m1 = (1/n)*af_log_abun.sum()
af_m2 = (1/n)*((af_log_abun)**2).sum()

c_na = np.exp(mode_na)
c_af = np.exp(mode_af)

def get_estimates(init, c):
    sigma = init[0]
    mu = init[1]
    a = -m1 + mu + ((np.sqrt(2/np.pi)*sigma*np.exp(-1*((np.log(c)-mu)**2)/(2*(sigma**2))))/erf((np.log(c)-mu)/(np.sqrt(2)*sigma)))
    b = -m2 + (sigma**2) + m1*mu + c*m1 - mu*c
    return np.array([a,b])

init = [1,1]
m1 = na_m1
m2 = na_m2
c = c_na
na_dist = fsolve(get_estimates, init, c)

m1 = af_m1
m2 = af_m2
c = c_af
af_dist = fsolve(get_estimates, init, c)

na_var = na_dist[0]**2
na_mu = na_dist[1]
af_var = af_dist[0]**2
af_mu = af_dist[1]

x_na = np.linspace(-8,-1,50)
y_na = np.array(1/(np.sqrt(2*np.pi*na_var)*np.exp(x_na))*np.exp(-1*(((x_na - na_mu)**2)/(2*na_var))))
#y_na = (np.sqrt(2)/(np.sqrt(np.pi*na_var)*np.exp(x_na)))*np.heaviside(np.exp(x_na) - c_na, 1)*(np.exp(-((x_na-na_mu)**2)/(2*na_var)))/(erf((np.log(c_na)-na_mu)/(np.sqrt(2)*np.sqrt(na_var))))
x_af = np.linspace(-8,-1,50)
y_af = np.array(1/(np.sqrt(2*np.pi*af_var)*np.exp(x_af))*np.exp(-1*(((x_af - af_mu)**2)/(2*af_var))))
#y_af = (np.sqrt(2)/(np.sqrt(np.pi*af_var)*np.exp(x_af)))*np.heaviside(np.exp(x_af) - c_af, 1)*(np.exp(-(x_af-af_mu)**2/(2*af_var)))/(erf((np.log(c_af)-af_mu)/(np.sqrt(2)*np.sqrt(af_var))))

rescaled_pred_na = (y_na - y_na.mean())/y_na.std()
rescaled_pred_af = (y_af - y_af.mean())/y_af.std()
hist_na, bin_edges_na = np.histogram(rescaled_pred_na, density=True, bins=25)
bins_mean_na = [0.5 * (bin_edges_na[i] + bin_edges_na[i+1]) for i in range(0, len(bin_edges_na)-1 )]
hist_af, bin_edges_af = np.histogram(rescaled_pred_af, density=True, bins=25)
bins_mean_af = [0.5 * (bin_edges_af[i] + bin_edges_af[i+1]) for i in range(0, len(bin_edges_af)-1 )]


hist_na_2, bin_edges_na_2 = np.histogram(standardized_na_abun, density=True, bins=25)
bins_mean_na_2 = [0.5 * (bin_edges_na_2[i] + bin_edges_na_2[i+1]) for i in range(0, len(bin_edges_na_2)-1 )]

hist_af_2, bin_edges_af_2 = np.histogram(standardized_africa_abun, density=True, bins=25)
bins_mean_af_2 = [0.5 * (bin_edges_af_2[i] + bin_edges_af_2[i+1]) for i in range(0, len(bin_edges_af_2)-1 )]

plt.clf()
plt.figure()
plt.scatter(bins_mean_na_2, hist_na_2, alpha=0.5, s=30, label = "North America")
plt.scatter(bins_mean_af_2, hist_af_2, alpha=0.5, s=30, label = "Africa")

plt.xlabel("Rescaled Log of Mean Abundances")
plt.ylabel("Density")
#axs[0, 0].set_yscale("log")
plt.legend()
plt.savefig('C:/Users/sarah/Garud Lab/plots/mad_densityhist.png', dpi=600)
