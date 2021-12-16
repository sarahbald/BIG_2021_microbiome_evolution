from __future__ import division
import math
import numpy
import sympy as sp
import scipy
import matplotlib.pyplot as plt



import calculate_species_relative_abundance
mean_abundance_dict = calculate_species_relative_abundance.load_species_mean_abun_dict()
af_abun_all = []
na_abun_all = []
for species_abundance in mean_abundance_dict.keys():
    af_abun_all.append(mean_abundance_dict[species_abundance]['Africa'])
    na_abun_all.append(mean_abundance_dict[species_abundance]['North America'])   
af_mad = numpy.array(af_abun_all)
na_mad = numpy.array(na_abun_all)

emp_mad = af_mad


# empirical mean abundance distribution
#emp_mad = numpy.asarray([0.00001, 0.00001, 0.0001, 0.001, 0.01, 0.1, 0.0001, 0.001, 0.000001, 0.000001, 0.000001, 0.000001,0.000001, 0.000001,0.000001 ,0.000001 ,0.000001 ])

def Klogn(emp_mad, c, mu0=-19,s0=5):
    # This function estimates the parameters (mu, s) of the lognormal distribution of K
    m1 = numpy.mean(numpy.log(emp_mad[emp_mad>c]))
    m2 = numpy.mean(numpy.log(emp_mad[emp_mad>c])**2)
    xmu = sp.symbols('xmu')
    xs = sp.symbols('xs')
    eq1 = - m1 + xmu + numpy.sqrt(2/math.pi)*xs*(sp.exp(-((numpy.log(c)-xmu)**2)/2/(xs**2))/(sp.erfc((numpy.log(c)-xmu)/numpy.sqrt(2)/xs)))
    eq2 = - m2 + xs**2 + m1*xmu+numpy.log(c)*m1-xmu*numpy.log(c)
    #print(eq1)
    sol = sp.nsolve([eq1,eq2],[xmu,xs],[mu0,s0])

    return(float(sol[0]),float(sol[1]))


#hist_na, bin_edges_na = numpy.histogram(na_mad, density=True, bins=numpy.logspace(-3,-0.5, 25))
#bins_mean_na = [0.5 * (bin_edges_na[i] + bin_edges_na[i+1]) for i in range(0, len(bin_edges_na)-1 )]

#hist_af, bin_edges_af = numpy.histogram(af_mad, density=True, bins=numpy.logspace(-3,-0.5, 25))
#bins_mean_af = [0.5 * (bin_edges_af[i] + bin_edges_af[i+1]) for i in range(0, len(bin_edges_af)-1 )]

mad_log10 = numpy.log10(emp_mad)
hist_mad, bin_edges_mad = numpy.histogram(mad_log10, density=True, bins=20)
bins_mean_mad = [0.5 * (bin_edges_mad[i] + bin_edges_mad[i+1]) for i in range(0, len(bin_edges_mad)-1 )]
prob_to_plot = [sum( (mad_log10>=bin_edges_mad[i]) & (mad_log10<=bin_edges_mad[i+1])  ) / len(mad_log10) for i in range(0, len(bin_edges_mad)-1 )]
bins_mean_plot = [numpy.exp(bins_mean_mad[i]) for i in range(len(bins_mean_mad))]

# range of abundances to plot
x=numpy.arange(-14,-1,step=0.1)
fig, ax = plt.subplots()

ax.scatter(bins_mean_plot, prob_to_plot, alpha=0.5, s=30, label = "North America")
#ax.scatter(bins_mean_af, hist_af, alpha=0.5, s=30, label = "Africa")
#c_values = [10e-8, 10e-7, 10e-6]
c_values = [10e-7]
#for value in numpy.linspace(numpy.quantile(emp_mad, 0.25), numpy.quantile(emp_mad, 0.5), 5):
for value in c_values:
    c = value
    (mu,sigma) = Klogn(emp_mad, c)
    ax.plot(numpy.exp(x), numpy.sqrt(2/math.pi)/sigma *numpy.exp(-(x-mu)**2 /2/(sigma**2))/scipy.special.erfc((numpy.log(c)-mu)/numpy.sqrt(2)/sigma),'-k')

#plot truncated lognormal fitted for K>log(c)
#ax.plot(numpy.log10(numpy.exp(x)),numpy.sqrt(2/math.pi)/sigma *numpy.exp(-(x-mu)**2 /2/(sigma**2))/scipy.special.erfc((numpy.log(c)-mu)/numpy.sqrt(2)/sigma),'-k')
#ax.plot(numpy.exp(x), numpy.sqrt(2/math.pi)/sigma *numpy.exp(-(x-mu)**2 /2/(sigma**2))/scipy.special.erfc((numpy.log(c)-mu)/numpy.sqrt(2)/sigma),'-k')
ax.set_xscale('log', basex=10)
ax.set_yscale('log', basey=10)
ax.set_xlim(10**(-3.5), 1)
plt.savefig('C:/Users/sarah/Garud Lab/plots/mad_estimates_na.png', dpi=600)
