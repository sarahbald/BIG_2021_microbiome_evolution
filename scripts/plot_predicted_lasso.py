#importance = numpy.abs(reg.coef_)
importance = numpy.abs(reg.coef_)

X_to_keep = X[:,importance != 0]

# keep the slopes that are not zero
kegg_ids_to_keep = kegg_ids[importance != 0]

observed_delta_prevalence_all = []
predicted_delta_prevalence_all = []

for k_idx, k in enumerate(kegg_ids_to_keep):

    k_column = X_to_keep[:,k_idx]

    observed_delta_prevalence_k = y[k_column>0]

    k_mean_delta_prevalence = numpy.mean(observed_delta_prevalence_k)

    predicted_delta_prevalence_k = predicted_delta_prev[k_column>0]

    observed_delta_prevalence_all.append(observed_delta_prevalence_k[0])
    predicted_delta_prevalence_all.append(predicted_delta_prevalence_k[0])

    print k, pathway_annotation_dict[k], importance[importance != 0][k_idx], k_mean_delta_prevalence



observed_delta_prevalence_all = numpy.asarray(observed_delta_prevalence_all)
predicted_delta_prevalence_all = numpy.asarray(predicted_delta_prevalence_all)



ss_tot = sum( (observed_delta_prevalence_all - numpy.mean(observed_delta_prevalence_all))**2 )
ss_res = sum( (observed_delta_prevalence_all - predicted_delta_prevalence_all)**2 )

r2 = 1 - (ss_res / ss_tot)

print(r2)



fig, ax = plt.subplots(figsize=(4,4))

ax.scatter(observed_delta_prevalence_all, predicted_delta_prevalence_all, c='royalblue', alpha=0.7)

min_ = min([min(observed_delta_prevalence_all), min(predicted_delta_prevalence_all)])
if min_ <0:
    min_ = min_*1.7
else:
    min_ = min_ * 0.8
max_ = max([max(observed_delta_prevalence_all), max(predicted_delta_prevalence_all)]) * 1.1


ax.set_xlim([min_, max_])
ax.set_ylim([min_, max_])

ax.plot([min_, max_],[min_, max_], ls='--', c='k')

ax.xaxis.set_tick_params(labelsize=7)
ax.yaxis.set_tick_params(labelsize=7)





fig.tight_layout()
#fig.subplots_adjust(hspace=0.2)
fig.savefig("%slasso_observed_predicted.png" % config.analysis_directory, format='png', bbox_inches = "tight", pad_inches = 0.4, dpi = 600)
plt.close()
