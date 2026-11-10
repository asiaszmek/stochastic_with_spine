import os
import sys
import h5py
import numpy as np
from scipy.fft import fft, fftfreq, fftshift
import matplotlib.pyplot as plt
import utility_functions as ut

plt.rcParams['text.usetex'] = True
colors =  {
    "1.2": 'tab:blue',
    "2.4": 'tab:purple',
    "6.0": 'tab:green'
}
names_dict = {
    "100%\n100%\n0%" :
    "model_RyRCaM_simple_SERCA_SOCE_tubes_diam_%s_um_10_um_dendrite.h5",
    "100%\n0\n0%":
    "model_noRyR_simple_SERCA_SOCE_tubes_diam_%s_um_10_um_dendrite.h5",
    "100%\n0\n100%%":
    "model_RyR_simple_SERCA_SOCE_tubes_diam_%s_um_2_um_dendrite.h5",
    "80%\n0\n100%":
    "model_RyR_simple_SERCA_0.8_PMCA_tubes_diam_%s_um_10_um_dendrite.h5",
    # "80%\n0\n100%%":
    # "model_RyR_simple_SERCA_SOCE_0.8_PMCA_tubes_diam_%s_um_2_um_dendrite.h5",
    "80%\n50%\n50%":
    "model_RyR_RyRCaM_0.8_PMCA_simple_SERCA_tubes_diam_%s_um_2_um_dendrite.h5",
    "80%\n100%\n100%":
    "model_2x_RyR_RyRCaM_0.8_PMCA_simple_SERCA_tubes_diam_%s_um_2_um_dendrite.h5",
    # "80%\n200%\n200%":
    # "model_4x_RyR_RyRCaM_0.8_PMCA_simple_SERCA_tubes_diam_%s_um_2_um_dendrite.h5",
}


dend_diam = ["1.2", "2.4", "6.0"]
output = "__main__"


def adjust_axes(ax):
    mini = min([min(x.get_ylim()) for x in ax])
    maxi = max([max(x.get_ylim()) for x in ax])
    for x in ax:
      
        x.set_ylim([mini, maxi])


if __name__ == "__main__":
    data_dir = os.path.join("..", "stacked_ER")
    x_labels = []
    x_labels_pmca = []
    fig_m_ca, ax_m_ca = plt.subplots(1, len(dend_diam),
                                           figsize=(len(dend_diam)*5, 5))

    for i, d in enumerate(dend_diam):
        means = []
        stds = []
        x_labels = []
        
        for key, fname in names_dict.items():
            path = os.path.join(data_dir, fname % d)
            my_file = h5py.File(path)
            grid_list = ut.get_grid_list(my_file)
            conc_list = []
            conc_std = []
            for trial in ["trial0", "trial1", "trial2", "trial3"]:
                data = ut.get_populations(my_file, trial=trial,
                                          output=output)
                specie_idx = ut.get_all_species(my_file,
                                                output=output).index("Ca")
                volume = sum(list(ut.region_volumes(my_file).values()))
                
                tot_conc =  ut.nano_molarity(data[:, :,
                                                  specie_idx].sum(axis=1),
                                         volume)
                conc_list.append(tot_conc.mean())
                conc_std.append(tot_conc.std())

            means.append(np.mean(conc_list))
            stds.append(sum([s**2 for s in conc_std])**0.5/2)
            x_labels.append(key)
    
 
        ax_m_ca[i].errorbar(x=x_labels, y=means, yerr=stds,
                        color=colors[d],
                        marker="o", linestyle="")
        ax_m_ca[i].set_title(r"diam %s  $\mathrm{\mu m}$" % d)
        if i:
            ax_m_ca[i].set_yticklabels([])
    ax_m_ca[0].set_ylabel(r"mean $\mathrm{Ca_i}$ (nM)",
                          fontsize=15)
    legend = "PMCA kcat\nRyR2CaM\nRyR2" #\nSOCE"
    adjust_axes(ax_m_ca)
    ax_m_ca[0].text(-2.5, min(ax_m_ca[0].get_ylim())
                    -(max(ax_m_ca[0].get_ylim())
                      -min(ax_m_ca[0].get_ylim()))*0.1698, legend,
                    horizontalalignment='left', fontsize=15)
    #ax_m_ca[0].set_xticklabels(x_labels, rotation=90)
    for ax in ax_m_ca:
        ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=15)

    fig_m_ca.savefig("mean_basal_ca.png", dpi=100,
                 bbox_inches="tight")
