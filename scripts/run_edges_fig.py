import os
import sys
import h5py

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
import utility_functions as ut


plt.rcParams['text.usetex'] = True
colors =  {
    "1.2": 'tab:blue',
    "2.4": 'tab:purple',
    "6.0": 'tab:green'
}

names_dict = {
   
    # "normal PMCA + no RyR2": "model_noRyR_simple_SERCA_SOCE_tubes_diam_%s_um_10_um_dendrite.h5",
    # "normal PMCA + RyR2": "model_RyR_simple_SERCA_SOCE_tubes_diam_%s_um_2_um_dendrite.h5",
    # "low PMCA, RyR2CaM": "model_RyRCaM_simple_SERCA_SOCE_0.8_PMCA_tubes_diam_%s_um_10_um_dendrite.h5",
    "80%\n100%\n0%":
    "model_RyRCaM_simple_SERCA_SOCE_0.8_PMCA_tubes_diam_%s_um_10_um_dendrite.h5",
    "80%\n50%\n50%": "model_RyR_RyRCaM_0.8_PMCA_simple_SERCA_tubes_diam_%s_um_2_um_dendrite.h5",
   "80%\n0\n100%": "model_RyR_simple_SERCA_0.8_PMCA_tubes_diam_%s_um_10_um_dendrite.h5",
    # "80%\n0\n100%\n100%":
    # "model_RyR_simple_SERCA_SOCE_0.8_PMCA_tubes_diam_%s_um_2_um_dendrite.h5",
    "80%\n100%\n100%": "model_2x_RyR_RyRCaM_0.8_PMCA_simple_SERCA_tubes_diam_%s_um_2_um_dendrite.h5",
     "80%\n200%\n200%": "model_4x_RyR_RyRCaM_0.8_PMCA_simple_SERCA_tubes_diam_%s_um_2_um_dendrite.h5"
    
}

dend_diam = ["1.2", "2.4", "6.0"]
output = "__main__"

        
def find_coords(P0, P1):
    return [[min(P0[0], P1[0]), min(P0[1], P1[1])],
            [max(P0[0], P1[0]), max(P0[1], P1[1])]]


def is_y(miny, maxy, old_miny, old_maxy):
    if miny <= old_maxy and miny >= old_miny :
        return True
    if maxy <= old_maxy and maxy >= old_miny:
        return True


def get_ca_conc(my_file, trial, output=output):
    try:
        data = ut.get_populations(my_file, trial=trial,
                                  output=output)
    except KeyError:
        return
                
    specie_idx = ut.get_all_species(my_file,
                                    output=output).index("Ca")

    volume = np.array([g[12] for g in grid_list])
    conc =  ut.nano_molarity(data[:, :, specie_idx],
                                         volume)
    new_conc = np.reshape(conc, (data.shape[0],
                                 102,
                                 len(volume)//102)).mean(axis=2).T
              
    return new_conc


def find_clusters(idx_array):
    clusters = [[[idx_array[0][0], idx_array[1][0]],
                 [idx_array[0][0], idx_array[1][0]]]]
    indices = zip(idx_array[0][1:], idx_array[1][1:])
    for i in range(len(idx_array[0][1:])):
        new_x, new_y = idx_array[0][i+1], idx_array[1][i+1]
        for  j, c in enumerate(clusters):
            if c[0][0]-10 <= new_x <= c[1][0]+10:
                if c[0][1]-4 <= new_y <= c[1][1]+4:
                    new_minx = min([c[0][0], c[1][0], new_x])
                    new_maxx = max([c[0][0], c[1][0], new_x])
                    new_miny = min([c[0][1], c[1][1], new_y])
                    new_maxy = max([c[0][1], c[1][1], new_y])
                    clusters[j] = [[new_minx, new_miny],
                                   [new_maxx, new_maxy]]
                    break
        else:
            clusters.append([[new_x, new_y], [new_x, new_y]])
    return clusters


def purge_clusters(clusters):
    new_clusters = []
    for i, c in enumerate(clusters):
        new_p1, new_p2 = c
        for j, nc in enumerate(new_clusters):
            p1, p2 = nc
            if (p1[0] <= new_p1[0] <= p2[0]
                or p1[0] <= new_p2[0] <= p2[0]):
                if (p1[1] <= new_p1[1] <= p2[1]
                    or p1[1]<= new_p2[1] <=p2[1]):
                    c_minx = min([p1[0], p2[0],
                                  new_p1[0], new_p2[0]])
                    c_maxx = max([p1[0], p2[0],
                                  new_p1[0], new_p2[0]])
                    c_miny = min([p1[1], p2[1],
                                  new_p1[1], new_p2[1]])
                    c_maxy = max([p1[1], p2[1],
                                  new_p1[1], new_p2[1]])
                    new_clusters[j] = [[c_minx, c_miny],
                                       [c_maxx, c_maxy]]
                    break
        else:   
            if  c[1][1] - c[0][1]>=5 and c[1][0] - c[0][0]>=5:
                new_clusters.append(c)
    if len(new_clusters):
        new_clusters = sorted(new_clusters, key=lambda x:x[0][0])
    return new_clusters


def get_width(p01, p02):
    return (p02[1] - p01[1])*0.2


def get_len(p01, p02):
    return (p02[0]-p01[0])/2+0.5


def max_amp(conc, p01, p02):
    return conc[p01[0]:p02[0]+1, p01[1]:p02[1]+1].max()


def peak_init(conc, p01):
    return conc[:, p01[1]].argmax()/2


def get_ipi(clusters):
    peak_dist = []
    for j, nc in enumerate(new_clusters):
        p1, p2 = nc
        if j:
            o_p1, o_p2 = new_clusters[j-1]
            new_t = (p2[1]+p1[1])/2*0.2
            old_t = (o_p2[1]+o_p1[1])/2*0.2
            peak_dist.append(abs(new_t - old_t))
    return peak_dist


def adjust_axes(ax):
    mini = min([min(x.get_ylim()) for x in ax])
    maxi = max([max(x.get_ylim()) for x in ax])
    for x in ax:
      
        x.set_ylim([mini, maxi])



if __name__ == "__main__":
    data_dir = os.path.join("..", "stacked_ER")
    x_labels = list(names_dict.keys())
    fig_m_peak_amp, ax_m_peak_amp = plt.subplots(1, 1,
                                                 figsize=(1*5,
                                                          5))
    fig_m_peak_width, ax_m_peak_width = plt.subplots(1, 1,
                                                     figsize=(1*5,
                                                              5))
    fig_m_peak_len, ax_m_peak_len = plt.subplots(1, 1,
                                                 figsize=(1*5,
                                                              5)) 
    fig_m_peak_dist, ax_m_peak_dist = plt.subplots(1, 1,
                                                   figsize=(1*5,
                                                            5))
    fig_m_peak_init, ax_m_peak_init = plt.subplots(1, len(dend_diam),
                                                   figsize=(len(dend_diam)*5,
                                                            5))
   
    ax_m_peak_init = ax_m_peak_init
    for i, d in enumerate(dend_diam):
        no_mean = []
        no_std = []
        amp_mean = []
        amp_std = []
        width_mean = []
        width_std = []
        len_mean = []
        len_std = []
        dist_mean = []
        dist_std = []
        xlabels = []
        for key, fname in names_dict.items():
            path = os.path.join(data_dir, fname % d)
            print(path)
            my_file = h5py.File(path)
            grid_list = ut.get_grid_list(my_file)
            peak_amp = []
            peak_width = []
            peak_len = []
            peak_dist = []
            peak_start = []
            for trial in ["trial0", "trial1", "trial2", "trial3", "trial4"]:
                new_conc = get_ca_conc(my_file, trial)
                if new_conc is None:
                    continue
                where = np.where(new_conc > 2.5*76)
              
                if not len(where[0]):
                    peak_dist.append(0)
                    peak_width.append(0)
                    peak_len.append(0)
                    peak_amp.append(0)
                    peak_start.append(0)
                    continue

                my_clusters = find_clusters(where)
                new_clusters = purge_clusters(my_clusters)
                peak_dist.append(len(new_clusters)/((new_conc.shape[1]-1)*0.2))
                                 
                for j, nc in enumerate(new_clusters):
                    p1, p2 = nc
                    peak_width.append(get_width(p1, p2))
                    peak_len.append(get_len(p1, p2))
                    peak_amp.append(new_conc[p1[0]:p2[0]+1,
                                             p1[1]:p2[1]+1].max())
                    peak_start.append(new_conc[:, p1[1]].argmax()/2)
                else:
                    peak_width.append(0)
                    peak_len.append(0)
                    peak_amp.append(0)
                    peak_start.append(0)
                    
            
            xlabels.append(key)
            
            peak_amp = np.array(peak_amp)
            peak_width = np.array(peak_width)
            peak_len = np.array(peak_len)
            peak_dist = np.array(peak_dist)
            xlabels_for_peak_start = [key]*len(peak_start)
            amp_mean.append(peak_amp.mean())
            amp_std.append(peak_amp.std()/len(peak_amp)**.5)
            width_mean.append(peak_width.mean())
            width_std.append( peak_width.std()/len(peak_width)**.5)
            len_mean.append(peak_len.mean())
            len_std.append( peak_len.std()/len(peak_len)**.5)
            dist_mean.append(peak_dist.mean())
            dist_std.append(peak_dist.std()/len(peak_dist)**.5)
            
            if len(peak_start):
                ax_m_peak_init[i].plot(xlabels_for_peak_start, peak_start,
                                       color=colors[d],
                                       marker="o", linestyle="")
       

        ax_m_peak_amp.errorbar(x=xlabels, y=amp_mean,
                               yerr=amp_std,
                               color=colors[d],
                               marker="o", linestyle="",
                               label=r"%s $\mathrm{\mu m}$ diam" % d)
        ax_m_peak_width.errorbar(x=xlabels, y=width_mean,
                                 yerr=width_std,
                                 color=colors[d],
                                 marker="o", linestyle="",
                                 label=r"%s $\mathrm{\mu m}$ diam" % d)
        ax_m_peak_len.errorbar(x=xlabels, y=len_mean,
                               yerr=len_std,
                               color=colors[d],
                               marker="o", linestyle="",
                               label=r"%s  $\mathrm{\mu m}$ diam" % d)
        ax_m_peak_dist.errorbar(x=xlabels, y=dist_mean,
                                yerr=dist_std,
                                color=colors[d],
                                marker="o", linestyle="",
                                label=r"%s  $\mathrm{\mu m}$ diam" % d)
                                
        if i:
            ax_m_peak_init[i].set_yticklabels([])

        ax_m_peak_init[i].set_title(r"diam %s  $\mathrm{\mu m}$" % d)
        ax_m_peak_init[i].set_ylabel(r"initiation location $(\mathrm{\mu m})$",
                                     fontsize=15)
        ax_m_peak_init[i].set_xticklabels(xlabels)
    ax_m_peak_amp.set_ylabel(r"$\mathrm{Ca_i}$ amplitude (nM)",
                             fontsize=15)
    ax_m_peak_width.set_ylabel(r"$\mathrm{Ca_i}$ peak width (s)",
                               fontsize=15)
    ax_m_peak_len.set_ylabel(r"$\mathrm{Ca_i}$ peak length (s)",
                             fontsize=15)
    ax_m_peak_dist.set_ylabel(r"$\mathrm{Ca_i}$ peak frequency (Hz)",
                              fontsize=15)
        
    ax_m_peak_amp.set_xticklabels(xlabels)
    ax_m_peak_width.set_xticklabels(xlabels)
    ax_m_peak_len.set_xticklabels(xlabels)
    ax_m_peak_dist.set_xticklabels(xlabels)
    legend = "PMCA kcat\nRyR2CaM\nRyR2"#\nSOCE"
    ax_m_peak_amp.legend()
    ax_m_peak_width.legend()
    ax_m_peak_len.legend()
    ax_m_peak_dist.legend()
    ax_m_peak_amp.tick_params(axis='x', labelsize=15)
    ax_m_peak_amp.tick_params(axis='y', labelsize=15)
    ax_m_peak_width.tick_params(axis='x', labelsize=15)
    ax_m_peak_width.tick_params(axis='y', labelsize=15)
    ax_m_peak_len.tick_params(axis='x', labelsize=15)
    ax_m_peak_len.tick_params(axis='y', labelsize=15)
    ax_m_peak_dist.tick_params(axis='x', labelsize=15)
    ax_m_peak_dist.tick_params(axis='y', labelsize=15)
    for ax in ax_m_peak_init:
        ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=15)
        
    adjust_axes(ax_m_peak_init)
    ax_m_peak_init[0].text(-1.25, min(ax_m_peak_init[0].get_ylim())
                    -(max(ax_m_peak_init[0].get_ylim())
                      -min(ax_m_peak_init[0].get_ylim()))*0.1698,
                           legend,
                           horizontalalignment='left', fontsize=15)
    ax_m_peak_amp.text(-1.25, min(ax_m_peak_amp.get_ylim())
                    -(max(ax_m_peak_amp.get_ylim())
                      -min(ax_m_peak_amp.get_ylim()))*0.1698,
                           legend,
                    horizontalalignment='left', fontsize=15)
    ax_m_peak_width.text(-1.25, min(ax_m_peak_width.get_ylim())
                    -(max(ax_m_peak_width.get_ylim())
                      -min(ax_m_peak_width.get_ylim()))*0.1698,
                           legend,
                    horizontalalignment='left', fontsize=15)
    ax_m_peak_len.text(-1.25, min(ax_m_peak_len.get_ylim())
                    -(max(ax_m_peak_len.get_ylim())
                      -min(ax_m_peak_len.get_ylim()))*0.1698,
                           legend,
                    horizontalalignment='left', fontsize=15)
    ax_m_peak_dist.text(-1.25, min(ax_m_peak_dist.get_ylim())
                    -(max(ax_m_peak_dist.get_ylim())
                      -min(ax_m_peak_dist.get_ylim()))*0.1698,
                           legend,
                    horizontalalignment='left', fontsize=15)
    

    fig_m_peak_amp.savefig("mean_peak_amp.png", dpi=100,
                 bbox_inches="tight")
    fig_m_peak_dist.savefig("mean_peak_dist.png", dpi=100,
                            bbox_inches="tight")
    fig_m_peak_width.savefig("mean_peak_width.png", dpi=100,
                            bbox_inches="tight")
    fig_m_peak_len.savefig("mean_peak_len.png", dpi=100,
                            bbox_inches="tight")
    fig_m_peak_init.savefig("peak_init.png", dpi=100,
                            bbox_inches="tight")
