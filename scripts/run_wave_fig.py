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
                                 len(volume)//102)).mean(axis=(1,2))
              
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
    fig_wave_f, ax_wave_f = plt.subplots(1, 1,
                                         figsize=(1*5, 5))
    for i, d in enumerate(dend_diam):
        f_mean = []
        f_std = []
        xlabels = []
        for key, fname in names_dict.items():
            path = os.path.join(data_dir, fname % d)
            my_file = h5py.File(path)
            grid_list = ut.get_grid_list(my_file)
            peak_no = []
            length = []
            for trial in ["trial0", "trial1", "trial2", "trial3", "trial4"]:
                new_conc = get_ca_conc(my_file, trial)
                if new_conc is None:
                    continue
                peak = 0
                where = np.where(new_conc > 2.5*76)
                if len(where[0]):
                    prev = where[0][0]
                    peak = 1
                    for idx in where[0][1:]:
                        if idx >= prev +1 and idx <= prev + 20:
                            prev = idx
                        else:
                            peak += 1
                            prev = idx
                peak_no.append(peak)
                length.append(len(new_conc)*0.2)
            xlabels.append(key)
            m_f = (np.array(peak_no)/np.array(length)).mean()
            std_f = (np.array(peak_no)/np.array(length)).std()/(len(peak_no))**0.5
            f_mean.append(m_f)
            f_std.append(std_f)

        ax_wave_f.errorbar(x=xlabels, y=f_mean,
                            yerr=f_std,
                            color=colors[d],
                            marker="o", linestyle="",
                            label=r"%s $\mathrm{\mu m}$ diam" % d)

                                
    ax_wave_f.set_ylabel(r"$\mathrm{Ca_i}$ wave frequency (Hz)",
                          fontsize=15)
        
    ax_wave_f.set_xticklabels(xlabels)
    ax_wave_f.tick_params(axis='x', labelsize=15)
    ax_wave_f.tick_params(axis='y', labelsize=15)
    legend = "PMCA kcat\nRyR2CaM\nRyR2"#\nSOCE"
    ax_wave_f.legend()
   
    ax_wave_f.text(-1.25, min(ax_wave_f.get_ylim())
                    -(max(ax_wave_f.get_ylim())
                      -min(ax_wave_f.get_ylim()))*0.1698,
                           legend,
                    horizontalalignment='left', fontsize=15)
    

    fig_wave_f.savefig("mean_wave_frequency.png", dpi=100,
                       bbox_inches="tight")
