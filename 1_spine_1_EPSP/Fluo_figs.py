import os
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
from lxml import etree
from scipy.constants import Avogadro
NA = Avogadro*1e-23


def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w


def Parser():
    parser = argparse.ArgumentParser(description='Generation of figures')
    parser.add_argument('input', nargs='+',
                        help='input files')
    parser.add_argument('--t_init', default=3000, type=float,
                        help='Stimulation initiation in ms')
    parser.add_argument('--specie', default='CaFluo4FF',
                        help='Dye-bound specie')
   

    return parser


def get_regions(my_file):
    grid_list = get_grid_list(my_file)
    return sorted(list(set([get_key(grid) for grid in grid_list])))


def get_key(cell):
    if cell[18]:
        return cell[15].decode('utf-8') + '_' + cell[18].decode('utf-8')
    return cell[15].decode('utf-8')


def get_region_indices(my_file):
    grid_list = get_grid_list(my_file)
    region_ind = {}
    for idx, cell in enumerate(grid_list):
        key = get_key(cell)
        if key not in region_ind:
            region_ind[key] = []
        region_ind[key].append(idx)
    return region_ind


def region_volumes(my_file):
    if isinstance(my_file, str):
        my_file = h5py.File(my_file)
    grid_list = get_grid_list(my_file)
    regions = get_regions(my_file)
    volumes = {}
    for region in regions:
        volumes[region] = 0
    for cell in grid_list:
        key = get_key(cell)
        volumes[key] += float(cell[12])

    return volumes


def get_grid_list(My_file):
    return np.array(My_file['model']['grid'])


def nano_molarity(N, V):
    return 10 * N / V / NA


def pico_sd(N, S):
    return 10 * N / S / NA


def get_populations(my_file, trial='trial0', output='__main__'):
    return np.array(my_file[trial]['output'][output]['population'])


def get_times(My_file, trial, output):
    return np.array(My_file[trial]['output'][output]['times'])


def get_all_species(My_file, output="__main__"):
    return [s.decode('utf-8') for s in My_file['model']['output'][output]['species']]


def sum_indices(my_file, region_list, spines=["_sa1[0]"]):
    reg_indices = get_region_indices(my_file)
    sum_indices = []
    for region in region_list:
        if region in reg_indices:
            sum_indices += reg_indices[region]
        else:
            for sp in spines:
                if region + sp in reg_indices:
                    sum_indices += reg_indices[region + sp]
    return sum_indices



def get_concentrations_region_list(my_file, my_list, trial, out, specie,
                                   spines=["_sa1[0]"]):

    grid_list = get_grid_list(my_file)
    species = get_all_species(my_file)
    idx = species.index(specie)
    idxs = sum_indices(my_file, my_list, spines)
    data = get_populations(my_file, trial=trial, output=out)
    numbers = data[:, idxs, idx].sum(axis=1)
    return numbers


def get_fluo_sig(signal, t_init, dt):
    min_len = min([len(dat) for dat in signal])
    out = np.array([data[:min_len] for data in signal]).mean(axis=0)
    basal = np.mean(out[int(1000/dt):int(t_init/dt)])
    out = (out-basal)/basal
    return out, min_len


if __name__ == "__main__":
    fnames = []
    args = Parser().parse_args()
    for name in args.input:
        fnames.append(name)
    if not fnames:
        sys.exit('Do specify at least one totals filename')
    t_init = args.t_init
    specie = args.specie
    for fname in fnames:
    
        data_spine = []
        data_dend = []
        my_file = h5py.File(fname, "r")

        for key in my_file.keys():
            if not key.startswith("trial"):
                continue
            print(key)
            spine = get_concentrations_region_list(my_file,["PSD",
                                                            "head", "neck"],
                                                   key, "__main__",
                                                   specie)
        
            dend = get_concentrations_region_list(my_file,["dend06"],
                                                  key, "__main__",
                                                  specie)

            data_dend.append(dend)
            data_spine.append(spine)
     
        

        time = get_times(my_file, key, "__main__")  - t_init

        dt = time[1] - time[0]
        out_spine, min_len = get_fluo_sig(data_spine, t_init, dt)
        out_dend, min_len = get_fluo_sig(data_dend, t_init, dt)
        fig, ax = plt.subplots(1, 1, figsize=(5, 5))
      
        ax.plot(time[:min_len]/1000, out_spine*100, "tab:green", label="Spine")
        ax.plot(time[:min_len]/1000, out_dend*100, "tab:blue",
                label="Dendrite")
        ax.set_xlabel("time (s)", fontsize=15)
        ax.set_ylabel("Relative fluorescence change", fontsize=15)
        ax.legend()
        output = "%s_%s_dR_R.png" % (fname.split(".h5")[0],
                                     specie)
        ax.tick_params(axis='x', labelsize=15)
        ax.tick_params(axis='y', labelsize=15)
        ax.set_title(fname)
        fig.savefig(output, format="png",
                    dpi=100, bbox_inches="tight")
    plt.show()
