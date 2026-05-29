import sys
import argparse

import h5py
import numpy as np
import matplotlib.pyplot as plt 
from lxml import etree

from scipy.constants import Avogadro
import utility_functions as utils


NA = Avogadro*1e-23
plt.rcParams['text.usetex'] = True
specie_dict = {
    "Ca": ["Ca"],
    "CaOut": ["CaOut"],
    "CaER": ["CaER"],
    "RyR2O": ["RyR2CaMO1", "RyR2CaMO2"],
    "RyR3O": ["RyR3CaMO1", "RyR3CaMO2"],
    "STIM_CaER": ["STIM_2CaER"],
    "Orai": ["OraiSTIM_4", "Orai2STIM_4", "Orai3STIM_4"],
    "Fura": ["Fura2Ca"],
 
}
multiplier = {
    "Ca": 1,
    "CaOut": 1,
    "CaER": 1,
    "RyRO1": 1,
    "RyRO2": 1,
    "STIM_2CaER": 1,
    "OraiSTIM_4": 1,
    "Orai2STIM_4": 2,
    "Orai3STIM_4": 3,
    "Fura2Ca": 1,
}
def Parser():
    parser = argparse.ArgumentParser(description='Generate figs of avg conc')
    parser.add_argument('input', nargs='+',
                        help='input h5 files')
    parser.add_argument('--species', default="Ca",
                        help='Ca, RyRO, CaER, CaOut, RyRO, Fura')
    parser.add_argument('--output', default="Ca",
                        help='Ca, __main__')

    return parser

    



if __name__ == '__main__':

    args = Parser().parse_args()
    specie = args.species
    output = args.output
    base = "dend"
    reg_list = [base, "dend01", "dend02", "dend03", "dend04", "dend05",
                "dend06", "dend07", "dend08", "dend09",]
    for i in range(10, 102, 1):
        reg_list.append("%s%d" %(base, i))
   
    figs, axes = [], []
    if len(sys.argv) == 1:
        sys.exit('No filename given')
    for fname in sys.argv[1:]:
        print(fname)
        my_file = h5py.File(fname, 'r')
        conc_dict = {}
        time_dict = {}
        vmin = 1200
        vmax = 0
        peak = []
        for trial in my_file.keys():
            if trial == "model":
                continue
            try:
                conc, voxels = utils.get_dynamics_in_region(my_file,
                                                            [specie],
                                                            reg_list, trial,
                                                            output)
                time = utils.get_times(my_file, trial, output)
            except KeyError:
                conc, voxels = utils.get_dynamics_in_region(my_file,
                                                            [specie],
                                                            reg_list, trial,
                                                            "__main__")
                time = utils.get_times(my_file, trial, "__main__")
                
            conc_dict[trial] = conc
            time_dict[trial] = time
            new_max = conc_dict[trial].max()
            new_min = conc_dict[trial].min()
            peak.append(vmax)
            if len(np.where(conc[3000:]>202.5)[0])>1:
                print(len(np.where(conc[3000:]>202.5)[1]))
                #print(set(np.where(conc[3000:]>202.5)[0][1]))
            if new_max > vmax:
                vmax = new_max
            if new_min < vmin:
                vmin = new_min

        print(np.mean(peak), np.std(peak))
        diam = fname.split("diam_")[-1][:3]

        for key in conc_dict:
            fig, ax = plt.subplots(1, 1)
            time = time_dict[key]
            im = ax.imshow(conc_dict[key].T, aspect="auto",
                           interpolation="none",
                           origin="lower", extent = [time[0]*1e-3,
                                                     time[-1]*1e-3,
                                                     voxels[0],
                                                     voxels[-1]],
                           cmap=plt.get_cmap("Reds"))
            ax.set_xlabel(r"time (s)", fontsize=14)
            ax.set_ylabel(r"dendrite $(\mathrm{\mu m})$", fontsize=14)
            fig.colorbar(im)
            
            ax.set_title(r"%s dynamics in %s $\mathrm{\mu m}$ dend" % (specie, diam),
                         fontsize=14)
            fig.savefig(fname[:-3]+"_"+key+".png", dpi=100,
                        bbox_inches="tight")
            plt.close()
    
    
 
                          
