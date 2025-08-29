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
    parser.add_argument('--scale', default="linear",
                        help='linear, log')

    return parser

    



if __name__ == '__main__':
    specie_list = ["Ca"]
    specie = "Ca"
    base = "dend"
    reg_list = [base, "dend01", "dend02", "dend03", "dend04", "dend05",
                "dend06", "dend07", "dend08", "dend09",]
    for i in range(10, 102, 1):
        reg_list.append("%s%d" %(base, i))
   
    figs, axes = [], []
    if len(sys.argv) == 1:
        sys.exit('No filename given')
    for fname in sys.argv[1:]:
        my_file = h5py.File(fname, 'r')
        conc_dict = {}
        time_dict = {}
        for trial in my_file.keys():
            if trial == "model":
                continue
            conc, voxels = utils.get_dynamics_in_region(my_file,
                                                        specie_list,
                                                        reg_list, trial, "__main__")
            conc_dict[trial] = conc
            time = utils.get_times(my_file, trial, "__main__")
            time_dict[trial] = time
        vmin = 0
        vmax = 1200
        diam = fname.split("diam_")[-1][:3]
        # for key in conc_dict:
        #     new_max = conc_dict[key].max()
        #     if new_max > vmax:
        #         vmax = new_max
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
    
    
    
 
                          
