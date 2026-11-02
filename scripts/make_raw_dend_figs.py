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
    "RyR2O": [ "RyR2CaMO2", "RyR2CaMO1"],
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
    fnames = []
    args = Parser().parse_args()
    for name in args.input:
        if name.endswith("h5"):
            fnames.append(name)
    if not fnames:
        sys.exit('Do specify at least one totals filename')
    chosen_specie = args.species
    if chosen_specie in ["Ca", "CaER", "CaOut", "RyR2O", "RyR3O", "Fura"]:
        output_name = "all"
    elif  chosen_specie in ["STIM_CaER", "Orai"]:
        output_name = "__main__"
    else:
        output_name = "all"

    try:
        specie_list = specie_dict[chosen_specie]
    except KeyError:
        specie_list = [chosen_specie]
    print(specie_list, output_name)
    base = "dend"
    reg_list = [base, "dend01", "dend02", "dend03", "dend04", "dend05",
                "dend06", "dend07", "dend08", "dend09",]
    for i in range(10, 102, 1):
        reg_list.append("%s%d" %(base, i))
   
    for j, fname in enumerate(fnames):
        my_file = h5py.File(fname)
        for trial in my_file.keys():
            if trial == "model":
                continue
        
            conc, vols = utils.get_dynamics_in_region(my_file, specie_list,
                                                      reg_list,
                                                      trial, output_name)

            time = utils.get_times(my_file, trial, output_name)
            print(conc.shape)
            dt = time[1]-time[0]
            time = time[int(3000//dt):len(time)//2]
            conc = conc[int(3000//dt):len(time)//2]
            fig, ax = plt.subplots(1, 1)
            im = ax.imshow(conc.T, aspect="auto",interpolation="none",
                          
                           origin="lower", extent = [time[0]*1e-3,
                                                     time[-1]*1e-3,
                                                     vols[0],
                                                     vols[-1]],
                           cmap=plt.get_cmap("Reds"))
            ax.set_xlabel(r"time (s)", fontsize=14)
            ax.set_ylabel(r"dendrite $(\mathrm{\mu m})$", fontsize=14)
            fig.colorbar(im)
            
            ax.set_title(fname,
                         fontsize=14)
            fig.savefig(fname[:-3]+"_"+trial+'_'+chosen_specie+".png", dpi=100,
                        bbox_inches="tight")
            plt.show()
    
    
    
 
                          
