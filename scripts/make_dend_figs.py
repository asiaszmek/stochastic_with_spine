import sys
import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt 
from lxml import etree
import utility_functions as utils
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
    "Orai2STIM_4": 2,
    "Orai3STIM_4": 3,
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
    if chosen_specie in ["Ca", "CaER", "CaOut", "RyRO", "Fura"]:
        output_name = "all"
    elif  chosen_specie in ["STIM_CaER", "Orai"]:
        output_name = "RyR_Orai"
    else:
        output_name = "__main__"

    try:
        specie_list = specie_dict[chosen_specie]
    except KeyError:
        specie_list = [chosen_specie]
    base = "dend"
    reg_list = [base, "dend01", "dend02", "dend03", "dend04", "dend05",
                "dend06", "dend07", "dend08", "dend09",]
    for i in range(10, 102, 1):
        reg_list.append("%s%d" %(base, i))
    fig1, ax1 = plt.subplots(1, len(fnames))
    im_list = []
    if len(fnames) == 1:
        ax1 = [ax1]
    for j, fname in enumerate(fnames):
        my_file = h5py.File(fname, 'r')
        conc_dict = {}
        time_dict = {}
        
        for trial in my_file.keys():
            print(trial)
            if trial == "model":
                continue
            try:
                conc, voxels = utils.get_dynamics_in_region(my_file,
                                                            specie_list,
                                                            reg_list, trial,
                                                            output_name)
                time = utils.get_times(my_file, trial, output_name)
            except KeyError:
                conc, voxels = utils.get_dynamics_in_region(my_file,
                                                            specie_list,
                                                            reg_list, trial,
                                                            "__main__")
                time = utils.get_times(my_file, trial, "__main__")
            conc_dict[trial] = conc
        
        
        
            time_dict[trial] = time
        vmax = max([conc.max() for conc in conc_dict.values()])
        vmin = min([conc.min() for conc in conc_dict.values()])
       
        lmin = min([len(conc) for conc in conc_dict.values()])
        
        shape2 = max([conc.shape[1] for conc in conc_dict.values()])
        conc_mean = np.zeros((lmin, shape2))
        for conc in conc_dict.values():
            conc_mean[:lmin, :] += conc[:lmin, :]
        conc_mean /= len(conc_dict)
        if chosen_specie == "Fura":
            mean = conc_mean[:,:int(3000/(time[1]-time[0]))].mean(axis=0)
            conc_mean = (conc_mean - mean)/mean

        # for i, key in enumerate(conc_dict):
        #     fig, ax = plt.subplots(1, 1)
        #     time = time_dict[key]
        #     im = ax.imshow(conc_dict[key].T, aspect="auto",
        #                    interpolation="none",
        #                    origin="lower", extent = [time[0]*1e-3,
        #                                              time[-1]*1e-3,
        #                                              voxels[0],
        #                                              voxels[-1]],
        #                    cmap=plt.get_cmap("Reds"), vmin=vmin, vmax=vmax)
        #     fig.colorbar(im)
        #     ax.set_title("%s trial %d %s" % (fname, i, specie))
        
        if args.scale == "log":
            im_list.append(np.log10(1e-9*conc_mean.T))
        elif args.scale == "linear":
            im_list.append(conc_mean.T)
    vmax = max([im.max() for im in im_list])
    vmin = min([im.min() for im in im_list])
    for j, x in enumerate(ax1):
        im = ax1[j].imshow(im_list[j], aspect="auto",
                          interpolation="none",
                          origin="lower", extent = [time[0]*1e-3,
                                                    time[-1]*1e-3,
                                                    voxels[0],
                                                    voxels[-1]],
                          cmap=plt.get_cmap("Reds"), vmin=vmin,
                          vmax=vmax)

        if "baloon" in fnames[j]:
            ax1[j].set_title("baloon ER mean %s" %  chosen_specie)
        elif "nc_tubes" in fnames[j]:
            ax1[j].set_title("RyR overexpression mean %s" %  chosen_specie)
        elif "tubes" in fnames[j]:
            ax1[j].set_title("tubes ER mean %s" %  chosen_specie)
            
    fig1.colorbar(im)
    plt.show()
                          
