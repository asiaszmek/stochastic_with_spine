import sys
import h5py
import numpy as np
import utility_functions as utils


length = 25000 #  sec
output = "__main__"

if __name__ == "__main__":
    for fname in sys.argv[1:]:
        print(fname)
        my_file = h5py.File(fname, 'r')
        conc_dict = {}
        time_dict = {}
        
        species = utils.get_all_species(my_file, output=output)
        regions = utils.get_regions(my_file)
        for trial in my_file.keys():
            if trial == "model":
                continue
            conc_dict[trial] = utils.get_concentrations(my_file, trial,
                                                        "__main__")
        times = utils.get_times(my_file, trial=trial, output=output)
        dt = times[1] - times[0]
        start = -int(length/dt)

        min_len = min([conc.shape[0] for conc in conc_dict.values()])
        temp_conc = np.zeros((min_len, len(regions), len(species)))
        for conc in conc_dict.values():
            temp_conc += conc[:min_len]
        mean_conc = temp_conc[start:].mean(axis=(0, 1))/len(conc_dict)
        for i, specie in enumerate(species):
            print(specie, mean_conc[i])
        
        
        
