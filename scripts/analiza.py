#!/usr/bin/env python
import h5py
import numpy as np
from lxml import etree
import sys
from scipy.constants import Avogadro
import utility_functions as utils


if __name__ == '__main__':
    if len(sys.argv) == 1:
        sys.exit('No filename given')
    for fname in sys.argv[1:]:
        my_file = h5py.File(fname, 'r')
        for trial in my_file.keys():
            
            if trial.startswith("trial"):
                try:
                    utils.save_concentrations(my_file, fname[:-3], 'all',
                                              trial=trial)
                except KeyError:
                    pass
                try:
                    utils.save_concentrations(my_file, fname[:-3], 'lowres',
                                              trial=trial)
                except KeyError:
                    pass
                try:
                    utils.save_concentrations(my_file, fname[:-3], 'fluo_dend',
                                              trial=trial)
                except KeyError:
                    pass
                try:
                    utils.save_concentrations(my_file, fname[:-3], 'fluo_PSD',
                                              trial=trial)
                except KeyError:
                    pass
                try:
                    utils.save_concentrations(my_file, fname[:-3], 'fluo_head',
                                              trial=trial)
                except KeyError:
                    pass
                try:
                    utils.save_concentrations(my_file, fname[:-3], 'fluo_neck',
                                              trial=trial)
                except KeyError:
                    pass

                utils.save_concentrations(my_file, fname[:-3], '__main__',
                                          trial=trial)
                

    print('Done')

