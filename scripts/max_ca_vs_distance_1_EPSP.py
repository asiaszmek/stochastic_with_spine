import os

import matplotlib.pyplot as plt

import utility_functions as utils


colors = {"1.2": 'tab:blue',
          "2.4": 'tab:purple',
          "6.0": 'tab:green'}

files_no_ER_ctrl = os.path.join("..", "1_spine_1_EPSP",
                           "model_RyR2CaM_%s_um_dend_spine_1_EPSP.h5")
files_ER_ctrl = os.path.join("..", "1_spine_1_EPSP",
                             "model_RyR2CaM_RyR3CaM_%s_um_dend_spine_1_EPSP.h5")
files_ER_big_spine = os.path.join("..", "1_spine_1_EPSP",
                                  "model_RyR2CaM_RyR3CaM_%s_um_dend_big_spine_1_EPSP.h5")
files_ER_big_spine_IC = os.path.join("..", "1_spine_1_EPSP",
                                  "model_RyR2CaM_RyR3CaM_%s_um_dend_big_spine_IC_1_EPSP.h5")

dend_diam = ["1.2", "2.4", "6.0"]
types = ["no spine ER", "spine ER", "big spine ER", "big spine RyR3CaM in spine head"]
markers = ["^", "o", "s", "*"]


if __name__ == '__main__':
    fig =  utils.make_distance_fig([files_no_ER_ctrl, files_ER_ctrl, files_ER_big_spine,
                                    files_ER_big_spine_IC],
                                   t_init=3000, dend_diam=dend_diam, what_species="Ca",
                                   output_name="all",
                                   colors=colors, types=types, markers=markers)
    plt.show()
    fig.savefig("Max_Ca_conc_1_EPSP_ctrl.png", dpi=100, bbox_inches="tight")
   
