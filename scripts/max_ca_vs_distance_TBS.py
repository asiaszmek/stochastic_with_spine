import os

import matplotlib.pyplot as plt

import utility_functions as utils


colors = {"1.2": 'tab:blue',
          "2.4": 'tab:purple',
          "6.0": 'tab:green'}

files_no_ER_ctrl = os.path.join("..", "TBS",
                           "model_RyR2CaM_%s_um_dend_spine_TBS_in_vitro.h5")
files_ER_ctrl = os.path.join("..", "TBS",
                             "model_RyR2CaM_RyR3CaM_%s_um_dend_spine_TBS_in_vitro.h5")
files_ER_big_spine = os.path.join("..", "TBS",
                                  "model_RyR2CaM_RyR3CaM_%s_um_dend_big_spine_TBS_in_vitro.h5")
dend_diam = ["1.2", "2.4", "6.0"]
types = ["no spine ER", "spine ER", "big spine ER"]
markers = ["^", "o", "s"]


if __name__ == '__main__':
    fig =  utils.make_distance_fig([files_no_ER_ctrl, files_ER_ctrl, files_ER_big_spine],
                                   t_init=3000, dend_diam=dend_diam, what_species="Ca",
                                   output_name="all",
                                   colors=colors, types=types, markers=markers)
   
    fig.savefig("Max_Ca_conc_TBS_in_vitro.png", dpi=100, bbox_inches="tight")
   
    plt.show()
