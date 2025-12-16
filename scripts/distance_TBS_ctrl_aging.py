import os

import matplotlib.pyplot as plt

import utility_functions as utils


colors = {"1.2": 'tab:blue',
          "2.4": 'tab:purple',
          "6.0": 'tab:green'}

files = [
         os.path.join("..", "TBS",
                      "model_RyR2CaM_%s_um_dend_spine_TBS_in_vivo.h5"),
         os.path.join("..", "TBS",
                      "model_RyR2CaM_RyR3CaM_%s_um_dend_spine_TBS_in_vivo.h5"),
         os.path.join("..", "TBS",
                      "model_RyR2CaM_%s_um_dend_spine_old_age_TBS_in_vivo.h5"),
         os.path.join("..", "TBS",
                      "model_RyR2CaM_RyR3CaM_%s_um_dend_spine_old_age_TBS_in_vivo.h5")]
    

dend_diam = ["1.2", "2.4"]
types = ["no spine ER in vivo", "RyR3CaM in neck in vivo", "old age no spine ER in vivo", "old age RyR3CaM in neck in vivo"]
markers = ["^", "o", "s", "d"]


if __name__ == '__main__':
    fig =  utils.make_distance_fig_compare_with_mean(files,
                                                     t_init=6000, stim_len=0, dend_diam=dend_diam,
                                                     what_species="Ca",
                                                     output_name="all",
                                                     colors=colors, types=types,
                                                     markers=markers)
   
    fig.savefig("TBS_distance_aging.png", dpi=100, bbox_inches="tight")
   
    plt.show()
