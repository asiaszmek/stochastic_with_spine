import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import utility_functions as utils

spines = {'sa1[0]': ['PSD_sa1[0]', 'head_sa1[0]', 'neck_sa1[0]']}
colors = {"1.2": 'tab:blue',
          "2.4": 'tab:purple',
          "6.0": 'tab:green'}


file_dend_ER_ctrl = os.path.join("..", "1_spine_different_Ca_amp",
                                  "model_RyR2CaM_1.2_um_dend_spine_%s_%d_stim.h5")
file_dend_ER_ctrl_IP3 = os.path.join("..", "1_spine_different_Ca_amp",
                                     "model_RyR2CaM_IP3_1.2_um_dend_spine_%s_%d_stim.h5")

file_spine_ER_dend_ER = os.path.join("..", "1_spine_different_Ca_amp",
                                     "model_RyR2CaM_RyR3CaM_1.2_um_dend_spine_%s_%d_stim.h5")
file_spine_ER_dend_ER_IP3 = os.path.join("..", "1_spine_different_Ca_amp",
                                     "model_RyR2CaM_RyR3CaM_IP3_1.2_um_dend_spine_%s_%d_stim.h5")

t_init = 3000
output = "Ca"
types = ["ER no SA", "ER no SA+IP3R", "ER+SA", "ER+SA+IP3R"]
markers = ["o", "o", "s", "s"]
fillstyles = ["none", "full", "none", "full"]
stim_dict = {
    "1EPSP": [1, 2, 3, 5, 8],
    "4EPSP_100Hz": [1, 2, 3, 5, 8],
}

specie = "Ca"

if __name__ == '__main__':
    print("ER_no_ER_IP3_aucdend_auchead_%s.png" %specie)
    directories = [ file_dend_ER_ctrl, file_dend_ER_ctrl_IP3, file_spine_ER_dend_ER, file_spine_ER_dend_ER_IP3]
    fig1, fig2, fig3 = utils.max_vs_auc_head_neck_dend(directories,
                                                       stim_dict, output,
                                                       types, markers,
                                                       fillstyles,
                                                       which_dend=["dend25"],
                                                       t_init=3000,
                                                       specie=specie)
       
    print("ER_no_ER_IP3R_aucdend_auchead_%s.png" %specie)
    fig1.savefig("ER_no_ER_IP3R_aucdend_auchead_%s.png" %specie,
                 dpi=100, bbox_inches="tight")
    fig2.savefig("ER_no_ER_IP3R_maxdend_maxhead_%s.png" %specie, dpi=100,
                 bbox_inches="tight")
    fig3.savefig("ER_no_ER_IP3R_aucneck_auchead_%s.png" %specie, dpi=100,
                 bbox_inches="tight")
