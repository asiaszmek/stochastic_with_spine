import os
import numpy as np
import matplotlib.pyplot as plt
import h5py
import utility_functions as utils

spines = {'sa1[0]': ['PSD_sa1[0]', 'head_sa1[0]', 'neck_sa1[0]']}
colors = {"1.2": 'tab:blue',
          "2.4": 'tab:purple',
          "6.0": 'tab:green'}

file_no_ER_ctrl = os.path.join("..", "1_spine_different_Ca_amp",
                               "model_1.2_um_dend_spine_%s_%d_stim.h5")
file_dend_ER_ctrl = os.path.join("..", "1_spine_different_Ca_amp",
                                 "model_RyR2CaM_1.2_um_dend_spine_%s_%d_stim.h5")
file_spine_ER_dend_ER = os.path.join("..", "1_spine_different_Ca_amp",
                                     "model_RyR2CaM_RyR3CaM_1.2_um_dend_spine_%s_%d_stim.h5")
file_spine_ER_no_dend_ER = os.path.join("..", "1_spine_different_Ca_amp",
                                        "model_RyR3CaM_1.2_um_dend_spine_%s_%d_stim.h5")
file_dend_ER_dend_Ca = os.path.join("..",
                                    "1_spine_different_Ca_amp_Ca_dend",
                                    "model_Ca_dend_RyR2CaM_1.2_um_dend_spine_%s_%d_stim.h5")
file_spine_ER_dend_ER_dend_Ca = os.path.join("..",
                                             "1_spine_different_Ca_amp_Ca_dend",
                                             "model_Ca_dend_RyR2CaM_RyR3CaM_1.2_um_dend_spine_%s_%d_stim.h5")
t_init = 3000
output = "Ca"

types = ["no ER", "ER no SA", "ER+SA", "SA"]
markers = ["o", "s", "o", "^"]
fillstyles = ["none", "none", "full", "none"]
stim_dict = {
    "1EPSP":  [1, 2, 3, 5, 8],
    "4EPSP_100Hz":  [1, 2, 3, 5, 8],
    # "100EPSP_100Hz": [1, 5, 10, 15, 20],
}

specie = "Ca"
if __name__ == '__main__':
    dirs = [file_no_ER_ctrl, file_dend_ER_ctrl, file_spine_ER_dend_ER, file_spine_ER_no_dend_ER]
    fig1, fig2 , fig3= utils.make_distance_figs(dirs, stim_dict, output, types, markers,
                                          fillstyles, length=42,
                                          t_init=t_init, max_ca=True, specie=specie)
    fig1.savefig("Spatial_spread_ER_no_ER_%s.png" % specie, dpi=100, bbox_inches="tight")
    fig2.savefig("Head_%s_vs_dend_%s_ER_no_ER.png" % (specie,
                                                      specie), dpi=100, bbox_inches="tight")
    fig3.savefig("Dend_%s_spatial_spread_ER_no_ER.png" % specie, dpi=100, bbox_inches="tight")
