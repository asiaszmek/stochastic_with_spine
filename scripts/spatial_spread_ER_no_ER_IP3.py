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
    "1EPSP": [1, 5, 10, 15, 20],
    "4EPSP_100Hz": [1, 5, 10, 15, 20],
    "100EPSP_100Hz": [1, 5, 10, 15, 20],
}
spine_dend = "dend11"
region_list = ["dend01", "dend02", "dend03", "dend04", "dend05", "dend06", "dend07","dend08", "dend08", "dend09", "dend10",
               "dend11", "dend12", "dend13", "dend14", "dend15", "dend16", "dend17",  "dend18", "dend18", "dend19", "dend20",
               "dend21"]

specie = "Ca"
if __name__ == '__main__':
    rollo = False
    dirs = [file_dend_ER_ctrl, file_dend_ER_ctrl_IP3, file_spine_ER_dend_ER, file_spine_ER_dend_ER_IP3]
    fig = utils.make_distance_figs(dirs, stim_dict, output, types, markers,
                                   fillstyles, length=42,
                                   t_init=t_init, max_ca=True, specie=specie)
    fig.savefig("Spatial_spread_ER_no_ER_IP3R_%s.png" % specie, dpi=100, bbox_inches="tight")
