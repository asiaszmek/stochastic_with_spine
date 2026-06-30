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
                                "model_1.2_um_dend_spine_%s_ms_%s.h5")
file_dend_ER_ctrl = os.path.join("..", "1_spine_different_Ca_amp",
                                  "model_RyR2CaM_1.2_um_dend_spine_%s_ms_%s.h5")
file_spine_ER_dend_ER = os.path.join("..", "1_spine_different_Ca_amp",
                                      "model_RyR2CaM_RyR3CaM_1.2_um_dend_spine_%s_ms_%s.h5")
file_spine_ER_no_dend_ER = os.path.join("..", "1_spine_different_Ca_amp",
                                         "model_RyR3CaM_1.2_um_dend_spine_%s_ms_%s.h5")
file_dend_ER_dend_Ca = os.path.join("..", "1_spine_different_Ca_amp_Ca_dend",
                                  "model_Ca_dend_RyR2CaM_1.2_um_dend_spine_%s_ms_%s.h5")
file_dend_ER_dend_Ca_old_age = os.path.join("..", "1_spine_different_Ca_amp_Ca_dend",
                                  "model_Ca_dend_RyR2CaM_old_age_1.2_um_dend_spine_%s_ms_%s.h5")
file_spine_ER_dend_ER_dend_Ca_old_age = os.path.join("..",
                                                     "1_spine_different_Ca_amp_Ca_dend",
                                                     "model_Ca_dend_RyR2CaM_RyR3CaM_old_age_1.2_um_dend_spine_%s_ms_%s.h5")
file_spine_ER_dend_ER_dend_Ca = os.path.join("..",
                                             "1_spine_different_Ca_amp_Ca_dend",
                                             "model_Ca_dend_RyR2CaM_RyR3CaM_1.2_um_dend_spine_%s_ms_%s.h5")

t_init = 3000
output = "Ca"
types = ["no SA+L-type", "no SA+L-type old age", "SA+L-type", "SA+L-type old age",]
markers = ["s", "s", "o", "o"]
fillstyles = ["full", "none", "full", "none"]
stim_dict = {
    "4": ["01750", "03500", "07000", "10500"],
    "40": ["0175", "0350", "0700", "1050"],
    "400": ["017.5", "035", "070", "105", "210"],
}
spine_dend = "dend11"
region_list = ["dend01", "dend02", "dend03", "dend04", "dend05", "dend06", "dend07","dend08", "dend08", "dend09", "dend10",
               "dend11", "dend12", "dend13", "dend14", "dend15", "dend16", "dend17",  "dend18", "dend18", "dend19", "dend20",
               "dend21"]
spine_idx = 21
if __name__ == '__main__':
    rollo = False
    dirs = [
        file_dend_ER_dend_Ca, file_dend_ER_dend_Ca_old_age,
        file_spine_ER_dend_ER_dend_Ca,
        file_spine_ER_dend_ER_dend_Ca_old_age,
    ]
    fig = utils.make_distance_figs(dirs, stim_dict, output, types, markers,
                                   fillstyles, length=42, spine_idx=21, t_init=t_init,
                                   max_ca=True)
    fig.savefig("Spatial_spread_Ltype_old_age_%s.png" % rollo, dpi=100,
                bbox_inches="tight")
