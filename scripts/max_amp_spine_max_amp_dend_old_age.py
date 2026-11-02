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
file_dend_ER_oa = os.path.join("..", "1_spine_different_Ca_amp",
                                  "model_RyR2CaM_old_age_1.2_um_dend_spine_%s_%d_stim.h5")
file_spine_ER_dend_ER = os.path.join("..", "1_spine_different_Ca_amp",
                                      "model_RyR2CaM_RyR3CaM_1.2_um_dend_spine_%s_%d_stim.h5")
file_spine_ER_dend_ER_oa = os.path.join("..", "1_spine_different_Ca_amp",
                                      "model_RyR2CaM_RyR3CaM_old_age_1.2_um_dend_spine_%s_%d_stim.h5")

file_spine_ER_no_dend_ER = os.path.join("..", "1_spine_different_Ca_amp",
                                         "model_RyR3CaM_1.2_um_dend_spine_%s_%d_stim.h5")
t_init = 3000
output = "Ca"
types = ["ER no SA ctrl", "ER no SA old age", "ER+SA ctrl", "ER+SA old age"]
markers = ["s", "s", "o", "o"]
fillstyles = ["full", "none", "full", "none"]
stim_dict = {
    "1EPSP": [1, 5, 10, 15, 20],
    "4EPSP_100Hz": [1, 5, 10, 15, 20],
}

if __name__ == '__main__':
    directories = [file_dend_ER_ctrl, file_dend_ER_oa, file_spine_ER_dend_ER, file_spine_ER_dend_ER_oa] 
    fig1, fig2, fig3 = utils.max_vs_auc_head_neck_dend(directories,
                                                       stim_dict, output,
                                                       types, markers,
                                                       fillstyles,
                                                       length=102,
                                                       which_dend=["dend25"],
                                                       t_init=3000)
       

    fig1.savefig("Old_age_no_Ltype_aucdend_auchead.png", dpi=100, bbox_inches="tight")
    fig2.savefig("Old_age_no_Ltype_maxdend_maxhead.png", dpi=100, bbox_inches="tight")
    fig3.savefig("Old_age_no_Ltype_aucneck_auchead.png", dpi=100, bbox_inches="tight")
