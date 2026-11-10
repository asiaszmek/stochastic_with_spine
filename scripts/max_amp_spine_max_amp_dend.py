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
t_init = 3000
output = "Ca"
types = ["no ER", "ER no SA", "ER+SA", "SA"]
markers = ["o", "s", "o", "^"]
fillstyles = ["none", "none", "full", "none"]
stim_dict = {
    "4": ["01750", "03500", "07000", "10500"],
    "40": ["0175", "0350", "0700", "1050"],
    "400": ["017.5", "035", "070", "105"],
}
spine_dend = "dend11"
region_list = ["dend01", "dend02", "dend03", "dend04", "dend05", "dend06", "dend07","dend08", "dend08", "dend09", "dend10",
               "dend11", "dend12", "dend13", "dend14", "dend15", "dend16", "dend17",  "dend18", "dend18", "dend19", "dend20",
               "dend21"]

if __name__ == '__main__':
    directories = [file_no_ER_ctrl, file_dend_ER_ctrl, file_spine_ER_dend_ER, file_spine_ER_no_dend_ER]
    fig1, fig2, fig3 = utils.max_vs_auc_head_neck_dend(directories,
                                                       stim_dict, output,
                                                       types, markers,
                                                       fillstyles,
                                                       length=42,
                                                       which_dend=["dend11"],
                                                       t_init=3000)
       

    fig1.savefig("ER_no_ER_aucdend_auchead.png", dpi=100, bbox_inches="tight")
    fig2.savefig("ER_no_ER_maxdend_maxhead.png", dpi=100, bbox_inches="tight")
    fig3.savefig("ER_no_ER_aucneck_auchead.png", dpi=100, bbox_inches="tight")
