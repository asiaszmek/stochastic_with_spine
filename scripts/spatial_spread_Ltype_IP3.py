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
file_dend_ER_IP3 = os.path.join("..", "1_spine_different_Ca_amp",
                                "model_RyR2CaM_IP3_1.2_um_dend_spine_%s_%d_stim.h5")
file_dend_ER_Ltype_IP3 = os.path.join("..", "1_spine_different_Ca_amp_Ca_dend",
                                      "model_Ca_dend_RyR2CaM_IP3_1.2_um_dend_spine_%s_%d_stim.h5")

file_dend_ER_Ltype = os.path.join("..", "1_spine_different_Ca_amp_Ca_dend",
                                  "model_Ca_dend_RyR2CaM_1.2_um_dend_spine_%s_%d_stim.h5")

file_spine_ER_dend_ER = os.path.join("..", "1_spine_different_Ca_amp",
                                     "model_RyR2CaM_RyR3CaM_1.2_um_dend_spine_%s_%d_stim.h5")
file_spine_ER_dend_ER_IP3 = os.path.join("..", "1_spine_different_Ca_amp",
                                         "model_RyR2CaM_RyR3CaM_IP3_1.2_um_dend_spine_%s_%d_stim.h5")
file_spine_ER_dend_ER_Ltype_IP3 = os.path.join("..",
                                               "1_spine_different_Ca_amp_Ca_dend",
                                               "model_Ca_dend_RyR2CaM_RyR3CaM_IP3_1.2_um_dend_spine_%s_%d_stim.h5")

file_spine_ER_dend_ER_Ltype = os.path.join("..",
                                           "1_spine_different_Ca_amp_Ca_dend",
                                           "model_Ca_dend_RyR2CaM_RyR3CaM_1.2_um_dend_spine_%s_%d_stim.h5")

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
types = [r"no SA no VGCC no $\mathrm{IP_3R}$",r"no SA + VGCC no $\mathrm{IP_3R}$", r"no SA + VGCC + $\mathrm{IP_3R}$", r"SA no VGCC no $\mathrm{IP_3R}$",r"SA + VGCC no $\mathrm{IP_3R}$", r"SA + VGCC + $\mathrm{IP_3R}$",]
markers = ["o", "^", "s",  "o", "^", "s"]
fillstyles = ["none", "none", "none", "full", "full", "full"]
stim_dict = {
    "1EPSP": [1, 2, 3, 5, 8],
    "4EPSP_100Hz": [1, 2, 3, 5, 8],
}

specie="Ca"
if __name__ == '__main__':
    colors = {file_dend_ER_ctrl: "tab:blue",
              file_dend_ER_Ltype: "tab:green",
              file_dend_ER_Ltype_IP3: "tab:purple",
              file_spine_ER_dend_ER: "tab:blue",
              file_spine_ER_dend_ER_Ltype: "tab:green",
              file_spine_ER_dend_ER_Ltype_IP3: "tab:purple"
    }
    dirs = [file_dend_ER_ctrl,file_dend_ER_Ltype, file_dend_ER_Ltype_IP3,
            file_spine_ER_dend_ER, file_spine_ER_dend_ER_Ltype,
            file_spine_ER_dend_ER_Ltype_IP3]
    fig1, fig2, fig3 = utils.make_distance_figs(dirs, stim_dict, output,
                                                types, markers, fillstyles,
                                                t_init=t_init, max_ca=True,
                                                colors=colors, specie=specie)
    fig1.savefig("Spatial_spread_Ltype_IP3_%s.png" % specie, dpi=100,
                 bbox_inches="tight")
    fig2.savefig("Head_%s_vs_dend_%s_spread_Ltype_IP3R.png" % (specie, specie),
                 dpi=100, bbox_inches="tight")
    fig3.savefig("Dend_%s_spatial_spread_Ltype_IP3.png" % specie, dpi=100,
                 bbox_inches="tight")
    
    #plt.show()
