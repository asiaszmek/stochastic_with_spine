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
file_dend_ER_IP3 = os.path.join("..", "1_spine_different_Ca_amp",
                                "model_RyR2CaM_IP3_1.2_um_dend_spine_%s_ms_%s.h5")

file_spine_ER_dend_ER = os.path.join("..", "1_spine_different_Ca_amp",
                                      "model_RyR2CaM_RyR3CaM_1.2_um_dend_spine_%s_ms_%s.h5")
file_spine_ER_dend_ER_IP3 = os.path.join("..", "1_spine_different_Ca_amp",
                                         "model_RyR2CaM_RyR3CaM_IP3_1.2_um_dend_spine_%s_ms_%s.h5")

file_spine_ER_no_dend_ER = os.path.join("..", "1_spine_different_Ca_amp",
                                         "model_RyR3CaM_1.2_um_dend_spine_%s_ms_%s.h5")
file_dend_ER_dend_Ca = os.path.join("..", "1_spine_different_Ca_amp_Ca_dend",
                                  "model_Ca_dend_RyR2CaM_1.2_um_dend_spine_%s_ms_%s.h5")
file_spine_ER_dend_ER_dend_Ca = os.path.join("..", "1_spine_different_Ca_amp_Ca_dend",
                                      "model_Ca_dend_RyR2CaM_RyR3CaM_1.2_um_dend_spine_%s_ms_%s.h5")
t_init = 3000
output = "Ca"
types = ["no SA no IP3 production", "no SA", "SA no IP3 production", "SA",]
markers = ["s", "s", "o", "o"]
fillstyles = ["none", "full", "none", "full"]
stim_dict = {
    "4": ["01750", "03500", "07000", "10500"],
    "40": ["0175", "0350", "0700", "1050"],
    "400": ["017.5", "035", "070", "105"],
}
spine_dend = "dend11"
region_list = ["dend01", "dend02", "dend03", "dend04", "dend05", "dend06", "dend07","dend08", "dend08", "dend09", "dend10",
               "dend11", "dend12", "dend13", "dend14", "dend15", "dend16", "dend17",  "dend18", "dend18", "dend19", "dend20",
               "dend21"]
spine_idx = 21
if __name__ == '__main__':
    for dur in stim_dict.keys():
        
        fig1, ax1 = plt.subplots(1, 1)
        
        
        period_spine = 1000+int(dur)
        period_dend = 500+int(dur)
        
        for i, fname in enumerate([file_dend_ER_ctrl,file_dend_ER_IP3, file_spine_ER_dend_ER,
                                   file_spine_ER_dend_ER_IP3]):
            auc_head_m = []
            auc_head_e = []
            spread_m = []
            spread_e = []
           

            for inj in stim_dict[dur]:
              
                print(fname % (dur, inj))
                if inj in ["400", "800"]:
                    continue

                try:
                    my_data = h5py.File(fname % (dur, inj))
                except FileNotFoundError:
                    print("Cound not open ", fname % (dur, inj))
                    continue
        
                dend_ca, time_dict = utils.get_conc(my_data, "Ca", region_list, output)
                head_ca, time_head_dict = utils.get_conc(my_data, "Ca", ['head'], output)
                
               
                dt = time_dict["trial0"][1]-time_dict["trial0"][0]
                t_start = int(t_init/dt)
               
                length = time_dict["trial0"][-1] - t_start
                nothing_spine = period_spine*dt*71
                nothing_dend = period_dend*dt*71
                trials = len(dend_ca.keys())
                auch = np.zeros((trials))
                spread = utils.get_distance(dend_ca, dt, t_init=3000, stim_len=int(dur), length=42)
                for j, trial in enumerate(dend_ca.keys()):
                    auch[j] = head_ca[trial][1, t_start:t_start+period_spine].sum()/nothing_spine
                    
                print(spread)
                auc_head_m.append(auch.mean())
                auc_head_e.append(auch.std()/trials**0.5)
                spread_m.append(spread.mean())
                spread_e.append(spread.std()/trials**0.5)
              
            print(spread_m)
            ax1.errorbar(x=auc_head_m, y=spread_m, xerr=auc_head_e, yerr=spread_e,label=types[i], marker=markers[i],
                         fillstyle=fillstyles[i], color="tab:blue", linewidth=0)

           
            
        ax1.set_ylabel("Spatial spread (um)")
        ax1.set_xlabel("Auc head Ca/Auc basal")

      

       

        ax1.legend()
        fig1.savefig("Spatial_spread_IP3_%s_ms.png" % dur, dpi=100, bbox_inches="tight")
        
    
    #plt.show()
