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
    for dur in stim_dict.keys():
        fig, ax = plt.subplots(1, 1)
        fig1, ax1 = plt.subplots(1, 1)
        fig_n, ax_n = plt.subplots(1, 1)
        fig1_n, ax1_n = plt.subplots(1, 1)
        fig2, ax2 = plt.subplots(1, 1)
        
        
        period_spine = 1000+int(dur)
        period_dend = 1000+int(dur)
        
        for i, fname in enumerate([file_no_ER_ctrl, file_dend_ER_ctrl, file_spine_ER_dend_ER, file_spine_ER_no_dend_ER]):
            auc_head_m = []
            auc_neck_m = []
            auc_dend_m = []
            max_head_m = []
            max_neck_m = []
            max_dend_m = []
            auc_head_e = []
            auc_neck_e = []
            auc_dend_e = []
            max_head_e = []
            max_neck_e = []
            max_dend_e = []
            
           

            for inj in stim_dict[dur]:
              
                print(fname % (dur, inj))
                try:
                    my_data = h5py.File(fname % (dur, inj))
                except FileNotFoundError:
                    print("Cound not open ", fname % (dur, inj))
                    continue
        
                dend_ca, time_dict = utils.get_conc(my_data, "Ca", ["dend11"], output)
                head_ca, time_head_dict = utils.get_conc(my_data, "Ca", ['head'], output)
                neck_ca, time_head_dict = utils.get_conc(my_data, "Ca", ['neck'], output)
               
                dt = time_dict["trial0"][1]-time_dict["trial0"][0]
                t_start = int(t_init/dt)
               
                length = time_dict["trial0"][-1] - t_start
                nothing_spine = period_spine*dt*71
                nothing_dend = period_dend*dt*71
                trials = len(dend_ca.keys())
                auch = np.zeros((trials))
                aucn = np.zeros((trials))
                aucd = np.zeros((trials))
                maxh = np.zeros((trials))
                maxn = np.zeros((trials))
                maxd = np.zeros((trials))
              

                for j, trial in enumerate(dend_ca.keys()):
                    maxh[j] = head_ca[trial][1, t_start:t_start+period_spine].max()/1000
                    maxn[j] = neck_ca[trial][:, t_start:t_start+period_spine].max()/1000
                    maxd[j] = dend_ca[trial][:, t_start:t_start+period_dend].mean(axis=0).max()
                    auch[j] = head_ca[trial][1, t_start:t_start+period_spine].sum()/nothing_spine
                    aucn[j] = neck_ca[trial][:, t_start:t_start+period_spine].sum()/nothing_spine
                    aucd[j] = dend_ca[trial][:, t_start:t_start+period_dend].sum()/2/nothing_dend
                    
             
                auc_head_m.append(auch.mean())
                auc_neck_m.append(aucn.mean())
                auc_dend_m.append(aucd.mean())
                max_head_m.append(maxh.mean())
                max_neck_m.append(maxn.mean())
                max_dend_m.append(maxd.mean())
                auc_head_e.append(auch.std()/trials**0.5)
                auc_neck_e.append(aucn.std()/trials**0.5)
                auc_dend_e.append(aucd.std()/trials**0.5)
                max_head_e.append(maxh.std()/trials**0.5)
                max_neck_e.append(maxn.std()/trials**0.5)
                max_dend_e.append(maxd.std()/trials**0.5)

              
                
            ax.errorbar(x=max_head_m, y=max_dend_m, xerr=max_head_e,
                    yerr=max_dend_e, label=types[i],
                    marker=markers[i], fillstyle=fillstyles[i],
                    color="tab:blue", linewidth=0)
            ax1.errorbar(x=auc_head_m, y=auc_dend_m, xerr=auc_head_e, yerr=auc_dend_e,label=types[i], marker=markers[i],
                         fillstyle=fillstyles[i], color="tab:blue", linewidth=0)
            ax_n.errorbar(x=auc_neck_m, y=auc_dend_m,xerr=auc_neck_e, yerr=auc_dend_e, label=types[i], marker=markers[i],
                          fillstyle=fillstyles[i], color="tab:blue", linewidth=0)
            ax1_n.errorbar(x=max_neck_m, y=max_dend_m,xerr=max_neck_e, yerr=max_dend_e, label=types[i], marker=markers[i], fillstyle=fillstyles[i], color="tab:blue",
                linewidth=0)
            ax2.errorbar(x=auc_head_m, y=auc_neck_m, xerr=auc_head_e,
                         yerr=auc_neck_e, label=types[i],
                         marker=markers[i], fillstyle=fillstyles[i],
                         color="tab:blue", linewidth=0)

           
            
        ax.set_ylabel("Max dend Ca (nM)")
        ax.set_xlabel("Max head Ca (uM)")
        ax2.set_ylabel("Auc neck Ca/Auc basal")
        ax2.set_xlabel("Auc head Ca/Auc basal")
        
        ax1.set_ylabel("Auc dend Ca/Auc basal")
        ax1.set_xlabel("Auc head Ca/Auc basal")

        ax1_n.set_ylabel("Max dend Ca (nM)")
        ax1_n.set_xlabel("Max neck Ca (uM)")
        ax_n.set_ylabel("Auc dend Ca/Auc basal")
        ax_n.set_xlabel("Auc neck Ca/Auc basal")
      

       

        ax.legend()
        ax1.legend()
        ax_n.legend()
        ax1_n.legend()

        fig.savefig("Just_spine_head_stim_dend_max_vs_head_max_%s_ms.png" % dur, dpi=100, bbox_inches="tight")
        fig2.savefig("Just_spine_head_stim_neck_auc_vs_head_auc_%s_ms.png" % dur, dpi=100, bbox_inches="tight")
        fig1.savefig("Just_spine_head_stim_dend_auc_vs_head_auc_%s_ms.png" % dur, dpi=100, bbox_inches="tight")
        fig_n.savefig("Just_spine_head_stim_dend_auc_vs_neck_auc_%s_ms.png" % dur, dpi=100, bbox_inches="tight")
        fig1_n.savefig("Just_spine_head_stim_dend_max_vs_neck_max_%s_ms.png" % dur, dpi=100, bbox_inches="tight")
        
    
    #plt.show()
