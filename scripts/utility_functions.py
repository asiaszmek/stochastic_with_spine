#!/usr/bin/env python
import os
import h5py
import numpy as np
from lxml import etree
import sys
from scipy.constants import Avogadro
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


plt.rcParams.update({
    "text.usetex": True,
    "font.family": "sans-serif",
    "font.sans-serif": "Helvetica",
})

plt.rcParams["text.latex.preamble"]+=r"\usepackage{sfmath}"
plt.rcParams["text.latex.preamble"]+=r"\usepackage{siunitx}"
plt.rcParams["text.latex.preamble"]+=r"\DeclareSIUnit{\molar}{M}"
plt.rcParams["text.latex.preamble"]+=r"\DeclareSIUnit{\Molar}{M}"

hatch_possibilities = ["/", "-", "+", "o"]
marker = ["d", "o", "v", "^"]
limit = 2.5
NA = Avogadro*1e-23
spine = ['PSD', 'head', 'neck']
t_init = 3000
window = 50

region_list = []
prefix = "dend"
ca_ctrl = 71
for i in range(1, 52):
    if i<10:
        region_list.append(prefix+"0"+str(i))
    else:
        region_list.append(prefix+str(i))


def pretty_axis(ax1):
    max_ylim = max([max(ax.get_ylim()) for ax in ax1])
    max_xlim = max([max(ax.get_xlim()) for ax in ax1])
    for x in ax1:
        x.set_ylim([-0.05, max_ylim+0.05])
        x.set_xlim([-0.05, max_xlim+0.05])
        

        

def get_array(conc_dict, specie):
    mini  = min([conc_dict[specie][key].shape[-1]
                 for key in conc_dict[specie].keys()])

    no = len(conc_dict[specie].keys())
    voxels = conc_dict[specie]["trial0"].shape[0]
    out = np.zeros((no, voxels, mini))
    for i, trial in enumerate(conc_dict[specie].keys()):
        out[i] = conc_dict[specie][trial][:,:mini]
    return out
    

def moving_average(x, w):
    return np.convolve(x, np.ones(w), 'valid') / w


def nano_molarity(N, V):
    return 10 * N / V / NA


def pico_sd(N, S):
    return 10 * N / S / NA


def get_length(My_file):
    if isinstance(My_file, str):
        try:
            my_file = h5py.File(My_file)
        except FileNotFoundError:
            return
    else:
        my_file = My_file    
    grid = get_grid_list(my_file)
    return grid[-1][3] - grid[0][0]
    
def get_grid_list(My_file):
    return np.array(My_file['model']['grid'])


def get_indx_no(My_file, spine_name, region_name):
    grid = get_grid_list(My_file)
    no = 0
    for g in grid:
        if g[-1].decode('utf-8') == spine_name and g[-4].decode('utf-8') == region_name:
            no += 1
    return no

def get_times(My_file, trial='trial0', output="__main__"):
    return np.array(My_file[trial]['output'][output]['times'])

def get_outputs(my_file):
    return my_file['model']['output'].keys()


def get_populations(my_file, trial='trial0', output='__main__'):
    return np.array(my_file[trial]['output'][output]['population'])


def get_all_species(My_file, output="__main__"):
    return [s.decode('utf-8') for s in My_file['model']['output'][output]['species']]


def get_all_anchored_species(root):
    all_species = []
    for son in root:
        if son.tag.endswith('ReactionScheme'):
            for grandson in son:
                if grandson.tag.endswith('Specie'):
                    if not float(grandson.get("kdiff")):
                        all_species.append(grandson.get('id'))
    return list(set(all_species))


def get_all_submembrane_species(my_file):
    root = etree.fromstring(my_file['model']['serialized_config'][0])
    all_anchored_species = get_all_anchored_species(root)
    anchored = []
    for son in root:
        if son.tag.endswith('InitialConditions'):
            for grandson in son:
                if grandson.tag.endswith("SurfaceDensitySet"):
                    for grandgrandson in grandson:
                        name = grandgrandson.get("specieID")
                        if name in all_anchored_species:
                            anchored.append(name)
    return list(set(anchored))

def get_output_regions(my_file):
    root = etree.fromstring(my_file['model']['serialized_config'][0])
    outputs = {}
    for son in root:
        if son.tag.endswith('OutputScheme'):
            for grandson in son:
                outputs[grandson.get("filename")] = grandson.get("region")
    return outputs
        
def get_key(cell):
    if cell[18]:
        return cell[15].decode('utf-8') + '_' + cell[18].decode('utf-8')
    return cell[15].decode('utf-8')

def region_volumes(my_file):
    if isinstance(my_file, str):
        my_file = h5py.File(my_file)
    grid_list = get_grid_list(my_file)
    regions = get_regions(my_file)
    volumes = {}
    for region in regions:
        volumes[region] = 0
    for cell in grid_list:
        key = get_key(cell)
        volumes[key] += float(cell[12])

    return volumes


def sum_volume(my_file, region_list):
    grid_list = get_grid_list(my_file)
    vol_sum = 0
    volumes = region_volumes(my_file)
    for region in region_list:
        if region in volumes:
            vol_sum += volumes[region]
    return vol_sum


def sum_indices(my_file, region_list):
    reg_indices = get_region_indices(my_file)
    sum_indices = []
    for region in region_list:
        if region in reg_indices:
            sum_indices += reg_indices[region]
    return sum_indices


def region_surface(grid_list, direction=0):
    submembrane_regions = []
    submembrane_regions_dict = {}
    for i, cell in enumerate(grid_list):
        if cell[17] == b'submembrane':
            new_name = cell[15].decode('utf-8')
            if  new_name not in submembrane_regions:
                submembrane_regions.append(new_name)
                submembrane_regions_dict[new_name] = []
            submembrane_regions_dict[new_name].append(i)
    surface = {}
    for key in submembrane_regions_dict:
        surface[key] = 0
        for cell_idx in submembrane_regions_dict[key]:
            if direction == 0:
                depth = grid_list[cell_idx][13]
                width = abs(grid_list[cell_idx][0] - grid_list[cell_idx][3])
                surface[key] += depth * width
            else:
                print('Unimplemented direction', direction)

    return surface


def get_region_indices(my_file):
    grid_list = get_grid_list(my_file)
    region_ind = {}
    for idx, cell in enumerate(grid_list):
        key = get_key(cell)
        if key not in region_ind:
            region_ind[key] = []
        region_ind[key].append(idx)
    return region_ind

def get_spines(regions):
    out = {}
    for region in regions:
        try:
            end = region.split("_")[1]
        except IndexError:
            continue
        if end == "":
            continue

        if end in out:
            out[end].append(region)
        else:
            out[end] = [region]
    return out


def get_regions(my_file):
    grid_list = get_grid_list(my_file)
    return sorted(list(set([get_key(grid) for grid in grid_list])))


def get_concentrations_region_list(my_file, my_list, trial, out):

    grid_list = get_grid_list(my_file)
    species = get_all_species(my_file)
    idxs = sum_indices(my_file, my_list)
    vol = sum_volume(my_file, my_list)
    data = get_populations(my_file, trial=trial, output=out)
    numbers = data[:, idxs, :].sum(axis=1)
    return nano_molarity(numbers, vol)


def get_concentrations(my_file, trial, out):
    grid_list = get_grid_list(my_file)
    data = get_populations(my_file, trial=trial, output=out)
    species = get_all_species(my_file, output=out)
    regions = get_regions(my_file)
    submembrane_species = get_all_submembrane_species(my_file)
    volume_dict = region_volumes(my_file)
    surface_dict = region_surface(grid_list)
    concentrations = np.zeros((data.shape[0], len(regions), len(species)))
    numbers = np.zeros_like(concentrations)
    region_indices = get_region_indices(my_file)

    for i, reg in enumerate(regions):
        # get numbers
        numbers[:, i, :] = data[:, region_indices[reg], :].sum(axis=1)
        
        if reg in surface_dict:
            for j, specie in enumerate(species):
                if specie in submembrane_species:
                    concentrations[:, i, j] = pico_sd(numbers[:, i, j],
                                                      surface_dict[reg])
                else:
                    concentrations[:, i, j] = nano_molarity(numbers[:, i, j],
                                                            volume_dict[reg])
        else:
            concentrations[:, i, :] = nano_molarity(numbers[:, i, :],
                                                    volume_dict[reg])

        
    return concentrations


def save_single_file(times, concentrations, species, fname):
    header = 'time'
    for specie in species:
        header += ' ' + specie
    what_to_save = np.zeros((concentrations.shape[0], len(species) + 1))
    what_to_save[:, 0] = times[:concentrations.shape[0]]
    what_to_save[:, 1:] = concentrations
    np.savetxt(fname, what_to_save, header=header, comments='')


def save_concentrations(my_file, fname_base, output, trial='trial0'):
    specialized_output = get_output_regions(my_file)
    regions = get_regions(my_file)
    times = get_times(my_file, trial=trial, output=output)
    species = get_all_species(my_file, output=output)
    if output == '__main__':
        add = ''
    else:
        add = output + '_'
    spines_dict = get_spines(regions)

    if output not in specialized_output or specialized_output[output] is None:
        concentrations = get_concentrations(my_file, trial, output)
        totals = get_concentrations_region_list(my_file, regions, trial, output)
        save_single_file(times, totals, species,
                         '%s_%s%s_%s.txt' % (fname_base, add, trial, 'total'))
        
        for spine_name in spines_dict.keys():
            spine_reg = spines_dict[spine_name]
            spine = get_concentrations_region_list(my_file, spine_reg,
                                                   trial, output)
            save_single_file(times, spine, species,
                             '%s_%s%s_%s.txt' % (fname_base, add,
                                                 trial, spine_name))
            print('%s_%s%s_%s.txt' % (fname_base, add,
                                      trial, spine_name))
        for i, region in enumerate(regions):
            fname = '%s_%s%s_%s.txt' % (fname_base, add, trial, region)
            print(fname)
            save_single_file(times, concentrations[:, i, :], species, fname)

    else:
        region = specialized_output[output]
        data = get_populations(my_file, trial=trial, output=output)
        start = 0
        if region in ["PSD", "head", "neck"]:
            for i, spine in enumerate(spines_dict.keys()):
                my_region = "%s_%s" % (region, spine)
                fname = '%s_%s_%s_%s.txt' % (fname_base, output, trial,
                                             my_region)
                idxs = get_indx_no(my_file, spine, region)

                volume = region_volumes(my_file)[my_region]
                concentrations = np.zeros((data.shape[0], len(species)))
                
                for j, specie in enumerate(species):
                    concentrations[:, j] = nano_molarity(data[:,
                                                              start:start+idxs,
                                                              j].sum(axis=1),
                                                         volume)
                save_single_file(times, concentrations, species, fname)
                print(fname)
                start += idxs


def get_dend_indices(grid, region=["dend"]):
    out = {}
    volumes = {}
    if not isinstance(region, list):
        region = [region]
    for i, line in enumerate(grid):
        if line[15].decode('utf-8') in region:
            pos = abs(np.round(line[0], 3))
            if pos in out:
                out[pos].append(i)
            else:
                out[pos] = [i]
            if pos in volumes:
                volumes[pos] += line[12]
            else:
                volumes[pos] = line[12]
    return out, volumes


def get_dynamics_in_region(my_file, specie, region, trial,
                           output="__main__"):
    if not isinstance(specie, list):
        specie = [specie]
    my_grid = get_grid_list(my_file)
    vox_ind, vols = get_dend_indices(my_grid, region=region)
    specie_list = get_all_species(my_file, output)
    population = get_populations(my_file, trial, output)
    specie_idx = []
    for sp in specie:
        specie_idx.append(specie_list.index(sp))
    voxel_list = sorted(vox_ind.keys())
    how_many_voxels = len(voxel_list)
    out = np.zeros((population.shape[0], how_many_voxels))
    for i, key in enumerate(voxel_list):
        volume = vols[key]
        for idx in specie_idx:
            h = population[:, vox_ind[key], idx].sum(axis=1)
            out[:, i] = out[:, i] + nano_molarity(h, volume)
    return out, voxel_list


def get_conc(my_file, specie, region_list, output):
    conc_dict = {}
    time_dict = {}
    for trial in my_file.keys():
        if trial == "model":
            continue
        try:
            pop, voxel = get_dynamics_in_region(my_file, specie,
                                                region_list,
                                                trial, output)
            conc_dict[trial] = pop.T
            time = get_times(my_file, trial, output)
            time_dict[trial] = time
        except IOError:
            print("Something wrong with", my_file)
            break
    return conc_dict, time_dict


def get_distance(conc_dict, dt, t_init, stim_len, spine_idx, rollo=True,
                 interval=700):
    
    decays = np.zeros((len(conc_dict), 1))
    shape = conc_dict["trial0"].shape[0]
    find_max = conc_dict["trial0"].argmax()
    full_shape = conc_dict["trial0"].shape
    
    for i, concentration in enumerate(conc_dict.values()):
        ca_conc = np.zeros((shape,))
        
        
        new_beg = int((t_init)/dt)
        indices = []
        for j in range(spine_idx, shape):
            
            try:
                new_idx = concentration[j, new_beg-100:new_beg + stim_len
                                        + interval].argmax()
            except ValueError:
                continue

            ca_conc[j] = concentration[j, new_beg-100+new_idx]

            if ca_conc[j] > limit*ca_ctrl and j==spine_idx:
                indices.append(j)
            elif ca_conc[j] > limit*ca_ctrl and j-1 in indices:
                    indices.append(j)
            elif ca_conc[j] > limit*ca_ctrl and j-2 in indices:
                    indices.append(j)
                    indices.append(j-1)
            else:
                break
            print(j, ca_conc[j], new_beg +new_idx -100-int((t_init)/dt))
            if rollo:
                new_beg = new_beg +new_idx -100
    
        new_beg = int((t_init)/dt)         
        for j in range(spine_idx-1, -1, -1):
            try:
                new_idx = concentration[j, new_beg-100:new_beg
                                        +stim_len+interval].argmax()
            except ValueError:
                continue
            ca_conc[j] = concentration[j, new_beg-100+new_idx]
         
            if ca_conc[j] > limit*ca_ctrl and j==spine_idx-1:
                indices.append(j)
                
            elif ca_conc[j] > limit*ca_ctrl and j+1 in indices:
                indices.append(j)
            elif ca_conc[j] > limit*ca_ctrl and j+2 in indices:
                indices.append(j)
                indices.append(j+1)
            else:
                break
            print(j, ca_conc[j], new_beg +new_idx -100-int((t_init)/dt))
            if rollo:
                new_beg = new_beg + new_idx -100
        decays[i] = len(indices)/2*0.5
        #print(indices)
    return decays



def max_vs_distance(conc, dt, t_init, spine_idx=49, length=51):
    out_conc = np.zeros((len(conc.values()), spine_idx))
    start = int(t_init/dt)
    for i, (key, concentration) in enumerate(conc.items()):
        for j in range(spine_idx):
            #print(spine_idx+j, spine_idx-j, out_conc.shape[1], len(concentration))
            out_conc[i, j] = max(max(concentration[spine_idx+j, start:]),
                                 max(concentration[spine_idx-j, start:]))
    
    return out_conc    


def get_max_basal(conc, dt, t_init, spine_idx, length):
    basal_start = int(1000/dt)
    start = int(t_init/dt)
    basal = []
    for i, concentration in enumerate(conc.values()):
        new_basal = np.zeros((spine_idx,))
        for j in range(spine_idx):
            new_basal[j] = max(max(concentration[spine_idx-j, basal_start:start]),
                               max(concentration[spine_idx+j, basal_start:start]))
        basal.append(new_basal)

       
    return np.array(basal)

def mean_vs_distance(conc, dt, t_init, period=500,  length=51):
    out_conc = np.zeros((len(conc.values()), spine_idx))
    start = int(t_init/dt)
    for i, (key, concentration) in enumerate(conc.items()):
        find_max = concentration.argmax()
        full_shape = concentration.shape
        spine_idx = np.unravel_index(find_max, full_shape)[0]
        for j in range(spine_idx):
            #print(spine_idx+j, spine_idx-j, out_conc.shape[1], len(concentration))
            out_conc[i, j] = (concentration[spine_idx+j,
                                            start:start+period].mean()+
                              concentration[spine_idx-j,
                                            start:start+period].mean())/2
    return out_conc    


def get_mean_basal(conc, dt, t_init, length):
    basal_start = int(1000/dt)
    start = int(t_init/dt)
    basal = []
    for i, concentration in enumerate(conc.values()):
        find_max = concentration.argmax()
        full_shape = concentration.shape
        spine_idx = np.unravel_index(find_max, full_shape)[0]
        new_basal = np.zeros((spine_idx,))
        for j in range(spine_idx):
            
            new_basal[j] = (np.mean(concentration[spine_idx-j, basal_start:start])+
                               np.mean(concentration[spine_idx+j, basal_start:start]))/2
        basal.append(new_basal)

       
    return np.array(basal)

def fit_exp(time, ca_conc, dt, duration=2000, t_init=3000, stim_len=3000, spatial=False):

    if not spatial:
        ca_conc_mean = ca_conc[:int(t_init/dt)].mean()
        ca_conc = ca_conc[int((t_init+stim_len)/dt):] - ca_conc_mean
        
        time = time[int((t_init+stim_len)/dt):] - t_init - stim_len
        
        try:
            ca_conc_log = ca_conc[:duration]
            new_time = time[:duration]
            
        except IndexError:
            return 0
        try:
            popt, pcov = curve_fit(lambda t, a, b, c: a*np.exp(-t/b)+c,
                                   new_time, ca_conc_log)
        except RuntimeError:
            return 0
        
        return popt[1]
    popt, pcov = curve_fit(lambda t, a, b, c: a*np.exp(-abs(t)/b)+c,
                           time, ca_conc-(ca_conc[0]+ca_conc[-1])/2)

    return popt[1]





def make_distance_figs(directories, stim_dict, output, types, markers,
                       fillstyles, length=42, spine_idx=21, t_init=3000,
                       xlabel=r"Auc head/Auc basal",
                       ylabel=r"Spatial spread ($\unit{\micro\metre}$))"):
    fig1, ax1 = plt.subplots(1, len(stim_dict.keys()),
                             figsize=(len(stim_dict.keys())*5, 5))
    for k, dur in enumerate(stim_dict.keys()):
        
        
        
        for i, fname in  enumerate(directories):
            auc_head_m = []
            auc_head_e = []
            spread_m = []
            spread_e = []
            for inj in stim_dict[dur]:
                print(fname % (dur, inj))
                try:
                    my_data = h5py.File(fname % (dur, inj))
                except FileNotFoundError:
                    print("Cound not open ", fname % (dur, inj))
                    continue
        
                dend_ca, time_dict = get_conc(my_data, "Ca", region_list, output)
                head_ca, time_head_dict = get_conc(my_data, "Ca", ['head'], output)
                dt = time_dict["trial0"][1]-time_dict["trial0"][0]
                t_start = int(t_init/dt)
               
                length = time_dict["trial0"][-1] - t_start
         
                trials = len(dend_ca.keys())
                auch = np.zeros((trials))
                for j, trial in enumerate(dend_ca.keys()):
                    auch[j] = head_ca[trial][1, t_start:].sum()/(71*len(head_ca[trial][1, t_start:]))

                spread = get_distance(dend_ca, dt, t_init=t_init, stim_len=int(dur),
                                      spine_idx=spine_idx)
                auc_head_m.append(auch.mean())
                auc_head_e.append(auch.std()/trials**0.5)
                spread_m.append(spread.mean())
                spread_e.append(spread.std()/trials**0.5)
                
            print(auc_head_m, spread_m, spread_e)
            ax1[k].errorbar(auc_head_m, spread_m, yerr=spread_e, xerr=auc_head_e,label=types[i], marker=markers[i],
                            fillstyle=fillstyles[i], color="tab:blue", linewidth=0)
        ax1[k].set_title("Stimulation %s (ms)"% dur)      

           
            
        ax1[k].set_ylabel(ylabel)
        ax1[k].set_xlabel(xlabel)
    pretty_axis(ax1)
    ax1[0].legend()
    return fig1


def max_vs_auc_head_neck_dend(directories, stim_dict, output, types, markers,
                       fillstyles, length=42, which_dend=["dend11"], t_init=3000):
    fig_dh, ax_dh = plt.subplots(1, len(stim_dict.keys()),
                                 figsize=(len(stim_dict.keys())*5, 5))
    fig_max_dh, ax_max_dh = plt.subplots(1, len(stim_dict.keys()),
                                          figsize=(len(stim_dict.keys())*5, 5))
    fig_nh, ax_nh = plt.subplots(1, len(stim_dict.keys()),
                                 figsize=(len(stim_dict.keys())*5, 5))
    for k, dur  in enumerate(stim_dict.keys()):
         
        for i, fname in enumerate(directories):
            auc_head_m = []
            auc_neck_m = []
            auc_dend_m = []
            max_head_m = []
            max_dend_m = []
            auc_head_e = []
            auc_neck_e = []
            auc_dend_e = []
            max_head_e = []
            max_dend_e = []
            for inj in stim_dict[dur]:
              
                print(fname % (dur, inj))
                try:
                    my_data = h5py.File(fname % (dur, inj))
                except FileNotFoundError:
                    print("Cound not open ", fname % (dur, inj))
                    continue
        
                dend_ca, time_dict = get_conc(my_data, "Ca", which_dend, output)
                head_ca, time_head_dict = get_conc(my_data, "Ca", ['head'], output)
                neck_ca, time_head_dict = get_conc(my_data, "Ca", ['neck'], output)
               
                dt = time_dict["trial0"][1]-time_dict["trial0"][0]
                t_start = int(t_init/dt)
               
                length = time_dict["trial0"][-1] - t_start
                
                trials = len(dend_ca.keys())
                auch = np.zeros((trials))
                aucn = np.zeros((trials))
                aucd = np.zeros((trials))
                maxh = np.zeros((trials))
                maxd = np.zeros((trials))
              

                for j, trial in enumerate(dend_ca.keys()):
                    try: 
                        maxh[j] = head_ca[trial][1, t_start:].max()/1000
                    except ValueError:
                        break
                    maxd[j] = dend_ca[trial][:, t_start:].mean(axis=0).max()
                    auch[j] = head_ca[trial][1, t_start:].sum()/71/len(head_ca[trial][1, t_start:])
                    aucn[j] = neck_ca[trial][:, t_start:].sum()/71/len(head_ca[trial][1, t_start:])
                    aucd[j] = dend_ca[trial][:, t_start:].sum()/2/71/len(head_ca[trial][1, t_start:])
                    
             
                auc_head_m.append(auch.mean())
                auc_neck_m.append(aucn.mean())
                auc_dend_m.append(aucd.mean())
                max_head_m.append(maxh.mean())
                max_dend_m.append(maxd.mean())
                auc_head_e.append(auch.std()/trials**0.5)
                auc_neck_e.append(aucn.std()/trials**0.5)
                auc_dend_e.append(aucd.std()/trials**0.5)
                max_head_e.append(maxh.std()/trials**0.5)
                max_dend_e.append(maxd.std()/trials**0.5)

              
                
            ax_max_dh[k].errorbar(x=max_head_m, y=max_dend_m, xerr=max_head_e,
                                  yerr=max_dend_e, label=types[i],
                                  marker=markers[i], fillstyle=fillstyles[i],
                                  color="tab:blue", linewidth=0)
            ax_dh[k].errorbar(x=auc_head_m, y=auc_dend_m, xerr=auc_head_e, yerr=auc_dend_e,label=types[i], marker=markers[i],
                              fillstyle=fillstyles[i], color="tab:blue", linewidth=0)
            ax_nh[k].errorbar(x=auc_head_m, y=auc_neck_m, xerr=auc_head_e,
                              yerr=auc_neck_e, label=types[i],
                              marker=markers[i], fillstyle=fillstyles[i],
                              color="tab:blue", linewidth=0)

           
        if not k:
            ax_max_dh[k].set_ylabel(r"Max dend Ca ($\unit{\nano\Molar}$)")
            ax_dh[k].set_ylabel("Auc dend Ca/Auc basal")
            ax_nh[k].set_ylabel("Auc neck Ca/Auc basal")                
        ax_max_dh[k].set_xlabel(r"Max head Ca ($\unit{\micro\Molar}$)")
        ax_dh[k].set_xlabel("Auc head Ca/Auc basal")
        ax_nh[k].set_xlabel("Auc head Ca/Auc basal")
        ax_max_dh[k].set_title("Stimulation %s (ms)"% dur)
        ax_dh[k].set_title("Stimulation %s (ms)"% dur)
        ax_nh[k].set_title("Stimulation %s (ms)"% dur)
    pretty_axis(ax_max_dh)
    pretty_axis(ax_dh)
    pretty_axis(ax_nh)
    ax_max_dh[0].legend()         
    ax_dh[0].legend()
    ax_nh[0].legend()
    return fig_dh, fig_max_dh, fig_nh
