#!/usr/bin/env python
import os
import h5py
import numpy as np
from lxml import etree
import sys
from scipy.constants import Avogadro
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



plt.rcParams['text.usetex'] = True
hatch_possibilities = ["/", "-", "+", "o"]
marker = ["d", "o", "v", "^"]
limit = 2.5
NA = Avogadro*1e-23
spine = ['PSD', 'head', 'neck']
t_init = 3000
window = 50

region_list = []
prefix = "dend"
for i in range(1, 52):
    if i<10:
        region_list.append(prefix+"0"+str(i))
    else:
        region_list.append(prefix+str(i))


     

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
    print(fname)
    np.savetxt(fname, what_to_save, header=header, comments='')


def save_concentrations(my_file, fname_base, output, trial='trial0'):
    regions = get_regions(my_file)
    times = get_times(my_file, trial=trial, output=output)
    species = get_all_species(my_file, output=output)
    concentrations = get_concentrations(my_file, trial, output)
    if output == '__main__':
        add = ''
    else:
        add = output + '_'
    for i, region in enumerate(regions):
        fname = '%s_%s%s_%s.txt' % (fname_base, add, trial, region)
        save_single_file(times, concentrations[:, i, :], species, fname)
    if len(regions) > 1:
        totals = get_concentrations_region_list(my_file, regions, trial, output)
        save_single_file(times, totals, species,
                         '%s_%s%s_%s.txt' % (fname_base, add, trial, 'total'))
        spines_dict = get_spines(regions)
        for spine_name in spines_dict.keys():
            spine_reg = spines_dict[spine_name]
            spine = get_concentrations_region_list(my_file, spine_reg,
                                                   trial, output)
            save_single_file(times, spine, species,
                             '%s_%s%s_%s.txt' % (fname_base, add,
                                                 trial, spine_name))


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


def get_conc(my_file, specie_list, region_list, output):
    if isinstance(specie_list, str):
        specie_list = [specie_list]
    conc_dict = {}
    time_dict = {}
    for specie in specie_list:
        conc_dict[specie] = {}
    for trial in my_file.keys():
        if trial == "model":
            continue
        try:
            for specie in specie_list:
                pop, voxel = get_dynamics_in_region(my_file, specie,
                                                    region_list,
                                                    trial, output)
                conc_dict[specie][trial] = pop.T
            time = get_times(my_file, trial, output)
            time_dict[trial] = time
        except IOError:
            print("Something wrong with", my_file)
            break
    return conc_dict, time_dict


def get_distance(conc_dict, dt, t_init=3000, length=102, spine_idx=49):
    decays = np.zeros((len(conc_dict), 1))
    shape = conc_dict["trial0"].shape[0]
    for i, concentration in enumerate(conc_dict.values()):
        ca_conc = np.zeros((shape,))
        ca_conc_mean = concentration[:, int(1000/dt):int(t_init/dt)].mean(axis=1)
        new_beg = int(t_init/dt)
        indices = []
        for j in range(spine_idx, length):
            try:
                new_idx = concentration[j, new_beg:].argmax()
            except ValueError:
                continue
            ca_conc[j] = concentration[j, new_beg+new_idx]
            
            if ca_conc[j] > limit*ca_conc_mean[j]:
                if not len(indices):
                    indices.append(j)
                elif j+1 in indices or j-1 in indices:
                    indices.append(j)
        new_beg = int(t_init/dt)         
        for j in range(spine_idx, -1, -1):
            try:
                new_idx = concentration[j, new_beg:].argmax()
            except ValueError:
                continue
            ca_conc[j] = concentration[j, new_beg+new_idx]
            if ca_conc[j] > limit*ca_conc_mean[j]:
                if not len(indices):
                    indices.append(j)
                elif j+1 in indices or j-1 in indices:
                    indices.append(j)
       
        decays[i] = len(indices)/2
    return decays



def max_vs_distance(conc, dt, t_init, spine_idx=49, length=51):
    out_conc = np.zeros((len(conc.values()), length-spine_idx-1))
    start = int(t_init/dt)
    for i, (key, concentration) in enumerate(conc.items()):
        for j in range(length-spine_idx-1):
            #print(spine_idx+j, spine_idx-j, out_conc.shape[1], len(concentration))
            out_conc[i, j] = max(max(concentration[spine_idx+j, start:]),
                                 max(concentration[spine_idx-j, start:]))
    
    return out_conc    


def get_max_basal(conc, dt, t_init, spine_idx, length):
    basal_start = int(1000/dt)
    start = int(t_init/dt)
    basal = []
    for i, concentration in enumerate(conc.values()):
        max_len = concentration.shape[0]-spine_idx+spine_idx-1
        new_basal = np.zeros((length-spine_idx-1,))
        for j in range(length-spine_idx-1):
            new_basal[j] = max(max(concentration[spine_idx-j, basal_start:start]),
                               max(concentration[spine_idx+j, basal_start:start]))
        basal.append(new_basal)

       
    return np.array(basal)


def make_distance_fig(files, t_init, dend_diam,
                      what_species, output_name, colors, types, markers):
    fig1, ax1 = plt.subplots(1, len(dend_diam), figsize=(len(dend_diam)*5, 5))
    for i, diam in enumerate(dend_diam):
        for j, fname in enumerate(files):
            new_fname = fname % (diam)
            try:
                my_file = h5py.File(new_fname)
                print(new_fname)
            except FileNotFoundError:
                print("File not found", new_fname)
                continue
            conc_dict, times_dict = get_conc(my_file,
                                             what_species,
                                             region_list,
                                             output_name)
            try:
                dt = times_dict["trial0"][1]-times_dict["trial0"][0]
            except KeyError:
                continue
            
            length = conc_dict["Ca"]["trial0"].shape[0]
            spine_idx = (length-1)//2
            grid = np.array(get_grid_list(my_file))
            dx = abs(grid[0][0]-grid[0][3])
            conc = max_vs_distance(conc_dict["Ca"],
                                   dt=dt, t_init=t_init, spine_idx=spine_idx, length=length)
            basal = get_max_basal(conc_dict["Ca"], dt=dt, t_init=t_init, spine_idx=spine_idx,
                              length=length)
            mean_basal = basal.mean(axis=0)
            distance = np.linspace(0, len(mean_basal)*dx, len(mean_basal))
            std_basal = basal.std(axis=0)/(basal.shape[0]**0.5)

            max_conc = conc.mean(axis=0)
            
            std_conc = conc.std(axis=0)/(len(conc)**0.5)
            ax1[i].errorbar(distance,
                            max_conc, yerr=std_conc, marker=markers[j], linestyle="",
                            color=colors[diam], fillstyle="full", label=types[j]+" after stim")
            ax1[i].errorbar(distance, mean_basal, yerr=std_basal, color="k",
                            fillstyle="none", marker=markers[j], label=types[j]+" before stim",
                            linestyle="")

        ax1[i].tick_params(axis='x', labelsize=15)
        ax1[i].tick_params(axis='y', labelsize=15)
    ax1[0].set_ylabel("$\mathrm{\max Ca_i\, (nM)}$", fontsize=15)
    ax1[0].set_xlabel(r"Distance from stim $\mu \mathrm{m}$", fontsize=15)
    mini = min([min(x.get_ylim()) for x in ax1])
    maxi = max([max(x.get_ylim()) for x in ax1])
    ax1[0].legend()
    for i, diam in enumerate(dend_diam):
        
        ax1[i].set_title(r"dend diam %s $\mathrm{\mu  m}$" % diam,
                         fontsize=15)
        ax1[i].set_ylim([mini, maxi])
        if i:
            ax1[i].set_yticks([])
    return fig1


def make_distance_fig_compare_with_mean(files, t_init, dend_diam,
                                        what_species, output_name, colors, types, markers):
    fig1, ax1 = plt.subplots(1, len(dend_diam), figsize=(len(dend_diam)*5, 5))
    for i, diam in enumerate(dend_diam):
        means = []
        stds = []
        paradigm = []
        for j, fname in enumerate(files):
            new_fname = fname % (diam)
            try:
                my_file = h5py.File(new_fname)
                print(new_fname)
            except FileNotFoundError:
                print("File not found", new_fname)
                continue
            conc_dict, times_dict = get_conc(my_file,
                                             what_species,
                                             region_list,
                                             output_name)
            try:
                dt = times_dict["trial0"][1]-times_dict["trial0"][0]
            except KeyError:
                continue
            length = conc_dict["Ca"]["trial0"].shape[0]
            spine_idx = (length-1)//2
            grid = np.array(get_grid_list(my_file))
            dx = abs(grid[0][0]-grid[0][3])
            
            distance = get_distance(conc_dict["Ca"], dt, t_init,
                                    length=length, spine_idx=spine_idx)*dx
            means.append(distance.mean())
            stds.append(distance.std()/(len(distance)**.5))
            paradigm.append(types[j])
        print(paradigm)
        ax1[i].errorbar(paradigm, means, yerr=stds, marker="s", linestyle="",
                        color=colors[diam], fillstyle="full")

        ax1[i].tick_params(axis='x', labelsize=15)
        ax1[i].tick_params(axis='y', labelsize=15)
        ax1[i].set_xticklabels(types, rotation=90)
    ax1[0].set_ylabel("Spatial extent $(\mathrm{\mu m})$",  fontsize=15)
    ax1[0].set_xlabel("Paradigm", fontsize=15)
    mini = min([min(x.get_ylim()) for x in ax1])
    maxi = max([max(x.get_ylim()) for x in ax1])
    ax1[0].legend()
    for i, diam in enumerate(dend_diam):
        
        ax1[i].set_title(r"dend diam %s $\mathrm{\mu  m}$" % diam,
                         fontsize=15)
        ax1[i].set_ylim([mini, maxi])
        if i:
            ax1[i].set_yticks([])
    return fig1



def fit_exp(time, ca_conc, dt, duration=2000, t_init=3000, stim_len=3000, spatial=False):
    if not spatial:
        ca_conc_mean = ca_conc[:int(t_init/dt)].mean()
        ca_conc = ca_conc[int((t_init+stim_len)/dt):] - ca_conc_mean
        
        time = time[int((t_init+stim_len)/dt):] - t_init - stim_len
        max_idx = int((t_init+stim_len/dt))+ca_conc[int((t_init+stim_len)/dt):].argmax()
        try:
            ca_conc_log = ca_conc[max_idx:duration+max_idx]
            new_time = time[max_idx:duration+max_idx]-time[max_idx]
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


def make_decay_fig(files, t_init, dend_diam,
                   what_species, output_name, colors, types, markers):
    stim_len = 3000
    fig1, ax1 = plt.subplots(1, len(dend_diam), figsize=(len(dend_diam)*5, 5))
    for i, diam in enumerate(dend_diam):
        means = []
        stds = []
        paradigm = []
        for j, fname in enumerate(files):
            new_fname = fname % (diam)
            try:
                my_file = h5py.File(new_fname)
                print(new_fname)
            except FileNotFoundError:
                print("File not found", new_fname)
                continue
            conc_dict, times_dict = get_conc(my_file,
                                             what_species,
                                             region_list=spine,
                                             output="all")
            try:
                dt = times_dict["trial0"][1]-times_dict["trial0"][0]
            except KeyError:
                continue
            t_decays1 = np.zeros((len(conc_dict["Ca"].keys())))
            for k, trial in enumerate(conc_dict["Ca"].keys()):
                
                ca = conc_dict["Ca"][trial].sum(axis=0)
                time = times_dict[trial]
                
                try:
                    t1 = fit_exp(time, ca, dt, t_init=t_init, stim_len=stim_len)
                except ValueError:
                    continue
                if t1 > 0 and t1 < 1000:
                    t_decays1[k] = t1
                    
            means.append(t_decays1.mean())
            stds.append(t_decays1.std()/(len(t_decays1)**.5))
            paradigm.append(types[j])
        print(means, stds, paradigm)
        ax1[i].errorbar(paradigm, means, yerr=stds, marker="s", linestyle="",
                        color=colors[diam], fillstyle="full")

        ax1[i].tick_params(axis='x', labelsize=15)
        ax1[i].tick_params(axis='y', labelsize=15)
        ax1[i].set_xticklabels(types, rotation=90)
    ax1[0].set_ylabel("Time decay $(\mathrm{ms})$",  fontsize=15)
    ax1[0].set_xlabel("Paradigm", fontsize=15)
    mini = min([min(x.get_ylim()) for x in ax1])
    maxi = max([max(x.get_ylim()) for x in ax1])
    ax1[0].legend()
    for i, diam in enumerate(dend_diam):
        
        ax1[i].set_title(r"dend diam %s $\mathrm{\mu  m}$" % diam,
                         fontsize=15)
        ax1[i].set_ylim([mini, maxi])
        if i:
            ax1[i].set_yticks([])
    return fig1

