import sys
import subprocess
from datetime import date
from lxml import etree
import h5py
import numpy as np
#import matplotlib.pyplot as plt



model_text = """<?xml version="1.0" encoding="UTF-8" standalone="yes"?>
<SDRun xmlns:xi="http://www.w3.org/2001/XInclude" xmlns="http://stochdiff.textensor.org">
    <xi:include href="../Rxn_small_neck.xml" />
    <xi:include href="small_comp.xml" />
    <xi:include href="IC_small_neck_%d_CaM_%d_RyR.xml" />
 
    <!--2D means the morphology is interpreted like a flatworm, 3D for
roundworms. The 2D case is good for testing as it is easy to visualize the
results (also, 3D may not work yet...)  -->
   
    <geometry>          2D           </geometry>
    <depth2D>           0.2          </depth2D>
    <distribution>      BINOMIAL     </distribution>
    <algorithm>         INDEPENDENT  </algorithm>
    <simulationSeed>    245         </simulationSeed>
    <outputQuantity>NUMBER</outputQuantity>

    <!-- run time for the calculation, milliseconds -->
    <runtime>100000</runtime>

    <!-- set the seed to get the same spines each time testing -->
    <spineSeed>123</spineSeed>

    <discretization>
       <defaultMaxElementSide>1</defaultMaxElementSide>
    </discretization>
    <tolerance>0.005</tolerance>

    <outputInterval>500</outputInterval>

    <calculation>GRID_ADAPTIVE</calculation>

</SDRun>
"""



IC_text = """<?xml version="1.0" encoding="utf-8"?>

<InitialConditions>
  <ConcentrationSet>
  
    <NanoMolarity specieID="Ca" value="75"/>
    <NanoMolarity specieID="CaOut" value="2000000"/>
    <!-- total calbindin = 160 uM -->
    <NanoMolarity specieID="Calbin" value="134888"/>
    <NanoMolarity specieID="CalbinC" value="15102"/>
    <NanoMolarity specieID="fixedbuffer" value="1996864"/>
    <NanoMolarity specieID="fixedbufferCa" value="3126"/>

    <NanoMolarity specieID="CaM" value="%f"/>
    <NanoMolarity specieID="RyR3"      value="%f"    />
  </ConcentrationSet>
  <SurfaceDensitySet>
    <PicoSD specieID="CaOutLeak" value="1830"/>
    <PicoSD specieID="ncx" value="2780"/>
    <PicoSD specieID="ncxCa" value="220"/>
    <PicoSD specieID="pmca1a" value="4419"/>
    <PicoSD specieID="pmca1aCa" value="58"/>
    <PicoSD specieID="CaMCa4pmca1a" value="17"/> 
    <PicoSD specieID="CaMCa4pmca1aCa" value="2"/>
    <PicoSD specieID="CaMpmca1a" value="147"/> 
  </SurfaceDensitySet>
</InitialConditions> 
"""





def get_key(cell):
    if cell[18]:
        return cell[15].decode('utf-8') + '_' + cell[18].decode('utf-8')
    return cell[15].decode('utf-8')

def get_regions(my_file):
    grid_list = get_grid_list(my_file)
    return sorted(list(set([get_key(grid) for grid in grid_list])))

def region_volumes(my_file):
    grid_list = get_grid_list(my_file)
    regions = get_regions(my_file)
    volumes = {}
    for region in regions:
        volumes[region] = 0
    for cell in grid_list:
        key = get_key(cell)
        volumes[key] += float(cell[12])
    return volumes

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

def get_output_regions(my_file):
    root = etree.fromstring(my_file['model']['serialized_config'][0])
    outputs = {}
    for son in root:
        if son.tag.endswith('OutputScheme'):
            for grandson in son:
                outputs[grandson.get("filename")] = grandson.get("region")
    return outputs

def sum_volume(my_file, region_list):
    grid_list = get_grid_list(my_file)
    vol_sum = 0
    volumes = region_volumes(my_file)
    for region in region_list:
        if region in volumes:
            vol_sum += volumes[region]
    return vol_sum

def get_all_open(data, species):
    length = len(data[:, 0, 0])//2
    if length %2:
        beg = 0
    else:
        beg = 1
    state = np.zeros(length)
    
    for specie in species:
        if "O" in specie:
            state +=  data[length+beg:, 0, species.index(specie)]
    
    count = len(np.where((state[1:] - state[0:-1])==1)[0])
    sum_times = state.sum()
    return sum_times, count


def get_numbers(my_file, output="__main__"):
    output_dict = get_output_regions(my_file)
    
    free_cams = []
    ryrscam = []
    means_ca = []
    for trial in my_file.keys():
        if trial == "model":
            continue

        times = get_times(my_file, trial=trial, output=output)
        species = get_all_species(my_file, output=output)
        data = get_populations(my_file, trial=trial, output=output)
        dt = times[1]-times[0]
        exp_len = len(times)//2
        vol = sum(region_volumes(my_file).values())
        means_ca.append(data[exp_len:, 0, species.index("Ca")].mean()*10/6.023/vol)
        free_cam = 0
        for specie in ["CaM", "CaMCa2C", "CaMCa2N", "CaMCa4"]:
            free_cam += data[exp_len:, 0, species.index(specie)]
        free_cams.append(free_cam.mean()*10/6.023/vol)
        
        ryrcam = np.zeros((5,))
        for idx, specie in enumerate(species):
            my_cam_no = 0
            if not specie.startswith("RyR3"):
                continue
            out = specie.split("_")
            for factor in out:
                if factor in ["Ca1", "Ca2", "Ca3", "Ca4"]:
                    continue
                if factor == "RyR3":
                    continue
                if "CaM" in factor:
                    my_cam_no += int(factor[0])
            number = data[-1, : , idx].sum()
        
        
            ryrcam[my_cam_no] += number
        ryrscam.append(ryrcam)
    print(ryrscam)
    return np.array(free_cams), np.array(ryrscam), np.array(means_ca)

if __name__ == "__main__":
  
    cam_conc_list = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
                     1500, 2000, 2500, 3000, 3500, 4000, 4500, 5000, 6000, 7000, 8000, 9000, 10000]
    ryr3_conc_list = [100, 150, 200, 250, 300, 350, 400, 450, 500, 550, 600, 650, 700, 750,
                      800, 850, 900, 950, 1000]
    output_4ryrcam = np.zeros((len(cam_conc_list), len(ryr3_conc_list)))
    mean_times = []
    free_cam_out = np.zeros((len(cam_conc_list), len(ryr3_conc_list), 20))
    means_ca_out = np.zeros((len(cam_conc_list), len(ryr3_conc_list), 20))
    for i, cam_conc in enumerate(cam_conc_list):
        for j, ryr3_conc in enumerate(ryr3_conc_list):
            
            IC_name = "IC_small_neck_%d_CaM_%d_RyR.xml" % (cam_conc, ryr3_conc)
            model_name = "CaM_%d_RyR_%d.xml" % (cam_conc, ryr3_conc)
            output_name = "CaM_%d_RyR_%d.h5" % (cam_conc, ryr3_conc)
            
            with open(IC_name, "w") as fic:
                fic.write(IC_text % (cam_conc, ryr3_conc))
            with open(model_name, "w") as fm:
                fm.write(model_text % (cam_conc, ryr3_conc))

            process = subprocess.run(["/usr/lib/jvm/java-8-openjdk-amd64/bin/java",
                                  "-jar",
                                  "/home/jszmek/new_neurord/neurord-3.2.3-all-deps.jar",
                                  "-Dneurord.trials=20", model_name],
                                 capture_output=True)
            print(process.returncode)
            if not process.returncode:
                my_file = h5py.File(output_name, 'r')
                free_cams, ryrscam, mean_ca = get_numbers(my_file, output="__main__")
            
                total_ryr = ryrscam.sum(axis=1)
                ryr_4cam = ryrscam[:, -1]
                print(total_ryr, ryr_4cam)
                output_4ryrcam[i, j] = np.mean(ryr_4cam/total_ryr)
                free_cam_out[i, j] = free_cams
                means_ca_out[i, j] = mean_ca
                print(cam_conc, ryr3_conc, output_4ryrcam[i, j], free_cam_out[i, j], mean_ca)
    # np.save("4CaM_bound_ryrs.npy", output_4ryrcam, cam_conc_list, ryr3_conc_list)
    # np.save("free_cam.npy", free_cam_out, cam_conc_list, ryr3_conc_list)
    # np.save("mean_ca.npy", means_ca_out, cam_conc_list, ryr3_conc_list)
    # fig, ax = plt.subplots(1, 1)
    # im = ax.imshow(output_4ryrcam, interpolation='none', aspect='auto', cmap="viridis", vmin=0, vmax=1,
    #                origin="lower")
    # ax.colorbar(im)
    # fig.savefig("4CaM_bound_ryrs.png", dpi=100,  transparent=False, pad_inches=.5,
    #             bbox_inches=None)
    
    # plt.show()
