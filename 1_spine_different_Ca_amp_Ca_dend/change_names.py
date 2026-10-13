import os
import glob
import subprocess


files = glob.glob("model*xml")
print(files)
for fname in files:
    if "400" in fname:
        new_fname = fname.replace("400", "100EPSP_100Hz")
        if "_017.5" in new_fname:
            new_fname = fname.replace("_017.5", "_1_stim")
        elif "_035" in new_fname:
            new_fname = fname.replace("_035", "_5_stim")
        elif "_070" in new_fname:
            new_fname = fname.replace("_070", "_10_stim")
        elif "_105" in new_fname:
            new_fname = fname.replace("_105", "_15_stim")
        elif "_200" in new_fname:
            new_fname = fname.replace("_200", "_20_stim")
    elif "40" in fname:
        new_fname = fname.replace("40", "4EPSP_100Hz")
        if "_0175" in new_fname:
            new_fname = fname.replace("_0175", "_1_stim")
        elif "_0350" in new_fname:
            new_fname = fname.replace("_0350", "_5_stim")
        elif "_0700" in new_fname:
            new_fname = fname.replace("_0700", "_10_stim")
        elif "_1050" in new_fname:
            new_fname = fname.replace("_1050", "_15_stim")
        elif "_2000" in new_fname:
            new_fname = fname.replace("_2000", "_20_stim")
    elif "4" in fname:
        new_fname = fname.replace("4", "1EPSP")
        if "_01750" in new_fname:
            new_fname = fname.replace("_01750", "_1_stim")
        elif "_03500" in new_fname:
            new_fname = fname.replace("_03500", "_5_stim")
        elif "_07000" in new_fname:
            new_fname = fname.replace("_07000", "_10_stim")
        elif "_10500" in new_fname:
            new_fname = fname.replace("_10500", "_15_stim")
        elif "_20000" in new_fname:
            new_fname = fname.replace("_20000", "_20_stim")

    print(fname, new_fname)
    subprocess.run(["git","mv",  fname,
                    new_fname],
                   capture_output=True)
