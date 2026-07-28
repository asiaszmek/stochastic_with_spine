import os
import glob
import subprocess


file_list = glob.glob("model_1.2*xml")
print(file_list)    

for fname in file_list:
  
    new_name = fname.replace("model_1.2", "model_Ca_dend_1.2")
    subprocess.run(["mv", fname,
                    new_name],
                   capture_output=True)
