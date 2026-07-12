import os
import subprocess
import math
import re

# 1. Define your specific (ptMin, ptMin_2) pairs
pt_pairs = [
    ## pt_lead = 2
    (2.0, 1.0),
    ## pt_lead = 5
    (5.0, 1.0),
    (5.0, 2.0),
    (5.0, 3.0),
    ## pt_lead = 6
    (6.0, 3.0),
    ## pt_lead = 10
    (10.0, 1.0),
    (10.0, 2.0),
    (10.0, 3.0),
    (10.0, 5.0),
    (10.0, 8.0),
    ## pt_lead = 20
    (20.0, 1.0),
    (20.0, 2.0),
    (20.0, 3.0),
    (20.0, 5.0),
    (20.0, 8.0),
    (20.0, 10.0),
    (20.0, 15.0),
    # pt_lead = pt_sublead
    (1.0, 1.0),
    (2.0, 2.0),
    (3.0, 3.0),
    (5.0, 5.0),
    (8.0, 8.0),
    (10.0, 10.0),
    (15.0, 15.0),
    (20.0, 20.0),

]

original_c_file = "Make2dMassHistos_M1M2_data_10SlicesEach_DDupdates.C"

# Read the original source code once to hold in memory
with open(original_c_file, "r") as f:
    original_code = f.read()

# Setup file counts (only need to calculate this once)
flist_name = "/home/awesole/forDDbarAnalysis/forTemplatesMethod/dataFiles/Run_2023_19Sept.list"
with open(flist_name) as f:
    flist = [line.strip() for line in f if line.strip()]

total_files = len(flist)
files_per_job = 900
njobs = math.ceil(total_files / files_per_job)

print(f"Total files: {total_files} | Files per job: {files_per_job} | SLURM array size: {njobs} jobs\n")

# 2. Iterate directly through your specific pairs
for pt1, pt2 in pt_pairs:
    print(f"--- Starting setup for ptMin={pt1}, ptMin_2={pt2} ---")
    str_pt1 = str(int(pt1))
    str_pt2 = str(int(pt2))
    
    # Modify the C++ code using regex to ensure we target the exact variable assignments
    modified_code = re.sub(r'const float ptMin\s*=\s*[\d\.]+;', f'const float ptMin = {pt1};', original_code)
    modified_code = re.sub(r'const float ptMin_2\s*=\s*[\d\.]+;', f'const float ptMin_2 = {pt2};', modified_code)
    
    temp_c_file = f"Make2dMassHistos_temp_{str_pt1}_{str_pt2}.C"
    exe_name = f"M1M2_data_10SlicesEach_{str_pt1}_{str_pt2}.exe"
    
    with open(temp_c_file, "w") as f:
        f.write(modified_code)
        
    print(f"Compiling {temp_c_file}...")
    compile_cmd = (
        f"g++ {temp_c_file} "
        "../../src/BDTHandler.cc ../../src/phiTools.cc "
        "-I../../include $(root-config --cflags --libs) "
        f"-O2 -Wall -o {exe_name}"
    )
    
    try:
        subprocess.run(compile_cmd, shell=True, check=True)
        print("Compilation successful.")
    except subprocess.CalledProcessError:
        print(f"Compilation failed for {temp_c_file}. exit\n")
        exit(1)
        
    # 4. Safely build directories using os.path.join
    base_dir = os.path.join("/scratch/negishi/awesole", f"DDbar_July_2026_pT_{str_pt1}_{str_pt2}/")
    
    os.makedirs(base_dir, exist_ok=True)
    for sub_dir in ["ROOT", "log", "err"]:
        os.makedirs(os.path.join(base_dir, sub_dir), exist_ok=True)
        
    # 5. Generate a unique SLURM Script for this pair
    slurm_script = f"run_M1M2_{str_pt1}_{str_pt2}.slurm"
    
    with open(slurm_script, "w") as f:
        f.write("#!/bin/bash\n")
        f.write(f"#SBATCH --job-name=M1M2_{str_pt1}_{str_pt2}\n")
        f.write(f"#SBATCH --output={base_dir}log/%A_%a.out\n")
        f.write(f"#SBATCH --error={base_dir}err/%A_%a.err\n")
        f.write("#SBATCH --time=04:00:00\n")
        f.write("#SBATCH --mem=16G\n")
        f.write("#SBATCH --cpus-per-task=1\n")
        f.write(f"#SBATCH --array=0-{njobs-1}\n\n")

        f.write("export X509_USER_PROXY=/home/awesole/myproxy\n\n")

        f.write(f"NFILE={files_per_job}\n")
        f.write("ISTART=$(( SLURM_ARRAY_TASK_ID * NFILE ))\n")
        f.write("IEND=$(( ISTART + NFILE ))\n\n")

        f.write(f"./{exe_name} $ISTART $IEND\n")
        
    # 6. Submit the SLURM job
    submit_cmd = f"sbatch -A cms -q standby -p cpu --time=04:00:00 {slurm_script}"
    print(f"Submitting SLURM array...")
    os.system(submit_cmd)
    
    print(f"Finished setup and submission for {str_pt1}, {str_pt2}.\n")

print("All jobs for the parameter pairs have been successfully submitted!")
