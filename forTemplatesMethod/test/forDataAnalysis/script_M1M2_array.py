import os
import subprocess
import math

print("Compiling C++ code...")

compile_cmd = (
    "g++ Make2dMassHistos_M1M2_data_10SlicesEach_DDupdates.C "
    "../../src/BDTHandler.cc ../../src/phiTools.cc "
    "-I../../include $(root-config --cflags --libs) "
    "-O2 -Wall -o M1M2_data_10SlicesEach_3_3_skipping.exe"
)

try:
    subprocess.run(compile_cmd, shell=True, check=True)
    print("Compilation successful.")
except subprocess.CalledProcessError:
    print("Compilation failed.")
    exit(1)

flist_name = "/home/awesole/forDDbarAnalysis/forTemplatesMethod/dataFiles/Run_2023_19Sept.list"
base_dir = "/scratch/negishi/awesole/DDbar_July_2026_pT_3_3/"
if not os.path.exists(base_dir):
    os.makedirs(base_dir)
for sub_dir in ["ROOT", "log", "err"]:
    os.makedirs(f"{base_dir}/{sub_dir}", exist_ok=True)

with open(flist_name) as f:
    flist = [line.strip() for line in f if line.strip()]

total_files = len(flist)
files_per_job = 600
njobs = math.ceil(total_files / files_per_job)

print(f"Total files: {total_files}")
print(f"Files per job: {files_per_job}")
print(f"SLURM array size: {njobs} jobs")

slurm_script = "run_M1M2.slurm"


with open(slurm_script, "w") as f:
    f.write("#!/bin/bash\n")
    f.write("#SBATCH --job-name=M1M2\n")
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

    f.write("./M1M2_data_10SlicesEach_3_3_skipping.exe $ISTART $IEND\n")

print(f"SLURM script created: {slurm_script}")

submit_cmd = f"sbatch -A cms -q standby -p cpu {slurm_script}"
print("Submitting jobs...")
os.system(submit_cmd)

print("All jobs submitted.")
