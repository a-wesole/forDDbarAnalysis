import os
import subprocess


#step 1 compile code
try:
    compile_command = [
        "g++",
        # "Make2dMassHistos_M1M2_data_10SlicesEach.C",
        "Make2dMassHistos_M1M2_data_10SlicesEach_DDupdates.C",
        "../../src/BDTHandler.cc",
        "../../src/phiTools.cc",
        "-I../../include",
        "$(root-config --cflags --libs)",
        "-Wall",
        "-O2",
        "-o",
        "./M1M2_data_10SlicesEach.exe"
    ]
    subprocess.run(" ".join(compile_command), shell=True, check=True)
    print("Compilation successful!")
    # exit(1)
except subprocess.CalledProcessError:
    print("Compilation failed. Please check your source code.")
    exit(1)

# Step 2: Run the executable or proceed with other steps
print("Running subsequent operations...")


class jobs_submission():
    def __init__(self):
        self.username = 'awesole'
    # Build job submission files 
    def build_job_files(self, job_name):
        dir_job_file = os.path.join(self.base_dir, 'job_files/')
        dir_job_out = os.path.join(self.base_dir, 'log/')
        dir_job_err = os.path.join(self.base_dir, 'err/')
        dir_job_output = self.base_dir + '/ROOT/'

        for d in [dir_job_file, dir_job_out, dir_job_err, dir_job_output]:
            if not os.path.isdir(d):
                os.makedirs(d)

        subfile = os.path.join(dir_job_file, job_name + '.sub')

        # SLURM job submission file
        with open(subfile, "w") as fsub:
            fsub.write('#!/bin/bash\n')
            fsub.write(f'#SBATCH --nodes=1\n')
            fsub.write(f'#SBATCH --mem=16G\n')
            fsub.write(f'#SBATCH -o {dir_job_out}{job_name}.out\n')
            fsub.write(f'#SBATCH -e {dir_job_err}{job_name}.err\n')
            fsub.write(f'#SBATCH --job-name="{job_name}"\n')
            fsub.write('#SBATCH --no-requeue\n')
            fsub.write('#SBATCH --export=ALL\n')
            fsub.write('export X509_USER_PROXY=/home/awesole/myproxy\n')
            fsub.write(f'./M1M2_data_10SlicesEach.exe {self.istart} {self.iend} \n')

        return subfile

    # Submit jobs
    def submit_jobs(self, base_dir, njobs, nfile_per_job, dry_run=False):
        self.base_dir = base_dir

        for i in range(njobs):
            self.istart = i*nfile_per_job
            self.iend = self.istart + nfile_per_job
            job_name = f"M1M2_f_{self.istart}_{self.iend}"

            subfile = self.build_job_files(job_name )

            if dry_run:
                print(f"Generated submission script: {subfile}")
            else:
                # os.system(f"sbatch -A physics -q standby -p cpu --time=01:00:00 {subfile}")
                os.system(f"sbatch -A cms -q standby -p cpu --time=01:00:00 {subfile}")



if __name__ == "__main__":
    base_dir = "/scratch/negishi/awesole/DDbar_Jan2026"
    flist_name = "/home/awesole/forDDbarAnalysis/forTemplatesMethod/dataFiles/Run_2023_19Sept.list"

    sub = jobs_submission()
    with open(flist_name) as f:
        flist = list(f)

    nfile_per_job = 1200
    njobs = int(len(flist) / nfile_per_job) + 1
    print(f"Number of files per job: {nfile_per_job}")
    print(f"Number of files in flist: {len(flist)}")
    print(f"Number of jobs: {njobs}")


    sub.submit_jobs(base_dir, njobs, nfile_per_job, dry_run=False)
