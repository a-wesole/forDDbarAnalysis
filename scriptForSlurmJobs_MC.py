import os
import subprocess

class jobs_submission():
    def __init__(self):
        self.username = 'awesole'

    # Build job submission files 
    def build_job_files(self, job_name, outfile, ptmin, ptmax, centmin, centmax):
        dir_job_file = self.base_dir + '/job_files/'
        dir_job_out = self.base_dir + '/log/'
        dir_job_err = self.base_dir + '/err/'
        dir_job_output = self.base_dir + '/ROOT/'

        if not os.path.isdir(dir_job_file):
            os.makedirs(dir_job_file)
        if not os.path.isdir(dir_job_out):
            os.makedirs(dir_job_out)
        if not os.path.isdir(dir_job_err):
            os.makedirs(dir_job_err)
        if not os.path.isdir(dir_job_output):
            os.makedirs(dir_job_output)

        subfile = dir_job_file + job_name + '.sub'

        # SLURM job submission file
        with open(subfile, "w") as fsub:
            fsub.write('#!/bin/bash\n')
            fsub.write(f'#SBATCH --nodes=1\n')
            fsub.write(f'#SBATCH -o {dir_job_out}{job_name}.out\n')
            fsub.write(f'#SBATCH -e {dir_job_err}{job_name}.err\n')
            fsub.write(f'#SBATCH --job-name="{job_name}"\n')
            fsub.write('#SBATCH --no-requeue\n')
            fsub.write('#SBATCH --export=ALL\n')
            fsub.write(f'./submitforMCloop.exe {outfile} {ptmin} {ptmax} {centmin} {centmax}\n')

        return subfile

    # Submit jobs
    def submit_jobs(self, base_dir, pt_ranges, cent_ranges, dry_run=False):
        self.base_dir = base_dir

        dir_job_output = os.path.join(base_dir, 'ROOT/')


        for ptmin, ptmax in pt_ranges:
            for centmin, centmax in cent_ranges:
                job_name = f"pt_{ptmin}_{ptmax}_cent_{centmin}_{centmax}"
                outfile = os.path.join(dir_job_output, f"mcOut_{job_name}.root")

                subfile = self.build_job_files(job_name, outfile, ptmin, ptmax, centmin, centmax)

                if dry_run:
                    print(f"Generated submission script: {subfile}")
                else:
                    os.system(f"sbatch --time=04:00:00 -A standby {subfile}")

if __name__ == "__main__":
    base_dir = "/scratch/bell/awesole/mcFilesForTraining"

    sub = jobs_submission()

    pt_ranges = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 8), (8, 10), (10, 15), (15, 20), (20, 30), (30, 50), (50, 100), (100, 1000)]
    cent_ranges = [(0, 10), (10, 30),(30, 50),(50, 90)]

    # Set dry_run to True to perform a dry run
    sub.submit_jobs(base_dir, pt_ranges, cent_ranges, dry_run=False)

