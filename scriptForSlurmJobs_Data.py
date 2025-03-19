import os
import subprocess

class jobs_submission():
    def __init__(self, flist_name):
        self.flist_name = flist_name
        self.max_njobs_per_round = 1500
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
            fsub.write('export X509_USER_PROXY=/home/awesole/myproxy\n')
            fsub.write(f'./updated_submitfordataloop.exe {ptmin} {ptmax} {centmin} {centmax} {self.istart} {self.iend} \n')
            fsub.close()

        return subfile

    # Submit jobs
    def submit_jobs(self, base_dir, pt_ranges, cent_ranges,  njobs, nfile_per_job, dry_run=False):
        self.base_dir = base_dir

        dir_job_output = os.path.join(base_dir, 'ROOT/')


        for ptmin, ptmax in pt_ranges:
            for centmin, centmax in cent_ranges:
                count = 0
                for i in range(njobs):
                    if count >= self.max_njobs_per_round: # hold until last round jobs are done
                        self.check_running_jobs()
                        count=0  # restart the counter for the next round
                    else:
                        count += 1

                    print('job: ', i)
                    self.istart = i*nfile_per_job
                    self.iend = self.istart + nfile_per_job
                    job_name = f"pt_{ptmin}_{ptmax}_cent_{centmin}_{centmax}_file_{self.istart}_{self.iend}"
                    outfile = os.path.join(dir_job_output, f"dataOut_{job_name}.root")
                    subfile = self.build_job_files(job_name, outfile, ptmin, ptmax, centmin, centmax)
                    if dry_run:
                        print(f"Generated submission script: {subfile}")
                    else:
                        os.system(f"sbatch --time=04:00:00 -A standby {subfile}")
                    
if __name__ == "__main__":
    base_dir = "/scratch/negishi/awesole/sidebands_complete_original_production"
    flist_name = "/home/awesole/forDDbarAnalysis/dataFiles/dataFilesTTrees.txt"

    sub = jobs_submission(flist_name)

    with open(flist_name) as f:
        flist = list(f)

    nfile_per_job = 1003
    njobs = int(len(flist) / nfile_per_job) 
    print(f"Number of files in flist: {len(flist)}")
    print(f"Number of jobs: {njobs}")

    pt_ranges = [(1, 2), (2, 3), (3, 4), (4, 5), (5, 6), (6, 8), (8, 10), (10, 15), (15, 20), (20, 30), (30, 50), (50, 100), (100, 1000)]
    cent_ranges = [(0, 10), (10, 30),(30, 50),(50, 90)]
    #pt_ranges = [(1, 2)]
    #cent_ranges = [(0, 10)]

    # Set dry_run to True to perform a dry run
    sub.submit_jobs(base_dir, pt_ranges, cent_ranges, njobs, nfile_per_job, dry_run=False) 

