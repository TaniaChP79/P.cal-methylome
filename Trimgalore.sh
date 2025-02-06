#!/bin/bash
 
#SBATCH --nodes=1                   # the number of nodes you want to reserve
#SBATCH --ntasks-per-node=1         # the number of tasks/processes per node
#SBATCH --cpus-per-task=36          # the number cpus per task
#SBATCH --partition=normal          # on which partition to submit the job
#SBATCH --time=24:00:00             # the max wallclock time (time limit your job will run)
 
#SBATCH --job-name=trimgalore         # the name of your job
#SBATCH --mail-type=ALL             # receive an email when your job starts, finishes normally or is aborted
#SBATCH --mail-user=tchavarr@uni-muenster.de # your mail address
 
# LOAD MODULES HERE IF REQUIRED

module load palma/2019a 
module load GCCcore/8.2.0 
module load TrimGalore/0.6.6	

# START THE APPLICATION

trim_galore /scratch/tmp/tchavarr/WGBSdata/raw/*/*_1.fq /scratch/tmp/tchavarr/WGBSdata/raw/*/_2.fq --paired --clip_R1 2 --clip_R2 2 --three_prime_clip_R1 2 --three_prime_clip_R2 2 --illumina -o /scratch/tmp/tchavarr/WGBSdata/clean/ 
