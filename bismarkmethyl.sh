#!/bin/bash
 
#SBATCH --nodes=1                   # the number of nodes you want to reserve
#SBATCH --ntasks-per-node=1         # the number of tasks/processes per node
#SBATCH --cpus-per-task=36          # the number cpus per task
#SBATCH --partition=normal          # on which partition to submit the job
#SBATCH --time=24:00:00             # the max wallclock time (time limit your job will run)
 
#SBATCH --job-name=bismarkmethy         # the name of your job
#SBATCH --mail-type=ALL             # receive an email when your job starts, finishes normally or is aborted
#SBATCH --mail-user=tchavarr@uni-muenster.de # your mail address
 
# LOAD MODULES HERE IF REQUIRED

module load palma/2021b  
module load GCC/11.2.0
module load OpenMPI/4.1.1
module load Bismark/0.23.1




# START THE APPLICATION

bismark_methylation_extractor -p --no_overlap --cytosine_report --genome_folder /scratch/tmp/tchavarr/Pogogenome/ /scratch/tmp/tchavarr/WGBSdata/alignment2/S40_R1_dedup_sort.bam -o /scratch/tmp/tchavarr/WGBSdata/methylation_calling/
bismark_methylation_extractor -p --no_overlap --cytosine_report --genome_folder /scratch/tmp/tchavarr/Pogogenome/ /scratch/tmp/tchavarr/WGBSdata/alignment2/S41_R1_dedup_sort.bam -o /scratch/tmp/tchavarr/WGBSdata/methylation_calling/