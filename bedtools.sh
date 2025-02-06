#!/bin/bash
 
#SBATCH --nodes=1                   # the number of nodes you want to reserve
#SBATCH --ntasks-per-node=1         # the number of tasks/processes per node
#SBATCH --cpus-per-task=36          # the number cpus per task
#SBATCH --partition=normal          # on which partition to submit the job
#SBATCH --time=24:00:00             # the max wallclock time (time limit your job will run)
 
#SBATCH --job-name=bedtools2         # the name of your job
#SBATCH --mail-type=ALL             # receive an email when your job starts, finishes normally or is aborted
#SBATCH --mail-user=tchavarr@uni-muenster.de # your mail address
 
# LOAD MODULES HERE IF REQUIRED

module load palma/2021a  
module load GCC/10.3.0
module load BEDTools/2.30.0

# START THE APPLICATION
bedtools intersect -a /scratch/tmp/tchavarr/pogohaplo/log_methylation_frequency_H1.bed -b /scratch/tmp/tchavarr/P.californicus_sorted.bed -wa -wb > /scratch/tmp/tchavarr/pogohaplo/methylation_exons_introns_H1.bed
