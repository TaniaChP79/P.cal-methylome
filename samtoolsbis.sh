#!/bin/bash
 
#SBATCH --nodes=1                   # the number of nodes you want to reserve
#SBATCH --ntasks-per-node=1         # the number of tasks/processes per node
#SBATCH --cpus-per-task=36          # the number cpus per task
#SBATCH --partition=normal          # on which partition to submit the job
#SBATCH --time=24:00:00             # the max wallclock time (time limit your job will run)
 
#SBATCH --job-name=samtoolsbis         # the name of your job
#SBATCH --mail-type=ALL             # receive an email when your job starts, finishes normally or is aborted
#SBATCH --mail-user=tchavarr@uni-muenster.de # your mail address
 
# LOAD MODULES HERE IF REQUIRED

module load palma/2022a
module load GCC/11.3.0
module load SAMtools/1.16.1


# START THE APPLICATION

samtools sort -n /scratch/tmp/tchavarr/WGBSdata/alignment/S40_R1_val_1_bismark_bt2_pe.bam -o /scratch/tmp/tchavarr/WGBSdata/alignment/S40_R1_sorted.bam
samtools fixmate -r -m /scratch/tmp/tchavarr/WGBSdata/alignment/S40_R1_sorted.bam /scratch/tmp/tchavarr/WGBSdata/alignment/S40_R1_fix.bam
samtools sort /scratch/tmp/tchavarr/WGBSdata/alignment/S40_R1_fix.bam -o /scratch/tmp/tchavarr/WGBSdata/alignment/S40_R1_fix_sort.bam
samtools markdup -r /scratch/tmp/tchavarr/WGBSdata/alignment/S40_R1_fix_sort.bam /scratch/tmp/tchavarr/WGBSdata/alignment/S40_R1_dedup.bam
samtools sort -n /scratch/tmp/tchavarr/WGBSdata/alignment/S40_R1_dedup.bam -o /scratch/tmp/tchavarr/WGBSdata/alignment2/S40_R1_dedup_sort.bam
 
