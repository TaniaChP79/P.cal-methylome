#!/bin/bash
 
#SBATCH --nodes=1                   # the number of nodes you want to reserve
#SBATCH --ntasks-per-node=1         # the number of tasks/processes per node
#SBATCH --cpus-per-task=36          # the number cpus per task
#SBATCH --partition=normal          # on which partition to submit the job
#SBATCH --time=24:00:00             # the max wallclock time (time limit your job will run)
#SBATCH --mem=90G                   # request RAM
#SBATCH --job-name=samtools         # the name of your job
#SBATCH --mail-type=ALL             # receive an email when your job starts, finishes normally or is aborted
#SBATCH --mail-user=tchavarr@uni-muenster.de # your mail address
 
# USAGE: sbatch samtools.sh

# LOAD MODULES HERE IF REQUIRED
# LOAD MODULES HERE IF REQUIRED

module load palma/2022a
module load GCC/11.3.0
module load SAMtools/1.16.1


# START THE APPLICATION

samtools sort -T tmp -o /scratch/tmp/tchavarr/pogoexp/output_C1.sorted.bam 
samtools index /scratch/tmp/tchavarr/pogoexp/output_C1.sorted.bam
