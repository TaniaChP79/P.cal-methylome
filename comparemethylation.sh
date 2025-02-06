#!/bin/bash
 
#SBATCH --nodes=1                   # the number of nodes you want to reserve
#SBATCH --ntasks-per-node=1         # the number of tasks/processes per node
#SBATCH --cpus-per-task=36          # the number cpus per task
#SBATCH --partition=normal          # on which partition to submit the job
#SBATCH --time=24:00:00             # the max wallclock time (time limit your job will run)
#SBATCH --mem=90G                   # request RAM
#SBATCH --job-name=compmethyl         # the name of your job
#SBATCH --mail-type=ALL             # receive an email when your job starts, finishes normally or is aborted
#SBATCH --mail-user=tchavarr@uni-muenster.de # your mail address
 
# USAGE: sbatch methylationfreq.sh

# LOAD MODULES HERE IF REQUIRED
module load palma/2022b
module load GCCcore/12.2.0
module load Python/3.10.8 

# START THE APPLICATION

python compare_methylation.py /scratch/tmp/tchavarr/bisulfite.ENCFF835NTC.example.tsv /scratch/tmp/tchavarr/methylation_frequency.tsv > /scratch/tmp/tchavarr/bisulfite_vs_nanopolish.tsv

