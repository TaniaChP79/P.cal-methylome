#!/bin/bash
 
#SBATCH --nodes=1                   # the number of nodes you want to reserve
#SBATCH --ntasks-per-node=1         # the number of tasks/processes per node
#SBATCH --cpus-per-task=8          # the number cpus per task
#SBATCH --partition=normal          # on which partition to submit the job
#SBATCH --time=24:00:00             # the max wallclock time (time limit your job will run)
#SBATCH --mem=92G                   # request RAM
#SBATCH --job-name=f5c        # the name of your job
#SBATCH --mail-type=ALL             # receive an email when your job starts, finishes normally or is aborted
#SBATCH --mail-user=tchavarr@uni-muenster.de # your mail address
 
# USAGE: sbatch nanopolishcha.sh

# LOAD MODULES HERE IF REQUIRED
# LOAD MODULES HERE IF REQUIRED
ml palma/2023a
ml gompi/2023a
ml slow5tools/1.1.0
ml f5c/1.4


# START THE APPLICATION

#convert fast5 files to slow5 files using 8 I/O processes
slow5tools f2s /scratch/tmp/tchavarr/PogoCNA3/Cpog_Pcq26b-4a-_20240111_Total/fast5_pass_total/ -d /scratch/tmp/tchavarr/PogoCNA3/blow5/ -p 8
#Merge all the slow5 files in to a single file using 8 threads
slow5tools merge /scratch/tmp/tchavarr/PogoCNA3/blow5/ -o /scratch/tmp/tchavarr/PogoCNA3/signals.blow5 -t8
#f5c index
f5c index -t 8  /scratch/tmp/tchavarr/PogoCNA3/PogoCNA3.fastq --slow5 /scratch/tmp/tchavarr/PogoCNA3/signals.blow5
#f5c methylation calling
f5c call-methylation -t 8 -r /scratch/tmp/tchavarr/PogoCNA3/PogoCNA3.fastq -g /scratch/tmp/tchavarr/Pcal3.1.clean.submitted.fa -b /scratch/tmp/tchavarr/PogoCNA3/outputCNA3.sorted.bam --slow5 /scratch/tmp/tchavarr/PogoCNA3/signals.blow5 > PogoCNA3meth.tsv



