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
 
# USAGE: sbatch f5c.sh

# LOAD MODULES HERE IF REQUIRED
# LOAD MODULES HERE IF REQUIRED
ml palma/2023a
ml gompi/2023a
ml slow5tools/1.1.0
ml f5c/1.4


# START THE APPLICATION

#f5c index
f5c index -t 8 --slow5 /scratch/tmp/tchavarr/PogoCNA3/signals.blow5 /scratch/tmp/tchavarr/PogoCNA3/PogoCNA3.fastq 
f5c index -d /scratch/tmp/tchavarr/PogoCNA3/Cpog_Pcq26b-4a-_20240111_Total/fast5_pass_total /scratch/tmp/tchavarr/PogoCNA3/PogoCNA3.fastq
#f5c methylation calling
f5c call-methylation -t 8 -b /scratch/tmp/tchavarr/PogoCNA3/outputCNA3.sorted.bam -g /scratch/tmp/tchavarr/Pcal3.1.clean.submitted.fa -r /scratch/tmp/tchavarr/PogoCNA3/PogoCNA3.fastq --slow5 /scratch/tmp/tchavarr/PogoCNA3/signals.blow5 > /scratch/tmp/tchavarr/PogoCNA3/CNA3meth.tsv
f5c call-methylation -b /scratch/tmp/tchavarr/PogoCNA3/outputCNA3.sorted.bam -g /scratch/tmp/tchavarr/Pcal3.1.clean.submitted.fa -r /scratch/tmp/tchavarr/PogoCNA3/PogoCNA3.fastq > /scratch/tmp/tchavarr/PogoCNA3/2meth.tsv # for fast5 --pore r10 required if R10.4.1
#methylation freq
f5c meth-freq -i /scratch/tmp/tchavarr/PogoCNA3/CNA3meth.tsv > /scratch/tmp/tchavarr/PogoCNA3/CNA3freq.tsv
f5c meth-freq -i /scratch/tmp/tchavarr/PogoCNA3/2meth.tsv > /scratch/tmp/tchavarr/PogoCNA3/2freq.tsv
