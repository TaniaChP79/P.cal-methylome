#!/bin/bash  
 
#SBATCH --partition=gpu2080
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --gres=gpu:8
#SBATCH --job-name=guppy_basecalling
#SBATCH --time=48:00:00  # Adjust the maximum runtime as needed
#SBATCH --output=guppy_basecalling.out
#SBATCH --error=guppy_basecalling.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=tchavarria@uni-muenster.de



# Load required modules

module load palma/2022a
module load ont-guppy/6.4.8-CUDA-11.7.0

# Run Guppy basecalling

guppy_basecaller --disable_pings \
  -i /scratch/tmp/tchavarr/PogoCA/fast5/ \
  -s /scratch/tmp/tchavarr/PogoCA/basecalling/ \
  -c dna_r9.4.1_450bps_fast.cfg \
  --min_qscore 7 \
  --recursive -x 'cuda:0' \
  --num_callers 16 \
  --cpu_threads_per_caller 1 \
  --gpu_runners_per_device 64 \
  --chunks_per_runner 512 \
  --chunk_size 1000 \
  --compress_fastq


