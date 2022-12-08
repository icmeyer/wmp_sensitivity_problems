#!/bin/bash
#SBATCH --job-name=4mwd_bwr  # Job name
#SBATCH --ntasks=12           # Number of tasks
# #SBATCH --mem=32gb          # Job memory request
#SBATCH --time=48:00:00       # Time limit hrs:min:sec

# Load modules
module purge
module load gcc-11.2.0-gcc-4.8.5-ytyiq4m
source /opt/intel/oneapi/setvars.sh
source ~/.bashrc

# Run openmc
python input.py
