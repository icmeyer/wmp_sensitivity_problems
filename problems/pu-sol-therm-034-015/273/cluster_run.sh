#!/bin/bash
#SBATCH --job-name=pu-sol  # Job name
#SBATCH --ntasks=12           # Number of tasks
# #SBATCH --mem=32gb          # Job memory request
#SBATCH -t 4-00:00            # Time limit: (D-HH:MM)

# Load modules
module purge
module load gcc-11.2.0-gcc-4.8.5-ytyiq4m
source /opt/intel/oneapi/setvars.sh
source ~/.bashrc

# Run openmc
python input.py
