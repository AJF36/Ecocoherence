#!/bin/bash
#SBATCH --cpus-per-task=20 # Number of CPUs
#SBATCH --mem=16GB # RAM
#SBATCH -t 7-00:00 # Maximum running time (D-HH:MM)
#SBATCH -o my_job_%j.out # Standard output file, %j is the <jobID>
#SBATCH -e my_job_%j.err # Standard error file, %j is the <jobID>
#SBATCH --mail-user adrian.jimenez@cnb.csic.es # E-mail address for notifications
#SBATCH --mail-type=FAIL,END # Notify successful or failed completion of job

# Get the directory where the script is located
#script_dir=$(dirname "$(realpath "$0")")

# Initialize Conda
source /home/ajimenez/miniconda3/etc/profile.d/conda.sh

# Activate the Conda environment
conda activate misosoup_cluster
# Check if the environment was activated successfully
if [[ $? -ne 0 ]]; then
    echo "Failed to activate Conda environment: smetana_gurobi"
    exit 1
fi

# Run misosoup
 
misosoup ./*xml --output misosoup_Agarose.yml --media media.yaml 
