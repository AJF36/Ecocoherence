#!/bin/bash
# Script to submit all 25 jobs to the computing cluster

# Adjust these variables according to your cluster system
# This example assumes SLURM, but you can modify for PBS, SGE, etc.

# Loop through the 25 jobs
for job_id in {1..4}; do
    echo "Submitting job $job_id of 25"
    
    # Create a job submission script
    cat > job_${job_id}.sh << EOF
#!/bin/bash
#SBATCH --cpus-per-task=5 # Number of CPUs
#SBATCH --mem=8GB # RAM
#SBATCH -t 0-160:00 # Maximum running time (D-HH:MM)
#SBATCH -o my_job_%j.out # Standard output file, %j is the <jobID>
#SBATCH -e my_job_%j.err # Standard error file, %j is the <jobID>
#SBATCH --mail-user adrian.jimenez@cnb.csic.es # E-mail address for notifications
#SBATCH --mail-type=FAIL,END # Notify successful or failed completion of job



#Activate conda
source /home/ajimenez/miniconda3/etc/profile.d/conda.sh
# Activate the Conda environment
conda activate r_scripts


# Run the R script with job_id as argument
Rscript creation_union_intersection_df.R ${job_id}
EOF

    # Submit the job
    sbatch job_${job_id}.sh
done

echo "All 25 jobs submitted successfully."
