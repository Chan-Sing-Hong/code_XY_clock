#!/bin/bash

#SBATCH --partition=edr1-al9_short        # Specify the partition
#SBATCH --output=tmp/#job_name#.out       # Standard output file
#SBATCH --job-name=#job_name#             # Set job name
#SBATCH --cpus-per-task=2                 # Number of CPU cores per task

# Set environment variables
export PATH="~/anaconda3/envs/cytnx_v1.0.1/bin:${PATH}"
mkdir -p tmp

# Log the start time
echo "Job started on $(date)"

# Run the Python script
python code_tmp_#job_name#.py

# Check if the script ran successfully
if [ $? -eq 0 ]; then
    echo "Job completed successfully on $(date)"
    rm code_tmp_#job_name#.py
    rm tmp/job_tmp_#job_name#.sh
else
    echo "Job failed on $(date)"
    exit 1
fi

echo -e "\nJob done.\n"

