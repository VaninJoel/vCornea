#!/bin/bash

#(CHANCE THIS TO YOUR PYTHON SCRIPT)
base_folder="/u/jvanin/Output_vCornea_Paper_version2/Output_20230920-164340"
num_jobs=103  # adjust as necessary Tatooine has 103 cores

# Find all Comb# folders
folders=( $(ls $base_folder | grep 'Comb') )

# Count the number of folders and calculate the number of folders each job should process
total_folders=${#folders[@]}
folders_per_job=$(( (total_folders + num_jobs - 1) / num_jobs ))  # ceiling division

# Submit jobs
for ((i=0; i<$num_jobs; i++)); do
    start_idx=$(( i * folders_per_job ))
    end_idx=$(( start_idx + folders_per_job - 1 ))

    # Extract the specific folders this job should process
    job_folders="${folders[@]:$start_idx:$folders_per_job}"

    sbatch --job-name=process_comb_data_$i <<EOL
#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH -J PlottingHPC
#SBATCH --mail-user=jvanin@iu.edu
#SBATCH --output=comb_data_out_$i.txt
#SBATCH --error=comb_data_err_$i.txt

OUTPUT_FILE="comb_data_out_${SLURM_JOBID}.txt"
ERROR_FILE="comb_data_err_${SLURM_JOBID}.err"

# Run the Python script (CHANGE THIS TO YOUR PYTHON SCRIPT)
python Plotting_Classifying_HPC_Parallel_v2.py $job_folders

# After Python script runs, send email
TEMP_FILE=\$(mktemp)

# Append output and errors to the temporary file
echo "----- OUTPUT -----" >> \$TEMP_FILE
cat "\$OUTPUT_FILE" >> \$TEMP_FILE
echo -e "\n\n----- ERRORS -----" >> \$TEMP_FILE
cat "\$ERROR_FILE" >> \$TEMP_FILE

# Send the merged contents in one email
mail -s "SLURM Job \${SLURM_JOBID} Results and Errors" jvanin@iu.edu < "\$TEMP_FILE"

# Remove the temporary file
rm -f \$TEMP_FILE
EOL

done