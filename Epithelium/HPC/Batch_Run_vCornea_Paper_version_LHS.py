from pathlib import Path
import os
import glob
import shutil
import subprocess
import datetime
import sys
import time
import pandas as pd


base_path = '/u/jvanin/'

#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
#                                                               |
#---- Time Scales ---
HOURtoMCS   = 10
DAYtoMCS    = 240
WEEKtoMCS   = 1680
MONTHtoMCS  = 7300
YEARtoMCS   = 87600

replicates = 3

var_dict = {
#---STEM---
    'InitSTEM_LambdaSurface'    : [2.0],
    'InitSTEM_TargetSurface'    : [18.0],

    'InitSTEM_LambdaVolume'     : [3.0],
    'InitSTEM_TargetVolume'     : [25.0],

    'DensitySTEM_HalfMaxValue'  : [120.0],  # 
    'EGF_STEM_HalfMaxValue'     : [0.8], #

    'InitSTEM_LambdaChemo'      : [100.0],
# Growth Scalars 
    'EGF_GrowthScalar_STEM'     : [1.0],
    'DensityGrowthScalar_STEM'  : [1.0],
#---BASAL---
    'InitBASAL_LambdaSurface'   : [2.0],
    'InitBASAL_TargetSurface'   : [20.0],

    'InitBASAL_LambdaVolume'    : [3.0],
    'InitBASAL_TargetVolume'    : [25.0],

    'InitBASAL_LambdaChemo'     : [100.0],
    'InitBASAL_Division'        : [1000.0],

    'DensityBASAL_HalfMaxValue' : [90.0], #[55.0, 60.0, 65.0][50.0, 55.0, 60.0, 65.0, 70.0]
    'EGF_BASAL_HalfMaxValue'    : [1.62], #[3.75, 3.5, 3.25]4.0[3.5, 4.0, 4.5]
# Growth Scalars 
    'EGF_GrowthScalar_BASAL'    : [1.0],
    'DensityGrowthScalar_BASAL' : [1.0],
#---WING---
    'InitWING_LambdaSurface'    : [5.0],
    'InitWING_TargetSurface'    : [23],

    'InitWING_LambdaVolume'     : [20.0],
    'InitWING_TargetVolume'     : [25.0],

    'InitWING_EGFLambdaChemo'   : [1000.0],
#---SUPERFICIAL---
    'InitSUPER_LambdaSurface'   : [5.0],
    'InitSUPER_TargetSurface'   : [25.0],

    'InitSUPER_LambdaVolume'    : [2.0],
    'InitSUPER_TargetVolume'    : [25.0],

    'EGF_SUPERDiffCoef'         : [0.089], # 1/2048 barrier | 1/1024 barrier | 1/512 barrier  

    'SloughProbability'         : [0.0],
# Death Scalars
    'DeathTimeScalar'           : [1],
    'DeathVolumeScalar'         : [1],
    'SloughScalar'              : [1],
# FIELDS
    'MovementBiasScreteAmount'  : [1],
    'MovementBiasUptake'        : [1],

    'EGF_ScreteAmount'          : [2],

    'EGF_FieldUptakeBASAL'      : [0.1],
    'EGF_FieldUptakeSTEM'       : [0.1],
    'EGF_FieldUptakeSuper'      : [0.1],
    'EGF_GlobalDecay'           : [0.05331901388], # lower, mean, upper
# WOUND
    'InjuryType'                : ["'A'"],
    'IsInjury'                  : [False],
    'InjuryTime'                : [500],
#---INJURY AREA---
    'InjuryX_Center'            : [150],
    'InjuryY_Center'            : [60],
    'InjuryRadius'              : [25],
# DEBUGGING
    'GrowthControl'             : [True],
    'MitosisControl'            : [True],
    'DeathControl'              : [True],
    'DifferentiationControl'    : [True],
#   PLOTS
    'CellCount'                 : [True],
    'XHypTracker'               : [False],
    'SloughTracker'             : [False],
    'PressureTracker'           : [False],
    'PressurePlot'              : [False],
    'VolumeTracker'             : [False],
    'EGF_SeenByCell'            : [True],
    'CenterBiasPlot'            : [False],
    'CenterBias'                : [False],
    'DivisionTracker'           : [False],
    'GrowthPlot'                : [True],
    'ThicknessPlot'             : [True],
    'VolumeSurfaceDetail'       : [False],
    'MitosisPlot'               : [False],

# TIME OF SIMULATION
    'SimTime'                   : [MONTHtoMCS*2],
}

#                                                               |
#||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

def estimate_average_time(sim_time, num_combinations, num_rep, mcs=7301, seconds=1292):
    # 14600 MCS took 2736 seconds.
    # Estimate average time for one combination
    avg_time_per_combination = (sim_time / mcs) * seconds
    # Calculate total time for all combinations and reps
    total_time = avg_time_per_combination * num_combinations * num_rep
    return total_time

def estimate_time_at_capacity(total_time, capacity_percentage, jobs_at_full_capacity=103):
    adjusted_jobs = (capacity_percentage / 100) * jobs_at_full_capacity
    return total_time / adjusted_jobs if adjusted_jobs != 0 else float('inf')

def calculate_total_combinations(dict):
    total_combinations = 1
    for key, values in dict.items():
        total_combinations *= len(values)
    return total_combinations


def generate_dicts_from_csv(csv_file_path):
    """
    Generate a list of dictionaries with parameter-value mappings from a CSV file.
    
    Args:
    - csv_file_path (str): Path to the CSV file containing the parameter combinations.
    
    Returns:
    - List of dictionaries where each dictionary represents a unique combination of parameter-value pairs.
    """
    # Read the CSV file
    df = pd.read_csv(csv_file_path)
    
    # Convert the DataFrame to a list of dictionaries
    combinations_dicts = df.to_dict(orient='records')
    
    return combinations_dicts

def _original_generate_lists(data, keys=None, values=None):
    """Original recursive logic to generate combinations from the provided dictionary."""
    if keys is None:
        keys = list(data.keys())
    if values is None:
        values = []

    if keys:
        current_key = keys.pop(0)
        current_values = data[current_key]
        new_values = []
        if values:
            for value in values:
                for current_value in current_values:
                    if isinstance(value, list):
                        new_values.append(value + [current_value])
                    else:
                        new_values.append([value, current_value])
        else:
            new_values = [[value] for value in current_values]
        return _original_generate_lists(data, keys, new_values)
    else:
        return values
    
def generate_lists(data, csv_file_path=None):
    """
    Generate a list of unique combinations from a dictionary or CSV file.
    
    Args:
    - data (dict): Dictionary containing parameter-value pairs.
    - csv_file_path (str, optional): Path to the CSV file containing parameter combinations.
    
    Returns:
    - List of lists where each inner list represents a unique combination of parameter values.
    """
    
    if csv_file_path:
        # Get the combinations as dictionaries
        combinations_dicts = generate_dicts_from_csv(csv_file_path)
        
        # Convert dictionaries to lists while preserving all keys from the original var_dict
        combinations = []
        for comb_dict in combinations_dicts:
            comb = []
            for key in data.keys():
                if key in comb_dict:
                    comb.append(comb_dict[key])
                else:
                    comb.append(data[key][0])
            combinations.append(comb)
    else:
        combinations = _original_generate_lists(data)
    
    return combinations

def generate_sbatch_string(cc3d_file, node=8, job_name='vCornea', out_file='simulation_out',
                            err_file='err_file'):
    sbatch_string = f"""#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks-per-node={node}    
#SBATCH -J {job_name}
#SBATCH -o {out_file}_%j.txt
#SBATCH -e {err_file}_%j.err
#SBATCH --mail-user=jvanin@iu.edu

export DISPLAY=:1
echo $DISPLAY
srun \
{base_path}.conda/envs/cc3d_441_310/bin/cc3d_runScript.sh -i {cc3d_file}

wait


OUTPUT_FILE="{out_file}_${{SLURM_JOBID}}.txt"
ERROR_FILE="{err_file}_${{SLURM_JOBID}}.err"
mail -s "SLURM Job ${{SLURM_JOBID}} Results" jvanin@iu.edu < "$OUTPUT_FILE"
mail -s "SLURM Job ${{SLURM_JOBID}} Errors" jvanin@iu.edu < "$ERROR_FILE"     
    """
    return sbatch_string

# Model inputs to be changed must be in a dictionary format with the same 
# name as the created BatchTest.py file

mult_dict = var_dict

# number of replicates
num_rep = replicates

# Model output frequency
model_out_freq = 1

# Output frequency of simulation data per simulation replica
out_freq = 120 # VTK frequence in MCS



# Define base output directory
base_output_folder = r'/u/jvanin/Output_vCornea_Paper_version3'

# Generate a timestamp
timestamp = datetime.datetime.now().strftime("%Y%m%d-%H%M%S")

# Append timestamp to base directory
sweep_output_folder = os.path.join(base_output_folder, f'Output_{timestamp}')

# Ensure directory exists
if not os.path.isdir(sweep_output_folder):
    Path(sweep_output_folder).mkdir(parents=True, exist_ok=True)

# TODO - Change this line with your CC3D simulation folder path 
simulation_folder = Path(r'/u/jvanin/Projects/vCornea_v_PaperHPC_version_3') 

files_toCopy = glob.glob(os.path.join(simulation_folder, "*.cc3d"))
files_toCopy.extend(glob.glob(os.path.join(simulation_folder, "*.piff")))
files_toCopy.extend(glob.glob(os.path.join(simulation_folder,'Simulation', "*.py")))
files_toCopy.extend(glob.glob(os.path.join(simulation_folder,'Simulation', "*.xml")))
files_toCopy.extend(glob.glob(os.path.join(simulation_folder,'screenshot_data', "*.json")))

output = generate_lists(var_dict, csv_file_path='/u/jvanin/lhs_samples_extended.csv')

# Calculate the total combinations
num_combinations = len(output)

total_time_seconds = estimate_average_time(var_dict['SimTime'][0], num_combinations, replicates)

days, remainder = divmod(total_time_seconds, 86400)  # 86400 seconds in a day
hours, remainder = divmod(remainder, 3600)
minutes, seconds = divmod(remainder, 60)
print('')
print(f"Total combinations possible: {num_combinations}")
print(f"Estimated time for completion if run sequentially: {days} days, {hours} hours, {minutes} minutes, and {seconds} seconds.")

capacities = [100, 75, 50, 25]
for cap in capacities:
    estimated_time = estimate_time_at_capacity(total_time_seconds, cap)
    e_days, e_remainder = divmod(estimated_time, 86400)
    e_hours, e_remainder = divmod(e_remainder, 3600)
    e_minutes, e_seconds = divmod(e_remainder, 60)
    print(f"At {cap}% capacity (approx. {cap/100*103:.0f} jobs in parallel), estimated time: {e_days} days, {e_hours} hours, {e_minutes} minutes, and {e_seconds:.2f} seconds.")

# Ask for user confirmation
confirmation = input("Do you want to proceed? (yes/no): ")

if confirmation.lower() != "yes":
    sys.exit("Operation aborted by user.")

# Create files and parameters 
for idx, comb in enumerate(output):        

        for ii in range(num_rep):

            rep_folder = os.path.join(sweep_output_folder, f'Comb{idx}',str(ii))
            Path(rep_folder).mkdir(parents=True, exist_ok=True)

            rep_folder_simulation = os.path.join(rep_folder, "Simulation")
            Path(rep_folder_simulation).mkdir(parents=True, exist_ok=True)
            

            for file in files_toCopy:

                if "Simulation" in file:
                    shutil.copy2(file,rep_folder_simulation)

                elif "screenshot_data" in file:
                    rep_folder_screenshot_data = os.path.join(rep_folder, "screenshot_data")
                    Path(rep_folder_screenshot_data).mkdir(parents=True, exist_ok=True)

                    shutil.copy2(file, rep_folder_screenshot_data)

                else:
                    shutil.copy2(file, rep_folder)

            cc3d_file = glob.glob(os.path.join(rep_folder, "*.cc3d"))[0]

            with open(os.path.join(rep_folder, "BatchCall.sh"),"w") as f:

                f.write(generate_sbatch_string(cc3d_file, node=1,
                                                job_name='simulation',
                                                  out_file='simulation_out',
                                                    err_file='err_file'))

            with open(os.path.join(rep_folder_simulation, "Parameters.py"),"w") as f:
                f.write(f'rep_folder = "{rep_folder}"\n')
                for val_idx, key in enumerate(var_dict.keys()):
                   f.write(f'{key}={comb[val_idx]}\n')

#  Call batch	
files_toCall = glob.glob(os.path.join(sweep_output_folder,"**", "BatchCall.sh"),recursive=True)
print(files_toCall)
processes = []
for files in files_toCall:
    os.chmod(files, 0o777)
    directory = os.path.dirname(files)
    p = subprocess.Popen(['sbatch', files], cwd=directory)
    processes.append(p)
    
