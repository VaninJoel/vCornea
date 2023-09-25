import os
import glob
import pandas as pd
import matplotlib.pyplot as plt
import argparse
from sklearn.linear_model import LinearRegression
import numpy as np
import shutil
from filelock import FileLock
import sys


class DataGatherer:
    def __init__(self, base_folder):       
            self.base_folder = base_folder

    def gather_csv_files_from_folders(self, folders_to_process):
        all_files = []
        for folder in folders_to_process:
            comb_folder = os.path.join(self.base_folder, folder)
            all_files.extend(glob.glob(os.path.join(comb_folder, "**", "*.csv"), recursive=True))
        return all_files
    
    def extract_params_from_path(self, path): 
        parts = path.split(os.path.sep)
        comb_part = [part for part in parts if part.startswith("Comb")]
        if not comb_part:
            print(f"Problematic path: {path}")
            raise ValueError("'Comb#' not found in path")
        comb_idx = comb_part[0][4:]
        replicate_idx = parts[parts.index(comb_part[0]) + 1]
        return comb_idx, replicate_idx
    
    def organize_files(self, csv_files):    
        organized = {}
        for file in csv_files:
            comb, _ = self.extract_params_from_path(file)
            base_name = os.path.basename(file)
            key = (comb, base_name)
            if key not in organized:
                organized[key] = []
            organized[key].append(file)
        return organized
   
class Classifier(DataGatherer):   
    CELL_COUNT_DIFF_THRESHOLD = {
        "Superficial": 7,
        "Wing": 55,
        "Basal": 12,
        "Stem": 7 }   
    LIMBUS_DIFF_THRESHOLD = 7
    PERIPHERY_DIFF_THRESHOLD = 7

    def __init__(self, organized_data):
        self.organized_data = organized_data

    @staticmethod
    def check_directionality(data, column, slope_threshold=0.01):
        X = np.arange(len(data)).reshape(-1, 1)
        y = data[column].values
        reg = LinearRegression().fit(X, y)
        slope = reg.coef_[0]        
        if slope > slope_threshold:
            return 'increasing', slope
        elif slope < -slope_threshold:
            return 'decreasing', slope
        else:
            return 'stable', slope
                
    def classify_data(self, files, data_type):
        flag_fail = False
        reasons = []                
        all_data = [pd.read_csv(f) for f in files]        
        combined_data = pd.concat(all_data).groupby(level=0).mean()

        if data_type == 'cell':
            columns = ['Superficial', 'Wing', 'Basal', 'Stem']
            initial_values = {col: combined_data[col].iloc[0] for col in columns}
            final_values = {col: combined_data[col].iloc[-1] for col in columns}
            diff_values = {col: final_values[col] - initial_values[col] for col in columns}
            trends = {col: self.check_directionality(combined_data, col) for col in columns}
            for col, (trend, slope) in trends.items():
                if trend in ['increasing', 'decreasing']:
                    reasons.append(f'{col} trend is {trend} ({slope})')
                    flag_fail = True
                else:
                    reasons.append(f'{col} trend is stable ({slope})')
            for col, diff in diff_values.items():
                if abs(diff) > self.CELL_COUNT_DIFF_THRESHOLD[col]:
                    reasons.append(f'Difference in {col} exceeded threshold {abs(diff)} > {self.CELL_COUNT_DIFF_THRESHOLD[col]}')
                    flag_fail = True
                else:
                    reasons.append(f'Difference in {col} is within threshold {abs(diff)} <= {self.CELL_COUNT_DIFF_THRESHOLD[col]}')
        elif data_type == 'thickness':
            initial_limbus = combined_data['Limbus Avg Thickness'].iloc[0]
            initial_periphery = combined_data['Periphery Avg Thickness'].iloc[0]
            final_limbus = combined_data['Limbus Avg Thickness'].iloc[-1]
            final_periphery = combined_data['Periphery Avg Thickness'].iloc[-1]
            diff_limbus = final_limbus - initial_limbus
            diff_periphery = final_periphery - initial_periphery
            limbus_trend, limbus_slope = self.check_directionality(combined_data, 'Limbus Avg Thickness')
            periphery_trend, periphery_slope = self.check_directionality(combined_data, 'Periphery Avg Thickness')
            if limbus_trend in ['increasing', 'decreasing']:
                reasons.append(f'Limbus trend is {limbus_trend} ({limbus_slope})')
                flag_fail = True
            else:
                reasons.append(f'Limbus trend is stable ({limbus_slope})')
            if periphery_trend in ['increasing', 'decreasing']:
                reasons.append(f'Periphery trend is {periphery_trend} ({periphery_slope})')
                flag_fail = True
            else:
                reasons.append(f'Periphery trend is stable ({periphery_slope})')
            if abs(diff_limbus) > self.LIMBUS_DIFF_THRESHOLD:
                reasons.append(f'Difference in Limbus exceeded threshold ({abs(diff_limbus)}) > {self.LIMBUS_DIFF_THRESHOLD}')
                flag_fail = True
            else:
                reasons.append(f'Difference in Limbus is within threshold ({abs(diff_limbus)}) <= {self.LIMBUS_DIFF_THRESHOLD}')
            if abs(diff_periphery) > self.PERIPHERY_DIFF_THRESHOLD:
                reasons.append(f'Difference in Periphery exceeded threshold ({abs(diff_periphery)}) > {self.PERIPHERY_DIFF_THRESHOLD}')
                flag_fail = True
            else:
                reasons.append(f'Difference in Periphery is within threshold ({abs(diff_periphery)}) <= {self.PERIPHERY_DIFF_THRESHOLD}')
        elif data_type == 'growth':
            pass
        else:
            raise ValueError("Invalid data type provided. Must be 'cell_count' or 'thickness'.", data_type) 
        
        comb_dir = os.path.dirname(os.path.dirname(os.path.dirname(files[0])))
        
        classification = 'Success' if not flag_fail else 'Failure'
        return classification, reasons, comb_dir
    
    def classify_comb_file_type(self, comb_id, file_type, files):             
        classification = self.classify_data(files, file_type)       
        return {
            'comb_id': comb_id,
            'file_type': file_type,
            'classification': classification
        }
    
    def classify_all(self):        
        results = []        
        for (comb_id, file_type), files in self.organized_data.items():            
            base_name_without_extension = os.path.splitext(file_type)[0]            
            actual_file_type = base_name_without_extension.split('_')[0]            
            comb_file_type_results = self.classify_comb_file_type(comb_id, actual_file_type, files)
            results.append(comb_file_type_results)                    
        return results
    
    def generate_combined_report(self, results):
        # Extract all unique comb_ids
        comb_ids = set(item['comb_id'] for item in results)

        # Prepare results dictionary
        classification_results = {}

        for cid in comb_ids:
            cell_data = [item for item in results if item['comb_id'] == cid and item['file_type'] == 'cell'][0]
            thickness_data = [item for item in results if item['comb_id'] == cid and item['file_type'] == 'thickness'][0]
            
            # Classify the comb_id based on the given criteria
            if 'Success' in cell_data['classification'][0] and 'Success' in thickness_data['classification'][0]:
                classification = 'Success'
            elif 'Success' in cell_data['classification'][0] or 'Success' in thickness_data['classification'][0]:
                classification = 'Partial'
            else:
                classification = 'Fail'
            
            # Build the report for the comb_id
            report = ["--- Classification Report for {} ---".format(cid)]
            report.append("\nCell Count Classification:")
            report.append(f"Result: {cell_data['classification'][0]}")
            for reason in cell_data['classification'][1]:
                report.append(f"- {reason}")
            report.append("\nThickness Classification:")
            report.append(f"Result: {thickness_data['classification'][0]}")
            for reason in thickness_data['classification'][1]:
                report.append(f"- {reason}")
            
            # Unique filename for each comb_id
            unique_filename = f"report_for_comb_{cid}.txt"
            report_path = os.path.join(cell_data['classification'][2], unique_filename)

            with open(report_path, 'w') as file:
                file.write('\n'.join(report))
            
            # Store classification and path in the dictionary
            classification_results[cid] = {
                'classification': classification,
                'report_path': report_path
            }
        return classification_results

class Plotter(DataGatherer):
    """
    The `Plotter` class provides functionality to visualize and save plots based on datasets of different parameters.
    
    Attributes:
    - organized_files: A dictionary of file basenames mapped to corresponding file lists.
    - output_folder: The directory where the generated plots will be saved.
    
    Methods:
    - plot_specific_parameters: Plots specific parameters from the provided data.
    - plot_data: Generates plots based on the parameters in the provided datasets, 
                 categorizing the data by 'thickness', 'growth', or other types. 
                 For each category, it plots individual data, computes the mean, 
                 and visualizes the mean values. The plots are then saved to the specified output directory.
    """
    def __init__(self, organized_files, output_folder):
        """
        Initialize the Plotter instance with organized files and output folder.

        Parameters:
        - organized_files (dict): Dictionary of files organized by their categorization.
        - output_folder (str): Directory where the plots will be saved.
        """
        self.organized_files = organized_files
        self.output_folder = output_folder
    
    def plot_specific_parameters(self, data, ax, y_columns, title):
        """
        Plot specified parameters on the given axis.

        Parameters:
        - data (pd.DataFrame): DataFrame containing data to plot.
        - ax (matplotlib.Axes): Axes object where the plot will be generated.
        - y_columns (list): List of columns to plot.
        - title (str): Title of the plot.
        """
        # Plot each specified column with an alpha value
        for col in y_columns:
            if col in data.columns:
                ax.plot(data["Time"], data[col], label=col, alpha=0.5)  # each replicate
        
        # Plot the mean of specified columns
        mean_data = data[y_columns].mean(axis=1)
        ax.plot(data["Time"], mean_data, 'k--', label='Mean', linewidth=2)  # mean line
        ax.set_title(title)
        ax.legend()
    
    def plot_data(self):
        """
        Plot the organized data based on the type of the file (thickness, growth, etc.).
        """
        try:
            colormap = plt.cm.tab10.colors

            # Loop through each category of organized files
            for (comb, file_basename), file_list in organized_files.items():
                all_data = []
                y_columns = []
                y_colors = {}

                # Decide on the plotting structure based on file type
                if 'thickness' in file_basename:
                    fig, axs = plt.subplots(3, 1, figsize=(10, 15))
                elif 'growth' in file_basename:
                    fig, axs = plt.subplots(4, 2, figsize=(15, 20))
                else:
                    plt.figure(figsize=(10, 6))

                # Loop through each file in the current category
                for file in file_list:
                    # print(file)            
                    data = pd.read_csv(file)
                    all_data.append(data)
                    x = data["Time"]
                    y_columns = [col for col in data.columns if col != "Time"]

                    # Assign colors to columns
                    for idx, col in enumerate(y_columns):
                        if col not in y_colors:
                            y_colors[col] = colormap[idx % len(colormap)]

                    # Plotting for thickness data
                    if 'thickness' in file_basename:
                        # Columns to be plotted for thickness data
                        limbus_params = ['Limbus Avg Thickness', 'Limbus Left Thickness', 'Limbus Right Thickness']
                        periphery_params = ['Periphery Avg Thickness', 'Periphery Left Thickness', 'Periphery Right Thickness']
                        
                        for ax, params in zip(axs, [y_columns, limbus_params, periphery_params]):
                            for col in params:
                                if col in y_colors:
                                    ax.plot(x, data[col], color=y_colors[col], alpha=0.2)

                    # Plotting for growth data
                    elif 'growth' in file_basename:
                        # Columns to be plotted for growth data
                        stem_params_group = [
                            ['Stem Average Total', 'Stem Average EGF', 'Stem Average Density'],
                            ['Stem Average EGF', 'Stem Average EGFseen', 'Stem Average Total'],
                            ['Stem Average Density', 'Stem Average Pressure', 'Stem Average Total'],
                            ['Stem Average TargetVolume', 'Stem Average Total']
                        ]
                        basal_params_group = [
                            [param.replace('Stem', 'Basal') for param in group]
                            for group in stem_params_group
                        ]

                        # Plotting the individual stem data with transparency
                        for data in all_data:
                            for i, stem_params in enumerate(stem_params_group):
                                ax = axs[i][0]  # Stem on the left column
                                for col in stem_params:
                                    ax.plot(x, data[col], color=y_colors[col], alpha=0.2)
                                if i in [1, 2, 3]:
                                    ax.set_yscale("log")

                        # Plotting the individual basal data with transparency
                        for data in all_data:
                            for i, basal_params in enumerate(basal_params_group):
                                ax = axs[i][1]  # Basal on the right column
                                for col in basal_params:
                                    ax.plot(x, data[col], color=y_colors[col], alpha=0.2)
                                if i in [1, 2, 3]:
                                    ax.set_yscale("log")

                    # Plotting for any other type of data
                    else:
                        for idx, col in enumerate(y_columns):
                            plt.plot(x, data[col], color=colormap[idx % len(colormap)], alpha=0.2)

                # Compute the mean data outside the file loop
                mean_data = pd.concat(all_data).groupby(level=0).mean()

                # Plotting mean data for 'thickness' files
                if 'thickness' in file_basename:
                    for ax, params, title in zip(axs, [y_columns, limbus_params, periphery_params], ["All Parameters", "Limbus", "Periphery"]):
                        for col in params:
                            if col in y_colors:
                                ax.plot(x, mean_data[col], label=col, color=y_colors[col])
                        ax.set_title(title)
                        ax.legend()
                    
                    plt.xlabel("Time (Hour)")
                    comb, rep = self.extract_params_from_path(file_list[0])
                    plt.suptitle(file_basename + f' Combination {comb}', y=1.02)
                    plt.tight_layout(h_pad=2.5)

                # Plotting mean data for 'growth' files
                elif 'growth' in file_basename:
                    # Plotting the mean data for 'Stem'
                    for i, stem_params in enumerate(stem_params_group):
                        ax = axs[i][0]  # Stem on the left column
                        for col in stem_params:
                            ax.plot(x, mean_data[col], label=col, color=y_colors[col])
                        ax.legend()
                        ax.set_title("Stem - " + ", ".join(stem_params))
                        if i in [1, 2, 3]:  # for second, third, and fourth rows
                            ax.set_yscale("log")

                    # Plotting the mean data for 'Basal'
                    for i, basal_params in enumerate(basal_params_group):
                        ax = axs[i][1]  # Basal on the right column
                        for col in basal_params:
                            ax.plot(x, mean_data[col], label=col, color=y_colors[col])
                        ax.legend()
                        ax.set_title("Basal - " + ", ".join(basal_params))
                        if i in [1, 2, 3]:  # for second, third, and fourth rows
                            ax.set_yscale("log")
                    plt.tight_layout()

                # Plotting mean data for any other type of files
                else:
                    for idx, col in enumerate(y_columns):
                        plt.plot(x, mean_data[col], label=col, color=colormap[idx % len(colormap)])

                    plt.xlabel("Time (Hour)")
                    comb, rep = self.extract_params_from_path(file_list[0])
                    plt.title(file_basename + f' Combination {comb}')
                    plt.legend() 

                # print(comb, file_basename)

                # Saving the plot to the specified directory
                plot_filename = file_basename + f' Combination {comb}.png'
                comb_output_folder = os.path.join(self.output_folder, f'Comb{comb}')
                # if not os.path.exists(comb_output_folder):  # Create the combination folder if it doesn't exist
                #     os.makedirs(comb_output_folder)
                plt.savefig(os.path.join(comb_output_folder, plot_filename), bbox_inches='tight')
            
        except FileNotFoundError:
            sys.stderr.write(f"Error: File {file} not found.")
        except pd.errors.EmptyDataError:
            sys.stderr.write(f"Error: File {file} contains no data.")
        except Exception as e:
            sys.stderr.write(f"An error occurred: {e}")
        finally:
            plt.close()    
            

class FileOrganizer:
    
    def move_comb_folders(self, classification_results):
        # Go through each comb_id and its results
        for cid, data in classification_results.items():
            classification = data['classification']
            comb_dir = os.path.dirname(data['report_path'])

            # Create the target directory based on classification
            target_dir = os.path.join(os.path.dirname(os.path.dirname(data['report_path'])), classification)
            if not os.path.exists(target_dir):
                with FileLock(f"{target_dir}.lock"):
                    # Make the directory if it doesn't exist after obtaining the lock
                    if not os.path.exists(target_dir):
                        os.makedirs(target_dir)
            # Compute the destination path
            destination_path = os.path.join(target_dir, os.path.basename(comb_dir))
            # Move the comb folder to the target directory
            with FileLock(f"{comb_dir}.lock"):
                if os.path.exists(comb_dir):
                    shutil.move(comb_dir, destination_path)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process specific Comb# folders.")
    parser.add_argument("folders", nargs='+', help="List of Comb# folders to process")
    args = parser.parse_args()    
    
    gatherer = DataGatherer(base_folder=r"C:\Users\joelv\OneDrive\Desktop\Output_20230918-122031")
    csv_files = gatherer.gather_csv_files_from_folders(args.folders)
    organized_files = gatherer.organize_files(csv_files)    
    
    # plotter = Plotter(organized_files=organized_files, output_folder=r"C:\Users\joelv\OneDrive\Desktop\Output_20230918-122031")
    # plotter.plot_data()    
    
    classifier = Classifier(organized_data=organized_files)
    results = classifier.classify_all()    
    classification = classifier.generate_combined_report(results)

    plotter = Plotter(organized_files=organized_files, output_folder=r"C:\Users\joelv\OneDrive\Desktop\Output_20230918-122031")
    plotter.plot_data() 

    # Organizing the folders
    FileOrganizer().move_comb_folders(classification)
    

