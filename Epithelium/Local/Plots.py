#!/usr/bin/env python
"""
Example script to run PlotManager from the command line.

Usage:
    python plot_manager.py /path/to/data 360 720 --out_dir /path/to/output
"""

import argparse
import sys
import os
from pathlib import Path

# Third-party libraries
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import skew, kurtosis, kstest, shapiro, probplot
import scipy.stats as stats

# -------------------------------------------------------------------
# Define the PlotManager class
# -------------------------------------------------------------------
class PlotManager:
    def __init__(self, Data_DIR, Stable_time=500, Max_time=720, Out_DIR=None):
        # Convert string paths to Path objects
        self.data_dir = Path(Data_DIR)
        self.out_dir = Path(Out_DIR) if Out_DIR else self.data_dir
        
        # Define color codes and other constants
        self.color_codes = ['#55ffff', '#0055ff', '#ffbe99', '#ff007f']
        self.stable_time = Stable_time  # e.g., 15*24 = 360 hours
        self.max_time = Max_time        # e.g., 720 hours

        # Collect and load data
        self.Param_file, self.Count_file, self.Thickness_file = self.collect_files(self.data_dir)
        self.memb_bins = self.load_memb_data()
        print("Membrane Bins:", self.memb_bins)

        self.count_data = self.data_collection(self.Count_file)
        self.thickness_data = self.data_collection(self.Thickness_file)
    
    def collect_files(self, directory):
        parameter_files = []
        count_files = []
        thickness_bins_files = []
        
        directory = Path(directory)
        for file in directory.rglob("*"):
            if file.is_file():
                if 'thickness_rep' in file.name:
                    thickness_bins_files.append(str(file))
                elif 'Parameters.py' in file.name:
                    parameter_files.append(str(file))
                elif "cell_count_" in file.name:
                    count_files.append(str(file))

        return parameter_files, count_files, thickness_bins_files

    def data_collection(self, file_paths):
        """
        Collect data from multiple files, supporting both CSV and Parquet formats.

        Args:
            file_paths (list of str): List of file paths

        Returns:
            dict: Dictionary containing DataFrames with keys based on path components
        """
        data_dict = {}
        
        for data_path in file_paths:
            file_extension = Path(data_path).suffix.lower()
            try:
                # Read file based on extension
                if file_extension == '.csv':
                    data = pd.read_csv(data_path)
                elif file_extension == '.parquet':
                    data = pd.read_parquet(data_path)
                else:
                    raise ValueError(
                        f"Unsupported file format: {file_extension}. "
                        "Only .csv and .parquet files are supported."
                    )
                
                # Create a dict key from parts of the path (adjust as needed)
                # Below, we are grabbing the "last-4" and "last-3" directory parts:
                # e.g., if path = .../some_dir/another_dir/file.csv
                # parts[-4], parts[-3] might be "some_dir_another_dir"
                key = f"{Path(data_path).parts[-4]}_{Path(data_path).parts[-3]}"
                data_dict[key] = data

            except Exception as e:
                print(f"Error reading file {data_path}: {str(e)}")
                continue
                
        return data_dict
    
    def load_memb_data(self):
        memb_data = {}
        bins_memb = None
        try:
            memb_file = self.data_dir.joinpath("membrane_height.csv")
            rows = pd.read_csv(memb_file, header=None)
            # Rows is a DataFrame with columns [0,1]
            for _, row in rows.iterrows():
                key = float(row[0])
                val = float(row[1])
                memb_data[key] = val
                
            bin_size = 200 / 10
            bin_indexes = list(range(10))
            bins_memb = {index: [] for index in bin_indexes}
            
            for key in memb_data.keys():
                bin_index_memb = int(key // bin_size)
                bins_memb[bin_index_memb].append(memb_data[key])
            
            for bin_index_memb in bins_memb.keys():
                if bins_memb[bin_index_memb]:
                    bins_memb[bin_index_memb] = np.mean(bins_memb[bin_index_memb])
                else:
                    bins_memb[bin_index_memb] = 0

        except Exception as e:
            print(f"Warning: Could not load membrane data: {str(e)}")
            bins_memb = None
        
        return bins_memb

    def plot_count(self, data_dict):
        # Merge all replicates
        combined_data = pd.concat(data_dict.values())
        combined_data = combined_data[combined_data['Time'] <= self.max_time]
        combined_stable_data = combined_data[combined_data['Time'] >= self.stable_time]
        mean_data = combined_data.groupby('Time').mean(numeric_only=True).reset_index()

        fig, axes = plt.subplots(4, 2, figsize=(14, 12))
        # Turn off unused subplots
        axes[0, 0].axis('off')
        axes[0, 1].axis('off')
        axes[3, 0].axis('off')
        axes[3, 1].axis('off')

        ax_longitudinal = plt.subplot2grid((4, 2), (0, 0), colspan=2)  
        ax_qq_combined = plt.subplot2grid((4, 2), (3, 0), colspan=2, rowspan=1)

        # Plot mean data lines
        # Assuming your data columns are ["Time", "ColA", "ColB", "ColC", "ColD", ...]
        # You might adjust based on actual column names
        for i, column in enumerate(mean_data.columns[1:5]):
            ax_longitudinal.plot(mean_data["Time"] / 24, mean_data[column],
                                 label=f'Mean {column}',
                                 color=self.color_codes[i], linewidth=2)

            # Determine row/col in the 2D subplot grid for histograms
            row = (i // 2) + 1
            col = i % 2
            ax_hist = axes[row, col]

            # Gather stable data for the current column
            stable_arrays = []
            for dkey, df in data_dict.items():
                stable_df = df[(df['Time'] >= self.stable_time) & (df['Time'] <= self.max_time)]
                stable_arrays.append(stable_df[column].values)

            all_stable_data = np.concatenate(stable_arrays)

            bin_edges = np.arange(int(all_stable_data.min()) - 0.5,
                                  int(all_stable_data.max()) + 1.5,
                                  1)
            # Plot each replicate histogram
            all_hists = []
            for _, df in data_dict.items():
                stable_df = df[(df['Time'] >= self.stable_time) & (df['Time'] <= self.max_time)]
                hist, _, _ = ax_hist.hist(stable_df[column],
                                          bins=bin_edges,
                                          color=self.color_codes[i],
                                          alpha=0.075)
                all_hists.append(hist)

            mean_hist = np.mean(all_hists, axis=0)  # not used, but you can keep if needed
            ax_hist.grid(True, linewidth=0.1)
            ax_hist.set_title(f'Variation of {column}-Cell Number')
            ax_hist.set_xlabel('Number of Cells per Time')
            ax_hist.set_ylabel('Frequency')

            # Print stats to console
            print(f"\n--- Stats for {column} ---")
            print(f"{column} Skewness: {skew(all_stable_data):.2f}")
            print(f"{column} Kurtosis: {kurtosis(all_stable_data):.2f}")

            ks_stat, ks_p = kstest(all_stable_data, 'norm',
                                   args=(np.mean(all_stable_data), np.std(all_stable_data)))
            print(f"{column} KS stat: {ks_stat:.2f}, KS p-value: {ks_p:.2e}")

            sw_stat, sw_p = shapiro(all_stable_data)
            print(f"{column} SW stat: {sw_stat:.2f}, SW p-value: {sw_p:.2e}")

            anderson_result = stats.anderson(all_stable_data, dist='norm')
            print(f"{column} A stat: {anderson_result.statistic:.2f}")
            print(f"{column} A Crit: {anderson_result.critical_values}")
            print(f"{column} A signi: {anderson_result.significance_level}")

            # Q-Q plot
            (osm, osr), (slope, intercept, r) = probplot(combined_stable_data[column], dist="norm", plot=ax_qq_combined)

            lines = ax_qq_combined.get_lines()
            # Each probplot() call typically adds 2 lines (data + theoretical line),
            # so we offset them by 2*i
            if len(lines) >= 2 * (i + 1):
                observed_line = lines[2 * i]
                expected_line = lines[2 * i + 1]

                observed_line.set_color(self.color_codes[i])
                observed_line.set_alpha(0.25)
                observed_line.set_label(column)

                expected_line.set_color('black')
                expected_line.set_linestyle('dashed')

                # Some stats for the Q-Q plot
                differences = np.abs(osm - osr)
                mean_diff = np.mean(differences)
                max_diff = np.max(differences)
                print(f"{column} Q-Q Slope: {slope:.2f}, Intercept: {intercept:.2f}, R: {r:.2f}")
                print(f"{column} Q-Q Mean diff: {mean_diff:.2f}, Max diff: {max_diff:.2f}")
        
        # Finalize Q-Q combined subplot
        ax_qq_combined.set_title('Combined Q-Q Plots')
        ax_qq_combined.grid(True, linewidth=0.1)
        ax_qq_combined.legend()

        # Plot replicate lines behind the mean lines
        for dkey, df in data_dict.items():
            df_filtered = df[df['Time'] <= self.max_time]
            for i, column in enumerate(df_filtered.columns[0:4]):
                ax_longitudinal.plot(df_filtered['Time'] / 24,
                                     df_filtered[column],
                                     color=self.color_codes[i],
                                     alpha=0.075, linewidth=2)

        ax_longitudinal.axvline(self.stable_time / 24, color='black',
                                linestyle='dashed', linewidth=1,
                                label="Stability Time")
        ax_longitudinal.legend(loc='center right')
        ax_longitudinal.grid(True, linewidth=0.1)
        ax_longitudinal.set_title('Longitudinal Cell Count')
        ax_longitudinal.set_xlabel('Time (Days)')
        ax_longitudinal.set_ylabel('Number of Cells')

        plt.tight_layout()
        plt.subplots_adjust(hspace=0.40)
        plt.show()

    def plot_thickness(self, data_dict, bin_memb):
        """
        Plot thickness-related data. 
        data_dict should be a dictionary of DataFrames, each with:
            - 'Time' column
            - 'Bin' column
            - 'Height' column (some measure of thickness)
        bin_memb: dictionary with membrane reference data 
        """
        if not data_dict or bin_memb is None:
            print("No thickness data or membrane bins to plot.")
            return

        unique_bins = set()
        for df in data_dict.values():
            unique_bins.update(df['Bin'].unique())
        unique_bins = sorted(unique_bins)

        color_palette = sns.color_palette("colorblind", len(unique_bins))
        bin_color_map = {bin_key: color_palette[i] for i, bin_key in enumerate(unique_bins)}

        bins_size = 200 // len(unique_bins)
        bin_x = {bin_key: [] for bin_key in unique_bins}

        for key in bin_x.keys():
            # Make an approximate horizontal coordinate
            bin_x[key] = [i for i in range(key * bins_size, int((key + 1) * bins_size))]
            bin_x[key] = np.mean(bin_x[key])

        fig, ax = plt.subplots(2, figsize=(14, 12))

        all_data = []
        for comb, df in data_dict.items():
            all_data.append(df)
            for bin_key in unique_bins:
                bin_data = df[df['Bin'] == bin_key]
                ax[0].plot(bin_data['Time'] / 24,
                           bin_data['Height'] * 2,  # if "Height" is half-thickness, multiply by 2
                           alpha=0.075, linewidth=1, color=bin_color_map[bin_key])

        combined_data = pd.concat(all_data)
        mean_data = combined_data.groupby(['Time', 'Bin']).mean(numeric_only=True).reset_index()
        std_data  = combined_data.groupby(['Time', 'Bin']).std(numeric_only=True).reset_index()

        bin_actual = []
        bin_relative = []

        for bin_key in unique_bins:
            bin_mean_data = mean_data[mean_data['Bin'] == bin_key]
            stable_bin_data = combined_data[
                (combined_data['Bin'] == bin_key) & 
                (combined_data['Time'] >= self.stable_time)
            ]
            ax[0].plot(bin_mean_data['Time'] / 24,
                       bin_mean_data['Height'] * 2,
                       label=f'X Coord. {bin_key*bins_size}-{int((bin_key+1)*bins_size)-1}',
                       linewidth=1, color=bin_color_map[bin_key])

            # Full thickness minus membrane reference
            thickness_actual = (stable_bin_data['Height'] - bin_memb[bin_key]) * 2
            thickness_relative = stable_bin_data['Height'] * 2

            # Mean thickness across stable time
            bin_actual.append(thickness_actual.mean())
            bin_relative.append(thickness_relative.mean())

            # Box plots: place them at x-coord = bin_x[bin_key] * 2 for clarity
            box_position = bin_x[bin_key] * 2

            ax[1].boxplot(thickness_actual,
                          positions=[box_position],
                          widths=6, patch_artist=True,
                          boxprops=dict(facecolor='tab:red', color='black', linewidth=0.5),
                          medianprops=dict(color='black', linewidth=0.5),
                          whiskerprops=dict(color='black', linewidth=0.5),
                          capprops=dict(color='black', linewidth=0.5),
                          flierprops=dict(marker='.', color='red', markersize=1, alpha=0.5))

            ax[1].boxplot(thickness_relative,
                          positions=[box_position],
                          widths=6, patch_artist=True,
                          boxprops=dict(facecolor='tab:blue', color='black', linewidth=0.5),
                          medianprops=dict(color='black', linewidth=0.5),
                          whiskerprops=dict(color='black', linewidth=0.5),
                          capprops=dict(color='black', linewidth=0.5),
                          flierprops=dict(marker='.', color='red', markersize=1, alpha=0.5))

        ax[0].axvline(self.stable_time / 24, color='black', linestyle='dashed',
                      linewidth=1, label="Stability Time")
        ax[0].set_title('Longitudinal Variation of Tissue Top Position Over Time')
        ax[0].set_xlabel('Time (Days)')
        ax[0].set_ylabel('Vertical Position (μm)')
        ax[0].legend(title='Bins', loc='lower right')
        ax[0].grid(True, linewidth=0.1)

        x_values = [val * 2 for val in bin_x.values()]
        ax[1].plot(x_values, bin_actual, color='tab:red', linestyle='dotted',
                   linewidth=1, label="Tissue Thickness")
        ax[1].plot(x_values, bin_relative, color='tab:blue', linestyle='dotted',
                   linewidth=1, label="Tissue Top Position")

        ax[1].set_title('Stable Tissue Top Position and Thickness Across Sections')
        ax[1].set_xlabel('Horizontal Position (μm)')
        ax[1].set_ylabel('Vertical Position / Thickness (μm)')
        ax[1].legend(title='Bins', loc='center right')
        ax[1].grid(True, linewidth=0.1)

        plt.tight_layout()
        plt.subplots_adjust(hspace=0.20)
        plt.show()

# -------------------------------------------------------------------
# Main CLI entry point
# -------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(
        description="Run PlotManager on a given simulation directory."
    )
    parser.add_argument("data_dir",
                        help="Path to the directory containing the simulation data.")
    parser.add_argument("stable_time", type=int,
                        help="Stability time in hours (e.g., 360).")
    parser.add_argument("max_time", type=int,
                        help="Maximum time in hours (e.g., 720).")
    parser.add_argument("--out_dir", default=None,
                        help="Optional output directory for plots. Defaults to data_dir if not specified.")

    args = parser.parse_args()

    # Create an instance of PlotManager
    pm = PlotManager(
        Data_DIR=args.data_dir,
        Stable_time=args.stable_time,
        Max_time=args.max_time,
        Out_DIR=args.out_dir
    )

    # Call plot methods
    # (Adjust as needed. If you only want to call one or the other, remove the extra call.)
    pm.plot_count(pm.count_data)
    pm.plot_thickness(pm.thickness_data, pm.memb_bins)


if __name__ == "__main__":
    main()
