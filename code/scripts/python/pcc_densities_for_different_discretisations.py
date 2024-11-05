import matplotlib
matplotlib.rcParams['font.family']='serif'
matplotlib.rcParams['mathtext.fontset']='cm'

import numpy as np
import matplotlib.pyplot as plt

import os
import csv

def load_csv_files_to_lists(directory):
    # Dictionary to store 2D lists, where keys are file names
    csv_data = {}

    # Loop through all files in the directory
    for filename in os.listdir(directory):
        filepath = os.path.join(directory, filename)
            
        # Read the CSV file and load it into a 2D list
        with open(filepath, mode='r', newline='') as file:
            reader = csv.reader(file)
            data = [row for row in reader]  # 2D list for the CSV content
            
            # Store the data in the dictionary
            csv_data[filename] = data
    
    return csv_data

dendritic_length = 500 #um

# Example usage:
# directory_path = '../../data/10_percent_transcription_rate/'  # Replace with the actual path
# directory_path = '../../data/normal/'  # Replace with the actual path
# directory_path = '../../data/0p1_percent_transcription_rate/'
directory_path = '../../data/0p1_percent_translation_rate/'
csv_files_data = load_csv_files_to_lists(directory_path)

fig, axs = plt.subplots(nrows=1, ncols=2, figsize=(10, 5))

# Output the loaded data (for demonstration)
pcc_soma_first_ds = []
n_ds = []
pcc_density_soma_first_ds = []
for filename, data in csv_files_data.items():
    print(f"Data from file: {filename}")
    n_ds_ = len(csv_files_data[filename][1]) - 3
    n_ds.append(n_ds_)
    # pcc_soma_first_ds.append(float(csv_files_data[filename][1][len(csv_files_data[filename][1]) - 3]))
    pcc_soma_first_ds.append(float(csv_files_data[filename][1][2]))
    pcc_density_soma_first_ds.append(float(csv_files_data[filename][1][2])*n_ds_/dendritic_length)

axs[0].set_title("PCC")
axs[0].plot(n_ds, pcc_soma_first_ds, linestyle='none', marker='*')
axs[0].grid()
axs[0].set_xlabel("Number of compartments")
axs[0].set_ylabel("PCC")

axs[1].set_title("PCC Density")
axs[1].plot(n_ds, pcc_density_soma_first_ds, linestyle='none', marker='*')
axs[1].grid()
axs[1].set_xlabel("Number of compartments")
axs[1].set_ylabel("PCC density $[\mu m^{-1}]$")


plt.tight_layout()
plt.show()
