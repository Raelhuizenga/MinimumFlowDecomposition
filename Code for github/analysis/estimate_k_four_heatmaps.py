import json
import os
from glob import glob
from collections import Counter
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# In seqwish container

def load_heatmap_data(json_dir):
    json_files = glob(os.path.join(json_dir, '*.json'))
    points = []

    for filepath in json_files:
        with open(filepath, 'r') as f:
            try:
                entry = json.load(f)
                points.append((entry['strain'], entry['found_paths']))
            except (json.JSONDecodeError, KeyError):
                print(f"Skipping malformed or incomplete JSON file: {filepath}")

    counter = Counter(points)
    df = pd.DataFrame.from_dict(counter, orient='index', columns=['count'])
    df.index = pd.MultiIndex.from_tuples(df.index, names=['strain', 'found_paths'])
    df = df.reset_index()
    return df

# Directories
dir1 = '/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/estimate_k/HIV_8000_dagify'
dir2 = '/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/estimate_k_no_constraints/HIV_8000_nodes_filtered'
# Load other two datasets (assuming directories exist)
dir3 = '/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/estimate_k/HCV_nodes_removed'
dir4 = '/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/estimate_k_no_constraints/HCV_nodes_removed'

# Load data
data1 = load_heatmap_data(dir1)
data2 = load_heatmap_data(dir2)
data3 = load_heatmap_data(dir3)
data4 = load_heatmap_data(dir4)

# Axis ranges for top and bottom rows
x_range_top = list(range(2, 6))
y_range_top = list(range(2, 9))

x_range_bottom = list(range(3, 8))
y_range_bottom = list(range(3, 11))

def create_matrix(df, x_range, y_range):
    mat = df.pivot(index='found_paths', columns='strain', values='count').fillna(0)
    mat = mat.reindex(index=sorted(y_range, reverse=True), columns=x_range, fill_value=0)
    return mat


# Create matrices
matrix1 = create_matrix(data1, x_range_top, y_range_top)
matrix2 = create_matrix(data2, x_range_top, y_range_top)
matrix3 = create_matrix(data3, x_range_bottom, y_range_bottom)
matrix4 = create_matrix(data4, x_range_bottom, y_range_bottom)

# Shared color scale
vmax = max(matrix1.values.max(), matrix2.values.max(), matrix3.values.max(), matrix4.values.max())

# Create 2x2 plot
fig, axs = plt.subplots(2, 2, figsize=(12, 12), sharex=False, sharey=False)

def plot_heatmap(matrix, ax, ylabel):
    sns.heatmap(matrix, annot=True, fmt=".0f", cmap='Greys', cbar=False,
                ax=ax, vmin=0, vmax=vmax, square=True, annot_kws={"size": 12})
    ax.set_xlabel('Number of haplotypes', fontsize=14)
    if ylabel:
        ax.set_ylabel('Number of found paths', fontsize=14)
    else:
        ax.set_ylabel('', fontsize=14)
    ax.tick_params(labelsize=12)

# Plot all 4 heatmaps with HCV on top, HIV on bottom
plot_heatmap(matrix3, axs[0, 0], True)   # HCV with subpaths
plot_heatmap(matrix4, axs[0, 1], False)  # HCV no subpaths
plot_heatmap(matrix1, axs[1, 0], True)   # HIV with subpaths
plot_heatmap(matrix2, axs[1, 1], False)  # HIV no subpaths


# Subfigure labels
fig.text(0.23, 0.93, '(a) HCV with subpaths     ', ha='left', fontsize=15)
fig.text(0.75, 0.93, '(b) HCV no subpaths', ha='right', fontsize=15)
fig.text(0.23, 0.46, '(c) HIV with subpaths     ', ha='left', fontsize=15)
fig.text(0.75, 0.46, '(d) HIV no subpaths', ha='right', fontsize=15)

# Layout
fig.subplots_adjust(wspace=-0.2, hspace=0.5)

# Save and show
plt.savefig("/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/graphs/heatmaps/heatmaps_comparison_four_HCV_top.png", dpi=300)
plt.show()
