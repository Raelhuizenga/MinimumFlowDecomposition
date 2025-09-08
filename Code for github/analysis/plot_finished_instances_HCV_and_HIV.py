import matplotlib.pyplot as plt
import numpy as np
import os

# HCV data
x_hcv = np.array([3, 4, 5, 6, 7])
mfd_hcv = np.array([80, 50, 50, 60, 20])
mfd_fast_hcv = np.array([80, 50, 20, 30, 20])
perfect_weights_hcv = [100, 100, 100, 100, 100]
k_means_weights_hcv = [100, 100, 70, 50, 30]

# HIV data (example values, replace with your actual data)
x_hiv = np.array([2, 3, 4, 5])
mfd_hiv = np.array([100, 87.5, 30, 25])
mfd_fast_hiv = np.array([80, 75, 20, 0])
perfect_weights_hiv = [100, 100, 100, 100]
k_means_weights_hiv = [90, 62.5, 90, 50]

bar_width = 0.25
offset = [-bar_width, 0, bar_width]
colors = ['salmon', '#8FB9E6', '#BFCFC2']

fig, axs = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

def plot_finished_instances(ax, x, mfd, mfd_fast, perfect_weights, k_means_weights, title):
    x_indices = np.arange(len(x))
    
    ax.bar(x_indices + offset[0], mfd_fast, width=bar_width,
           color=colors[0], hatch='///', edgecolor='black', label='MFD (correct)')
    
    mfd_top = mfd - mfd_fast
    mfd_top[mfd_top < 0] = 0
    ax.bar(x_indices + offset[0], mfd_top, width=bar_width,
           bottom=mfd_fast, color=colors[0], edgecolor='black', label='MFD')
    
    ax.bar(x_indices + offset[1], k_means_weights, width=bar_width,
           label='k-means weights', color=colors[1], edgecolor='black')
    ax.bar(x_indices + offset[2], perfect_weights, width=bar_width,
           label='perfect weights', color=colors[2], edgecolor='black')
    
    ax.set_xlabel('Number of haplotypes', fontsize=14)
    ax.set_title(title, fontsize=15)
    ax.set_xticks(x_indices)
    ax.set_xticklabels(x)
    ax.tick_params(axis='both', labelsize=12)
    ax.set_ylim(0, 105)
    ax.grid(True, axis='y', linestyle='--', alpha=0.5)

# Plot HCV data
plot_finished_instances(axs[0], x_hcv, mfd_hcv, mfd_fast_hcv, perfect_weights_hcv, k_means_weights_hcv, 'Finished HCV samples per method')
axs[0].set_ylabel('Finished samples (%)', fontsize=14)

# Plot HIV data
plot_finished_instances(axs[1], x_hiv, mfd_hiv, mfd_fast_hiv, perfect_weights_hiv, k_means_weights_hiv, 'Finished HIV samples per method')

# Make y-axis labels and ticks visible on the right plot (HIV)
axs[1].yaxis.set_ticks_position('both')
axs[1].tick_params(axis='y', labelleft=True)

# Adjust layout BEFORE adding legend
plt.tight_layout(rect=[0, 0, 0.85, 1])

# Add legend on the right side of the figure
handles, labels = axs[0].get_legend_handles_labels()
fig.legend(handles, labels, title="Method", bbox_to_anchor=(0.85, 0.9), loc='upper left', fontsize=12, title_fontsize=14)

# Save with bbox_inches='tight' so legend is not cut off
output_dir = '/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/graphs/finished_graphs'
os.makedirs(output_dir, exist_ok=True)
output_path = os.path.join(output_dir, 'finished_instances_barplot_color_weights_striped_side_by_side_colors adjusted.png')
plt.savefig(output_path, dpi=300, bbox_inches='tight')

plt.show()
