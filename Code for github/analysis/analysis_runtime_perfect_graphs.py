import sys
import glob
import json
import matplotlib.pyplot as plt
import numpy as np
from collections import defaultdict
import os

def main():
    # t, g_length = read_results_from_file_length("output/solutions/genomelength/*.json")
    # plot_genomelength_vs_time(t, g_length)

    # t, theoretical, haps = read_results_from_file()
    # plot_haps_vs_time(t, theoretical, haps)

    result_by_strain = read_results_by_strain()
    # plot_runtime_vs_nodes_by_strains(result_by_strain)
    plot_runtime_vs_length(result_by_strain)


def plot_runtime_vs_nodes(results_by_strain):
    print(results_by_strain.keys())
    nodes, runtimes, _ = zip(*sorted(results_by_strain.get(3, results_by_strain.get("3", []))))
    plt.scatter(nodes, runtimes, marker='o', color='black')
    plt.xlabel("Number of nodes")
    plt.ylabel("Runtime (s)")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/haplotype_gridsearch/length/runtime_vs_nodes_3_strains.png")
    plt.show()


def plot_runtime_vs_length(results_by_strain):
    plt.rcParams.update({
        'font.family': 'DejaVu Sans',
        'font.size': 12,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11
    })
    plt.figure(figsize=(8, 5)) 

    print(results_by_strain.keys())
    nodes, runtimes, length = zip(*sorted(results_by_strain.get(3, results_by_strain.get("3", []))))
    plt.scatter(length, runtimes, marker='o', color='black', alpha=0.5)

    # Fit a linear line: y = m*x + b
    slope, intercept = np.polyfit(length, runtimes, 1)
    fit_line = np.poly1d((slope, intercept))
    x_vals = np.linspace(min(length), max(length), 100)
    plt.plot(x_vals, fit_line(x_vals), color='red', label=f"Fit: y={slope:.5f}x+{intercept:.2f}")

    plt.xlabel("Length haplotypes")
    plt.ylabel("Runtime (s)")
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.legend()
    plt.tight_layout()
    plt.savefig("/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/haplotype_gridsearch/length/runtime_vs_length_more_points.png")
    plt.show()


def read_results_by_strain():
    solution_files = glob.glob("/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/haplotype_gridsearch/length/*.json")
    # solution_files = glob.glob("/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition/output_simulations/haplotypes_grid_search/flowpackage/*.json")
    # solution_files = glob.glob("/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition/output_simulations/grid_search/*.json")

    results_by_strain = {}

    for file in solution_files:
        filename = os.path.basename(file)
        length = filename.strip(".json").split("_")[2][1::]
        with open(file, "r") as f:
            data = json.load(f)
            time = data["time"]
            nodes = data["num_vertices"]
            num_strains = len(data["paths"])
            if num_strains == 0:
                continue
            results_by_strain.setdefault(num_strains, []).append((nodes, time, int(length)))
    return results_by_strain


def plot_runtime_vs_strains_with_error(results_by_strain):
    # Use default clean sans-serif font
    plt.rcParams.update({
        'font.family': 'DejaVu Sans',
        'font.size': 12,
        'axes.titlesize': 14,
        'axes.labelsize': 12,
        'xtick.labelsize': 11,
        'ytick.labelsize': 11
    })

    strain_counts = sorted(results_by_strain.keys())
    medians = []
    lower_errors = []
    upper_errors = []

    for strain_count in strain_counts:
        runtimes = [runtime for _, runtime in results_by_strain[strain_count]]
        runtimes = np.array(runtimes)
        median = np.median(runtimes)
        q25 = np.percentile(runtimes, 25)
        q75 = np.percentile(runtimes, 75)

        medians.append(median)
        lower_errors.append(median - q25)
        upper_errors.append(q75 - median)

    # Combine lower and upper errors into format expected by errorbar
    asymmetric_error = [lower_errors, upper_errors]

    plt.figure(figsize=(8, 5))
    plt.errorbar(
        strain_counts,
        medians,
        yerr=asymmetric_error,
        fmt='o-',
        capsize=5,
        linewidth=2,
        markersize=6,
        color='black',
        ecolor='gray',
        elinewidth=1
    )

        # Add star annotations
    for strain, label in [(8, "*"), (9, "**")]:
        if strain in strain_counts:
            idx = strain_counts.index(strain)
            plt.annotate(
                label,
                xy=(strain, medians[idx]),
                xytext=(-10, 10),
                textcoords='offset points',
                ha='center',
                fontsize=14,
                color='black',
                weight='bold'
            )

    plt.yscale('log')
    plt.xlabel("Number of haplotypes", fontsize=16)
    plt.ylabel("Median runtime (s) (with IQR)", fontsize=16)
    plt.title("The runtime increases with the number of haplotypes", fontsize=20, pad=15)
    plt.grid(True, linestyle='--', alpha=0.6)
    plt.ylim(bottom=0)
    plt.tick_params(labelsize=12)
    plt.xticks(strain_counts)  # Ensure integer ticks on x-axis
    # plt.savefig("/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition/output_simulations/plots/runtime_vs_strains_errorbars.png")
    plt.savefig("/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/haplotype_gridsearch/runtime_vs_strains_error_haplotypes_log.png")
    plt.show()



def plot_runtime_vs_nodes_by_strains(results_by_strain):
    """
    Plot runtime vs number of nodes for different strain counts.

    results_by_strain: dict
        Key = strain count (e.g., 3, 4, 5)
        Value = list of tuples (number_of_nodes, runtime)
    """
    for strain_count, values in sorted(results_by_strain.items()):
        nodes, runtimes = zip(*sorted(values))  # Sort by node count
        plt.plot(nodes, runtimes, marker='o', label=f"{strain_count} strains")

    plt.xlabel("Number of nodes")
    plt.ylabel("Runtime (s)")
    plt.title("The runtime of a MFD depends on the complexity of the graph")
    plt.legend(title="Strains")
    plt.grid(True)
    plt.tight_layout()
    plt.savefig("/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/haplotype_gridsearch/length/runtime_vs_nodes_by_strain.png")
    plt.show()


def read_results_from_file():
    solution_files = glob.glob("/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition/output_simulations/grid_search/*.json")  # Change path if needed

    times = []
    theoretical_haps = []
    num_haps = []
    num_haps_theoretical = 2

    # Read each file and extract the time
    for file in solution_files:
        with open(file, "r") as f:
            data = json.load(f)
            times.append(data["time"])
            num_haps.append(len(data["weights"]))
            theoretical_haps.append(num_haps_theoretical)
            num_haps_theoretical += 1
    print(times)
    return times, theoretical_haps, num_haps



def read_results_from_file_length(file_path):
    solution_files = glob.glob(file_path)  # Change path if needed

    times = []
    genome_lengths = []
    genome_length = 10

    # Read each file and extract the time
    for file in solution_files:
        with open(file, "r") as f:
            data = json.load(f)
            times.append(data["time"])
            genome_lengths.append(genome_length)
            genome_length += 50
    return times, genome_lengths

def plot_genomelength_vs_time(times, genome_lengths):

    plt.plot(genome_lengths, times)
    plt.xlabel("Genome length")
    plt.ylabel("Time (s)")
    plt.title("Time vs genome length")
    plt.show()
    plt.savefig('output/plot_genome_length.png')


def plot_haps_vs_time(times, theoretical_haps, num_haps):
    plt.plot(theoretical_haps, times)
    plt.xlabel("Number of haplotypes")
    plt.ylabel("Time (s)")
    plt.title("Time vs theortical number of haplotypes")
    plt.show()
    plt.savefig('output/plot_theorethical_haps.png')

    plt.scatter(num_haps, times, label="Real")
    plt.xlabel("Number of haplotypes")
    plt.ylabel("Time (s)")
    plt.title("Time vs found number of haplotypes")
    plt.show()
    plt.savefig('output/plot_found_haps.png')

if __name__ == '__main__':
    sys.exit(main())