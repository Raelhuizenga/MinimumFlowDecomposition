import os
import matplotlib.pyplot as plt
import re
import seaborn as sns
import pandas as pd
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from collections import defaultdict, Counter
from scipy.stats import wilcoxon


# Run in seqwish container


# samples = [(3,0), (3,1), (3,2), (3,4), (3,5), (3,6), (3,7), (3,9), (4,0), (4,2), (4,5), (4,6), (4,7), (4,8), (5,1), (5,2), (5,5), (5,6), (5,7), (5,8), (6,0), (6,1), (6,2), (6,4), (6,5), (6,6), (6,9), (7,6), (7,8)]
# samples = [(3,0),  (4,5), (5,1), (7,6)]
# samples = [(3,0), (3,6), (3,7), (3,9)]
samples = [ (3,0), (3,1), (3,2), (3,4), (3,5), (3,6), (3,7), (3,9), (4,0), (4,5), (4,6), (4,7), (4,8), (5,1), (5,6), (6,1), (6,4), (6,5), (7,6), (7,8)]

# directory_solutions = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/min_cover_HCV_new"
# directory_ground_truth = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition/Data/simulation_data/simulation_data/art_reads/150bp/random_abundances"
# directory_vg_flow = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition/output_VG_flow/HCV"

directory_solutions = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/min_cover_strict_constraints/nodes_removed"
directory_ground_truth = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition/Data/simulation_data/simulation_data/art_reads/150bp/random_abundances"
directory_vg_flow = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition/output_VG_flow/HCV"


def load_fasta(filepath):
    # Even if there's no .fasta extension, we force the parser to treat it as FASTA
    with open(filepath, "r") as f:
        return list(SeqIO.parse(f, "fasta"))
    


def read_fasta_vg_flow(fasta_path):
    """
    Parses a FASTA file where the header includes a weight like '947x'
    and returns:
    - a list of SeqRecords
    - a dictionary {record.description: weight}
    """
    records = []
    id_to_weight = {}

    for record in SeqIO.parse(fasta_path, "fasta"):
        header = record.description
        match = re.search(r'(\d+)x', header)
        if match:
            weight = int(match.group(1))
        else:
            raise ValueError(f"Could not find weight like '947x' in header: {header}")
        
        records.append(record)
        id_to_weight[record.description] = weight

    return records, id_to_weight


def load_ground_truth_abundance(folder_path):
    """
    Parses the ground truth file (non-FASTA) into SeqRecords with full descriptions
    and a separate abundance dict. Returns (list_of_seq_records, abundances_dict)
    """
    gt_file = os.path.join(folder_path, "ground_truth_abundances.txt")
    if not os.path.isfile(gt_file):
        raise FileNotFoundError(f"Ground truth file not found: {gt_file}")

    seq_records = []
    abundances = {}
    current_id = None
    current_seq = []
    current_description = ""

    with open(gt_file, "r") as f:
        for line in f:
            if line.startswith(">"):
                if current_id is not None:
                    # Store the previous record
                    record = SeqRecord(Seq("".join(current_seq)), id=current_id, description=current_description)
                    seq_records.append(record)
                header_line = line.strip()[1:]  # remove '>'
                current_description = header_line
                match = re.match(r"(\S+), weight: ([\d\.]+)x", header_line)
                if not match:
                    continue
                current_id, weight = match.groups()
                abundances[current_id] = float(weight)
                current_seq = []
            else:
                current_seq.append(line.strip())

        # Don't forget the last one
        if current_id is not None:
            record = SeqRecord(Seq("".join(current_seq)), id=current_id, description=current_description)
            seq_records.append(record)

    return seq_records, abundances


def align_and_assign(reconstructed_seqs, gt_seqs):
    assignments = {}
    aligner = Align.PairwiseAligner()
    aligner.match_score = 1
    for rec in reconstructed_seqs:
        best_match = None
        best_similarity = -1.0
        for gt in gt_seqs:
            score_forward = aligner.score(gt.seq, rec.seq)
            score_reverse = aligner.score(gt.seq, rec.seq.reverse_complement())
            score = max(score_forward, score_reverse)
            if score > best_similarity:
                best_similarity = score
                best_match = gt
        if best_similarity < 0.8 * len(rec.seq):  # Treshold similarity
            print(f"Low alignment score for {rec.description}, best match: {best_match.id}, score: {best_similarity}")
        else:
            assignments[rec.description] = best_match.id
    return assignments


def get_abundances_from_headers(seqs):
    # Extracts weight from header formatted like: >id, weight: 681.0x
    abundances = {}
    for seq in seqs:
        header = seq.description
        match = re.search(r"weight:\s*([\d.]+)x", header)
        if match:
            abundance = float(match.group(1))
        else:
            print('weight not found in fasta file')
            abundance = 1.0  # fallback
        abundances[seq.description] = abundance
    return abundances


def plot_relative_errors(gt_norm, rec_norm, strain, seed, vg_norm=None):
    """
    Plots relative error vs ground truth abundance for a single sample.
    Includes optional VG-flow results.
    """
    def compute_rel_errors(gt_norm, other_norm):
        x_vals, y_vals = [], []
        for desc in gt_norm:
            gt_val = gt_norm[desc]
            rec_val = other_norm.get(desc, 0.0)
            mean_val = (gt_val + rec_val) / 2
            rel_error = abs(gt_val - rec_val) / mean_val if mean_val > 0 else 0.0
            x_vals.append(gt_val)
            y_vals.append(rel_error)
        return x_vals, y_vals

    rec_x, rec_y = compute_rel_errors(gt_norm, rec_norm)
    vg_x, vg_y = compute_rel_errors(gt_norm, vg_norm) if vg_norm else ([], [])

    plt.figure(figsize=(6, 4))
    plt.scatter(rec_x, rec_y, color='steelblue', label='Flowpaths', alpha=0.6)
    if vg_norm:
        plt.scatter(vg_x, vg_y, color='darkorange', label='VG-flow', alpha=0.6)

    plt.xlabel("Ground truth abundance")
    plt.ylabel("Relative error")
    plt.title(f"HCV {strain}-strain mixture")
    plt.ylim(0, 0.025)
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(f"./flow_assembly/run_flowpaths/output/abundance_estimations_visualizations/relative_error_strain{strain}_seed{seed}.png")
    plt.close()


def plot_relative_errors_all_seeds(strain_data, output_dir="./flow_assembly/run_flowpaths/output/abundance_estimations_visualizations"):
    """
    Plots relative error vs ground truth abundance for all seeds of each strain.
    """
    os.makedirs(output_dir, exist_ok=True)

    # y_limit = max(strain_data.values()['rec_y'], strain_data.values()['vg_y'])

    for strain, data in strain_data.items():
        plt.figure(figsize=(6, 4))
        plt.scatter(data['rec_x'], data['rec_y'], color='steelblue', label='Flowpaths', alpha=0.6)
        if data['vg_x']:
            plt.scatter(data['vg_x'], data['vg_y'], color='darkorange', label='VG-flow', alpha=0.6)

        plt.xlabel("Ground truth abundance")
        plt.ylabel("Relative error")
        plt.title(f"HCV {strain}-strain mixture (all seeds)")
        plt.ylim(0, 0.025)
        plt.legend()
        plt.grid(True)
        plt.tight_layout()
        output_path = os.path.join(output_dir, f"relative_error_strain{strain}_ALL_SEEDS.png")
        plt.savefig(output_path)
        plt.close()


def plot_l1_norm_boxplots(samples, l1_errors_per_method, output_path="boxplot_l1_norms.png"):
    """
    Creates a boxplot of L1 norm errors per number of strains.
    
    Parameters:
    - samples: list of (strain, seed) tuples
    - l1_errors_per_method: dict of the form:
        {
            'flowpaths': {(strain, seed): L1_error, ...},
            'vgflow': {(strain, seed): L1_error, ...}
        }
    - output_path: where to save the plot
    """

    grouped_errors = defaultdict(lambda: {'flowpaths': [], 'vgflow': []})

    for (strain, seed) in samples:
        l1_flow = l1_errors_per_method['flowpaths'].get((strain, seed))
        l1_vg = l1_errors_per_method['vgflow'].get((strain, seed))
        if l1_flow is not None and l1_vg is not None:
            grouped_errors[strain]['flowpaths'].append(l1_flow)
            grouped_errors[strain]['vgflow'].append(l1_vg)

    # Prepare boxplot data
    labels = []
    box_data = []

    for strain in sorted(grouped_errors.keys()):
        fp = grouped_errors[strain]['flowpaths']
        vg = grouped_errors[strain]['vgflow']
        if fp and vg:
            labels.extend([f"{strain} strains\nFlowpaths", f"{strain} strains\nVG-flow"])
            box_data.extend([fp, vg])

    # Plot
    plt.figure(figsize=(10, 6))
    bplot = plt.boxplot(box_data, patch_artist=True, labels=labels, widths=0.6)

    # Color Flowpaths and VG-flow differently
    colors = ['steelblue', 'darkorange']
    for i, patch in enumerate(bplot['boxes']):
        patch.set_facecolor(colors[i % 2])

    plt.ylabel("L1 norm (|est. - ground truth|)")
    plt.title("L1 norm per method grouped by number of strains")
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()


def plot_abundance_error_boxplots_per_strain_count(strain_data, output_path):
    """
    Create a boxplot of relative errors grouped by number of ground truth strains,
    with Flowpaths and VG-flow comparison, and Wilcoxon test annotations.
    """

    grouped_data = defaultdict(lambda: {'Flowpaths': [], 'VG-flow': []})
    for strain, data in strain_data.items():
        num_strains = len(data['gt_abunds'])
        grouped_data[num_strains]['Flowpaths'].extend(data['rec_y'])
        grouped_data[num_strains]['VG-flow'].extend(data['vg_y'])

    fig, ax = plt.subplots(figsize=(10, 6))
    box_data = []
    box_positions = []
    xtick_positions = []
    xtick_labels = []
    significance_annotations = []

    offset = 0.2
    colors = ['skyblue', 'salmon']
    index = 0

    for i, strain_count in enumerate(sorted(grouped_data.keys())):
        flowpaths_errors = grouped_data[strain_count]['Flowpaths']
        vgflow_errors = grouped_data[strain_count]['VG-flow']

        if not flowpaths_errors or not vgflow_errors:
            continue

        # Add data
        box_data.append(flowpaths_errors)
        box_positions.append(i - offset)
        box_data.append(vgflow_errors)
        box_positions.append(i + offset)

        # X-tick at the center between two bars
        xtick_positions.append(i)
        xtick_labels.append(f"{strain_count} strains")

        # Wilcoxon test
        min_len = min(len(flowpaths_errors), len(vgflow_errors))
        stat, p_value = wilcoxon(flowpaths_errors[:min_len], vgflow_errors[:min_len])
        if p_value < 0.001:
            annotation = '***'
        elif p_value < 0.01:
            annotation = '**'
        elif p_value < 0.05:
            annotation = '*'
        else:
            annotation = 'ns'

        # Get y position for annotation
        y_max = max(max(flowpaths_errors), max(vgflow_errors))
        significance_annotations.append((i, y_max + 0.01, annotation))

    # Plot
    bplot = ax.boxplot(box_data, positions=box_positions, widths=0.35, patch_artist=True)
    for patch, color in zip(bplot['boxes'], colors * len(grouped_data)):
        patch.set_facecolor(color)

    # X-axis
    ax.set_xticks(xtick_positions)
    ax.set_xticklabels(xtick_labels)
    ax.set_ylabel("Absolute error (estimated - true)")
    ax.set_title("Abundance estimation error per strain count (Wilcoxon significance)")
    ax.grid(axis='y', linestyle='--', alpha=0.7)

    # Legend
    handles = [plt.Line2D([0], [0], color=color, lw=10) for color in colors]
    ax.legend(handles, ['Flowpaths', 'VG-flow'], loc='upper right')

    # Significance annotations
    for x, y, annotation in significance_annotations:
        ax.text(x, y, annotation, ha='center', va='bottom', fontsize=12, fontweight='bold')

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()



strain_data = defaultdict(lambda: {'rec_x': [], 'rec_y': [], 'vg_x': [], 'vg_y': []})
l1_errors = {'flowpaths': {}, 'vgflow': {}}

for strain, seed in samples:
    print(f"\nProcessing strain {strain}, seed {seed}")
    
    # Full paths to files
    solution_path = os.path.join(directory_solutions, f"solution_{strain}_{seed}")
    gt_folder = os.path.join(directory_ground_truth, f"{strain}-strain_seed_{seed}")
    vg_flow_solution_path = os.path.join(directory_vg_flow, f"haps_{strain}_{seed}.fasta")
    
    try:
        reconstructed = load_fasta(solution_path)
        ground_truth, gt_abundances = load_ground_truth_abundance(gt_folder)
        vg_sol, vg_abundances = read_fasta_vg_flow(vg_flow_solution_path)
    except FileNotFoundError as e:
        print(f"{e}, skipping.")
        continue

    if not reconstructed or not ground_truth:
        print("Empty file(s), skipping...")
        continue

    rec_abundances = get_abundances_from_headers(reconstructed)

    assignments = align_and_assign(reconstructed_seqs=reconstructed, gt_seqs=ground_truth)
    assignments_vg = align_and_assign(reconstructed_seqs=vg_sol, gt_seqs=ground_truth)
    print("Assignment counts:", Counter(assignments.values()))


    rec_sum_by_gt = defaultdict(float)
    for rec_desc, gt_desc in assignments.items():
        rec_sum_by_gt[gt_desc] += rec_abundances.get(rec_desc, 0)

    vg_sum_by_gt = defaultdict(float)
    for vg_desc, gt_desc in assignments_vg.items():
        vg_sum_by_gt[gt_desc] += vg_abundances.get(vg_desc, 0)

    assigned_gt_descs = set(rec_sum_by_gt.keys()).union(vg_sum_by_gt.keys())

    # Normalize GT abundances
    gt_total = sum(gt_abundances.get(desc, 0) for desc in assigned_gt_descs)
    if gt_total == 0:
        print("Ground truth total abundance is 0, skipping normalization...")
        continue

    gt_norm = {desc: gt_abundances.get(desc, 0) / gt_total for desc in assigned_gt_descs}

    # Normalize reconstructed abundances
    rec_total = sum(rec_sum_by_gt.values())
    rec_norm = {desc: rec_sum_by_gt[desc] / rec_total for desc in assigned_gt_descs}

    vg_total = sum(vg_sum_by_gt.values())
    vg_norm = {desc: vg_sum_by_gt[desc] / vg_total for desc in assigned_gt_descs}

    def compute_rel_errors(gt_norm, other_norm):
        x_vals, y_vals = [], []
        for desc in gt_norm:
            gt_val = gt_norm[desc]
            rec_val = other_norm.get(desc, 0.0)
            mean_val = (gt_val + rec_val) / 2
            rel_error = abs(gt_val - rec_val) / mean_val if mean_val > 0 else 0.0
            x_vals.append(gt_val)
            y_vals.append(rel_error)
        return x_vals, y_vals

    rec_x, rec_y = compute_rel_errors(gt_norm, rec_norm)
    vg_x, vg_y = compute_rel_errors(gt_norm, vg_norm)

    # Add to per-strain data for aggregate plot
    strain_data[strain]['rec_x'].extend(rec_x)
    strain_data[strain]['rec_y'].extend(rec_y)
    strain_data[strain]['vg_x'].extend(vg_x)
    strain_data[strain]['vg_y'].extend(vg_y)

    # Collect L1 norms
    def compute_l1_norm(gt_norm, est_norm):
        return sum(abs(gt_norm[desc] - est_norm.get(desc, 0.0)) for desc in gt_norm)
    
    # Store data for all-seed plotting and statistical analysis
    strain_data[strain]['gt_abunds'] = [gt_norm[desc] for desc in assigned_gt_descs]
    strain_data[strain]['rec_x'].extend(list(gt_norm.values()))
    strain_data[strain]['rec_y'].extend([abs(gt_norm[desc] - rec_norm.get(desc, 0)) / ((gt_norm[desc] + rec_norm.get(desc, 0)) / 2) for desc in assigned_gt_descs])
    strain_data[strain]['vg_x'].extend(list(gt_norm.values()))
    strain_data[strain]['vg_y'].extend([abs(gt_norm[desc] - vg_norm.get(desc, 0)) / ((gt_norm[desc] + vg_norm.get(desc, 0)) / 2) for desc in assigned_gt_descs])


    l1_errors.setdefault('flowpaths', {})[(strain, seed)] = compute_l1_norm(gt_norm, rec_norm)
    l1_errors.setdefault('vgflow', {})[(strain, seed)] = compute_l1_norm(gt_norm, vg_norm)



    # plot_relative_errors(gt_norm, rec_norm, strain, seed, vg_norm=vg_norm)

    print("GT_description\t\tTrue_freq\tRecon_freq")
    for desc in sorted(assigned_gt_descs):
        print(f"{desc}\t{gt_norm.get(desc, 0.0):.4f}\t\t{vg_norm.get(desc, 0.0):.4f}")

# To store per-point absolute errors for boxplot
per_point_errors = dict()

for strain, seed in samples:
    key = (strain, seed)
    if not os.path.isfile(os.path.join(directory_solutions, f"solution_{strain}_{seed}")):
        continue

    # Make sure normalized abundances exist from previous processing
    if not 'gt_norm' in locals() or not 'rec_norm' in locals() or not 'vg_norm' in locals():
        continue

    abs_err_rec = [abs(gt_norm[desc] - rec_norm.get(desc, 0.0)) for desc in gt_norm]
    abs_err_vg = [abs(gt_norm[desc] - vg_norm.get(desc, 0.0)) for desc in gt_norm]

    per_point_errors[key] = {
        'flowpaths': abs_err_rec,
        'vgflow': abs_err_vg
    }


# Compute Wilcoxon signed-rank test on relative errors across all strains/seeds
all_flowpaths_errors = []
all_vgflow_errors = []

for strain, data in strain_data.items():
    # Make sure to only use paired values where both methods have a corresponding GT abundance
    paired_errors = list(zip(data['rec_y'], data['vg_y']))
    for rec_err, vg_err in paired_errors:
        all_flowpaths_errors.append(rec_err)
        all_vgflow_errors.append(vg_err)

# Perform Wilcoxon test
if len(all_flowpaths_errors) != len(all_vgflow_errors):
    print("Mismatch in number of errors between Flowpaths and VG-flow. Cannot perform Wilcoxon test.")
else:
    stat, p_value = wilcoxon(all_flowpaths_errors, all_vgflow_errors)
    print(f"\nWilcoxon signed-rank test result:")
    print(f"Statistic: {stat:.4f}, p-value: {p_value:.4e}")
    if p_value < 0.05:
        print("Result is statistically significant (p < 0.05).")
    else:
        print("Result is not statistically significant (p >= 0.05).")


# Now plot the boxplot
plot_abundance_error_boxplots_per_strain_count(strain_data, output_path="./flow_assembly/run_flowpaths/output/abundance_estimations_visualizations/boxplot_l1_norms_with_significance_nodes_removed.png")


# plot_relative_errors_all_seeds(strain_data)
# plot_l1_norm_boxplots(samples, l1_errors, output_path="./flow_assembly/run_flowpaths/output/abundance_estimations_visualizations/l1_norm_boxplot.png")


       


