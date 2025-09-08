import os
import re
from collections import defaultdict
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import os

# Helper functions
def load_fasta(path):
    with open(path, "r") as f:
        return list(SeqIO.parse(f, "fasta"))

def load_fasta_with_weight(path):
    records = []
    weights = {}
    for record in SeqIO.parse(path, "fasta"):
        match = re.search(r"(\d+)x", record.description)
        if not match:
            raise ValueError(f"No weight in header: {record.description}")
        weights[record.description] = int(match.group(1))
        records.append(record)
    return records, weights

def load_gt_abundances(folder):
    gt_file = os.path.join(folder, "ground_truth_abundances.txt")
    seqs, abund = [], {}
    with open(gt_file) as f:
        for line in f:
            if line.startswith(">"):
                if 'current_id' in locals():
                    seqs.append(SeqRecord(Seq("".join(seq_buf)), id=current_id, description=desc))
                match = re.match(r"(\S+), weight: ([\d\.]+)x", line.strip()[1:])
                if match:
                    current_id, w = match.groups()
                    abund[current_id] = float(w)
                    desc = line.strip()[1:]
                    seq_buf = []
            else:
                seq_buf.append(line.strip())
        if 'current_id' in locals():
            seqs.append(SeqRecord(Seq("".join(seq_buf)), id=current_id, description=desc))
    return seqs, abund

def assign_recon_to_gt(recon, gt):
    aligner = Align.PairwiseAligner()
    aligner.match_score = 1
    mapping = {}
    for r in recon:
        best, best_score = None, -1
        for g in gt:
            score = max(aligner.score(g.seq, r.seq), aligner.score(g.seq, r.seq.reverse_complement()))
            if score > best_score:
                best, best_score = g, score
        if best_score >= 0.8 * len(r.seq):
            mapping[r.description] = best.id
    return mapping

def extract_abundances(seqs):
    abund = {}
    for s in seqs:
        match = re.search(r"weight:\s*([\d.]+)x", s.description)
        abund[s.description] = float(match.group(1)) if match else 1.0
    return abund


def create_dataframe(samples, dir_solutions, dir_vg, dir_virusvg, dir_gt, csv_path): 
    data = []

    for strain, seed in samples:
        print(f"Processing {strain}, seed {seed}")
        sol_path = os.path.join(dir_solutions, f"solution_{strain}_{seed}")
        vg_path = os.path.join(dir_vg, f"haps_{strain}_{seed}.fasta")
        virusvg_path = os.path.join(dir_virusvg, f"{strain}-strain_seed_{seed}", f"haps_{strain}_{seed}.fasta")
        gt_path = os.path.join(dir_gt, f"{strain}-strain_seed_{seed}")

        try:
            rec_seqs = load_fasta(sol_path)
            vg_seqs, vg_weights = load_fasta_with_weight(vg_path)
            virus_seqs, virus_weights = load_fasta_with_weight(virusvg_path)
            gt_seqs, gt_abund = load_gt_abundances(gt_path)
        except Exception as e:
            print(f"Error loading files for {strain}-{seed}: {e}")
            continue

        # Assignments
        rec_assign = assign_recon_to_gt(rec_seqs, gt_seqs)
        vg_assign = assign_recon_to_gt(vg_seqs, gt_seqs)
        virus_assign = assign_recon_to_gt(virus_seqs, gt_seqs)

        # Extract raw abundances
        rec_abund = extract_abundances(rec_seqs)
        rec_total = defaultdict(float)
        for desc, gt_id in rec_assign.items():
            rec_total[gt_id] += rec_abund.get(desc, 0)

        vg_total = defaultdict(float)
        for desc, gt_id in vg_assign.items():
            vg_total[gt_id] += vg_weights.get(desc, 0)

        virus_total = defaultdict(float)
        for desc, gt_id in virus_assign.items():
            virus_total[gt_id] += virus_weights.get(desc, 0)

        # Union of all assigned ground truth IDs
        assigned = set(rec_total) | set(vg_total) | set(virus_total)
        gt_sum = sum(gt_abund.get(d, 0) for d in assigned)
        if gt_sum == 0:
            print(f"GT abundance sum is 0 for {strain}-{seed}. Skipping...")
            continue

        # Normalize ground truth
        gt_norm = {d: gt_abund.get(d, 0)/gt_sum for d in assigned}

        # Normalize estimates for each method
        rec_sum = sum(rec_total.values())
        vg_sum = sum(vg_total.values())
        virus_sum = sum(virus_total.values())

        rec_norm = {d: rec_total[d]/rec_sum if rec_sum > 0 else 0 for d in assigned}
        vg_norm = {d: vg_total[d]/vg_sum if vg_sum > 0 else 0 for d in assigned}
        virus_norm = {d: virus_total[d]/virus_sum if virus_sum > 0 else 0 for d in assigned}

        for gt_id in assigned:
            data.append({
                "strain": strain,
                "seed": seed,
                "gt_id": gt_id,
                "ground_truth_abundance": gt_norm.get(gt_id, 0),
                "estimated_abundance": rec_norm.get(gt_id, 0),
                "method": "MFD"
            })
            data.append({
                "strain": strain,
                "seed": seed,
                "gt_id": gt_id,
                "ground_truth_abundance": gt_norm.get(gt_id, 0),
                "estimated_abundance": vg_norm.get(gt_id, 0),
                "method": "vg-flow"
            })
            data.append({
                "strain": strain,
                "seed": seed,
                "gt_id": gt_id,
                "ground_truth_abundance": gt_norm.get(gt_id, 0),
                "estimated_abundance": virus_norm.get(gt_id, 0),
                "method": "virus-vg"
            })

    # Create and export dataframe
    df = pd.DataFrame(data)
    df.to_csv(csv_path, index=False)


def plot_abundance_scatter_by_method_and_dataset(df_hiv, df_hcv, output_dir):
    method_colors = {
        "MFD": "salmon",
        "vg-flow": "skyblue",
        "virus-vg": "goldenrod"
    }

    os.makedirs(output_dir, exist_ok=True)

    combined_df = pd.concat([df_hiv, df_hcv], ignore_index=True)
    methods = ["MFD", "vg-flow", "virus-vg"]
    datasets = ["HIV", "HCV"]

    for dataset in datasets:
        for method in methods:
            df_sub = combined_df[(combined_df["method"] == method) & (combined_df["dataset"] == dataset)]

            if df_sub.empty:
                print(f"No data for {method} on {dataset}. Skipping.")
                continue

            # Compute R²
            y_true = df_sub["ground_truth_abundance"]
            y_pred = df_sub["estimated_abundance"]
            ss_res = ((y_true - y_pred) ** 2).sum()
            ss_tot = ((y_true - y_true.mean()) ** 2).sum()
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')

            # Plot
            plt.figure(figsize=(6, 6))
            plt.plot([0, 1], [0, 1], 'k--', label="Perfect estimation", alpha=0.6)

            sns.scatterplot(
                data=df_sub,
                x="ground_truth_abundance",
                y="estimated_abundance",
                color=method_colors[method],
                alpha=0.9
            )

            plt.xlim(0, 1)
            plt.ylim(0, 1)
            plt.xlabel("Ground truth abundance", fontsize=12)
            plt.ylabel("Estimated abundance", fontsize=12)
            plt.title(f"{method} on {dataset}", fontsize=14)
            plt.grid(True)

            # Annotate R²
            plt.text(0.05, 0.95, f"$R^2$ = {r2:.3f}", transform=plt.gca().transAxes,
                     fontsize=12, verticalalignment='top')

            plt.tight_layout()
            output_path = os.path.join(output_dir, f"abundance_scatter_{method}_{dataset}.png")
            plt.savefig(output_path)
            plt.close()
            print(f"Saved: {output_path}")


def plot_abundance_subplots(df_hiv, df_hcv, output_path):

    method_colors = {
        "MFD": "salmon",
        "vg-flow": "skyblue",
        "virus-vg": "goldenrod"
    }

    methods = ["MFD", "vg-flow", "virus-vg"]
    datasets = ["HIV", "HCV"]

    # Merge data
    df_combined = pd.concat([df_hiv, df_hcv], ignore_index=True)

    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(18, 10), sharex=True, sharey=True)

    for row_idx, dataset in enumerate(datasets):
        for col_idx, method in enumerate(methods):
            ax = axes[row_idx, col_idx]
            df_sub = df_combined[
                (df_combined["method"] == method) &
                (df_combined["dataset"] == dataset)
            ]

            if df_sub.empty:
                ax.set_visible(False)
                continue

            # Compute R²
            y_true = df_sub["ground_truth_abundance"]
            y_pred = df_sub["estimated_abundance"]
            ss_res = ((y_true - y_pred) ** 2).sum()
            ss_tot = ((y_true - y_true.mean()) ** 2).sum()
            r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')

            # Plot
            ax.plot([0, 1], [0, 1], 'k--', alpha=0.6)
            sns.scatterplot(
                data=df_sub,
                x="ground_truth_abundance",
                y="estimated_abundance",
                color=method_colors[method],
                ax=ax,
                alpha=0.9
            )

            ax.set_xlim(0, 1)
            ax.set_ylim(0, 1)
            ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
            ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
            ax.tick_params(axis='both', labelsize=14)
            ax.set_title(f"{method} ({dataset})", fontsize=15)
            ax.grid(True)

            # Axis labels
            if row_idx == 1:
                ax.set_xlabel("Ground truth abundance", fontsize=14)
            else:
                ax.set_xlabel("")
            if col_idx == 0:
                ax.set_ylabel("Estimated abundance", fontsize=14)
            else:
                ax.set_ylabel("")

            # R² annotation
            ax.text(0.05, 0.95, f"$R^2$ = {r2:.3f}", transform=ax.transAxes,
                    fontsize=13, verticalalignment='top')

    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()
    print(f"Saved combined plot to: {output_path}")


# ---------------------------
# HIV dataset
# ---------------------------
samples_hiv = [(2, 2), (2, 3), (2,4), (2,5), (2,6), (2,7), (2,8), (2,9),
               (3, 10), (3, 13), (3,15), (3,16), (3,18), (3,19),
               (4, 29)]

dir_solutions_hiv = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/min_cover/HIV_8000_dagify"
dir_gt_hiv = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition/Data/simulation_data/simulation_data/art_reads/150bp/HIV_8000"
dir_vg_hiv = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition/output_VG_flow/HIV_dagify"
dir_virusvg_hiv = "/tudelft.net/staff-umbrella/FlowDecomposition/virus-vg/output/HIV"
csv_hiv = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/abundance_estimations_visualizations/supplement/abundance_HIV.csv"

if not os.path.exists(csv_hiv):
    create_dataframe(samples_hiv, dir_solutions_hiv, dir_vg_hiv, dir_virusvg_hiv, dir_gt_hiv, csv_hiv)
df_hiv = pd.read_csv(csv_hiv)
df_hiv["dataset"] = "HIV"

# ---------------------------
# HCV dataset
# ---------------------------
samples_hcv = [(3, 0), (3, 1), (3, 2), (3, 4), (3, 5), (3, 6), (3, 7), (3, 9),
               (4, 0), (4, 5), (4, 6), (4, 7), (4, 8), (5, 1), (5, 6),
               (6, 1), (6, 4), (6, 5), (7, 6), (7, 8)]

dir_solutions_hcv = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/min_cover_strict_constraints/nodes_removed"
dir_gt_hcv = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition/Data/simulation_data/simulation_data/art_reads/150bp/random_abundances"
dir_vg_hcv = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition/output_VG_flow/HCV_nodes_removed"
dir_virusvg_hcv = "/tudelft.net/staff-umbrella/FlowDecomposition/virus-vg/output/HCV"
csv_hcv = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/abundance_estimations_visualizations/supplement/abundance_HCV.csv"

if not os.path.exists(csv_hcv):
    create_dataframe(samples_hcv, dir_solutions_hcv, dir_vg_hcv, dir_virusvg_hcv, dir_gt_hcv, csv_hcv)
df_hcv = pd.read_csv(csv_hcv)
df_hcv["dataset"] = "HCV"

output_dir = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/abundance_estimations_visualizations/supplement/abundance_subplots_combined.png"
plot_abundance_subplots(df_hiv, df_hcv, output_dir)



# # Manual R² computation
# y_true = df["ground_truth_abundance"]
# y_pred = df["estimated_abundance"]
# ss_res = ((y_true - y_pred) ** 2).sum()
# ss_tot = ((y_true - y_true.mean()) ** 2).sum()
# r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')

# # Plot
# plt.figure(figsize=(8, 8))
# plt.plot([0, 0.8], [0, 0.8], 'k-', label="Perfect estimation", alpha=0.5)
# sns.scatterplot(data=df, x="ground_truth_abundance", y="estimated_abundance", hue="method", alpha=0.8)
# plt.ylim(bottom=0)
# plt.xlim(left=0)
# plt.xlabel("Ground truth abundance", fontsize=14)
# plt.ylabel("Estimated abundance (MFD)", fontsize=14)
# plt.title("Abundance: ground truth vs MFD estimate", fontsize=15)
# plt.xticks(fontsize=14)
# plt.yticks(fontsize=14)
# plt.legend()
# plt.grid(True)

# # Annotate R² in the top-left corner
# plt.text(0.05, 0.95, f"$R^2$ = {r2:.3f}", transform=plt.gca().transAxes, fontsize=14, verticalalignment='top')

# plt.tight_layout()

# # Save to file
# output_path = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/abundance_estimations_visualizations/supplement/abundance_estimates_MFD_with_line.png"
# plt.savefig(output_path)
# plt.close()
# print(f"Plot saved to: {output_path}")

