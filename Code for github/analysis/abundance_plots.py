import os
import re
from collections import defaultdict
from Bio import SeqIO, Align
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from scipy.stats import wilcoxon
from itertools import combinations

# In seqwish container

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


def compute_l1(gt_norm, est_norm):
    return sum(abs(gt_norm[desc] - est_norm.get(desc, 0)) for desc in gt_norm)


def create_l1_errors_dataframe(samples, dir_solutions, dir_vg, dir_virusvg, dir_gt, csv_path):
    # Main data structure
    l1_errors = defaultdict(dict)

    # Main loop
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

        # Align and assign
        rec_assign = assign_recon_to_gt(rec_seqs, gt_seqs)
        vg_assign = assign_recon_to_gt(vg_seqs, gt_seqs)
        virus_assign = assign_recon_to_gt(virus_seqs, gt_seqs)

        # Aggregate and normalize
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

        assigned = set(rec_total) | set(vg_total) | set(virus_total)
        gt_sum = sum(gt_abund.get(d, 0) for d in assigned)
        if gt_sum == 0:
            print("GT abundance sum is 0. Skipping...")
            continue

        gt_norm = {d: gt_abund.get(d, 0)/gt_sum for d in assigned}
        rec_norm = {d: rec_total[d]/sum(rec_total.values()) for d in assigned}
        vg_norm = {d: vg_total[d]/sum(vg_total.values()) for d in assigned}
        virus_norm = {d: virus_total[d]/sum(virus_total.values()) for d in assigned}

        # Store L1 errors
        l1_errors['MFD'][(strain, seed)] = compute_l1(gt_norm, rec_norm)
        l1_errors['vg-flow'][(strain, seed)] = compute_l1(gt_norm, vg_norm)
        l1_errors['virus-vg'][(strain, seed)] = compute_l1(gt_norm, virus_norm)
    # Convert to DataFrame
    rows = []
    for method, errors in l1_errors.items():
        for (strain, seed), error in errors.items():
            rows.append({
                "strain_count": strain,
                "seed": seed,
                "method": method,
                "error": error
            })
    df_l1 = pd.DataFrame(rows)

    # Save to CSV
    df_l1.to_csv(csv_path, index=False)

    # Output example
    print("\nSample L1 errors (first 5):")
    for method in l1_errors:
        print(f"{method}:")
        for k, v in list(l1_errors[method].items())[:5]:
            print(f"  {k}: {v:.4f}")

def plot_l1_error_per_sample(l1_errors_df, output_path):
    """
    Plots L1 error per sample:
    - X-axis: number of strains (aligned samples)
    - Y-axis: L1 error
    - Methods shown with side-by-side points (no jitter)
    - Custom colors
    """
    # Method-specific offsets and colors
    method_offsets = {
        "MFD": -0.2,
        "vg-flow": 0,
        "virus-vg": 0.2
    }
    method_colors = {
        "MFD": "salmon",
        "vg-flow": "skyblue",
        "virus-vg": "goldenrod"  
    }

    plt.figure(figsize=(10, 6))

    for method in l1_errors_df["method"].unique():
        df_m = l1_errors_df[l1_errors_df["method"] == method].copy()
        df_m["x"] = df_m["strain_count"] + method_offsets[method]
        plt.scatter(
            df_m["x"],
            df_m["error"],
            label=method,
            edgecolors=method_colors[method],   # Colored border
            facecolors='none',                  # Transparent inside
            linewidths=1.5,
            s=100,
            marker='o'
        )
        
    # Count datapoints per strain
    strain_counts = l1_errors_df.groupby("strain_count")["seed"].nunique().to_dict()

    strain_levels = sorted(l1_errors_df["strain_count"].unique())
    xtick_labels = [
        f"{s}\n(n={strain_counts.get(s, 0)})" for s in strain_levels
    ]
    plt.xticks(strain_levels, xtick_labels, fontsize=12)
    plt.yscale("log")
    plt.xlabel("Number of haplotypes", fontsize=14)
    plt.ylabel("Abundance estimation error (L1 norm)", fontsize=14)
    plt.title("HCV per-sample L1 error by method", fontsize=15)
    plt.xticks(sorted(l1_errors_df["strain_count"].unique()))
    plt.grid(axis="y", linestyle="--", alpha=0.5)
    plt.ylim(bottom=0)
    plt.ylim(top=0.6)
    plt.legend(title="Method", bbox_to_anchor=(1.05, 1), loc="upper left")
    plt.tight_layout()
    plt.savefig(output_path)
    plt.close()

def plot_l1_error_per_sample_side_by_side(df_combined, output_path):
    """
    Plots HIV and HCV L1 errors side by side with shared legend.
    """
    method_offsets = {
        "MFD": -0.2,
        "vg-flow": 0.0,
        "virus-vg": 0.2
    }
    method_colors = {
        "MFD": "salmon",
        "vg-flow": "skyblue",
        "virus-vg": "goldenrod"
    }

    datasets = ["HIV", "HCV"]
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)

    for ax, dataset in zip(axes, datasets):
        df = df_combined[df_combined["dataset"] == dataset]
        strain_counts = df.groupby("strain_count")["seed"].nunique().to_dict()
        strain_levels = sorted(df["strain_count"].unique())
        xtick_labels = [
            f"{s}\n(n={strain_counts.get(s, 0)})" for s in strain_levels
        ]

        for method in df["method"].unique():
            df_m = df[df["method"] == method].copy()
            df_m["x"] = df_m["strain_count"] + method_offsets[method]
            ax.scatter(
                df_m["x"],
                df_m["error"],
                label=method,
                edgecolors=method_colors[method],
                facecolors='none',
                linewidths=1.5,
                s=100,
                marker='o'
            )

        ax.set_title(f"{dataset}", fontsize=15)
        ax.set_xlabel("Number of haplotypes", fontsize=13)
        ax.set_xticks(strain_levels)
        ax.set_xticklabels(xtick_labels, fontsize=11)
        ax.grid(axis="y", linestyle="--", alpha=0.5)

    axes[0].set_ylabel("Abundance estimation error (L1 norm)", fontsize=13)
    axes[0].set_yscale("log")
    axes[0].set_ylim(bottom=0, top=0.6)
    axes[1].legend(title="Method", bbox_to_anchor=(1.05, 1), loc="upper left")


    plt.tight_layout()
    plt.subplots_adjust(right=0.85)
    plt.savefig(output_path)
    plt.close()

from scipy.stats import wilcoxon
from itertools import combinations


def run_wilcoxon_tests(df_combined):
    """
    Performs pairwise two-sided Wilcoxon signed-rank tests on L1 errors
    across methods, separately for each dataset (HIV, HCV).
    
    Assumes df_combined contains:
      - 'dataset'
      - 'strain_count'
      - 'seed'
      - 'method'
      - 'error'
    
    Returns a DataFrame with test results.
    """
    results = []

    for dataset in df_combined["dataset"].unique():
        df_dataset = df_combined[df_combined["dataset"] == dataset]
        methods = df_dataset["method"].unique()

        for m1, m2 in combinations(methods, 2):
            df_m1 = df_dataset[df_dataset["method"] == m1]
            df_m2 = df_dataset[df_dataset["method"] == m2]

            # Merge on (strain_count, seed) to get paired samples
            merged = pd.merge(
                df_m1,
                df_m2,
                on=["strain_count", "seed"],
                suffixes=(f"_{m1}", f"_{m2}")
            )

            if len(merged) < 5:
                print(f"Skipping {dataset}: too few paired samples between {m1} and {m2}")
                continue

            if m1 == "MFD":
                alternative = "less"
            elif m2 == "MFD":
                alternative = "greater"
            else:
                alternative = "two-sided"

            # Perform two-sided Wilcoxon test
            stat, p = wilcoxon(
                merged[f"error_{m1}"],
                merged[f"error_{m2}"],
                alternative=alternative
            )

            results.append({
                "dataset": dataset,
                "method_1": m1,
                "method_2": m2,
                "n_samples": len(merged),
                "statistic": stat,
                "p_value": p
            })

    return pd.DataFrame(results)



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
csv_hiv = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/abundance_estimations_visualizations/sample_error/l1_errors_HIV.csv"

if not os.path.exists(csv_hiv):
    create_l1_errors_dataframe(samples_hiv, dir_solutions_hiv, dir_vg_hiv, dir_virusvg_hiv, dir_gt_hiv, csv_hiv)
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
csv_hcv = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/abundance_estimations_visualizations/sample_error/l1_errors_HCV.csv"

if not os.path.exists(csv_hcv):
    create_l1_errors_dataframe(samples_hcv, dir_solutions_hcv, dir_vg_hcv, dir_virusvg_hcv, dir_gt_hcv, csv_hcv)
df_hcv = pd.read_csv(csv_hcv)
df_hcv["dataset"] = "HCV"



# plot_l1_error_per_sample(
#     df_hcv,
#     output_path="/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/abundance_estimations_visualizations/sample_error/l1_error_per_sample_HCV_log_scale.png"
# )

# Combine and plot
# df_combined = pd.concat([df_hiv, df_hcv])
# plot_l1_error_per_sample_side_by_side(
#     df_combined,
#     output_path="/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/abundance_estimations_visualizations/sample_error/l1_error_HIV_HCV_side_by_side.png"
# )

df_combined = pd.concat([df_hiv, df_hcv])
df_stats = run_wilcoxon_tests(df_combined)
print(df_stats)
