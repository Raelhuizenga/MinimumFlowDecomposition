import os
import re
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import wilcoxon

BASE_DIR = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/quast_results"

METHOD_DIRS_HCV = {
    "MFD": f"{BASE_DIR}/HCV_filtered_flowpaths_constraints",
    "VG-flow": f"{BASE_DIR}/HCV_vg-flow",
    "Virus-VG": f"{BASE_DIR}/HCV_virus-vg",
    "SAVAGE": f"{BASE_DIR}/savage_HCV"
}

METHOD_DIRS_HIV = {
    "MFD": f"{BASE_DIR}/HIV_MFD",
    "VG-flow": f"{BASE_DIR}/HIV_vg-flow",
    "Virus-VG": f"{BASE_DIR}/HIV_virus-vg",
    "SAVAGE": f"{BASE_DIR}/savage_HIV"
}

included_combinations_HCV = [
    (3,0), (3,1), (3,3), (3,4), (3,5), (3,6), (3,7), (3,9),
    (4,0), (4,5), (4,6), (4,7), (4,8),
    (5,1), (5,6),
    (6,1), (6,4), (6,5),
    (7,6), (7,8)
]

included_combinations_HIV = [
    (2,2), (2,3), (2,4), (2,5), (2,6), (2,7), (2,8), (2,9),
    (3,10), (3,13), (3,15), (3,16), (3,18), (3,19),
    (4,20), (4,29),
]

FIELDS_ERROR = {
    "# N's per 100 kbp": float,
    "# mismatches per 100 kbp": float,
    "# indels per 100 kbp": float,
}

FIELDS_GENOME_FRAC = {
    "Genome fraction (%)": float,
}

# Set the metric here: either "Error rate" or "Genome fraction (%)"
METRIC = "Genome fraction (%)"  # Change to "Error rate" if you want the original behavior
# METRIC = "Error rate"

def extract_metric(base_dir, included_combinations, metric=METRIC):
    data = []
    if not os.path.isdir(base_dir):
        print(f"Warning: base dir does not exist: {base_dir}")
        return data

    for subdir in sorted(os.listdir(base_dir)):
        if not subdir.startswith("strain_") or "_seed_" not in subdir:
            continue
        m = re.match(r"strain_(\d+)_seed_(\d+)", subdir)
        if not m:
            continue
        s, seed = int(m.group(1)), int(m.group(2))
        if (s, seed) not in included_combinations:
            continue

        report_path = os.path.join(base_dir, subdir, "report.txt")
        if not os.path.isfile(report_path):
            continue

        stats = {}
        with open(report_path, "r") as fh:
            for line in fh:
                line = line.strip()
                if metric == "Error rate":
                    for field, dtype in FIELDS_ERROR.items():
                        if line.startswith(field):
                            token = line.split()[-1].replace(",", "")
                            try:
                                stats[field] = dtype(token)
                            except Exception:
                                mnum = re.search(r"[-+]?\d*\.\d+|\d+", line)
                                if mnum:
                                    stats[field] = dtype(mnum.group(0))
                                else:
                                    stats[field] = None
                elif metric == "Genome fraction (%)":
                    for field, dtype in FIELDS_GENOME_FRAC.items():
                        if line.startswith(field):
                            token = line.split()[-1].replace(",", "")
                            try:
                                stats[field] = dtype(token)
                            except Exception:
                                mnum = re.search(r"[-+]?\d*\.\d+|\d+", line)
                                if mnum:
                                    stats[field] = dtype(mnum.group(0))
                                else:
                                    stats[field] = None

        if metric == "Error rate":
            if all(stats.get(k) is not None for k in FIELDS_ERROR):
                err = (stats["# N's per 100 kbp"]
                       + stats["# mismatches per 100 kbp"]
                       + stats["# indels per 100 kbp"]) / 1000.0
                data.append((s, seed, err))
        elif metric == "Genome fraction (%)":
            val = stats.get("Genome fraction (%)")
            if val is not None:
                data.append((s, seed, val))
    return data


def collect_df(method_dirs, included_combinations, dataset_name, metric=METRIC):
    rows = []
    for method, base_dir in method_dirs.items():
        data = extract_metric(base_dir, included_combinations, metric)
        for s, seed, val in data:
            rows.append({
                "Method": method,
                metric: val,
                "Dataset": dataset_name,
                "Strain": s,
                "Seed": seed,
            })
    return pd.DataFrame(rows)


df_hcv = collect_df(METHOD_DIRS_HCV, included_combinations_HCV, "HCV", metric=METRIC)
df_hiv = collect_df(METHOD_DIRS_HIV, included_combinations_HIV, "HIV", metric=METRIC)
df = pd.concat([df_hcv, df_hiv], ignore_index=True)

df_hcv_wide = df_hcv.pivot_table(
    index=["Strain", "Seed"],
    columns="Method",
    values=METRIC
).reset_index()

df_hiv_wide = df_hiv.pivot_table(
    index=["Strain", "Seed"],
    columns="Method",
    values=METRIC
).reset_index()

# Plotting config
method_order = ["MFD", "VG-flow", "Virus-VG", "SAVAGE"]
method_colors = {
    "MFD": "salmon",
    "VG-flow": "skyblue",
    "Virus-VG": "goldenrod",
    "SAVAGE": "mediumseagreen"
}

fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
for ax, dataset in zip(axes, ["HCV", "HIV"]):
    sub = df[df["Dataset"] == dataset]

    # If plotting genome fraction, exclude SAVAGE
    if METRIC == "Genome fraction (%)":
        sub = sub[sub["Method"] != "SAVAGE"]

    if sub.empty:
        ax.text(0.5, 0.5, "No data", ha="center", va="center")
        ax.set_title(dataset)
        continue

    sns.boxplot(
        data=sub,
        x="Method",
        y=METRIC,  # Use dynamic metric variable for y axis
        order=[m for m in method_order if m != "SAVAGE" or METRIC != "Genome fraction (%)"],
        palette={k: v for k, v in method_colors.items() if k != "SAVAGE" or METRIC != "Genome fraction (%)"},
        showmeans=True,
        meanprops={"marker": "o", "markerfacecolor": "black", "markeredgecolor": "black"},
        ax=ax
    )
    
    ax.grid(True, linestyle='--', alpha=0.5)
    
    # Set yscale and labels accordingly
    if METRIC == "Error rate":
        ax.set_yscale("symlog", linthresh=1e-6)
        ax.set_ylabel("Error rate" if dataset == "HCV" else "")
    else:
        ax.set_yscale("linear")
        ax.set_ylabel("Genome fraction (%)" if dataset == "HCV" else "")

    ax.set_title(f"{dataset} {METRIC}")
    ax.set_xlabel("")
    ax.tick_params(axis="x", rotation=30)

plt.tight_layout()
out_path = f"/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/graphs/box_plot_{METRIC.lower().replace(' ', '_')}/box_plot_hcv_hiv.png"
os.makedirs(os.path.dirname(out_path), exist_ok=True)
fig.savefig(out_path, dpi=300, bbox_inches="tight")
print(f"Saved figure to {out_path}")
plt.show()


def wilcoxon_test_mfd_vs_other(df_wide, other_method, metric=METRIC):
    paired = df_wide.dropna(subset=["MFD", other_method])
    stat, pvalue = wilcoxon(
        paired["MFD"],
        paired[other_method],
        alternative='less' if metric == "Error rate" else 'greater'  # For genome fraction, test if MFD > others
    )
    print(f"Wilcoxon signed-rank test MFD vs {other_method} ({metric}): stat={stat:.3f}, p-value={pvalue:.4f}")


print(f"HCV dataset tests for metric: {METRIC}")
for method in ["VG-flow", "Virus-VG", "SAVAGE"]:
    wilcoxon_test_mfd_vs_other(df_hcv_wide, method, metric=METRIC)

print(f"\nHIV dataset tests for metric: {METRIC}")
for method in ["VG-flow", "Virus-VG", "SAVAGE"]:
    wilcoxon_test_mfd_vs_other(df_hiv_wide, method, metric=METRIC)
