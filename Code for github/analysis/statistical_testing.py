import os
import re
import pandas as pd
from scipy.stats import wilcoxon

# === Configuration ===
BASE_DIR = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/quast_results"

METHOD_DIRS = {
    "Perfect weights": f"{BASE_DIR}/HIV_weights_perfect",
    "k-means weights": f"{BASE_DIR}/HIV_k-means_weights",
    "VG-flow": f"{BASE_DIR}/HIV_vg-flow"
}

included_combinations_HCV = [
    (3,0), (3,1), (3,2), (3,3), (3,4), (3,5), (3,6), (3,7), (3,8), (3,9),
    (4,0), (4,1), (4,2), (4,3), (4,4), (4,5), (4,6), (4,7), (4,8), (4,9),
    (5,0), (5,1), (5,4), (5,5), (5,6), (5,7), (5,8),
    (6,1), (6,4), (6,5), (6,6), (6,9),
    (7,6), (7,8), (7,9),
]

included_combinations_HIV = [
    (2,0), (2,1), (2,2), (2,3), (2,4), (2,5), (2,6), (2,7), (2,8), (2,9),
    (3,12), (3,13), (3,14), (3,18), (3,19),
    (4,21), (4,22), (4,23), (4,24), (4,25), (4,26), (4,27), (4,28), (4,29),
    (5,32), (5,33), (5,36), (5,37),
]

FIELDS_ERROR = {
    "# N's per 100 kbp": float,
    "# mismatches per 100 kbp": float,
    "# indels per 100 kbp": float,
}

FIELDS_GENOME_FRAC = {
    "Genome fraction (%)": float,
}

# === Extraction functions ===
def extract_metric(base_dir, included_combinations, metric):
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


def collect_df(method_dirs, included_combinations, metric):
    rows = []
    for method, base_dir in method_dirs.items():
        data = extract_metric(base_dir, included_combinations, metric)
        for s, seed, val in data:
            rows.append({
                "Method": method,
                metric: val,
                "Strain": s,
                "Seed": seed,
            })
    return pd.DataFrame(rows)


def wilcoxon_test(df_wide, other_method, metric):
    paired = df_wide.dropna(subset=["Perfect weights", other_method])
    if paired.empty:
        print(f"No paired data for {other_method}.")
        return

    stat, pvalue = wilcoxon(
        paired["Perfect weights"],
        paired[other_method],
        alternative='less' if metric == "Error rate" else 'greater'
    )
    print(f"Wilcoxon test Perfect weights vs {other_method} ({metric}): stat={stat:.3f}, p-value={pvalue:.4g}")


if __name__ == "__main__":
    for metric in ["Error rate", "Genome fraction (%)"]:
        print("=" * 60)
        print(f"Testing metric: {metric}")
        df = collect_df(METHOD_DIRS, included_combinations_HIV, metric)

        if df.empty:
            print("No data found.")
            continue

        df_wide = df.pivot_table(
            index=["Strain", "Seed"],
            columns="Method",
            values=metric
        ).reset_index()

        for method in ["k-means weights", "VG-flow"]:
            wilcoxon_test(df_wide, method, metric)
