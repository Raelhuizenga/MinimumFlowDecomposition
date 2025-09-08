import os
import re
import json
import numpy as np
from collections import defaultdict


# All combinations HCV
# included_combinations = [
#     (3,0), (3,1), (3,2), (3,3), (3,4), (3,5), (3,6), (3,7), (3,8), (3,9), 
#     (4,0), (4,1), (4,2), (4,3), (4,4), (4,5), (4,6), (4,7), (4,8), (4,9), 
#     (5,0), (5,1), (5,2), (5,3), (5,4), (5,5), (5,6), (5,7), (5,8), (5,9), 
#     (6,0), (6,1), (6,2), (6,3), (6,4), (6,5), (6,6), (6,7), (6,8), (6,9), 
#     (7,0), (7,1), (7,2), (7,3), (7,4), (7,5), (7,6), (7,7), (7,8), (7,9)
#                          ]

# # combinations 
# included_combinations = [
#     (3,0), (3,1), (3,2), (3,3), (3,4), (3,5), (3,6), (3,7), (3,8), (3,9), 
#     (4,0), (4,1), (4,2), (4,3), (4,4), (4,5), (4,6), (4,7), (4,8), (4,9), 
#     (5,0), (5,1), (5,2), (5,3), (5,4), (5,5), (5,6), (5,7), (5,8), (5,9), 
#     (6,0), (6,1), (6,2), (6,3), (6,4), (6,5), (6,6), (6,7), (6,8), (6,9), 
#     (7,0), (7,1), (7,2), (7,3), (7,4), (7,5), (7,6), (7,7), (7,8), (7,9)
#                          ]

# # Slow combinations HCV
included_combinations = [
   (5,2), (5,7), (5,8),
   (6,0), (6,2), (6,6)
]

# SLow combinations HIV 
# (4,20) and (4,22) are not in here because Virus-VG did not finish
# included_combinations = [
#    (2,0), (2,1),
#    (3,14), 
#    (5,32), (5,33)
# ]

# k-means clustering combinations HCV
# included_combinations = [
#     (3,0), (3,1), (3,2), (3,3), (3,4), (3,5), (3,6), (3,7), (3,8), (3,9),
#     (4,0), (4,1), (4,2), (4,3), (4,4), (4,5), (4,6), (4,7), (4,8), (4,9), 
#     (5,0), (5,1), (5,4), (5,5), (5,6), (5,7), (5,8),
#     (6,1), (6,4), (6,5), (6,6), (6,9),
#     (7,6), (7,8), (7,9)
#                          ]

# Combinations HIV that finished fast
# included_combinations = [
#     (2,2), (2,3), (2,4), (2,5), (2,6), (2,7), (2,8), (2,9),
#     (3, 10), (3, 13), (3, 15), (3, 16), (3, 18), (3, 19),
#     (4,20), (4, 29),
#                          ]

# k-means clustering combinations HIV
# included_combinations = [
#     (2,0), (2,1), (2,2), (2,3), (2,5), (2,6), (2,7), (2,8), (2,9),
#     (3, 12), (3, 13), (3, 14), (3, 16), (3, 18), (3, 19),
#     (4, 21), (4, 22), (4, 23), (4, 24), (4, 25), (4, 26), (4, 27), (4, 28), (4,29),
#     (5, 32), (5, 33), (5, 36), (5, 37)
#                          ]

# Directory containing strain_seed subdirectories
BASE_DIR = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/quast_results/HCV_vg-flow_constraints"

FIELDS = {
    "Genome fraction (%)": float,
    "Duplication ratio": float,
    "# N's per 100 kbp": float,
    "# mismatches per 100 kbp": float,
    "# indels per 100 kbp": float,
    "NG50": int,
    "# contigs": int,
}

strain_stats = defaultdict(list)

for subdir in os.listdir(BASE_DIR):
    if not subdir.startswith("strain_") or "_seed_" not in subdir:
        continue

    match = re.match(r"strain_(\d+)_seed_(\d+)", subdir)
    if not match:
        continue

    strain_id = int(match.group(1))
    seed_id = int(match.group(2))

    if (strain_id, seed_id) not in included_combinations:
        continue  # skip excluded combinations

    report_path = os.path.join(BASE_DIR, subdir, "report.txt")
    if not os.path.isfile(report_path):
        continue

    stats = {}
    with open(report_path, "r") as f:
        for line in f:
            line = line.strip()
            for field, dtype in FIELDS.items():
                if line.startswith(field):
                    try:
                        value = dtype(line.split()[-1])
                        stats[field] = value
                    except:
                        pass

    if all(k in stats for k in ["# N's per 100 kbp", "# mismatches per 100 kbp", "# indels per 100 kbp"]):
        error_rate = (
            stats["# N's per 100 kbp"]
            + stats["# mismatches per 100 kbp"]
            + stats["# indels per 100 kbp"]
        ) / 1000
    else:
        error_rate = None

    strain_stats[strain_id].append({
        "Genome fraction (%)": stats.get("Genome fraction (%)"),
        "Duplication ratio": stats.get("Duplication ratio"),
        "Error rate": error_rate,
        "NG50": stats.get("NG50"),
        "# contigs": stats.get("# contigs")
    })

result = {}

for strain_id, records in sorted(strain_stats.items()):
    stat_values = defaultdict(list)
    for record in records:
        for k, v in record.items():
            if v is not None:
                stat_values[k].append(v)

    strain_result = {}
    for k in ["Genome fraction (%)", "Duplication ratio", "Error rate", "NG50", "# contigs"]:
        values = stat_values.get(k, [])
        if values:
            strain_result[k] = round(np.mean(values), 3)
            strain_result[f"{k}_std"] = round(np.std(values), 3)
        else:
            strain_result[k] = None
            strain_result[f"{k}_std"] = None

    # Only count seeds for that strain
    n_total = sum(1 for s, _ in included_combinations if s == strain_id)
    strain_result["Seeds found"] = f"{len(records)}/{n_total}"

    result[f"strain_{strain_id}"] = strain_result

result_file = "strain_quast_stats_subset_slow_instances.json"

with open(f"{BASE_DIR}/{result_file}", "w") as f:
    json.dump(result, f, indent=4)

print(f"Summary written to {result_file}")
