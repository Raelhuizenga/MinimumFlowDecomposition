import os
import re
from pathlib import Path

def filter_fasta_by_weight(input_dir, output_dir, min_weight=1.0):
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    fasta_files = list(input_dir.glob("solution*")) 

    weight_pattern = re.compile(r"weight:\s*([\d.]+)x")

    for fasta_file in fasta_files:
        output_lines = []
        keep = False
        with open(fasta_file, "r") as f:
            for line in f:
                if line.startswith(">"):
                    match = weight_pattern.search(line)
                    if match and float(match.group(1)) >= min_weight:
                        keep = True
                        output_lines.append(line)
                    else:
                        keep = False
                else:
                    if keep:
                        output_lines.append(line)

        # Write the filtered data only if it's non-empty
        if output_lines:
            out_path = output_dir / fasta_file.name
            with open(out_path, "w") as f_out:
                f_out.writelines(output_lines)
            print(f"Filtered: {fasta_file.name} : {out_path.name}")

if __name__ == "__main__":
    input_directory = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/k-means_clustering/small_weight"     # CHANGE THIS
    output_directory = "/tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths/output/k-means_clustering/small_weight/filter_out_small_weight"      # CHANGE THIS
    filter_fasta_by_weight(input_directory, output_directory)