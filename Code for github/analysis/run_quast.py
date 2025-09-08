import os
import re
import subprocess
from pathlib import Path
import shutil
import time

# In flow_paths container:
# apptainer shell \
#   --bind /tudelft.net/staff-umbrella/FlowDecomposition/flowpaths:/mnt/flowpaths \
#   --bind /tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition/Data:/mnt/data \
#   --bind /tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths:/mnt/run_flowpaths \
#   containers/flow_paths/flow-paths.sif


# Paths inside the container
# output_dir = Path("flow_assembly/run_flowpaths/output/min_cover_strict_constraints/nodes_removed")  # Your solution files
output_dir = Path("/mnt/run_flowpaths/output/min_cover_strict_constraints/nodes_removed")  # Your solution files
# output_dir = Path("/mnt/output_HCV")
# output_dir = Path("/mnt/data/simulation_data/simulation_data/contigs/HIV_8000")  # Your solution files
quast_output_base = Path("/mnt/run_flowpaths/quast_results/HCV_vg-flow_constraints")  # Final results location
# quast_output_base = Path("flow_assembly/run_flowpaths/quast_results/HCV_filtered_strict_constraints")  # Final results location
# quast_output_base = Path("flow_assembly/run_flowpaths/quast_results/HIV_weights_perfect")  # Final results location
reference_base = Path("flow_assembly/AssemblyFlowDecomposition/Data/simulation_data/simulation_data/art_reads/150bp/random_abundances")  # Ground truth location

# Ensure result base exists
quast_output_base.mkdir(parents=True, exist_ok=True)

# Pattern to match output files like solution_6_1
# pattern = re.compile(r"haps_(\d+)_(\d+)$")
# pattern = re.compile(r"(\d+)-strain_seed_(\d+)$")
pattern = re.compile(r"solution_(\d+)_(\d+)$")

for file in output_dir.iterdir():
    if file.is_file() or file.is_dir():
        match = pattern.match(file.stem)
        if not match:
            print('no match')
            continue
        strain, seed = match.groups()
        if file.is_dir():
            # strain_file = f'{file}/filtered_contigs.fasta'
            strain_file = f'{file}/haps_{strain}_{seed}.fasta'
        else:
            strain_file = file

        
        # if strain != '5' and seed != '2' and seed != '5' and seed != '7' and seed != '8':
        #     continue
        reference_file = reference_base / f"{strain}-strain_seed_{seed}/ground_truth_abundances.txt"

        if not reference_file.exists():
            print(f"Reference file not found for strain {strain} seed {seed}, skipping.")
            continue

        print(f"Running QUAST for strain {strain}, seed {seed}...")

        quast_output = Path(f"{quast_output_base}/strain_{strain}_seed_{seed}")

        subprocess.run([
            "quast.py", "-m", "500",
            "-R", str(reference_file),
            "-o", str(quast_output),
            str(strain_file)
        ], check=True)

