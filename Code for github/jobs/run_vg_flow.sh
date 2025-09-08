#!/bin/sh
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=8GB
#SBATCH --mail-type=END
#SBATCH --output=output/slurm_files/slurm_%j.out
#SBATCH --error=output/slurm_files/slurm_%j.err

echo "Starting batch job"

# Set container and base paths
CONTAINER="/tudelft.net/staff-umbrella/FlowDecomposition/containers/vg-flow/vg-flow-image.sif"
BASE_DIR="/tudelft.net/staff-umbrella/FlowDecomposition"
SCRIPT="/mnt/flow_assembly/AssemblyFlowDecomposition/scripts/vg-flow.py"

# Bind the base directory into the container as /mnt
for strain in $(seq 3 7); do
  for seed in $(seq 0 9); do
    echo "Processing strain $strain, seed $seed"

    DATA_DIR="/mnt/flow_assembly/AssemblyFlowDecomposition/Data/simulation_data/simulation_data/seqwish_graphs/random_abundances/wrong_nodes_removed/${strain}-strain_seed_${seed}"
    OUTPUT_DIR="/mnt/flow_assembly/AssemblyFlowDecomposition/output_VG_flow/HCV_nodes_removed"

    m=$(echo "scale=0; 750 * 0.005 * $strain / 1" | bc)
    c=$(echo "750 * 0.01 * $strain" | bc)

    cd "$OUTPUT_DIR"

    srun apptainer exec \
      --bind "$BASE_DIR:/mnt" \
      "$CONTAINER" \
      python "$SCRIPT" -m "$m" -c "$c" -t 8 \
        -o "$OUTPUT_DIR/haps_${strain}_${seed}.fasta" \
        "$DATA_DIR/node_abundance_remapped.txt" \
        "$DATA_DIR/mod_graph_remapped.gfa"

    # Exit check for this combination
    if [ $? -eq 0 ]; then
      echo "Strain $strain seed $seed finished successfully."
    else
      echo "Strain $strain seed $seed failed."
    fi

  done
done

echo "Batch job completed."
