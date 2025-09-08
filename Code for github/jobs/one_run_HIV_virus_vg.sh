#!/bin/sh
#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --mail-type=END
#SBATCH --output=output/slurm_files/slurm_%j.out
#SBATCH --error=output/slurm_files/slurm_%j.err
#SBATCH --exclude=3dgi1,3dgi2,gpu01,gpu02,gpu03,gpu04,gpu05,gpu06,gpu07,gpu08,gpu09,gpu10,gpu11,gpu12,gpu14,gpu15,gpu16,gpu17,gpu18,gpu19,gpu20,gpu21,gpu22,gpu23,gpu24,gpu25,gpu26,gpu27,gpu28,gpu29,gpu30,gpu31,gpu32,gpu33,gpu34,gpu35,influ2,influ4,influ5,influ6,cor1


echo "Starting job"

# Constants
CONTAINER="/tudelft.net/staff-umbrella/FlowDecomposition/containers/virus-vg/virus-vg.sif"
SCRIPT="/mnt/virus-vg/jbaaijens-virus-vg-69a05f3e74f2/scripts/optimize_strains.py"
REMAPPER="/mnt/virus-vg/jbaaijens-virus-vg-69a05f3e74f2/scripts/remap_gfa_ids.py"

# Set strain and seed values (single job version)
strain=3
seed=16

# Set paths for HIV instance
DATA_DIR="/mnt/flow_assembly/AssemblyFlowDecomposition/Data/simulation_data/simulation_data/seqwish_graphs/HIV_8000_dagify/${strain}-strain_seed_${seed}"
OUTPUT_DIR="/tudelft.net/staff-umbrella/FlowDecomposition/virus-vg/output/HIV/${strain}-strain_seed_${seed}"

# Compute m and c values
m=$(echo "scale=0; 750 * 0.005 * $strain / 1" | bc)
c=$(echo "750 * 0.01 * $strain" | bc)

# Create output dir
mkdir -p "$OUTPUT_DIR"

# === STEP 1: Run ID remapping ===
echo "Remapping IDs..."
srun apptainer exec \
  --bind /tudelft.net/staff-umbrella/FlowDecomposition:/mnt \
  "$CONTAINER" \
  python "$REMAPPER" \
    --abundance-in "${DATA_DIR}/final_node_abundances.txt" \
    --abundance-out "${DATA_DIR}/node_abundance_remapped.txt" \
    --gfa-in "${DATA_DIR}/final_graph.gfa" \
    --gfa-out "${DATA_DIR}/mod_graph_remapped.gfa"

# === STEP 2: Run optimization ===
echo "Running optimization..."
cd "$OUTPUT_DIR"

srun apptainer exec \
  --bind /tudelft.net/staff-umbrella/FlowDecomposition:/mnt \
  "$CONTAINER" \
  python "$SCRIPT" \
    -m "$m" -c "$c" -t 8 -o "haps_${strain}_${seed}.fasta" \
    "${DATA_DIR}/node_abundance_remapped.txt" \
    "${DATA_DIR}/mod_graph_remapped.gfa"

# Check job success
if [ $? -eq 0 ]; then
  echo "Script finished successfully."
else
  echo "Script failed."
  exit 1
fi
