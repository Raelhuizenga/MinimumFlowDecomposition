#!/bin/sh
#SBATCH --partition=general
#SBATCH --qos=medium
#SBATCH --time=24:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=64GB
#SBATCH --mail-type=END
#SBATCH --output=output/slurm_files/slurm_k_%j.out
#SBATCH --error=output/slurm_files/slurm_k_%j.err
#SBATCH --exclude=3dgi1,3dgi2,gpu01,gpu02,gpu03,gpu04,gpu05,gpu06,gpu07,gpu08,gpu09,gpu10,gpu11,gpu12,gpu14,gpu15,gpu16,gpu17,gpu18,gpu19,gpu20,gpu21,gpu22,gpu23,gpu24,gpu25,gpu26,gpu27,gpu28,gpu29,gpu30,gpu31,gpu32,gpu33,gpu34,gpu35,influ2,influ4,influ5,influ6,cor1

# Check for input arguments
if [ $# -ne 2 ]; then
  echo "Usage: sbatch $0 <strain> <seed>"
  exit 1
fi

STRAIN=$1
SEED=$2

echo "strain: ${STRAIN} seed: ${SEED}"

# Set paths
CONTAINER="/tudelft.net/staff-umbrella/FlowDecomposition/containers/flow_paths/flow-paths.sif"
SCRIPT="/mnt/run_flowpaths/scripts/weights_superset.py"
GRAPH_DIR="/mnt/data/simulation_data/simulation_data/seqwish_graphs/HIV_8000_dagify"
K_DIR="/mnt/run_flowpaths/output/estimate_k/HIV_8000_dagify"
OUT_DIR='/mnt/run_flowpaths/output/k-means_clustering/zero_weight_HIV'

# Run the script with parameters inside the container
srun apptainer exec \
  --bind /tudelft.net/staff-umbrella/FlowDecomposition/flowpaths:/mnt/flowpaths \
  --bind /tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition/Data:/mnt/data \
  --bind /tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths:/mnt/run_flowpaths \
  "$CONTAINER" \
  python "$SCRIPT" --strain "$STRAIN" --seed "$SEED" --graph_dir "$GRAPH_DIR" --k_dir "$K_DIR" --out_dir "$OUT_DIR"

# Exit handling
if [ $? -eq 0 ]; then
  echo "Script finished successfully."
else
  echo "Script failed."
  exit 1
fi