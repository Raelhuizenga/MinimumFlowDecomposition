#!/bin/sh
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16GB
#SBATCH --mail-type=END
#SBATCH --output=output/slurm_files/slurm_%j.out
#SBATCH --error=output/slurm_files/slurm_%j.err
#SBATCH --exclude=3dgi1,3dgi2,gpu01,gpu02,gpu03,gpu04,gpu05,gpu06,gpu07,gpu08,gpu09,gpu10,gpu11,gpu12,gpu14,gpu15,gpu16,gpu17,gpu18,gpu19,gpu20,gpu21,gpu22,gpu23,gpu24,gpu25,gpu26,gpu27,gpu28,gpu29,gpu30,gpu31,gpu32,gpu33,gpu34,influ2,influ4,influ5,influ6,cor1

echo "Starting batch job"

# Set paths
CONTAINER="/tudelft.net/staff-umbrella/FlowDecomposition/containers/flow_paths/flow-paths.sif"
SCRIPT="/mnt/run_flowpaths/scripts/exact_flow_decomposition.py"
DATA_DIR="/mnt/data/simulation_data/gridsearch_length"
OUTPUT_DIR="/mnt/run_flowpaths/output/haplotype_gridsearch/length"

for length in 7500 8000 8500 9000 9500 10000 10500 11000; do
  for repeat in 0 1 2 3 4 5 6 7 8 9; do
    echo "Running simulation for $length length and instance $repeat"

    apptainer exec \
      --bind /tudelft.net/staff-umbrella/FlowDecomposition/flowpaths:/mnt/flowpaths \
      --bind /tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition/Data:/mnt/data \
      --bind /tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths:/mnt/run_flowpaths \
      "$CONTAINER" \
      python "$SCRIPT" --strain 3 \
                       --seed "$repeat" \
                       --length "$length" \
                       --graph_directory "$DATA_DIR" \
                       --output_directory "$OUTPUT_DIR"

    if [ $? -eq 0 ]; then
      echo "Run finished successfully."
    else
      echo "Run failed. Exiting."
    fi
  done
done

echo "All simulations completed."