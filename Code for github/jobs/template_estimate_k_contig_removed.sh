#!/bin/sh
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --time=0:10:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=8GB
#SBATCH --output=output/slurm_files/slurm_strain__STRAIN___seed__SEED___%j.out
#SBATCH --error=output/slurm_files/slurm_strain__STRAIN___seed__SEED___%j.err
#SBATCH --job-name=strain__STRAIN___seed__SEED__

echo "Running strain=__STRAIN__ seed=__SEED__"


CONTAINER="/tudelft.net/staff-umbrella/FlowDecomposition/containers/flow_paths/flow-paths.sif"
SCRIPT="/mnt/run_flowpaths/scripts/path_cover.py"
OUTPUT_DIR="/mnt/run_flowpaths/output/estimate_k/HCV_contigs_removed"
INPUT_DIR="/mnt/data/simulation_data/simulation_data/seqwish_graphs/random_abundances/wrong_contigs_removed"

srun apptainer exec \
  --bind /tudelft.net/staff-umbrella/FlowDecomposition/flowpaths:/mnt/flowpaths \
  --bind /tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/AssemblyFlowDecomposition/Data:/mnt/data \
  --bind /tudelft.net/staff-umbrella/FlowDecomposition/flow_assembly/run_flowpaths:/mnt/run_flowpaths \
  "$CONTAINER" \
    /usr/bin/time -v -o "$OUTPUT_DIR/strain__STRAIN___seed__SEED___.time" -a \
  python "$SCRIPT" \
    --strain __STRAIN__ \
    --seed __SEED__ \
    --subpaths True \
    --vgflow_constraints False \
    --output_dir "$OUTPUT_DIR" \
    --input_dir "$INPUT_DIR"
