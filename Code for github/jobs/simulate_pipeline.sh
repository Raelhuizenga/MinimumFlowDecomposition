#!/bin/sh
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=12GB
#SBATCH --output=output/slurm_files/slurm_%j.out
#SBATCH --error=output/slurm_files/slurm_%j.err

STRAIN=$1
SEED=$2
echo "Running pipeline for strain $STRAIN"

# ========= SETUP COMMON PATHS =========
BASE="/tudelft.net/staff-umbrella/FlowDecomposition"
AFD="$BASE/flow_assembly/AssemblyFlowDecomposition"
DATA="$AFD/Data/simulation_data"
SCRIPT_DIR="$AFD/scripts"

# ========= STEP 1: SIMULATE READS =========
echo "Starting read simulation..."
export APPTAINER_ROOT="$BASE/containers/art_read_simulation"
export APPTAINER_NAME="art-read-simulator.sif"
RAW_READS_DIR="$DATA/simulation_data/art_reads/150bp/HIV_8000/${STRAIN}-strain_seed_${SEED}"
STRAINS_DIR="$DATA/strains/HIV_strains_8000_region"
mkdir -p "$RAW_READS_DIR"

srun --ntasks=1 --cpus-per-task=16 --exclusive --job-name=sim_${STRAIN} \
    apptainer exec \
    -B "$BASE:$BASE" \
    -B "$HOME:$HOME" \
    "$APPTAINER_ROOT/$APPTAINER_NAME" \
    python "$SCRIPT_DIR/simulate_reads.py" \
        -i "$STRAINS_DIR" \
        -o "$RAW_READS_DIR" \
        -n "$STRAIN" \
        -s "$SEED"
        
if [ $? -ne 0 ]; then
    echo "Read simulation failed for strain $STRAIN."
    continue
fi

echo "Starting PEAR merging..."
PEAR_BIN="$BASE/pear"
PEAR_OUTPUT_DIR="$DATA/simulation_data/pear_merged/150bp/HIV_8000/${STRAIN}-strain_seed_${SEED}"
mkdir -p "$PEAR_OUTPUT_DIR"

srun --ntasks=1 --cpus-per-task=16 --exclusive --job-name=pear_${STRAIN} \
    "$PEAR_BIN" \
    -f "${RAW_READS_DIR}/all_reads_1.fq" \
    -r "${RAW_READS_DIR}/all_reads_2.fq" \
    -o "${PEAR_OUTPUT_DIR}/pear_read"

if [ $? -ne 0 ]; then
    echo "PEAR merging failed for strain $STRAIN."
    continue
fi

echo "Starting SAVAGE assembly..."
export APPTAINER_ROOT="$BASE/containers/SAVAGE"
export APPTAINER_NAME="savage-image.sif"

SAVAGE_JOB_DIR="$AFD/output/${SLURM_JOB_ID}/savage_${STRAIN}"
mkdir -p "$SAVAGE_JOB_DIR"
cd "$SAVAGE_JOB_DIR"

srun --ntasks=1 --cpus-per-task=16 --exclusive --job-name=savage_${STRAIN} \
    apptainer exec \
    -B "$BASE:$BASE" \
    -B "$HOME:$HOME" \
    "$APPTAINER_ROOT/$APPTAINER_NAME" \
    savage \
        --split "$STRAIN" \
        -s "$PEAR_OUTPUT_DIR/pear_read.assembled.fastq" \
        -p1 "$PEAR_OUTPUT_DIR/pear_read.unassembled.forward.fastq" \
        -p2 "$PEAR_OUTPUT_DIR/pear_read.unassembled.reverse.fastq" \
        -t 16

if [ $? -eq 0 ]; then
    echo "SAVAGE finished successfully for strain $STRAIN."
    CONTIGS_DIR="$DATA/simulation_data/contigs/HIV_8000/${STRAIN}-strain_seed_${SEED}"
    mkdir -p "$CONTIGS_DIR"
    cp "$SAVAGE_JOB_DIR/contigs_stage_c.fasta" "$CONTIGS_DIR/contigs_stage_c.fasta"
    if [ $? -eq 0 ]; then
    echo "Contigs copied for strain $STRAIN."
    else
    echo "Failed to copy contigs for strain $STRAIN."
    fi
else
    echo "SAVAGE failed for strain $STRAIN."
fi