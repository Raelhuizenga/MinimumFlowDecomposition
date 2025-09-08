#!/bin/sh
#SBATCH --partition=general
#SBATCH --qos=short
#SBATCH --time=4:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16GB
#SBATCH --output=output/slurm_files/slurm_%j.out
#SBATCH --error=output/slurm_files/slurm_%j.err

# Define container folder and name
export APPTAINER_ROOT="/tudelft.net/staff-umbrella/FlowDecomposition/containers/seqwish_graph"
export APPTAINER_NAME="seqwish.sif"

# Define working directories
WORKDIR="/tudelft.net/staff-umbrella/FlowDecomposition"
SCRIPT_DIR="$WORKDIR/flow_assembly/AssemblyFlowDecomposition/scripts"
DATA_DIR="$WORKDIR/flow_assembly/AssemblyFlowDecomposition/Data/simulation_data"

for STRAIN in 2 3 4 5; do
  if [ "$STRAIN" -eq 2 ]; then
    SEED_START=0
    SEED_END=9
  elif [ "$STRAIN" -eq 3 ]; then
    SEED_START=10
    SEED_END=19
  elif [ "$STRAIN" -eq 4 ]; then
    SEED_START=20
    SEED_END=29
  elif [ "$STRAIN" -eq 5 ]; then
    SEED_START=30
    SEED_END=39
  fi

  for SEED in $(seq $SEED_START $SEED_END); do
    OUTPUT_DIR="$DATA_DIR/simulation_data/seqwish_graphs/HIV_8000_dagify/${STRAIN}-strain_seed_${SEED}"

    if [ -d "$OUTPUT_DIR" ]; then
      echo "Skipping NUM_STRAINS=$STRAIN SEED=$SEED: Output already exists at $OUTPUT_DIR"
      continue
    fi

    echo "Running pipeline for NUM_STRAINS=$STRAIN SEED=$SEED"

    READ1="$DATA_DIR/simulation_data/art_reads/150bp/HIV_8000/${STRAIN}-strain_seed_${SEED}/all_reads_1.fq"
    READ2="$DATA_DIR/simulation_data/art_reads/150bp/HIV_8000/${STRAIN}-strain_seed_${SEED}/all_reads_2.fq"
    CONTIGS="$DATA_DIR/simulation_data/contigs/HIV_8000/${STRAIN}-strain_seed_${SEED}"

    # export TMPDIR="$CONTIGS/tmp"
    # mkdir -p "$TMPDIR"
    # cd "$CONTIGS"

    # # Index contigs
    # apptainer exec -B "$WORKDIR:$WORKDIR" -B "$TMPDIR:$TMPDIR" "$APPTAINER_ROOT/$APPTAINER_NAME" \
    #   kallisto index -i "$CONTIGS/contigs.idx" "$CONTIGS/contigs_stage_c.fasta"

    # for i in {1..5}; do
    #   if [ -f "$CONTIGS/contigs.idx" ]; then
    #     break
    #   else
    #     echo "Waiting for kallisto index to be written..."
    #     sleep 2
    #   fi
    # done

    # # Quantify abundance
    # apptainer exec -B "$WORKDIR:$WORKDIR" -B "$TMPDIR:$TMPDIR" "$APPTAINER_ROOT/$APPTAINER_NAME" \
    #   kallisto quant -i "$CONTIGS/contigs.idx" -o "$CONTIGS/kallisto_out" -b 100 "$READ1" "$READ2"

    # # Filter contigs
    # apptainer exec -B "$WORKDIR:$WORKDIR" "$APPTAINER_ROOT/$APPTAINER_NAME" \
    #   python "$SCRIPT_DIR/filter_contigs_by_abundance.py" \
    #   --contigs "$CONTIGS/contigs_stage_c.fasta" \
    #   --abundance "$CONTIGS/kallisto_out/abundance.tsv" \
    #   --output "$CONTIGS/filtered_contigs.fasta" \
    #   --threshold "2.0"

    # if [ $? -eq 0 ]; then
    #   echo "Kallisto filtering finished successfully."
    # fi

    mkdir -p "$OUTPUT_DIR"

    # Build graph
    srun apptainer exec \
      -B "$WORKDIR:$WORKDIR" \
      "$APPTAINER_ROOT/$APPTAINER_NAME" \
      python "$SCRIPT_DIR/build_variation_graph.py" \
      "$CONTIGS/filtered_contigs.fasta" \
      "$READ1" \
      "$READ2" \
      "$OUTPUT_DIR" \
      -vg "$WORKDIR/new_vg"

    # rm -rf "$TMPDIR"

    if [ $? -eq 0 ]; then
      echo "Graph building finished successfully for NUM_STRAINS=$STRAIN SEED=$SEED"
    fi

  done
done
