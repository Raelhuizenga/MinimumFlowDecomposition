import argparse
import subprocess
from pathlib import Path
import random
import shutil
import numpy as np

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input_directory', type=str, required=True, help='The directory where the fasta files of the strains are located')
    parser.add_argument('-o', '--output_directory', type=str, required=True, help='The directory where the reads should be written to')
    parser.add_argument('-n', '--num_strains', type=int, default=10, help='The number of strains to simulate')
    parser.add_argument('-s', '--seed', type=int, help='Random seed for reproducibility')
    parser.add_argument('-strains', '--strains', type=int, nargs='+', help='List of indices of strains to use')
    parser.add_argument('-f', '--frequency', type=float, default=0, help='Frequency minor strain in mixture')
    args = parser.parse_args()

    generate_reads(args.input_directory, args.output_directory, args.num_strains, args.seed, args.strains, args.frequency)

def generate_reads(input_directory, output_directory, num_strains, seed, strains, frequency):
    """
    Generate reads from the reference genomes and save them to the output directory.
    Also stores the ground truth abundances in a text file.
    """
    if seed == None:
        seed = random.randint(0, 2**32 - 1)
    else:
        random.seed(seed)

    input_path = Path(input_directory)
    output_path = Path(output_directory)
    output_path.mkdir(parents=True, exist_ok=True)

    fasta_files = sorted(input_path.glob("*.fasta"))
    if len(fasta_files) < num_strains:
        raise ValueError(f"Requested {num_strains} strains, but only found {len(fasta_files)} fasta files.")

    if not strains:
        selected_strains = random.sample(fasta_files, num_strains)
    else:
        selected_strains = [fasta_files[i] for i in strains]
    abundances = generate_constrained_exponential(num_strains, seed)

    if frequency != 0:
        abundances = np.array([frequency*1500, (1-frequency)*1500])

    read_files_1 = []
    read_files_2 = []

    for i, (strain_file, abundance) in enumerate(zip(selected_strains, abundances)):
        prefix = output_path / f"reads_{i}"
        simulate_reads(strain_file, abundance, prefix, seed + i)
        read_files_1.append(f"{prefix}1.fq")
        read_files_2.append(f"{prefix}2.fq")

    # Store ground truth
    ground_truth_path = output_path / "ground_truth_abundances.txt"
    write_ground_truth_fasta(selected_strains, abundances, ground_truth_path)

    all_reads_1 = output_path / "all_reads_1.fq"
    all_reads_2 = output_path / "all_reads_2.fq"
    # Merge reads into single files
    merge_reads(read_files_1, all_reads_1)
    merge_reads(read_files_2, all_reads_2)

    # Trim reads using cutadapt
    # trimmed_1 = output_path / "all_reads_trimmed_1.fq"
    # trimmed_2 = output_path / "all_reads_trimmed_2.fq"
    # trim_reads_with_cutadapt(all_reads_1, all_reads_2, trimmed_1, trimmed_2)


def write_ground_truth_fasta(strain_files, abundances, output_path):
    """
    Write a FASTA file where each sequence is labeled with its abundance.

    Parameters:
    - strain_files: list of Path objects pointing to the selected FASTA files
    - abundances: list or array of abundance values corresponding to each strain
    - output_path: Path object for the output FASTA file
    """
    output_path = Path(output_path)
    with open(output_path, "w") as out_fasta:
        for i, (strain_file, abundance) in enumerate(zip(strain_files, abundances)):
            with open(strain_file, "r") as f:
                lines = f.readlines()
                sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
            out_fasta.write(f">{strain_file.stem}, weight: {abundance:.1f}x\n{sequence}\n")


# def simulate_reads_old(strain_file, abundance, output_prefix, seed):
#     """
#     Simulate reads from a strain and save them to paired FASTQ files using ART.
#     """
#     cmd = (
#         f'art_illumina '
#         f'-ss MSv3 -sam '
#         f'-i "{strain_file}" '
#         f'-p -l 250 -f {abundance} '
#         f'-m 600 -s 10 '
#         f'-rs {seed} '
#         f'-o "{output_prefix}"'
#     )
#     subprocess.check_call(cmd, shell=True)

def simulate_reads(strain_file, abundance, output_prefix, seed):
    """
    Simulate reads from a strain and save them to paired FASTQ files using ART.
    """
    cmd = (
        f'art_illumina '
        f'-ss HS25 '
        f'-i "{strain_file}" '
        f'-p -l 150 -f {abundance} '
        f'-m 300 -s 10 '
        f'-rs {seed} '
        f'-o "{output_prefix}"'
    )
    subprocess.check_call(cmd, shell=True)

def merge_reads(input_files, output_file):
    """
    Concatenate multiple FASTQ files into a single output file.
    """
    with open(output_file, 'wb') as wfd:
        for f in input_files:
            with open(f, 'rb') as fd:
                shutil.copyfileobj(fd, wfd)


# def trim_reads_with_cutadapt(read1_path, read2_path, output1_path, output2_path):
#     """
#     Use cutadapt to trim reads.
#     """
#     cmd = (
#         f'cutadapt -u 50 -U 50 -q 20,20 '
#         f'-o "{output1_path}" -p "{output2_path}" '
#         f'"{read1_path}" "{read2_path}"'
#     )
#     subprocess.check_call(cmd, shell=True)


def generate_constrained_exponential(n, seed, min_val=0.05, average_coverage=750):
    np.random.seed(seed)  # Set the random seed for reproducibility
    if min_val * n > 1:
        raise ValueError("Minimum value too large to satisfy sum constraint.")
    
    remaining = 1.0 - min_val * n
    raw = np.random.exponential(1, n)
    raw_normalized = raw / raw.sum()
    scaled = raw_normalized * remaining
    final = scaled + min_val
    return np.round(final * average_coverage * n, 0)

if __name__ == "__main__":
    main()
