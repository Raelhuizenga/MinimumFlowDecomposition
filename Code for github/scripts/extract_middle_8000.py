import os
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from pathlib import Path

def extract_middle_8000_from_folder(input_folder, output_folder):
    input_folder = Path(input_folder)
    output_folder = Path(output_folder)
    output_folder.mkdir(parents=True, exist_ok=True)

    num_written = 0

    for fasta_file in input_folder.glob("*"):
        if not fasta_file.is_file():
            continue

        try:
            record = next(SeqIO.parse(fasta_file, "fasta"))
        except Exception as e:
            print(f"Skipping {fasta_file.name}: error parsing FASTA ({e})")
            continue

        seq_len = len(record.seq)
        if seq_len < 8000:
            print(f"Skipping {fasta_file.name}: sequence too short ({seq_len} bp)")
            continue

        start = (seq_len - 8000) // 2
        end = start + 8000
        middle_seq = record.seq[start:end]

        new_record = SeqRecord(
            middle_seq,
            id=record.id,
            description=f"{record.description} | middle 8000bp ({start}-{end})"
        )

        # Construct output filename with _8000.fasta suffix
        output_file = output_folder / (fasta_file.stem + "_8000.fasta")
        SeqIO.write(new_record, output_file, "fasta")
        num_written += 1

    print(f"Extracted and saved middle 8000 bp for {num_written} sequences in '{output_folder}'.")

if __name__ == "__main__":
    input_folder = "flow_assembly/AssemblyFlowDecomposition/Data/simulation_data/strains/HIV_strains"
    output_folder = "flow_assembly/AssemblyFlowDecomposition/Data/simulation_data/strains/HIV_strains_8000_region"
    extract_middle_8000_from_folder(input_folder, output_folder)
