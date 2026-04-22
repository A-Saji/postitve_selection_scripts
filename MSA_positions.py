#!/usr/bin/env python3
import csv

def read_fasta(path):
    seqs = {}
    current = None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                seqs[current] = []
            else:
                seqs[current].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def main():
    msa_file = "transAln_cds_300_GRMZM2G005155.aa_cleanali.fasta"   # ← change if needed
    out_csv = "MSA_positions_table.csv"

    print(f"Reading MSA: {msa_file}")

    seqs = read_fasta(msa_file)
    species = list(seqs.keys())

    # All sequences must be aligned
    aln_len = len(next(iter(seqs.values())))
    for s in species:
        if len(seqs[s]) != aln_len:
            raise ValueError(f"Sequence {s} has different alignment length.")

    print(f"Alignment length = {aln_len}")
    print(f"Number of sequences = {len(species)}")

    # Prepare CSV output
    with open(out_csv, "w", newline="") as f:
        writer = csv.writer(f)
        header = ["Position"] + species
        writer.writerow(header)

        for pos in range(aln_len):
            row = [pos + 1]  # 1-based
            for sp in species:
                row.append(seqs[sp][pos])
            writer.writerow(row)

    print(f"\nWrote output: {out_csv}")


if __name__ == "__main__":
    main()
