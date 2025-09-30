#!/usr/bin/env python3

import pandas as pd
import argparse

# Argument parser setup
parser = argparse.ArgumentParser(description="Generate a dummy GTF file from a Trinity count matrix.")
parser.add_argument("input_file", help="Path to the Trinity matrix file (e.g., trinity_matrix.txt)")
parser.add_argument("--output_file", default="dummy.gtf", help="Output GTF filename (default: dummy.gtf)")
args = parser.parse_args()

# Load Trinity matrix
df = pd.read_csv(args.input_file, sep="\t", index_col=0)
transcript_ids = df.index.to_list()

# Write fixed GTF
with open(args.output_file, "w") as gtf:
    for tx in transcript_ids:
        # Set both gene_id and transcript_id to the full transcript ID
        line = (
            f"{tx}\tTrinity\texon\t1\t1000\t.\t+\t.\t"
            f'gene_id "{tx}"; transcript_id "{tx}";\n'
        )
        gtf.write(line)

print(f"✅ Fixed dummy GTF saved as: {args.output_file}")