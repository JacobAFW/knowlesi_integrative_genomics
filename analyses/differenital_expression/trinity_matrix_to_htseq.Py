import pandas as pd
import os
import argparse

# Set up command-line argument parser
parser = argparse.ArgumentParser(description="Convert Trinity matrix to per-sample HTSeq-format files.")
parser.add_argument("input_file", help="Path to the Trinity matrix file (e.g., trinity_matrix.txt)")
parser.add_argument("--output_dir", default="ht_seq", help="Directory to save HTSeq files (default: ht_seq)")
args = parser.parse_args()

# Load Trinity matrix
df = pd.read_csv(args.input_file, sep="\t", index_col=0)

# Create output directory if it doesn't exist
os.makedirs(args.output_dir, exist_ok=True)

# Loop through each sample (column)
for sample in df.columns:
    sample_counts = df[[sample]].copy()
    sample_counts.reset_index(inplace=True)
    sample_counts.columns = ["gene_id", "count"]
    
    # Round to nearest integer for HTSeq/DE compatibility
    sample_counts["count"] = sample_counts["count"].round().astype(int)

    # Construct output file path
    output_path = os.path.join(args.output_dir, f"{sample}.txt")

    # Save in HTSeq format
    sample_counts.to_csv(output_path, sep="\t", index=False, header=False)
    print(f"Generated: {output_path}")