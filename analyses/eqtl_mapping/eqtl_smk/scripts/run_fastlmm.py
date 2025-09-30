#!/usr/bin/env python

import argparse
from fastlmm.association import single_snp

parser = argparse.ArgumentParser()
parser.add_argument("--pheno-gene", required=True)
parser.add_argument("--expr-matrix", required=True)
parser.add_argument("--covar", required=True)
parser.add_argument("--plink-prefix", required=True)
parser.add_argument("--out", required=True)
parser.add_argument("--ld-bed-covar", required=True)
args = parser.parse_args()

# Just extract the correct phenotype to a new file
import pandas as pd

# Read and filter pheno
expr_df = pd.read_csv(args.expr_matrix, sep="\t")
pheno_df = expr_df[["FID", "IID", args.pheno_gene]].dropna()
pheno_df.columns = ["FID", "IID", "pheno"]
pheno_path = args.out + ".pheno.tmp"
pheno_df.to_csv(pheno_path, sep="\t", index=False, header=False)
ld_bed_covar= args.ld_bed_covar

# Run FastLMM
results = single_snp(
    test_snps=args.plink_prefix,
    pheno=pheno_path,
    covar=args.covar,
    K0=ld_bed_covar,
    count_A1=False  # suppress warning
)

results.to_csv(args.out, sep="\t", index=False)