# Import tools
import numpy as np # type: ignore
from fastlmm.association import single_snp # type: ignore
from fastlmm.association import single_snp_linreg # type: ignore
from fastlmm.util import example_file # type: ignore
import pandas as pd # type: ignore
import os

# Define data paths for fast LMM function arguments
bed_path = 'fastLMM.bed'
pheno_path = 'fastLMM_clean.fam'
ld_bed = 'fastLMM_ld.bed'

bed_path_covar = 'fastLMM_covar.bed'
pheno_path_covar = 'fastLMM_covar_clean.fam'
ld_bed_covar = 'fastLMM_covar_ld.bed'
covar_path_covar = 'multi_column_covar.fam'

################################################################################################################

# No covariates
print("# Running GWAS with No Covariates")
## Lin-reg
print("# Model: Linear Regression without Covariates")
gwas_results = single_snp_linreg(bed_path, pheno_path, count_A1=True)
gwas_results.to_csv('testing.tsv', sep = '\t')
## Distance matrix
print("# Model: Mixed Model (Distance Matrix) without Covariates")
gwas_results = single_snp(bed_path, pheno_path, K0=ld_bed, count_A1=True).sort_values('PValue')
gwas_results.to_csv('testing_K0_ld.tsv', sep = '\t')

# Covariates
print("# Running GWAS with Covariates")
## Lin-reg
print("# Model: Linear Regression with Covariates")
gwas_results = single_snp_linreg(bed_path_covar, pheno_path_covar, covar=covar_path_covar, count_A1=True)
gwas_results.to_csv('testing_covar.tsv', sep = '\t')
## Distance matrix
print("# Model: Mixed Model (Distance Matrix) with Covariates")
gwas_results =  single_snp(bed_path_covar, pheno_path_covar, K0=ld_bed_covar, covar=covar_path_covar, count_A1=True).sort_values('PValue')
gwas_results.to_csv('testing_K0_ld_covar.tsv', sep = '\t')
