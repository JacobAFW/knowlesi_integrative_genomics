# Integrative Genomics of Plasmodium knowlesi

This repository contains all analysis code, pipelines, and supporting materials for our study:

**Integrative genomics of _Plasmodium knowlesi_ reveals parasite-intrinsic regulators of severe human malaria** (Westaway et al., 2026)

We present the first integrated genomic and transcriptomic analysis of clinical _P. knowlesi_ infections, combining GWAS, bulk RNA-seq, de novo transcriptome assembly, and eQTL mapping across over 500 human infections from Malaysia. Our multi-omic approach reveals distinct parasite programs associated with parasite burden, disease severity, and host-intrinsic factors — including a previously unrecognised role of host sex in modulating parasite gene expression. 

A key finding is the transcriptional divergence between stress-response programs (linked to severe disease) and invasion/evasion pathways (linked to parasite load and transmission potential). We identify genetically regulated expression patterns in conserved and variant gene families (e.g., SICAvar, kir), and describe potential transmission from low-density infections using intraerythrocytic deconvolution.

---

### 🧬 Analyses included

- [x] **Whole-genome variant calling** from Illumina WGS using GATK + bcftools consensus
- [x] **Genome-wide association studies (GWAS)** for parasitemia and severity phenotypes
- [x] **Transcriptome assembly** (Trinity-based) from bulk field RNA-seq
- [x] **Transcript abundance quantification** and differential expression analysis (DESeq2, edgeR, voom)
- [x] **Functional annotation** (Trinotate, GO, KEGG, Pfam, Infernal)
- [x] **Stage deconvolution** using Scaden and single-cell references from the Malaria Cell Atlas
- [x] **eQTL mapping** using FastLMM, including cis and trans scans, with integration of GWAS signals
- [x] **Targeted exploratory analysis** of gene families (SICAvar, kir, invasion genes, stress regulators)
- [x] **GO enrichment analysis** using clusterProfiler

### 📂 Directory Structure

```plaintext
├── README.md
├── pipelines/
├── analyses/
│   ├── gwas/
│   ├── rna_seq_qc/
│   ├── annotation/
│   ├── idc_deconvolution/
│   ├── differential_expression/
│   ├── go_enrichment/
│   ├── targeted_expression/
│   └── eqtl_mapping/
```

---

### 🔁 Reproducibility & transparency statement

All analyses were designed with transparency in mind, using modular Snakemake workflows and version-controlled scripts.

However, the repository contains more than 18,000 lines of code and comments, spanning pipelines for variant detection, transcriptomic assembly and quantification, differential expression, eQTL mapping, functional annotation, and figure generation.

So, if you are interested in replicating or adapting any part of this workflow for your own research, feel free to reach out — I am happy to help.

---

### 📄 License

This repository is licensed under the [MIT License](./LICENSE).