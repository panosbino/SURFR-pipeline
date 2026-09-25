# Sequence-based Unbiased Reference-Free small-RNA sequencing (SURFR) pipeline

This repository provides scripts for running the **SURFR pipeline**, including tools to identify cancer-specific small RNAs as described in *Kalogeropoulos et al., 2025 (preprint)*.  

## Software Requirements

We recommend running the SURFR pipeline on a computational cluster with a job scheduling system such as **SLURM** for analysis of very large datasets such us the ones used in the *Kalogeropoulos et al., 2025 (preprint)* manuscript. For smaller datasets the pipeline can be also run locally.   

The following software is required:  
- [Samtools](http://www.htslib.org/)  
- [pigz](https://zlib.net/pigz/)  
- [miRTrace](https://github.com/friedlanderlab/mirtrace/)  
- [KMC](https://github.com/refresh-bio/KMC)  
- [R](https://www.r-project.org/)  

## Pipeline Overview

<img width="230" height="600" alt="SURFR pipeline schematic" src="https://github.com/user-attachments/assets/d2706f55-0a40-4750-92d8-daf084176b24" />

---

### 1. `MergeAndConvertBAMs.sh`

- Uses **Samtools** to merge BAM files for each condition downloaded from the Genomic Data Commons (GDC).
- Merged files are converted to FASTQ format and compressed with **pigz**.  
- If your input files are already in FASTQ format, the conversion step can be skipped.  

---

### 2. `QCandFasta.sh`

- Runs quality control on the FASTQ files with **miRTrace**.  
- Converts FASTQ to uncollapsed FASTA format.  

---

### 3. `countKmers.sh`

- Uses **KMC** to generate k-mers from the FASTA files.  
- Filters out k-mers with fewer than 30 counts in the cancer condition to reduce noise.  
- Produces a k-mer count table by left-joining cancer and adjacent sample tables.  

---

### 4. `FindCancerSpecificRNAs.r`

Identifies cancer-specific sequences for each cohort using the following criteria:  
- **Cancer counts > 200**  
- **Cancer enrichment > 40**  
- **Adjacent counts < 100**  

The workflow continues as follows:  
1. Cancer-specific sequences are identified independently in TCGA and CPTAC cohorts.  
2. Overlapping sequences are intersected across the two cohorts.  
3. Expression of overlapping k-mers is quantified in an independent cohort of non-cancer samples.  
4. k-mers with 200 or more counts in controls are removed.  
5. Overlapping and offset k-mers are merged using **dekupl-mergeTags** (k = 17, stranded, minimum overlap 15 nt).  

---

## Testing

Two tests run inside the SURFR container, from the repository root:

```bash
# 1. Synthetic known-answer test of Step 4 (environment and code run end to end)
python3 tests/make_test_data.py test_data
Rscript FindCancerSpecificRNAs.r TEST test_data/merged_tables test_data/SRA_kmer_counts.txt \
    test_data/metadata test_data/results /opt/dekupl/bin/mergeTags
# expected: 1 merged sequence, TTGCGTGAACCTGAGAATGAGCAGT

# 2. Regression test against the published LUAD result
Rscript tests/test_step4_LUAD_published.R /opt/dekupl/bin/mergeTags
# expected: PASS, 73 sequences and representative k-mers identical to Supplementary Table 2
```

## Notes on the preprint text

The code in this repository is the reference for how the published results were produced. Where the preprint text is ambiguous or differs, the following applies:

- **Enrichment** is the ratio of raw pooled k-mer counts (cancer / adjacent normal) in each cohort. Although the Methods equation defines the number of tumour and adjacent samples, no normalisation by sample number or library size is applied.
- **Non-cancer filter:** k-mers with 200 or more counts in the non-cancer cohort are removed. No published candidate has exactly 200 counts (the maximum is 199), so this boundary does not affect the published results.
- **Merging:** dekupl-mergeTags is run with k = 17, stranded merging and a minimum overlap of 15 nt (its defaults). Each merged sequence is represented by the k-mer with the lowest non-cancer count, and the count columns reported for a sequence are those of this representative k-mer, not totals over the sequence.
- **LUAD** yields 73 sequences. One of them (AAAATGCTCAGACTCCT) appears in Supplementary Table 2 under LUSC (oncRNA-8), so the table lists 72 sequences under LUAD.
- **samtools:** the container uses samtools 1.23.1 instead of 1.20 cited in the Methods. The release notes between these versions list no change that alters the reads produced by `samtools merge` and `samtools fastq` as used here.

