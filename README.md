# Sequence-based Unbiased Reference-Free small-RNA sequencing (SURFR) pipeline

This repository provides scripts for running the **SURFR pipeline**, including tools to identify cancer-specific small RNAs as described in *Kalogeropoulos et al., 2025 (preprint)*.  

## Software Requirements

We recommend running the SURFR pipeline on a computational cluster with a job scheduling system such as **SLURM**.  

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
- ⚠️ If your input files are already in FASTQ format, the conversion step can be skipped.  

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
4. k-mers with >200 counts in controls are removed.  
5. Overlapping and offset k-mers are merged using **dekupl-mergeTags**.  

---
