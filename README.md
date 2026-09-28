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
- [dekupl-mergeTags](https://github.com/Transipedia/dekupl-mergeTags)  

All of these, at the versions used here, are provided by the container described below.

## Container

The `Dockerfile` builds an image with all tools at fixed versions (R 4.4.1 with a CRAN snapshot of 2024-10-30; exact R package versions are listed in `R_package_versions.tsv`). Build it for x86-64 and convert it to a Singularity image (SIF) for use on the cluster:

```bash
# on your own machine
docker buildx build --platform linux/amd64 --tag surfr_pipeline:latest --load .
docker save surfr_pipeline:latest -o surfr_pipeline_docker.tar
# copy the .tar to the cluster, then on the cluster
singularity build surfr_pipeline.sif docker-archive://$PWD/surfr_pipeline_docker.tar
```

Run commands in the container with:

```bash
SINGULARITY_TMPDIR=/tmp singularity exec --no-mount bind-paths -B /cfs/klemming surfr_pipeline.sif <command>
```

`run_SURFR_pipeline.sh` already does this through its `SING_EXEC` setting. Both options were needed on the Dardel cluster (PDC, KTH) and may need adjusting elsewhere:

- `SINGULARITY_TMPDIR=/tmp` keeps Singularity's temporary mount point off the Lustre file system, where FUSE mounts of the image are not allowed.
- `--no-mount bind-paths` skips the site-configured bind mounts into `/etc`. With them, Singularity's underlay mode leaves the container's `/etc` empty for this image, and R fails to start (`libblas.so.3` not found).
- `-B /cfs/klemming` makes the project file system available inside the container; replace it with the file system holding your data.

## Pipeline Overview

<img width="230" height="600" alt="SURFR pipeline schematic" src="https://github.com/user-attachments/assets/d2706f55-0a40-4750-92d8-daf084176b24" />

---

### 1. `MergeAndConvertBAMs.sh`

- Uses **Samtools** to merge BAM files for each condition downloaded from the Genomic Data Commons (GDC).
- Merged files are converted to FASTQ format and compressed with **pigz**.  
- If your input files are already in FASTQ format, the conversion step can be skipped.  
- Expected input layout: `<projPath>/Data/<project>/<dataset>/bams/cancer_bams/*.bam` and `.../bams/adjacent_bams/*.bam`. The scripts used for the preprint called the second folder `healthy_bams`; rename it if your data follows that layout.  

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

#### Non-cancer (SRA) k-mer table

Step 4 needs a table of 17-mer counts in the non-cancer public cohort: no header, two whitespace-separated columns (k-mer, count). The cohort consists of the 2,602 public sequencing runs (SRR, ERR and DRR accessions) listed in Supplementary Table 1 (`Supplemental Tables_1.xlsx`).

**This table is not included in the repository, and neither is the code that generated the published version.** To regenerate it, process these runs with the same miRTrace quality control as Step 2 and count strand-specific 17-mers with KMC as in Step 3 (non-canonical k-mers, `-b`).

The table must contain counts for **all** k-mers that pass the TCGA and CPTAC filters, not only for those that pass the non-cancer filter. k-mers absent from the table are given a count of 0 and therefore pass the filter, so an incomplete table silently lets through k-mers that should be removed. `tests/data/sums_SRA_filtered_LUAD.txt` is such a filtered subset (the 156 LUAD k-mers that passed) and is meant only for the regression test.

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

## Verification

What has been checked for this release, and what has not:

**Verified**
- The container builds and runs on a Mac with Docker and on the Dardel cluster with Singularity.
- The non-cancer filter and the dekupl-mergeTags merge in `FindCancerSpecificRNAs.r` reproduce the published LUAD result exactly: the same 73 sequences with the same representative k-mers as Supplementary Table 2 (`tests/test_step4_LUAD_published.R`).
- `FindCancerSpecificRNAs.r` runs end to end, including all plots, on synthetic data with a known answer (`tests/make_test_data.py`).
- The SLURM job chain in `run_SURFR_pipeline.sh` was checked with `--dry-run`.

**Not verified**
- Steps 1–3 were not rerun on the TCGA and CPTAC data for this release. Their commands match the scripts used for the preprint, apart from changes that only affect execution (paths, module loading, overwriting outputs on rerun).
- The enrichment filters and the TCGA–CPTAC intersection in Step 4 were tested on synthetic data only, not rerun on the published data.
- The non-cancer k-mer table is not included (see above), so Step 4 cannot currently be run on the real data from this repository alone.

## Notes on the preprint text

The code in this repository is the reference for how the published results were produced. Where the preprint text is ambiguous or differs, the following applies:

- **Enrichment** is the ratio of raw pooled k-mer counts (cancer / adjacent normal) in each cohort. Although the Methods equation defines the number of tumour and adjacent samples, no normalisation by sample number or library size is applied.
- **Non-cancer filter:** k-mers with 200 or more counts in the non-cancer cohort are removed. No published candidate has exactly 200 counts (the maximum is 199), so this boundary does not affect the published results.
- **Merging:** dekupl-mergeTags is run with k = 17, stranded merging and a minimum overlap of 15 nt (its defaults). Each merged sequence is represented by the k-mer with the lowest non-cancer count, and the count columns reported for a sequence are those of this representative k-mer, not totals over the sequence.
- **LUAD** yields 73 sequences. One of them (AAAATGCTCAGACTCCT) appears in Supplementary Table 2 under LUSC (oncRNA-8), so the table lists 72 sequences under LUAD.
- **samtools:** the container uses samtools 1.23.1 instead of 1.20 cited in the Methods. The release notes between these versions list no change that alters the reads produced by `samtools merge` and `samtools fastq` as used here.

