#!/usr/bin/env python3
"""
SURFR v1 — synthetic known-answer test data for Step 4 (FindCancerSpecificRNAs.r).

Creates tiny inputs in the exact formats Step 4 expects, with a known answer:
  * SEQ_PASS (25 nt) -> its 9 overlapping 17-mers pass every filter in both
    cohorts and have low control counts. Expected final output: ONE merged
    sequence equal to SEQ_PASS.
  * SEQ_SRA (25 nt)  -> passes both cohorts but is abundant in the controls,
    so its 9 17-mers must be removed by the SRA filter (>= 200).
  * 2,000 random background 17-mers with cancer/adjacent ratio <= 2, so none
    pass the enrichment filter.

This tests that the environment and code run end to end. It does NOT test
the scientific result; that needs the real Step 3 tables.

Usage: python3 make_test_data.py <output_dir>
"""
import os
import random
import sys

K = 17
PROJECT = "TEST"
SEQ_PASS = "TTGCGTGAACCTGAGAATGAGCAGT"  # 25 nt
SEQ_SRA = "GGATCCAAGTCTTAGCCATGCATTA"   # 25 nt

out = sys.argv[1] if len(sys.argv) > 1 else "test_data"
rng = random.Random(42)  # fixed seed: identical data on every run


def kmers(seq):
    return [seq[i:i + K] for i in range(len(seq) - K + 1)]


pass_kmers, sra_kmers = kmers(SEQ_PASS), kmers(SEQ_SRA)
planted = set(pass_kmers) | set(sra_kmers)

background = set()
while len(background) < 2000:
    km = "".join(rng.choice("ACGT") for _ in range(K))
    if km not in planted:
        background.add(km)
background = sorted(background)

merged_dir = os.path.join(out, "merged_tables")
meta_dir = os.path.join(out, "metadata")
os.makedirs(merged_dir, exist_ok=True)
os.makedirs(meta_dir, exist_ok=True)


def write_merged(cohort, planted_rows):
    # Same format as countKmers.sh output: kmer<TAB>cancer<TAB>adjacent, no header,
    # sorted by k-mer, cancer counts >= 30.
    rows = []
    for km in background:
        c = rng.randint(30, 5000)
        a = int(c * rng.uniform(0.5, 2.0))  # ratio <= 2: must fail enrichment
        rows.append((km, c, a))
    rows += planted_rows
    rows.sort()
    path = os.path.join(merged_dir, f"all_{PROJECT}_{cohort}_{K}mers_merged.txt")
    with open(path, "w") as fh:
        for km, c, a in rows:
            fh.write(f"{km}\t{c}\t{a}\n")


# TCGA: enrichment 1000/5 = 200; CPTAC: adjacent 0 -> Inf. Both pass.
write_merged("TCGA", [(km, 1000, 5) for km in pass_kmers + sra_kmers])
write_merged("CPTAC", [(km, 800, 0) for km in pass_kmers + sra_kmers])

# SRA control table: kmer<TAB>count, no header (format assumed by Step 4).
with open(os.path.join(out, "SRA_kmer_counts.txt"), "w") as fh:
    for km in pass_kmers:
        fh.write(f"{km}\t10\n")
    for km in sra_kmers:
        fh.write(f"{km}\t5000\n")
    for km in background[:500]:
        fh.write(f"{km}\t{rng.randint(0, 300)}\n")

# Metadata: TCGA is ';'-separated, CPTAC is space-separated (as Step 4 reads them).
with open(os.path.join(meta_dir, f"clean_metadata_long_TCGA_{PROJECT}.txt"), "w") as fh:
    fh.write("Sample_ID;Sample_Type\n")
    for i in range(20):
        fh.write(f"T{i};{'Tumor' if i < 18 else 'Normal'}\n")
with open(os.path.join(meta_dir, f"clean_metadata_long_CPTAC_{PROJECT}.txt"), "w") as fh:
    fh.write("Sample_ID Sample_Type\n")
    for i in range(20):
        fh.write(f"C{i} {'Tumor' if i < 10 else 'Normal'}\n")

print(f"Test data written to {out}/")
print(f"Expected: 18 k-mers intersected, 9 after SRA filter, 1 merged sequence = {SEQ_PASS}")
