# Normal vs Cancer miRNA Bulk Analysis

This project is designed to compare miRNA expression between normal tissue and tumor samples across multiple cancer types using TCGA-style sample barcodes and DeepBase-derived expression data.

## Goal

The main goal is to identify miRNAs that are differentially expressed between normal samples and cancer samples for each cancer type. The analysis produces one summary table per cancer type, reporting:

- mean miRNA expression in normal samples
- mean miRNA expression in cancer samples
- number of normal samples used
- number of cancer samples used
- p-value from a statistical comparison
- significance label (e.g. *, **, ***, ****, or ns)

## How it works

The script reads files named like:

- metafile_LUAD.sampleID.csv
- metafile_BRCA.sampleID.csv
- metafile_UVM.sampleID.csv

Each metafile contains miRNA RPM values across samples. The script:

1. Identifies RPM columns in the DeepBase metadata files.
2. Parses the TCGA barcode from each sample to classify it as either:
   - normal (sample type 11)
   - cancer/tumor (sample types 01 or 02)
3. For each miRNA, calculates mean expression in normal and cancer groups.
4. Compares groups using a Welch two-sample t-test.
5. Writes a CSV file for each cancer type summarizing the results.

## Input data

This workflow expects DeepBase miRNA expression files in the same folder as the script. The dataset is organized by cancer type, with metafile files corresponding to each disease cohort.

The project also includes output CSV files such as:

- LUAD.csv
- BRCA.csv
- UVM.csv
- etc.

These files are generated from the corresponding metafile inputs and contain a per-miRNA comparison between normal and tumor samples.

## Script

The main script is:

- reader_normalVsTumor.py

It can be run from the project directory with:

```bash
python reader_normalVsTumor.py
```

or with an explicit path:

```bash
python reader_normalVsTumor.py --path "/Users/darmirador/Desktop/Python Scripts/normalVsCancer_miRNA_bulk"
```

## Output format

Each cancer-specific output CSV has the following columns:

- miRNA_ID
- normal_mean
- cancer_mean
- n_normal
- n_cancer
- p_value
- significance_level

The significance column uses standard asterisk notation:

- **** for p < 0.0001
- *** for p < 0.001
- ** for p < 0.01
- * for p < 0.05
- ns for non-significant results

## Why this analysis matters

This type of comparison helps identify miRNAs whose expression is consistently altered in tumor tissue compared with normal tissue. Such miRNAs may be useful as:

- biomarkers for cancer detection
- markers of disease subtype or progression
- candidates for further mechanistic or validation studies

## Notes

- This is a bulk-expression comparison and is intended for exploratory analysis.
- Statistical testing is performed only when both normal and cancer groups have enough valid numeric values.
- The project is focused on miRNA expression profiles from TCGA-associated cancer datasets and DeepBase-derived RPM values.
