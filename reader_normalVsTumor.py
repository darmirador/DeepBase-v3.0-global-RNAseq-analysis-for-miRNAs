#!/usr/bin/env python3
import os
import csv
import argparse
from pathlib import Path
import numpy as np
from scipy import stats


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Extract RPM values for all miRNAs from DeepBase metafile CSVs "
            "and classify samples as normal vs cancer, generating one summary "
            "file per cancer type."
        )
    )
    parser.add_argument(
        "-p", "--path",
        default="/Users/darmirador/Desktop/Python Scripts/normalVsCancer_miRNA_bulk",
        help="Directory containing DeepBase miRNA metafile CSVs."
    )
    return parser.parse_args()


def cancer_type_from_filename(filename: str) -> str:
    """
    Extract cancer type from filenames such as:
        metafile_UVM.sampleID.csv -> UVM
        metafile_LUAD.csv -> LUAD
    """
    base = os.path.basename(filename)

    if base.startswith("metafile"):
        base = base[len("metafile"):]
    if base.startswith("_"):
        base = base[1:]

    core = os.path.splitext(base)[0]
    core = core.split(".")[0]

    return core


def significance_from_p(p: float) -> str:
    if np.isnan(p):
        return "NA"
    elif p < 0.0001:
        return "****"
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


def clean_numeric(values):
    cleaned = []
    for x in values:
        try:
            cleaned.append(float(x))
        except (ValueError, TypeError):
            pass
    return np.array(cleaned, dtype=float)


def classify_sample_type(barcode: str):
    """
    TCGA sample type classification from barcode.
    01, 02 = cancer/tumor
    11     = normal
    """
    try:
        tcga_parts = barcode.split("-")
        sample_type = tcga_parts[3][:2]
    except (IndexError, AttributeError):
        return None

    if sample_type in ("01", "02"):
        return "cancer"
    elif sample_type == "11":
        return "normal"
    else:
        return None


def main():
    args = parse_args()
    base_dir = Path(args.path)

    if not base_dir.exists():
        raise FileNotFoundError(f"Directory not found: {base_dir}")

    os.chdir(base_dir)

    meta_files = sorted(
        f for f in os.listdir(base_dir)
        if f.startswith("metafile") and f.endswith(".csv")
    )

    if not meta_files:
        raise FileNotFoundError(f"No metafile CSVs found in: {base_dir}")

    for filename in meta_files:
        file_path = base_dir / filename
        cancer_type = cancer_type_from_filename(filename)
        output_path = base_dir / f"{cancer_type}.csv"

        print(f"Processing {filename} -> {output_path.name}")

        with file_path.open("r", newline="") as csv_file:
            csv_reader = csv.reader(csv_file)

            try:
                header = next(csv_reader)
            except StopIteration:
                print(f"[WARNING] Empty file: {filename}; skipping.")
                continue

            rpm_indices = []
            rpm_barcodes = []

            for i, col in enumerate(header):
                if col.endswith("rpm"):
                    rpm_indices.append(i)
                    rpm_barcodes.append(col)

            if not rpm_indices:
                print(f"[WARNING] No RPM columns found in {filename}; skipping.")
                continue

            sample_classes = [classify_sample_type(barcode) for barcode in rpm_barcodes]

            with output_path.open("w", newline="") as summary_fh:
                summary_writer = csv.writer(summary_fh)
                summary_writer.writerow([
                    "miRNA_ID",
                    "normal_mean",
                    "cancer_mean",
                    "n_normal",
                    "n_cancer",
                    "p_value",
                    "significance_level",
                ])

                for line in csv_reader:
                    if not line:
                        continue

                    mirna_id = line[0].strip()
                    if not mirna_id:
                        continue

                    normal_vals_raw = []
                    cancer_vals_raw = []

                    for idx, sample_class in zip(rpm_indices, sample_classes):
                        if idx >= len(line):
                            continue

                        value = line[idx]

                        if sample_class == "normal":
                            normal_vals_raw.append(value)
                        elif sample_class == "cancer":
                            cancer_vals_raw.append(value)

                    normal_vals = clean_numeric(normal_vals_raw)
                    cancer_vals = clean_numeric(cancer_vals_raw)

                    if len(normal_vals) > 1 and len(cancer_vals) > 1:
                        normal_mean = float(np.mean(normal_vals))
                        cancer_mean = float(np.mean(cancer_vals))
                        _, p_val = stats.ttest_ind(
                            normal_vals,
                            cancer_vals,
                            equal_var=False,
                            nan_policy="omit"
                        )
                        sig_level = significance_from_p(p_val)
                    else:
                        normal_mean = np.nan
                        cancer_mean = np.nan
                        p_val = np.nan
                        sig_level = "NA"

                    summary_writer.writerow([
                        mirna_id,
                        normal_mean,
                        cancer_mean,
                        len(normal_vals),
                        len(cancer_vals),
                        p_val,
                        sig_level,
                    ])

    print("Done.")


if __name__ == "__main__":
    main()