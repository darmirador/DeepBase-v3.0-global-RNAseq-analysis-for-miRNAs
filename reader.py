#!/usr/bin/env python3
import os
import csv
import argparse
from pathlib import Path
import numpy as np
from scipy import stats


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract RPM values for a specific miRNA from DeepBase metafiles "
                    "and classify samples as normal vs cancer."
    )
    parser.add_argument(
        "-p", "--path",
        default="/Users/darmirador/Desktop/Python Scripts/normalVsCancer_miRNA",
        help="Directory containing DeepBase miRNA metafile CSVs."
    )
    parser.add_argument(
        "-m", "--mirna",
        required=True,
        help="miRNA name to extract, e.g. 'hsa-mir-490'."
    )
    return parser.parse_args()


def cancer_type_from_filename(filename: str) -> str:
    base = os.path.basename(filename)
    if base.startswith("metafile"):
        base = base[len("metafile"):]
    if base.startswith("_"):
        base = base[1:]
    return os.path.splitext(base)[0]


def significance_from_p(p: float) -> str:
    if p < 0.0001:
        return "****"
    elif p < 0.001:
        return "***"
    elif p < 0.01:
        return "**"
    elif p < 0.05:
        return "*"
    else:
        return "ns"


def main():
    args = parse_args()
    pathname = args.path
    mirnaName = args.mirna

    base_dir = Path(pathname)
    os.chdir(pathname)

    meta_files = [f for f in os.listdir(pathname) if f.startswith("metafile")]

    # SUMMARY FILE (per miRNA)
    summary_path = base_dir / f"{mirnaName}.csv"
    with summary_path.open("w", newline="") as summary_fh:
        summary_writer = csv.writer(summary_fh)
        summary_writer.writerow(
            [
                "cancer_type",
                "normal-mean",
                "cancer-mean",
                "n_normal",
                "n_cancer",
                "p-value",
                "significance_level",
            ]
        )

        # PROCESS EACH METAFILE
        for filename in meta_files:

            # --- READ & EXTRACT DATA --- #
            tcgaId = []
            expressionLevel = []
            tcgaIdSubset = []
            expressionLevelSubset = []

            with open(filename, "r") as csvFile:
                csvReader = csv.reader(csvFile)
                tcgaId = next(csvReader)

                for line in csvReader:
                    if line[0] == mirnaName:
                        expressionLevel = line
                        break

            if not expressionLevel:
                print(f"[WARNING] {mirnaName} not in {filename}; skipping.")
                continue

            # RPM columns only
            for i in range(len(tcgaId)):
                if tcgaId[i].endswith("rpm"):
                    tcgaIdSubset.append(tcgaId[i])
                    expressionLevelSubset.append(expressionLevel[i])

            class tissueType:
                def __init__(self):
                    self.sampleId = []
                    self.normalizedTranscripts = []

            normalSample = tissueType()
            cancerSample = tissueType()

            counter = 0
            for barcode in tcgaIdSubset:
                tcgaString = barcode.split("-")
                sampleTypeString = tcgaString[3][:2]

                if sampleTypeString in ("01", "02"):
                    cancerSample.sampleId.append(barcode[:28])
                    cancerSample.normalizedTranscripts.append(expressionLevelSubset[counter])
                elif sampleTypeString == "11":
                    normalSample.sampleId.append(barcode[:28])
                    normalSample.normalizedTranscripts.append(expressionLevelSubset[counter])
                counter += 1

            # --- PER-CANCER OUTPUT CSV --- #
            cancer_type = cancer_type_from_filename(filename)
            out_name = f"{mirnaName}_{cancer_type}.csv"    # << NEW FORMAT HERE
            out_path = base_dir / out_name

            with out_path.open("w") as csvOutput:
                csvOutput.write("tissue,barcode,rpm\n")

            with out_path.open("a") as csvOutput:
                for i in range(len(normalSample.sampleId)):
                    csvOutput.write(
                        f"normal,{normalSample.sampleId[i]},{normalSample.normalizedTranscripts[i]}\n"
                    )
                for i in range(len(cancerSample.sampleId)):
                    csvOutput.write(
                        f"cancer,{cancerSample.sampleId[i]},{cancerSample.normalizedTranscripts[i]}\n"
                    )

            print(f"Saved per-cancer file: {out_name}")

            # --- SUMMARY STATS --- #
            def clean(arr):
                vals = []
                for x in arr:
                    try:
                        vals.append(float(x))
                    except:
                        pass
                return np.array(vals)

            normal_vals = clean(normalSample.normalizedTranscripts)
            cancer_vals = clean(cancerSample.normalizedTranscripts)

            if len(normal_vals) > 1 and len(cancer_vals) > 1:
                normal_mean = float(np.mean(normal_vals))
                cancer_mean = float(np.mean(cancer_vals))
                t_stat, p_val = stats.ttest_ind(normal_vals, cancer_vals, equal_var=False)
                sig_level = significance_from_p(p_val)
            else:
                normal_mean = np.nan
                cancer_mean = np.nan
                p_val = np.nan
                sig_level = "NA"

            summary_writer.writerow(
                [
                    cancer_type,
                    normal_mean,
                    cancer_mean,
                    len(normal_vals),
                    len(cancer_vals),
                    p_val,
                    sig_level,
                ]
            )

    print(f"Done generating summary file: {mirnaName}.csv")


if __name__ == "__main__":
    main()
