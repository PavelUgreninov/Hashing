import csv
import os
import argparse
import pandas as pd

def load_hashes_from_csv(filepath):

    hashes = set()
    with open(filepath, "r", encoding="utf-8") as f:
        next(f, None)
        reader = csv.reader(f, skipinitialspace=True)
        for row in reader:
            if row:
                hashes.add(row[0])
    return hashes


def calculate_jaccard(set1: set, set2: set):
    """Calculates the Jaccard Similarity index between two sets."""
    intersection = len(set1.intersection(set2))
    union = len(set1.union(set2))
    return intersection / union if union > 0 else 0.0


def run_paired_inference(output_file, ref_files, test_files):
    if not ref_files or not test_files:
        print("Error: Both reference files and test files must be provided.")
        return

    ref_labels = [os.path.basename(f) for f in ref_files]
    test_labels = [os.path.basename(f) for f in test_files]

    heatmap = pd.DataFrame(index=test_labels, columns=ref_labels, dtype=float)

    print(f"Computing Jaccard for {len(test_files)}x{len(ref_files)} matrix...")

    try:
        test_cache = {
            os.path.basename(f): load_hashes_from_csv(f) for f in test_files
        }
        ref_cache = {
            os.path.basename(f): load_hashes_from_csv(f) for f in ref_files
        }
    except FileNotFoundError as e:
        print(f"Error reading file: {e}")
        return

    for test_name in test_labels:
        for ref_name in ref_labels:
            heatmap.loc[test_name, ref_name] = calculate_jaccard(
                test_cache[test_name],
                ref_cache[ref_name]
            )

    heatmap.to_csv(output_file)
    print(f"Stats saved to {output_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Jaccard distance matrix calculation')

    parser.add_argument('-r', '--references', nargs='+', required=True,
                        help='List of reference CSV files')
    parser.add_argument('-t', '--tests', nargs='+', required=True,
                        help='List of test CSV files')
    parser.add_argument('-o', '--output_name', required=True,
                        help='Output CSV filepath')

    args = parser.parse_args()

    run_paired_inference(
        output_file=args.output_name,
        ref_files=args.references,
        test_files=args.tests
    )