import subprocess
import time
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob
import os
import argparse
import pandas as pd


def run_stage(command):
    """Executes a shell command and returns the wall-clock duration."""
    start = time.perf_counter()
    subprocess.run(command, shell=True, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
    return time.perf_counter() - start


def main():
    parser = argparse.ArgumentParser(description='K-mer Benchmarking Suite for GitHub Repo')
    parser.add_argument('-k', '--kmer', type=int, required=True, help='K-mer size')
    parser.add_argument('-H', '--hash_funcs', nargs='+', required=True,
                        help='List of hash functions to test (e.g., python_hash minimizer_hash)')
    parser.add_argument('-r', '--ref_pattern', required=True, help='Pattern for Ref files')
    parser.add_argument('-t', '--test_pattern', required=True, help='Pattern for Test files')
    parser.add_argument('-o', '--output_dir', default="./bench_results", help='Output directory')
    parser.add_argument('-n', '--reps', type=int, default=5, help='Repetitions for inference')
    parser.add_argument('--extra_args', type=str, default="", help='Additional args like "-m 26 1"')

    args = parser.parse_args()
    os.makedirs(args.output_dir, exist_ok=True)
    os.makedirs("./temp_data", exist_ok=True)

    ref_files = sorted(glob.glob(args.ref_pattern))
    test_files = sorted(glob.glob(args.test_pattern))

    if not ref_files or not test_files:
        print("Error: Files not found.")
        return

    barplot_data = []

    heatmap_results = []

    for h_func in args.hash_funcs:
        print(f"\n--- Benchmarking Hash Type: {h_func} ---")

        table_cache = {}
        build_times = []

        def get_table_with_timing(fpath):
            if fpath in table_cache:
                return table_cache[fpath]

            base = os.path.basename(fpath).split('.')[0]
            kmer_file = f"./temp_data/{base}.kmers"
            table_file = f"./temp_data/{base}.table"

            run_stage(f"python3 kmer_split.py -k {args.kmer} -f {fpath} -o ./temp_data")

            cmd = f"python3 kmer2hash.py -H {h_func} -f {kmer_file} {args.extra_args} -o ./temp_data"
            b_time = run_stage(cmd)
            build_times.append(b_time)

            table_cache[fpath] = table_file
            return table_file

        compare_times = []

        for r_path in ref_files:
            r_table = get_table_with_timing(r_path)
            r_name = os.path.basename(r_path)

            for t_path in test_files:
                t_table = get_table_with_timing(t_path)
                t_name = os.path.basename(t_path)

                print(f"  Inference: {r_name} vs {t_name}")
                pair_times = []
                for _ in range(args.reps):
                    c_time = run_stage(f"python3 kmer_infer.py {r_table} {t_table}")
                    pair_times.append(c_time)

                avg_pair_time = np.mean(pair_times)
                compare_times.append(avg_pair_time)

                # Store for the heatmap (only for the current hash function)
                heatmap_results.append({
                    "HashFunction": h_func,
                    "Reference": r_name,
                    "Test": t_name,
                    "Time": avg_pair_time
                })

        # Log average Build and Compare times for this Hash Type
        barplot_data.append({"HashType": h_func, "Phase": "Build Time", "Time": np.mean(build_times)})
        barplot_data.append({"HashType": h_func, "Phase": "Compare Time", "Time": np.mean(compare_times)})

    # --- 1. GENERATE GROUPED BARPLOT ---
    bench_df = pd.DataFrame(barplot_data)
    plt.figure(figsize=(10, 6))
    sns.set_style("whitegrid")
    ax = sns.barplot(data=bench_df, x='HashType', y='Time', hue='Phase', palette='viridis')

    # Annotate bars
    for p in ax.patches:
        if p.get_height() > 0:
            ax.annotate(f'{p.get_height():.4f}', (p.get_x() + p.get_width() / 2., p.get_height()),
                        ha='center', va='center', xytext=(0, 10), textcoords='offset points', fontsize=8)

    plt.title(f"Build vs Compare Performance (K={args.kmer})")
    plt.ylabel("Time (seconds)")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "benchmark_barplot.png"))

    # --- 2. GENERATE HEATMAP (For the last hash function in the list) ---
    h_df = pd.DataFrame(heatmap_results)
    last_h = args.hash_funcs[-1]
    pivot_df = h_df[h_df['HashFunction'] == last_h].pivot(index="Test", columns="Reference", values="Time")

    plt.figure(figsize=(12, 8))
    sns.heatmap(pivot_df, annot=True, fmt=".4f", cmap='magma')
    plt.title(f"Inference Time Heatmap (Rows: Test, Cols: Ref)\nHash: {last_h}")
    plt.tight_layout()
    plt.savefig(os.path.join(args.output_dir, "inference_heatmap.png"))

    # --- 3. SAVE LOG (.out) ---
    with open(os.path.join(args.output_dir, "summary.out"), "w") as f:
        f.write("K-mer Hashing Benchmark Results\n")
        f.write("=" * 30 + "\n")
        f.write(bench_df.to_string(index=False))

    print(f"\nDone! Outputs saved to {args.output_dir}")


if __name__ == "__main__":
    main()