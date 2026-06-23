import os
import argparse
import functools
import pandas as pd
import hashing
from hashing import HashTable, GeneralizedHash

def kmer2hash(filename, hash_func):
    """ Transform a .txt file with kmers into a hash table"""
    ht = HashTable()
    gen_id = os.path.basename(filename).replace('_kmers.txt', '')

    if not os.path.exists(filename):
        print(f"Warning: File {filename} not found. Skipping.")
        return None

    with open(filename, 'r') as f:
        next(f)
        for line in f:
            kmer = line.strip()
            if not kmer:
                continue
            kmer_hash = GeneralizedHash(kmer, hash_func)
            ht.add(gen_id, kmer_hash)
    ht.sort()

    return ht


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Read k-mer files, create hash tables and write down them into csv files.',
    )

    parser.add_argument('-f', '--kmer_files', nargs='+', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-H', '--hash_function', required=True, help='Hash function to use, provide any of the following: murmur, encode2bit, polynomial_hash, hash_family')
    parser.add_argument('-s', '--seed', nargs='+', type=int, required=False, default=None, help='Enable masked mode, provide a list of positions to exclude from hashing')
    parser.add_argument('-m', '--minimizer', nargs=2, type=int, required=False, default=None, help='Enable minimizer modem requires two integers: the len of a minimizer and the number of minimizers that will be used for hashing')
    parser.add_argument('-p', '--polyhash', type=int, required=False, default=None, help='Turn on hash family mode,take an integer as the number of hash functions to create')
    parser.add_argument('-n', '--iterations', type=int, required=False, default=1, help='Number of iterations for speed benchmarking')
    parser.add_argument('-S', '--global_seed', type=int, default=None, required=False, help='Seed for reproducibility across separate runs')

    args = parser.parse_args()
    if args.global_seed is not None:
        import random
        random.seed(args.global_seed)

    if not os.path.exists(args.output):
        os.makedirs(args.output)

    if args.polyhash is not None:
        functions = list(hashing.universal_hashing(args.polyhash))

        if args.seed is not None:
            base = lambda seq: hashing.hash_family_masked(seq, args.seed, functions)
        else:
            base = lambda seq: hashing.hash_family(seq, functions)

    elif args.seed is not None:
        func_name = args.hash_function if args.hash_function == 'rolling_hash' else args.hash_function + '_masked'
        raw_hash_func = getattr(hashing, func_name)
        base = functools.partial(raw_hash_func, skip_positions=args.seed)

    else:
        base = getattr(hashing, args.hash_function)

    if args.minimizer is not None:
        w_size, num_var = args.minimizer

        if args.polyhash is not None:
            minimizer_base = lambda seq: hashing.hash_family(seq, functions)
        else:
            minimizer_base = getattr(hashing, args.hash_function)

        if args.seed is not None:

            selected_hash_func = lambda seq: hashing.minimizer_hash_masked(
                seq, base_hash=minimizer_base, skip_positions=args.seed, w=w_size, var=num_var
            )
        else:
            selected_hash_func = lambda seq: hashing.minimizer_hash(
                seq, base_hash=minimizer_base, w=w_size, var=num_var
            )
    else:
        selected_hash_func = base

    iters = args.iterations if args.iterations else 1

    for i in range(iters):
        if iters > 1:
            print(f"\n--- Iteration {i + 1}/{iters} ---")

        for input_file in args.kmer_files:
            print(f"Processing: {input_file}...")

            hash_table_obj = kmer2hash(input_file, selected_hash_func)

            if hash_table_obj:
                base_name = os.path.basename(input_file).rsplit('_kmers.txt', 1)[0]

                suffix = f"_n{i + 1}" if iters > 1 else ""
                output_file = os.path.join(args.output, f"{base_name}{suffix}.csv")

                df = pd.DataFrame(list(hash_table_obj.table.items()), columns=['Hash', 'OrganismID'])

                with open(output_file, 'w') as f:
                    f.write(f"# Hash function: {args.hash_function}\n")
                    if iters > 1:
                        f.write(f"# Iteration index: {i + 1}\n")
                    df.to_csv(f, index=False)

                print(f"  Saved to {output_file}")