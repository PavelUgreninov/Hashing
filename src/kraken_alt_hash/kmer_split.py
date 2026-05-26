import argparse
import gzip
import os

def rc(seq):
    """ Create the reverse complement of seq """
    DNA_TRANS = str.maketrans("ACGT", "TGCA")
    return seq.translate(DNA_TRANS)[::-1]


def split_kmers(fasta_file, k):
    """ Return the generator of lexicographically smaller kmers """
    with gzip.open(fasta_file, "rt") as f:
        next(f)
        window = ""
        for line in f:
            line_seq = line.strip()
            if not line_seq or line_seq.startswith(">"):
                continue

            window += line_seq

            while len(window) >= k:
                kmer = window[:k]
                yield min(kmer, rc(kmer))
                window = window[1:]


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Split fasta files into canonical k-mers'
    )

    parser.add_argument('-f', '--fasta_files', nargs='+', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-k', '--kmer_size', type=int, required=True)

    args = parser.parse_args()

    kmer_size = args.kmer_size

    for fasta_path in args.fasta_files:

        # Check that fasta files exist
        if not os.path.exists(fasta_path):
            print(f"Warning: File {fasta_path} not found. Skipping.")
            continue

        print(f"Processing: {fasta_path}...")

        # Extract filename and combine with the output directory
        base_name = os.path.basename(fasta_path).split('.fna.gz')[0]
        output_name = os.path.join(args.output, f"{base_name}_kmers.txt")

        # Check that the output directory exists, otherwise, create it
        if not os.path.exists(args.output):
            os.makedirs(args.output)

        # Write down the k-mer size and every k-mers each on a separate line
        with open(output_name, 'w') as f:
            f.write(f"K-mer size: {kmer_size}\n")
            for kmer in split_kmers(fasta_path, kmer_size):
                f.write(kmer + '\n')

        print(f"  Saved to {output_name}")