import argparse
import gzip
import os


def rc(seq):
    """ Create the reverse complement of seq """
    DNA_TRANS = str.maketrans("ACGT", "TGCA")
    return seq.translate(DNA_TRANS)[::-1]


def spaced_offsets(k, gap=None):
    """
    Return the list of positions (offsets) sampled from a window.

    gap = (n_bases, gap_size):
      (1, 1) -> 1010101...
      (2, 1) -> 110110110...
    k is the number of sampled nucleotides (the weight of the pattern),
    so the window span is larger than k when gaps are used.
    If gap is None, a contiguous k-mer is used.
    """
    if gap is None:
        return list(range(k))

    n_bases, gap_size = gap
    if n_bases < 1 or gap_size < 0:
        raise ValueError("-g needs: bases >= 1 and gap size >= 0")

    unit = n_bases + gap_size
    offsets = []
    pos = 0
    while len(offsets) < k:
        if pos % unit < n_bases:
            offsets.append(pos)
        pos += 1
    return offsets


def split_kmers(fasta_file, k, gap=None):
    """ Return the generator of lexicographically smaller (spaced) kmers """
    offsets = spaced_offsets(k, gap)
    span = offsets[-1] + 1  # length of the window covered by the pattern

    with gzip.open(fasta_file, "rt") as f:
        window = ""
        for line in f:
            line_seq = line.strip().upper()
            if not line_seq:
                continue
            if line_seq.startswith(">"):
                window = ""  # new record: don't build k-mers across records
                continue

            window += line_seq

            start = 0
            while start + span <= len(window):
                kmer = "".join(window[start + i] for i in offsets)
                yield min(kmer, rc(kmer))
                start += 1

            window = window[start:]  # keep the last span-1 characters


if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        description='Split fasta files into canonical (optionally spaced) k-mers'
    )

    parser.add_argument('-f', '--fasta_files', nargs='+', required=True)
    parser.add_argument('-o', '--output', required=True)
    parser.add_argument('-k', '--kmer_size', type=int, required=True,
                        help='Number of sampled nucleotides in a k-mer')
    parser.add_argument('-g', '--gap', nargs=2, type=int,
                        metavar=('BASES', 'GAP'), default=None,
                        help='Spaced pattern: BASES nucleotides followed by a '
                             'gap of GAP positions. E.g. "-g 1 1" -> 1010101..., '
                             '"-g 2 1" -> 110110110...')

    args = parser.parse_args()

    kmer_size = args.kmer_size
    gap = tuple(args.gap) if args.gap else None

    # Validate pattern early
    offsets = spaced_offsets(kmer_size, gap)

    for fasta_path in args.fasta_files:

        if not os.path.exists(fasta_path):
            print(f"Warning: File {fasta_path} not found. Skipping.")
            continue

        print(f"Processing: {fasta_path}...")

        base_name = os.path.basename(fasta_path).split('.fna.gz')[0]
        output_name = os.path.join(args.output, f"{base_name}_kmers.txt")

        if not os.path.exists(args.output):
            os.makedirs(args.output)

        with open(output_name, 'w') as f:
            f.write(f"K-mer size: {kmer_size}\n")
            if gap:
                pattern = "".join("1" if i in set(offsets) else "0"
                                  for i in range(offsets[-1] + 1))
                f.write(f"Spaced pattern: {pattern}\n")
            for kmer in split_kmers(fasta_path, kmer_size, gap):
                f.write(kmer + '\n')

        print(f"  Saved to {output_name}")
