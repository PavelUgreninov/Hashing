# Hashing Module

This package provides a comprehensive framework for applying and evaluating different hashing approaches to genomic data. It enables the analysis of Jaccard similarity between resulting hash tables, offering an intuitive workflow for benchmarking hashing functions across various configurations.

The module facilitates the validation of hashing strategies in genomic contexts by enabling researchers to:
- **Benchmark** different hashing functions and parameters
- **Assess** the discriminative power of hash-based representations across species
- **Generate** similarity matrices that inform the development of novel approaches
- **Create** training data for annotating metagenomic classifiers

## Scripts Overview

### 1. split.py
Split genome into canonical k-mers of the defined length and write them as a simple `.txt` file.

**Usage:**
```bash
python split.py -k <kmer_size> -f <.fna.gz files> -o <output directory>
```

**Parameters:**
| Parameter | Description |
|-----------|-------------|
| `-k` | K-mer size |
| `-f` | Input `.fna.gz` files |
| `-o` | Output directory |

---

### 2. kmer2hash.py
Apply different hashing functions to k-mers, enabling masking or usage of minimizers. Outputs resulting hash tables as `.csv` files.

**Usage:**
```bash
python kmer2hash.py -H <hash function> -s <space-separated list of skip-positions> -m <minimizer_length #_of_minimizers> -f <.txt files> -o <output directory>
```

**Parameters:**
| Parameter | Description |
|-----------|-------------|
| `-H` | Hashing function (e.g., murmur3, encode2bit, polynomial_hash, hash_family) |
| `-s` | Space-separated list of skip positions |
| `-m` | Minimizer length and number of minimizers |
| `-f` | Input `.txt` files |
| `-p` | Number of universal hash functions (works only with "hashing_family" hashing function) |
| `-S` | Specifies random seed for creating universal hash functions (works only with "hashing_family" hashing function) |
| `-o` | Output directory |

---

### 3. kmer_infer.py
Estimate Jaccard similarity for provided reference and test files. The matrix of Jaccard similarities is written to a `.csv` file.

**Usage:**
```bash
python kmer_infer.py -r <reference tables> -t <test tables> -o <output.csv>
```

**Parameters:**
| Parameter | Description |
|-----------|-------------|
| `-r` | Reference tables |
| `-t` | Test tables |
| `-o` | Output `.csv` file |

---

## Example Workflow

1. Split genomes into k-mers:
   ```bash
   python split.py -k 30 -f genome.fna.gz -o kmers/
   ```

2. Generate hash tables:
   ```bash
   python kmer2hash.py -H murmur3 -s 7 14 21 -m 25 1 -f kmers/*.txt -o hash_tables/
   ```

3. Calculate Jaccard similarity:
   ```bash
   python kmer_infer.py -r hash_tables/reference.csv -t hash_tables/test.csv -o results.csv
   ```

## License

Licensed under the [MIT License](LICENSE)
