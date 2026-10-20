from kraken_alt_hash import hashing
from functools import partial
import pytest

from kraken_alt_hash.hashing import encode2bit, GeneralizedHash, HashTable


class TestHashing:
    """Check the initialization and for equality of two GeneralizedHash objects"""
    seq1 = "ACTGGTCAGCAGCTACGACGACACGACTACGACGG"
    seq2 = "ACTGGTCAGCAACTACTTACATCGACTACGACTTG"

    @pytest.mark.parametrize("function", [hashing.murmur3, partial(hashing.murmur3_masked, skip_positions=(7,)),
                                          hashing.encode2bit, partial(hashing.encode2bit_masked, skip_positions=(7,)),
                                          hashing.polynomial_hash, partial(hashing.polynomial_hash_masked, skip_positions=(7,)),
                                          partial(hashing.hash_family, funcs=hashing.universal_hashing(5)),
                                          partial(hashing.hash_family_masked, funcs=hashing.universal_hashing(5), skip_positions=(7,))])
    def test_generalized_hash(self, function):
        if function == hashing.hash_family or function == hashing.hash_family_masked:
            self.seq1 = hashing.encode2bit(self.seq1)
            self.seq2 = hashing.encode2bit(self.seq2)

        hash1 = hashing.GeneralizedHash(self.seq1, function)
        hash2 = hashing.GeneralizedHash(self.seq2, function)

        # test initialization
        assert isinstance(hash1, GeneralizedHash) and isinstance(hash2, GeneralizedHash)
        assert hasattr(hash1, "value") and  hasattr(hash2, "value")
        assert hasattr(hash1, "seq") and hasattr(hash2, "seq")
        assert hasattr(hash1, "hash_func") and hasattr(hash2, "hash_func")

        #test for equality of hashes
        result = hash1.value == hash2.value
        assert not result


class TestHashtable:
    seq1 = "ACTGGTCAGCAGCTACGACGACACGACTACGACGG"
    seq2 = "ACTGGTCAGCAACTACTTACATCGACTACGACTTG"
    seq3 = "ACTGGTCAGCAGTTACGACGACACGACTACGACGG"
    seq4 = "ACTGACCAGCAACTACTTACATCGACTACGACTTG"

    def test_hash_table(self):
        ht = hashing.HashTable()
        ht2 = hashing.HashTable()

        #check initialization
        assert isinstance(ht, HashTable)
        assert hasattr(ht, "table")

        #check that len method works correctly
        assert len(ht) == 0

        #check that add method works correctly
        ht.add('seq1', hashing.GeneralizedHash(self.seq1, hashing.murmur3))
        ht.add('seq2', hashing.GeneralizedHash(self.seq3, hashing.murmur3))
        assert len(ht) == 2
        ht2.add('seq2', hashing.GeneralizedHash(self.seq2, hashing.murmur3))
        ht2.add('seq4', hashing.GeneralizedHash(self.seq4, hashing.murmur3))
        assert len(ht2) == 2

        # check that merge works
        ht.merge(ht2)
        assert len(ht) == 4

        #check that sort works
        ht.sort()
        assert list(ht.table.keys()) == sorted(ht.table.keys())



# ---------------------------------------------------------------- regression tests for the fixes
import gzip
import random
from kraken_alt_hash import kmer_split


def _write_fna(path, records):
    with gzip.open(path, "wt") as f:
        for i, seq in enumerate(records):
            f.write(f">r{i}\n")
            for j in range(0, len(seq), 60):
                f.write(seq[j:j + 60] + "\n")


def _random_dna(n, seed=0):
    rnd = random.Random(seed)
    return "".join(rnd.choice("ACGT") for _ in range(n))


class TestFixes:
    def test_generalized_hash_eq_compares_other(self):
        a = hashing.GeneralizedHash("x", lambda s: (1, 2, 3))
        b = hashing.GeneralizedHash("y", lambda s: (7, 8, 9))
        assert not (a == b)          # used to be True (min(self) == min(self))
        assert a == hashing.GeneralizedHash("z", lambda s: (1, 5, 6))

    @pytest.mark.parametrize("k,gap", [(21, None), (21, (1, 1)), (21, (2, 1)), (22, (3, 1)), (20, (2, 1))])
    def test_spaced_kmers_strand_symmetric(self, tmp_path, k, gap):
        g = _random_dna(3000)
        _write_fna(tmp_path / "a.fna.gz", [g])
        _write_fna(tmp_path / "b.fna.gz", [kmer_split.rc(g)])
        A = set(kmer_split.split_kmers(str(tmp_path / "a.fna.gz"), k, gap))
        B = set(kmer_split.split_kmers(str(tmp_path / "b.fna.gz"), k, gap))
        assert A == B and len(A) > 0

    def test_ambiguous_bases_skipped(self, tmp_path):
        seq = "ACGTACGTAC" + "N" + "GTACGTACGT"
        _write_fna(tmp_path / "n.fna.gz", [seq])
        kmers = list(kmer_split.split_kmers(str(tmp_path / "n.fna.gz"), 5))
        assert all(set(x) <= set("ACGT") for x in kmers)
        assert len(kmers) == 2 * (10 - 5 + 1)   # windows overlapping the N are dropped

    def test_records_do_not_overlap(self, tmp_path):
        _write_fna(tmp_path / "r.fna.gz", ["ACGTA", "CGTAC"])
        assert len(list(kmer_split.split_kmers(str(tmp_path / "r.fna.gz"), 5))) == 2

    @pytest.mark.parametrize("func", [hashing.murmur3_masked, hashing.encode2bit_masked, hashing.polynomial_hash_masked])
    def test_masked_hash_strand_symmetric(self, func):
        s = _random_dna(30, seed=3)
        skip = (3, 7, 12)
        assert func(s, skip_positions=skip) == func(kmer_split.rc(s), skip_positions=skip)

    def test_masked_minimizer_runs_and_orders(self):
        s = _random_dna(30, seed=5)
        for order in ("lex", "hash"):
            v1 = hashing.minimizer_hash(s, 10, hashing.encode2bit, var=2, order=order)
            v2 = hashing.minimizer_hash_masked(s, 10, [2, 5], hashing.encode2bit, var=2, order=order)
            assert len(v1) == 2 and len(v2) == 2
        with pytest.raises(ValueError):
            hashing.minimizer_hash(s, 10, hashing.encode2bit, order="nope")
