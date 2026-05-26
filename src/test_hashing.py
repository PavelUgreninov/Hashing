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

