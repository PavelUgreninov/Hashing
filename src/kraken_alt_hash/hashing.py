from collections import defaultdict
from collections.abc import Iterable
from typing import Callable
import random
import mmh3
import numpy as np

class GeneralizedHash:
    """ Create an object of generalized hash """

    def __init__(self, seq: str, hash_func: Callable):
        self.seq = seq
        self.hash_func = hash_func
        self.value = hash_func(seq)

    def __eq__(self, other):
        if not isinstance(other, GeneralizedHash):
            return False
        if isinstance(self.value, Iterable) and isinstance(other.value, Iterable):
            return min(self.value) == min(other.value)
        return self.value == other.value

            # correct! for example, for list - compare all / min / etc depend on self.hash_func
    def __len__(self):
        return len(self.value) if hasattr(self.value, "__len__") else None

class HashTable:
    """ Create an object of generalized hash table """

    def __init__(self, probing = None):
        self.table = defaultdict(list)

    def add(self, gen_id, generalized_hash):
            self.table[generalized_hash.value].append(gen_id)

    def __len__(self):
        return len(self.table)

    def sort(self):
        self.table = dict(sorted(self.table.items()))

    def merge(self, other):
        for key, gen_ids in other.table.items():
            self.table[key].extend(gen_ids)



def rc(seq):
    """ Create the reverse complement of seq """
    DNA_TRANS = str.maketrans("ACGT", "TGCA")
    return seq.translate(DNA_TRANS)[::-1]

def masked_canonical(seq, skip_positions):
    """
    Remove skip_positions from seq and from its reverse complement and return
    the smaller result. This makes a masked k-mer a true spaced seed that is
    identical for both strands (masking an already canonicalised k-mer gives a
    pattern that flips from k-mer to k-mer).
    """
    skip = set(skip_positions)
    fwd = "".join(c for i, c in enumerate(seq) if i not in skip)
    rev = "".join(c for i, c in enumerate(rc(seq)) if i not in skip)
    return fwd if fwd <= rev else rev


lut = np.zeros(256, dtype=np.uint64)
lut[ord('A')] = 0
lut[ord('C')] = 1
lut[ord('G')] = 2
lut[ord('T')] = 3


def encode2bit(seq):
    arr = lut[np.frombuffer(seq.encode('ascii'), dtype=np.uint8)]
    k_effective = len(arr)

    if k_effective == 0: return 0

    shifts = np.arange(k_effective - 1, -1, -1, dtype=np.uint64)
    powers = np.left_shift(1, 2 * shifts)  # Equivalent to 4**shifts

    return int(np.dot(arr, powers))

def encode2bit_masked(seq, skip_positions=None):
    if skip_positions:
        seq = masked_canonical(seq, skip_positions)
    return encode2bit(seq)

lut_14 = np.zeros(256, dtype=np.uint64)
lut_14[ord('A')], lut_14[ord('C')], lut_14[ord('G')], lut_14[ord('T')] = 1, 2, 3, 4


def polynomial_hash(seq, p=5, m=(1 << 61) - 1):
    arr = lut_14[np.frombuffer(seq.encode('ascii'), dtype=np.uint8)]
    hash_val = 0
    for val in reversed(arr.tolist()):
        hash_val = (hash_val * p + val) % m

    return hash_val


def polynomial_hash_masked(seq, skip_positions=None, p=5, m=(1 << 61) - 1):
    if skip_positions:
        seq = masked_canonical(seq, skip_positions)
    return polynomial_hash(seq, p=p, m=m)

def murmur3(data, seed=0):
    res = mmh3.hash(data)
    return res

def murmur3_masked(data, skip_positions, seed=0):
    return mmh3.hash(masked_canonical(data, skip_positions))


def _window_order(order):
    """Sort key for choosing minimizers: plain lexicographic (biased towards
    A-rich l-mers) or the murmur3 hash of the l-mer (pseudo-random order)."""
    if order == "lex":
        return None
    if order == "hash":
        return lambda s: mmh3.hash(s, signed=False)
    raise ValueError("order must be 'lex' or 'hash'")


def minimizer_hash(seq: str, w: int, base_hash, var=1, order="lex"):
    num_windows = len(seq) - w + 1
    if num_windows <= 0:
        return None if var == 1 else tuple()

    windows = [seq[i:i + w] for i in range(num_windows)]
    windows.sort(key=_window_order(order))

    if var == 1:
        return base_hash(windows[0])

    chosen = windows[:var] if isinstance(var, int) else windows

    return tuple(base_hash(kmer) for kmer in chosen)

def minimizer_hash_masked(seq: str, w: int, skip_positions: list, base_hash, var=1, order="lex"):
    num_windows = len(seq) - w + 1
    if num_windows <= 0:
        return None if var == 1 else tuple()
    windows = [seq[i:i + w] for i in range(num_windows)]
    windows.sort(key=_window_order(order))

    # Select minimizers
    if var == 1:
        chosen = [windows[0]]
    else:
        chosen = windows[:var] if isinstance(var, int) else windows

    results = []

    for kmer in chosen:
        results.append(base_hash(masked_canonical(kmer, skip_positions)))

    return results[0] if var == 1 else tuple(results)

def universal_hashing(num_functions, m=2 ** 32):
    p = (1 << 61) - 1
    hash_functions = []

    for _ in range(num_functions):
        a = random.randint(1, p - 1)
        if a % 2 == 0:
            a += 1
        b = random.randint(0, p - 1)
        h = lambda x, a=a, b=b: ((a * x + b) % p) % m
        hash_functions.append(h)

    return hash_functions

def hash_family(seq, funcs):
    seq = encode2bit(seq)
    return tuple(f(seq) for f in funcs)

def hash_family_masked(seq, skip_positions, funcs:list):
    seq = encode2bit_masked(seq, skip_positions)
    return tuple(f(seq) for f in funcs)
