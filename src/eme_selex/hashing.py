import mmh3
import itertools

from Bio import Seq
from eme_selex.sequence import generate_all_kmers
from eme_selex.sequence import canonical

def hash_kmer(kmer, seed=21):
    """Calculates murmurhash with keyword argument seed
    (int: default=21). 
    """
    rc_kmer = Seq.reverse_complement(kmer)

    # determine whether original k-mer or reverse complement is lesser
    if kmer < rc_kmer:
        canonical_kmer = kmer
    else:
        canonical_kmer = rc_kmer

    # calculate murmurhash
    mhash = mmh3.hash64(canonical_kmer, seed)[0]
    if mhash < 0: mhash += 2**64

    # done
    return mhash

def build_hashmaps(k=5):
    """Create hashmaps of canonical kmers and murmur hash
    (mhash)

    Keyword Args:
        k (int): default=5 length of the kmer

    Returns:
        dict: hashmap with kmer as key and mhash as value
        dict: hashmap with mhash as value and tuple of
            canonical and non-canonical kmer as value
        dict: hashmap with kmer as key and its canonical
            kmer as value
    """
    kmer_mhash, mhash_kmer, kmer_canonical = {}, {}, {}
    for i in generate_all_kmers(k=k):
        kmer_mhash[i] = hash_kmer(i)
    
    for i, j in kmer_mhash.items():
        canonical_kmer, other_kmer = canonical(i)
        mhash_kmer[j] = (canonical_kmer, other_kmer)
        kmer_canonical[i] = canonical_kmer

    return kmer_mhash, mhash_kmer, kmer_canonical