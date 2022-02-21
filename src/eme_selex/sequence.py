import itertools
import random

from Bio import Seq
from collections import Counter
from eme_selex.readers import file_handle

mismatch_orderings = ["".join(i) for i in itertools.permutations(['A', 'T', 'G', 'C'], 4)]

def canonical(kmer):
    """Returns the lesser of the original k-mer and its 
    reverse complement as a tuple

    Args:
        kmer (str): Sequence kmer

    Returns:
        str: canonical kmer
        str: reverse complement of the canonical kmer
    """
    rc_kmer = Seq.reverse_complement(kmer)

    return min(kmer, rc_kmer), max(kmer, rc_kmer)

def sequence_kmers(sequence, kmer_canonical, k):
    """Returns a set of canonical kmers of length k from
    the sequence

    Args:
        sequence (str): Sequencing read
        kmer_canonical (dict): Hashmap of kmer as key
            and its canonical kmer as the value
        k (int): Length of kmer

    Returns:
        set: Seq of canonical kmers present in the read
    """
    sequence_kmer_set = set()
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]        
        sequence_kmer_set.add(kmer_canonical[kmer])

    return sequence_kmer_set

def count_kmers(sequence_file, kmer_mhash, mhash_kmer, kmer_canonical, k=5, txt=False):
    """Counts kmers in the sequencing reads and returns a Counter object

    Args:
        sequence_file (str): File containing sequencing
            reads, can be a txt, fasta or fastq file
        kmer_mhash (dict): hashmap with kmer as key and 
            mhash as value
        mhash_kmer (dict): hashmap with mhash as value 
            and tuple of canonical and non-canonical 
            kmer as value
        kmer_canonical (dict): hashmap with kmer as key
            and its canonical kmer as value

    Kwargs:
        k (int): length of kmer (default=5)
        txt (bool): Whether the sequene_file is a txt file
            (default=False)

    Returns:
        counts (collections.Counter): kmer occurernce (at
            most one per read) for all kmers
        n (int): number of sequences in the sequence_file
    """
    counts = Counter()
    n = 0
    for seq in file_handle(sequence_file, txt=txt):
        n += 1
        for kmer in sequence_kmers(seq, kmer_canonical, k):
            mhash = kmer_mhash[kmer]
            if mhash_kmer[mhash][0] == mhash_kmer[mhash][1]:
                counts[mhash_kmer[mhash][0]] += 1
            else:
                counts[mhash_kmer[mhash][0]] += 1
                counts[mhash_kmer[mhash][1]] += 1
    return counts, n

def generate(s):
    """Generate all possible mismatches of a kmer
    """
    mismatch_kmers = set()
    range_N = list(range(len(s)))
    random.shuffle(range_N)
    letters = random.choice(mismatch_orderings)

    for index in range_N:
        nucleotides = ['A', 'T', 'G', 'C']
        nucleotides.remove(s[index])
        for nucleotide in nucleotides:
            mismatched_s = list(s)
            mismatched_s[index] = nucleotide
            mismatch_kmers.add("".join(mismatched_s))
    return mismatch_kmers

def generate_all_kmers(k=5):
    """Returns a generator with all the possible kmers of
    length k

    Keyword Args:
        k (int): default=5 length of the kmer

    Yields:
        str: kmer
    """
    for i in itertools.product(['A','T','G','C'], repeat=k):
        yield ''.join(i)

def hamming_distance(top_kmer, other_kmer):
    """Returns the hamming distance between two kmers. Because we use the canonical representation of a kmer, this function calculates hamming distances between all combinations of kmer and its reverse complement and returns the minimum value

    Args:
        kmer (str): Most occuring kmer
        kmer (str): Another kmer

    Returns:
        int: Hamming distance between the two kmers
    """
    one = sum(c1 != c2 for c1, c2 in zip(top_kmer,
        other_kmer))
    two = sum(c1 != c2 for c1, c2 in zip(top_kmer,
        Seq.reverse_complement(other_kmer)))
    three = sum(c1 != c2 for c1, c2 in zip(
        Seq.reverse_complement(top_kmer),
        Seq.reverse_complement(other_kmer)))
    four = sum(c1 != c2 for c1, c2 in zip(
        Seq.reverse_complement(top_kmer),
        other_kmer))
    return min(one, two, three, four)