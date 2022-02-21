import sys

from collections import defaultdict
from collections import Counter

from eme_selex.sequence import generate
from eme_selex.sequence import count_kmers
from eme_selex.hashing import build_hashmaps

def calculate_pfm_models(counts, top=50):
    """Calculate PWM models from the counts of all kmers

    To generate a PWM model of a particular kmer, iterate over all the possible mismatches of that kmer and the kmer itself. While iterating, keep adding the nucleotide frequencies at every position from the counts

    Args:
        counts (collections.Counter): kmer occurernce (at
        most one per read) for all kmers

    Kwargs:
        top (int): No. of top PWM models to calculate based
        on the seed kmer occurence (default=50)

    Returns:
        collections.defaultdict(collections.defaultdict(
            collections.Counter)): A nested data structure
            with keys as kmer and values as a defaultdict
            of Counter containing nucleotide frequencies at
            every position along the length of the kmer 

    """
    pfm_models = defaultdict()
    for kmer, count in counts.most_common()[:top]:
        pfm_models[kmer] = defaultdict(Counter)
        mismatch_kmers = generate(kmer)
        mismatch_kmers.add(kmer)
        for mismatch_kmer in mismatch_kmers:
            for pos, nucleotide in enumerate(mismatch_kmer,
                start=1):
                pfm_models[kmer][pos][nucleotide] += counts[mismatch_kmer]

    return pfm_models

def kmer_fraction_from_file(sequence_file, k=5, top=50,
    txt=False):
    """Iterate over the sequencing reads file and calculate
    kmer occurences (at most one per read) and subsequent
    PWM models from the top kmers

    Args:
        sequence_file (str): File containing sequencing
            reads, can be a txt, fasta or fastq file

    Kwargs:
        k (int): length of kmer (default=5)
        top (int): No. of top PWM models to calculate based
        on the seed kmer occurence (default=50)
        txt (bool): Whether the sequene_file is a txt file
            (default=False)

    Returns:
        counts (collections.Counter): kmer occurence (at most
            one per read) for all kmers
        fraction (dict): fraction of reads
            containing the kmer for all kmers
        collections.defaultdict(collections.defaultdict(
            collections.Counter)): A nested data structure
            with keys as kmer and values as a defaultdict
            of Counter containing nucleotide frequencies at
            every position along the length of the kmer
    """    
    if k > 10:
        try:            
            sys.exit()
        finally:
            pass        
            
    kmer_mhash, mhash_kmer, kmer_canonical = build_hashmaps(k=k)

    fraction = {}

    counts, n = count_kmers(sequence_file, kmer_mhash, mhash_kmer, kmer_canonical, k=k, txt=txt)
    for kmer, count in counts.items():
        fraction[kmer] = float(count / n)
    
    pfm_models = calculate_pfm_models(counts, top=top)
    
    return counts, fraction, pfm_models