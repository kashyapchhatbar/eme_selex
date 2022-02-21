import pytest
import os
import random

from eme_selex.sequence import canonical
from eme_selex.sequence import sequence_kmers
from eme_selex.sequence import generate
from eme_selex.sequence import count_kmers
from eme_selex.hashing import build_hashmaps
from testfixtures import TempDirectory

data = b"""ATGCA
ATTGC
GCCAT
TATAT
GATTA
ACGTA
GCATG
CATGA
"""

data2 = b"""GGATCCAAGCTT
AAGCTTAAGCTT
GCATGCGCATGC
AAGCTTGGATCC
GCATGGACGTAC
"""

def test_canonical_upper_case():
	assert canonical("ATATT") == ("AATAT", "ATATT")

def test_canonical_lower_case():
	assert canonical("atatt") == ("aatat", "atatt")

def test_sequence_kmers():
	example_sequence = "ATATTA"
	kmer_canonical = {"ATA": "ATA",	"TAT": "ATA", "ATT": "AAT", "TTA": "TAA"}
	assert sequence_kmers(example_sequence, kmer_canonical, 3) == set(["ATA", "AAT", "TAA"])

def test_generate():
	example_kmer = "ATA"
	assert generate(example_kmer) == set(["TTA", "CTA", "GTA", "AAA", "AGA", "ACA", "ATG", "ATC", "ATT"])

def test_count_kmers():
	a = TempDirectory()
	a.write('test.txt', data)	
	kmer_mhash, mhash_kmer, kmer_canonical = build_hashmaps(k=3)
	counts, _ = count_kmers(os.path.join(a.path, 'test.txt'), 
		kmer_mhash, mhash_kmer, kmer_canonical, k=3, txt=True)
	a.cleanup()
	assert counts['ATT'] == 2
	assert counts['AAT'] == 2

def test_count_kmers_6():
	a = TempDirectory()
	a.write('test2.txt', data2)	
	kmer_mhash, mhash_kmer, kmer_canonical = build_hashmaps(k=6)
	counts, _ = count_kmers(os.path.join(a.path, 'test2.txt'), 
		kmer_mhash, mhash_kmer, kmer_canonical, k=6, txt=True)
	a.cleanup()
	assert counts['GGATCC'] == 2
	assert counts['AAGCTT'] == 3