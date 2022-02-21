import gzip
import functools

from Bio import SeqIO

def txt(func):    
    @functools.wraps(func)
    def return_handle(*args, **kwargs):
        if args[0].endswith("gz"): 
            file_handle = gzip.open(args[0], "rt")
        else:
            file_handle = open(args[0], "rt")
        return func(file_handle)        
    return return_handle

def fastx(func):
    @functools.wraps(func)
    def return_handle(*args, **kwargs):
        if args[0].endswith("gz"):
            gzip_handle = gzip.open(args[0], "rt")
            if args[0].endswith("fasta.gz") or args[0].endswith("fa.gz"):
                file_handle = SeqIO.parse(gzip_handle,
                    "fasta")
            elif args[0].endswith("fastq.gz") or args[0].endswith("fq.gz"):
                file_handle = SeqIO.parse(gzip_handle,
                    "fastq")
        elif args[0].endswith("fa") or args[0].endswith("fasta"):
            _file_handle = open(args[0], "rt")
            file_handle = SeqIO.parse(_file_handle, "fasta")
        elif args[0].endswith("fastq") or args[0].endswith("fq"):
            _file_handle = open(args[0], "rt")
            file_handle = SeqIO.parse(_file_handle,
                "fastq")
        return func(file_handle)
    return return_handle