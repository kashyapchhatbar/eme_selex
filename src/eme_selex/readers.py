from eme_selex.decorators import txt
from eme_selex.decorators import fastx

@txt
def read_txt(txt_file):
    """Iterates over a text file. If the file name ends in
    "gz", the decorator automatically uses the gzip module.

    Args:
        txt_file (str): Sequencing reads in text format

    Yields:
        str: One sequencing read at a time (uppercase)
    """    
    for line in txt_file:
        yield line.strip().upper()
    txt_file.close()

@fastx
def read_fastx(fastx_file):
    """Iterates over a fasta/q file using Biopython's
    "SeqIO.parse" iterator. If the file name ends in
    "gz", the decorator automatically uses the gzip module.

    Args:
        fastx_file (str): Sequencing reads in fasta/q format

    Yields:
        Bio.SeqRecord.SeqRecord: A sequence with annotation
    """
    for rec in fastx_file:
        yield str(rec.seq)


def file_handle(sequence_file, txt=False):
    """Returns the relevant file handle based on the
    filename

    Args:
        sequence_file (str): File containing sequencing
            reads, can be a txt, fasta or fastq file

    Kwargs:
        txt (bool): Whether the sequene_file is a txt file
            (default=False)
    """
    try:    
        if txt:
            handle = read_txt(sequence_file)
        else:
            handle = read_fastx(sequence_file)

        return handle
    except Exception:
        console.print_exception()
        raise