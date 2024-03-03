#!/usr/bin/env python3

from Bio import Entrez, SeqIO
from io import StringIO # Needed to "transform" a string into a file

Entrez.email = 'juan.arboleda2@udea.edu.co'

def get_seq(id, start, stop, ret_type='fasta'):
    '''Get sequence as FASTA string'''

    handle = Entrez.efetch(db='nuccore', id=id, seq_start=start, seq_stop=stop, rettype=ret_type)
    seq = handle.read()
    handle.close()

    return seq

def get_seq_record(id, start, stop, ret_type='fasta'):
    '''Get a SeqRecord from nuccore'''

    seq = get_seq(id, start, stop, ret_type)

    if ret_type == 'fasta':
        return SeqIO.read(StringIO(seq), 'fasta')
    elif ret_type == 'gb':
        return SeqIO.read(StringIO(seq), 'gb')

