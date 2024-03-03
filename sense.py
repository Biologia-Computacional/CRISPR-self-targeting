#!/usr/bin/env python3

import pandas as pd
from Bio import Entrez, SeqIO
from sys import stderr # For printing to stderr
from utils import *

Entrez.email = 'juan.arboleda2@udea.edu.co'

def get_strand(seq_record):
    '''Get sequence strand for genes.
    If the region is not a gene returns None.
    
    Parameters:
    seq_record: A Bio.SeqRecord.SeqRecord obtained by importing a Genbank file
    '''

    for feature in seq_record.features:
        if feature.type == 'gene':
            return feature.strand
    return None

def get_ps_sense(df):
    '''Gets information about Protospacer sense of transcription
    from NCBI nuccore.'''

    results_list = []
    n = 0 # Counter
    total = len(df)

    # Get CDS info
    for i, start, end in zip(df['Refseq ID'], df['Proto-spacer Start'], df['Proto-spacer End']):
        n += 1
        if (n % 200) == 0:
            print('%2d%% complete.' % (n/total*100))

        try:
            cds = get_seq_record(i, start, end, 'gb')
            strand = get_strand(cds)
            results_list.append(strand)
        except Exception as e:
            print(e, 'Refseq ID:', i, 'Start:', start, 'End:', end, file=stderr)
            results_list.append('Error: '+str(e))
            continue

    return results_list

def get_orientation(seq_1, seq_2) -> int:
    '''Returns
    ---

    0: Sequences do not match in any orientation.
    1: Sequences are in the same strand and direction.
    2: Sequences are in the same strand but different direction.
    3: Sequences are in a different strand but same direction.
    4: Sequences are in a different strand and direction.'''

    # Ensure all sequences are uppercase
    seq_1 = seq_1.upper()
    seq_2 = seq_2.upper()

    if seq_1 == seq_2:
        return 1
    if seq_1 == seq_2[::-1]:
        return 2
    if seq_1 == seq_2.reverse_complement():
        return 3
    if seq_1 == seq_2.complement():
        return 4
    else:
        return 0

def complement(c, td_1, td_2):
    '''Determine if RNA transcribed from sequences seq_1 and seq_2
    are complementary.

    Arguments:
     - c: Case (0,1,2,3 or 4) according to get_orientation
     - td_1: Transcription direction of sequence seq_1
     - td_2: Transcription direction of sequence seq_2'''

    assert (type(td_1) == int) or (type(td_1) == float), 'td_1 is not int or float'
    assert (type(td_2) == int) or (type(td_2) == float), 'td_2 is not int or float'

    if (c == 2) or (c == 4):
        return 'Not complementary'
    if (c == 1):
        if ((td_1 == 1) and (td_2 == -1)) or ((td_1 == -1) and (td_2 == 1)):
            return 'Complementary'
        else:
            return 'Complementary to non-coding strand'
    if (c == 3):
        if ((td_1 == 1) and (td_2 == 1)) or ((td_1 == -1) and (td_2 == -1)):
            return 'Complementary'
        else:
            return 'Complementary to non-coding strand'
    else:
        return 'No info'

def complement_from_seq(seq_1, seq_2, td_1, td_2):
    '''Determine if RNA transcribed from sequences seq_1 and seq_2
    are complementary.

    Arguments:
     - td_1: Transcription direction of sequence seq_1
     - td_2: Transcription direction of sequence seq_2'''

    c = compare_seqs(seq_1, seq_2)
    return complement(c, td_1, td_2)

def orientations(df):
    '''Process dataframe df and gets the corresponding orientation of
    spacer with respect to the protospacer (0,1,2,3 or 4) according to
    get_orientation.'''

    result_list = []
    n = 0 # Counter
    total = len(df)

    for seq_id, ps_start, ps_end, s_start, s_end in zip(df['Refseq ID'],
                                                   df['Proto-spacer Start'], df['Proto-spacer End'],
                                                   df['Spacer Start'], df['Spacer End']):
        n += 1
        if (n % 200) == 0:
            print('%2d%% complete.' % (n/total*100))

        try:
            ps_rec = get_seq_record(seq_id, ps_start, ps_end, 'fasta')
            s_rec = get_seq_record(seq_id, s_start, s_end, 'fasta')
            result_list.append(get_orientation(ps_rec.seq, s_rec.seq))
        except Exception as e:
            print(e, 'Refseq ID:', seq_id, 'Start:', s_start, 'End:', s_end, file=stderr)
            result_list.append('Error: '+str(e))
            continue

    return result_list

def spacers2fasta(df, file, offset=150):
    '''Gets the sequences of spacers within an offset and save it to file in
    FASTA fromat.'''

    n = 0 # Counter
    total = len(df)
    
    with open(file, 'a') as output_handle:
        for seqid, start, end in zip(df['Refseq ID'], df['Spacer Start'], df['Spacer End']):
            n += 1
            if (n % 200) == 0:
                print('%2d%% complete.' % (n/total*100))
            try:
                spacer = get_seq_record(seqid, start-offset, end+offset)
                spacer.id = seqid + ':' + str(start) + '-' + str(end)
                SeqIO.write(spacer, output_handle, "fasta")
            except Exception as e:
                print(e, 'Refseq ID:', seqid, 'Start:', start, 'End:', end, file=stderr)
                continue

def download_spacers(df, file, offset=150):
    '''Get sequences of spacers within an offset and save it to file in
    FASTA format'''
    
    ids = list(df['Refseq ID'])
    starts = [i-offset for i in df['Spacer Start']]
    ends = [i+offset for i in df['Spacer End']]
    
    try:
        spacers = get_seq(ids, starts, ends, 'fasta')
    except Exception as e:
        print(e, file=stderr)
    
    with open(file, 'w') as out_handle:
        out_handle.write(spacers)
    print('Created file:', file)

if __name__ == '__main__':
    # Load Self-Targeting database
    df = pd.read_csv('self-target-proteins.tsv', sep='\t')
    
    # Protospacer sense
    if 'Protospacer_sense' not in df.columns:
        print('Running get_ps_sense(df)...')
        df['Protospacer_sense'] = get_ps_sense(df)

        # Save table to file
        print('Saving results to file...')
        df.to_csv('self-target-proteins.tsv', index=False, sep='\t')

    # Relative orientation of protospacers and spacers
    if 'Orientations' not in df.columns:
        print('Running orientations(df)...')
        df['Orientations'] = orientations(df)

        # Save table to file
        print('Saving results to file...')
        df.to_csv('self-target-proteins.tsv', index=False, sep='\t')

    # Download spacer context sequences
    print('Downloading spacers...')
    spacers2fasta(df, 'spacers.fasta', 200)
