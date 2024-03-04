#!/usr/bin/env python3

'''
Get the result from CRISPRDirection v2.0 and save it to a TSV table
for later use.

Juan Camilo Arboleda Rivera
email: juan.arboleda2@udea.edu.co
'''

import pandas as pd
import parse

def create_table():
    # Create new dataframe
    df = pd.DataFrame(columns=['RefSeq ID', 'spacer_start', 'spacer_end',
                               'direction', 'confidence'])

    with open('results.txt') as file:
        for line in file:
            if line.startswith('>'):
                record = parse.search('>{seq_id}:{start:d}-{end:d} ', line)

            if line.startswith('# Final direction:'):
                direction = parse.search('# Final direction:         {direction} [{:g},{}   Confidence: {conf}]', line)
                df.loc[len(df)] = [record['seq_id'], record['start'],
                                   record['end'], direction['direction'],
                                   direction['conf']]

    return df

if __name__ == '__main__':
    print('Running create_table()...')
    df = create_table()
    df.to_csv('sts_transc_direction.tsv', index=False, sep='\t')
