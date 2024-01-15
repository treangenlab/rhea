## function for counting the number of queries with >80% alignment
## used to estimate the percent of reads included in curated MAGs

import argparse
import pandas as pd


parser = argparse.ArgumentParser(description='')
parser.add_argument('minimap_SAM', help='Path to minimap sam output')
args = parser.parse_args()
COLUMN_HEADERS = ['query', 'q_len', 'flag', 'align_len', 'direction', 'ref', '1', '2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12']

df = pd.read_csv(args.minimap_SAM, sep='\t', names=COLUMN_HEADERS)
align_percent_df = df[['query', 'q_len', 'align_len']]
align_percent_df['percent_align'] = align_percent_df['align_len']/align_percent_df['q_len']
align_percent_df_grouped = align_percent_df.groupby('query').max()
count_greater_than_80p = (align_percent_df_grouped['percent_align'] > 0.8).sum()
print(count_greater_than_80p)