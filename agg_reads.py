import argparse
import sys
import pandas as pd

exe_parser = argparse.ArgumentParser()
exe_parser.add_argument('tsv_file', type=str, help='<input file>')

args = exe_parser.parse_args()

if args.tsv_file:
    tsv_file = args.tsv_file
else:
    sys.exit('Exiting. No tsv file specified.')

tsv_file = '/Users/jfloresu/analysis/tim_explorer/ERR594413.tsv'
tsv_file_agg = tsv_file.replace('.tsv','.agg.tsv')

df = pd.read_csv(tsv_file,sep='\t',names=['read','genome','id','gene'])
del df['genome']

df['id'] = df.loc[:,'id'].round(0)
df['dataset'] = df['read'].apply(lambda x: x.split('.')[0])

gdf = df.groupby(['dataset','gene','id']).agg('count').reset_index()

gdf.to_csv(tsv_file_agg, header=False, sep='\t', index=False)
