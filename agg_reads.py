#!/usr/bin/env python

import argparse
import sys
import pandas as pd

exe_parser = argparse.ArgumentParser()
exe_parser.add_argument('tsv_file', type=str, help='<input file>')
exe_parser.add_argument('--mode', type=str, help='output mode: cov - output tsv for cov plot \n reg - regular tsv output for map.')
exe_parser.add_argument('--metadata', type=str, help="TARA Oceans metadata file, only required for cov mode.")

args = exe_parser.parse_args()

if args.tsv_file:
	tsv_file = args.tsv_file
else:
	sys.exit('Exiting. No tsv file specified.')

if args.mode in {'cov','reg'}:
	mode = args.mode
else:
	sys.exit('Exiting. No tsv file mode specified.')

if args.mode == "reg":
	if args.metadata == None:
		sys.exit('Exiting. No metadata tsv file specified.')

df = pd.read_csv(tsv_file,sep='\t',names=['read','genome','id','gene','aln_start','aln_end','dataset'])
del df['genome']

df['id'] = df.loc[:,'id'].round(0)

if mode == 'reg':
	del df['aln_start']
	del df['aln_end']
	
	tsv_file_output = tsv_file.replace('.tsv','.agg.tsv')	

	gdf = df.groupby(['dataset','gene','id']).agg('count').reset_index()
	gdf.to_csv(tsv_file_output, header=False, sep='\t', index=False)
else:
	tsv_file_output = tsv_file.replace('.tsv','.cov.tsv')
	metadata_df = pd.read_csv(args.metadata)
	metadata_df.set_index('ENA-RUN',inplace=True)
	metadata_df = metadata_df.rename(columns={'Sampling Station':'Station', "Longitude Start":"Longitude", "Latitude Start":"Latitude"})

	gdf = df.groupby(['dataset','gene','id','aln_start','aln_end']).agg('count').reset_index()
	
	gdf.loc[:,'ENA-RUN'] = gdf.loc[:,'dataset'].copy()
	gdf.set_index('ENA-RUN', inplace=True)
	gdf = pd.concat([gdf,metadata_df[["Station","Longitude","Latitude","Depth","Protocol Label"]]],axis=1,join_axes=[gdf.index])

	viral_criterion_fraction = gdf['Protocol Label'].map(lambda x: x in {'GIRUS_NUC-dry_W0.1-0.22','VIRUS_NUC-DNA-Fe(20L)_W<-0.22','VIRUS_NUC-Fe_W<-0.22'})
	gdf.loc[viral_criterion_fraction,'Fraction'] = 'VIRUS'
	gdf.loc[viral_criterion_fraction == False,'Fraction'] = 'BACT'

	del gdf["Protocol Label"]
    
	gdf.reset_index(drop=True, inplace=True)
    
	del gdf['dataset']
	
	gdf.to_csv(tsv_file_output, index=False, header=False, sep='\t')
