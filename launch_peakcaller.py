import numpy as np 
import pickle
import argparse 
import subprocess
import glob
import os
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('chr_map_pickle',type=str)
parser.add_argument('hnrnpk_for_path',type=str)
args = parser.parse_args()

chr_map = pickle.load(open(args.chr_map_pickle,'rb'))
hnrnpk_for = pd.read_csv(args.hnrnpk_for_path,header=None,sep='\t')

for i,row in hnrnpk_for.iterrows():
    peak_range = set(range(row[1],row[2]))
    for query,vals in chr_map[row[0]].items():
        strand = vals[0]['strand']
        if strand=='-':
            pass
        elif strand=='+':  
            query_range = set(vals[0]['seq_range'])
            if peak_range.issubset(query_range):
                if len(vals)>1:
                    if peak_range.intersection(set(vals[1]['seq_range'])):
                        vals[1]['peaks'] += 1
                    else:
                        vals[0]['peaks'] += 1
                else:
                    vals[0]['peaks']+=1


pickle.dump(chr_map,open('for_hnrnpk_chrmap_out.p','wb'))