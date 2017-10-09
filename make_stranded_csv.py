#!/usr/bin/env python
'''
Generates metrics related to strandeness from a directory of <.ct> files (connectivity tables).
'''
#Imports
from collections import Counter
import argparse
import glob
import os

#Functions
def get_all_strandedness(CT_directory='CT'):
    '''Returns a nested dictionary, True (Stranded) or False (Single Stranded)'''
    os.chdir(CT_directory)
    all_data = {}
    for ct_fyle in glob.glob('*.ct'):
        with open(ct_fyle,'r') as f:
             big_key = f.readline().strip().split()[1]
             for line in f:
                 pos, base, neg, pos, pair, nat = line.strip().split()
                 all_data.setdefault(big_key, {})[pos] = False if pair == '0' else True
    os.chdir('..')
    return all_data

def write_out_stranded_info(info_dict,outfyle='derp.csv',suffix=None):
    '''Writes out the info'''
    with open(outfyle,'w') as g:
        if suffix == None:
            header = ','.join(['transcript','double_bases','double_percent','single_bases','single_percent'])+'\n'
        else:
            header = ','.join(['transcript','double_bases'+'_'+suffix,'double_percent'+'_'+suffix,'single_bases'+'_'+suffix,'single_percent'+'_'+suffix])+'\n'
        g.write(header)
        for k, v in info_dict.items():
            stranded_counts = Counter(v.values())
            ds,ss,total = float(stranded_counts[True]),float(stranded_counts[False]),float(sum(stranded_counts.values()))
            try:
                x_values = ','.join([k]+[str(x) for x in [ds,float(ds)/total,ss,float(ss)/total]])+'\n'
            except:
                x_values = ','.join([k]+['NA']*4)+'\n'
            g.write(x_values)

def main():
    parser = argparse.ArgumentParser(description='Creates a <.csv> of the amount of strandeness per transcript from a directory of <.ct> files')
    parser.add_argument('directory',type=str,help='Path to CT directory')
    parser.add_argument('-name',default='strands.csv', help='Specify output file name')
    parser.add_argument('-suffix',default=None, help='Append given suffix to data column names')
    args = parser.parse_args()
    
    #Read in the data
    data = get_all_strandedness(args.directory)
    
    #Write data
    write_out_stranded_info(data,args.name,args.suffix)

if __name__ == '__main__':
    main()
