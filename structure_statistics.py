#!/usr/bin/env python2

'''
Generates numerical summaries of batch MFE foldings. Depending on the mode it either...
R)Summarizes multiple directories of MFE <.ct> files into a single <.csv>.
C)Runs 'scorer' between multiple directories of MFE <.ct> files, generating a pairwise PPV
for the same oligo in different folds.
'''

#Imports
import os
import argparse
from sf2libs.connectivity_table import collect_connectivity_tables
from sf2libs.structure_io import check_extension

#Functions
def collect_connectivity_dirs(list_of_directories,offset):
    '''Applys the CT_Report table to directories, gathers information'''
    mass_data = {}
    for directory in list_of_directories:
        raw_dict = collect_connectivity_tables(directory)
        #We want to collapse entries around transcript between directories
        #Some transcript nomenclatures have underscores in them, thus an offset can work around it
        melted_dict = dict([('_'.join(k.split('_')[:1+offset]),v) for k, v in raw_dict.items()])  
        mass_data[os.path.basename(directory)] = melted_dict
    return mass_data

def write_out_raw(raw_data,out_fyle='derp.csv'):
    '''Writes out raw data'''
    vertical_keys = sorted(list(set.union(*map(set, raw_data.values()))))
    horizontal_keys = sorted(raw_data.keys())
    cols_type = ['double_stranded','single_stranded','DeltaG']
    cols_labels = [item+'_'+affix for item in horizontal_keys for affix in cols_type]
    with open(out_fyle,'w') as g:
        header = ','.join(['transcript']+cols_labels)
        g.write(header+'\n')
        for v in vertical_keys:
            garble_line = [[raw_data[h][v].perc_stranded,raw_data[h][v].perc_unstranded,raw_data[h][v].energy] if v in raw_data[h] else ['NA']*3 for h in horizontal_keys]
            flat_line = [str(item) for sublist in garble_line for item in sublist]
            outline = ','.join([v]+flat_line)
            g.write(outline+'\n')

def main():
    parser = argparse.ArgumentParser(description='Summarizes or compares MFE <.ct> files')
    parser.add_argument('-d',type=str,help='CT directories', nargs='+')
    parser.add_argument('-mode',type=str.upper,default= None,choices = ['R','C'],help='Raw statistics or Comparative Analysis')
    parser.add_argument('-name',type=str,default = None, help = 'Output file name')
    parser.add_argument('-offset',type=int,default = 0, help = 'Number of Underscores in Transcript Names')
    args = parser.parse_args()
    
    if args.mode =='R':
        data = collect_connectivity_dirs(args.d,args.offset)
        default_name = '_'.join(sorted(data.keys())+['statisics'])+'.csv'
        out_name = default_name if args.name == None else check_extension(args.outname,'.csv')
        write_out_raw(data,out_name)
    
    if args.mode == 'C':
        print 'Mode Not Yet Implemented, Call Again'

if __name__ == '__main__':
    main()
