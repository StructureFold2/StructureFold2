#!/usr/bin/env python

'''
Takes two directories of CT files, reports gains and losses of stems between two conditions.
Still experimental. Only works on MFE CT files. 
'''

#Imports
import os
import glob
import argparse
from ct_tools import *
from collections import Counter

#Functions
def delta_loops(directory1,directory2):
    '''Extracts all loops from two directories, returns losses and gains of loops for each parallel transcript'''
    all_losses,all_gains = Counter(),Counter()
    control = collect_connectivity_tables(directory1)
    experimental = collect_connectivity_tables(directory2)
    common_keys = set(control.keys()).intersection(set(experimental.keys()))
    for key in common_keys:
        con,exp = control[key],experimental[key]
        con2,exp2 = Counter(con.primary_loops()),Counter(exp.primary_loops())
        gains,losses = exp2 - con2, con2 - exp2
        all_losses += losses
        all_gains += gains
    return all_losses,all_gains

def delta_stems(directory1,directory2):
    '''Extracts all stems from two directories, returns losses and gains of stems for each parallel transcript'''
    all_losses,all_gains = Counter(),Counter()
    control = collect_connectivity_tables(directory1)
    experimental = collect_connectivity_tables(directory2)
    common_keys = set(control.keys()).intersection(set(experimental.keys()))
    for key in common_keys:
        con,exp = control[key],experimental[key]
        con2,exp2 = Counter(con.primary_stems()),Counter(exp.primary_stems())
        gains,losses = exp2 - con2, con2 - exp2
        all_losses += losses
        all_gains += gains
    return all_losses,all_gains

def sequence_dinucleotide_freq(seq):
    '''Simple dinucelotide'''
    return Counter([seq[n-1:n+1] for n in range(1,len(seq))])

def write_fasta(info,outfyle='out.fasta',LW=80):
    '''Writes out the <.fasta> file, names are just transcript+step'''
    with open(outfyle,'w') as g:
        for name,seq in sorted(info.items()):
            g.write('>' + name + '\n')
            for i in xrange(0,len(seq),LW):
                g.write(seq[i:i+LW] + '\n')

def components_to_fasta_dict(acounter,prefix='Component',threshold=10):
    '''Takes a counter, makes a fancy out dict for writing (MEMEs!)'''
    new = {}
    over_thresh = dict([(x,b) for x,b in acounter.items() if len(x)>= threshold])
    for i, (k, v) in enumerate(sorted(over_thresh.items()),1):
        for sub_index in range(1,v+1):
            new['_'.join([prefix,str(i),str(sub_index)])] = k
    return new

def main():
    parser = argparse.ArgumentParser(description='Pulls loop/stem gains/losses between to CT directories')
    parser.add_argument('control',type=str,help='Control Directory')
    parser.add_argument('experimental',type=str,help='Experimental Directory')
    parser.add_argument('-flen',type=int,default=10,help='Feature Minimum Length')
    parser.add_argument('-prefix',type=str,default='CT_Compare', help='Prefix for out files')
    args = parser.parse_args()
    #Workflow
    lost_loops,gained_loops = delta_loops(args.control,args.experimental)
    lost_stems,gained_stems = delta_stems(args.control,args.experimental)
    data_sets = [lost_loops,gained_loops,lost_stems,gained_stems]
    suffixes = ['lost_loops','gained_loops','lost_stems','gained_stems']
    for data, suffix in zip(data_sets,suffixes):
        temp = components_to_fasta_dict(data,suffix,args.flen)
        write_fasta(temp,'_'.join([args.prefix,suffix,str(args.flen)+'min'])+'.fa')

if __name__ == '__main__':
    main()
