#!/usr/bin/env python2

#Imports
import argparse
from itertools import islice
from collections import Counter
import glob

comp_dict = {"A":"T","T":"A","C":"G","G":"C","N":"N"}

def report_lig(infyle,outfyle):
    '''Checks .fastqs for patterns'''
    first_base_counts,total_counts = Counter(),Counter()
    with open(infyle, 'r') as f, open(outfyle, 'w') as g:
        while True:
            next_n_lines = list(islice(f, 4))
            if not next_n_lines:
                break
            read_hash,seq,plus,qual = [n.strip() for n in next_n_lines]
            total_counter = Counter(seq)
            total_counts = total_counts + total_counter
            first_base_counts[seq[0]]+=1
        g.write(','.join(['nt','count','position'])+'\n')
        for k,v in first_base_counts.items():
            g.write(','.join([comp_dict[k],str(v),'position_1'])+'\n')
        for k,v in total_counts.items():
            g.write(','.join([k,str(v),'total_counts'])+'\n')

def main():
    parser = argparse.ArgumentParser(description='Creates <.csv> of total nt counts and the complement nt counts at position 1 from all reads in <.fastq> files in the directory')
    parser.add_argument('-single',default = None, help = 'Operate on this single file, rather than the directory')
    args = parser.parse_args()
    
    if args.single != None:
        outname = args.single.replace('.fastq','_ligation_bias.csv')
        report_lig(args.single, outname)
    else:
        for fyle in sorted(glob.glob('*.fastq')):
            outname = fyle.replace('.fastq','_ligation_bias.csv')
            report_lig(fyle, outname)

if __name__ == '__main__':
    main()
