#!/usr/bin/env python

#Imports
import glob
from collections import Counter
import argparse
from itertools import islice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def read_in_fasta(fasta_fyle):
    '''Reads in a fasta file to a dictionary'''
    fasta_dict = {}
    fasta_sequences,fasta_dict = SeqIO.parse(open(fasta_fyle),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

def read_in_rtsc(rtsc_fyle):
    '''Reads in a reactivity file to a dictionary'''
    information = {}
    with open(rtsc_fyle,'r') as f:
        while True:
            next_n_lines = list(islice(f, 3))
            if not next_n_lines:
                break
            transcript,stops,empty_line = [n.strip() for n in next_n_lines]
            information[transcript] = [int(x) for x in stops.split('\t')]
    return information

def coverage_dictionary_generic(sequence_dictionary,rtstop_dictionary,specificity='AC'):
    '''Calculates coverage on each transcript, returns a dictionary of transcript to coverage value with user defined specificity'''
    coverage_dict = {}
    for transcript,stops in rtstop_dictionary.items():
        seq = sequence_dictionary[transcript].upper()[:-1]
        seq_counter = Counter(seq)
        target_counts = sum([seq_counter[letter] for letter in specificity])
        matched = zip(list(seq),stops[1:])
        target_hits = sum([derp[1] for derp in matched if derp[0] in specificity])
        try:
            coverage_dict[transcript] = float(target_hits)/target_counts
        except ZeroDivisionError:
            coverage_dict[transcript] = 0
    return coverage_dict

def batch_coverage_dictionary(fasta_fyle,batch_specificity):
    '''generates a dictionary for each file'''
    grande_dictionary = {}
    sequence_dictionary = read_in_fasta(fasta_fyle)
    for rtsc_fyle in sorted(glob.glob('*.rtsc')):
        rtsc_dict = read_in_rtsc(rtsc_fyle)
        grande_dictionary[rtsc_fyle.replace('.rtsc','')] = coverage_dictionary_generic(sequence_dictionary,rtsc_dict,batch_specificity)
    return grande_dictionary

def write_single_coverage(coverage_dict,infyle):
    '''Writes out in csv format'''
    outfile = infyle.replace('.rtsc','')+'_coverage.csv'
    header = ','.join(['transcript',infyle.replace('.rtsc','')+'_coverage'])
    with open(outfile,'w') as g:
        g.write(header+'\n')
        for transcript, value in sorted(coverage_dict.items(), key=lambda x: x[1],reverse=True):
            g.write(','.join([transcript,str(value)])+'\n')

def write_batch_coverage(batch_coverage_dictionary):
    '''Write out all the information based on alphabetical order'''
    outfile = '_'.join(sorted(batch_coverage_dictionary.keys()))+'_coverage.csv'
    header = ','.join(['transcript']+[name+'_coverage' for name in sorted(batch_coverage_dictionary.keys())])
    all_keys = [z.keys() for z in batch_coverage_dictionary.values()]
    key_set = set([j for i in all_keys for j in i])
    with open(outfile,'w') as g:
        g.write(header+'\n')
        for transcript in key_set:
            outline = [str(batch_coverage_dictionary[name][transcript]) if transcript in batch_coverage_dictionary[name] else 'NA' for name in sorted(batch_coverage_dictionary.keys())]
            g.write(','.join([transcript]+outline)+'\n')

def main():
    parser = argparse.ArgumentParser(description='Creates <.csv> of stop coverages from <.rtsc> files in the directory')
    parser.add_argument("index",type=str,help="<.fasta> file used to generate the <.rtsc>")
    parser.add_argument('-single',default = None, help = 'Operate on this single file, rather than the directory')
    parser.add_argument('-bases',type=str,default='AC', help='[default = AC] Coverage Specificity')
    args = parser.parse_args()
    if args.single != None:
        input_fasta = read_in_fasta(args.index)
        input_rtsc = read_in_rtsc(args.single)
        single_data = coverage_dictionary_generic(input_fasta,input_rtsc,args.bases)
        write_single_coverage(single_data,args.single)
    else:
        batch_coverage = batch_coverage_dictionary(args.index,args.bases)
        write_batch_coverage(batch_coverage)
        

if __name__ == '__main__':
    main()

