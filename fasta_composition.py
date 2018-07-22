#!/usr/bin/env python

'''
This script will generate a convenient <.csv> detailing the nucleotide composition of 
sequences in fasta format. This <.csv> may be imported and merged with other <.csv>s to build 
a comprehensive data sheet.
'''

#Imports
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import argparse
import glob

#Functions
def read_in_fasta(afasta):
    '''Reads in a fasta file to a dictionary'''
    fasta_dict = {}
    fasta_sequences,fasta_dict = SeqIO.parse(open(afasta),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

def generate_comp(fasta_dictonary):
    '''generates statistics for a fasta dictionary'''
    out = {}
    for name,sequence in fasta_dictonary.items():
        nucleotide_count = Counter(sequence)
        total_nucleotides = sum(nucleotide_count.values())
        GC_content = (nucleotide_count['G']+nucleotide_count['C'])/float(total_nucleotides)
        AT_content = (nucleotide_count['A']+nucleotide_count['T'])/float(total_nucleotides)
        AC_content = (nucleotide_count['A']+nucleotide_count['C'])/float(total_nucleotides)
        A_content = (nucleotide_count['A'])/float(total_nucleotides)
        C_content = (nucleotide_count['C'])/float(total_nucleotides)
        G_content = (nucleotide_count['G'])/float(total_nucleotides)
        T_content = (nucleotide_count['T'])/float(total_nucleotides)
        entire_length = len(sequence)
        out[name] = [name,GC_content,AT_content,AC_content,A_content,C_content,G_content,T_content,entire_length]
    return out
        
def write_out_csv(composition_dictionary,out_fyle='fasta_comp.csv'):
    '''Writes out the data as a csv. This csv may be combined with others based on transcript'''
    with open(out_fyle,'w') as g:
        #Write Header
        g.write(','.join(['transcript','GC_content','AT_content','AC_content',
                          'A_content','C_content','G_content','T_content','length'])+'\n')
        #Write Data
        for data in composition_dictionary.values():
            g.write(','.join([str(x) for x in data])+'\n')

def main():
    parser = argparse.ArgumentParser(description='Creates <.csv> of sequence compositions from <.fasta> files in the directory')
    parser.add_argument('-single',default = None, help = 'Operate on this single file')
    parser.add_argument('-suffix',type=str,default='composition', help='[default = composition] <.csv> file suffix')
    args = parser.parse_args()
    if args.single != None:
        sequences = read_in_fasta(args.single)
        composition = generate_comp(sequences)
        new_name = args.single.split('.')[0]+'_'+args.suffix+'.csv'
        write_out_csv(composition,new_name)
    else:
        all_fastas = [derp for herp in [glob.glob(x) for x in ['*.fa','*.fasta','*.fas']] for derp in herp]
        for fyle in all_fastas:
            sequences = read_in_fasta(fyle)
            composition = generate_comp(sequences)
            new_name = fyle.split('.')[0]+'_'+args.suffix+'.csv'
            write_out_csv(composition,new_name)


if __name__ == '__main__':
    main()
