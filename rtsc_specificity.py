#!/usr/bin/env python

'''
By default, this will compare every <.rtsc> in the directory to the -index fasta file, and generate a small report.
A custom set of files can instead be given to the program.
'''

#Imports
import glob
import argparse
from collections import Counter
from itertools import islice
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

#Functions
def read_in_rtsc(rtsc_fyle):
    '''Reads a <.rtsc> file into a dictionary'''
    information = {}
    with open(rtsc_fyle,'r') as f:
        while True:
            next_n_lines = list(islice(f, 3))
            if not next_n_lines:
                break
            transcript,stops,empty_line = [n.strip() for n in next_n_lines]
            information[transcript] = [int(x) for x in stops.split('\t')]
    return information

def read_in_fasta(fasta_fyle):
    '''Reads in a fasta file to a dictionary'''
    fasta_dict = {}
    fasta_sequences,fasta_dict = SeqIO.parse(open(fasta_fyle),'fasta'),{}
    for fasta in fasta_sequences:
        fasta_dict[fasta.id] = str(fasta.seq)
    return fasta_dict

def rtsc_to_specificity(rtsc_fyle,fasta_index):
    '''Returns a Counter object of the frequencies of the nucleotide specificities'''
    Xcounter = Counter()
    rtsc_dictionary = read_in_rtsc(rtsc_fyle)
    for transcript, stops in rtsc_dictionary.items():
        seq = fasta_index[transcript].upper()[:-1]
        matched_values = zip(list(seq),stops[1:])
        for pair in matched_values:
            Xcounter[pair[0]]+=int(pair[1])
    return Xcounter
            
def write_specificity_data(infod,outfyle='specificity.csv'):
    '''Dumps to a .csv'''
    #x_names = [x.strip('.rtsc') for x in infod.keys()]
    x_names = [x.replace('.rtsc','') for x in infod.keys()]
    x_labs_1 = [y+'_count' for y in x_names]
    x_labs_2 = [z+'_specificity' for z in x_names]
    x_header = ','.join(['base']+sorted(x_labs_1+x_labs_2))
    A_totals = ['A']+[(infod[q]['A'],infod[q]['A']/float(sum(infod[q].values()))) for q in sorted(infod.keys())]
    T_totals = ['T']+[(infod[q]['T'],infod[q]['T']/float(sum(infod[q].values()))) for q in sorted(infod.keys())]
    G_totals = ['G']+[(infod[q]['G'],infod[q]['G']/float(sum(infod[q].values()))) for q in sorted(infod.keys())]
    C_totals = ['C']+[(infod[q]['C'],infod[q]['C']/float(sum(infod[q].values()))) for q in sorted(infod.keys())]
    x_data = [A_totals,T_totals,G_totals,C_totals]
    with open(outfyle,'w') as g:
        g.write(x_header+'\n')
        for data in x_data:
            g.write(','.join(str(y) for x in data for y in x)+'\n')

def main():
    parser = argparse.ArgumentParser(description='Analyzes native/reagent nucleotide stop specificity')
    parser.add_argument("-index",type=str,help="<.fasta> file used to generate the <.rtsc>")
    parser.add_argument('-rtsc',default = None, help='Operate on specific <.rtsc>', nargs='+')
    parser.add_argument('-name',default = None, help='Specify output file name')
    args = parser.parse_args()
    
    #Read in fasta sequences
    input_fasta = read_in_fasta(args.index)

    #Batch Mode
    if args.rtsc == None:
        #Generate Out Name
        #default_name = '_'.join(sorted([x.strip('.rtsc') for x in glob.glob('*.rtsc')]))+'_specificity.csv'
        default_name = '_'.join(sorted([x.replace('.rtsc','') for x in glob.glob('*.rtsc')]))+'_specificity.csv'
        out_name = default_name if args.name == None else args.name
        collected_info = {}
        for fyle in glob.glob('*.rtsc'):
            collected_info[fyle] = rtsc_to_specificity(fyle,input_fasta)
        write_specificity_data(collected_info,out_name)
    
    #Specific Mode
    if args.rtsc != None:
        #Generate Out Name
        #default_name = '_'.join(sorted([x.strip('.rtsc') for x in args.rtsc]))+'_specificity.csv'
        default_name = '_'.join(sorted([x.replace('.rtsc','') for x in args.rtsc]))+'_specificity.csv'
        out_name = default_name if args.name == None else args.name
        collected_info = {}
        for fyle in args.rtsc:
            collected_info[fyle] = rtsc_to_specificity(fyle,input_fasta)
        write_specificity_data(collected_info,out_name)



if __name__ == '__main__': 
    main()

